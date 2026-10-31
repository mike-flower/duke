# ==============================================================================
# Duke - alignment to reference flanks
# ==============================================================================
# Each reference contains one block of >=5 Ns marking the repeat region. It is
# split into a left and a right flank, and reads are aligned to each flank
# separately with minimap2, so a read aligns whatever length its repeat is.
#
# Everything for one sample happens in one worker (align -> best alignment per
# flank -> both flanks on the same strand -> left before right -> strand
# correction -> manual exclusions), which keeps the per-read tables small and
# parallelises across samples.
#
# Coordinates: for reads on the reverse strand minimap2 reports query
# coordinates on the reverse complement. Those are exactly the coordinates of
# the strand-corrected read, so the read is reverse-complemented and the
# coordinates are kept as reported.
# ==============================================================================

parse_reference_sequence <- function(ref_path) {
  ref_sequence <- Biostrings::readDNAStringSet(ref_path, format = "fasta")
  if (length(ref_sequence) != 1) stop("Reference must contain exactly one sequence: ", ref_path)
  ref_name <- sub("\\s.*$", "", names(ref_sequence)[[1]])
  full <- as.character(ref_sequence[[1]])
  blocks <- gregexpr("[Nn]{5,}", full)[[1]]
  if (blocks[1] == -1) stop("Reference must contain at least 5 consecutive Ns marking the repeat region: ", ref_path)
  if (length(blocks) > 1) stop("Reference must contain exactly one N block: ", ref_path,
                               "\n  For tracts that are far apart, give each its own reference in the tracts table.")
  n_start <- blocks[1]; n_end <- n_start + attr(blocks, "match.length")[1] - 1
  left  <- substr(full, 1, n_start - 1)
  right <- substr(full, n_end + 1, nchar(full))
  ref_sequences <- Biostrings::DNAStringSet(setNames(c(left, right), c(paste0(ref_name, "_left"), paste0(ref_name, "_right"))))
  list(
    ref_id = ref_name, path = ref_path, full = full, n_start = n_start, n_end = n_end,
    sequences = ref_sequences,
    summary = data.frame(ref_id = ref_name, flank = c("left", "right"), rname = names(ref_sequences),
                         length_bp = Biostrings::width(ref_sequences),
                         first_bases = substr(c(left, right), 1, 20),
                         last_bases = c(substr(left, nchar(left) - 19, nchar(left)), substr(right, nchar(right) - 19, nchar(right))),
                         stringsAsFactors = FALSE)
  )
}

# Leading soft/hard clip length from a CIGAR string (vectorised)
.leading_clip <- function(cigar) {
  m <- regmatches(cigar, regexpr("^([0-9]+[SH])+", cigar))
  out <- integer(length(cigar))
  has <- grepl("^([0-9]+[SH])+", cigar)
  out[has] <- vapply(m, function(x) sum(as.integer(regmatches(x, gregexpr("[0-9]+", x))[[1]])), integer(1))
  out
}

# Align one sample's reads to one flank; returns one row per alignment record
.align_to_flank <- function(query_path, ref_seq, ref_name, stem, minimap2_args, keep_bam, bam_dir) {
  ref_path <- paste0(stem, "-", ref_name, "-ref.fasta")
  Biostrings::writeXStringSet(ref_seq, ref_path)
  sam_path <- paste0(stem, "-", ref_name, ".sam")
  args <- strsplit(minimap2_args %||% "-x sr -w 1", "\\s+")[[1]]
  args[args == "-ax"] <- "-x"
  args <- args[!args %in% c("-a", "--MD", "--cs")]
  status <-system2("minimap2", c("-a", args, ref_path, query_path), stdout = sam_path, stderr = FALSE)
  if (!file.exists(sam_path) || file.size(sam_path) == 0) {
    unlink(c(ref_path, sam_path)); stop("minimap2 produced no output for ", basename(stem), " / ", ref_name)
  }
  bam_path <- paste0(stem, "-", ref_name, ".bam")
  system2("samtools", c("sort", "-o", bam_path, sam_path), stdout = FALSE, stderr = FALSE)
  system2("samtools", c("index", bam_path), stdout = FALSE, stderr = FALSE)
  bam <- Rsamtools::scanBam(bam_path, param = Rsamtools::ScanBamParam(
    what = c("qname", "flag", "rname", "strand", "pos", "mapq", "cigar"), tag = "NM"))[[1]]
  if (keep_bam && !is.null(bam_dir)) {
    dir.create(bam_dir, showWarnings = FALSE, recursive = TRUE)
    file.copy(c(bam_path, paste0(bam_path, ".bai")), bam_dir, overwrite = TRUE)
  }
  unlink(c(ref_path, sam_path, bam_path, paste0(bam_path, ".bai")))
  keep <- !is.na(bam$rname) & bitwAnd(bam$flag, 4L) == 0
  if (!any(keep)) return(data.table::data.table())
  cig <- bam$cigar[keep]
  qstart <- .leading_clip(cig) + 1L
  data.table::data.table(
    qname = bam$qname[keep],
    strand = as.character(bam$strand[keep]),
    rstart = bam$pos[keep],
    rend = bam$pos[keep] + GenomicAlignments::cigarWidthAlongReferenceSpace(cig) - 1L,
    qstart = qstart,
    qend = qstart + GenomicAlignments::cigarWidthAlongQuerySpace(cig, after.soft.clipping = TRUE) - 1L,
    mapq = bam$mapq[keep],
    nm = bam$tag$NM[keep] %||% NA_integer_,
    secondary = bitwAnd(bam$flag[keep], 256L) > 0
  )
}

# Full alignment step for one sample and one reference
process_sample_alignment <- function(sample, reads, ref, out_dir, minimap2_args,
                                     keep_bam = FALSE, bam_dir = NULL, exclusions = NULL) {
  cache <- file.path(out_dir, paste0("align_v31_", ref$ref_id, "_", sample, ".RData"))
  if (file.exists(cache) && safe_load_cache(cache, environment())) return(result)

  n_input <- length(reads)
  empty <- list(table = NULL, funnel = data.frame(file_stem = sample, ref_id = ref$ref_id, input = n_input,
                aligned_left = 0L, aligned_right = 0L, both_flanks = 0L, same_strand = 0L,
                correct_order = 0L, manually_excluded = 0L, final = 0L),
                strand = data.frame(file_stem = sample, ref_id = ref$ref_id, reverse_left = NA_real_, reverse_right = NA_real_,
                                    reverse_both = NA_real_, reverse_final = NA_real_))
  if (n_input == 0) return(empty)

  stem <- file.path(out_dir, paste0(sample, "-", ref$ref_id))
  long_reads <- max(Biostrings::width(reads)) > 200000
  if (is(reads, "QualityScaledXStringSet") && !long_reads) {
    query_path <- paste0(stem, ".fastq")
    Biostrings::writeQualityScaledXStringSet(reads, query_path)
  } else {
    query_path <- paste0(stem, ".fasta")
    plain <- Biostrings::DNAStringSet(reads); names(plain) <- names(reads)
    Biostrings::writeXStringSet(plain, query_path)
  }

  aln <- lapply(c(left = 1, right = 2), function(i) {
    .align_to_flank(query_path, ref$sequences[i], names(ref$sequences)[i], stem, minimap2_args,
                    keep_bam, if (keep_bam) bam_dir else NULL)
  })
  unlink(query_path)

  # Best alignment per read and flank: highest MAPQ, then longest on the reference
  best <- lapply(aln, function(a) {
    if (nrow(a) == 0) return(a)
    a[, rlen := rend - rstart + 1L]
    data.table::setorder(a, qname, -mapq, -rlen, secondary)
    a[!duplicated(qname)]
  })
  n_left <- nrow(best$left); n_right <- nrow(best$right)
  if (n_left == 0 || n_right == 0) { empty$funnel$aligned_left <- n_left; empty$funnel$aligned_right <- n_right; return(empty) }

  both <- merge(best$left, best$right, by = "qname", suffixes = c("_l", "_r"))
  n_both <- nrow(both)
  # Share of reads on the reverse strand at each step, so a strand imbalance in the
  # final reads can be traced to the step that introduced it
  rev_left <- mean(best$left$strand == "-"); rev_right <- mean(best$right$strand == "-")
  rev_both <- if (n_both) mean(both$strand_l == "-") else NA_real_
  both <- both[strand_l == strand_r]
  n_strand <- nrow(both)
  both <- both[qstart_l < qstart_r]
  n_order <- nrow(both)

  n_excl <- 0L
  if (!is.null(exclusions) && sample %in% names(exclusions)) {
    drop <- both$qname %in% exclusions[[sample]]
    n_excl <- sum(drop); both <- both[!drop]
  }

  idx <- match(both$qname, names(reads))
  seqs <- as.character(reads[idx])
  rev <- both$strand_l == "-"
  if (any(rev)) seqs[rev] <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seqs[rev])))
  left_len <- Biostrings::width(ref$sequences)[1]

  tab <- data.table::data.table(
    file_stem = sample, qname = both$qname, ref_id = ref$ref_id,
    strand = ifelse(rev, "+c", "+"),
    query_sequence = seqs, query_length = nchar(seqs),
    left_qstart = both$qstart_l, left_qend = both$qend_l,
    right_qstart = both$qstart_r, right_qend = both$qend_r,
    left_rstart = both$rstart_l, left_rend = both$rend_l,
    right_rstart = both$rstart_r, right_rend = both$rend_r,
    left_mapq = both$mapq_l, right_mapq = both$mapq_r,
    left_nm = both$nm_l, right_nm = both$nm_r,
    # How far each flank alignment stops short of the repeat-adjacent end of the flank.
    # Non-zero values mean the aligner placed the flank boundary away from the N block.
    left_end_shortfall = left_len - both$rend_l,
    right_start_shortfall = both$rstart_r - 1L
  )
  result <- list(table = tab, funnel = data.frame(
    file_stem = sample, ref_id = ref$ref_id, input = n_input, aligned_left = n_left, aligned_right = n_right,
    both_flanks = n_both, same_strand = n_strand, correct_order = n_order,
    manually_excluded = n_excl, final = nrow(tab)),
    strand = data.frame(file_stem = sample, ref_id = ref$ref_id, reverse_left = rev_left, reverse_right = rev_right,
                        reverse_both = rev_both, reverse_final = if (nrow(tab)) mean(rev) else NA_real_))
  save(result, file = cache)
  result
}

# ------------------------------------------------------------------------------
# Manual read exclusions: named list file_stem -> read names
# ------------------------------------------------------------------------------
load_manual_read_exclusions <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path)) return(NULL)
  if (!file.exists(path)) stop("Manual read exclusion file not found: ", path)
  ext <- tolower(tools::file_ext(path))
  x <- if (ext == "csv") utils::read.csv(path, stringsAsFactors = FALSE)
       else if (ext %in% c("tsv", "txt")) utils::read.delim(path, stringsAsFactors = FALSE)
       else {
         sheets <- readxl::excel_sheets(path)
         as.data.frame(readxl::read_excel(path, sheet = if ("Exclusions" %in% sheets) "Exclusions" else 1))
       }
  if (!all(c("file_name", "read_name") %in% names(x))) stop("Exclusion file needs columns file_name and read_name")
  x$file_name <- sub("\\.(bam|fastq|fq|fasta|fa)(\\.gz)?$", "", x$file_name)
  split(as.character(x$read_name), x$file_name)
}
