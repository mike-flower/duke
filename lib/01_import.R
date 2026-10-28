# ==============================================================================
# Duke Pipeline - Import Functions
# ==============================================================================
# Functions for importing FASTQ, FASTA, and BAM files
# ==============================================================================

# Function to import fastq including quality data
fastq_import <- function(path) {
  readfastq_file <- readFastq(path)
  read_names <- as.character(ShortRead::id(readfastq_file))
  # Strip tags from read names by keeping only the first field before the first space
  clean_read_names <- sapply(strsplit(read_names, " "), `[`, 1)
  quality_scaled_dnastringset <- 
      QualityScaledDNAStringSet(sread(readfastq_file), as(quality(readfastq_file), "PhredQuality"))
  names(quality_scaled_dnastringset) <- clean_read_names
  return(quality_scaled_dnastringset)
}

# Function to import fasta
fasta_import <- function(path) {
  readDNAStringSet(path, format = "fasta")
}

# ------------------------------------------------------------------------------
# BAM reading via samtools
# ------------------------------------------------------------------------------
# BAMs are read with `samtools view` rather than Rsamtools::scanBam. Rsamtools
# prints "[samopen] no @SQ lines in the header" from C for every unaligned BAM
# (all PacBio HiFi BAMs), three times per file, which R cannot silence; samtools'
# stderr can simply be discarded. Secondary and supplementary records are
# skipped (-F 0x900); the sequence is the stored SEQ, as scanBam returned it.
# Returns data.table(qname, flag, seq, qual, np, rq, ec, zm, bq); absent tags NA.
.bam_awk <- paste0("'BEGIN{FS=OFS=\"\\t\"} {np=rq=ec=zm=bq=\"NA\"; ",
                   "for(i=12;i<=NF;i++){t=substr($i,1,2); v=substr($i,6); ",
                   "if(t==\"np\")np=v; else if(t==\"rq\")rq=v; else if(t==\"ec\")ec=v; ",
                   "else if(t==\"zm\")zm=v; else if(t==\"bq\")bq=v}; ",
                   "print $1,$2,$10,$11,np,rq,ec,zm,bq}'")

read_bam_records <- function(path, n = NULL) {
  cmd <- paste("samtools view -F 0x900", shQuote(path), "2>/dev/null |",
               if (!is.null(n)) paste("head -n", as.integer(n), "|") else "", "awk", .bam_awk)
  cols <- c("qname", "flag", "seq", "qual", "np", "rq", "ec", "zm", "bq")
  out <- tryCatch(data.table::fread(cmd = cmd, sep = "\t", header = FALSE, quote = "", na.strings = "NA",
                                    colClasses = c("character", "integer", "character", "character",
                                                   "integer", "numeric", "numeric", "integer", "integer"),
                                    col.names = cols, showProgress = FALSE),
                  error = function(e) NULL)
  if (is.null(out) || !nrow(out)) {
    if (is.null(n) && Sys.which("samtools") == "") stop("samtools is not on the PATH (needed to read BAM files)")
    out <- data.table::data.table(qname = character(), flag = integer(), seq = character(), qual = character(),
                                  np = integer(), rq = numeric(), ec = numeric(), zm = integer(), bq = integer())
  }
  out
}

bam_import <- function(path) {
  r <- read_bam_records(path)
  seqs <- Biostrings::DNAStringSet(r$seq)
  if (!nrow(r) || all(r$qual == "*")) { names(seqs) <- r$qname; return(seqs) }
  out <- Biostrings::QualityScaledDNAStringSet(seqs, Biostrings::PhredQuality(r$qual))
  names(out) <- r$qname
  out
}

# Function to import sequencing files
import_sequencing <- function(path, format) {
  if (format %in% c("fastq", "fastq.gz", "fq", "fq.gz")) {
    return(fastq_import(path))
  } else if (format == "fasta" || format == "fa") {
    return(fasta_import(path))
  } else if (format == "bam") {
    return(bam_import(path))
  } else {
    stop("Unsupported file format")
  }
}

# ------------------------------------------------------------------------------
# PacBio HiFi detection and per-read aux-tag extraction (v2.3.0)
# ------------------------------------------------------------------------------

# detect_pacbio: Cheap auto-detection of PacBio HiFi BAMs.
#   Primary check: the np:i: tag on the first record (one-record scan).
#   Secondary check (only on failure): @PG header lines with PN matching
#   ccs / skera / lima (case-insensitive).
#   Returns FALSE for non-BAM formats, missing files, and any error.
#
#   Usage:
#     is_pacbio <- detect_pacbio(path = "/path/to/sample.bam", format = "bam")
detect_pacbio <- function(path, format) {

  # Non-BAM inputs are never PacBio
  if (!is.character(format) || length(format) != 1 || tolower(format) != "bam") {
    return(FALSE)
  }
  if (!is.character(path) || length(path) != 1 || !file.exists(path)) {
    return(FALSE)
  }

  tryCatch({
    # Primary: the np:i: tag on the first record
    first <- read_bam_records(path, n = 1)
    if (nrow(first) && !is.na(first$np[1])) return(TRUE)
    # Secondary: @PG header lines for ccs / skera / lima
    hdr <- suppressWarnings(system(paste("samtools view -H", shQuote(path), "2>/dev/null"), intern = TRUE))
    pn <- unlist(regmatches(hdr[startsWith(hdr, "@PG")], gregexpr("PN:[^\t]+", hdr[startsWith(hdr, "@PG")])))
    any(grepl("^PN:(ccs|skera|lima)$", pn, ignore.case = TRUE))
  }, error = function(e) FALSE)
}


# extract_pacbio_metadata: Single-pass extraction of PacBio HiFi per-read aux
#   tags from a BAM. Returns a data.frame with one row per read and columns
#   read_name, np, rq, ec, zm, bq. Absent tags are NA.
#
#   Usage:
#     md <- extract_pacbio_metadata(path = "/path/to/sample.bam")
extract_pacbio_metadata <- function(path) {
  r <- read_bam_records(path)
  data.frame(read_name = r$qname, np = r$np, rq = r$rq, ec = r$ec, zm = r$zm, bq = r$bq,
             stringsAsFactors = FALSE)
}
