# ==============================================================================
# Duke Pipeline - Alignment Functions
# ==============================================================================
# Functions for reference processing and minimap2 alignment
# ==============================================================================

# Parse reference sequence - split at NNNNN into left/right flanks
parse_reference_sequence <- function(ref_path) {
  # Import reference sequence
  ref_sequence <- readDNAStringSet(ref_path, format = "fasta")
  ref_name <- names(ref_sequence)[[1]]
  
  # Find 'NNNNN' split point
  split_point <- stringr::str_locate(as.character(ref_sequence), "N{5,}")[1,]
  
  if (is.na(split_point[1])) {
    stop("Reference sequence must contain at least 5 consecutive Ns (NNNNN) to mark the repeat region")
  }
  
  # Split reference sequence
  ref_sequence_left <- substr(ref_sequence, 1, split_point[1] - 1)
  ref_sequence_right <- substr(ref_sequence, split_point[2] + 1, nchar(ref_sequence))
  
  # Create named DNAStringSet
  ref_sequences <- DNAStringSet(setNames(c(ref_sequence_left, ref_sequence_right),
                                         c(paste0(ref_name, "_left"), paste0(ref_name, "_right"))))
  
  # Summarise reference sequences
  ref_summary <- tibble::tibble(
    rname = names(ref_sequences),
    ref_id = sub("_(left|right)$", "", names(ref_sequences)),
    ref_flank = sub("^.*_(left|right)$", "\\1", names(ref_sequences)),
    ref_length = width(ref_sequences),
    ref_sequence = as.character(ref_sequences))
  
  return(list(
    sequences = ref_sequences,
    summary = ref_summary
  ))
}


# Run minimap2 alignment for one sample
minimap2_align_sample <- function(sample_name,
                                  sequences,
                                  ref_sequences,
                                  out_dir,
                                  minimap2_args = "-x map-ont",
                                  resume = FALSE) {
  
  # Extract DNA stringset for this sample
  dna_stringset <- sequences[[sample_name]]
  
  # Exit if empty
  if (is.null(dna_stringset) || length(dna_stringset) == 0) {
    return(NULL)
  }
  
  # Check for cached result
  alignment_info_path <- file.path(out_dir, paste0("align_", sample_name, ".RData"))
  if (resume && file.exists(alignment_info_path)) {
    load(alignment_info_path)
    return(alignment_info)
  }
  
  # Vector of ref names
  ref_names <- names(ref_sequences)
  
  # Align to each reference flank
  alignment_list <- lapply(ref_names, function(ref_name) {
    
    # Define temporary file stem
    temp_stem <- file.path(out_dir, paste0(sample_name, "-", ref_name))
    
    # Check if any sequences exceed FASTQ write limit (200,000 characters)
    # This is a limitation of Biostrings::writeQualityScaledXStringSet
    max_seq_length <- max(width(dna_stringset))
    use_fastq <- is(dna_stringset, "QualityScaledXStringSet") && max_seq_length <= 200000
    
    # Warn if switching to FASTA due to long reads
    if (is(dna_stringset, "QualityScaledXStringSet") && max_seq_length > 200000) {
      warning("Sample ", sample_name, " has very long reads (max: ", max_seq_length, 
              " bp). Switching to FASTA format (quality scores will be lost).")
    }
    
    # Export query sequences
    if (use_fastq) {
      query_path <- paste0(temp_stem, ".fastq")
      writeQualityScaledXStringSet(dna_stringset, filepath = query_path)
    } else {
      query_path <- paste0(temp_stem, ".fasta")
      # Convert to DNAStringSet if it's QualityScaledXStringSet (drops quality scores)
      if (is(dna_stringset, "QualityScaledXStringSet")) {
        dna_stringset_plain <- DNAStringSet(as.character(dna_stringset))
        names(dna_stringset_plain) <- names(dna_stringset)
        writeXStringSet(dna_stringset_plain, filepath = query_path)
      } else {
        writeXStringSet(dna_stringset, filepath = query_path)
      }
    }
    
    # Export reference
    ref_path <- paste0(temp_stem, "-ref.fasta")
    writeXStringSet(ref_sequences[ref_name], ref_path)
    
    # Parse minimap2 args
    args_user <- if (is.na(minimap2_args) || minimap2_args == "") {
      c("-x", "map-ont")
    } else {
      strsplit(minimap2_args, "\\s+")[[1]]
    }
    
    # Handle -ax -> -x substitution
    if ("-ax" %in% args_user) {
      ax_index <- which(args_user == "-ax")
      args_user[ax_index] <- "-x"
    }
    
    # Required arguments
    args_required <- c("-a", "--MD", "--cs")
    args_user <- args_user[!args_user %in% args_required]
    
    # Full minimap2 command
    minimap2_command <- c(args_required, args_user, ref_path, query_path)
    
    # SAM and BAM paths
    sam_path <- paste0(temp_stem, ".sam")
    bam_path <- paste0(temp_stem, ".bam")
    
    # Run minimap2
    minimap2_result <- system2("minimap2", args = minimap2_command, 
                               stdout = sam_path, stderr = NULL)
    
    # Check if SAM file was created
    if (!file.exists(sam_path) || file.size(sam_path) == 0) {
      warning("minimap2 failed or produced empty output for ", sample_name, "-", ref_name)
      # Clean up
      temp_files <- c(query_path, ref_path, sam_path)
      file.remove(temp_files[file.exists(temp_files)])
      return(NULL)
    }
    
    # Sort SAM to BAM (using old syntax that works on macOS)
    system(paste("samtools sort", sam_path, ">", bam_path))
    
    # Check if BAM was created
    if (!file.exists(bam_path) || file.size(bam_path) == 0) {
      warning("samtools sort failed for ", sample_name, "-", ref_name)
      # Clean up
      temp_files <- c(query_path, ref_path, sam_path, bam_path)
      file.remove(temp_files[file.exists(temp_files)])
      return(NULL)
    }
    
    # Index BAM
    system(paste("samtools index", bam_path))
    
    # Import BAM
    tryCatch({
      param <- ScanBamParam(what = scanBamWhat(), tag = c("NM", "MD"))
      bam_data <- scanBam(bam_path, param = param)
      alignments <- bam_data[[1]]
      
      # Expand tags to match length
      if (!is.null(alignments$tag)) {
        for (tag_name in names(alignments$tag)) {
          if (length(alignments$tag[[tag_name]]) < length(alignments$qname)) {
            alignments$tag[[tag_name]] <- c(
              alignments$tag[[tag_name]], 
              rep(NA, length(alignments$qname) - length(alignments$tag[[tag_name]]))
            )
          }
        }
      }
      
      # Ensure NM tag exists
      if (is.null(alignments$tag$NM)) {
        alignments$tag$NM <- rep(NA, length(alignments$qname))
      }
      
      # Convert to data frame
      alignment_df <- as.data.frame(alignments) %>%
        dplyr::rename(ref_start = pos, 
                      mismatches = tag.NM,
                      aligned_subsequence = seq)
      
      # Clean up temp files
      temp_files <- c(query_path, ref_path, sam_path, bam_path, paste0(bam_path, ".bai"))
      file.remove(temp_files[file.exists(temp_files)])
      
      return(alignment_df)
      
    }, error = function(e) {
      warning("Failed to read BAM file for ", sample_name, "-", ref_name, ": ", e$message)
      # Clean up
      temp_files <- c(query_path, ref_path, sam_path, bam_path, paste0(bam_path, ".bai"))
      file.remove(temp_files[file.exists(temp_files)])
      return(NULL)
    })
  })
  
  # Combine all alignments (remove NULLs from failed references)
  alignment_list <- alignment_list[!sapply(alignment_list, is.null)]
  
  # If all alignments failed, return empty tibble
  if (length(alignment_list) == 0) {
    alignment_info <- tibble::tibble(
      qname = character(), rname = character(), strand = character(),
      ref_start = integer(), mapq = integer(), mismatches = integer()
    )
    save(alignment_info, file = alignment_info_path)
    return(alignment_info)
  }
  
  # Combine all alignments
  alignment_info <- data.table::rbindlist(alignment_list, fill = TRUE)
  
  # Remove unaligned reads
  alignment_info <- alignment_info %>%
    dplyr::filter(!is.na(rname))
  
  # If no alignments, return empty tibble
  if (nrow(alignment_info) == 0) {
    alignment_info <- tibble::tibble(
      qname = character(), rname = character(), strand = character(),
      ref_start = integer(), mapq = integer(), mismatches = integer()
    )
    save(alignment_info, file = alignment_info_path)
    return(alignment_info)
  }
  
  # Add query sequences
  query_seqs <- data.frame(
    qname = names(dna_stringset),
    query_sequence = as.character(dna_stringset),
    stringsAsFactors = FALSE
  )
  alignment_info <- alignment_info %>%
    left_join(query_seqs, by = "qname")
  
  # Calculate alignment metrics from CIGAR
  cigar_rle <- cigarToRleList(alignment_info$cigar)
  
  # Alignment length on reference
  alignment_info$alignment_length_on_ref <- vapply(cigar_rle, function(rle) {
    sum(runLength(rle)[runValue(rle) %in% c("M", "D", "N", "=", "X")])
  }, integer(1))
  
  # Calculate ref_end (standard genomic coordinates: ref_start < ref_end always)
  alignment_info$ref_end <- alignment_info$ref_start + alignment_info$alignment_length_on_ref - 1
  
  # Insertion/deletion counts
  alignment_info$insertions <- vapply(cigar_rle, function(rle) {
    sum(runLength(rle)[runValue(rle) == "I"])
  }, integer(1))
  alignment_info$deletions <- vapply(cigar_rle, function(rle) {
    sum(runLength(rle)[runValue(rle) == "D"])
  }, integer(1))
  
  # Alignment length on query
  alignment_info$alignment_length_on_query <- vapply(cigar_rle, function(rle) {
    sum(runLength(rle)[runValue(rle) %in% c("M", "I", "=", "X")])
  }, integer(1))
  
  # Soft-clip at start
  alignment_info$softclip_start <- vapply(cigar_rle, function(rle) {
    if (length(rle) > 0 && runValue(rle)[1] == "S") {
      as.integer(runLength(rle)[1])
    } else {
      0L
    }
  }, integer(1))
  
  # Query coordinates (simple CIGAR-based initially)
  alignment_info$query_start <- alignment_info$softclip_start + 1
  alignment_info$query_end <- alignment_info$query_start + alignment_info$alignment_length_on_query - 1
  
  # Clean up intermediate columns and reorder
  alignment_info <- alignment_info %>%
    dplyr::select(-softclip_start) %>%
    dplyr::select(qname, rname, strand, 
                  ref_start, ref_end, alignment_length_on_ref,
                  query_start, query_end, alignment_length_on_query,
                  mapq, mismatches, insertions, deletions,
                  query_sequence, aligned_subsequence, cigar) %>%
    dplyr::arrange(qname, desc(mapq))
  
  # Save alignment info
  save(alignment_info, file = alignment_info_path)
  
  return(alignment_info)
}
