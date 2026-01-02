# ==============================================================================
# Duke Pipeline - Consensus Functions
# ==============================================================================
# Functions for generating consensus sequences from clustered reads
# ==============================================================================

# Generate consensus sequence from a DNAStringSet
consensus_function <- function(dna_stringset, threads = 1, consensus_threshold = 0.5, sample_size = NA) {
  
  # Optional downsample
  if (!is.na(sample_size)) {
    dna_stringset <- downsample(dna_stringset, sample_size)
  }
  
  # Read count
  read_count <- length(dna_stringset)
  
  # Default outputs
  alignment <- NULL
  consensus_clean <- NA
  
  # Generate consensus sequence if there is more than one sequence
  if (read_count > 1) {
    # Suppress DECIPHER progress bars
    alignment <- suppressMessages(AlignSeqs(dna_stringset, processors = threads, verbose = FALSE))
    consensus <- suppressMessages(ConsensusSequence(alignment, 
                                   threshold = consensus_threshold, 
                                   ambiguity = FALSE,
                                   noConsensusChar = "N",
                                   includeNonLetters = FALSE,
                                   includeTerminalGaps = FALSE)) 
    # threshold sets how much disagreement is allowed at each position.
    # - threshold = 0.05 → strict: ≥95% of reads must agree to call a base.
    # - threshold = 0.5  → majority rule: ≥50% agreement required.
    # Higher values make consensus more permissive, lower values more conservative.
    
    # ambiguity controls how input ambiguity codes (e.g. "R" = A or G) are interpreted.
    # - FALSE = treat them as literal characters (e.g. "R" is counted as "R").
    # - TRUE  = split them between possible bases (e.g. "R" = 0.5 A + 0.5 G).
    
    # Clean the consensus string:
    # - Remove alignment gaps ("-")
    # - Replace remaining ambiguity characters with "N"
    consensus_clean <- gsub("-", "", as.character(consensus))
    consensus_clean <- gsub("[^ACGT]", "N", consensus_clean)
    
  } else if (read_count == 1) {
    consensus_clean <- as.character(dna_stringset[[1]])
  } else {
    consensus_clean <- NA
  }
  
  # Output
  return(list(consensus = consensus_clean,
              input_seqs = dna_stringset,
              alignment = alignment))
}

# Call variants between query and reference
variant_call <- function(query_dnastring, ref_dnastring) {
  
  # Pairwise alignment
  my_alignment <- pairwiseAlignment(pattern = query_dnastring,
                                    subject = ref_dnastring,
                                    substitutionMatrix = NULL,
                                    gapOpening = 5,
                                    gapExtension = 2,
                                    type = "global")
  
  # Extract alignment
  ref_bases <- strsplit(as.character(alignedSubject(my_alignment)), "")[[1]]
  query_bases <- strsplit(as.character(alignedPattern(my_alignment)), "")[[1]]
  
  # Table of variants
  variant_df <- tibble(position = seq_along(ref_bases),
                       ref = ref_bases,
                       query = query_bases) %>%
    mutate(
      variant_type = case_when(
        ref == query ~ "match",
        ref == "-" ~ "insertion",
        query == "-" ~ "deletion",
        TRUE ~ "mismatch"
      )
    ) %>%
    filter(variant_type != "match")  # Only keep variants
  
  # Output
  return(variant_df)
}

# Call variants for all consensus sequences against reference
call_variants_for_clusters <- function(consensus_summary, reference_sequence) {
  
  if (is.null(consensus_summary) || nrow(consensus_summary) == 0) {
    return(NULL)
  }
  
  # Convert reference to DNAString
  ref_dna <- DNAString(reference_sequence)
  
  # Call variants for each consensus
  variant_results <- list()
  
  for (i in 1:nrow(consensus_summary)) {
    sample_name <- consensus_summary$file_stem[i]
    cluster_num <- consensus_summary$cluster_number[i]
    consensus_seq <- consensus_summary$consensus_full[i]
    
    if (!is.na(consensus_seq) && consensus_seq != "") {
      query_dna <- DNAString(consensus_seq)
      variants <- variant_call(query_dna, ref_dna)
      
      if (nrow(variants) > 0) {
        variants$file_stem <- sample_name
        variants$cluster_number <- cluster_num
        variant_results[[length(variant_results) + 1]] <- variants
      }
    }
  }
  
  # Combine all variants
  if (length(variant_results) > 0) {
    all_variants <- bind_rows(variant_results)
    return(all_variants)
  } else {
    return(NULL)
  }
}

# Convert alignment to matrix format
alignment_to_matrix <- function(alignment, include_names = TRUE) {
  if (is.null(alignment) || length(alignment) == 0) return(NULL)
  
  # Convert to matrix of characters
  matrix_df <- alignment %>%
    as.character() %>%
    strsplit("") %>%
    do.call(rbind, .) %>%
    as.data.frame(stringsAsFactors = FALSE)
  
  # Add column names for positions
  colnames(matrix_df) <- paste0("pos", seq_len(ncol(matrix_df)))
  
  # Optionally add read names
  if (include_names) {
    matrix_df$read_name <- names(alignment)
    matrix_df <- dplyr::relocate(matrix_df, read_name)
  }
  
  return(matrix_df)
}

# Export sequences to FASTA
write_fasta <- function(named_sequences, sequence_type, output_path) {
  
  # Replace NA or empty sequences with a placeholder
  placeholder <- "-"
  named_sequences[is.na(named_sequences) | named_sequences == ""] <- placeholder
  
  # Clean sequences to remove invalid characters
  if (sequence_type == "DNA") {
    named_sequences <- setNames(gsub("[^ACGTN\\-]", "-", named_sequences), 
                                names(named_sequences))
    sequence_stringset <- DNAStringSet(named_sequences)
  } else if (sequence_type == "protein") {
    named_sequences <- setNames(gsub("[^ACDEFGHIKLMNPQRSTVWY\\-]", "-", named_sequences), 
                                names(named_sequences))
    sequence_stringset <- AAStringSet(named_sequences)
  } else {
    stop("Invalid sequence type: should be either 'DNA' or 'protein'")
  }
  
  # Export fasta
  writeXStringSet(sequence_stringset, 
                  filepath = output_path, 
                  format = "fasta")
}

# Generate summary of variants per cluster
summarise_variants <- function(variant_df) {
  
  if (is.null(variant_df) || nrow(variant_df) == 0) {
    return(NULL)
  }
  
  variant_summary <- variant_df %>%
    group_by(file_stem, cluster_number, variant_type) %>%
    summarise(
      count = n(),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = variant_type, 
                values_from = count, 
                values_fill = 0)
  
  return(variant_summary)
}
