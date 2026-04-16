# ==============================================================================
# Duke Pipeline - Alignment Processing Functions
# ==============================================================================
# Functions for post-alignment processing, filtering, and quality control
# Used primarily in Module 2 (Alignment and Processing)
# ==============================================================================

# ==============================================================================
# Strand Correction
# ==============================================================================

# Reverse complement a DNA sequence
reverse_complement <- function(sequence) {
  if (is.na(sequence) || sequence == "") {
    return(sequence)
  }
  as.character(reverseComplement(DNAString(sequence)))
}


# Correct strand orientation for one sample
correct_strand_orientation <- function(sample_name,
                                       alignment_list,
                                       out_dir,
                                       resume = FALSE) {
  
  # Check cache
  cache_path <- file.path(out_dir, paste0("strand_correct_", sample_name, ".RData"))
  if (resume && file.exists(cache_path)) {
    load(cache_path)
    return(alignment_corrected)
  }
  
  # Extract alignment data
  alignment_info <- alignment_list[[sample_name]]
  
  # Exit if NULL
  if (is.null(alignment_info)) {
    return(NULL)
  }
  
  # Exit if empty
  if (nrow(alignment_info) == 0) {
    return(alignment_info)
  }
  
  # Add query_length if not present
  if (!"query_length" %in% names(alignment_info)) {
    alignment_info <- alignment_info %>%
      dplyr::mutate(query_length = nchar(query_sequence))
  }
  
  # Correct sequence and coordinates for reverse strand reads (vectorized)
  alignment_corrected <- alignment_info

  # Convert strand from factor to character first — assigning "+c" to a factor
  # without that level would silently produce NAs
  alignment_corrected$strand <- as.character(alignment_corrected$strand)

  rev_idx <- which(alignment_corrected$strand == "-")

  if (length(rev_idx) > 0) {
    # Reverse complement query_sequence for all reverse-strand reads in one C-level call
    alignment_corrected$query_sequence[rev_idx] <-
      as.character(reverseComplement(DNAStringSet(alignment_corrected$query_sequence[rev_idx])))

    # Reverse complement aligned_subsequence where it exists
    has_subseq <- rev_idx[!is.na(alignment_corrected$aligned_subsequence[rev_idx])]
    if (length(has_subseq) > 0) {
      alignment_corrected$aligned_subsequence[has_subseq] <-
        as.character(reverseComplement(DNAStringSet(alignment_corrected$aligned_subsequence[has_subseq])))
    }

    # For reverse strand: minimap2 coordinates are already in RC space, so
    # they remain valid for the RC'd sequence — no coordinate transformation needed.

    # Mark corrected reads
    alignment_corrected$strand[rev_idx] <- "+c"
  }
  
  # Save cache
  save(alignment_corrected, file = cache_path)
  
  return(alignment_corrected)
}


# ==============================================================================
# Alignment Filtering
# ==============================================================================

# Filter alignments for reads spanning both flanks
filter_dual_flank_alignments <- function(alignment_info, ref_summary) {
  
  # Exit if NULL or empty
  if (is.null(alignment_info) || nrow(alignment_info) == 0) {
    return(alignment_info)
  }
  
  # Add ref_id and ref_flank if not present
  if (!"ref_id" %in% names(alignment_info)) {
    alignment_info <- alignment_info %>%
      left_join(ref_summary %>% dplyr::select(rname, ref_id, ref_flank),
                by = "rname")
  }
  
  # Count unique reference sequences per read per ref_id
  ref_pair_count <- alignment_info %>%
    dplyr::filter(!is.na(rname)) %>%
    group_by(file_stem, qname, ref_id) %>%
    dplyr::summarise(n_flanks = n_distinct(rname), .groups = "drop")
  
  # Keep only reads that align to both flanks
  reads_both_flanks <- ref_pair_count %>%
    dplyr::filter(n_flanks == 2) %>%
    dplyr::select(file_stem, qname, ref_id)
  
  # Filter alignment data
  alignment_filtered <- alignment_info %>%
    semi_join(reads_both_flanks, by = c("file_stem", "qname", "ref_id"))
  
  # Keep best alignment per reference flank
  # FIXED: Use explicit filtering and slice_max with with_ties=FALSE
  alignment_filtered <- alignment_filtered %>%
    group_by(file_stem, qname, rname) %>%
    dplyr::filter(mapq == max(mapq)) %>%  # Keep only alignments with max MAPQ
    dplyr::slice_max(alignment_length_on_ref, n = 1, with_ties = FALSE) %>%  # Pick longest, break ties
    ungroup()
  
  # Ensure both flanks are on same strand
  alignment_filtered <- alignment_filtered %>%
    group_by(file_stem, qname, ref_id) %>%
    dplyr::filter(n_distinct(strand) == 1) %>%
    ungroup()
  
  return(alignment_filtered)
}


# ==============================================================================
# Validate Correct Reference Order (AFTER strand correction)
# ==============================================================================
# After strand correction, ALL reads should have left→right reference order
# This validation checks for the biologically expected order and removes
# any reads with backwards or unusual patterns (spurious alignments)

# Validate that alignments follow correct left→right order after correction
validate_correct_reference_order <- function(alignment_corrected_list) {
  
  # Merge all alignment data
  alignment_merge <- data.table::rbindlist(alignment_corrected_list, fill = TRUE)
  
  # Exit if empty
  if (nrow(alignment_merge) == 0) {
    return(list(
      alignments = alignment_corrected_list,
      validation_summary = tibble::tibble()
    ))
  }
  
  # Ensure ref_flank exists
  if (!"ref_flank" %in% names(alignment_merge)) {
    stop("ref_flank column not found in alignment data")
  }
  
  # For each read, check if flanks are in correct left→right order
  # After strand correction with coordinates preserved:
  # - ALL reads (both original sense and corrected antisense) should have left→right
  # - Minimap2's RC-space coordinates become real-space after we RC the sequence
  # - So corrected antisense reads ALSO show left→right order
  read_order_check <- alignment_merge %>%
    group_by(file_stem, qname, ref_id) %>%
    dplyr::filter(n() == 2) %>%  # Must have both flanks
    arrange(query_start, .by_group = TRUE) %>%
    dplyr::mutate(
      flank_position = row_number(),
      first_flank = ref_flank[flank_position == 1],
      second_flank = ref_flank[flank_position == 2]
    ) %>%
    dplyr::slice(1) %>%  # One row per read for checking
    ungroup() %>%
    dplyr::mutate(
      # After proper strand correction, only left→right is valid
      correct_order = (first_flank == "left" & second_flank == "right")
    )
  
  # Summarize validation results
  validation_summary <- read_order_check %>%
    group_by(file_stem, correct_order) %>%
    dplyr::summarise(n_reads = n(), .groups = "drop") %>%
    pivot_wider(names_from = correct_order, 
                values_from = n_reads,
                values_fill = 0,
                names_prefix = "order_")
  
  # Add column names if missing
  if (!"order_TRUE" %in% names(validation_summary)) {
    validation_summary$order_TRUE <- 0
  }
  if (!"order_FALSE" %in% names(validation_summary)) {
    validation_summary$order_FALSE <- 0
  }
  
  validation_summary <- validation_summary %>%
    dplyr::rename(correct_order = order_TRUE,
                  incorrect_order = order_FALSE) %>%
    dplyr::mutate(
      total = correct_order + incorrect_order,
      retention_rate = correct_order / total
    )
  
  # Keep only reads with correct order
  reads_correct_order <- read_order_check %>%
    dplyr::filter(correct_order) %>%
    dplyr::select(file_stem, qname, ref_id)
  
  # Filter alignment data
  alignment_validated <- alignment_merge %>%
    semi_join(reads_correct_order, by = c("file_stem", "qname", "ref_id"))
  
  # Split back into list by file_stem
  alignment_validated_list <- alignment_validated %>%
    split(f = as.factor(.$file_stem))
  
  # Ensure list has same names as input (preserve order)
  alignment_validated_list <- alignment_validated_list[names(alignment_corrected_list)]
  
  return(list(
    alignments = alignment_validated_list,
    validation_summary = validation_summary
  ))
}


# ==============================================================================
# Sequence Splitting
# ==============================================================================

# Split sequences into pre/left/mid/right/post segments
split_sequences_enhanced <- function(alignment_info) {
  
  # Exit if NULL
  if (is.null(alignment_info)) {
    return(NULL)
  }
  
  # Exit if empty
  if (nrow(alignment_info) == 0) {
    return(alignment_info %>%
             dplyr::mutate(pre = NA_character_,
                           left = NA_character_,
                           mid = NA_character_,
                           right = NA_character_,
                           post = NA_character_))
  }
  
  # Ensure query_length exists
  if (!"query_length" %in% names(alignment_info)) {
    alignment_info <- alignment_info %>%
      dplyr::mutate(query_length = nchar(query_sequence))
  }
  
  # Reshape alignments wider so each read has one row
  # Group by file_stem, qname, ref_id to get both flanks on one row
  alignment_wide <- alignment_info %>%
    dplyr::select(file_stem, qname, rname, ref_id, ref_flank, strand, 
                  ref_start, ref_end, alignment_length_on_ref, 
                  query_start, query_end, alignment_length_on_query, 
                  query_sequence, query_length,
                  mapq, mismatches, insertions, deletions, cigar) %>%
    group_by(file_stem, qname, ref_id) %>%
    dplyr::arrange(query_start, query_end, .by_group = TRUE) %>%
    dplyr::mutate(alignment_number = row_number()) %>%
    ungroup() %>%
    pivot_wider(
      names_from = alignment_number,
      values_from = -all_of(c("file_stem", "qname", "ref_id", 
                              "query_sequence", "query_length")),
      names_sep = "_"
    )
  
  # Select and order columns
  alignment_wide <- alignment_wide %>%
    dplyr::select(file_stem, qname, ref_id, query_sequence, query_length,
                  matches("_1$"), matches("_2$")) %>%
    dplyr::select(-starts_with("alignment_number_"))
  
  # Exit if no sequences after reshaping
  if (nrow(alignment_wide) == 0) {
    return(alignment_wide %>%
             dplyr::mutate(pre = NA_character_,
                           left = NA_character_,
                           mid = NA_character_,
                           right = NA_character_,
                           post = NA_character_))
  }
  
  # Split sequences into segments
  alignment_split <- alignment_wide %>%
    dplyr::mutate(
      # Pre: before left flank alignment
      pre = case_when(
        query_start_1 > 1 ~ str_sub(query_sequence, 1, query_start_1 - 1),
        TRUE ~ NA_character_
      ),
      
      # Left: aligned to left flank
      left = case_when(
        !is.na(query_start_1) & !is.na(query_end_1) ~ 
          str_sub(query_sequence, query_start_1, query_end_1),
        TRUE ~ NA_character_
      ),
      
      # Mid: between flanks (contains repeat region)
      mid = case_when(
        !is.na(query_end_1) & !is.na(query_start_2) & 
          query_start_2 > query_end_1 ~ 
          str_sub(query_sequence, query_end_1 + 1, query_start_2 - 1),
        TRUE ~ NA_character_
      ),
      
      # Right: aligned to right flank
      right = case_when(
        !is.na(query_start_2) & !is.na(query_end_2) ~ 
          str_sub(query_sequence, query_start_2, query_end_2),
        TRUE ~ NA_character_
      ),
      
      # Post: after right flank alignment
      post = case_when(
        !is.na(query_end_2) & query_end_2 < query_length ~ 
          str_sub(query_sequence, query_end_2 + 1, query_length),
        TRUE ~ NA_character_
      )
    )
  
  # Handle empty strings (convert to NA)
  alignment_split <- alignment_split %>%
    dplyr::mutate(
      pre = ifelse(pre == "", NA_character_, pre),
      left = ifelse(left == "", NA_character_, left),
      mid = ifelse(mid == "", NA_character_, mid),
      right = ifelse(right == "", NA_character_, right),
      post = ifelse(post == "", NA_character_, post)
    )
  
  return(alignment_split)
}


# ==============================================================================
# Manual Read Exclusion
# ==============================================================================

# Load manual read exclusions from file
load_manual_read_exclusions <- function(exclusion_path) {
  
  # Check file exists
  if (!file.exists(exclusion_path)) {
    warning("Manual read exclusion file not found: ", exclusion_path)
    return(NULL)
  }
  
  # Determine file type and load
  file_ext <- tolower(tools::file_ext(exclusion_path))
  
  if (file_ext %in% c("txt", "tsv")) {
    # Read as tab-delimited or plain text
    exclusion_data <- tryCatch({
      read.delim(exclusion_path, header = TRUE, stringsAsFactors = FALSE)
    }, error = function(e) {
      warning("Failed to read exclusion file: ", e$message)
      return(NULL)
    })
    
  } else if (file_ext == "csv") {
    # Read as CSV
    exclusion_data <- tryCatch({
      read.csv(exclusion_path, header = TRUE, stringsAsFactors = FALSE)
    }, error = function(e) {
      warning("Failed to read exclusion file: ", e$message)
      return(NULL)
    })
    
  } else if (file_ext %in% c("xlsx", "xls")) {
    # Read Excel file
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      warning("openxlsx package required to read Excel files")
      return(NULL)
    }
    
    exclusion_data <- tryCatch({
      openxlsx::read.xlsx(exclusion_path, sheet = 1)
    }, error = function(e) {
      warning("Failed to read exclusion file: ", e$message)
      return(NULL)
    })
    
  } else {
    warning("Unsupported file format for read exclusions: ", file_ext)
    return(NULL)
  }
  
  # Check for required columns
  if (is.null(exclusion_data)) {
    return(NULL)
  }
  
  if (!all(c("file_name", "read_name") %in% names(exclusion_data))) {
    warning("Exclusion file must contain 'file_name' and 'read_name' columns")
    return(NULL)
  }
  
  # Convert to list by file_name
  exclusion_list <- exclusion_data %>%
    group_by(file_name) %>%
    dplyr::summarise(reads = list(unique(read_name)), .groups = "drop") %>%
    deframe()
  
  # Convert to named list
  exclusion_list <- setNames(exclusion_list, names(exclusion_list))
  
  return(exclusion_list)
}


# Apply manual read exclusions to alignment data
#
# Exclusions are scoped by BOTH file_stem AND qname. This means:
#   - The 'file_name' column in the exclusion file must match the pipeline's
#     file_stem exactly (i.e. filename without extension, including any _R1/_R2
#     suffix for paired-end data).
#   - A read name listed under sample A will never affect sample B, even if
#     both samples happen to contain a read with the same qname.
#
# alignment_info: combined data frame for ALL samples (has a file_stem column)
# exclusion_list: named list; names = file_stem values, elements = read name vectors
# sample_name:    the file_stem to process in this call
apply_manual_exclusions <- function(alignment_info, exclusion_list, sample_name) {
  
  # Exit if no exclusion list
  if (is.null(exclusion_list)) {
    return(alignment_info)
  }
  
  # Exit if NULL or empty alignment data
  if (is.null(alignment_info) || nrow(alignment_info) == 0) {
    return(alignment_info)
  }
  
  # Check if sample has exclusions
  if (!sample_name %in% names(exclusion_list)) {
    return(alignment_info)
  }
  
  # Get reads to exclude for this sample
  reads_to_exclude <- exclusion_list[[sample_name]]
  
  # Filter out excluded reads, scoped to this sample only.
  # The file_stem guard ensures read names are matched within the correct
  # sample and cannot accidentally remove reads from other samples.
  alignment_filtered <- alignment_info %>%
    dplyr::filter(!(file_stem == sample_name & qname %in% reads_to_exclude))
  
  # Report how many were excluded
  n_excluded <- nrow(alignment_info) - nrow(alignment_filtered)
  if (n_excluded > 0) {
    message("Excluded ", n_excluded, " reads from ", sample_name)
  }
  
  return(alignment_filtered)
}
