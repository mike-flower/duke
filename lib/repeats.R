# ==============================================================================
# Duke Pipeline - Repeat Detection Functions (Vectorized)
# ==============================================================================
# Optimized functions for finding and counting repeat tracts
# ==============================================================================

# Find repeat tracts in a single sequence
find_repeat_tracts <- function(sequence,
                               pattern,
                               min_repeats = 3,
                               max_mismatch = 0,
                               start_perfect_repeats = 1,
                               end_perfect_repeats = 1,
                               max_gap = 3,
                               max_tract_gap = NULL,
                               return_option = "all") {
  
  # Exit if there is no sequence
  if (is.na(sequence) || sequence == "") {
    return(list(start = NA, end = NA, repeat_count = NA))
  }
  
  # Find repeat units throughout the sequence
  matches <- matchPattern(pattern = pattern, 
                          subject = DNAString(sequence), 
                          max.mismatch = max_mismatch)
  
  # Exit if no matches found
  if (length(matches) == 0) {
    return(list(start = NA, end = NA, repeat_count = NA))
  }
  
  # Check if only one match is found
  if (length(matches) == 1) {
    if (min_repeats > 1) {
      return(list(start = NA, end = NA, repeat_count = NA))
    } else {
      return(list(start = start(matches)[1], end = end(matches)[1], repeat_count = 1))
    }
  }
  
  # Initiate lists to store repeat tract information
  tracts <- list()
  current_tract <- list(start = start(matches)[1], 
                        end = end(matches)[1], 
                        repeat_count = 1)
  
  # Analyse the matches to form repeat tracts
  for (i in 2:length(matches)) {
    
    # Measure gap between this match and the next
    gap_between_matches <- start(matches)[i] - end(matches)[i - 1] - 1
    
    # Use the gap to determine whether to end the current repeat tract
    if (gap_between_matches <= max_gap) {
      
      # Update repeat tract end
      current_tract$end <- end(matches)[i]
      
      # Increment repeat count
      current_tract$repeat_count <- current_tract$repeat_count + 1
      
    } else {
      
      # Store the current tract if it is long enough
      if (current_tract$repeat_count >= min_repeats) {
        tracts[[length(tracts) + 1]] <- current_tract
      }
      
      # Initiate a new tract starting at the current match
      current_tract <- list(start = start(matches)[i], 
                            end = end(matches)[i], 
                            repeat_count = 1)
    }
  }
  
  # Store the current tract if it meets the minimum repeat requirement
  if (current_tract$repeat_count >= min_repeats) {
    tracts[[length(tracts) + 1]] <- current_tract
  }
  
  # Optionally combine tracts that are within the tract_gap distance
  if (!is.null(max_tract_gap) && length(tracts) > 1) {
    
    combined_tracts <- list()
    current_tract <- tracts[[1]]
    
    for (i in 2:length(tracts)) {
      
      gap_between_tracts <- tracts[[i]]$start - current_tract$end - 1
      
      if (gap_between_tracts <= max_tract_gap) {
        current_tract$end <- tracts[[i]]$end
        current_tract$repeat_count <- current_tract$repeat_count + tracts[[i]]$repeat_count
      } else {
        combined_tracts[[length(combined_tracts) + 1]] <- current_tract
        current_tract <- tracts[[i]]
      }
    }
    
    combined_tracts[[length(combined_tracts) + 1]] <- current_tract
    tracts <- combined_tracts
  }
  
  # Trim the start and end points based on perfect repeat constraints
  trimmed_tracts <- lapply(tracts, function(tract) {
    
    # Trim forward from the start for x perfect repeats
    if (start_perfect_repeats > 0) {
      for (pos in seq(tract$start, tract$end - (start_perfect_repeats * nchar(pattern)) + 1)) {
        if (identical(substring(sequence, pos, pos + (start_perfect_repeats * nchar(pattern)) - 1), 
                      paste(rep(pattern, start_perfect_repeats), collapse = ""))) {
          tract$start <- pos
          break
        }
      }
    }
    
    # Trim backward from the end for y perfect repeats
    if (end_perfect_repeats > 0) {
      for (pos in seq(tract$end, tract$start + (end_perfect_repeats * nchar(pattern)) - 1, by = -1)) {
        if (identical(substring(sequence, pos - (end_perfect_repeats * nchar(pattern)) + 1, pos), 
                      paste(rep(pattern, end_perfect_repeats), collapse = ""))) {
          tract$end <- pos
          break
        }
      }
    }
    
    return(tract)
  })
  
  # Exit if no valid tracts
  if (length(trimmed_tracts) == 0) {
    return(list(start = NA, end = NA, repeat_count = NA))
  }
  
  # Select which tract(s) to return
  if (return_option == "first") {
    return(trimmed_tracts[[1]])
  } else if (return_option == "longest") {
    longest <- which.max(sapply(trimmed_tracts, function(x) x$end - x$start + 1))
    return(trimmed_tracts[[longest]])
  } else {
    return(trimmed_tracts)
  }
}

# ==============================================================================
# VECTORIZED: Count repeats for a vector of sequences
# ==============================================================================
# Calculates three different repeat count methods:
#   1. repeat_count_full: tract_length ÷ pattern_length (robust, includes errors)
#   2. repeat_count_match: count of exact matches (stringent, excludes errors)
#   3. repeat_count_tracts: tract structure like "3,2,2" (for QC visualization)
# ==============================================================================

count_repeats_vectorized <- function(sequences, pattern, round_digits = NA) {
  
  # Handle NA or empty sequences
  is_valid <- !is.na(sequences) & sequences != ""
  
  # Initialize output vectors
  n <- length(sequences)
  repeat_count_full <- rep(NA_real_, n)
  repeat_count_match <- rep(NA_integer_, n)
  repeat_count_tracts <- rep(NA_character_, n)
  
  if (any(is_valid)) {
    valid_seqs <- sequences[is_valid]
    
    # Vectorized: Calculate full repeat count
    repeat_count_full[is_valid] <- nchar(valid_seqs) / nchar(pattern)
    if (!is.na(round_digits)) {
      repeat_count_full[is_valid] <- round(repeat_count_full[is_valid], digits = round_digits)
    }
    
    # Vectorized: Count exact matches
    repeat_count_match[is_valid] <- str_count(valid_seqs, pattern)
    
    # For tract counting, we need to iterate (can't fully vectorize this part)
    # This calculates tract structure (e.g., "3,2,2" for three tracts)
    # ANY gap (even 1bp) breaks tract continuity
    for (i in which(is_valid)) {
      matches <- str_locate_all(sequences[i], pattern)[[1]]
      
      if (nrow(matches) > 0) {
        # Calculate gaps between consecutive matches
        gaps <- matches[-1, 1] - matches[-nrow(matches), 2] - 1
        
        # Any gap > 0 breaks the tract (strict continuity)
        stretch_starts <- c(1, which(gaps > 0) + 1)
        stretch_ends <- c(which(gaps > 0), nrow(matches))
        repeat_count_tracts[i] <- paste(stretch_ends - stretch_starts + 1, collapse = ",")
      } else {
        repeat_count_tracts[i] <- "0"
      }
    }
  }
  
  # Return as data frame for easy merging
  return(data.frame(
    repeat_count_full = repeat_count_full,
    repeat_count_match = repeat_count_match,
    repeat_count_tracts = repeat_count_tracts,
    stringsAsFactors = FALSE
  ))
}

# Apply find_repeat_tracts to a data frame - IMPROVED VERSION
apply_find_repeat_tracts <- function(alignment_info,
                                     sequence_column,
                                     pattern,
                                     min_repeats = 3,
                                     max_mismatch = 0,
                                     start_perfect_repeats = 1,
                                     end_perfect_repeats = 1,
                                     max_gap = 3,
                                     max_tract_gap = NULL,
                                     return_option = "all") {
  
  if (is.null(alignment_info) || nrow(alignment_info) == 0) return(NULL)
  
  # Store original read names for verification
  original_qnames <- alignment_info$qname
  
  # Use pbapply for progress bar without rowwise()
  results <- pblapply(seq_len(nrow(alignment_info)), function(i) {
    seq_value <- alignment_info[[sequence_column]][i]
    split_result <- find_repeat_tracts(
      sequence = seq_value,
      pattern = pattern,
      min_repeats = min_repeats,
      max_mismatch = max_mismatch,
      start_perfect_repeats = start_perfect_repeats,
      end_perfect_repeats = end_perfect_repeats,
      max_gap = max_gap,
      max_tract_gap = max_tract_gap,
      return_option = return_option
    )
    
    # Extract sequence segments
    if (!is.na(split_result$start)) {
      seq1 <- if (split_result$start > 1) {
        substr(seq_value, 1, split_result$start - 1)
      } else { NA_character_ }
      
      seq2 <- substr(seq_value, split_result$start, split_result$end)
      
      seq3 <- if (split_result$end < nchar(seq_value)) {
        substr(seq_value, split_result$end + 1, nchar(seq_value))
      } else { NA_character_ }
    } else {
      seq1 <- seq_value
      seq2 <- NA_character_
      seq3 <- NA_character_
    }
    
    return(list(
      pattern_start = split_result$start,
      pattern_end = split_result$end,
      pattern_count = split_result$repeat_count,
      seq1 = seq1,
      seq2 = seq2,
      seq3 = seq3
    ))
  })
  
  # Convert results to data frame
  results_df <- data.frame(
    pattern_start = sapply(results, `[[`, "pattern_start"),
    pattern_end = sapply(results, `[[`, "pattern_end"),
    pattern_count = sapply(results, `[[`, "pattern_count"),
    stringsAsFactors = FALSE
  )
  
  # Add sequence columns
  seq1_col <- paste0(sequence_column, ".1")
  seq2_col <- paste0(sequence_column, ".2")
  seq3_col <- paste0(sequence_column, ".3")
  
  results_df[[seq1_col]] <- sapply(results, `[[`, "seq1")
  results_df[[seq2_col]] <- sapply(results, `[[`, "seq2")
  results_df[[seq3_col]] <- sapply(results, `[[`, "seq3")
  
  # Combine with original data
  alignment_info <- cbind(alignment_info, results_df)
  
  # VERIFICATION: Check read names still match
  if (!all(alignment_info$qname == original_qnames)) {
    stop("CRITICAL ERROR: Read order has been scrambled!")
  }
  
  # Rename pattern columns
  names(alignment_info)[names(alignment_info) == "pattern_start"] <- 
    paste0(sequence_column, "_pattern_start")
  names(alignment_info)[names(alignment_info) == "pattern_end"] <- 
    paste0(sequence_column, "_pattern_end")
  names(alignment_info)[names(alignment_info) == "pattern_count"] <- 
    paste0(sequence_column, "_pattern_count")
  
  # Remove original sequence column
  alignment_info[[sequence_column]] <- NULL
  
  return(alignment_info)
}

# ==============================================================================
# Apply count_repeats to alignment data - VECTORIZED VERSION
# ==============================================================================
# Calculates all three repeat counting methods and adds them as columns:
#   - repeat_count_full: length-based (most robust for ONT data)
#   - repeat_count_match: exact match counting (stringent)
#   - repeat_count_tracts: tract structure string (for QC only)
# 
# The primary method is selected in duke_run.R via repeat_count_method parameter
# ==============================================================================

apply_count_repeats <- function(alignment_info, string_column, pattern, round_digits = 0) {
  
  if (is.null(alignment_info) || nrow(alignment_info) == 0) return(NULL)
  
  # Vectorized counting (much faster!)
  count_results <- count_repeats_vectorized(
    sequences = alignment_info[[string_column]],
    pattern = pattern,
    round_digits = round_digits
  )
  
  # Add results as new columns
  alignment_info$repeat_count_full <- count_results$repeat_count_full
  alignment_info$repeat_count_match <- count_results$repeat_count_match
  alignment_info$repeat_count_tracts <- count_results$repeat_count_tracts
  
  # Relocate new columns after the sequence column
  col_idx <- which(names(alignment_info) == string_column)
  alignment_info <- alignment_info[, c(
    seq_len(col_idx),
    which(names(alignment_info) %in% c("repeat_count_full", "repeat_count_match", "repeat_count_tracts")),
    setdiff(seq_len(ncol(alignment_info)), 
            c(seq_len(col_idx), 
              which(names(alignment_info) %in% c("repeat_count_full", "repeat_count_match", "repeat_count_tracts"))))
  )]
  
  return(alignment_info)
}

# Helper function: replace NA with blank
na_to_blank <- function(string) {
  ifelse(is.na(string), "", string)
}

# ==============================================================================
# Per-sample wrapper for repeat tract finding (for use with apply_fn)
# ==============================================================================
find_repeat_tracts_sample <- function(file_stem,
                                      alignment_split,
                                      sequence_column = "mid",
                                      pattern,
                                      min_repeats = 3,
                                      max_mismatch = 0,
                                      start_perfect_repeats = 1,
                                      end_perfect_repeats = 1,
                                      max_gap = 3,
                                      max_tract_gap = NULL,
                                      return_option = "all",
                                      out_dir = NULL,
                                      resume = FALSE) {
  
  # Check for cached result
  if (!is.null(out_dir) && resume) {
    cache_file <- file.path(out_dir, paste0("repeat_tracts_", file_stem, ".RData"))
    if (file.exists(cache_file)) {
      load(cache_file)  # Loads 'sample_repeats'
      return(sample_repeats)
    }
  }
  
  # Get data for this sample
  sample_data <- alignment_split[alignment_split$file_stem == file_stem, ]
  
  if (nrow(sample_data) == 0) {
    warning("No data for sample: ", file_stem)
    return(NULL)
  }
  
  # Process this sample
  sample_repeats <- apply_find_repeat_tracts(
    alignment_info = sample_data,
    sequence_column = sequence_column,
    pattern = pattern,
    min_repeats = min_repeats,
    max_mismatch = max_mismatch,
    start_perfect_repeats = start_perfect_repeats,
    end_perfect_repeats = end_perfect_repeats,
    max_gap = max_gap,
    max_tract_gap = max_tract_gap,
    return_option = return_option
  )
  
  # Cache the result
  if (!is.null(out_dir)) {
    cache_file <- file.path(out_dir, paste0("repeat_tracts_", file_stem, ".RData"))
    save(sample_repeats, file = cache_file)
  }
  
  return(sample_repeats)
}
