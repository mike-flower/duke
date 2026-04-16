# ==============================================================================
# Range Analysis Functions
# ==============================================================================
# Functions for analyzing repeat distributions within defined ranges
# ==============================================================================

#' Parse Analysis Ranges from Bracket Notation
#'
#' Parse range definitions from settings Excel (e.g., "[0-35][36-1000]")
#'
#' @param range_string Character string with range definition
#'
#' @return Named list of numeric vectors (each vector has [min, max])
#'
#' @details
#' Supports two formats:
#' - Bracketed: "[0-35][36-1000]" → list(range1 = c(0, 35), range2 = c(36, 1000))
#' - Simple: "0-35" → list(range1 = c(0, 35))
#' 
#' NA values are preserved and handled later (replaced with min/max of data)
#'
#' @export
parse_analysis_ranges <- function(range_string) {
  
  # Check for bracket notation
  has_brackets <- grepl("\\[", range_string)
  
  if (has_brackets) {
    # Extract bracketed ranges: [0-35][36-1000]
    ranges <- stringr::str_extract_all(range_string, "\\[([^\\[\\]]+)\\]")[[1]]
    ranges <- stringr::str_remove_all(ranges, "\\[|\\]")
    
    # Convert to list of numeric vectors
    range_list <- lapply(ranges, function(x) {
      parts <- as.numeric(unlist(strsplit(x, "-")))
      if (length(parts) != 2) stop("Each range must have format 'min-max'")
      parts
    })
    
    names(range_list) <- paste0("range", seq_along(range_list))
    
  } else {
    # Simple format: 0-35
    parts <- as.numeric(unlist(strsplit(range_string, "-")))
    if (length(parts) != 2) stop("Range must have format 'min-max'")
    range_list <- list(range1 = parts)
  }
  
  return(range_list)
}


#' Parse Bracketed Parameters (Generic)
#'
#' Parse any parameter that may be bracketed (floor, max_peaks, etc.)
#'
#' @param param_string Character string with parameter definition
#'
#' @return Named numeric vector of values (names are "range1", "range2", etc.)
#'
#' @details
#' Handles both formats:
#' - Bracketed: "[1][3]" → c(range1=1, range2=3)
#' - Simple: "1" → c(range1=1)
#'
#' @export
parse_bracketed_param <- function(param_string) {
  
  has_brackets <- grepl("\\[", param_string)
  
  if (has_brackets) {
    # Extract bracketed values: [1][3]
    values <- stringr::str_extract_all(param_string, "\\[([^\\[\\]]+)\\]")[[1]]
    values <- as.numeric(stringr::str_remove_all(values, "\\[|\\]"))
  } else {
    # Single value
    values <- as.numeric(param_string)
  }
  
  # Add names
  names(values) <- paste0("range", seq_along(values))
  
  return(values)
}


#' Find Modal Peaks in Repeat Distribution
#'
#' Detect local maxima in repeat length frequency distribution
#'
#' @param repeat_counts Vector of repeat counts
#' @param range_min Minimum repeat length to consider
#' @param range_max Maximum repeat length to consider  
#' @param floor Frequency threshold (absolute or proportional)
#' @param max_peaks Maximum number of peaks to return
#' @param span Window size for peak detection (must be odd)
#'
#' @return Data frame with detected peaks
#'
#' @export
find_modal_peaks <- function(repeat_counts, range_min, range_max, 
                             floor = NA, max_peaks = 2, span = 3) {
  
  # Filter to range
  repeat_counts <- repeat_counts[repeat_counts >= range_min & repeat_counts <= range_max]
  
  if (length(repeat_counts) == 0) {
    return(data.frame(
      modal_length = NA,
      modal_frequency = NA,
      peak_rank = NA
    ))
  }
  
  # Frequency table
  freq_table <- as.data.frame(table(repeat_counts))
  names(freq_table) <- c("length", "freq")
  freq_table$length <- as.numeric(as.character(freq_table$length))
  freq_table <- freq_table[order(freq_table$length), ]
  
  # Apply floor filter
  if (!is.na(floor)) {
    max_freq <- max(freq_table$freq)
    floor_val <- if (floor < 1) floor * max_freq else floor
    freq_table <- freq_table[freq_table$freq >= floor_val, ]
  }
  
  if (nrow(freq_table) == 0) {
    return(data.frame(
      modal_length = NA,
      modal_frequency = NA,
      peak_rank = NA
    ))
  }
  
  # Find local maxima using pracma::findpeaks
  if (nrow(freq_table) < span) {
    # Not enough points for peak detection - return global maximum
    max_idx <- which.max(freq_table$freq)
    return(data.frame(
      modal_length = freq_table$length[max_idx],
      modal_frequency = freq_table$freq[max_idx],
      peak_rank = 1
    ))
  }
  
  peaks <- pracma::findpeaks(freq_table$freq, 
                              minpeakheight = 0,
                              minpeakdistance = 1,
                              npeaks = max_peaks,
                              sortstr = TRUE)
  
  if (is.null(peaks)) {
    # No peaks found - return global maximum
    max_idx <- which.max(freq_table$freq)
    return(data.frame(
      modal_length = freq_table$length[max_idx],
      modal_frequency = freq_table$freq[max_idx],
      peak_rank = 1
    ))
  }
  
  # Extract peak information
  if (is.matrix(peaks)) {
    peak_heights <- peaks[, 1]
    peak_positions <- peaks[, 2]
  } else {
    peak_heights <- peaks[1]
    peak_positions <- peaks[2]
  }
  
  # Create results
  data.frame(
    modal_length = freq_table$length[peak_positions],
    modal_frequency = peak_heights,
    peak_rank = seq_along(peak_heights)
  )
}


#' Calculate Distribution Summary Statistics
#'
#' Comprehensive summary of repeat length distribution
#'
#' Returns central tendency, spread, shape, and tail metrics
#'
#' Metrics include:
#' - Central: modal_length, mean_length, median_length
#' - Spread: sd, sem, cv, iqr, percentiles
#' - Shape: skewness, kurtosis
#' - Distribution metrics: Gini, entropy, tail ratios
#'
#' @export
distribution_summary_function <- function(repeat_counts, floor = NA) {
  
  # Remove NAs
  repeat_counts <- as.numeric(na.omit(repeat_counts))
  n_reads_total <- length(repeat_counts)
  
  # Default result (all NA)
  default_result <- list(
    modal_length = NA, modal_frequency = NA,
    mean_length = NA, mean_to_modal_difference = NA,
    median_length = NA, min_length = NA, max_length = NA, range_length = NA, 
    sd = NA, sem = NA, coefficient_of_variation = NA,
    interquartile_range = NA, percentile_5th = NA, percentile_95th = NA,
    skewness = NA, kurtosis = NA, gini_coefficient = NA,
    n_reads_analysed = NA, 
    reads_below_modal = NA, reads_above_modal = NA,
    contraction_proportion = NA, expansion_proportion = NA, 
    tail_balance = NA,
    contraction_to_modal_ratio = NA, expansion_to_modal_ratio = NA,
    tail_modal_balance = NA,
    entropy_normalised = NA, peak_sharpness = NA,
    floor_threshold = floor, n_reads_total = n_reads_total, n_reads_used = NA
  )
  
  if (n_reads_total == 0) return(default_result)
  
  # Single value case
  if (n_reads_total == 1) {
    default_result$modal_length <- repeat_counts[1]
    default_result$modal_frequency <- 1
    default_result$mean_length <- repeat_counts[1]
    default_result$median_length <- repeat_counts[1]
    default_result$n_reads_analysed <- 1
    default_result$n_reads_used <- 1
    return(default_result)
  }
  
  # Create frequency table
  freq_table <- as.data.frame(table(repeat_counts))
  names(freq_table) <- c("length", "freq")
  freq_table$length <- as.numeric(as.character(freq_table$length))
  freq_table <- freq_table[order(freq_table$length), ]
  
  x <- freq_table$length
  y <- freq_table$freq
  
  # Apply floor filter
  max_freq <- max(y)
  floor_threshold <- if (!is.na(floor) && floor < 1) floor * max_freq else floor
  
  if (!is.na(floor_threshold)) {
    keep_idx <- y >= floor_threshold
    x <- x[keep_idx]
    y <- y[keep_idx]
  }
  
  n_reads_used <- sum(y)
  
  if (length(x) == 0) {
    default_result$floor_threshold <- floor_threshold
    default_result$n_reads_used <- 0
    return(default_result)
  }
  
  # Reconstruct a per-read vector from the floor-filtered frequency table.
  # n_reads_total (above) preserves the original read count for sequencing QC.
  # All descriptive statistics below are computed on this filtered population so
  # that mean, SD, SEM, and spread metrics are consistent with reads_below_modal
  # and reads_above_modal, which are also derived from the filtered data.
  filtered_counts <- rep(x, y)
  
  # Mode
  mode_idx <- which.max(y)
  mode_rpt <- x[mode_idx]
  mode_y <- y[mode_idx]
  
  # Central tendency
  mean_val <- mean(filtered_counts)
  median_val <- median(filtered_counts)
  
  # Spread
  sd_val <- sd(filtered_counts)
  sem_val <- sd_val / sqrt(n_reads_used)
  cv_val <- if (mean_val != 0) sd_val / mean_val else NA
  min_val <- min(filtered_counts)
  max_val <- max(filtered_counts)
  range_val <- max_val - min_val
  iqr_val <- IQR(filtered_counts)
  perc5 <- quantile(filtered_counts, 0.05, na.rm = TRUE)
  perc95 <- quantile(filtered_counts, 0.95, na.rm = TRUE)
  
  # Shape
  skew_val <- if (requireNamespace("moments", quietly = TRUE)) {
    moments::skewness(filtered_counts)
  } else NA
  
  kurt_val <- if (requireNamespace("moments", quietly = TRUE)) {
    moments::kurtosis(filtered_counts)
  } else NA
  
  # Gini coefficient
  gini_val <- NA
  if (requireNamespace("ineq", quietly = TRUE)) {
    gini_val <- ineq::Gini(y)
  }
  
  # Read count metrics (renamed from AUC)
  n_reads_analyzed <- sum(y)
  below_mode_idx <- x < mode_rpt
  above_mode_idx <- x > mode_rpt
  
  reads_below_modal <- if (any(below_mode_idx)) sum(y[below_mode_idx]) else 0
  reads_above_modal <- if (any(above_mode_idx)) sum(y[above_mode_idx]) else 0
  
  # Tail proportions (relative to total reads)
  contraction_proportion <- if (n_reads_analyzed > 0) reads_below_modal / n_reads_analyzed else NA
  expansion_proportion <- if (n_reads_analyzed > 0) reads_above_modal / n_reads_analyzed else NA
  tail_balance <- if (!is.na(contraction_proportion) && !is.na(expansion_proportion)) {
    (expansion_proportion - contraction_proportion) / (expansion_proportion + contraction_proportion)
  } else NA
  
  # Tail ratios (relative to modal peak)
  contraction_to_modal_ratio <- if (mode_y > 0) reads_below_modal / mode_y else NA
  expansion_to_modal_ratio <- if (mode_y > 0) reads_above_modal / mode_y else NA
  tail_modal_balance <- if (!is.na(contraction_to_modal_ratio) && !is.na(expansion_to_modal_ratio)) {
    (expansion_to_modal_ratio - contraction_to_modal_ratio) / (expansion_to_modal_ratio + contraction_to_modal_ratio)
  } else NA
  
  # Entropy
  p <- y / sum(y)
  entropy <- -sum(p * log(p + 1e-10))
  max_entropy <- log(length(y))
  entropy_normalized <- if (max_entropy > 0) entropy / max_entropy else NA
  
  # Peak sharpness
  peak_sharpness <- mode_y / mean(y)
  
  # Compile results
  list(
    modal_length = mode_rpt,
    modal_frequency = mode_y,
    mean_length = mean_val,
    mean_to_modal_difference = mean_val - mode_rpt,
    median_length = median_val,
    min_length = min_val,
    max_length = max_val,
    range_length = range_val,
    sd = sd_val,
    sem = sem_val,
    coefficient_of_variation = cv_val,
    interquartile_range = iqr_val,
    percentile_5th = as.numeric(perc5),
    percentile_95th = as.numeric(perc95),
    skewness = skew_val,
    kurtosis = kurt_val,
    gini_coefficient = gini_val,
    n_reads_analysed = n_reads_analyzed,
    reads_below_modal = reads_below_modal,
    reads_above_modal = reads_above_modal,
    contraction_proportion = contraction_proportion,
    expansion_proportion = expansion_proportion,
    tail_balance = tail_balance,
    contraction_to_modal_ratio = contraction_to_modal_ratio,
    expansion_to_modal_ratio = expansion_to_modal_ratio,
    tail_modal_balance = tail_modal_balance,
    entropy_normalised = entropy_normalized,
    peak_sharpness = peak_sharpness,
    floor_threshold = floor,
    n_reads_total = n_reads_total,
    n_reads_used = n_reads_used
  )
}


#' Calculate Instability Metrics Relative to a Setpoint
#'
#' Calculates instability from raw read distribution
#'
#' @param numbers Vector of repeat counts (raw read data)
#' @param setpoint Reference repeat length for comparison
#' @param floor Frequency threshold (absolute or proportional)
#' @return Data frame with instability metrics
#'
#' @export
calculate_instability_metrics <- function(numbers, setpoint, floor = NA) {
  
  # Validate and clean input
  if (!is.numeric(numbers)) stop("Input 'numbers' must be numeric.")
  numbers <- as.numeric(na.omit(numbers))
  n_input <- length(numbers)
  
  # Default result
  default_result <- data.frame(
    setpoint_repeat_length = setpoint,
    instability_index = NA,
    expansion_index = NA,
    contraction_index = NA,
    reads_at_setpoint = NA,
    reads_below_setpoint = NA,
    reads_above_setpoint = NA,
    contraction_ratio = NA,
    expansion_ratio = NA,
    total_tail_ratio = NA,
    expansion_vs_contraction_balance = NA,
    floor_threshold = NA,
    n_reads_total = n_input,
    n_reads_used = NA
  )
  
  # Exit early if no input
  if (n_input == 0) return(default_result)
  
  # Frequency table
  freq_table <- data.frame(table(numbers)) %>%
    dplyr::rename(number = numbers, freq = Freq) %>%
    dplyr::mutate(number = as.numeric(as.character(number))) %>%
    dplyr::arrange(number)
  
  # Filter out numbers with frequency below the floor
  max_freq <- max(freq_table$freq)
  floor_threshold <- if (!is.na(floor) && floor < 1) floor * max_freq else floor
  if (!is.na(floor_threshold)) {
    freq_table <- freq_table %>% 
      dplyr::filter(freq >= floor_threshold)
  }
  n_pass_filter <- sum(freq_table$freq)
  if (nrow(freq_table) == 0) {
    default_result$floor_threshold <- floor_threshold
    default_result$n_reads_used <- 0
    return(default_result)
  }
  
  # Filter the original number vector
  filtered_numbers <- numbers[numbers %in% freq_table$number]
  
  # Calculate distance from setpoint
  distance_from_setpoint <- filtered_numbers - setpoint
  
  # Data frame of distances
  distance_df <- data.frame(
    number = filtered_numbers,
    distance_from_setpoint = distance_from_setpoint
  ) %>%
    dplyr::count(number, distance_from_setpoint, name = "freq")
  
  # Frequency-weighted mean function
  frequency_weighted_mean <- function(data, value_column, weight_column) {
    sum(data[[value_column]] * data[[weight_column]], na.rm = TRUE) / 
      sum(data[[weight_column]], na.rm = TRUE)
  }
  
  # Calculate instability indices
  instability_index <- frequency_weighted_mean(distance_df, "distance_from_setpoint", "freq")
  expansion_index <- frequency_weighted_mean(
    dplyr::filter(distance_df, number > setpoint), 
    "distance_from_setpoint", "freq"
  )
  contraction_index <- frequency_weighted_mean(
    dplyr::filter(distance_df, number < setpoint), 
    "distance_from_setpoint", "freq"
  )
  
  # Split into three regions (using rounded setpoint for discrete counts)
  setpoint_round <- round(setpoint, 0)
  distance_setpoint <- dplyr::filter(distance_df, number == setpoint_round)
  distance_left <- dplyr::filter(distance_df, number < setpoint_round)
  distance_right <- dplyr::filter(distance_df, number > setpoint_round)
  
  # Count reads in each region
  reads_at_setpoint <- sum(distance_setpoint$freq, na.rm = TRUE)
  reads_below_setpoint <- sum(distance_left$freq, na.rm = TRUE)
  reads_above_setpoint <- sum(distance_right$freq, na.rm = TRUE)
  
  # Calculate tail ratios
  contraction_ratio <- ifelse(reads_at_setpoint > 0, 
                              reads_below_setpoint / reads_at_setpoint, NA)
  expansion_ratio <- ifelse(reads_at_setpoint > 0, 
                           reads_above_setpoint / reads_at_setpoint, NA)
  total_tail_ratio <- ifelse(reads_at_setpoint > 0, 
                             (reads_below_setpoint + reads_above_setpoint) / reads_at_setpoint,
                             NA)
  
  # Balance ratio: expansion vs contraction
  expansion_vs_contraction_balance <- ifelse(reads_below_setpoint > 0,
                                            expansion_ratio / contraction_ratio,
                                            ifelse(reads_above_setpoint > 0, Inf, NA))
  
  # Output summary
  instability_summary <- data.frame(
    setpoint_repeat_length = setpoint,
    instability_index = instability_index,
    expansion_index = expansion_index,
    contraction_index = contraction_index,
    reads_at_setpoint = reads_at_setpoint,
    reads_below_setpoint = reads_below_setpoint,
    reads_above_setpoint = reads_above_setpoint,
    contraction_ratio = contraction_ratio,
    expansion_ratio = expansion_ratio,
    total_tail_ratio = total_tail_ratio,
    expansion_vs_contraction_balance = expansion_vs_contraction_balance,
    floor_threshold = floor_threshold,
    n_reads_total = n_input,
    n_reads_used = n_pass_filter
  ) %>%
    dplyr::mutate(across(everything(), ~ ifelse(is.nan(.), NA, .)))
  
  return(instability_summary)
}


#' Calculate Group Control Setpoint
#' 
#' @param control_data Data frame with distribution summary for control samples
#' @param metric Column name to use as setpoint (e.g., "modal_length", "mean_length", "median_length")
#' @param aggregation_method How to aggregate control samples: "mean", "median", "trimmed_mean"
#' @return List with setpoint, SD, SEM, and count
#'
#' @export
calculate_control_setpoint <- function(control_data, metric, aggregation_method = "mean") {
  
  values <- control_data[[metric]]
  values <- values[!is.na(values)]
  
  if (length(values) == 0) {
    return(list(
      setpoint = NA_real_,
      sd = NA_real_,
      sem = NA_real_,
      n = 0
    ))
  }
  
  setpoint <- switch(aggregation_method,
    "mean" = mean(values, na.rm = TRUE),
    "median" = median(values, na.rm = TRUE),
    "trimmed_mean" = mean(values, trim = 0.1, na.rm = TRUE),  # 10% trimmed mean
    stop("Unknown aggregation method: ", aggregation_method)
  )
  
  # Calculate SD and SEM (useful for z-scores)
  sd_val <- sd(values, na.rm = TRUE)
  sem_val <- sd_val / sqrt(length(values))
  
  return(list(
    setpoint = setpoint,
    sd = sd_val,
    sem = sem_val,
    n = length(values)
  ))
}
