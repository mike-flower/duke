# ==============================================================================
# Duke Pipeline - Clustering Functions
# ==============================================================================
# Functions for clustering reads by repeat length or haplotype
# ==============================================================================

# Cluster reads by repeat length using GMM or k-means
cluster_by_repeat_length <- function(values, cluster_max = 10, cluster_downsample = 10000) {
  
  # Identify non-NA values
  non_na_index <- which(!is.na(values))
  non_na_values <- values[non_na_index]
  
  # Exit early if not enough variation — return all reads in cluster 1
  # Named by position (not by value) to match the output format of the full path
  if (length(unique(non_na_values)) < 2) {
    return(setNames(rep(1L, length(values)), seq_along(values)))
  }
  
  # Report dataset size
  n_values <- length(non_na_values)
  message("  Clustering ", n_values, " repeat length values")
  
  # Subsample for cluster optimisation if dataset is large
  use_subsample <- n_values > cluster_downsample
  if (use_subsample) {
    message("  Large dataset detected - using ", cluster_downsample, " subsampled values for cluster optimisation")
    set.seed(123)
    subsample_idx <- sample(length(non_na_values), cluster_downsample)
    optimisation_values <- non_na_values[subsample_idx]
  } else {
    optimisation_values <- non_na_values
  }
  
  # Auto-detect optimal number of clusters using GMM or silhouette
  
  # Attempt GMM (suppress verbose output with capture.output to avoid connection leaks)
  gmm <- tryCatch({
    capture.output(
      result <- Mclust(optimisation_values, G = 1:cluster_max, verbose = FALSE)
    )
    result
  },
    error = function(e) {
      message("  GMM failed: ", conditionMessage(e), ". Falling back to silhouette-based k-means.")
      return(NULL)
    }
  )
    
  if (!is.null(gmm)) {
    # GMM succeeded - use it to determine optimal k
    optimal_k <- gmm$G
    if (optimal_k == 1) {
      warning("GMM selected a single cluster (G = 1). Check for possible multimodality.")
    }
    message("  GMM selected ", optimal_k, " clusters")
    
    # Now cluster the FULL dataset with the optimal k
    if (optimal_k == 1) {
      cluster_ids <- rep(1, length(non_na_values))
    } else {
      # Use k-means on full dataset with optimal k from GMM
      set.seed(123)
      km_full <- kmeans(non_na_values, centers = optimal_k, nstart = 10)
      cluster_ids <- km_full$cluster
    }
    
  } else {
    # Fall back to silhouette-based k-means
    k_max <- min(cluster_max, length(unique(optimisation_values)) - 1)
    if (k_max < 2) {
      cluster_n <- 1
      cluster_ids <- rep(1, length(non_na_values))
    } else {
      # Compute silhouette on subsampled data
      message("  Computing silhouette widths for k=2 to k=", k_max)
      sil_widths <- sapply(2:k_max, function(k) {
        km <- kmeans(optimisation_values, centers = k, nstart = 10)
        ss <- silhouette(km$cluster, dist(optimisation_values))
        mean(ss[, "sil_width"])
      })
      optimal_k <- which.max(sil_widths) + 1  # +1 because index 1 = 2 clusters
      message("  Silhouette analysis selected ", optimal_k, " clusters")
      
      # Now cluster the FULL dataset with optimal k
      set.seed(123)
      km_full <- kmeans(non_na_values, centers = optimal_k, nstart = 10)
      cluster_ids <- km_full$cluster
    }
  }
  
  # Assign to full-length vector
  cluster_assignments <- rep(NA_integer_, length(values))
  cluster_assignments[non_na_index] <- cluster_ids
  
  # Table of values and cluster assignment
  cluster_df <- tibble(value = values, cluster = cluster_assignments)
  
  # Order clusters by mean value
  cluster_means <- cluster_df %>%
    dplyr::filter(!is.na(value)) %>%
    group_by(cluster) %>%
    dplyr::summarise(mean_value = mean(value), .groups = "drop") %>%
    dplyr::arrange(mean_value) %>%
    dplyr::mutate(cluster_new = row_number())
  
  # Reassign cluster numbers by rank order of mean
  cluster_df <- cluster_df %>%
    left_join(cluster_means, by = "cluster") %>%
    dplyr::mutate(cluster = cluster_new) %>%
    dplyr::select(-cluster_new, -mean_value)
  
  # Return named vector
  return(setNames(cluster_df$cluster, cluster_df$value))
}

# Cluster reads by haplotype (variants in flanking regions)
# Works with mixed short-read and long-read data
cluster_by_haplotype <- function(sequences, 
                                method = "levenshtein",
                                max_length_diff = NULL,
                                trim_to_length = NA,
                                cluster_downsample = 5000) {  # Downsample for speed
  
  # Remove NA sequences
  valid_idx <- which(!is.na(sequences) & sequences != "")
  valid_seqs <- sequences[valid_idx]
  
  if (length(valid_seqs) == 0) {
    return(rep(NA_integer_, length(sequences)))
  }
  
  # Check sequence lengths
  seq_lengths <- nchar(valid_seqs)
  median_length <- median(seq_lengths)
  
  # Report length distribution
  message("  Sequence length range: ", min(seq_lengths), "-", max(seq_lengths), " bp (median: ", median_length, ")")
  
  # Optional: filter outlier lengths only if max_length_diff is specified
  # Default is NULL = no filtering (handles mixed short/long read data)
  if (!is.null(max_length_diff)) {
    length_filter <- abs(seq_lengths - median_length) <= max_length_diff
    if (sum(!length_filter) > 0) {
      message("  Excluding ", sum(!length_filter), " reads with lengths >", max_length_diff, 
              "bp from median")
      message("  (median: ", median_length, " bp, excluded range: ", 
              paste(range(seq_lengths[!length_filter]), collapse = "-"), " bp)")
    }
    
    # Filter to selected sequences
    normal_idx <- valid_idx[length_filter]
    normal_seqs <- valid_seqs[length_filter]
  } else {
    # No filtering - use all sequences
    normal_idx <- valid_idx
    normal_seqs <- valid_seqs
    # Don't say "using all" yet - we might subsample later!
  }
  
  if (length(normal_seqs) < 2) {
    # Not enough sequences
    cluster_assignments <- rep(1L, length(sequences))
    return(cluster_assignments)
  }
  
  # Get sequence lengths for determining trim strategy
  seq_lengths_filtered <- nchar(normal_seqs)
  
  # Check if lengths are still variable
  still_variable <- length(unique(seq_lengths_filtered)) > 1
  
  # Force levenshtein if variable lengths and user requested hamming
  if (still_variable && method == "hamming") {
    length_range <- range(seq_lengths_filtered)
    if (diff(length_range) <= 5) {
      # Small variation (≤5bp) - probably just alignment edges
      message("  Note: Small length variation (", length_range[1], "-", length_range[2], " bp)")
      message("  Switching from Hamming to Levenshtein distance")
    } else {
      # Larger variation - mixed read types or variable flanks
      message("  Note: Variable sequence lengths detected (", length_range[1], "-", length_range[2], " bp)")
      message("  This may indicate mixed short/long read data")
      message("  Switching from Hamming to Levenshtein distance")
    }
    method <- "levenshtein"
  }
  
  # Convert to DNAStringSet
  dna_set <- DNAStringSet(normal_seqs)
  
  # Subsampling strategy for large datasets
  use_subsample <- !is.null(cluster_downsample) && !is.na(cluster_downsample) && 
                   length(normal_seqs) > cluster_downsample
  
  if (use_subsample) {
    # Sample representative reads for clustering
    set.seed(123)
    sample_idx <- sample(length(normal_seqs), cluster_downsample)
    sample_seqs <- normal_seqs[sample_idx]
    
    message("  Large dataset detected (", length(normal_seqs), " sequences)")
    message("  Subsampling ", cluster_downsample, " representatives for initial clustering...")
    
    # Calculate distance matrix on sample only
    dna_sample <- DNAStringSet(sample_seqs)
    message("  Computing ", method, " distances for ", length(sample_seqs), " sampled sequences...")
    
    if (method == "hamming") {
      dist_matrix_sample <- stringDist(dna_sample, method = "hamming")
    } else {
      dist_matrix_sample <- stringDist(dna_sample, method = "levenshtein")
    }
    
    message("  Distance calculation complete for sample")
    
    # Hierarchical clustering on sample
    hc_sample <- hclust(dist_matrix_sample, method = "complete")
    
    # Cut tree to identify clusters
    sample_clusters <- cutree(hc_sample, h = max(dist_matrix_sample) * 0.2)
    n_clusters <- length(unique(sample_clusters))
    
    message("  Identified ", n_clusters, " haplotype group(s) in sample")
    
    # Calculate centroid (consensus) for each cluster.
    # APPROXIMATION: the centroid is the first sequence in each cluster rather
    # than a true consensus. This is a speed trade-off — computing a full
    # multiple-sequence alignment for every cluster on every large sample would
    # be prohibitively slow. In practice the approximation works well when
    # within-cluster sequences are highly similar (typical for flanking regions
    # of a single haplotype). It can misassign reads that sit near a cluster
    # boundary, particularly with shorter or more variable sequences.
    message("  Computing cluster centroids (first-sequence approximation)...")
    centroids <- lapply(1:n_clusters, function(i) {
      cluster_seqs <- sample_seqs[sample_clusters == i]
      cluster_seqs[1]
    })
    
    # Assign remaining sequences to nearest centroid
    message("  Assigning ", length(normal_seqs) - cluster_downsample, " remaining sequences to clusters...")
    
    remaining_idx <- setdiff(1:length(normal_seqs), sample_idx)
    remaining_assignments <- sapply(remaining_idx, function(i) {
      seq_i <- normal_seqs[i]
      # Calculate distance to each centroid
      distances <- sapply(centroids, function(centroid) {
        if (method == "hamming" && nchar(seq_i) == nchar(centroid)) {
          sum(strsplit(seq_i, "")[[1]] != strsplit(centroid, "")[[1]])
        } else {
          # Use simple levenshtein approximation
          adist(seq_i, centroid)[1,1]
        }
      })
      which.min(distances)
    })
    
    # Combine sample and remaining assignments
    all_assignments <- integer(length(normal_seqs))
    all_assignments[sample_idx] <- sample_clusters
    all_assignments[remaining_idx] <- remaining_assignments
    
    cluster_ids <- all_assignments
    
  } else {
    # Standard approach: compute full distance matrix
    message("  Using all ", length(normal_seqs), " sequences (no subsampling)")
    message("  Computing ", method, " distances for ", length(normal_seqs), " sequences...")
    
    if (method == "hamming") {
      dist_matrix <- stringDist(dna_set, method = "hamming")
    } else {
      dist_matrix <- stringDist(dna_set, method = "levenshtein")
    }
    
    message("  Distance calculation complete")
    
    # Hierarchical clustering
    hc <- hclust(dist_matrix, method = "complete")
    
    # Cut tree to identify clusters
    cluster_ids <- cutree(hc, h = max(dist_matrix) * 0.2)
  }
  
  
  # For subsampled data, we already have cluster_ids
  # For full data, determine optimal k using silhouette
  if (!use_subsample) {
    max_k <- min(10, length(normal_seqs) - 1)
    if (max_k < 2) {
      cluster_ids <- rep(1, length(normal_seqs))
    } else {
      message("  Determining optimal number of clusters...")
      sil_widths <- sapply(2:max_k, function(k) {
        clusters <- cutree(hc, k = k)
        ss <- silhouette(clusters, dist_matrix)
        mean(ss[, "sil_width"])
      })
      optimal_k <- which.max(sil_widths) + 1
      cluster_ids <- cutree(hc, k = optimal_k)
      
      message("  Detected ", optimal_k, " haplotype group(s) based on ", 
              ifelse(method == "hamming", "SNPs", "SNPs and indels"))
    }
  }
  
  # cluster_ids now contains assignments for all normal_seqs
  n_clusters <- length(unique(cluster_ids))
  message("  Final result: ", n_clusters, " haplotype cluster(s)")
  
  # Assign to full vector
  cluster_assignments <- rep(NA_integer_, length(sequences))
  cluster_assignments[normal_idx] <- cluster_ids
  
  return(cluster_assignments)
}

# Cluster by combination with flexible ordering
# cluster_order can be c("haplotype", "repeat") or c("repeat", "haplotype")
cluster_by_combination <- function(sequences, repeat_values, 
                                   cluster_order = c("repeat", "haplotype"),
                                   haplotype_method = "levenshtein",
                                   haplotype_cluster_max = 10,
                                   repeat_cluster_max = 20,
                                   haplotype_cluster_downsample = 5000,
                                   repeat_cluster_downsample = 10000) {
  
  # Validate cluster_order
  valid_orders <- list(
    c("repeat", "haplotype"),
    c("haplotype", "repeat")
  )
  
  if (!any(sapply(valid_orders, function(x) identical(sort(cluster_order), sort(x))))) {
    stop("cluster_order must be c('repeat', 'haplotype') or c('haplotype', 'repeat')")
  }
  
  # Determine primary and secondary clustering
  primary <- cluster_order[1]
  secondary <- cluster_order[2]
  
  message("Clustering order: ", primary, " → ", secondary)
  
  # Step 1: Primary clustering
  if (primary == "repeat") {
    # Cluster by repeat length first
    primary_clusters <- cluster_by_repeat_length(
      values = repeat_values,
      cluster_max = repeat_cluster_max,
      cluster_downsample = repeat_cluster_downsample
    )
    primary_data <- repeat_values
    secondary_data <- sequences
    secondary_max <- haplotype_cluster_max
    message("  Primary: Found ", length(unique(primary_clusters[!is.na(primary_clusters)])), 
            " repeat cluster(s)")
    
  } else {  # primary == "haplotype"
    # Cluster by haplotype first
    primary_clusters <- cluster_by_haplotype(
      sequences = sequences,
      method = haplotype_method,
      max_length_diff = NULL,
      trim_to_length = NA,
      cluster_downsample = haplotype_cluster_downsample
    )
    primary_data <- sequences
    secondary_data <- repeat_values
    secondary_max <- repeat_cluster_max
    message("  Primary: Found ", length(unique(primary_clusters[!is.na(primary_clusters)])), 
            " haplotype cluster(s)")
  }
  
  # Step 2: Within each primary cluster, do secondary clustering
  unique_primary <- unique(primary_clusters[!is.na(primary_clusters)])
  
  final_clusters <- rep(NA_integer_, length(sequences))
  repeat_cluster_labels <- rep(NA_character_, length(sequences))
  haplotype_cluster_labels <- rep(NA_character_, length(sequences))
  cluster_counter <- 1
  
  for (p in unique_primary) {
    # Get reads in this primary cluster
    primary_idx <- which(primary_clusters == p)
    
    # Secondary clustering
    if (secondary == "repeat") {
      # Cluster by repeat within this primary cluster
      secondary_clusters <- cluster_by_repeat_length(
        values = repeat_values[primary_idx],
        cluster_max = secondary_max,
        cluster_downsample = repeat_cluster_downsample
      )
      
    } else {  # secondary == "haplotype"
      # Cluster by haplotype within this primary cluster
      secondary_clusters <- cluster_by_haplotype(
        sequences = sequences[primary_idx],
        method = haplotype_method,
        max_length_diff = NULL,
        trim_to_length = NA,
        cluster_downsample = haplotype_cluster_downsample
      )
    }
    
    # Assign unique cluster numbers and track R/H labels
    unique_secondary <- unique(secondary_clusters[!is.na(secondary_clusters)])
    for (s in unique_secondary) {
      final_idx <- primary_idx[which(secondary_clusters == s)]
      final_clusters[final_idx] <- cluster_counter
      
      # Assign R and H labels based on clustering order
      if (primary == "repeat") {
        repeat_cluster_labels[final_idx] <- paste0("R", p)
        haplotype_cluster_labels[final_idx] <- paste0("H", s)
      } else {
        haplotype_cluster_labels[final_idx] <- paste0("H", p)
        repeat_cluster_labels[final_idx] <- paste0("R", s)
      }
      
      cluster_counter <- cluster_counter + 1
    }
  }
  
  message("  Final: ", cluster_counter - 1, " total cluster(s)")
  
  # Return list with cluster assignments and metadata
  return(list(
    cluster_number = final_clusters,
    repeat_cluster = repeat_cluster_labels,
    haplotype_cluster = haplotype_cluster_labels,
    n_repeat_clusters = length(unique(repeat_cluster_labels[!is.na(repeat_cluster_labels)])),
    n_haplotype_clusters = length(unique(haplotype_cluster_labels[!is.na(haplotype_cluster_labels)]))
  ))
}

# Wrapper function to apply clustering based on user choice
apply_clustering <- function(alignment_data, 
                             cluster_by = "repeat",  # Can be vector or single string
                             repeat_column = "repeat_count",
                             sequence_column = "sequence",
                             left_column = "left",
                             right_column = "right",
                             haplotype_cluster_max = 10,
                             repeat_cluster_max = 20,
                             repeat_cluster_downsample = 10000,
                             haplotype_cluster_downsample = 5000,
                             haplotype_method = "levenshtein",
                             haplotype_region = "both",
                             haplotype_trim_length = NA) {
  
  # Parse cluster_by parameter (can be vector or single value)
  if (length(cluster_by) == 0 || (length(cluster_by) == 1 && cluster_by == "none")) {
    # No clustering
    message("No clustering - all reads in single group")
    cluster_assignments <- rep(1L, nrow(alignment_data))
    return(cluster_assignments)
    
  } else if (length(cluster_by) == 1) {
    # Single method clustering
    
    if (cluster_by == "repeat" || cluster_by == "repeat_length") {
      # Accept both "repeat" and "repeat_length" for backward compatibility
      message("Strategy: repeat_length")
      cluster_assignments <- cluster_by_repeat_length(
        values = alignment_data[[repeat_column]],
        cluster_max = repeat_cluster_max,
        cluster_downsample = repeat_cluster_downsample
      )
      
    } else if (cluster_by == "haplotype") {
      message("Strategy: haplotype")
    # Cluster by haplotype (flanking sequences)
    
    # Choose which region to use for haplotype
    # If trimming, keep bases ADJACENT to repeat (most informative, most reliable)
    
    # Determine trim length (can be numeric, "auto", or NA)
    if (!is.null(haplotype_trim_length) && !is.na(haplotype_trim_length)) {
      
      if (haplotype_trim_length == "auto") {
        # Automatic trimming: use modal (most common) length
        
        # Get lengths of left and right flanks
        left_lengths <- nchar(alignment_data[[left_column]][!is.na(alignment_data[[left_column]])])
        right_lengths <- nchar(alignment_data[[right_column]][!is.na(alignment_data[[right_column]])])
        
        # Find modal length for each flank
        get_modal_length <- function(lengths) {
          length_table <- table(lengths)
          as.numeric(names(length_table)[which.max(length_table)])
        }
        
        trim_left <- get_modal_length(left_lengths)
        trim_right <- get_modal_length(right_lengths)
        
        message("  Auto-trim: Using modal lengths (left: ", trim_left, " bp, right: ", trim_right, " bp)")
        
        # Store as list for asymmetric trimming
        trim_per_flank <- list(left = trim_left, right = trim_right)
        
      } else {
        # Manual trim length specified
        trim_per_flank <- ceiling(as.numeric(haplotype_trim_length) / 2)  # Split between flanks
      }
      
    } else {
      trim_per_flank <- NA
    }
    
    if (haplotype_region == "left") {
      haplotype_seqs <- alignment_data[[left_column]]
      # Trim left flank: keep LAST N bases (closest to repeat)
      if (!is.na(trim_per_flank)[1]) {
        trim_length <- if (is.list(trim_per_flank)) trim_per_flank$left else trim_per_flank
        haplotype_seqs <- sapply(haplotype_seqs, function(seq) {
          if (!is.na(seq) && nchar(seq) > trim_length) {
            substr(seq, nchar(seq) - trim_length + 1, nchar(seq))
          } else {
            seq
          }
        })
      }
      
    } else if (haplotype_region == "right") {
      haplotype_seqs <- alignment_data[[right_column]]
      # Trim right flank: keep FIRST N bases (closest to repeat)
      if (!is.na(trim_per_flank)[1]) {
        trim_length <- if (is.list(trim_per_flank)) trim_per_flank$right else trim_per_flank
        haplotype_seqs <- sapply(haplotype_seqs, function(seq) {
          if (!is.na(seq) && nchar(seq) > trim_length) {
            substr(seq, 1, trim_length)
          } else {
            seq
          }
        })
      }
      
    } else if (haplotype_region == "both") {
      # Trim BEFORE concatenating
      # Left: keep LAST N bases, Right: keep FIRST N bases
      left_seqs <- alignment_data[[left_column]]
      right_seqs <- alignment_data[[right_column]]
      
      if (!is.na(trim_per_flank)[1]) {
        # Get trim lengths (may be different for left and right in auto mode)
        trim_left <- if (is.list(trim_per_flank)) trim_per_flank$left else trim_per_flank
        trim_right <- if (is.list(trim_per_flank)) trim_per_flank$right else trim_per_flank
        
        # Trim left flank: keep END (closest to repeat)
        left_seqs <- sapply(left_seqs, function(seq) {
          if (!is.na(seq) && nchar(seq) > trim_left) {
            substr(seq, nchar(seq) - trim_left + 1, nchar(seq))
          } else {
            seq
          }
        })
        # Trim right flank: keep START (closest to repeat)
        right_seqs <- sapply(right_seqs, function(seq) {
          if (!is.na(seq) && nchar(seq) > trim_right) {
            substr(seq, 1, trim_right)
          } else {
            seq
          }
        })
        
        if (is.list(trim_per_flank)) {
          message("  Trimmed flanks (left: ", trim_left, " bp, right: ", trim_right, " bp)")
        } else {
          message("  Trimmed flanks to ", trim_per_flank, " bp each (keeping bases adjacent to repeat)")
        }
      }
      
      # Now concatenate trimmed flanks
      haplotype_seqs <- paste0(left_seqs, right_seqs)
      
    } else {
      stop("haplotype_region must be 'left', 'right', or 'both'")
    }
    
    cluster_assignments <- cluster_by_haplotype(
      sequences = haplotype_seqs,
      method = haplotype_method,
      max_length_diff = NULL,
      trim_to_length = haplotype_trim_length,
      cluster_downsample = haplotype_cluster_downsample
    )
    return(cluster_assignments)
      
    } else {
      stop("cluster_by must be 'repeat', 'haplotype', c('repeat', 'haplotype'), or 'none'")
    }
    
  } else if (length(cluster_by) == 2) {
    # Two-step clustering - order specified by vector
    
    # Validate input
    valid_methods <- all(cluster_by %in% c("repeat", "haplotype", "repeat_length"))
    unique_methods <- length(unique(gsub("repeat_length", "repeat", cluster_by))) == 2
    
    if (!valid_methods || !unique_methods) {
      stop("cluster_by must contain exactly 'repeat' and 'haplotype' in desired order")
    }
    
    # Normalize "repeat_length" to "repeat"
    cluster_order <- gsub("repeat_length", "repeat", cluster_by)
    
    message("Strategy: ", paste(cluster_order, collapse = " → "))
    message("  Haplotype region: ", haplotype_region)
    message("  Haplotype method: ", haplotype_method)
    # First haplotype, then repeat length within haplotype
    
    # Choose haplotype region and trim if requested
    if (!is.null(haplotype_trim_length) && !is.na(haplotype_trim_length)) {
      
      if (haplotype_trim_length == "auto") {
        # Automatic trimming: use modal (most common) length
        left_lengths <- nchar(alignment_data[[left_column]][!is.na(alignment_data[[left_column]])])
        right_lengths <- nchar(alignment_data[[right_column]][!is.na(alignment_data[[right_column]])])
        
        get_modal_length <- function(lengths) {
          length_table <- table(lengths)
          as.numeric(names(length_table)[which.max(length_table)])
        }
        
        trim_left <- get_modal_length(left_lengths)
        trim_right <- get_modal_length(right_lengths)
        
        message("  Auto-trim: Using modal lengths (left: ", trim_left, " bp, right: ", trim_right, " bp)")
        trim_per_flank <- list(left = trim_left, right = trim_right)
        
      } else {
        trim_per_flank <- ceiling(as.numeric(haplotype_trim_length) / 2)
      }
      
    } else {
      trim_per_flank <- NA
    }
    
    if (haplotype_region == "left") {
      haplotype_seqs <- alignment_data[[left_column]]
      if (!is.na(trim_per_flank)[1]) {
        trim_length <- if (is.list(trim_per_flank)) trim_per_flank$left else trim_per_flank
        haplotype_seqs <- sapply(haplotype_seqs, function(seq) {
          if (!is.na(seq) && nchar(seq) > trim_length) {
            substr(seq, nchar(seq) - trim_length + 1, nchar(seq))
          } else {
            seq
          }
        })
      }
    } else if (haplotype_region == "right") {
      haplotype_seqs <- alignment_data[[right_column]]
      if (!is.na(trim_per_flank)[1]) {
        trim_length <- if (is.list(trim_per_flank)) trim_per_flank$right else trim_per_flank
        haplotype_seqs <- sapply(haplotype_seqs, function(seq) {
          if (!is.na(seq) && nchar(seq) > trim_length) {
            substr(seq, 1, trim_length)
          } else {
            seq
          }
        })
      }
    } else {
      left_seqs <- alignment_data[[left_column]]
      right_seqs <- alignment_data[[right_column]]
      
      if (!is.na(trim_per_flank)[1]) {
        trim_left <- if (is.list(trim_per_flank)) trim_per_flank$left else trim_per_flank
        trim_right <- if (is.list(trim_per_flank)) trim_per_flank$right else trim_per_flank
        
        left_seqs <- sapply(left_seqs, function(seq) {
          if (!is.na(seq) && nchar(seq) > trim_left) {
            substr(seq, nchar(seq) - trim_left + 1, nchar(seq))
          } else {
            seq
          }
        })
        right_seqs <- sapply(right_seqs, function(seq) {
          if (!is.na(seq) && nchar(seq) > trim_right) {
            substr(seq, 1, trim_right)
          } else {
            seq
          }
        })
        
        if (is.list(trim_per_flank)) {
          message("  Trimmed flanks (left: ", trim_left, " bp, right: ", trim_right, " bp)")
        } else {
          message("  Trimmed flanks to ", trim_per_flank, " bp each (keeping bases adjacent to repeat)")
        }
      }
      
      haplotype_seqs <- paste0(left_seqs, right_seqs)
    }
    
    # Use flexible ordering with cluster_by_combination
    cluster_result <- cluster_by_combination(
      sequences = haplotype_seqs,
      repeat_values = alignment_data[[repeat_column]],
      cluster_order = cluster_order,
      haplotype_method = haplotype_method,
      haplotype_cluster_max = haplotype_cluster_max,
      repeat_cluster_max = repeat_cluster_max,
      haplotype_cluster_downsample = haplotype_cluster_downsample,
      repeat_cluster_downsample = repeat_cluster_downsample
    )
    
    # cluster_by_combination returns a list with metadata
    return(cluster_result)
    
  } else {
    stop("cluster_by must be length 1 ('repeat'/'haplotype'/'none') or length 2 (c('repeat', 'haplotype'))")
  }
}
