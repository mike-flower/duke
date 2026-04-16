# ==============================================================================
# Waterfall Plotting Functions
# ==============================================================================
# Functions for creating waterfall plots visualizing read structure
# ==============================================================================

#' Create Waterfall Plot
#'
#' Generate a waterfall plot showing read structure with left flank, repeat region, and right flank
#'
#' @param data Data frame with columns: left, mid, right (sequences), repeat_length, cluster_number
#' @param sample_name Character string for plot title
#' @param downsample Integer or NA. If integer, randomly sample this many reads. Default NA (no downsampling)
#' @param dna_colours Named vector of colours for bases (A, C, T, G, -, N)
#' @param repeat_count_col Character string naming the column with repeat counts for y-axis labels. Default "repeat_length"
#' @param repeat_count_method Character string naming the method used (e.g., "repeat_count_full", "repeat_count_match") for better y-axis label
#' @param y_axis_labels Character or integer. Controls y-axis labeling density:
#'   - "auto" = dynamic based on read count (default)
#'   - "all" = label every read
#'   - Integer (e.g. 10, 20) = target number of labels
#' @param total_reads Integer. Total reads before filtering (for subtitle). If NA, not shown.
#'
#' @return ggplot object
#'
#' @details
#' Creates a tile plot where:
#' - X-axis: Position (centered on repeat region, 0 = start of repeat)
#' - Y-axis: Read order (sorted by repeat length, longest at bottom)
#' - Fill: Base colour (A/C/T/G)
#' - Segments: Left flank (right-aligned) | Repeat (left-aligned) | Right flank
#'
#' @examples
#' p <- create_waterfall_plot(
#'   data = clustered_reads,
#'   sample_name = "blood",
#'   downsample = 5000,
#'   dna_colours = c("A" = "darkgreen", "C" = "blue", "T" = "red", "G" = "black"),
#'   y_axis_labels = "auto",  # or 20, or "all"
#'   total_reads = 13757
#' )
#'
#' @export
create_waterfall_plot <- function(data, 
                                   sample_name, 
                                   downsample = NA, 
                                   dna_colours = c("A" = "darkgreen", "C" = "blue", 
                                                   "T" = "red", "G" = "black", 
                                                   "-" = "lightgrey", "N" = "grey"),
                                   repeat_count_col = "repeat_length",
                                  repeat_count_method = NULL,  # NEW: for better label
                                   y_axis_labels = "auto",
                                   total_reads = NA) {
  
  # Input validation
  required_cols <- c("left", "mid", "right", repeat_count_col, "cluster_number")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Store total reads before downsampling (for subtitle)
  if (is.na(total_reads)) {
    total_reads <- nrow(data)
  }
  
  # Downsample if needed
  if (!is.na(downsample) && nrow(data) > downsample) {
    set.seed(123)  # Reproducible sampling
    data <- data %>% slice_sample(n = downsample)
    message("  Downsampled to ", downsample, " reads for ", sample_name)
  }
  
  # Sort by repeat length (descending - longest at bottom)
  data <- data %>%
    arrange(desc(!!sym(repeat_count_col)), cluster_number) %>%
    mutate(read_order = row_number())
  
  # Calculate y-axis breaks based on y_axis_labels parameter
  if (is.character(y_axis_labels)) {
    if (y_axis_labels == "all") {
      # Label every read
      y_breaks_calc <- seq(1, nrow(data), by = 1)
    } else if (y_axis_labels == "auto") {
      # Dynamic based on read count (from original Duke)
      n_reads <- nrow(data)
      if (n_reads <= 30) {
        step <- 1
      } else if (n_reads <= 50) {
        step <- 3
      } else if (n_reads <= 100) {
        step <- 5
      } else if (n_reads <= 200) {
        step <- 10
      } else if (n_reads <= 500) {
        step <- 20
      } else if (n_reads <= 1000) {
        step <- 50
      } else if (n_reads <= 5000) {
        step <- 100
      } else {
        step <- 500
      }
      y_breaks_calc <- seq(1, nrow(data), by = step)
    } else {
      # Default to pretty breaks if unrecognized string
      y_breaks_calc <- pretty_breaks(n = 10)
    }
  } else if (is.numeric(y_axis_labels)) {
    # User specified number of breaks
    y_breaks_calc <- pretty_breaks(n = y_axis_labels)
  } else {
    # Default
    y_breaks_calc <- pretty_breaks(n = 10)
  }
  
  # Convert sequences to long format for plotting
  plot_data_list <- list()
  
  for (i in 1:nrow(data)) {
    read_row <- data[i, ]
    
    # Left flank (right-aligned, negative positions)
    left_seq <- strsplit(as.character(read_row$left), "")[[1]]
    left_len <- length(left_seq)
    
    if (left_len > 0) {
      left_df <- data.frame(
        read_order = i,
        position = seq(-left_len, -1),
        base = left_seq,
        segment = "left",
        cluster = read_row$cluster_number,
        repeat_length = read_row[[repeat_count_col]],
        stringsAsFactors = FALSE
      )
      plot_data_list[[length(plot_data_list) + 1]] <- left_df
    }
    
    # Mid (repeat region, starts at position 0)
    mid_seq <- strsplit(as.character(read_row$mid), "")[[1]]
    mid_len <- length(mid_seq)
    
    if (mid_len > 0) {
      mid_df <- data.frame(
        read_order = i,
        position = seq(0, mid_len - 1),
        base = mid_seq,
        segment = "mid",
        cluster = read_row$cluster_number,
        repeat_length = read_row[[repeat_count_col]],
        stringsAsFactors = FALSE
      )
      plot_data_list[[length(plot_data_list) + 1]] <- mid_df
    }
    
    # Right flank (continues from mid)
    right_seq <- strsplit(as.character(read_row$right), "")[[1]]
    right_len <- length(right_seq)
    
    if (right_len > 0) {
      right_df <- data.frame(
        read_order = i,
        position = seq(mid_len, mid_len + right_len - 1),
        base = right_seq,
        segment = "right",
        cluster = read_row$cluster_number,
        repeat_length = read_row[[repeat_count_col]],
        stringsAsFactors = FALSE
      )
      plot_data_list[[length(plot_data_list) + 1]] <- right_df
    }
  }
  
  # Combine all data
  if (length(plot_data_list) == 0) {
    stop("No sequence data to plot for ", sample_name)
  }
  plot_data <- rbindlist(plot_data_list)
  
  # Create plot
  p <- ggplot(plot_data, aes(x = position, y = read_order, fill = base)) +
    geom_tile(width = 1, height = 1, colour = NA) +
    scale_fill_manual(values = dna_colours, name = "Base") +
    scale_y_continuous(
      breaks = y_breaks_calc,
      labels = function(x) {
        # Show repeat length instead of read number
        sapply(x, function(idx) {
          if (is.na(idx) || idx < 1 || idx > nrow(data)) {
            ""
          } else {
            round(data[[repeat_count_col]][idx])
          }
        })
      }
    ) +
    labs(
      title = paste("Waterfall plot -", sample_name),
      subtitle = paste(nrow(data), "reads plotted from", total_reads, "total"),
      x = "Position (bp, centered on repeat region)",
      y = if (!is.null(repeat_count_method)) {
        # Use the full method name as selected by user
        paste0("Repeat length (", repeat_count_method, ")")
      } else {
        # Fallback to column name
        paste0("Repeat length (", repeat_count_col, ")")
      }
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "right"
    )
  
  # Add segment annotations and dividers
  # Calculate positions for labels
  left_positions <- plot_data$position[plot_data$segment == "left"]
  mid_positions <- plot_data$position[plot_data$segment == "mid"]
  right_positions <- plot_data$position[plot_data$segment == "right"]
  
  # Add vertical line at repeat start
  p <- p + geom_vline(xintercept = 0, linetype = "dashed", colour = "black", alpha = 0.5)
  
  # Add segment labels if segments exist
  label_y <- nrow(data) * 1.05
  
  if (length(left_positions) > 0) {
    left_center <- min(left_positions) + (max(left_positions) - min(left_positions)) / 2
    p <- p + annotate("text", x = left_center, y = label_y, 
                     label = "Left flank", size = 4, fontface = "bold")
  }
  
  if (length(mid_positions) > 0) {
    mid_center <- min(mid_positions) + (max(mid_positions) - min(mid_positions)) / 2
    p <- p + annotate("text", x = mid_center, y = label_y, 
                     label = "Repeat", size = 4, fontface = "bold", colour = "red")
  }
  
  if (length(right_positions) > 0) {
    right_center <- min(right_positions) + (max(right_positions) - min(right_positions)) / 2
    p <- p + annotate("text", x = right_center, y = label_y, 
                     label = "Right flank", size = 4, fontface = "bold")
  }
  
  return(p)
}


#' Prepare Waterfall Data
#'
#' Prepare alignment data for waterfall plotting.
#' Flank length filtering is applied upstream in Module 3 if enabled —
#' this function does not apply any outlier removal.
#'
#' @param alignment_data Data frame with alignment information
#' @param repeat_count_method Character string naming the repeat count column to use
#' @param group_vars Character vector of grouping variables (unused, kept for compatibility)
#'
#' @return Data frame ready for waterfall plotting
#'
#' @export
prepare_waterfall_data <- function(alignment_data,
                                   repeat_count_method = "repeat_count_full",
                                   group_vars = c("file_stem", "cluster_number")) {
  
  data <- alignment_data %>%
    dplyr::filter(!is.na(mid), !is.na(cluster_number)) %>%
    dplyr::mutate(
      left_length   = nchar(left),
      right_length  = nchar(right),
      mid_length    = nchar(mid),
      repeat_length = !!sym(repeat_count_method)
    )
  
  message("Waterfall data prepared: ", format(nrow(data), big.mark = ","),
          " reads across ", length(unique(data$file_stem)), " samples")
  
  return(data)
}
