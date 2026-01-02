# Repeat Distribution Visualisation Functions
# For Module 7: Repeat Distribution Visualisation

#' Plot Repeat Frequency Histogram
#'
#' @param data Data frame with repeat counts
#' @param repeat_col Column name containing repeat lengths
#' @param binwidth Histogram bin width
#' @param manual_points Data frame with points to overlay (x, y, name)
#' @param manual_boxes Data frame with boxes to overlay (xmin, xmax, name)
#' @param range_colors Named vector of colors for ranges
#' @param title Plot title
#' @param base_text_size Base text size
#'
#' @return ggplot object
plot_frequency_histogram <- function(data,
                                     repeat_col,
                                     binwidth = 1,
                                     manual_points = NULL,
                                     manual_boxes = NULL,
                                     range_colors = NULL,
                                     title = NULL,
                                     base_text_size = 12) {
  
  # Rename for point labels
  point_rename <- c("mode" = "Modal", "mean" = "Mean")
  
  # Base plot
  p <- ggplot(data, aes(x = .data[[repeat_col]])) +
    geom_histogram(binwidth = binwidth,
                   fill = "lightgrey",
                   colour = "black",
                   alpha = 0.7)
  
  # Add range boxes if provided
  if (!is.null(manual_boxes)) {
    p <- p +
      geom_rect(data = manual_boxes,
                aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf,
                    fill = name),
                inherit.aes = FALSE,
                colour = NA,
                alpha = 0.1)
    
    # Apply consistent colors if provided
    if (!is.null(range_colors)) {
      p <- p + scale_fill_manual(values = range_colors, name = "Analysis Range")
    }
  }
  
  # Add points if provided
  if (!is.null(manual_points)) {
    # Define metric shapes (distinct and clear)
    metric_shapes <- c("modal" = 16,      # Circle (filled)
                       "mean" = 17,       # Triangle (filled)
                       "median" = 15)     # Square (filled)
    metric_labels <- c("modal" = "Modal", "mean" = "Mean", "median" = "Median")
    
    # Filter to only metrics present in data
    present_metrics <- unique(manual_points$name)
    metric_shapes_filtered <- metric_shapes[names(metric_shapes) %in% present_metrics]
    metric_labels_filtered <- metric_labels[names(metric_labels) %in% present_metrics]
    
    # If manual_points has 'range' column and range_colors provided, color by range
    # Otherwise, use black for all points
    if ("range" %in% names(manual_points) && !is.null(range_colors)) {
      p <- p +
        geom_point(data = manual_points,
                   aes(x = x, y = y, 
                       shape = name,
                       colour = range),
                   size = 4) +
        geom_text_repel(data = manual_points,
                        aes(x = x, y = y, 
                            label = paste0(metric_labels[name], ": ", round(x, 1)),
                            colour = range),
                        size = base_text_size / 3.5,
                        show.legend = FALSE) +
        scale_shape_manual(values = metric_shapes_filtered,
                          labels = metric_labels_filtered,
                          name = "Metric") +
        scale_colour_manual(values = range_colors,
                           name = "Analysis Range")
    } else {
      # No range info - use shape for metric, black color
      p <- p +
        geom_point(data = manual_points,
                   aes(x = x, y = y, shape = name),
                   size = 4,
                   colour = "black") +
        geom_text_repel(data = manual_points,
                        aes(x = x, y = y, 
                            label = paste0(metric_labels[name], ": ", round(x, 1))),
                        size = base_text_size / 3.5,
                        colour = "black",
                        show.legend = FALSE) +
        scale_shape_manual(values = metric_shapes_filtered,
                          labels = metric_labels_filtered,
                          name = "Metric")
    }
  }
  
  # Formatting
  p <- p +
    scale_x_continuous(breaks = pretty_breaks(),
                       guide = guide_axis(angle = 45),
                       limits = c(-(0.5 * binwidth), NA)) +
    scale_y_continuous(breaks = pretty_breaks(),
                       expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Repeat Length",
         y = "Number of Reads",
         fill = "Analysis Range",
         title = title) +
    theme_bw(base_size = base_text_size) +
    guides(colour = guide_legend(override.aes = list(size = 3)))
  
  return(p)
}


#' Plot Vertical Scatter/Violin
#'
#' Samples on x-axis, repeat length on y-axis
#'
#' @param data Data frame with repeat counts
#' @param x_var Column for x-axis (typically file_stem)
#' @param y_var Column for y-axis (repeat length)
#' @param points Optional data frame with modal/mean points
#' @param title Plot title
#' @param base_text_size Base text size
#'
#' @return ggplot object
plot_scatter_vertical <- function(data,
                                  x_var,
                                  y_var,
                                  points = NULL,
                                  title = NULL,
                                  base_text_size = 12) {
  
  # Base plot
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.25,
                fill = "lightgrey", colour = "black") +
    geom_jitter(height = 0, width = 0.25, size = 0.5, alpha = 0.3)
  
  # Add points if provided
  if (!is.null(points)) {
    # Define metric colors and labels
    metric_colors <- c("modal" = "red", "mean" = "blue", "median" = "purple")
    metric_labels <- c("modal" = "Modal", "mean" = "Mean", "median" = "Median")
    
    # Filter to only metrics present in data
    present_metrics <- unique(points$name)
    metric_colors_filtered <- metric_colors[names(metric_colors) %in% present_metrics]
    metric_labels_filtered <- metric_labels[names(metric_labels) %in% present_metrics]
    
    p <- p +
      geom_point(data = points,
                 aes(colour = name),
                 size = 3) +
      geom_text_repel(data = points,
                      aes(label = round(.data[[y_var]], 1),
                          colour = name),
                      size = base_text_size / 3.5,
                      angle = 45) +
      scale_colour_manual(values = metric_colors_filtered,
                         labels = metric_labels_filtered,
                         name = "Metric")
  }
  
  # Formatting
  p <- p +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_y_continuous(breaks = pretty_breaks(n = 10)) +
    labs(x = "Sample",
         y = "Repeat Length",
         title = title) +
    theme_bw(base_size = base_text_size)
  
  return(p)
}


#' Plot Horizontal Scatter/Violin
#'
#' Repeat length on x-axis, single sample
#'
#' @param data Data frame with repeat counts
#' @param x_var Column for x-axis (repeat length)
#' @param points Optional data frame with modal/mean points
#' @param boxes Optional data frame with range boxes
#' @param range_colors Named vector of colors for ranges
#' @param title Plot title
#' @param base_text_size Base text size
#'
#' @return ggplot object
plot_scatter_horizontal <- function(data,
                                    x_var,
                                    points = NULL,
                                    boxes = NULL,
                                    range_colors = NULL,
                                    title = NULL,
                                    base_text_size = 12) {
  
  # Base plot
  p <- ggplot(data, aes(x = .data[[x_var]], y = "")) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.25,
                fill = "lightgrey", colour = "black") +
    # geom_jitter(width = 0, height = 0.1, size = 1, alpha = 0.3)
    geom_jitter(width = 0, size = 1, alpha = 0.3)
  
  # Add boxes if provided
  if (!is.null(boxes)) {
    p <- p +
      geom_rect(data = boxes,
                aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                    fill = name),
                inherit.aes = FALSE,
                colour = NA,
                alpha = 0.1)
    
    # Apply consistent colors if provided
    if (!is.null(range_colors)) {
      p <- p + scale_fill_manual(values = range_colors, name = "Analysis Range")
    } else {
      p <- p + scale_fill_discrete(name = "Analysis Range")
    }
  }
  
  # Add new colour scale for points
  if (!is.null(points) && !is.null(boxes)) {
    p <- p + new_scale_colour()
  }
  
  # Add points if provided
  if (!is.null(points)) {
    # Define metric colors and labels
    metric_colors <- c("modal" = "red", "mean" = "blue", "median" = "purple")
    metric_labels <- c("modal" = "Modal", "mean" = "Mean", "median" = "Median")
    
    # Filter to only metrics present in data
    present_metrics <- unique(points$name)
    metric_colors_filtered <- metric_colors[names(metric_colors) %in% present_metrics]
    metric_labels_filtered <- metric_labels[names(metric_labels) %in% present_metrics]
    
    p <- p +
      geom_point(data = points,
                 aes(x = x, y = 1, colour = name),
                 size = 3,
                 inherit.aes = FALSE) +
      geom_text_repel(data = points,
                      aes(x = x, y = 1, 
                          label = paste0(name, ": ", round(x, 1)),
                          colour = name),
                      size = base_text_size / 3.5,
                      inherit.aes = FALSE,
                      show.legend = FALSE,
                      direction = "y",
                      nudge_y = 0.2) +
      scale_colour_manual(values = metric_colors_filtered,
                         labels = metric_labels_filtered,
                         name = "Metric")
  }
  
  # Formatting
  p <- p +
    scale_x_continuous(breaks = pretty_breaks(n = 10),
                       guide = guide_axis(angle = 45)) +
    labs(x = "Repeat Length",
         y = NULL,
         title = title) +
    theme_bw(base_size = base_text_size) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  return(p)
}
