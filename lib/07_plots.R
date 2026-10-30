# ==============================================================================
# Duke - repeat-length figures (Module 7)
# ==============================================================================

# Wrap a long label for a facet strip, breaking after "-", "_", "." or space
# where possible, so that names such as "compare_gt-QS5---ipsc-cell" stay legible.
wrap_label <- function(x, width = 14) {
  vapply(x, function(s) {
    tok <- regmatches(s, gregexpr("[^-_. ]*[-_. ]*", s))[[1]]; tok <- tok[nzchar(tok)]
    lines <- character(0); cur <- ""
    for (t in tok) {
      while (nchar(t) > width) { if (nzchar(cur)) { lines <- c(lines, cur); cur <- "" }; lines <- c(lines, substr(t, 1, width)); t <- substring(t, width + 1) }
      if (nchar(cur) + nchar(t) > width && nzchar(cur)) { lines <- c(lines, cur); cur <- t } else cur <- paste0(cur, t)
    }
    paste(c(lines, if (nzchar(cur)) cur), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

# Heatmap: one column per sample, one row per length bin, fill = share of the
# sample's reads (or molecules). Samples are grouped into panels by group.
# Lengths above the cohort's 99.9th percentile (rare outliers) are not drawn, so the
# alleles fill the panel; the colour is on a log scale from 0.01%, so a minor allele
# or a thin expansion tail is as visible as the main peak, and lengths with no reads
# stay white.
length_heatmap <- function(d, title, binwidth = 1, ylab = "Repeat length", unit = "reads") {
  d <- data.table::as.data.table(d)
  ymax <- stats::quantile(d$length, 0.999, na.rm = TRUE)
  n_out <- sum(d$length > ymax, na.rm = TRUE)
  d[, bin := floor(length / binwidth) * binwidth]
  h <- d[, .N, by = .(sample_label, group_label, bin)][, prop := N / sum(N), by = sample_label][bin <= ymax]
  ggplot2::ggplot(h, ggplot2::aes(sample_label, bin, fill = prop)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(option = "plasma", labels = scales::label_percent(), name = "Share\n(log scale)",
                                  trans = "log10", limits = c(1e-4, max(h$prop)), oob = scales::squish) +
    ggplot2::facet_grid(~ group_label, scales = "free_x", space = "free_x", labeller = ggplot2::as_labeller(wrap_label)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = scales::pretty_breaks(10)) +
    ggplot2::labs(x = NULL, y = ylab, title = title,
                  subtitle = sprintf("Lengths above %s (the top 0.1%%, %s %s) not drawn", round(ymax), format(n_out, big.mark = ","), unit)) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
                   panel.spacing = grid::unit(2, "pt"), strip.text = ggplot2::element_text(size = 7),
                   panel.grid = ggplot2::element_blank())
}

# One histogram panel: reads (grey bars) and, with UMIs, molecules (teal line
# scaled to the same total); modal / median / mean marked and labelled.
length_panel <- function(reads, molecules = NULL, marks = NULL, title = NULL, xlab = "Repeat length",
                         binwidth = 1, sqrt_y = FALSE) {
  p <- ggplot2::ggplot(data.frame(length = reads), ggplot2::aes(length)) +
    ggplot2::geom_histogram(binwidth = binwidth, fill = "grey60", colour = NA, boundary = -0.5 * binwidth)
  if (!is.null(molecules) && length(molecules)) {
    b <- floor(molecules / binwidth + 0.5)
    mh <- data.frame(length = seq(min(b) - 1, max(b) + 1) * binwidth)
    mh$Freq <- tabulate(b - min(b) + 2, nbins = nrow(mh))
    mh$scaled <- mh$Freq * length(reads) / length(molecules)
    p <- p + ggplot2::geom_step(data = mh, ggplot2::aes(length, scaled), colour = "#00838f", direction = "mid", linewidth = 0.5)
  }
  if (!is.null(marks) && nrow(marks)) {
    cols <- c(modal = "#e41a1c", median = "#984ea3", mean = "#ff7f00")
    ymax <- max(table(round(reads / binwidth)))
    marks$y <- ymax * c(modal = 0.95, median = 0.8, mean = 0.65)[as.character(marks$stat)]
    p <- p + ggplot2::geom_vline(data = marks, ggplot2::aes(xintercept = value, colour = stat), linewidth = 0.4, linetype = "dashed") +
      ggplot2::geom_label(data = marks, ggplot2::aes(x = value, y = y, label = sprintf("%s %.1f", stat, value), colour = stat),
                          inherit.aes = FALSE, size = 2.6, hjust = -0.05, label.size = 0, fill = "white", alpha = 0.8) +
      ggplot2::scale_colour_manual(values = cols, guide = "none")
  }
  # The x-axis stops at the 99.9th percentile of the reads, so a few very long reads do
  # not squash the distribution into the left of the panel
  xr <- c(min(reads, na.rm = TRUE), stats::quantile(reads, 0.999, na.rm = TRUE)) + c(-2, 2) * binwidth
  p + ggplot2::coord_cartesian(xlim = xr) +
    (if (sqrt_y) ggplot2::scale_y_sqrt(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.05)))
       else ggplot2::scale_y_continuous(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.05)))) +
    ggplot2::labs(x = xlab, y = paste0("Reads", if (sqrt_y) " (square-root scale)" else "",
                                       if (!is.null(molecules) && length(molecules)) "\nmolecules scaled to the same total" else ""),
                  title = title) +
    ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
}
