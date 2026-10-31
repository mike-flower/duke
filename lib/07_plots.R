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

# Length distribution of reads and molecules drawn in separate rows, each on its own
# y-axis: bars for the share of the unit at each length, and a line for a kernel
# density (bandwidth `bw` units; by default 0.5 for short tracts, 1 otherwise) scaled to the same bars, so that a sparse
# distribution (molecules) is read from its smooth shape rather than from the
# bin-to-bin noise of a few thousand values. `d` has columns length and unit.
dist_rows <- function(d, title = NULL, xlab = "Repeat length", binwidth = 1, bw = NULL, xlim = NULL,
                      unit_cols = c(Reads = "grey45", Molecules = "#00838f"), marks = NULL) {
  d <- data.table::as.data.table(d)[is.finite(length)]
  d[, unit := factor(unit, levels = intersect(names(unit_cols), unique(unit)))]
  # A short tract (a few units, e.g. CCG) would be flattened by a 1-unit bandwidth
  if (is.null(bw)) bw <- if (stats::median(d$length) < 30) 0.5 else 1
  if (is.null(xlim)) xlim <- c(min(d$length), stats::quantile(d$length, 0.999)) + c(-2, 2) * binwidth
  h <- d[, .N, by = .(unit, bin = floor(length / binwidth + 0.5) * binwidth)][, share := N / sum(N), by = unit]
  den <- d[, {
    if (.N >= 2) { k <- stats::density(length, bw = bw, from = xlim[1], to = xlim[2], n = 1024)
                   .(x = k$x, y = k$y * binwidth) } else .(x = numeric(0), y = numeric(0))
  }, by = unit]
  n_lab <- d[, .(label = sprintf("%s (%s)", unit, format(.N, big.mark = ","))), by = unit]
  lab <- stats::setNames(n_lab$label, as.character(n_lab$unit))
  p <- ggplot2::ggplot(h, ggplot2::aes(bin, share)) +
    ggplot2::geom_col(ggplot2::aes(fill = unit), width = binwidth, alpha = 0.35) +
    ggplot2::geom_line(data = den, ggplot2::aes(x, y, colour = unit), linewidth = 0.7) +
    ggplot2::facet_wrap(~ unit, ncol = 1, scales = "free_y", labeller = ggplot2::as_labeller(lab)) +
    ggplot2::scale_fill_manual(values = unit_cols, guide = "none") + ggplot2::scale_colour_manual(values = unit_cols, guide = "none") +
    ggplot2::scale_y_continuous(labels = scales::label_percent(), expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::labs(x = xlab, y = "Share at each length", title = title) + ggplot2::theme_bw()
  if (!is.null(marks) && nrow(marks)) {
    # marks: unit, stat, value - drawn per row, so read and molecule statistics are each beside their own data
    marks <- data.table::as.data.table(marks)[, unit := factor(unit, levels = levels(d$unit))]
    ytop <- h[, .(top = max(share)), by = unit]
    marks <- ytop[marks, on = "unit"][, y := top * c(modal = 0.95, median = 0.78, mean = 0.61)[as.character(stat)]]
    cols <- c(modal = "#e41a1c", median = "#984ea3", mean = "#ff7f00")
    p <- p +
      ggplot2::geom_vline(data = marks, ggplot2::aes(xintercept = value), colour = cols[as.character(marks$stat)], linewidth = 0.4, linetype = "dashed") +
      ggplot2::geom_label(data = marks, ggplot2::aes(x = value, y = y, label = sprintf("%s %.1f", stat, value)), colour = cols[as.character(marks$stat)],
                          inherit.aes = FALSE, size = 2.6, hjust = -0.05, label.size = 0, fill = "white", alpha = 0.85)
  }
  p
}

# Ridge plot of several samples (or one): each sample's read length distribution as a
# grey density (bandwidth 1 unit, scaled to its own height), the molecule density as a
# teal line on the same row when there are molecules, every read (up to dots_max per
# sample, drawn at random) as a dot just below the ridge, and the reads' modal length
# as a red dashed line with its value. `rd` and `md` have columns label and length;
# label is a factor whose first level is drawn at the top.
ridge_dots_plot <- function(rd, md = NULL, title = NULL, xlab = "Repeat length", xlim = NULL, dots_max = 2000L, subtitle = NULL) {
  rd <- data.table::as.data.table(rd)[is.finite(length)]
  if (!nrow(rd)) return(NULL)
  L <- levels(droplevels(rd$label)); if (is.null(L)) L <- unique(as.character(rd$label))
  idx <- stats::setNames(rev(seq_along(L)), L)
  if (is.null(xlim)) xlim <- c(min(rd$length), stats::quantile(rd$length, 0.999)) + c(-2, 2)
  dens <- function(d) d[, {
    if (.N >= 2 && diff(range(length)) > 0) { k <- stats::density(length, bw = 1, from = xlim[1], to = xlim[2], n = 512); .(x = k$x, y = k$y / max(k$y))
    } else .(x = rep(length[1], 2), y = c(0, 1))
  }, by = label][, base := idx[as.character(label)]][, top := base + 0.85 * y]
  rdn <- dens(rd)
  set.seed(123)
  # dots jittered sideways by up to 0.35 units, so reads of one whole length do not stack into a bar
  dots <- rd[, .SD[sample(.N, min(.N, dots_max))], by = label][, `:=`(y = idx[as.character(label)] - 0.06 - stats::runif(.N, 0, 0.26),
                                                                      length = length + stats::runif(.N, -0.35, 0.35))]
  modal <- rd[, .(modal = { tb <- table(round(length)); as.numeric(names(tb)[which.max(tb)]) }), by = label][, base := idx[as.character(label)]]
  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = rdn, ggplot2::aes(x, ymin = base, ymax = top, group = label), fill = "grey82", colour = "grey25", linewidth = 0.3) +
    ggplot2::geom_point(data = dots, ggplot2::aes(length, y), size = 0.35, alpha = 0.3, colour = "black", stroke = 0)
  if (!is.null(md) && nrow(md)) {
    md <- data.table::as.data.table(md)[is.finite(length) & label %in% L]
    if (nrow(md)) p <- p + ggplot2::geom_line(data = dens(md), ggplot2::aes(x, top, group = label, colour = "Molecules"), linewidth = 0.5) +
      ggplot2::scale_colour_manual(values = c(Molecules = "#00838f"), name = NULL)
  }
  p + ggplot2::geom_segment(data = modal, ggplot2::aes(x = modal, xend = modal, y = base, yend = base + 0.9), colour = "red", linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_text(data = modal, ggplot2::aes(x = modal, y = base + 0.9, label = round(modal)), colour = "red", size = 2.6, hjust = -0.15, vjust = 1) +
    ggplot2::scale_y_continuous(breaks = unname(idx), labels = names(idx), expand = ggplot2::expansion(add = c(0.45, 0.25))) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(10)) +
    ggplot2::labs(x = xlab, y = NULL, title = title,
                  subtitle = subtitle %||% sprintf("Grey = reads (density, scaled to each sample's height)%s\nDots = %s reads; red = the reads' modal length",
                                                   if (!is.null(md) && nrow(md)) ", teal = molecules" else "",
                                                   if (rd[, max(.N), by = label][, max(V1)] > dots_max) paste("up to", format(dots_max, big.mark = ","), "random") else "all")) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank())
}

# Histograms of several samples: one row per sample, reads and (with UMIs) molecules
# in separate columns, each panel on its own y-axis as the share of that sample's
# reads (or molecules) at each length.
hist_rows_plot <- function(rd, md = NULL, title = NULL, xlab = "Repeat length", xlim = NULL, binwidth = 1) {
  d <- rbind(data.table::as.data.table(rd)[, .(label, length, unit = "Reads")],
             if (!is.null(md) && nrow(md)) data.table::as.data.table(md)[, .(label, length, unit = "Molecules")])
  d <- d[is.finite(length)]
  if (!nrow(d)) return(NULL)
  d[, unit := factor(unit, levels = c("Reads", "Molecules"))]
  if (is.null(xlim)) xlim <- c(min(d$length), stats::quantile(d$length, 0.999)) + c(-2, 2) * binwidth
  h <- d[, .N, by = .(label, unit, bin = floor(length / binwidth + 0.5) * binwidth)][, share := N / sum(N), by = .(label, unit)]
  h[, panel := factor(paste0(label, " - ", unit), levels = as.vector(t(outer(levels(factor(h$label)), levels(h$unit), paste, sep = " - "))))]
  ggplot2::ggplot(h, ggplot2::aes(bin, share, fill = unit)) + ggplot2::geom_col(width = binwidth) +
    ggplot2::facet_wrap(~ panel, ncol = data.table::uniqueN(h$unit), scales = "free_y", drop = TRUE) +
    ggplot2::scale_fill_manual(values = c(Reads = "grey45", Molecules = "#00838f"), guide = "none") +
    ggplot2::scale_y_continuous(labels = scales::label_percent(), n.breaks = 3) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::labs(x = xlab, y = "Share at each length", title = title,
                  subtitle = "One row per sample; each panel is the share of that sample's reads (grey) or molecules (teal) at each length, on its own scale") +
    ggplot2::theme_bw() + ggplot2::theme(strip.text = ggplot2::element_text(size = 7, margin = ggplot2::margin(1, 1, 1, 1)),
                                         panel.spacing.y = grid::unit(2, "pt"), axis.text.y = ggplot2::element_text(size = 6))
}
