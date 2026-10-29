# ==============================================================================
# Duke - range analysis and instability metrics
# ==============================================================================
# Definitions (x = repeat lengths of the reads or molecules in a range,
# s = setpoint, d = x - s, n = number of reads or molecules):
#
#   instability_index  mean(d)                         (Lee et al. 2010; signed)
#   expansion_index    sum(d[d > 0]) / n               (shares the n denominator,
#   contraction_index  sum(d[d < 0]) / n                 so II = EI + CI)
#   mean_expansion_size    mean(d[d > 0])  - how large a typical expansion is
#   mean_contraction_size  mean(d[d < 0])  - (these were Duke's EI/CI before v3)
#   prop_at / prop_above / prop_below   |d| <= 0.5, d > 0.5, d < -0.5
#   prop_above_k / prop_below_k        d >= k, d <= -k   (k = expansion_thresholds)
#   expansion_ratio / contraction_ratio  reads above (below) / reads at the setpoint
#       Kept for continuity. Its denominator is a single bin, typically ~6% of
#       reads in an expanded allele, so it is noisy; prop_above_k is the stable
#       alternative.
#   modal_diff / median_diff / mean_diff  sample statistic - s
#
# Setpoints: sample_relative uses the sample's own control_setpoint_metric;
# control_relative uses the group control (or a manual control length from the
# settings file, which takes precedence when given).
# ==============================================================================

parse_bracketed <- function(x) {
  x <- as.character(x)
  if (is.na(x) || !nzchar(trimws(x))) return(NULL)
  if (grepl("\\[", x)) {
    v <- regmatches(x, gregexpr("\\[([^]]*)\\]", x))[[1]]
    gsub("\\[|\\]", "", v)
  } else trimws(x)
}

parse_range_bounds <- function(x) {
  b <- strsplit(x, "-", fixed = TRUE)[[1]]
  if (length(b) != 2) stop("Range '", x, "' must be min-max, e.g. 36-NA")
  b <- trimws(b)
  open <- b == "" | toupper(b) == "NA"
  v <- suppressWarnings(as.numeric(b))
  if (any(is.na(v) & !open)) stop("Range '", x, "' has a bound that is not a number or NA")
  c(lo = ifelse(open[1], -Inf, v[1]), hi = ifelse(open[2], Inf, v[2]))
}

load_settings <- function(path, samples_in_data) {
  if (is.na(path) || !file.exists(path)) stop("Settings file not found: ", path)
  ext <- tolower(tools::file_ext(path))
  s <- if (ext == "csv") utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
       else as.data.frame(readxl::read_excel(path, sheet = 1))
  if (!"file_name" %in% names(s)) stop("Settings file needs a file_name column")
  s$file_stem <- sub("\\.(bam|fastq|fq|fasta|fa)(\\.gz)?$", "", s$file_name)
  for (col in c("group", "time", "group_control_sample", "exclude")) if (!col %in% names(s)) s[[col]] <- NA
  s$group <- dplyr::na_if(trimws(as.character(s$group)), "")
  s$time <- suppressWarnings(as.numeric(s$time))
  s$is_control <- truthy(s$group_control_sample)
  s$is_excluded <- truthy(s$exclude)
  s <- s[!duplicated(s$file_stem), , drop = FALSE]
  if (!any(s$file_stem %in% samples_in_data)) {
    ex <- function(v) paste(utils::head(v, 3), collapse = ", ")
    stop("None of the ", nrow(s), " file names in the settings file match the ", length(samples_in_data),
         " samples in dir_data.\n  Settings file_name (first 3): ", ex(s$file_name),
         "\n  Samples in the data (first 3): ", ex(samples_in_data),
         "\n  file_name must be the input file name, with or without its extension.")
  }
  list(all = s, used = s[s$file_stem %in% samples_in_data & !s$is_excluded, , drop = FALSE],
       excluded = s[s$file_stem %in% samples_in_data & s$is_excluded, , drop = FALSE],
       missing = setdiff(samples_in_data, s$file_stem))
}

# Long range table: one row per sample x measure x range.
# For the primary measure the columns are analysis_ranges / floor / max_peaks /
# manual_control_repeat_length. Another measure may have its own columns with a
# suffix (analysis_ranges_CCG, floor_CCG, ...), in which case its ranges are on
# its own scale. Without them it uses the primary measure's ranges: the reads
# are selected by their primary length (e.g. the expanded CAG allele) and the
# metrics are calculated on the other measure.
build_range_table <- function(settings, measures, primary) {
  rows <- list()
  for (m in measures) {
    get_col <- function(nm) if (nm %in% names(settings)) settings[[nm]] else rep(NA, nrow(settings))
    ar0 <- get_col("analysis_ranges"); fl0 <- get_col("floor"); mp0 <- get_col("max_peaks")
    arm <- get_col(paste0("analysis_ranges_", m)); flm <- get_col(paste0("floor_", m)); mpm <- get_col(paste0("max_peaks_", m))
    mc <- if (m == primary) get_col("manual_control_repeat_length") else get_col(paste0("manual_control_repeat_length_", m))
    for (i in seq_len(nrow(settings))) {
      # a measure's own ranges when its column is filled for this sample, else the primary ranges
      own <- m != primary && !is.null(parse_bracketed(arm[i]))
      ar <- if (own) arm else ar0; fl <- if (own) flm else fl0; mp <- if (own) mpm else mp0
      r <- parse_bracketed(ar[i])
      if (is.null(r)) next
      f <- parse_bracketed(fl[i]); k <- parse_bracketed(mp[i]); c <- parse_bracketed(mc[i])
      n <- length(r)
      pad <- function(v, default) { v <- suppressWarnings(as.numeric(v)); if (!length(v)) rep(default, n) else if (length(v) != n) stop("Sample ", settings$file_stem[i], ": bracketed columns for ", m, " have different numbers of ranges") else v }
      b <- t(vapply(r, parse_range_bounds, numeric(2)))
      width <- nchar(as.character(n))
      rows[[length(rows) + 1]] <- data.frame(
        file_stem = settings$file_stem[i], measure = m, range_measure = if (own || m == primary) m else primary,
        range_name = paste0("range", formatC(seq_len(n), width = max(2, width), flag = "0")),
        range_label = r, lo = b[, "lo"], hi = b[, "hi"],
        floor = pad(f, NA_real_), max_peaks = pad(k, 1), manual_control = pad(c, NA_real_),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) stop("No analysis ranges could be read from the settings file")
  out <- do.call(rbind, rows)
  flo <- is.finite(out$lo); fhi <- is.finite(out$hi)
  out$range_text <- ifelse(!flo & !fhi, "all",
                    ifelse(!flo, paste0("<=", out$hi),
                    ifelse(!fhi, paste0(">=", out$lo), paste0(out$lo, "-", out$hi))))
  out$range_text <- ifelse(out$range_measure == out$measure, out$range_text,
                           paste(out$range_measure, out$range_text))
  out
}

# Remove repeat lengths whose bin holds fewer than `floor` reads (or, when
# floor < 1, less than that fraction of the modal bin). Note this trims the
# sparse expansion tail, increasingly so at lower depth; floor = 1 keeps all.
apply_floor <- function(x, floor) {
  x <- x[is.finite(x)]
  if (!length(x) || is.na(floor) || floor <= 0) return(x)
  tb <- table(x)
  thr <- if (floor < 1) floor * max(tb) else floor
  keep <- as.numeric(names(tb))[tb >= thr]
  x[x %in% keep]
}

.skew <- function(x) { m <- mean(x); s <- sqrt(mean((x - m)^2)); if (s == 0) NA_real_ else mean(((x - m) / s)^3) }
.kurt <- function(x) { m <- mean(x); s <- sqrt(mean((x - m)^2)); if (s == 0) NA_real_ else mean(((x - m) / s)^4) }
.mode <- function(x) { tb <- table(x); as.numeric(names(tb)[which.max(tb)]) }

distribution_stats <- function(x) {
  n <- length(x)
  if (n == 0) return(data.frame(n_used = 0L, modal_length = NA_real_, modal_frequency = NA_real_, mean_length = NA_real_,
                                median_length = NA_real_, sd = NA_real_, iqr = NA_real_, p05 = NA_real_, p95 = NA_real_,
                                min_length = NA_real_, max_length = NA_real_, cv = NA_real_, skewness = NA_real_,
                                kurtosis = NA_real_, prop_below_modal = NA_real_, prop_above_modal = NA_real_))
  mo <- .mode(x)
  data.frame(n_used = n, modal_length = mo, modal_frequency = sum(x == mo), mean_length = mean(x),
             median_length = stats::median(x), sd = if (n > 1) stats::sd(x) else NA_real_,
             iqr = stats::IQR(x), p05 = unname(stats::quantile(x, 0.05)), p95 = unname(stats::quantile(x, 0.95)),
             min_length = min(x), max_length = max(x),
             cv = if (mean(x) != 0 && n > 1) stats::sd(x) / mean(x) else NA_real_,
             skewness = if (n > 2) .skew(x) else NA_real_, kurtosis = if (n > 3) .kurt(x) else NA_real_,
             prop_below_modal = mean(x < mo), prop_above_modal = mean(x > mo))
}

instability_stats <- function(x, s, thresholds = c(5, 10)) {
  n <- length(x)
  base <- data.frame(setpoint = s, n_used = n, instability_index = NA_real_, expansion_index = NA_real_,
                     contraction_index = NA_real_, mean_expansion_size = NA_real_, mean_contraction_size = NA_real_,
                     prop_at = NA_real_, prop_above = NA_real_, prop_below = NA_real_,
                     reads_at = NA_integer_, reads_above = NA_integer_, reads_below = NA_integer_,
                     expansion_ratio = NA_real_, contraction_ratio = NA_real_)
  for (k in thresholds) { base[[paste0("prop_above_", k)]] <- NA_real_; base[[paste0("prop_below_", k)]] <- NA_real_ }
  if (n == 0 || is.na(s)) return(base)
  d <- x - s
  at <- abs(d) <= 0.5; up <- d > 0.5; dn <- d < -0.5
  base$instability_index <- mean(d)
  base$expansion_index <- sum(d[d > 0]) / n
  base$contraction_index <- sum(d[d < 0]) / n
  base$mean_expansion_size <- if (any(d > 0)) mean(d[d > 0]) else NA_real_
  base$mean_contraction_size <- if (any(d < 0)) mean(d[d < 0]) else NA_real_
  base$prop_at <- mean(at); base$prop_above <- mean(up); base$prop_below <- mean(dn)
  base$reads_at <- sum(at); base$reads_above <- sum(up); base$reads_below <- sum(dn)
  base$expansion_ratio <- if (sum(at) > 0) sum(up) / sum(at) else NA_real_
  base$contraction_ratio <- if (sum(at) > 0) sum(dn) / sum(at) else NA_real_
  for (k in thresholds) { base[[paste0("prop_above_", k)]] <- mean(d >= k); base[[paste0("prop_below_", k)]] <- mean(d <= -k) }
  base
}

# Modal peaks inside one range (local maxima of the frequency table)
find_modal_peaks <- function(x, max_peaks = 1, span = 3) {
  if (!length(x)) return(data.frame(peak_rank = NA_integer_, peak_length = NA_real_, peak_frequency = NA_real_))
  tb <- table(x); len <- as.numeric(names(tb)); fr <- as.numeric(tb)
  if (length(fr) < span || length(fr) < 3) {
    i <- which.max(fr); return(data.frame(peak_rank = 1L, peak_length = len[i], peak_frequency = fr[i]))
  }
  pk <- pracma::findpeaks(fr, minpeakdistance = 1, npeaks = max(1, max_peaks), sortstr = TRUE)
  if (is.null(pk)) { i <- which.max(fr); return(data.frame(peak_rank = 1L, peak_length = len[i], peak_frequency = fr[i])) }
  pk <- matrix(pk, ncol = 4)
  data.frame(peak_rank = seq_len(nrow(pk)), peak_length = len[pk[, 2]], peak_frequency = pk[, 1])
}

aggregate_controls <- function(values, method) {
  values <- values[is.finite(values)]
  if (!length(values)) return(NA_real_)
  switch(method, mean = mean(values), median = stats::median(values), trimmed_mean = mean(values, trim = 0.1),
         stop("Unknown control_aggregation_method: ", method))
}
