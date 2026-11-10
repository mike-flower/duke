# ==============================================================================
# Duke - general utilities
# ==============================================================================
# Excel export, threading, timing, safe cache loading, path sanitising,
# parameter hashing, house-style tables and plot saving.
# ==============================================================================

# The run's parameter list as seen from the caller: only a list called params counts (a
# function of that name elsewhere on the search path is skipped), NULL when there is none
.run_params <- function(env) get0("params", envir = env, mode = "list", inherits = TRUE)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

# ------------------------------------------------------------------------------
# Excel export
# ------------------------------------------------------------------------------
# openxlsx (4.2.5) links every sheet to a drawing part (xl/drawings/drawing1.xml) that it does
# not write. Excel ignores the broken link, but Python (openpyxl, pandas) refuses the file. The
# link is removed after saving; if the zip tool is missing the file is left as openxlsx wrote it.
drop_missing_drawings <- function(path) {
  if (!nzchar(Sys.which("zip"))) return(invisible(FALSE))
  path <- normalizePath(path)
  tmp <- tempfile("xlsx_"); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  if (is.null(tryCatch(utils::unzip(path, exdir = tmp), error = function(e) NULL))) return(invisible(FALSE))
  ws <- file.path(tmp, "xl", "worksheets")
  changed <- FALSE
  for (r in list.files(file.path(ws, "_rels"), pattern = "\\.rels$", full.names = TRUE)) {
    x <- paste(readLines(r, warn = FALSE), collapse = "")
    rel <- regmatches(x, gregexpr("<Relationship [^>]*/>", x))[[1]]
    tgt <- sub('.*Target="([^"]+)".*', "\\1", rel)
    gone <- rel[!grepl('TargetMode="External"', rel) & !file.exists(file.path(ws, tgt))]
    if (!length(gone)) next
    for (g in gone) x <- sub(g, "", x, fixed = TRUE)
    writeLines(x, r)
    sh <- file.path(ws, sub("\\.rels$", "", basename(r)))
    if (file.exists(sh)) {
      y <- paste(readLines(sh, warn = FALSE), collapse = "\n")
      for (id in sub('.*Id="([^"]+)".*', "\\1", gone)) y <- gsub(sprintf('<(drawing|legacyDrawing) r:id="%s"/>', id), "", y)
      writeLines(y, sh)
    }
    changed <- TRUE
  }
  if (!changed) return(invisible(FALSE))
  out <- tempfile(fileext = ".xlsx")
  wd <- setwd(tmp); on.exit(setwd(wd), add = TRUE)
  fl <- list.files(".", recursive = TRUE, all.files = TRUE)
  ok <- utils::zip(out, c("[Content_Types].xml", setdiff(fl, "[Content_Types].xml")), flags = "-q -X") == 0
  if (ok && file.exists(out)) file.copy(out, path, overwrite = TRUE)
  invisible(ok)
}

write_excel_function <- function(named_list, path, overwrite_sheet = TRUE) {
  named_list <- Filter(function(x) !is.null(x) && is.data.frame(x), named_list)
  if (length(named_list) == 0) return(invisible(NULL))
  # Re-saving a workbook opened with openxlsx::loadWorkbook leaves links to drawing parts
  # that are not in the file, which other readers (Python openpyxl/pandas) refuse to open.
  # So the existing sheets are read back as plain tables and the workbook is rebuilt.
  wb <- openxlsx::createWorkbook()
  if (file.exists(path)) {
    for (sh in openxlsx::getSheetNames(path)) {
      if (overwrite_sheet && sh %in% vapply(names(named_list), function(n) if (nchar(n) > 31) paste0(substr(n, 1, 29), "..") else n, "")) next
      old <- tryCatch(openxlsx::read.xlsx(path, sheet = sh, colNames = TRUE, check.names = FALSE, sep.names = " "),
                      error = function(e) NULL)
      # A sheet with headers but no rows comes back empty: keep its header row
      if (is.null(old) || !ncol(old)) old <- tryCatch(openxlsx::read.xlsx(path, sheet = sh, colNames = FALSE, rows = 1),
                                                      error = function(e) NULL)
      openxlsx::addWorksheet(wb, sh)
      if (!is.null(old)) openxlsx::writeData(wb, sh, old, colNames = !is.null(names(old)) && !all(grepl("^X\\d+$", names(old))))
    }
  }
  truncate_name <- function(name, max_length = 31) {
    if (nchar(name) > max_length) paste0(substr(name, 1, max_length - 2), "..") else name
  }
  for (table_name in names(named_list)) {
    export_name <- truncate_name(table_name)
    if (export_name %in% names(wb)) {
      if (overwrite_sheet) {
        openxlsx::removeWorksheet(wb, export_name)
      } else {
        suffix <- 1
        while (paste0(export_name, "_", suffix) %in% names(wb)) suffix <- suffix + 1
        export_name <- paste0(export_name, "_", suffix)
      }
    }
    openxlsx::addWorksheet(wb, export_name)
    openxlsx::writeData(wb, export_name, as.data.frame(named_list[[table_name]]))
  }
  openxlsx::saveWorkbook(wb, file = path, overwrite = TRUE)
  drop_missing_drawings(path)
  invisible(path)
}

# ------------------------------------------------------------------------------
# Threading wrapper. Uses params$threads from the calling environment.
# ------------------------------------------------------------------------------
apply_fn <- function(X, FUN, ..., threads = NULL) {
  threads <- threads %||% tryCatch(.run_params(parent.frame())$threads,
                                   error = function(e) 1)
  threads <- threads %||% 1
  lost <- character(0)
  if (threads > 1 && length(X) > 1) {
    # A worker process that is killed (usually for running out of memory) returns NULL for
    # every item it held, with only a warning; catch that warning so the run stops instead
    # of silently dropping those samples
    res <- withCallingHandlers(
      pbmcapply::pbmclapply(X, FUN, ..., mc.cores = min(threads, length(X))),
      warning = function(w) {
        if (grepl("did not deliver|core.*(killed|died)", conditionMessage(w))) {
          lost <<- c(lost, conditionMessage(w)); invokeRestart("muffleWarning")
        }
      })
  } else {
    res <- pbapply::pblapply(X, FUN, ...)
  }
  res <- extract_apply_fn_result(res)
  if (length(lost) || any(grepl("did not deliver", attr(res, "warning") %||% "")))
    stop("A worker process stopped without returning its results (usually it ran out of memory), so ",
         "some samples would be missing. Give the job more memory or lower `threads`.\n  ",
         paste(unique(lost), collapse = "\n  "))
  # Surface worker errors rather than letting them become silent NULLs
  errs <- vapply(res, function(r) inherits(r, "try-error"), logical(1))
  if (any(errs)) stop("Worker failed for ", sum(errs), " item(s); first error:\n  ",
                      as.character(res[[which(errs)[1]]]))
  if (!is.null(names(X))) names(res) <- names(X) else if (is.character(X)) names(res) <- X
  res
}

extract_apply_fn_result <- function(result) {
  if (is.list(result) && "value" %in% names(result) && "warning" %in% names(result)) {
    if (!is.null(result$warning)) {
      if (grepl("did not deliver", result$warning$message))
        stop("A worker process stopped without returning its results (usually it ran out of memory), so ",
             "some samples would be missing. Give the job more memory or lower `threads`.\n  ", result$warning$message)
      warning(result$warning$message)
    }
    result <- result$value
  }
  result
}

# ------------------------------------------------------------------------------
# Timing
# ------------------------------------------------------------------------------
format_size <- function(bytes) {
  bytes <- as.numeric(bytes)
  if (is.na(bytes)) return(NA_character_)
  if (bytes < 1024) sprintf("%.0f B", bytes)
  else if (bytes < 1024^2) sprintf("%.1f KB", bytes / 1024)
  else if (bytes < 1024^3) sprintf("%.1f MB", bytes / 1024^2)
  else sprintf("%.2f GB", bytes / 1024^3)
}

format_elapsed <- function(secs) {
  secs <- round(as.numeric(secs))
  if (secs < 3600) sprintf("%02d:%02d", secs %/% 60, secs %% 60)
  else sprintf("%02d:%02d:%02d", secs %/% 3600, (secs %% 3600) %/% 60, secs %% 60)
}

build_timing_table <- function(timings) {
  if (length(timings) == 0) timings <- c(None = 0)
  elapsed_secs <- as.numeric(timings)
  total <- sum(elapsed_secs)
  data.frame(
    section    = c(names(timings), "Total"),
    secs       = c(round(elapsed_secs, 1), round(total, 1)),
    elapsed    = c(vapply(elapsed_secs, format_elapsed, ""), format_elapsed(total)),
    cumulative = c(vapply(cumsum(elapsed_secs), format_elapsed, ""), format_elapsed(total)),
    stringsAsFactors = FALSE
  )
}

build_manifest <- function(output_list, descriptions) {
  data.frame(
    object = names(output_list),
    description = unname(descriptions[names(output_list)] %||% NA_character_),
    class = vapply(output_list, function(x) if (is.null(x)) "NULL" else class(x)[1], ""),
    size = vapply(output_list, function(x) format_size(if (is.null(x)) 0 else as.numeric(object.size(x))), ""),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# Safe cache loading
# ------------------------------------------------------------------------------
safe_load_cache <- function(path, envir = parent.frame()) {
  if (!file.exists(path)) return(FALSE)
  if (file.size(path) == 0) { file.remove(path); return(FALSE) }
  tryCatch({ load(path, envir = envir); TRUE },
           error = function(e) { warning("Cache file is corrupt, will recompute: ", basename(path))
                                 file.remove(path); FALSE })
}

safe_load_module <- function(path, module_num, envir = parent.frame()) {
  if (!file.exists(path)) {
    stop("Module ", module_num, " output not found. Has Module ", module_num, " been run?\n  Expected: ", path)
  }
  if (file.size(path) == 0) {
    stop("Module ", module_num, " output is empty (zero bytes) - the previous run was likely ",
         "interrupted during save.\n  Delete and rerun Module ", module_num, ":\n  ", path)
  }
  tryCatch(load(path, envir = envir),
           error = function(e) stop("Module ", module_num, " output is corrupt. Delete and rerun Module ",
                                    module_num, ":\n  ", path, "\n  Original error: ", e$message))
}

# ------------------------------------------------------------------------------
# Parameter hashing for resume. A module's results are reused only when the
# parameters that affect it are unchanged; per-sample caches live in a temp
# sub-directory named by the hash, so changed settings never pick up old caches.
# ------------------------------------------------------------------------------
params_hash <- function(params, keys) {
  sub <- params[intersect(keys, names(params))]
  sub <- sub[order(names(sub))]
  substr(digest::digest(sub, algo = "md5"), 1, 10)
}

module_temp_dir <- function(dir_out, module_tag, hash) {
  d <- file.path(dir_out, "temp", paste0(module_tag, "_", hash))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}

# ------------------------------------------------------------------------------
# Flag parsing. Settings flags may be TRUE/FALSE, 1/0, yes/no or blank.
# ------------------------------------------------------------------------------
truthy <- function(x) {
  v <- tolower(trimws(as.character(x)))
  !is.na(v) & v %in% c("true", "t", "1", "yes", "y", "x")
}

# ------------------------------------------------------------------------------
# Filesystem-safe labels (display labels are never changed)
# ------------------------------------------------------------------------------
safe_path <- function(x, replacement = "-") {
  x <- as.character(x)
  # Characters that cannot be in a file name, and spaces, become the replacement; a run
  # of them becomes one ("a /  b" gives a-b). Separators already in the name are kept,
  # so a sample called bc1002--bc1050 keeps its two dashes.
  x <- gsub("[/\\\\:*?\"<>|#%[:cntrl:][:space:]]+", "\001", x, perl = TRUE)   # # and % break links in the HTML reports
  x <- gsub(paste0("\\Q", replacement, "\\E*\001[\001", replacement, "]*"), "\001", x, perl = TRUE)
  x <- gsub("\001", replacement, x, fixed = TRUE)
  x <- gsub(paste0("^\\Q", replacement, "\\E+|\\Q", replacement, "\\E+$"), "", x, perl = TRUE)
  x[is.na(x) | x == ""] <- "unnamed"
  x
}

warn_if_sanitised <- function(original, what = "label") {
  original <- unique(as.character(original))
  safe <- safe_path(original)
  changed <- original != safe
  if (any(changed)) {
    message("Note: ", sum(changed), " ", what, "(s) contain characters unsafe for filenames; ",
            "sanitised for output paths (display labels unchanged)")
  }
  dup <- unique(safe[duplicated(safe)])
  if (length(dup) > 0) {
    stop("Sanitising ", what, " names for output paths produced collisions: ",
         paste(vapply(dup, function(d) paste(original[safe == d], collapse = " / "), ""), collapse = "; "))
  }
  invisible(safe)
}

# ------------------------------------------------------------------------------
# House-style tables
# ------------------------------------------------------------------------------
dt_full <- function(df, caption = NULL, digits = NULL, page_length = 25) {
  df <- as.data.frame(df)
  dt <- DT::datatable(
    df, caption = caption, extensions = "Buttons",
    options = list(pageLength = page_length, scrollX = TRUE, scrollY = "500px",
                   dom = "Bfrtip", buttons = c("copy", "csv", "excel")),
    rownames = FALSE, filter = "top", class = "cell-border stripe compact",
    style = "bootstrap")
  dt <- DT::formatStyle(dt, columns = colnames(df), fontSize = "85%")
  if (!is.null(digits)) {
    num <- names(df)[vapply(df, function(x) is.numeric(x) && !is.integer(x), logical(1))]
    if (length(num)) dt <- DT::formatRound(dt, columns = num, digits = digits)
  }
  dt
}

dt_compact <- function(df, caption = NULL, digits = NULL) {
  df <- as.data.frame(df)
  dt <- DT::datatable(df, caption = caption,
                      options = list(dom = "t", pageLength = max(10, nrow(df)), scrollX = TRUE),
                      rownames = FALSE, class = "cell-border stripe compact")
  dt <- DT::formatStyle(dt, columns = colnames(df), fontSize = "85%")
  if (!is.null(digits)) {
    num <- names(df)[vapply(df, function(x) is.numeric(x) && !is.integer(x), logical(1))]
    if (length(num)) dt <- DT::formatRound(dt, columns = num, digits = digits)
  }
  dt
}

# ------------------------------------------------------------------------------
# Plot saving. Every plot is assigned, printed (when inline) and saved.
# ------------------------------------------------------------------------------
save_plot <- function(p, path, width = 10, height = 6, dpi = NULL) {
  prm <- .run_params(parent.frame())
  dpi <- dpi %||% prm$plot_dpi %||% 150
  # A plot whose type (the start of its file name) is not in `plots` is not saved
  ty <- plot_type_of(path)
  if (!is.na(ty) && !is.null(prm$plot_types) && !ty %in% prm$plot_types) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  # ragg when it is installed (several times faster, and needs no X11), otherwise the cairo
  # PNG device: the default device on an HPC node without X11 can take a minute or more
  # over one large figure
  if (requireNamespace("ragg", quietly = TRUE))
    ggplot2::ggsave(path, p, width = min(width, 49), height = min(height, 49),
                    dpi = dpi, bg = "white", limitsize = FALSE, device = ragg::agg_png)
  else
    ggplot2::ggsave(path, p, width = min(width, 49), height = min(height, 49),
                    dpi = dpi, bg = "white", limitsize = FALSE, device = "png", type = "cairo")
  invisible(path)
}

# Sample labels in figures: the full file stem, wrapped onto as many lines as it needs so a
# long name fits beside a row or in a panel strip. Lines break after ".", "_", "-" or a space
# (not inside a number such as 0.01uM) when one falls in the second half of the line, and
# are cut at `width` characters otherwise. Only line breaks are added, so every sample keeps
# a label as unique as its file stem, and the label is the same whichever samples are in
# the run (the short labels of v3.3.17 and earlier were not, and fell back to an index).
wrap_stem <- function(x, width = 30) {
  vapply(as.character(x), function(s) {
    if (is.na(s)) return(NA_character_)
    lines <- character(0)
    while (nchar(s) > width) {
      ch <- strsplit(substr(s, 1, width + 1), "")[[1]]
      ok <- which(ch[seq_len(width)] %in% c(".", "_", "-", " "))
      # a "." between two digits is part of a number
      ok <- ok[!(ch[ok] == "." & ok > 1 & grepl("[0-9]", ch[pmax(ok - 1, 1)]) & grepl("[0-9]", ch[pmin(ok + 1, length(ch))]))]
      # break after the last separator of a run (e.g. after "--", not between the dashes)
      ok <- ok[!(ok < length(ch) & ch[pmin(ok + 1, length(ch))] %in% c(".", "_", "-", " "))]
      cut <- if (length(ok) && max(ok) >= width / 2) max(ok) else width
      lines <- c(lines, substr(s, 1, cut)); s <- substring(s, cut + 1)
    }
    paste(c(lines, s), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

# Each sample's label for figures, named by file stem: lab <- sample_labels(stems); lab[stem].
# 45 characters a line suits the figures with one row per sample for the whole run (Modules 1,
# 3 and 5): a PacBio-style name takes two lines, so a run of 150 samples keeps a 7 pt font
sample_labels <- function(x, width = 45) { u <- unique(as.character(x)); stats::setNames(wrap_stem(u, width), u) }

# Height (inches) and label font size (pt) for a figure with one row per sample: each row
# is at least `per_row` inches and tall enough for its label's lines; when that would make
# the figure taller than `hi`, the font is made smaller (down to `min_size`) before the
# height is capped. `extra` is the height of everything but the rows (title, axis, legend).
label_rows <- function(lab, extra = 2, per_row = 0.14, size = 7, hi = 45, min_size = 4) {
  n <- length(unique(lab)); if (!n) return(list(height = extra + per_row, size = size))
  nl <- max(1, lengths(strsplit(as.character(lab), "\n", fixed = TRUE)))
  need <- function(sz) extra + n * max(per_row, nl * sz / 72 * 1.25)
  while (need(size) > hi && size > min_size) size <- size - 0.5
  list(height = min(need(size), hi), size = size)
}

# PacBio read quality (rq, predicted accuracy) as a Phred-scaled QV, capped at 60
rq_to_qv <- function(rq) pmin(-10 * log10(pmax(1 - rq, 1e-6)), 60)

# Label for a percentage axis that carries the % sign on the ticks
# (accuracy NULL lets scales choose enough decimals to keep the ticks distinct)
pct_axis <- function(accuracy = NULL) scales::label_percent(accuracy = accuracy)

# Consistent summary layers: box + every point + mean diamond, with the mean
# written in red beside the box (clear of the points). percent = TRUE writes a
# proportion as a percentage with `digits` decimals. jitter = 0 keeps points on the
# category line (needed when lines join paired points).
box_points_mean <- function(point_colour = "grey30", digits = 1, label = TRUE, percent = FALSE, jitter = 0.15, fmt = NULL) {
  lay <- list(
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "grey92", colour = "grey40", width = 0.5),
    ggplot2::geom_jitter(width = jitter, height = 0, size = 1.8, alpha = 0.7, colour = point_colour),
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 18, size = 4, colour = "red")
  )
  if (is.null(fmt)) fmt <- if (percent) function(v) scales::percent(v, accuracy = 10^-digits) else function(v) formatC(v, format = "f", digits = digits)
  if (label) lay <- c(lay, list(ggplot2::stat_summary(
    fun = mean, geom = "text", colour = "red", size = 3.2, hjust = 0, position = ggplot2::position_nudge(x = 0.3),
    ggplot2::aes(label = fmt(ggplot2::after_stat(y))))))
  lay
}

log_step <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# A local file path for a markdown image link: %, # and spaces percent-encoded, so pandoc finds
# the file however the sample or group is named
md_path <- function(path) {
  x <- normalizePath(path)
  x <- gsub("%", "%25", x, fixed = TRUE); x <- gsub("#", "%23", x, fixed = TRUE); gsub(" ", "%20", x, fixed = TRUE)
}

# A group name for display: leading and trailing dashes, underscores and spaces removed (a
# group made by joining settings columns, some blank, can end in "----"). The group itself is
# unchanged, since it is the key for group controls and ranges.
display_group <- function(x) gsub("^[-_ ]+|[-_ ]+$", "", as.character(x))


# ------------------------------------------------------------------------------
# Plot catalogue: every plot Duke saves, by type
# ------------------------------------------------------------------------------
# The type is the first part of the file name and, for types with more than one
# file, the folder they are saved in. `plots` selects types (default "all": every
# type with default = TRUE); `plot_units` selects the reads / molecules versions
# of the types that can show both (units = TRUE). scope: run (one figure for the
# run), sample (one per sample), range, group, measure.
plot_catalogue <- function() {
  x <- read.csv(text = "
type,module,scope,units,default,description
reads_per_sample,1,run,FALSE,TRUE,Reads per sample at each import step
read_length_distribution,1,run,FALSE,TRUE,Read length distribution pooled and per sample
quality_vs_length,1,run,FALSE,TRUE,Mean base quality against read length per sample
pacbio_np_rq,1,run,FALSE,TRUE,PacBio pass count and read quality per sample
pacbio_np_distribution,1,run,FALSE,TRUE,PacBio pass count of every read: share at each pass count with the cap marked; short and long reads compared
kinnex_segments,1,run,FALSE,TRUE,Kinnex reads: segments per concatenated read and pass count and read quality by the segment's start position
pacbio_np_vs_length,1,run,FALSE,TRUE,PacBio pass count against read length: median and middle 50% and 80% of reads per length band; share at the cap
adapter_trimming,1,run,FALSE,TRUE,Adapter trimming outcome
read_accounting,2,run,FALSE,TRUE,Reads kept at each alignment step
flank_boundary_shortfall,2,run,FALSE,TRUE,Flank alignments that stop short of the repeat
mapping_quality,2,run,FALSE,TRUE,Mapping quality per flank
strand_by_step,2,run,FALSE,TRUE,Share of forward-strand reads at each step
alignment,2,sample,FALSE,TRUE,Each read's flank alignments (one per sample and reference)
detection_rate,3,run,FALSE,TRUE,Reads with each tract found
flank_length,3,run,FALSE,TRUE,Flank length between primer and tract per sample
strand_by_length,3,run,FALSE,TRUE,Strand balance by tract length
count_method_agreement,3,run,FALSE,TRUE,How often the counting methods differ by tract length
count_method_gap,3,run,FALSE,TRUE,How far the counting methods differ by tract length
break_positions,3,run,FALSE,TRUE,Where broken tracts are broken
between_tract_sequences,3,run,FALSE,TRUE,Sequences between adjacent tracts
extraction_outcome,4,run,FALSE,TRUE,UMI extraction outcome per sample
umi_length,4,run,FALSE,TRUE,UMI length per sample
umi_neighbours,4,run,FALSE,TRUE,UMIs one mismatch from another against chance
length_rule_scorecard,4,run,FALSE,TRUE,Molecule length rules scored
length_rule_comparison,4,run,FALSE,TRUE,Molecule length rules compared
shortened_flag_by_length,4,run,FALSE,TRUE,Molecules flagged as shortened by length
shortened_flag_by_depth_band,4,run,FALSE,TRUE,Molecules flagged as shortened by length and depth
shortened_filter_distribution,4,run,FALSE,TRUE,Shortened-read share per molecule
shortened_flag_per_sample,4,run,FALSE,TRUE,Molecules flagged as shortened per sample
reads_per_molecule,4,run,FALSE,TRUE,Reads per molecule
molecules_per_sample,4,run,FALSE,TRUE,Molecules found and kept per sample
molecules_and_depth,4,run,FALSE,TRUE,Molecules against read depth per sample
within_molecule_error,4,run,FALSE,TRUE,Read length minus molecule length by molecule length
within_molecule_error_hist,4,run,FALSE,TRUE,Read length minus molecule length (histogram)
reads_far_from_molecule,4,run,FALSE,TRUE,Reads far from their molecule's length per sample
reads_vs_molecules_cohort,4,run,FALSE,TRUE,Read and molecule length distributions for the run
reads_vs_molecules,4,sample,FALSE,TRUE,Read and molecule length distributions (one per sample)
allele_shares,5,run,FALSE,TRUE,Each sample's reads and molecules split between its alleles
clusters_overview,5,run,FALSE,TRUE,Clusters in every sample
clusters,5,sample,FALSE,TRUE,One sample's lengths by cluster
variants_by_sample,5,run,FALSE,TRUE,Variants in alleles 1 and 2 of every sample
alleles,5,sample,FALSE,TRUE,One sample's alleles and variants
waterfall,5,sample,FALSE,TRUE,Reads stacked by length with each base coloured
waterfall_cluster,5,sample,FALSE,FALSE,Waterfall of each cluster of each sample
experiment_map,6,run,FALSE,TRUE,Samples by group and time
counts_per_range,6,run,TRUE,TRUE,Reads and molecules in each range
control_setpoints,6,run,TRUE,TRUE,Control samples and group setpoints
timecourse,6,range,TRUE,TRUE,Headline metrics against time (one per measure and range)
instability_index_by_range,6,run,TRUE,TRUE,Instability index in every range
pair_trend,6,range,TRUE,TRUE,Change in one tract against change in another (one per range)
pair_correlation,6,range,TRUE,TRUE,Per-sample correlation between two tracts (one per range)
pair_scatter,6,sample,TRUE,TRUE,One sample's lengths of two tracts per read or molecule (every range)
pair_map,6,measure,TRUE,TRUE,Two tracts' lengths against each other for every read or molecule (one per tract pair; a panel per group)
pair_map_group,6,group,TRUE,TRUE,One group's samples: two tracts' lengths against each other (one per group and tract pair)
heatmap,7,measure,TRUE,TRUE,Every sample's length distribution (one per measure)
group_density,7,group,TRUE,TRUE,A group's samples as ridges (all lengths and each range)
group_histogram,7,group,TRUE,TRUE,A group's samples as histograms (all lengths and each range)
sample_ranges,7,sample,TRUE,TRUE,One sample's length in every range
sample_card,7,sample,TRUE,TRUE,One sample's every measure (all lengths and the ranges drawn in the report)
sample_density,7,sample,TRUE,FALSE,One sample's ridge for all lengths and each range (one file each)
sample_histogram,7,sample,TRUE,FALSE,One sample's histogram for all lengths and each range (one file each)
", stringsAsFactors = FALSE, strip.white = TRUE)
  x$units <- as.logical(x$units); x$default <- as.logical(x$default)
  x
}

# The plot types this run saves, from `plots` ("all", or a list of types; "-type"
# removes a type, e.g. "all,-waterfall,sample_density") and the older switches.
resolve_plot_types <- function(p) {
  cat_ <- plot_catalogue()
  req <- p$plots %||% "all"
  if (!length(req)) req <- "all"
  add <- req[!startsWith(req, "-")]; drop <- sub("^-", "", req[startsWith(req, "-")])
  renamed <- c(clusters_per_sample = "allele_shares")   # plot types renamed since they were added
  for (old in intersect(names(renamed), c(add, drop))) {
    message("Plot type '", old, "' is now called '", renamed[[old]], "'")
    add[add == old] <- renamed[[old]]; drop[drop == old] <- renamed[[old]]
  }
  bad <- setdiff(c(add, drop), c("all", "none", cat_$type))
  if (length(bad)) stop("Unknown plot type(s) in plots: ", paste(bad, collapse = ", "), ". Known types: ", paste(cat_$type, collapse = ", "))
  out <- unique(c(if ("all" %in% add) cat_$type[cat_$default], setdiff(add, c("all", "none"))))
  # Older switches, so existing job scripts keep their meaning
  if (isFALSE(p$waterfall)) drop <- c(drop, "waterfall", "waterfall_cluster")
  if (isTRUE(p$waterfall_per_cluster)) out <- c(out, "waterfall_cluster")
  if (isFALSE(p$sample_cards)) drop <- c(drop, "sample_card")
  if (isTRUE(p$plot_sample_ranges)) out <- c(out, "sample_density", "sample_histogram")
  if (isFALSE(p$visualise_alignment)) drop <- c(drop, "alignment")
  if (isFALSE(p$plot_per_sample)) drop <- c(drop, cat_$type[cat_$scope == "sample"])
  setdiff(unique(out), drop)
}

want_plot <- function(type, p = NULL) {
  # knitr also evaluates chunk options once where params can be some package's function,
  # so anything but Duke's parameter list is ignored
  if (!is.list(p)) p <- .run_params(parent.frame())
  if (is.null(p)) return(TRUE)
  type %in% (p$plot_types %||% resolve_plot_types(p))
}

# The unit versions to draw for a plot that can show reads and molecules: any of
# "combined" (both in one figure), "reads", "molecules". Without molecules only
# reads can be drawn. can_combine = FALSE (e.g. one heatmap per unit) turns
# "combined" into one figure per unit.
plot_units_for <- function(p, have_molecules, can_combine = TRUE) {
  u <- p$plot_units %||% "combined"
  if (!have_molecules) return(if (length(u)) "reads" else character(0))
  if (!can_combine && "combined" %in% u) u <- c(setdiff(u, "combined"), "reads", "molecules")
  intersect(c("combined", "reads", "molecules"), u)
}
# Which units' data a version shows
unit_members <- function(u) if (identical(u, "combined")) c("reads", "molecules") else u

# File path for a saved plot: {type}_{part}_..._{unit}.png, in a folder named after
# the type when the type has a file per sample, range, group or measure
plot_path <- function(plots_dir, type, ..., unit = NULL, folder = NULL) {
  parts <- unlist(list(...)); parts <- parts[!is.na(parts) & nzchar(parts)]
  scope <- plot_catalogue()$scope[match(type, plot_catalogue()$type)]
  if (is.null(folder)) folder <- scope %in% c("sample", "range", "group", "measure")
  name <- paste(c(type, safe_path(parts), unit), collapse = "_")
  if (!isTRUE(folder)) return(file.path(plots_dir, paste0(name, ".png")))
  # A type with many files is split into subfolders, so each folder holds one set to
  # flick through: first by unit (combined / reads / molecules), then, for a type with
  # one file per sample or group in each range, by range (all-lengths, range01, ...), and
  # for a type with one file per range, by measure (or tract pair). The file name still
  # carries every part, so a copied file keeps its meaning.
  sub <- if (!is.null(unit)) unit else character(0)
  if (isTRUE(scope %in% c("sample", "group"))) {
    rg <- parts[grepl("^(range[0-9]+|all-lengths)$", parts)]
    if (length(rg)) sub <- c(sub, rg[1])
  }
  if (isTRUE(scope == "range") && length(parts)) sub <- c(sub, safe_path(parts[1]))
  do.call(file.path, as.list(c(plots_dir, type, sub, paste0(name, ".png"))))
}

# The plot type a file name starts with (the longest catalogue type that matches), or NA
.plot_types_cache <- new.env()
plot_type_of <- function(path) {
  if (is.null(.plot_types_cache$t)) { t <- plot_catalogue()$type; .plot_types_cache$t <- t[order(-nchar(t))] }
  t <- .plot_types_cache$t
  vapply(sub("\\.png$", "", basename(path)), function(b) { m <- t[b == t | startsWith(b, paste0(t, "_"))]; if (length(m)) m[1] else NA_character_ }, "", USE.NAMES = FALSE)
}

# ------------------------------------------------------------------------------
# Kinnex (skera) segments
# ------------------------------------------------------------------------------
# skera names each segment <movie>/<zmw>/ccs[/fwd|rev]/<start>_<end>: where in bp the segment
# sits in the concatenated HiFi read it was cut from. Returns, for every name, the concatenated
# read it belongs to (read_id = movie/zmw) and the segment's start and end, or NULL when most
# names are not Kinnex segments. Names are split once on fixed characters (fast on millions).
kinnex_positions <- function(qname) {
  if (!length(qname)) return(NULL)
  ok <- grepl("/ccs(/[a-z]+)?/[0-9]+_[0-9]+$", qname)
  if (mean(ok) < 0.5) return(NULL)
  parts <- data.table::tstrsplit(qname, "/", fixed = TRUE)
  lastf <- if (length(parts) > 4) sub(".*/", "", qname) else parts[[length(parts)]]
  se <- data.table::tstrsplit(lastf, "_", fixed = TRUE)
  d <- data.table::data.table(read_id = paste(parts[[1]], parts[[2]], sep = "/"),
                              seg_start = suppressWarnings(as.numeric(se[[1]])), seg_end = suppressWarnings(as.numeric(se[[2]])))
  d[!ok, `:=`(read_id = NA_character_, seg_start = NA_real_, seg_end = NA_real_)]
  d[]
}

