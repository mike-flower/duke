# ==============================================================================
# Duke - general utilities
# ==============================================================================
# Excel export, threading, timing, safe cache loading, path sanitising,
# parameter hashing, house-style tables and plot saving.
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

# ------------------------------------------------------------------------------
# Excel export
# ------------------------------------------------------------------------------
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
  invisible(path)
}

# ------------------------------------------------------------------------------
# Threading wrapper. Uses params$threads from the calling environment.
# ------------------------------------------------------------------------------
apply_fn <- function(X, FUN, ..., threads = NULL) {
  threads <- threads %||% tryCatch(get("params", envir = parent.frame())$threads,
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
  x <- gsub("[/\\\\:*?\"<>|[:cntrl:]]", replacement, x, perl = TRUE)
  x <- gsub("[[:space:]]+", replacement, x)
  x <- gsub(paste0("^", replacement, "+|", replacement, "+$"), "", x)
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
  dpi <- dpi %||% tryCatch(get("params", envir = parent.frame())$plot_dpi, error = function(e) 150)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(path, p, width = min(width, 49), height = min(height, 49),
                  dpi = dpi, bg = "white", limitsize = FALSE)
  invisible(path)
}

# Shorter sample labels for plot facets: drop the prefix every name shares (cut back
# to the last '.', '_' or '-' so a label never starts mid-word), e.g.
# "m84277_..._s1.hifi_reads.bc2083.skera.demux.kxF01--kxR13" -> "kxF01--kxR13"
short_labels <- function(x) {
  u <- unique(as.character(x))
  if (length(u) < 2) return(setNames(u, u))
  ch <- strsplit(u, "")
  n <- min(lengths(ch)); k <- 0
  while (k < n && length(unique(vapply(ch, `[`, "", k + 1))) == 1) k <- k + 1
  pre <- substr(u[1], 1, k)
  cut <- max(c(0, gregexpr("[._-]", pre)[[1]]))
  out <- substring(u, cut + 1)
  if (any(!nzchar(out)) || anyDuplicated(out)) out <- u
  setNames(out, u)
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
