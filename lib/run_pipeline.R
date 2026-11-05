# ==============================================================================
# Duke - pipeline runner (used by scripts/duke_run.R and scripts/duke_cli.R)
# ==============================================================================

duke_version <- "3.3.13"

.module_files <- c(
  `1` = "01_import_qc.Rmd", `2` = "02_alignment.Rmd", `3` = "03_repeat_detection.Rmd",
  `4` = "04_umi.Rmd", `5` = "05_allele_structure.Rmd", `6` = "06_instability.Rmd",
  `7` = "07_repeat_lengths.Rmd")
.module_data <- c(
  `1` = "01_import_qc.RData", `2` = "02_alignment.RData", `3` = "03_repeat_detection.RData",
  `4` = "04_umi.RData", `5` = "05_allele_structure.RData", `6` = "06_instability.RData",
  `7` = "07_repeat_lengths.RData")
.module_names <- c(`1` = "Import and QC", `2` = "Alignment", `3` = "Repeat detection", `4` = "UMI",
                   `5` = "Allele structure", `6` = "Instability", `7` = "Repeat lengths")

# Parameters that change each module's results. A module is re-run when these,
# or anything upstream of it, change.
.module_keys <- list(
  `1` = c("dir_data", "data_signature", "import_patterns", "import_recursive", "r1_pattern", "r2_pattern", "select_one_of_pair",
          "downsample", "np_min", "check_duplicate_readnames", "trim", "path_trim_patterns",
          "trim_max_mismatch", "trim_with_indels"),
  `2` = c("minimap2_args", "path_manual_exclusions", "ref_signature", "keep_bam", "visualise_alignment_downsample"),
  `3` = c("tracts_table", "composites_table", "rm_flank_length_outliers", "flank_iqr_multiplier", "export_read_counts"),
  `4` = c("umi", "umi_flank_5", "umi_flank_3", "umi_flank_max_mismatch", "umi_search_window", "umi_len_min",
          "umi_len_max", "umi_cluster_max_dist", "umi_min_reads", "umi_collapse_stat", "umi_top_cluster_k",
          "umi_top_cluster_w", "umi_top_cluster_wfrac", "umi_top_cluster_min_share", "umi_shortened_filter", "umi_shortened_min_units", "umi_shortened_min_share"),
  # (not settings_signature: Module 4 uses the settings file only to group one figure, and
  # adding or editing the settings must not re-run UMI collapsing and everything after it)
  `5` = c("cluster", "cluster_by", "cluster_unit", "repeat_cluster_max", "repeat_cluster_min_prop", "repeat_cluster_downsample",
          "haplotype_cluster_max", "haplotype_cluster_downsample", "haplotype_region", "haplotype_trim_length",
          "consensus_threshold", "consensus_downsample", "consensus_parts", "call_variants", "allele_minor_share", "allele_support_reads",
          "allele_bam", "allele_bam_reads", "allele_bam_method", "allele_bam_minimap2", "waterfall", "waterfall_downsample",
          "waterfall_per_cluster", "waterfall_flank_bp"),
  `6` = c("settings_signature", "range_peak_span", "range_min_n", "group_control", "control_sample_selection",
          "control_setpoint_metric", "control_aggregation_method", "expansion_thresholds", "cohort_plot_by", "time_label",
          "headline_metrics", "headline_setpoint_type", "plot_ranges", "measure_pairs", "pair_max_points"),
  `7` = c("sample_cards", "repeat_histogram_binwidth", "settings_signature", "plot_every_range", "plot_sample_ranges", "plot_dots_max", "plot_dots_unit"))
.module_upstream <- list(`1` = character(0), `2` = "1", `3` = "2", `4` = "3", `5` = c("3", "4"),
                         `6` = c("3", "4"), `7` = c("3", "4", "6"))

# Names, sizes and modification times of the input files, so that adding or
# replacing a file in dir_data re-runs Module 1
.dir_signature <- function(dir, recursive = TRUE) {
  if (is.na(dir) || !dir.exists(dir)) return(NA_character_)
  f <- list.files(dir, recursive = recursive, full.names = TRUE)
  if (!length(f)) return(NA_character_)
  i <- file.info(f)
  digest::digest(paste(basename(f), i$size, as.numeric(i$mtime)), algo = "md5")
}

.file_signature <- function(paths) {
  paths <- paths[!is.na(paths) & file.exists(paths)]
  if (!length(paths)) return(NA_character_)
  paste(basename(paths), unname(tools::md5sum(paths)), collapse = ";")
}

module_hashes <- function(p) {
  h <- list()
  for (m in names(.module_keys)) {
    up <- unlist(h[.module_upstream[[m]]])
    h[[m]] <- substr(digest::digest(list(up, p[intersect(.module_keys[[m]], names(p))]), algo = "md5"), 1, 10)
  }
  # Which plots a module saves changes only that module: a suffix for its plot choice is
  # added to its hash (not passed downstream), and only when the choice differs from the
  # default, so results made before the plots setting existed are still resumed
  cat_ <- plot_catalogue()
  def <- resolve_plot_types(list(plots = "all"))
  for (m in names(h)) {
    ty <- cat_$type[cat_$module == as.integer(m)]
    mine <- sort(intersect(if (is.null(p$plot_types)) def else p$plot_types, ty)); dflt <- sort(intersect(def, ty))
    units_changed <- any(cat_$units[cat_$module == as.integer(m)]) && !identical(sort(p$plot_units %||% "combined"), "combined")
    if (!identical(mine, dflt) || units_changed)
      h[[m]] <- paste0(h[[m]], "-p", substr(digest::digest(list(mine, if (units_changed) sort(p$plot_units)), algo = "md5"), 1, 6))
  }
  h
}

# ------------------------------------------------------------------------------
# Reusing modules from an earlier run (reuse_from, reuse_modules)
# ------------------------------------------------------------------------------
# Copies a module's results from an earlier output folder into this run's, so the
# module is skipped instead of re-run: its saved data (module_data/), its output
# folder, its report and its working caches (temp/). A module is reused only when
# its data were made with the same settings as this run (the same data part of the
# module hash, which takes in the settings of every module upstream), so a reused
# module is the same as a re-run one would be, apart from anything the newer Duke
# version changed in that module's own figures and reports. The plot choice (plots,
# plot_units) is not part of that check: a reused module keeps the earlier run's
# figures. reuse_modules names the modules to take (default: every module whose
# data still match); a named module that cannot be reused stops the run before
# anything is done, so a long job does not start on a wrong assumption.
.data_hash <- function(h) sub("-p.*$", "", h)
.hash_version <- function(path) {
  side <- paste0(path, ".hash")
  if (!file.exists(side)) return(NA_character_)
  l <- readLines(side, warn = FALSE)
  v <- grep("^duke ", l, value = TRUE)
  v <- if (length(v)) sub("^duke ", "", v[1]) else NA_character_
  if (identical(v, "unknown")) NA_character_ else v
}
# The version that last drew a module's report and figures, when that was a redraw
.hash_drawn <- function(path) {
  side <- paste0(path, ".hash")
  if (!file.exists(side)) return(NA_character_)
  v <- grep("^redrawn by duke ", readLines(side, warn = FALSE), value = TRUE)
  if (length(v)) sub("^redrawn by duke ", "", v[1]) else NA_character_
}
# A module's saved data are linked rather than copied when both folders are on the same file
# system: the file is not duplicated, and the copy takes seconds instead of minutes. Duke never
# writes into a saved data file (a module that runs again deletes its file first and writes a
# new one), so the earlier run's file cannot be changed through the link.
.link_or_copy <- function(from, to) {
  unlink(to)
  ok <- suppressWarnings(file.link(from, to))
  if (!isTRUE(ok)) ok <- file.copy(from, to, overwrite = TRUE, copy.date = TRUE)
  isTRUE(ok)
}
reuse_earlier_run <- function(p, hashes, results_dir) {
  src <- normalizePath(p$reuse_from, mustWork = FALSE)
  if (!dir.exists(src)) stop("reuse_from: folder not found: ", src)
  if (identical(src, normalizePath(p$dir_out, mustWork = FALSE))) stop("reuse_from is the output folder itself; give the earlier run's folder")
  named <- length(p$reuse_modules) > 0
  want <- as.character(if (named) p$reuse_modules else sort(unique(p$run_modules)))
  bad <- setdiff(want, names(.module_data)); if (length(bad)) stop("reuse_modules: no module ", paste(bad, collapse = ", "))
  # Modules this run does not run at all have nothing to reuse
  off <- c(if (!isTRUE(p$umi)) "4", if (!isTRUE(p$cluster) && !isTRUE(p$waterfall)) "5", if (is.na(p$path_settings) || !file.exists(p$path_settings)) "6")
  want <- setdiff(want, off)
  cat("\nReusing from: ", src, "\n", sep = "")
  plan <- list()
  for (m in want) {
    from <- file.path(src, "module_data", .module_data[[m]])
    here <- file.path(results_dir, .module_data[[m]])
    why <- NULL
    if (identical(.data_hash(stored_hash(here)), .data_hash(hashes[[m]]))) { cat(sprintf("  Module %s: already here with the current settings\n", m)); next }
    if (!file.exists(from)) why <- "not in that folder"
    else if (!identical(.data_hash(stored_hash(from)), .data_hash(hashes[[m]])))
      why <- "made with different settings (this module's or an upstream module's), so it would not match this run"
    else if (!named && !identical(.hash_version(from), duke_version))
      why <- sprintf("made by Duke %s; its settings match, but this version may have changed its outputs. Name it in reuse_modules to reuse it anyway",
                     if (is.na(.hash_version(from))) "before 3.3.2" else .hash_version(from))
    if (!is.null(why)) {
      if (named) {
        msg <- sprintf("reuse_modules: Module %s (%s) cannot be reused from %s: %s. Leave it out of reuse_modules to run it. Nothing has been copied or run.", m, .module_names[[m]], src, why)
        cat("\nSTOPPED:", msg, "\n")
        stop(msg, call. = FALSE)
      }
      cat(sprintf("  Module %s: not reused (%s); it will run\n", m, why)); next
    }
    plan[[m]] <- from
  }
  for (m in names(plan)) {
    stem <- tools::file_path_sans_ext(.module_files[[m]])
    unlink(c(file.path(p$dir_out, stem), file.path(p$dir_out, paste0(stem, ".html"))), recursive = TRUE)
    if (dir.exists(file.path(src, stem))) file.copy(file.path(src, stem), p$dir_out, recursive = TRUE, copy.date = TRUE)
    if (file.exists(file.path(src, paste0(stem, ".html")))) file.copy(file.path(src, paste0(stem, ".html")), p$dir_out, overwrite = TRUE, copy.date = TRUE)
    tmp <- list.files(file.path(src, "temp"), pattern = paste0("^m", m, "_", .data_hash(hashes[[m]])), full.names = TRUE)
    if (length(tmp)) { dir.create(file.path(p$dir_out, "temp"), showWarnings = FALSE); file.copy(tmp, file.path(p$dir_out, "temp"), recursive = TRUE, copy.date = TRUE) }
    if (!.link_or_copy(plan[[m]], file.path(results_dir, .module_data[[m]]))) stop("reuse_from: could not copy ", plan[[m]])
    v <- .hash_version(plan[[m]]); dv <- .hash_drawn(plan[[m]])
    # The current hash, so the module is skipped; the version lines keep the version that made
    # the results and, after a redraw, the version that drew the report and figures
    writeLines(c(hashes[[m]], paste("duke", if (is.na(v)) "unknown" else v), if (!is.na(dv)) paste("redrawn by duke", dv), paste("reused from", src)),
               paste0(file.path(results_dir, .module_data[[m]]), ".hash"))
    vtxt <- if (is.na(v)) "before 3.3.2" else v
    drawn <- if (is.na(dv)) v else dv
    note <- if (m %in% as.character(p$redraw)) c(if (is.na(v) || v != duke_version) sprintf("made by Duke %s", vtxt), "its report and figures will be redrawn") else
      c(if (!is.na(dv)) sprintf("results made by Duke %s, report and figures drawn by Duke %s", vtxt, dv)
        else if (is.na(v) || v != duke_version) sprintf("made by Duke %s, so its report and figures are that version's", vtxt),
        if (!identical(stored_hash(plan[[m]]), hashes[[m]])) "its figures follow the earlier run's plots / plot_units")
    cat(sprintf("  Module %s (%s): reused%s\n", m, .module_names[[m]], if (length(note)) paste0("; ", paste(note, collapse = "; ")) else ""))
  }
  invisible(names(plan))
}

# The hash a module's results were made with. A small sidecar file
# (<module>.RData.hash, written after each successful run) avoids loading
# large .RData files just to read it; it is only trusted when it is newer
# than the .RData.
stored_hash <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  side <- paste0(path, ".hash")
  if (file.exists(side) && file.mtime(side) >= file.mtime(path)) return(readLines(side, n = 1, warn = FALSE))
  e <- new.env()
  ok <- tryCatch({ load(path, envir = e); TRUE }, error = function(err) FALSE)
  if (!ok) return(NA_character_)
  obj <- mget(ls(e), envir = e)[[1]]
  obj$hash %||% NA_character_
}

# Warnings raised while a module renders (outside what its report shows) are written to
# the log as they happen, each different message once with the module it came from, and
# counted at the end of the run. Without this R prints only "There were N warnings" when
# the run ends, with no message and no module.
.warn_log <- new.env()
with_warning_log <- function(where, expr) {
  withCallingHandlers(expr, warning = function(w) {
    msg <- gsub("\\s+", " ", trimws(conditionMessage(w)))
    cl <- conditionCall(w); cl <- if (is.null(cl)) "" else paste(deparse(cl, nlines = 1L), collapse = "")
    if (startsWith(cl, "withCallingHandlers(")) cl <- ""
    key <- paste(where, msg, sep = "\r")
    n <- .warn_log$counts[[key]] %||% 0L
    if (n == 0L) {
      cat(sprintf("  Warning (%s): %s%s\n", where, msg, if (nzchar(cl)) sprintf("  [in %s]", substr(cl, 1, 120)) else ""))
      .warn_log$first <- rbind(.warn_log$first, data.frame(where = where, message = msg, call = cl, stringsAsFactors = FALSE))
    }
    .warn_log$counts[[key]] <- n + 1L
    invokeRestart("muffleWarning")
  })
}

run_duke <- function(user_params, duke_root = getwd()) {
  .warn_log$counts <- list(); .warn_log$first <- NULL
  duke_root <- normalizePath(duke_root)
  setwd(duke_root)
  # Reports use en dashes and other non-ASCII text; under a non-UTF-8 locale (e.g. "C" on
  # some nodes) R writes them as <U+2013>
  if (!isTRUE(l10n_info()$`UTF-8`))
    for (loc in c("C.UTF-8", "en_GB.UTF-8", "en_US.UTF-8")) if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
  p <- resolve_params(user_params)
  # A settings path that was given but does not exist is almost always a typo: stop now
  # rather than run every other module and silently skip Module 6
  if (!is.na(p$path_settings) && !file.exists(p$path_settings))
    stop("Settings file not found: ", p$path_settings, "\n  Check the path, or leave path_settings unset to run without Module 6.")
  dir.create(p$dir_out, showWarnings = FALSE, recursive = TRUE)

  p$tracts_table <- resolve_tracts(p)
  p$composites_table <- resolve_composites(p, p$tracts_table)
  p$ref_signature <- .file_signature(unique(p$tracts_table$path_ref))
  p$settings_signature <- .file_signature(p$path_settings)
  p$data_signature <- .dir_signature(p$dir_data, p$import_recursive)
  p$duke_version <- duke_version
  p$duke_root <- duke_root
  hashes <- module_hashes(p)

  # Logs go in logs/ at the Duke root, one folder per run
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  log_dir <- file.path(duke_root, "logs", stamp)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(log_dir, paste0("duke_run-", stamp, ".log"))
  params_file <- file.path(log_dir, paste0("params-", stamp, ".R"))
  writeLines(c(paste("# Duke", duke_version, "run parameters,", format(Sys.time())), "params <- "), params_file)
  utils::capture.output(dput(p[setdiff(names(p), c("tracts_table", "composites_table"))]), file = params_file, append = TRUE)
  utils::write.csv(p$tracts_table, file.path(log_dir, "tracts.csv"), row.names = FALSE)
  log_con <- file(log_file, open = "wt")
  sink(log_con, type = "output", split = TRUE); sink(log_con, type = "message")
  on.exit({ sink(type = "message"); sink(type = "output"); close(log_con) }, add = TRUE)

  start <- Sys.time()
  cat("\n=================================================================\n")
  cat(" DUKE", duke_version, "\n")
  cat("=================================================================\n")
  cat("Started:  ", format(start), "\nData:     ", p$dir_data, "\nOutput:   ", p$dir_out,
      "\nThreads:  ", p$threads, "\nModules:  ", paste(p$run_modules, collapse = ", "), "\nLog:      ", log_file, "\n")
  cat("Tracts:   ", paste0(p$tracts_table$tract, " (", p$tracts_table$motif, ", ", basename(p$tracts_table$path_ref), ")", collapse = "; "), "\n")
  if (nrow(p$composites_table)) cat("Combined: ", paste0(p$composites_table$name, " = ", p$composites_table$tracts, " (", p$composites_table$mode, ")", collapse = "; "), "\n")

  results_dir <- file.path(p$dir_out, "module_data")
  dir.create(results_dir, showWarnings = FALSE)
  reused <- if (!is.na(p$reuse_from %||% NA)) with_warning_log("reuse", reuse_earlier_run(p, hashes, results_dir)) else character(0)
  # redraw needs saved results made with the current settings: checked before any module
  # runs, so a long job does not stop hours in
  for (m in setdiff(as.character(p$redraw), "7")) {
    if (m == "4" && !isTRUE(p$umi)) next
    out_path <- file.path(results_dir, .module_data[[m]])
    if (!identical(.data_hash(stored_hash(out_path)), .data_hash(hashes[[m]]))) {
      msg <- paste0("redraw: Module ", m, " (", .module_names[[m]], ") has no saved results made with the current settings in ", results_dir,
                    if (file.exists(out_path)) " (they were made with different settings)" else "",
                    ". Leave it out of redraw to re-run it, or take its results from an earlier run with reuse_from and reuse_modules. Nothing has been run.")
      cat("\nSTOPPED:", msg, "\n"); stop(msg, call. = FALSE)
    }
  }
  # Modules whose results in dir_out belong to this run (run now, resumed or reused); the
  # run summary and index use only these, not files an earlier run left behind
  current <- character(0)
  for (m in as.character(sort(unique(p$run_modules)))) {
    cat("\n-----------------------------------------------------------------\n")
    cat("Module ", m, ": ", .module_names[[m]], "\n", sep = "")
    redraw_now <- m %in% as.character(p$redraw)
    if (m %in% reused && !redraw_now) { cat("Reused from", p$reuse_from, "\n"); current <- c(current, m); next }
    if (m == "4" && !isTRUE(p$umi)) { cat("Skipped (umi = FALSE)\n"); next }
    if (m == "5" && !isTRUE(p$cluster) && !isTRUE(p$waterfall)) { cat("Skipped (cluster = FALSE and waterfall = FALSE)\n"); next }
    if (m == "6" && (is.na(p$path_settings) || !file.exists(p$path_settings))) {
      cat("Skipped (no settings file - set path_settings to calculate instability metrics)\n"); next }
    out_path <- file.path(results_dir, .module_data[[m]])
    # redraw: the report and figures are drawn again from the saved results, which are
    # kept as they are (Module 7 only draws, so it simply re-runs)
    in_place <- redraw_now && m != "7"
    if (in_place) {
      cat(if (m %in% reused) "Reused from " else "", if (m %in% reused) p$reuse_from else "", if (m %in% reused) "; " else "",
          "redrawing the report and figures from the saved results\n", sep = "")
    } else if (redraw_now) {
      cat("Redrawn by re-running (this module only draws)\n")
    } else if (isTRUE(p$resume) && identical(stored_hash(out_path), hashes[[m]])) {
      cat("Skipped - results found with identical parameters\n"); current <- c(current, m); next
    }
    if (!redraw_now && file.exists(out_path)) cat("Parameters or upstream results changed - re-running\n")
    # Every upstream module must hold results made with the current parameters,
    # otherwise this module would be built on stale data and later resumed as current
    for (u in .module_upstream[[m]]) {
      if (u == "4" && !isTRUE(p$umi)) next
      if (u == "6" && (is.na(p$path_settings) || !file.exists(p$path_settings))) next
      # (the data part only: a module whose plot choice changed still holds current data)
      if (!identical(sub("-p.*$", "", stored_hash(file.path(results_dir, .module_data[[u]]))), sub("-p.*$", "", hashes[[u]])))
        stop("Module ", m, " needs Module ", u, " results made with the current parameters, but they are missing or out of date. ",
             "Add ", u, " to run_modules (resume will reuse whatever is still current).")
    }
    p$module_hash <- hashes[[m]]
    p$redraw_module <- in_place
    t0 <- Sys.time()
    # Clear the previous result first: a module saves its data before its report is
    # finished, so if this run stopped part-way the old hash must not let a resume skip it
    # (a redraw keeps the saved results, and they are still current if it stops part-way)
    data_version <- .hash_version(out_path)
    if (!in_place) unlink(c(out_path, paste0(out_path, ".hash")))
    # Each run renders in its own working folder. By default rmarkdown writes its working
    # files (the .knit.md and the figure folder) next to the Rmd, i.e. in the shared repo,
    # so two runs rendering the same module at once could delete each other's figures.
    work_dir <- file.path(p$dir_out, "temp", paste0("render_", m))
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    with_warning_log(paste("Module", m), rmarkdown::render(file.path(duke_root, "modules", .module_files[[m]]), output_dir = p$dir_out,
                      intermediates_dir = work_dir,
                      envir = list2env(list(params = p), parent = globalenv()),
                      knit_root_dir = duke_root, quiet = TRUE))
    unlink(c(work_dir, file.path(p$dir_out, paste0(tools::file_path_sans_ext(.module_files[[m]]), "_files"))), recursive = TRUE)
    # After a redraw the version line keeps the version that made the results; the hash
    # takes the current plot choice, as the figures now follow it
    if (in_place) writeLines(c(hashes[[m]], paste("duke", if (is.na(data_version)) "unknown" else data_version), paste("redrawn by duke", duke_version)), paste0(out_path, ".hash"))
    else if (file.exists(out_path)) writeLines(c(hashes[[m]], paste("duke", duke_version)), paste0(out_path, ".hash"))
    current <- c(current, m)
    cat("Done in", format_elapsed(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "\n")
  }

  # Modules already current before this run's list (e.g. run_modules = 6:7 after a full run)
  for (m in setdiff(names(.module_data), current))
    if (identical(.data_hash(stored_hash(file.path(results_dir, .module_data[[m]]))), .data_hash(hashes[[m]]))) current <- c(current, m)
  p$current_modules <- current
  rs <- tryCatch(with_warning_log("run summary", write_run_summary(p)), error = function(e) { cat("\nRun summary not written:", conditionMessage(e), "\n"); NULL })
  if (!is.null(rs)) cat("\nRun summary: run_summary.xlsx (one row per sample, from every module)\n")
  ix <- tryCatch(with_warning_log("index", write_run_index(p, rs)), error = function(e) { cat("\n00_index.html not written:", conditionMessage(e), "\n"); NULL })
  if (!is.null(ix)) cat("Index:       00_index.html (links every report, table and figure)\n")

  if (isTRUE(p$remove_temp)) {
    td <- file.path(p$dir_out, "temp")
    if (dir.exists(td)) { unlink(td, recursive = TRUE); cat("\nRemoved temp/\n") }
  }
  if (length(.warn_log$counts)) {
    w <- .warn_log$first
    w$times <- unlist(.warn_log$counts[paste(w$where, w$message, sep = "\r")])
    utils::write.table(w, file.path(log_dir, paste0("warnings-", stamp, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("\nWarnings: %d (%d different), each shown above where it happened; listed in %s\n", sum(w$times), nrow(w), file.path(log_dir, paste0("warnings-", stamp, ".tsv"))))
    for (i in seq_len(nrow(w))) cat(sprintf("  %s, %d x: %s\n", w$where[i], w$times[i], substr(w$message[i], 1, 200)))
  } else cat("\nWarnings: none\n")
  cat("\n=================================================================\n")
  cat("Finished:", format(Sys.time()), " Total:", format_elapsed(as.numeric(difftime(Sys.time(), start, units = "secs"))), "\n")
  cat("Reports:  ", paste(list.files(p$dir_out, pattern = "\\.html$"), collapse = ", "), "\n")
  invisible(p)
}
