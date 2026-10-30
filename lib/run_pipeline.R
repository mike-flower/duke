# ==============================================================================
# Duke - pipeline runner (used by scripts/duke_run.R and scripts/duke_cli.R)
# ==============================================================================

duke_version <- "3.0.10"

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
  `2` = c("minimap2_args", "path_manual_exclusions", "ref_signature"),
  `3` = c("tracts_table", "composites_table", "rm_flank_length_outliers", "flank_iqr_multiplier"),
  `4` = c("umi", "umi_flank_5", "umi_flank_3", "umi_flank_max_mismatch", "umi_search_window", "umi_len_min",
          "umi_len_max", "umi_cluster_max_dist", "umi_min_reads", "umi_collapse_stat", "umi_top_cluster_k",
          "umi_top_cluster_w", "umi_top_cluster_wfrac", "umi_top_cluster_min_share", "umi_long_read_filter", "umi_long_read_units", "umi_long_read_share"),
  `5` = c("cluster", "cluster_by", "repeat_cluster_max", "repeat_cluster_min_prop", "repeat_cluster_downsample",
          "haplotype_cluster_max", "haplotype_cluster_downsample", "haplotype_region", "haplotype_trim_length",
          "consensus_threshold", "consensus_downsample", "call_variants", "waterfall", "waterfall_downsample",
          "waterfall_per_cluster", "waterfall_flank_bp"),
  `6` = c("settings_signature", "range_peak_span", "group_control", "control_sample_selection",
          "control_setpoint_metric", "control_aggregation_method", "use_manual_controls", "expansion_thresholds", "cohort_plot_by",
          "headline_metrics", "headline_setpoint_type", "plot_ranges", "measure_pairs", "pair_max_points"),
  `7` = c("sample_cards", "repeat_histogram_binwidth", "settings_signature"))
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
  h
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

run_duke <- function(user_params, duke_root = getwd()) {
  duke_root <- normalizePath(duke_root)
  setwd(duke_root)
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
  for (m in as.character(sort(unique(p$run_modules)))) {
    cat("\n-----------------------------------------------------------------\n")
    cat("Module ", m, ": ", .module_names[[m]], "\n", sep = "")
    if (m == "4" && !isTRUE(p$umi)) { cat("Skipped (umi = FALSE)\n"); next }
    if (m == "5" && !isTRUE(p$cluster) && !isTRUE(p$waterfall)) { cat("Skipped (cluster = FALSE and waterfall = FALSE)\n"); next }
    if (m == "6" && (is.na(p$path_settings) || !file.exists(p$path_settings))) {
      cat("Skipped (no settings file - set path_settings to calculate instability metrics)\n"); next }
    out_path <- file.path(results_dir, .module_data[[m]])
    if (isTRUE(p$resume) && identical(stored_hash(out_path), hashes[[m]])) {
      cat("Skipped - results found with identical parameters\n"); next
    }
    if (file.exists(out_path)) cat("Parameters or upstream results changed - re-running\n")
    # Every upstream module must hold results made with the current parameters,
    # otherwise this module would be built on stale data and later resumed as current
    for (u in .module_upstream[[m]]) {
      if (u == "4" && !isTRUE(p$umi)) next
      if (u == "6" && (is.na(p$path_settings) || !file.exists(p$path_settings))) next
      if (!identical(stored_hash(file.path(results_dir, .module_data[[u]])), hashes[[u]]))
        stop("Module ", m, " needs Module ", u, " results made with the current parameters, but they are missing or out of date. ",
             "Add ", u, " to run_modules (resume will reuse whatever is still current).")
    }
    p$module_hash <- hashes[[m]]
    t0 <- Sys.time()
    rmarkdown::render(file.path(duke_root, "modules", .module_files[[m]]), output_dir = p$dir_out,
                      envir = list2env(list(params = p), parent = globalenv()),
                      knit_root_dir = duke_root, quiet = TRUE)
    if (file.exists(out_path)) writeLines(hashes[[m]], paste0(out_path, ".hash"))
    cat("Done in", format_elapsed(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "\n")
  }

  if (isTRUE(p$remove_temp)) {
    td <- file.path(p$dir_out, "temp")
    if (dir.exists(td)) { unlink(td, recursive = TRUE); cat("\nRemoved temp/\n") }
  }
  cat("\n=================================================================\n")
  cat("Finished:", format(Sys.time()), " Total:", format_elapsed(as.numeric(difftime(Sys.time(), start, units = "secs"))), "\n")
  cat("Reports:  ", paste(list.files(p$dir_out, pattern = "\\.html$"), collapse = ", "), "\n")
  invisible(p)
}
