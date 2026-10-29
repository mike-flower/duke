#!/usr/bin/env Rscript
# ==============================================================================
# Duke - command-line interface
# ==============================================================================
# Usually called through the ./duke wrapper:
#   ./duke --dir_data DIR --dir_out DIR --path_ref REF.fasta [options]
#   ./duke --help                 every option, grouped by module, with defaults
#   ./duke ... --dry_run          print the resolved parameters and tracts, then stop
#
# Options are generated from lib/params.R, so every parameter in duke_run.R is
# available here. Vectors are comma-separated (--run_modules 3,6,7;
# --composites "CAG_CCG_sum=CAG+CCG:sum,CAG_CCG_span=CAG+CCG:span").
# Logical options take TRUE / FALSE; a bare flag (--umi) means TRUE.
# ==============================================================================

duke_root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(dirname(normalizePath(f))) else getwd()
})
setwd(duke_root)
source("lib/load_all.R")

spec <- duke_param_spec()
names(spec) <- vapply(spec, `[[`, "", "name")

show_help <- function() {
  cat("Duke", duke_version, "- repeat-length and somatic instability pipeline\n\n")
  cat("Usage: ./duke --dir_data DIR --dir_out DIR --path_ref REF.fasta [options]\n\n")
  cat("  --help        Show this help\n  --dry_run     Resolve and print the parameters and tracts, then stop\n")
  show_default <- function(d) {
    if (is.null(d) || length(d) == 0) return("(none)")
    if (length(d) == 1 && is.na(d)) return("NA")
    if (is.numeric(d) && length(d) > 2 && all(diff(d) == 1)) return(paste0(min(d), ":", max(d)))
    paste(d, collapse = ",")
  }
  grp <- ""
  for (s in spec) {
    if (s$group != grp) { cat("\n", s$group, "\n", sep = ""); grp <- s$group }
    cat(sprintf("  --%-30s %s\n", s$name, paste0("[", s$type, "] default: ", show_default(s$default))))
    cat(paste0(strwrap(s$help, 76, indent = 6, exdent = 6), collapse = "\n"), "\n")
  }
  cat("\nSee README.md for the tracts table, composites, UMI design and metric definitions.\n")
}

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || any(args %in% c("--help", "-h"))) { show_help(); quit(status = 0) }

user <- list(); dry_run <- FALSE; i <- 1
while (i <= length(args)) {
  a <- args[i]
  if (!startsWith(a, "--")) stop("Unexpected argument: ", a, " (options start with --)")
  kv <- sub("^--", "", a)
  if (grepl("=", kv)) { key <- sub("=.*", "", kv); val <- sub("^[^=]*=", "", kv); i <- i + 1 }
  else {
    key <- kv
    nxt <- if (i < length(args)) args[i + 1] else NA
    if (key == "dry_run") { dry_run <- TRUE; i <- i + 1; next }
    if (!key %in% names(spec)) stop("Unknown option --", key, ". See ./duke --help")
    if (spec[[key]]$type == "lgl" && (is.na(nxt) || startsWith(nxt, "--"))) { val <- "TRUE"; i <- i + 1 }
    else { if (is.na(nxt)) stop("--", key, " needs a value"); val <- nxt; i <- i + 2 }
  }
  if (key == "dry_run") { dry_run <- truthy(val); next }
  if (!key %in% names(spec)) stop("Unknown option --", key, ". See ./duke --help")
  type <- spec[[key]]$type
  user[[key]] <- switch(type,
    lgl = { v <- toupper(val); if (!v %in% c("TRUE", "FALSE", "T", "F", "1", "0")) stop("--", key, " must be TRUE or FALSE"); v %in% c("TRUE", "T", "1") },
    intvec = unlist(lapply(trimws(strsplit(val, ",")[[1]]), function(x) {
      if (grepl("^[0-9]+:[0-9]+$", x)) { b <- as.integer(strsplit(x, ":")[[1]]); seq(b[1], b[2]) }
      else if (grepl("^-?[0-9]+$", x)) as.integer(x) else stop("--", key, " expects whole numbers, e.g. 1:3,6")
    })),
    val)
}

if (dry_run) {
  p <- resolve_params(user)
  tr <- resolve_tracts(p); comp <- resolve_composites(p, tr)
  cat("Duke", duke_version, "dry run: parameters that differ from the defaults\n\n")
  base <- resolve_params(list(dir_data = p$dir_data, dir_out = p$dir_out))
  for (n in names(p)) if (n %in% c("dir_data", "dir_out") || !identical(p[[n]], base[[n]]))
    cat(sprintf("  %-28s %s\n", n, paste(p[[n]], collapse = ", ")))
  cat("\nTracts (5' to 3'):\n"); print(tr[, intersect(c("tract", "motif", "path_ref", "ref_group", "count_method", "na_handling"), names(tr))], row.names = FALSE)
  if (!is.null(comp) && nrow(comp)) { cat("\nComposite measures:\n"); print(comp, row.names = FALSE) }
  cat("\nModules to run:", paste(p$run_modules, collapse = ", "), "\n")
  paths <- c(dir_data = p$dir_data, path_settings = p$path_settings, path_tracts = p$path_tracts,
             path_trim_patterns = p$path_trim_patterns, path_manual_exclusions = p$path_manual_exclusions,
             setNames(unique(tr$path_ref), rep("path_ref", length(unique(tr$path_ref)))))
  paths <- paths[!is.na(paths)]
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    cat("\nNot found - a real run would stop:\n"); cat(sprintf("  %-22s %s\n", names(missing), missing), sep = "")
    quit(status = 1)
  }
  cat("\nAll input paths found.\n")
  quit(status = 0)
}

status <- tryCatch({ run_duke(user, duke_root); 0 }, error = function(e) { message("\nDuke stopped: ", conditionMessage(e)); 1 })
quit(status = status)
