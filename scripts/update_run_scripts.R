# ==============================================================================
# Duke - write the parameter lists in scripts/duke_run.R and duke_run_local.R
# ==============================================================================
# Rebuilds each script's `params <- list(...)` from lib/params.R, so every parameter
# is listed with its current default and description. Each script keeps its own
# dir_data, dir_out, path_ref, path_settings, path_trim_patterns and threads.
# Run from the Duke folder after adding or changing a parameter:
#
#   Rscript scripts/update_run_scripts.R
# ==============================================================================
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
setwd(if (length(file_arg)) normalizePath(file.path(dirname(file_arg), "..")) else getwd())
source("lib/params.R")
spec <- duke_param_spec()
lit <- function(x, type) {
  if (length(x) == 0) return(switch(type, chrvec = "character(0)", intvec = "integer(0)", numvec = "numeric(0)", "character(0)"))
  if (type %in% c("int", "intvec")) {
    if (all(is.na(x))) return("NA_integer_")
    if (length(x) > 2 && all(diff(x) == 1)) return(sprintf("%d:%d", x[1], x[length(x)]))
    v <- paste0(x, "L"); return(if (length(v) == 1) v else paste0("c(", paste(v, collapse = ", "), ")"))
  }
  if (type %in% c("num", "numvec")) { if (all(is.na(x))) return("NA_real_"); v <- format(x, scientific = FALSE, trim = TRUE); return(if (length(v) == 1) v else paste0("c(", paste(v, collapse = ", "), ")")) }
  if (type == "lgl") return(if (is.na(x)) "NA" else as.character(x))
  if (all(is.na(x))) return("NA_character_")
  paste(deparse(x, width.cutoff = 500L), collapse = "")
}
for (f in c("scripts/duke_run.R", "scripts/duke_run_local.R")) {
  L <- readLines(f); a <- grep("^params <- list\\(", L); b <- grep("^\\)$", L); b <- b[b > a][1]
  old <- L[(a + 1):(b - 1)]
  keep <- c("dir_data", "dir_out", "path_ref", "path_settings", "path_trim_patterns", "threads")
  own <- list(); for (k in keep) { i <- grep(paste0("^\\s+", k, " = "), old); if (length(i)) own[[k]] <- sub(",\\s*(#.*)?$", "", sub(paste0("^\\s+", k, " = "), "", old[i[1]])) }
  out <- character(0); g0 <- ""
  for (r in spec) {
    if (r$group != g0) { g0 <- r$group; out <- c(out, if (length(out)) "", sprintf("  # --- %s %s", g0, strrep("-", max(3, 74 - nchar(g0))))) }
    out <- c(out, paste0("  # ", strwrap(r$help, 86, exdent = 0)))
    val <- if (!is.null(own[[r$name]])) own[[r$name]] else lit(r$default, r$type)
    out <- c(out, sprintf("  %s = %s,", r$name, val))
  }
  out[length(out)] <- sub(",$", "", out[length(out)])
  L2 <- c(L[1:a], out, L[b:length(L)])
  writeLines(L2, f)
  cat(f, ": ", length(spec), " parameters; kept own values for ", paste(names(own), collapse = ", "), "\n", sep = "")
}
