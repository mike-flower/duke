# ==============================================================================
# Duke - write the parameter tables in README.md from lib/params.R
# ==============================================================================
# Each table sits between two markers in README.md:
#   <!-- params:<group>:start -->  ...  <!-- params:<group>:end -->
# where <group> is the group name used in lib/params.R (Paths, Import,
# Alignment, Repeats, UMI, Allele structure, Instability, Repeat lengths,
# Runtime). Run after adding or changing a parameter:
#
#   Rscript scripts/update_readme_params.R
#
# The script stops if a group has no markers in README.md, so no parameter can
# be left undocumented.
# ==============================================================================

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
duke_root <- if (length(file_arg)) normalizePath(file.path(dirname(file_arg), "..")) else getwd()
source(file.path(duke_root, "lib", "params.R"))

type_label <- c(chr = "text", int = "whole number", num = "number", lgl = "TRUE / FALSE",
                chrvec = "text list", intvec = "whole-number list", numvec = "number list")

format_default <- function(x, type) {
  if (length(x) == 0) return("–")
  if (type == "chr" && is.na(x[1])) return("–")
  if (length(x) == 1 && is.na(x)) return("`NA`")
  if (type == "intvec" && length(x) > 2 && all(diff(x) == 1)) return(sprintf("`%d:%d`", x[1], x[length(x)]))
  paste0("`", paste(as.character(x), collapse = ", "), "`")
}

format_help <- function(h) {
  h <- gsub("'([^']+)'", "`\\1`", h)   # 'longest' -> `longest`
  gsub("\\|", "\\\\|", h)
}

make_table <- function(rows) {
  c("| Parameter | Default | Format | What it does |",
    "|---|---|---|---|",
    vapply(rows, function(r) sprintf("| `%s` | %s | %s | %s |", r$name, format_default(r$default, r$type),
                                     type_label[[r$type]], format_help(r$help)), ""))
}

spec <- duke_param_spec()
groups <- unique(vapply(spec, `[[`, "", "group"))
readme_path <- file.path(duke_root, "README.md")
readme <- readLines(readme_path, warn = FALSE)

for (g in groups) {
  start <- which(readme == sprintf("<!-- params:%s:start -->", g))
  end <- which(readme == sprintf("<!-- params:%s:end -->", g))
  if (length(start) != 1 || length(end) != 1 || end < start)
    stop("README.md needs one pair of markers for parameter group '", g, "'")
  rows <- Filter(function(r) r$group == g, spec)
  readme <- c(readme[seq_len(start)], "", make_table(rows), "", readme[end:length(readme)])
}

writeLines(readme, readme_path)
cat("Updated", length(groups), "parameter tables (", length(spec), "parameters ) in", readme_path, "\n")
