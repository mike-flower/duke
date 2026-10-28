# ==============================================================================
# Duke - load all function libraries
# ==============================================================================
# Works from the Duke root, from modules/, or when sourced by path.
# ==============================================================================

.duke_lib_dir <- local({
  this <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(this) && file.exists(this)) dirname(normalizePath(this))
  else if (dir.exists("lib")) "lib"
  else if (dir.exists("../lib")) "../lib"
  else stop("Cannot locate lib/. Run from the Duke root directory.")
})

for (.f in c("00_utils.R", "params.R", "01_import.R", "01_sequence_qc.R", "02_alignment.R",
             "03_repeats.R", "04_umi.R", "05_structure.R", "06_instability.R", "07_plots.R",
             "run_pipeline.R")) {
  source(file.path(.duke_lib_dir, .f))
}
rm(.f)
