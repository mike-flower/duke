# ==============================================================================
# Duke Pipeline - Load All Functions
# ==============================================================================
# Master script to load all function libraries
# Files are prefixed with module numbers (00-07) for easy identification
# ==============================================================================

# Intelligently find the lib directory
# Works from: duke root, modules/ subdirectory, or via pipeline rendering
find_lib_dir <- function() {
  # Try to detect where this script is located
  this_script <- tryCatch({
    # Method 1: sys.frame (works when sourced)
    sys.frame(1)$ofile
  }, error = function(e) NULL)
  
  if (!is.null(this_script) && file.exists(this_script)) {
    # We found the script location - use that directory
    return(dirname(normalizePath(this_script)))
  }
  
  # Method 2: Check current working directory and common locations
  if (dir.exists("lib")) {
    return("lib")
  } else if (dir.exists("../lib")) {
    return("../lib")
  } else {
    stop("Cannot locate lib directory. Ensure you're running from duke root or modules/ subdirectory.")
  }
}

script_dir <- find_lib_dir()

# Load function files in order (00-07 prefix indicates primary module)
source(file.path(script_dir, "00_utils.R"))                   # General utilities (all modules)
source(file.path(script_dir, "01_import.R"))                  # File import (Module 1)
source(file.path(script_dir, "01_sequence_qc.R"))             # Quality control (Module 1)
source(file.path(script_dir, "02_alignment.R"))               # Alignment functions (Module 2)
source(file.path(script_dir, "02_alignment_processing.R"))    # Alignment utilities (Module 2)
source(file.path(script_dir, "03_repeats.R"))                 # Repeat detection (Module 3)
source(file.path(script_dir, "04_clustering.R"))              # Clustering (Module 4)
source(file.path(script_dir, "04_consensus.R"))               # Consensus sequences (Module 4)
source(file.path(script_dir, "05_waterfall.R"))               # Waterfall plots (Module 5)
source(file.path(script_dir, "06_range_analysis.R"))          # Range analysis (Module 6)
source(file.path(script_dir, "07_visualisation.R"))           # Repeat visualization (Module 7)

message("All function libraries loaded successfully")
message("  00_utils, 01_import, 01_sequence_qc")
message("  02_alignment, 02_alignment_processing")
message("  03_repeats, 04_clustering, 04_consensus")
message("  05_waterfall, 06_range_analysis, 07_visualisation")
