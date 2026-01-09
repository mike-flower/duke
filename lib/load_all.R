# ==============================================================================
# Duke Pipeline - Load All Functions
# ==============================================================================
# Master script to load all function libraries
# Files are prefixed with module numbers (00-07) for easy identification
# ==============================================================================

# Get the directory where this script is located
script_dir <- if (exists("params") && !is.null(params$default_wd)) {
  file.path(params$default_wd, "lib")
} else {
  "lib"
}

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
