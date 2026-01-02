# ==============================================================================
# Duke Pipeline - Load All Functions
# ==============================================================================
# Master script to load all function libraries
# ==============================================================================

# Get the directory where this script is located
script_dir <- if (exists("params") && !is.null(params$default_wd)) {
  file.path(params$default_wd, "lib")
} else {
  "lib"
}

# Load function files in order
source(file.path(script_dir, "utils.R"))                # General utilities (Module 1)
source(file.path(script_dir, "import.R"))               # File import (Module 1)
source(file.path(script_dir, "sequence_qc.R"))          # Quality control (Module 1)
source(file.path(script_dir, "alignment.R"))            # Alignment functions (Module 2)
source(file.path(script_dir, "alignment_processing.R")) # Alignment utilities (Module 2)
source(file.path(script_dir, "repeats.R"))              # Repeat detection (Module 3)
source(file.path(script_dir, "clustering.R"))           # Clustering (Module 4)
source(file.path(script_dir, "consensus.R"))            # Consensus sequences (Module 4)
source(file.path(script_dir, "waterfall.R"))            # Waterfall plots (Module 5)
source(file.path(script_dir, "range_analysis.R"))       # Range analysis (Module 6)
source(file.path(script_dir, "visualisation.R"))        # Repeat visualisation (Module 7)

message("All function libraries loaded successfully")
