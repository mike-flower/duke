# Duke Pipeline 2.0.1

A modular pipeline for amplicon sequencing analysis with comprehensive repeat length characterisation and instability metrics.

**Latest Updates (v2.0.1):**
- ✨ **New Command-Line Interface** - Run Duke without editing files
- ✅ **Optional Trimming** - Enable/disable adapter trimming with `--trim` flag
- 🔧 **Enhanced Documentation** - Comprehensive parameter explanations with examples
- 📦 **Organized Library Files** - Clear module number prefixes (00-07)
- 🐛 **Bug Fixes** - Resume detection, parameter naming consistency

---

## Quick Start

### Three Ways to Run Duke

**1. Command-Line Interface (NEW!)** - No file editing required:
```bash
./duke --dir_data ~/data --dir_out ~/results --path_ref ~/ref.fasta
```

**2. Script-Based** - Full control with duke_run.R:
```bash
# Edit duke_run.R, then:
Rscript duke_run.R
```

**3. Interactive** - RStudio development:
```r
source("duke_run_local.R")
```

---

## Table of Contents

- [Installation](#installation)
- [File Structure](#file-structure)
- [Run Analysis](#run-analysis)
  - [Command-Line Interface](#1-command-line-interface-new)
  - [Script-Based](#2-script-based-duke_runr)
  - [Interactive RStudio](#3-interactive-rstudio)
- [HPC Deployment](#hpc-deployment)
- [Parameters](#parameters)
  - [New in 2.0.1](#new-parameters-in-201)
  - [Essential Parameters](#essential-parameters)
  - [Complete Reference](#complete-parameter-reference)
- [Module Overview](#module-overview)
- [Output Structure](#output-structure)
- [Troubleshooting](#troubleshooting)

---

## Installation

```r
# Core packages
install.packages(c("tidyverse", "data.table", "rmarkdown", "knitr"))

# Visualization
install.packages(c("ggplot2", "ggrepel", "ggridges", "RColorBrewer", "cowplot"))

# Tables and reports
install.packages(c("openxlsx", "DT", "kableExtra", "htmltools"))

# Analysis
install.packages(c("mclust", "cluster", "ineq", "moments", "pracma"))

# Parallel processing
install.packages(c("pbapply", "pbmcapply"))

# Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biostrings", "ShortRead", "GenomicAlignments", 
                       "Rsamtools", "DECIPHER"))
```

**HPC Setup:** Point R to your existing package library:
```bash
export R_LIBS_USER=~/R/library
```

---

## File Structure

```
duke_pipeline/
├── duke_cli.R                    # CLI script (NEW!)
├── duke                          # CLI wrapper (NEW!)
├── duke_run.R                    # HPC runner
├── duke_run_local.R              # Local runner
├── duke_myriad_cli.sh            # Myriad CLI job script (NEW!)
├── duke_kathleen_cli.sh          # Kathleen CLI job script (NEW!)
├── duke_myriad.sh                # Myriad job script
├── duke_kathleen.sh              # Kathleen job script
├── 01_import_and_qc.Rmd through 07_repeat_visualisation.Rmd
└── lib/                          # Function libraries
    ├── load_all.R
    ├── 00_utils.R                # General utilities
    ├── 01_import.R               # Module 1
    ├── 01_sequence_qc.R          # Module 1
    ├── 02_alignment.R            # Module 2
    ├── 02_alignment_processing.R # Module 2
    ├── 03_repeats.R              # Module 3
    ├── 04_clustering.R           # Module 4
    ├── 04_consensus.R            # Module 4
    ├── 05_waterfall.R            # Module 5
    ├── 06_range_analysis.R       # Module 6
    └── 07_visualisation.R        # Module 7
```

**Note:** Library files prefixed with module numbers for easy identification

---

## Run Analysis

### 1. Command-Line Interface (NEW!)

#### Quick Start - Minimal Example
```bash
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/output \
  --path_ref /path/to/reference.fasta
```

#### Complete Test Run Example

**Prerequisites on HPC:**
```bash
# Load required modules
module load r/recommended
module load samtools/1.11/gnu-4.9.2

# Point R to your package library
export R_LIBS_USER=~/R/library

# Add minimap2 to PATH
export PATH=$HOME/Scratch/bin/minimap2:$PATH
```

**Run Duke with test dataset:**
```bash
# Navigate to Duke directory
cd /home/skgtmdf/Scratch/bin/duke

# Run test analysis with 3 threads
./duke \
  --dir_data /home/skgtmdf/Scratch/data/2025.12.17_pb_test/data \
  --dir_out /home/skgtmdf/Scratch/data/2025.12.17_pb_test/result_duke \
  --path_ref /home/skgtmdf/Scratch/refs/HTTset20/HTTset20.fasta \
  --path_trim_patterns /home/skgtmdf/Scratch/refs/adapters/adapters.csv \
  --path_settings /home/skgtmdf/Scratch/data/2025.12.17_pb_test/settings/settings_duke.xlsx \
  --threads 3 \
  --resume TRUE \
  --remove_intermediate TRUE \
  --cleanup_temp FALSE
```

**Notes:**
- Use `--threads 3` for login node testing (never use high thread counts on login nodes)
- For production runs, submit to queue with job scripts (see [HPC Deployment](#hpc-deployment))
- Results will be in `/home/skgtmdf/Scratch/data/2025.12.17_pb_test/result_duke/`

#### With Adapter Trimming
```bash
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/output \
  --path_ref /path/to/reference.fasta \
  --path_trim_patterns /path/to/adapters.csv
```

#### Without Adapter Trimming
```bash
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/output \
  --path_ref /path/to/reference.fasta \
  --trim FALSE
```

#### See All Options
```bash
./duke --help
```

**Key Features:**
- ✅ 54+ configurable parameters
- ✅ Resume enabled by default (skip completed modules)
- ✅ Comprehensive help with examples
- ✅ No file editing required

**⚠️ Important:** Don't run Duke directly on login nodes for production datasets. Use job scripts via `qsub` or `sbatch`.

---

#### HPC Cluster Test Run (Myriad/Kathleen)

**For running test data on compute nodes via job scripts:**

Edit `duke_myriad_cli.sh` or `duke_kathleen_cli.sh` and update the `./duke` command:

```bash
./duke \
  --dir_data /home/skgtmdf/Scratch/data/2025.12.17_pb_test/data \
  --dir_out /home/skgtmdf/Scratch/data/2025.12.17_pb_test/result_duke \
  --path_ref /home/skgtmdf/Scratch/refs/HTTset20/HTTset20.fasta \
  --path_trim_patterns /home/skgtmdf/Scratch/refs/adapters/adapters.csv \
  --path_settings /home/skgtmdf/Scratch/data/2025.12.17_pb_test/settings/settings_duke.xlsx \
  --threads 36 \
  --resume TRUE \
  --remove_intermediate TRUE \
  --cleanup_temp FALSE
```

Then submit:
```bash
# Myriad (36 cores)
qsub duke_myriad_cli.sh

# Kathleen (80 cores - change --threads to 80 in script)
qsub duke_kathleen_cli.sh
```

**Monitor the job:**
```bash
# Check status
qstat -u $USER

# Watch log
tail -f logs/duke_<JOB_ID>.out
```

---

### 2. Script-Based (duke_run.R)

**Local:**
```r
# Edit duke_run_local.R
params$dir_data <- "/path/to/data"
params$dir_out <- "/path/to/output"
params$path_ref <- "/path/to/reference.fasta"
params$trim <- FALSE  # Optional: disable trimming

# Run
source("duke_run_local.R")
```

**HPC:**
```bash
# Edit duke_run.R with HPC paths
nano ~/Scratch/bin/duke/duke_run.R

# Submit
qsub duke_myriad.sh
```

---

### 3. Interactive RStudio

```r
# Open duke_run_local.R in RStudio
# Edit parameters
# Run line-by-line or source
```

---

## HPC Deployment

### Cluster Selection

| Cluster | Best For | Cores | Scheduler |
|---------|----------|-------|-----------|
| **Myriad** | ≤500 samples | 1-36 | SGE (`-pe smp`) |
| **Kathleen** | 1000+ samples | 80-160 | SGE (`-pe mpi`) |

### Job Script Setup

**IMPORTANT:** All HPC job scripts must include:

```bash
# After module loads, add these lines:
export R_LIBS_USER=~/R/library
export PATH=$HOME/Scratch/bin/minimap2:$PATH
```

### Memory Optimization (Kathleen)

For large datasets, request more cores than you use:

```bash
#$ -pe mpi 160     # Request 160 cores
#$ -l mem=2G       # 2GB/core = 320GB total

# In duke_run.R or CLI:
threads = 80       # Use only 80 threads
```

**Why:** Provides memory buffer (320GB / 80 = 4GB per thread) to prevent fork failures.

### Job Submission

**Traditional (Edit duke_run.R):**
```bash
# Myriad
qsub duke_myriad.sh

# Kathleen  
qsub duke_kathleen.sh
```

**CLI (Edit job script parameters):**
```bash
# Myriad
qsub duke_myriad_cli.sh

# Kathleen
qsub duke_kathleen_cli.sh
```

**Performance:** ~15-18 hours for 1000 samples (Kathleen 80 cores)

---

## Parameters

### New Parameters in 2.0.1

```r
# Optional adapter trimming (NEW!)
trim = TRUE                          # Enable/disable trimming
                                    # Set FALSE to skip trimming

# Renamed for consistency
visualise_alignment_downsample = 1000  # Was: visualise_alignment_n_reads
                                      # Max reads to plot (NA = all)
```

### Essential Parameters

```r
# Required
dir_data = "/path/to/data"           # Input directory
dir_out = "/path/to/output"          # Output directory
path_ref = "/path/to/reference.fasta" # Reference FASTA
path_settings = "/path/to/settings.xlsx" # Settings file

# Trimming (if trim = TRUE)
path_trim_patterns = "/path/to/adapters.csv"

# Runtime
threads = 2                          # CPU cores
resume = TRUE                        # Skip completed modules (default)
```

### Complete Parameter Reference

#### File Paths
- `dir_data` - Input directory **(required)**
- `dir_out` - Output directory **(required)**
- `path_ref` - Reference FASTA **(required)**
- `path_trim_patterns` - Adapter file (required if `trim=TRUE`)
- `path_settings` - Settings Excel **(required for Module 6)**

#### Import Options
- `import_patterns` - File extensions to import
- `downsample` - Limit reads per sample (NA = all)
- `select_one_of_pair` - For paired-end: "R1", "R2", or NA

#### Adapter Trimming
- `trim` - Enable/disable (default: TRUE)
- `trim_max_mismatch` - Max errors (default: 3)
- `trim_with_indels` - Allow indels (default: TRUE)

#### Alignment
- `minimap2_args` - minimap2 settings (default: "-t 2 -x sr -w 1")
- `visualise_alignment` - Generate plots (default: TRUE)
- `visualise_alignment_downsample` - Max reads to plot (default: 1000)

#### Repeat Detection
All parameters have comprehensive documentation with examples in `duke_run.R` or via `./duke --help`:

- `rpt_pattern` - Repeat motif (default: "CAG")
- `rpt_min_repeats` - Minimum count (default: 2)
- `rpt_max_mismatch` - Error tolerance (default: 0)
- `rpt_start_perfect_repeats` - Perfect at start (default: 2)
- `rpt_end_perfect_repeats` - Perfect at end (default: 2)
- `rpt_max_gap` - Max gap within tract (default: 6)
  - Example: "CAG-TA-CAG" = 2bp gap
- `rpt_max_tract_gap` - Max gap between tracts (default: 18)
- `rpt_return_option` - "longest" or "all"

#### Clustering
- `cluster` - Enable clustering (default: TRUE)
- `cluster_by` - Method: "repeat", "haplotype", or both
- `haplotype_cluster_max` - Max haplotypes (default: 10)
- `repeat_cluster_max` - Max repeat lengths (default: 20)

#### Module Control
- `waterfall` - Generate waterfall plots (default: TRUE)
- `repeat_histogram` - Generate histograms (default: TRUE)
- `repeat_scatter` - Generate scatter plots (default: TRUE)

#### Runtime
- `threads` - CPU cores (default: 80 HPC / 2 local)
- `resume` - Skip completed modules (default: TRUE)
- `remove_intermediate` - Free RAM (default: TRUE)
- `cleanup_temp` - Delete temp files (default: FALSE)
- `run_modules` - Which modules to run (default: c(1:7))
- `log_dir` - Log directory (default: "logs")

---

## Module Overview

### Module 1: Import and QC
- Imports FASTQ/FASTA/BAM files
- **Optionally trims adapters** (controlled by `trim` parameter)
- Removes duplicates
- QC metrics and plots

### Module 2: Alignment
- minimap2 alignment to reference
- Strand correction
- Coverage visualization

### Module 3: Repeat Detection
- Identifies repeat tracts
- Tolerates sequencing errors
- Multiple counting methods

### Module 4: Allele Calling
- Clusters by repeat/haplotype
- Consensus sequences
- Variant calling (VCF export)

### Module 5: Waterfall Plots
- Visual read inspection
- Per-sample and per-cluster plots
- Configurable downsampling

### Module 6: Range Analysis
- Modal peak detection
- Instability metrics
- Control comparisons
- Excel export

### Module 7: Repeat Visualization
- Frequency histograms
- Scatter/violin plots
- Publication-ready figures

---

## Output Structure

```
result_duke/
├── 01_import_and_qc.html             # HTML reports
├── 02_alignment.html
├── ...
├── 07_repeat_visualisation.html
│
├── module_data/                      # Cached results (for resume)
│   ├── 01_import_qc_results.RData
│   ├── 02_alignment_results.RData
│   └── ...
│
├── 01_import_qc/                     # Module outputs
│   ├── qc.xlsx
│   └── plots/
├── 02_alignment/
│   └── plots/coverage_by_sample/
├── ...
└── 07_repeat_visualisation/
    └── plots/
```

---

## Troubleshooting

### Common Issues

**"Rscript: command not found" (CLI on HPC)**
```bash
# Load R module first
module load r/recommended

# Also set library path
export R_LIBS_USER=~/R/library
```

**"Error: there is no package called 'openxlsx'"**
```bash
# Point R to your existing package library
export R_LIBS_USER=~/R/library

# Or install packages (see Installation section)
```

**"Error: --path_ref is required"**
```bash
# path_ref is now required (no default)
./duke --path_ref /path/to/reference.fasta ...
```

**"Cannot allocate memory" (Kathleen)**
```bash
# Request more cores than you use:
#$ -pe mpi 160
#$ -l mem=2G
# threads = 80 in duke_run.R
```

**"visualise_alignment_n_reads not found"**
```r
# Renamed in 2.0.1:
visualise_alignment_downsample = 1000  # New name
```

**"path_trim_patterns required"**
```bash
# Either provide it:
./duke --path_trim_patterns /path/to/adapters.csv ...

# Or disable trimming:
./duke --trim FALSE ...
```

**Resume not working (CLI)**
```bash
# Resume is enabled by default (changed in 2.0.1)
# To force re-run:
./duke --resume FALSE ...
```

**Want to re-run specific modules**
```bash
# Delete module output:
rm result_duke/module_data/03_repeat_detection_results.RData

# Or run specific modules only:
./duke --run_modules 3,4,5 ...
```

---

## Version History

### 2.0.1 (January 2025) - Current
- ✨ **NEW:** Command-line interface
- ✅ **NEW:** `trim` parameter (optional adapter trimming)
- 🔧 **RENAMED:** `visualise_alignment_downsample` (was: `_n_reads`)
- 🐛 **FIXED:** CLI resume detection
- 🐛 **FIXED:** Resume default now TRUE (consistent)
- 📚 **ENHANCED:** Comprehensive repeat parameter documentation
- 📦 **ORGANIZED:** Library files with module prefixes (00-07)
- 🔧 **UPDATED:** Conditional trimming in Module 1
- ❌ **REMOVED:** `visualise_alignment_corrected` (didn't exist)

### 2.0.0
- Modular architecture (7 modules)
- Resume capability
- Cleanup options

---

## Citation

[Citation information to be added]

---

## Contact

For issues or questions:
- GitHub Issues: [repository URL]
- Email: [contact email]

---

**Duke Pipeline 2.0.1 - Ready for production!** 🎯
