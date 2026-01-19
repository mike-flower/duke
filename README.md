# Duke Pipeline 2.1.0

A modular pipeline for amplicon sequencing analysis with comprehensive repeat length characterisation and instability metrics.

**Latest Updates (v2.1.0):**
- 🗂️ **Reorganised Directory Structure** - Cleaner root with `modules/` and `scripts/` subdirectories
- ✨ **Command-Line Interface** - Run Duke without editing files
- ✅ **Optional Trimming** - Enable/disable adapter trimming with `--trim` flag
- 🔧 **Enhanced Documentation** - Comprehensive parameter explanations with examples
- 📦 **Organised Library Files** - Clear module number prefixes (00-07)
- 🐛 **Bug Fixes** - Resume detection, parameter naming consistency, knit directory handling
- 📊 **Improved Logging** - End time display and consistent HH:MM:SS duration format

---

## Quick Start

### Three Ways to Run Duke

**1. Command-Line Interface (RECOMMENDED)** - No file editing required:
```bash
./duke --dir_data ~/data --dir_out ~/results --path_ref ~/ref.fasta
```

**2. Script-Based** - Full control with duke_run.R:
```bash
# Edit scripts/duke_run.R, then:
Rscript scripts/duke_run.R
```

**3. Interactive** - RStudio development:
```r
source("scripts/duke_run.R")
```

---

## Table of Contents

- [Installation](#installation)
- [File Structure](#file-structure)
- [Run Analysis](#run-analysis)
  - [Command-Line Interface](#1-command-line-interface-recommended)
  - [Script-Based](#2-script-based-duke_runr)
  - [Interactive RStudio](#3-interactive-rstudio)
- [HPC Deployment](#hpc-deployment)
  - [Job Submission](#job-submission)
  - [Monitoring Jobs](#monitoring-jobs)
  - [Critical Memory Requirements](#critical-memory-requirements)
  - [Resource Recommendations](#resource-recommendations)
- [Parameters](#parameters)
- [Module Overview](#module-overview)
- [Output Structure](#output-structure)
- [Troubleshooting](#troubleshooting)

---

## Installation

```r
# Core packages
install.packages(c("tidyverse", "data.table", "rmarkdown", "knitr"))

# Visualisation
install.packages(c("ggplot2", "ggrepel", "ggridges", "RColorBrewer", "cowplot", "ggnewscale"))

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
duke/
├── duke                          # CLI wrapper executable
├── README.md                     # This file
├── lib/                          # Function libraries (UNCHANGED)
│   ├── load_all.R
│   ├── 00_utils.R                # General utilities
│   ├── 01_import.R               # Module 1
│   ├── 01_sequence_qc.R          # Module 1
│   ├── 02_alignment.R            # Module 2
│   ├── 02_alignment_processing.R # Module 2
│   ├── 03_repeats.R              # Module 3
│   ├── 04_clustering.R           # Module 4
│   ├── 04_consensus.R            # Module 4
│   ├── 05_waterfall.R            # Module 5
│   ├── 06_range_analysis.R       # Module 6
│   └── 07_visualisation.R        # Module 7
├── modules/                      # Module Rmd files
│   ├── 01_import_and_qc.Rmd
│   ├── 02_alignment.Rmd
│   ├── 03_repeat_detection.Rmd
│   ├── 04_allele_calling.Rmd
│   ├── 05_waterfall.Rmd
│   ├── 06_range_analysis.Rmd
│   └── 07_repeat_visualisation.Rmd
├── scripts/                      # R and shell scripts
│   ├── duke_run.R                # Script-based runner
│   ├── duke_cli.R                # CLI script
│   ├── duke_myriad.sh            # Myriad (SGE) job script
│   ├── duke_kathleen.sh          # Old Kathleen (SGE) job script
│   └── duke_kathleen_slurm.sh    # New Kathleen (Slurm) job script
├── www/                          # Example files
│   ├── read_exclusions_example.xlsx
│   └── settings_example.xlsx
├── logs/                         # Job output logs
├── demo/                         # Demo datasets
└── archive/                      # Old versions
```

**Note:** Library files prefixed with module numbers for easy identification (00-07)

---

## Run Analysis

### 1. Command-Line Interface (RECOMMENDED)

#### Quick Start - Minimal Example
```bash
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/output \
  --path_ref /path/to/reference.fasta
```

#### Local Test Run Example
```bash
# Navigate to Duke directory
cd ~/Scratch/bin/duke

# Run Duke with test dataset
./duke \
  --dir_data demo/2025.12.17_pb_test/data \
  --dir_out demo/2025.12.17_pb_test/result_duke \
  --path_ref ~/refs/HTTset20/HTTset20.fasta \
  --path_trim_patterns ~/refs/adapters/adapters.csv \
  --threads 3 \
  --resume TRUE \
  --remove_intermediate TRUE \
  --cleanup_temp FALSE
```

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

---

### 2. Script-Based (duke_run.R)

**Local:**
```r
# Edit scripts/duke_run.R
params$dir_data <- "demo/2025.12.17_pb_test/data"
params$dir_out <- "demo/2025.12.17_pb_test/result_duke"
params$path_ref <- "~/refs/HTTset20/HTTset20.fasta"
params$trim <- FALSE  # Optional: disable trimming

# Run
source("scripts/duke_run.R")
```

**HPC:**
```bash
# Edit scripts/duke_run.R with HPC paths
nano ~/Scratch/bin/duke/scripts/duke_run.R

# Submit via job script
qsub scripts/duke_myriad.sh
```

---

### 3. Interactive RStudio

```r
# Open scripts/duke_run.R in RStudio
# Edit parameters
# Run line-by-line or source
```

---

## HPC Deployment

### Cluster Selection

| Cluster | Best For | Cores | Scheduler | Script |
|---------|----------|-------|-----------|--------|
| **Myriad** | ≤500 samples | 1-36 | SGE (`-pe smp`) | `duke_myriad.sh` |
| **Old Kathleen** | 1000+ samples | 80-160 | SGE (`-pe mpi`) | `duke_kathleen.sh` |
| **New Kathleen** | 1000+ samples | 80-160 | Slurm | `duke_kathleen_slurm.sh` |

---

### Job Submission

#### Myriad (SGE)
```bash
# Make executable (one-time)
chmod +x scripts/duke_myriad.sh

# Edit parameters in the script
nano scripts/duke_myriad.sh

# Submit
qsub scripts/duke_myriad.sh
```

#### Old Kathleen (SGE)
```bash
# Make executable (one-time)
chmod +x scripts/duke_kathleen.sh

# Edit parameters in the script
nano scripts/duke_kathleen.sh

# Submit
qsub scripts/duke_kathleen.sh
```

#### New Kathleen (Slurm)
```bash
# No chmod needed for Slurm

# Edit parameters in the script
nano scripts/duke_kathleen_slurm.sh

# Submit
sbatch scripts/duke_kathleen_slurm.sh
```

---

### Monitoring Jobs

#### Check Job Status

**SGE (Myriad, Old Kathleen):**
```bash
# Check your jobs
qstat -u $USER

# Watch queue (updates every 5 seconds)
watch -n 5 'qstat -u $USER'

# Detailed job info
qstat -j <JOB_ID>
```

**Slurm (New Kathleen):**
```bash
# Check your jobs
squeue -u $USER

# Watch queue (updates every 5 seconds)
watch -n 5 'squeue -u $USER'

# Detailed job info
scontrol show job <JOB_ID>
```

#### Monitor Log Files

```bash
# Watch Duke output log
tail -f logs/duke_<JOB_ID>.out

# Watch Duke error log
tail -f logs/duke_<JOB_ID>.err
```

#### Check Recent File Activity

Useful for checking if Duke is actively processing files:

```bash
# Files modified in last 2 hours (sorted)
watch -n 10 "find . -type f -newermt '2 hours ago' -print | sed 's|^\./||' | sort"

# Files modified in last 1 minute (sorted) - for active monitoring
watch -n 10 "find . -type f -newermt '1 minute ago' -print | sed 's|^\./||' | sort"

# One-time check (no watch)
find . -type f -newermt '30 minutes ago' -print | sed 's|^\./||' | sort
```

These commands help you verify that:
- Modules are progressing (new HTML/RData files appearing)
- Log files are being updated
- Results are being written to output directory

**Tip:** Adjust the time window based on dataset size:
- Small datasets (< 50 samples): `'1 minute ago'`
- Medium datasets (50-500 samples): `'10 minutes ago'`
- Large datasets (500+ samples): `'1 hour ago'`

#### Cancel Jobs

**SGE:**
```bash
qdel <JOB_ID>
```

**Slurm:**
```bash
scancel <JOB_ID>
```

---

### Critical Memory Requirements

**⚠️ IMPORTANT: Duke requires minimum 4GB RAM per core for Module 4 (Clustering)**

#### Empirical Testing Results (3 samples, Myriad)

| Memory/Core | Total RAM | Runtime | Status |
|-------------|-----------|---------|--------|
| **2GB** | 6GB | >4 hours | ❌ **TOO SLOW** - Stuck at clustering |
| **4GB** | 12GB | ~52 min | ✅ **MINIMUM VIABLE** |
| **8GB** | 24GB | ~28 min | ✅ **RECOMMENDED** |
| **16GB** | 48GB | ~33 min | ✅ Good (diminishing returns) |
| **32GB** | 96GB | ~38 min | ✅ Good (diminishing returns) |
| **64GB** | 192GB | ~24 min | ✅ Good (diminishing returns) |

**Key Finding:** Module 4 clustering uses parallel processing (`pbmclapply`) which forks R processes. With only 2GB/core, memory pressure causes severe slowdowns. **Minimum 4GB/core required, 8GB/core recommended.**

---

### Resource Recommendations

#### Test Runs (1-10 samples)
```bash
#$ -pe smp 3
#$ -l mem=8G        # 8GB/core recommended
#$ -l tmpfs=10G
#$ -l h_rt=2:00:00
```

#### Small Jobs (10-50 samples)
```bash
#$ -pe smp 12
#$ -l mem=8G        # 8GB/core for safe parallel clustering
#$ -l tmpfs=50G
#$ -l h_rt=6:00:00
```

#### Medium Jobs (50-200 samples) - Myriad
```bash
#$ -pe smp 36       # Maximum for Myriad
#$ -l mem=8G        # 288GB total (36 × 8GB)
#$ -l tmpfs=100G
#$ -l h_rt=24:00:00
```

#### Large Jobs (200-1000 samples) - Kathleen
```bash
#$ -pe mpi 80       # Old Kathleen
#$ -l mem=4G        # 320GB total (80 × 4GB)
#$ -l h_rt=36:00:00
```

**Or for Slurm (New Kathleen):**
```bash
#SBATCH --nodes=2
#SBATCH --ntasks=80
#SBATCH --mem-per-cpu=4G
#SBATCH --time=36:00:00
```

#### Very Large Jobs (1000+ samples) - Kathleen
```bash
#$ -pe mpi 160      # Request 160 cores
#$ -l mem=2G        # 320GB total
#$ -l h_rt=48:00:00
# Use --threads 80 in Duke (not 160!)
```

**Key Notes:**
- **Always use 4GB+ per core** for production runs
- Myriad maximum: 36 cores
- Kathleen: Request cores in multiples of 40
- Match `--threads` parameter to core count (except memory workaround above)

---

## Parameters

### New Parameters in 2.1.0

All parameters remain the same as v2.0.1. The reorganisation only changed directory structure, not functionality.

### Essential Parameters

```r
# Required
dir_data = "/path/to/data"           # Input directory
dir_out = "/path/to/output"          # Output directory
path_ref = "/path/to/reference.fasta" # Reference FASTA

# Trimming (if trim = TRUE)
path_trim_patterns = "/path/to/adapters.csv"

# Runtime
threads = 2                          # CPU cores
resume = TRUE                        # Skip completed modules (default)
```

### Complete Parameter Reference

See `./duke --help` for comprehensive documentation of all 54+ parameters, including:

#### File Paths
- `dir_data` - Input directory **(required)**
- `dir_out` - Output directory **(required)**
- `path_ref` - Reference FASTA **(required)**
- `path_trim_patterns` - Adapter file (required if `trim=TRUE`)
- `path_settings` - Settings Excel (required for Module 6)

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
All parameters have comprehensive documentation with examples in `scripts/duke_run.R` or via `./duke --help`:

- `rpt_pattern` - Repeat motif (default: "CAG")
- `rpt_min_repeats` - Minimum count (default: 2)
- `rpt_max_mismatch` - Error tolerance (default: 0)
- `rpt_start_perfect_repeats` - Perfect at start (default: 2)
- `rpt_end_perfect_repeats` - Perfect at end (default: 2)
- `rpt_max_gap` - Max gap within tract (default: 6)
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
- Coverage visualisation

### Module 3: Repeat Detection
- Identifies repeat tracts
- Tolerates sequencing errors
- Multiple counting methods

### Module 4: Allele Calling
- **⚠️ Most memory-intensive module** - requires 4GB+ per core
- Clusters by repeat/haplotype using parallel processing
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

### Module 7: Repeat Visualisation
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

**"Error: cannot open the connection" when running Rmd**
```bash
# CAUSE: Rmd files moved to modules/ subdirectory
# SOLUTION: Ensure using v2.1.0 scripts with knit_root_dir fix

# Verify you have the updated scripts:
grep "knit_root_dir" scripts/duke_run.R
grep "knit_root_dir" scripts/duke_cli.R

# Should both return lines with: knit_root_dir = getwd()
```

**Job stuck at Module 4 clustering (very slow)**
```bash
# CAUSE: Insufficient memory (2GB/core)
# SOLUTION: Use minimum 4GB/core, 8GB recommended

# Update job script:
#$ -l mem=8G

# Resubmit
qsub scripts/duke_myriad.sh
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
# threads = 80 in duke command
```

**"duke: command not found"**
```bash
# Make duke executable
chmod +x duke

# Run from duke directory
cd ~/Scratch/bin/duke
./duke --help
```

**"scripts/duke_cli.R not found"**
```bash
# CAUSE: Using old duke wrapper or wrong directory
# SOLUTION: Ensure duke wrapper points to scripts/ and run from duke root

# Check duke wrapper:
head duke

# Should show: Rscript "$SCRIPT_DIR/scripts/duke_cli.R" "$@"
```

**"Error: there is no package called 'openxlsx'"**
```bash
# Point R to your existing package library
export R_LIBS_USER=~/R/library

# Or install packages (see Installation section)
```

**Want to re-run specific modules**
```bash
# Delete module output:
rm result_duke/module_data/03_repeat_detection_results.RData

# Or run specific modules only:
./duke --run_modules 3,4,5 ...
```

**Files not in modules/ or scripts/ directories**
```bash
# CAUSE: Haven't completed reorganisation
# SOLUTION: See REORGANISATION_INSTRUCTIONS.md

# Directory should look like:
# duke/
# ├── lib/
# ├── modules/     <- All .Rmd files here
# ├── scripts/     <- All .R and .sh files here
# └── duke         <- Wrapper in root
```

---

## Version History

### 2.1.0 (January 2026) - Current
- 🗂️ **REORGANISED:** Clean directory structure with `modules/` and `scripts/` subdirectories
- 🐛 **FIXED:** `knit_root_dir` handling for Rmd files in subdirectories
- 📚 **ENHANCED:** Comprehensive reorganisation documentation
- 📊 **ADDED:** Job monitoring commands in README

### 2.0.1 (January 2025)
- ✨ **NEW:** Command-line interface
- ✅ **NEW:** `trim` parameter (optional adapter trimming)
- 🔧 **RENAMED:** `visualise_alignment_downsample` (was: `_n_reads`)
- 🐛 **FIXED:** CLI resume detection
- 🐛 **FIXED:** Resume default now TRUE (consistent)
- 📚 **ENHANCED:** Comprehensive repeat parameter documentation
- 📦 **ORGANISED:** Library files with module prefixes (00-07)
- 📊 **IMPROVED:** Logging format

### 2.0.0
- Modular architecture (7 modules)
- Resume capability
- Cleanup options

---

## Citation

If you use Duke Pipeline in your research, please cite:

Hölbling, B.V. et al. A multimodal screening platform for endogenous dipeptide repeat proteins in C9orf72 patient iPSC-neurons. *Cell Reports* (2025). https://www.sciencedirect.com/science/article/pii/S2211124725004668

---

## Contact

**Michael Flower**  
Senior Clinical Research Fellow  
Department of Neurodegenerative Disease  
UCL Queen Square Institute of Neurology  
London, UK  
WC1N 3BG

- Email: michael.flower@ucl.ac.uk  
- UCL Profile: https://profiles.ucl.ac.uk/45681-michael-flower
- ORCID: https://orcid.org/0000-0001-5568-6239
- GitHub: https://github.com/mike-flower/duke

For issues or questions:
- GitHub Issues: https://github.com/mike-flower/duke/issues
- Email: michael.flower@ucl.ac.uk

---

**Duke Pipeline 2.1.0 - Clean, organised, and production-ready!** 🎯
