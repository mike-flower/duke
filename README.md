# Duke Pipeline 2.0

A modular pipeline for amplicon sequencing analysis with comprehensive repeat length characterisation and instability metrics.

## Table of Contents

**Quick Start:**
- [Installation](#installation)
- [File Structure](#file-structure)
- [Run Analysis](#run-analysis)
- [HPC Deployment](#hpc-deployment-eg-ucl-myriad)
  - [Quick Start for First-Time Users](#quick-start-for-first-time-hpc-users)
  - [Initial Setup](#initial-setup-on-myriad)
  - [Updating Duke from GitHub](#updating-duke-from-github)
  - [Installing R Dependencies](#installing-r-dependencies-on-hpc)
  - [Configure Duke for HPC](#configure-duke-for-hpc)
  - [Submit and Monitor Jobs](#submit-and-monitor-jobs)
- [Review Results](#review-results)

**Configuration:**
- [Essential Parameters](#essential-parameters)
- [Settings File Format](#settings-file-format)
- [Complete Parameter Reference](#complete-parameter-reference)

**Understanding Duke:**
- [Module Overview](#module-overview)
- [Control Sample Analysis](#control-sample-analysis)
- [Analysis Ranges](#analysis-ranges)
- [Output Structure](#output-structure)

**Reference:**
- [Troubleshooting](#troubleshooting)
- [Advanced Configuration](#advanced-configuration)
- [Package Dependencies](#package-dependencies)

---

## Features

- **7 Analysis Modules**: From import to advanced visualisation
- **Flexible Parameters**: Extensive customisation options  
- **Multiple Outputs**: HTML reports, plots, Excel tables, VCF files
- **Repeat Analysis**: Modal peaks, instability metrics, distribution analysis
- **Visualisation**: Waterfall plots, histograms, scatter/violin plots
- **Resume Capability**: Skip completed modules automatically
- **Auto Cleanup**: Optional removal of intermediate files

---

## Quick Start

### Installation

```r
# Core tidyverse and data manipulation
install.packages(c("tidyverse", "data.table", "dplyr", "tidyr", "stringr", 
                   "tibble", "readxl", "plyr"))

# Plotting and visualisation
install.packages(c("ggplot2", "ggrepel", "ggridges", "ggnewscale", "scales", 
                   "RColorBrewer", "viridisLite", "cowplot"))

# Reporting and tables
install.packages(c("rmarkdown", "knitr", "kableExtra", "DT", "htmltools", 
                   "openxlsx"))

# Statistical analysis
install.packages(c("mclust", "cluster", "ineq", "moments", "pracma"))

# Parallel processing
install.packages(c("pbapply", "pbmcapply"))

# Utilities
install.packages(c("rstudioapi"))

# Install Bioconductor manager and packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biostrings", "ShortRead", "GenomicAlignments", 
                       "Rsamtools", "DECIPHER"))
```

**Note**: The `tidyverse` package includes `ggplot2`, `dplyr`, `tidyr`, `stringr`, and `tibble`, but they are listed explicitly above for clarity.

### File Structure

```
duke_pipeline/
├── duke_run.R                    # Main runner script
├── 01_import_and_qc.Rmd          # Module 1
├── 02_alignment.Rmd              # Module 2
├── 03_repeat_detection.Rmd       # Module 3
├── 04_allele_calling.Rmd         # Module 4
├── 05_waterfall.Rmd              # Module 5
├── 06_range_analysis.Rmd         # Module 6
├── 07_repeat_visualisation.Rmd   # Module 7
└── lib/                          # Function libraries
    ├── load_all.R
    ├── utils.R
    ├── import.R
    ├── sequence_qc.R
    ├── alignment.R
    ├── alignment_processing.R
    ├── repeats.R
    ├── clustering.R
    ├── consensus.R
    ├── waterfall.R
    ├── range_analysis.R
    └── visualisation.R
```

### Run Analysis

**From R/RStudio:**
```r
# Edit duke_run.R with your paths
params$dir_data <- "/path/to/your/fastq/files"
params$dir_out <- "/path/to/output/directory"
params$path_ref <- "/path/to/reference.fasta"
params$path_settings <- "/path/to/settings.xlsx"  # Required for Module 6

# Run the pipeline
source("duke_run.R")
```

**From command line:**
```bash
# Navigate to pipeline directory
cd /path/to/duke_pipeline

# Run the pipeline
Rscript duke_run.R

# For long runs (background)
nohup Rscript duke_run.R &
```

All runs are automatically logged to `logs/TIMESTAMP/` directories.

---

## HPC Deployment (e.g. UCL Myriad)

Duke can be deployed on HPC clusters for large datasets. This section covers setup, dependency installation, and job submission.

### Quick Start for First-Time HPC Users

**Complete workflow from scratch to running job:**

```bash
# 1. Connect to Myriad
ssh your_username@myriad.rc.ucl.ac.uk

# 2. Verify prerequisites are available
module purge
module load r/recommended
R --version  # Should show R 4.2.0+
module load samtools/1.11/gnu-4.9.2
samtools --version  # Should show samtools 1.11

# 3. Set up directories and Duke pipeline
cd ~/Scratch/bin
git clone https://github.com/your_repo/duke_pipeline.git duke
cd duke
chmod +x duke_myriad.sh
mkdir -p logs

# To update Duke later (if already cloned):
# git stash                 # Save local changes
# git pull origin main      # Get updates
# git stash pop             # Restore local changes

# 4. Set up R environment
echo 'export R_LIBS_USER=~/R/library' >> ~/.bashrc
mkdir -p ~/R/library
source ~/.bashrc

# 5. Install minimap2
cd ~/Scratch/bin
git clone https://github.com/lh3/minimap2
cd minimap2 && make
echo 'export PATH=$HOME/Scratch/bin/minimap2:$PATH' >> ~/.bashrc
source ~/.bashrc
minimap2 --version  # Test it works

# 6. Install R packages (see Installing R Dependencies section below)
qrsh -l h_rt=2:00:00,mem=8G
module load r/recommended
R
# ... follow R installation instructions ...
quit(save = "no")
exit  # Exit interactive session

# 7. Configure Duke for your data
cd ~/Scratch/bin/duke
nano duke_run_myriad.R  # Edit paths to your data

# 8. Submit job
qsub duke_myriad.sh

# 9. Monitor progress
qstat -u $USER  # Check if job is running
tail -f logs/duke_job_<JOB_ID>.out  # Watch live output
```

That's it! For detailed explanations of each step, see sections below.

---

### Initial Setup on Myriad

**1. Connect to Myriad and verify prerequisites:**
```bash
# SSH to Myriad
ssh your_username@myriad.rc.ucl.ac.uk

# Verify R is available
module purge
module load r/recommended
R --version
# Should show R version 4.2.0 or higher

# Verify samtools is available
module load samtools/1.11/gnu-4.9.2
samtools --version
# Should show samtools 1.11

# Exit R if you started it
# (press Ctrl+D or type q())
```

**2. Transfer Duke to HPC:**
```bash
# Option A: Transfer from local machine
# (run this from your LOCAL machine, not Myriad)
scp -r duke_pipeline/ your_username@myriad.rc.ucl.ac.uk:~/Scratch/bin/

# Option B: Clone from GitHub repository
# (run this ON Myriad after SSH)
cd ~/Scratch/bin
git clone https://github.com/your_repo/duke_pipeline.git duke
cd duke
```

**3. Make job script executable:**
```bash
cd ~/Scratch/bin/duke
chmod +x duke_myriad.sh

# Verify it's executable
ls -lh duke_myriad.sh
# Should show -rwxr-xr-x (note the 'x' for executable)
```

**4. Set up R library path:**
```bash
# Add to ~/.bashrc for persistence
echo 'export R_LIBS_USER=~/R/library' >> ~/.bashrc
source ~/.bashrc

# Create R library directory
mkdir -p ~/R/library
```

**5. Install minimap2 (not available as module on Myriad):**
```bash
# Navigate to bin directory
cd ~/Scratch/bin

# Download and compile minimap2
git clone https://github.com/lh3/minimap2
cd minimap2
make

# Add to PATH (add to ~/.bashrc for persistence)
echo 'export PATH=$HOME/Scratch/bin/minimap2:$PATH' >> ~/.bashrc
source ~/.bashrc

# Test installation
minimap2 --version
# Should show: 2.xx-rxxx
```

### Updating Duke from GitHub

When bug fixes or new features are released, update your Duke installation:

**If you have NO local modifications:**
```bash
cd ~/Scratch/bin/duke
git pull origin main
```

**If you have local modifications (e.g., edited duke_run_myriad.R):**
```bash
cd ~/Scratch/bin/duke

# Save your local changes temporarily
git stash

# Pull latest updates
git pull origin main

# Restore your local changes
git stash pop

# If there are conflicts, git will notify you
# Edit the conflicting files to resolve, then:
git add <conflicting_file>
git stash drop
```

**Recommended workflow to avoid conflicts:**

Keep your configuration separate from the Duke codebase:
```bash
cd ~/Scratch/bin/duke

# Make a copy of your configuration
cp duke_run_myriad.R ~/duke_run_myriad_MY_CONFIG.R

# Now you can safely pull updates without stashing
git pull origin main

# After updates, copy your config back
cp ~/duke_run_myriad_MY_CONFIG.R duke_run_myriad.R
```

**Check what version you have:**
```bash
cd ~/Scratch/bin/duke
git log --oneline -5  # Show last 5 commits
head -20 README.md | grep "Version"  # Check version in README
```

**Updating after bug fixes (e.g., v2.0.1):**
```bash
cd ~/Scratch/bin/duke

# Check current status
git status

# If you have uncommitted changes to config files only:
git stash push duke_run_myriad.R -m "My config"

# Pull the bug fixes
git pull origin main

# The following files should update automatically:
# - lib/consensus.R (Module 4 fix)
# - 06_range_analysis.Rmd (Module 6 fix)  
# - README.md (updated documentation)

# Restore your config
git stash pop

# Verify you have the fixes
grep -n "includeNonLetters" lib/consensus.R
# Should NOT appear (bug was removed)

grep -n "openxlsx" 06_range_analysis.Rmd  
# Should appear on line 23 (replaced writexl)
```

### Installing R Dependencies on HPC

Load R and install packages in an interactive session:

```bash
# Request interactive node with more resources for installation
qrsh -l h_rt=2:00:00,mem=8G

# Load R
module load r/recommended
R
```

Then in R:
```r
# Set library path
.libPaths("~/R/library")

# Install CRAN packages
install.packages(c("tidyverse", "data.table", "dplyr", "tidyr", "stringr", 
                   "tibble", "readxl", "plyr", "ggplot2", "ggrepel", 
                   "ggridges", "ggnewscale", "scales", "RColorBrewer", 
                   "viridisLite", "cowplot", "rmarkdown", "knitr", 
                   "kableExtra", "DT", "htmltools", "openxlsx", "mclust", 
                   "cluster", "ineq", "moments", "pracma", "pbapply", 
                   "pbmcapply"),
                repos = "https://cloud.r-project.org")

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install(c("Biostrings", "ShortRead", "GenomicAlignments", 
                       "Rsamtools", "DECIPHER"))

# Verify installation
required_cran <- c("tidyverse", "data.table", "dplyr", "tidyr", "stringr", 
                   "tibble", "readxl", "plyr", "ggplot2", "ggrepel", 
                   "ggridges", "ggnewscale", "scales", "RColorBrewer", 
                   "viridisLite", "cowplot", "rmarkdown", "knitr", 
                   "kableExtra", "DT", "htmltools", "openxlsx", "mclust", 
                   "cluster", "ineq", "moments", "pracma", "pbapply", 
                   "pbmcapply")

required_bioc <- c("Biostrings", "ShortRead", "GenomicAlignments", 
                   "Rsamtools", "DECIPHER")

missing_cran <- required_cran[!sapply(required_cran, requireNamespace, quietly = TRUE)]
missing_bioc <- required_bioc[!sapply(required_bioc, requireNamespace, quietly = TRUE)]

if (length(missing_cran) == 0 && length(missing_bioc) == 0) {
  cat("✓ All packages installed successfully!\n")
} else {
  cat("Missing packages:\n")
  print(c(missing_cran, missing_bioc))
}

quit(save = "no")
```

### Configure Duke for HPC

**1. Create/edit `duke_run_myriad.R`:**
```r
# Copy and modify duke_run.R for Myriad paths
params$dir_data <- "~/Scratch/data/your_project/fastq_files"
params$dir_out <- "~/Scratch/data/your_project/result_duke"
params$path_ref <- "~/Scratch/data/your_project/reference.fasta"
params$path_settings <- "~/Scratch/data/your_project/settings.xlsx"

# Match threads to job script
params$threads <- 3

# Enable for large datasets
params$remove_intermediate <- TRUE
params$cleanup_temp <- FALSE  # Keep for debugging, enable later
```

**2. Create job submission script (`duke_myriad.sh`):**

See included `duke_myriad.sh` template. Key parameters:
```bash
#$ -N duke_pipeline          # Job name
#$ -wd ~/Scratch/bin/duke    # Working directory
#$ -l h_rt=24:00:00          # 24 hour walltime
#$ -pe smp 3                 # 3 CPU cores (match params$threads)
#$ -l mem=8G                 # 8GB per CPU (24GB total)
#$ -l tmpfs=50G              # Temp space for BAM files
```

**3. Create logs directory:**
```bash
cd ~/Scratch/bin/duke
mkdir -p logs
```

**4. Verify setup is complete:**

Before submitting your first job, check that everything is ready:

```bash
# Checklist - run these commands and verify output
module purge
module load r/recommended && R --version          # ✓ R 4.2.0+
module load samtools/1.11/gnu-4.9.2 && samtools --version  # ✓ samtools 1.11
minimap2 --version                                 # ✓ minimap2 2.xx
ls -lh ~/Scratch/bin/duke/duke_myriad.sh          # ✓ -rwxr-xr-x (executable)
ls -d ~/Scratch/bin/duke/logs                     # ✓ Directory exists
ls ~/R/library | head -5                          # ✓ R packages installed

# Test R package installation
R --quiet --no-save << 'EOF'
.libPaths("~/R/library")
required <- c("Biostrings", "DECIPHER", "dplyr", "ggplot2", "openxlsx")
missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) == 0) {
  cat("✓ All key packages available\n")
} else {
  cat("✗ Missing packages:", paste(missing, collapse = ", "), "\n")
}
EOF
```

If all checks pass (✓), you're ready to submit jobs!

### Submit and Monitor Jobs

**Submit job:**
```bash
cd ~/Scratch/bin/duke
qsub duke_myriad.sh
```

You'll see output like:
```
Your job 1234567 ("duke_pb_HTT_6c_24h") has been submitted
```

Note the **job ID** (1234567 in this example) - you'll need it for monitoring.

**Monitor job status:**
```bash
# Check if job is running/queued
qstat -u $USER

# Example output:
# job-ID  prior   name       user   state  submit/start at     queue
# 1234567 0.50500 duke_pb_HT skgtmd r      01/06/2026 10:30:00 main.q@node-x-y-z

# State codes:
# qw = queued, waiting
# r  = running
# Eqw = error state
# (empty) = completed

# Get detailed job information
qstat -j 1234567

# Watch output in real-time (replace 1234567 with your job ID)
tail -f logs/duke_job_1234567.out

# Check for errors
tail -f logs/duke_job_1234567.err

# Check how long job has been running
qstat -u $USER | grep duke
```

**While job is running:**
```bash
# Duke also creates its own log with timestamps
# Find the latest log directory
ls -lt ~/Scratch/data/your_project/result_duke/logs/

# Watch Duke's internal log
tail -f ~/Scratch/data/your_project/result_duke/logs/20260106_103000/20260106_103000_duke_run.log
```

**After job completes (disappears from qstat):**
```bash
# Check exit status in job output
tail logs/duke_job_1234567.out
# Look for "Duke Pipeline Job Completed" message

# Check for errors
cat logs/duke_job_1234567.err

# Verify results exist
ls -lh ~/Scratch/data/your_project/result_duke/

# Check module completion (should see 7 .RData files if all modules ran)
ls -lh ~/Scratch/data/your_project/result_duke/module_data/*.RData

# Check HTML reports (7 files)
ls -lh ~/Scratch/data/your_project/result_duke/*.html

# Download results to local machine (run from LOCAL machine):
scp -r your_username@myriad.rc.ucl.ac.uk:~/Scratch/data/your_project/result_duke ./
```

**Cancel job if needed:**
```bash
qdel 1234567  # Replace with your job ID
```

### HPC Resource Guidelines

Match resources to your dataset size:

| Dataset | Samples | Reads/Sample | Walltime | CPUs | Memory | Temp |
|---------|---------|--------------|----------|------|--------|------|
| Small | 1-5 | <50K | 12h | 2 | 4GB | 20GB |
| Medium | 5-20 | 50-200K | 24h | 4 | 8GB | 50GB |
| Large | 20-50 | 200K-1M | 48h | 8 | 16GB | 100GB |
| Very Large | 50+ | >1M | 72h | 16 | 32GB | 200GB |

**Important:** Match `-pe smp N` in job script to `params$threads = N` in R script!

### HPC Troubleshooting

**Job fails immediately:**
```bash
cat logs/duke_job_<JOB_ID>.err
# Check module loading and file paths
```

**Job runs but Duke fails:**
```bash
# Check Duke log
cat logs/<TIMESTAMP>/<TIMESTAMP>_duke_run.log

# Check intermediate files
ls -lh ~/Scratch/data/your_project/result_duke/module_data/
```

**Out of memory:**
- Increase: `-l mem=16G` in job script
- Enable: `remove_intermediate = TRUE` in R script

**Out of temp space:**
- Increase: `-l tmpfs=100G` in job script

**Resume from failure:**
```r
# Set in duke_run_myriad.R
params$resume <- TRUE

# Delete problematic cache files
rm ~/Scratch/data/your_project/result_duke/temp/*.RData

# Resubmit
qsub duke_myriad.sh
```

---

### Review Results

Check `dir_out` for:
- **HTML reports**: `01_import_and_qc.html` through `07_repeat_visualisation.html`
- **Module directories**: Organised by analysis type
- **Excel files**: QC metrics, alignments, variants, range analysis
- **VCF files**: `04_allele_calling/variants/`

---

## Essential Parameters

Duke has **55 configurable parameters**. Here are the most important:

### Required Parameters (4)

| Parameter | Description | Example |
|-----------|-------------|---------|
| `dir_data` | Input directory with FASTQ/FASTA/BAM files | `/path/to/data` |
| `dir_out` | Output directory for all results | `/path/to/results` |
| `path_ref` | Reference FASTA with NNNNN separator | `/path/to/reference.fasta` |
| `path_settings` | Settings Excel/CSV (required for Module 6-7) | `/path/to/settings.xlsx` |

### Commonly Modified Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `threads` | `4` | CPU cores for parallel processing |
| `resume` | `TRUE` | Skip completed modules |
| `downsample` | `NA` | Subsample to N reads (use `NA` for all) |
| `rpt_pattern` | `"CAG"` | DNA motif to detect |
| `cluster_by` | `"repeat"` | Clustering strategy |
| `group_control` | `TRUE` | Enable control comparison |
| `waterfall` | `TRUE` | Generate waterfall plots |
| `remove_intermediate` | `TRUE` | Free RAM during run |
| `cleanup_temp` | `FALSE` | Delete temp files after |

### Reference Format ⚠️ IMPORTANT

The reference FASTA file **must** use uppercase `NNNNN` (5 consecutive Ns) to mark the repeat region location:

```
>reference_name
AGCTGATCGATCG...FLANKING_SEQUENCE_LEFT...
NNNNN
...FLANKING_SEQUENCE_RIGHT...GATCGATCGAT
```

**Example HTT reference:**
```
>HTT_CAG_repeat_amplicon
GGCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTC
NNNNN
CAACAGCCGCCACCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGC
```

---

## Settings File Format

Duke requires a settings file (Excel or CSV) defining sample-specific parameters for Module 6-7.

**Supported formats:** `.xlsx`, `.xls`, `.csv`

**Example provided:** `settings_example.xlsx`

### Required Columns

| Column | Description | Example |
|--------|-------------|---------|
| `file_name` | **Complete** FASTQ filename (Duke extracts stem automatically) | `sample_001.fastq.gz` |
| `analysis_ranges` | Range definitions using square brackets | `[0-35]` or `[0-35][36-200]` |
| `floor` | Peak detection frequency threshold(s), one per range | `[3]` or `[3][10]` |
| `max_peaks` | Maximum peaks per range | `[3]` or `[3][2]` |

### Optional Columns

| Column | Description | Example |
|--------|-------------|---------|
| `group` | Sample grouping for comparative analysis | `patient1`, `control` |
| `group_control_sample` | Any non-blank value = control sample | `1`, `TRUE`, or blank |
| `time` | Timepoint for temporal analysis | `0`, `6`, `12` |
| `manual_control_repeat_length` | Manual control setpoint(s), one per range | `[24][120]` or blank |
| `exclude` | Any non-blank value = exclude from analysis | `1`, `TRUE`, or blank |

### Understanding the `floor` Parameter

The `floor` parameter controls **frequency filtering** for peak detection. It determines the minimum number of times a repeat length must appear to be considered for analysis.

**How it works:**

- **If `floor >= 1`** (e.g., 1, 3, 10): Treated as an **absolute count** (minimum number of reads)
  - `floor = 3` means "keep only repeat lengths appearing ≥3 times"
  - Removes low-frequency noise (singletons, doubletons)
  
- **If `floor < 1`** (e.g., 0.05, 0.1): Treated as a **proportion** of maximum frequency
  - `floor = 0.05` means "keep only repeat lengths appearing ≥5% as often as the most common length"
  - Automatically scales with sequencing depth

**Common values:**

| Value | Effect | Use Case |
|-------|--------|----------|
| `[0]` | No filtering, keep all data | Maximum sensitivity, noisy |
| `[1]` | Keep all (functionally same as 0) | Minimal filtering |
| `[3]` | **Recommended default** | Removes technical noise, good balance |
| `[5]` or `[10]` | Higher stringency | Deep sequencing, focus on major peaks |
| `[0.05]` | Keep ≥5% of max frequency | Adaptive filtering, scales with depth |

**Example for diploid Huntington's disease:**
```
analysis_ranges: [0-35][36-200]
floor:           [3][10]
```
- Normal range (0-35): Lower threshold (3) to detect both alleles
- Expanded range (36-200): Higher threshold (10) as expanded alleles may be less frequent

**What gets filtered:**
- If you have a repeat length appearing only 1-2 times (likely sequencing errors)
- With `floor = [3]`, these low-frequency lengths are excluded from peak detection and statistics
- The `n_reads_total` vs `n_reads_used` columns in output show how many reads were filtered

### Analysis Ranges Syntax

Use **square bracket notation**:

**Single range:**
```
[0-35]
```

**Multiple ranges:**
```
[0-35][36-200]
```

**Open-ended range:**
```
[36-NA]    # NA = no upper limit
```

### Multi-Range Parameters

When you define **multiple analysis ranges**, these parameters must have **matching numbers of values**:

| Parameter | Format | Example for 2 ranges |
|-----------|--------|---------------------|
| `analysis_ranges` | `[min-max][min-max]` | `[0-35][36-200]` |
| `floor` | `[val1][val2]` | `[3][10]` |
| `max_peaks` | `[val1][val2]` | `[3][2]` |
| `manual_control_repeat_length` | `[val1][val2]` or blank | `[24][120]` |

**Example row:**
```
file_name: sample_001.fastq.gz
analysis_ranges: [0-35][36-200]
floor: [3][10]
max_peaks: [3][2]
group: patient1
group_control_sample: 
```

### Important Notes

- **File name matching**: Duke uses the **complete filename**, not a separate file_stem column. Duke automatically extracts the stem by removing extensions.
- **Control designation**: **Any non-blank value** in `group_control_sample` marks it as a control (`1`, `TRUE`, `yes`, etc. all work)
- **Exclusion**: **Any non-blank value** in `exclude` column excludes the sample
- **Range count consistency**: `floor` and `max_peaks` **must** have the same number of values as `analysis_ranges`

---

## Module Overview

| Module | Name | Outputs |
|--------|------|---------|
| 1 | Import & QC | Read counts, quality metrics, filtering stats |
| 2 | Alignment | Alignment to reference, MAPQ, coverage, strand |
| 3 | Repeat Detection | Repeat length calling, distributions |
| 4 | Allele Calling | Clustering, consensus sequences, variants (VCF) |
| 5 | Waterfall Plots | Per-cluster visualisations |
| 6 | Range Analysis | Modal peaks, instability metrics, control comparisons |
| 7 | Distribution Visualisation | Histograms, scatter/violin plots |

### Which Modules Use the Settings File?

**Only Module 6 (Range Analysis)** requires the settings file specified in `path_settings`.

| Module | Uses Settings? | Rerun After Settings Change? |
|--------|---------------|------------------------------|
| Module 1: Import & QC | ❌ No | ❌ No - skip with `run_module_1 = FALSE` |
| Module 2: Alignment | ❌ No | ❌ No - skip with `run_module_2 = FALSE` |
| Module 3: Repeat Detection | ❌ No | ❌ No - skip with `run_module_3 = FALSE` |
| Module 4: Allele Calling | ❌ No | ❌ No - skip with `run_module_4 = FALSE` |
| Module 5: Waterfall Plots | ❌ No | ❌ No - skip with `run_module_5 = FALSE` |
| **Module 6: Range Analysis** | **✅ Yes** | **✅ Yes - must rerun** |
| Module 7: Visualisation | ❌ No* | ✅ Yes - uses Module 6 outputs |

*Module 7 doesn't read settings directly, but visualises Module 6 outputs.

---

## Control Sample Analysis

Duke quantifies somatic instability by comparing samples to germline/control baselines.

### Setup

**1. Enable in duke_run.R:**
```r
group_control = TRUE
control_sample_selection = "flagged"  # or "earliest" or "all"
control_setpoint_metric = "median_length"  # or "modal_length" or "mean_length"
control_aggregation_method = "median"  # or "mean" or "trimmed_mean"
```

**2. Mark controls in settings file:**
Add `group_control_sample = 1` for baseline samples (blood, germline, t0).

### Outputs

**Instability metrics** (two types per sample):
- **Sample-relative:** Compare to sample's own baseline (always calculated)
- **Control-relative:** Compare to group control baseline (only if `group_control = TRUE`)
  - Includes z-scores and statistical measures
  - Detects somatic changes

**Visualisations:**
- Group control setpoints table
- Violin plots of control distributions  
- Instability index by group
- Expansion/contraction ratios

**Example HD workflow:** Blood as control (germline) → Brain samples (test somatic expansion) → Z-scores quantify deviation

---

## Analysis Ranges

Analysis ranges focus analysis on specific repeat length windows.

### Benefits

- **Biological:** Separate normal/pathogenic repeats, clinical severity categories
- **Technical:** Filter artifacts, focus on high-confidence calls  
- **Statistical:** Independent peak detection and metrics per range

### Syntax

Defined in settings file using square brackets:
```
[0-35]              # Single range
[0-35][36-200]      # Two ranges
[0-35][36-NA]       # Open-ended (NA = no upper limit)
```

### Outputs Per Range

- **Module 6:** Modal peaks, distribution metrics (mean/median/modal), spread (SD/IQR/CV), instability metrics, ridge plots
- **Module 7:** Histograms, scatter plots, combined visualisations

### Automatic "Full" Range

Duke always creates a "Full" range with all data, regardless of defined ranges. Use this to identify artifacts excluded from focused ranges.

### Best Practices

- Keep ranges consistent across samples for comparability
- Use biologically meaningful boundaries (e.g., disease thresholds)
- Typical: 2-4 ranges (more = more plots = longer runtime)

**Example HD ranges:**
```
[0-35][36-200]     # Normal vs. expanded
```

---

## Output Structure

```
dir_out/
├── 01_import_and_qc.html
├── 01_import_qc/
│   ├── plots/
│   │   ├── quality_vs_length.png                      # Combined plot
│   │   └── quality_vs_length_by_sample/               # Individual plots
│   │       └── {sample}-quality_vs_length.png
│   └── qc/
│       └── qc.xlsx
│
├── 02_alignment.html
├── 02_alignment/
│   ├── plots/
│   │   ├── alignment_length_by_flank.png              # Combined plots
│   │   ├── mapq_distribution.png
│   │   ├── repeat_segment_length.png
│   │   ├── alignment_length_by_flank_by_sample/       # Individual plots
│   │   │   └── {sample}-alignment_length_by_flank.png
│   │   ├── mapq_distribution_by_sample/
│   │   │   └── {sample}-mapq_distribution.png
│   │   ├── coverage_by_sample/
│   │   │   └── {sample}-coverage.png
│   │   └── repeat_segment_length_by_sample/
│   │       └── {sample}-repeat_segment_length.png
│   └── alignment/
│       └── alignment_qc.xlsx
│
├── 03_repeat_detection.html
├── 03_repeat_detection/
│   ├── plots/
│   │   ├── repeat_histogram.png                       # Combined plots
│   │   ├── repeat_violin.png
│   │   ├── repeat_density.png
│   │   ├── repeat_histogram_by_sample/                # Individual plots
│   │   │   └── {sample}-repeat_histogram.png
│   │   ├── repeat_violin_by_sample/
│   │   │   └── {sample}-repeat_violin.png
│   │   └── repeat_density_by_sample/
│   │       └── {sample}-repeat_density.png
│   └── repeat_detection/
│       └── repeat_summaries.xlsx
│
├── 04_allele_calling.html
├── 04_allele_calling/
│   ├── plots/
│   │   ├── scatter.png                                # Combined plots
│   │   ├── violin.png
│   │   ├── scatter_by_sample/                         # Individual plots
│   │   │   └── {sample}-scatter.png
│   │   └── violin_by_sample/
│   │       └── {sample}-violin.png
│   ├── consensus/
│   │   ├── consensus_sequences.fasta
│   │   └── consensus_summary.xlsx
│   └── variants/
│       ├── all_samples_variants.vcf
│       ├── variant_summary.xlsx
│       └── vcf_by_sample/
│           └── {sample}_variants.vcf
│
├── 05_waterfall.html
├── 05_waterfall/
│   └── plots/
│       ├── waterfall_{sample}.png
│       └── waterfall_by_cluster/
│           └── {sample}_cluster{N}-waterfall.png
│
├── 06_range_analysis.html
├── 06_range_analysis/
│   ├── plots/
│   │   ├── 01_modal_peaks.png                         # 20 main plots (01-20)
│   │   ├── 02_modal_vs_mean.png
│   │   ├── ...
│   │   ├── 19_read_proportions_by_region.png
│   │   ├── 20_expansion_contraction_balance.png
│   │   ├── 01_modal_peaks_by_sample/                  # By-sample plots
│   │   │   └── {sample}-modal_peaks.png
│   │   ├── 02_modal_vs_mean_by_range/                 # By-range plots
│   │   │   └── modal_vs_mean-{range}.png
│   │   ├── 05_tail_balance_by_range/
│   │   │   └── tail_balance-{range}.png
│   │   ├── 21_density_by_group_and_time/              # Density plots
│   │   │   ├── Full/
│   │   │   │   ├── density.png
│   │   │   │   └── density-{group}.png
│   │   │   └── {range}/
│   │   │       ├── density-{range}.png
│   │   │       └── density-{group}-{range}.png
│   │   ├── distribution_by_range/
│   │   │   └── metrics-{range}.png
│   │   ├── distribution_by_sample/
│   │   │   └── {sample}-distribution.png
│   │   ├── instability_by_range/
│   │   │   ├── instability_index-{range}.png
│   │   │   └── read_counts-{range}.png
│   │   └── instability_by_sample/
│   │       └── {sample}-instability.png
│   └── range_analysis/
│       └── range_analysis_results.xlsx
│
├── 07_repeat_visualisation.html
├── 07_repeat_visualisation/
│   └── plots/
│       ├── scatter.png                                # Combined scatter plot
│       ├── histograms_by_sample/
│       │   ├── {sample}-histogram-full.png            # Full range
│       │   └── {sample}-histogram-{range}.png         # Per range
│       ├── scatter_by_range/
│       │   └── scatter-{range}.png
│       └── scatter_by_sample/
│           └── {sample}-scatter.png
│
├── module_data/
│   └── *.RData (for resume functionality)
└── temp/
    └── (intermediate files, removed if cleanup_temp = TRUE)
```

**Plot Organization Notes:**
- Combined plots are in the root `plots/` directory
- Individual sample plots are in `*_by_sample/` subdirectories
- Individual range plots are in `*_by_range/` subdirectories
- Filenames follow pattern: `{identifier}-{plot_type}.png`
- Module 6 has 21 numbered plot groups (01-21)

---

## Complete Parameter Reference

Duke has **55 configurable parameters** organised by module. All are specified in `duke_run.R`.

### Essential Paths

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `dir_data` | string | - | ✅ | Directory containing input FASTQ/FASTA/BAM files |
| `dir_out` | string | - | ✅ | Output directory for all results |
| `path_ref` | string | - | ✅ | Reference FASTA file (must contain NNNNN separator for repeat region) |
| `path_settings` | string | - | ✅ (Module 6-7) | Settings Excel/CSV file defining analysis ranges and groupings |
| `path_trim_patterns` | string | - | Optional | CSV file with adapter/primer sequences for trimming |
| `path_manual_exclusions` | string | `NA` | Optional | Excel/CSV file listing specific read names to exclude from analysis |

### Module 1: Import & QC

**File Import:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `import_patterns` | character vector | `c("\\.fastq$", "\\.fastq.gz$", "\\.fasta$", "\\.fa$", "\\.bam$")` | File extensions to import |
| `import_recursive` | logical | `TRUE` | Search subdirectories recursively |
| `r1_pattern` | string | `"_R1_"` | Pattern to identify R1 files (paired-end) |
| `r2_pattern` | string | `"_R2_"` | Pattern to identify R2 files (paired-end) |
| `select_one_of_pair` | string | `"R1"` | For paired-end: "R1", "R2", or `NA` to keep both |
| `downsample` | numeric | `NA` | Subsample to N reads per sample |

**Adapter Trimming:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `trim_max_mismatch` | numeric | `3` | Maximum mismatches in adapter matching (0 = perfect, 3 = ~10% error) |
| `trim_with_indels` | logical | `TRUE` | Allow insertions/deletions in adapter matching |

### Module 2: Alignment

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `minimap2_args` | string | `"-t 2 -x sr -w 1"` | minimap2 arguments (flags -a, --MD, --cs added automatically) |
| `visualise_alignment` | logical | `TRUE` | Generate raw alignment plots |
| `visualise_alignment_corrected` | logical | `TRUE` | Generate plots after strand correction |

**Common minimap2 presets:**
- `"-t 2 -x sr -w 1"` - Amplicon sequencing (recommended)
- `"-t 2 -x map-ont -k 15"` - Long genomic ONT reads
- `"-t 2 -x map-hifi"` - PacBio HiFi genomic data

### Module 3: Repeat Detection

**Pattern Matching:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `rpt_pattern` | string | `"CAG"` | DNA motif to search for |
| `rpt_min_repeats` | numeric | `2` | Minimum consecutive repeats |
| `rpt_max_mismatch` | numeric | `0` | Maximum mismatches per repeat unit |
| `rpt_start_perfect_repeats` | numeric | `2` | Minimum perfect repeats at tract start |
| `rpt_end_perfect_repeats` | numeric | `2` | Minimum perfect repeats at tract end |
| `rpt_max_gap` | numeric | `6` | Maximum gap within a tract |
| `rpt_max_tract_gap` | numeric | `18` | Maximum gap between tracts to merge |
| `rpt_return_option` | string | `"longest"` | Which tract: "longest", "first", or "all" |

**Repeat Counting:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `repeat_count_method` | string | `"repeat_count_full"` | `"repeat_count_full"` (tract length ÷ pattern, recommended) or `"repeat_count_match"` (exact matches) |
| `na_repeat_handling` | string | `"convert_to_zero"` | `"convert_to_zero"`, `"filter"`, or `"flag_only"` |

### Module 4: Allele Calling & Clustering

**Clustering Strategy:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cluster` | logical | `TRUE` | Enable clustering |
| `cluster_by` | string/vector | `"repeat"` | `"repeat"`, `"haplotype"`, `c("repeat", "haplotype")`, or `"none"` |
| `repeat_cluster_max` | numeric | `20` | Maximum repeat-based clusters |
| `haplotype_cluster_max` | numeric | `10` | Maximum haplotype-based clusters |

**Haplotype Clustering:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `haplotype_region` | string | `"both"` | Flanking region: `"left"`, `"right"`, or `"both"` |
| `haplotype_method` | string | `"levenshtein"` | Distance: `"levenshtein"` or `"hamming"` |
| `haplotype_trim_length` | string/numeric | `"auto"` | `NA` (no trim), `"auto"` (modal length), or numeric (bp) |
| `haplotype_subsample` | numeric | `250` | Subsample to N reads for distance matrix |

**Consensus & Variants:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cluster_consensus` | logical | `TRUE` | Generate consensus sequences |
| `consensus_threshold` | numeric | `0.5` | Minimum base frequency (0.5 = majority) |
| `consensus_downsample` | numeric | `50` | Maximum reads per consensus |
| `call_variants` | logical | `TRUE` | Call variants (generates VCF files) |

### Module 5: Waterfall Plots

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `waterfall` | logical | `TRUE` | Generate waterfall plots |
| `waterfall_downsample` | numeric | `1000` | Maximum reads to plot per sample |
| `waterfall_rm_flank_length_outliers` | logical | `TRUE` | Remove reads with unusual flank lengths (IQR method) |
| `waterfall_y_axis_labels` | string/numeric | `"auto"` | Y-axis label density: `"auto"`, `"all"`, or integer |
| `waterfall_per_cluster` | logical | `TRUE` | Generate separate waterfall per cluster |
| `dna_colours` | named vector | See code | Hex colours for DNA bases (A, C, T, G, -, N) |

### Module 6: Range Analysis

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `path_settings` | string | - | ✅ **Required** - Settings file with analysis ranges |
| `range_peak_span` | numeric | `3` | Peak detection smoothing window (odd integer) |
| `group_control` | logical | `TRUE` | Enable control comparison analysis |
| `control_sample_selection` | string | `"flagged"` | `"flagged"`, `"earliest"`, or `"all"` |
| `control_setpoint_metric` | string | `"median_length"` | `"modal_length"`, `"mean_length"`, or `"median_length"` |
| `control_aggregation_method` | string | `"median"` | `"mean"`, `"median"`, or `"trimmed_mean"` |

### Module 7: Distribution Visualisation

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `repeat_histogram` | logical | `TRUE` | Generate frequency histograms |
| `repeat_histogram_binwidth` | numeric | `1` | Histogram bin width in repeat units |
| `repeat_scatter` | logical | `TRUE` | Generate scatter/violin plots |
| `repeat_distribution_metrics` | character vector | `c("modal_length", "mean_length", "median_length")` | Metrics to overlay on plots |

### Performance & System

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `threads` | numeric | `4` | CPU cores for parallel processing |
| `resume` | logical | `TRUE` | Skip completed modules using cached files |
| `remove_intermediate` | logical | `TRUE` | Delete large objects from memory during run (frees RAM) |
| `cleanup_temp` | logical | `FALSE` | Delete temp/ directory after completion (frees disk space) |

---

## Troubleshooting

### Issue: "Cannot find reference file"
**Solution**: Check `path_ref` is absolute path and file exists.

### Issue: "No reads pass QC"
**Solution**: Check:
- Adapter sequences in `path_trim_patterns`
- Reference sequence matches your amplicon

### Issue: "Module fails to render"
**Solution**: 
- Check all required packages installed
- Look in module HTML for specific error
- Try `resume = FALSE` to re-run from scratch

### Issue: "Out of memory"
**Solution**:
- Reduce `threads` parameter
- Enable `downsample` for testing
- Use `remove_intermediate = TRUE`
- Use `cleanup_temp = TRUE`

### Issue: "Alignments have low MAPQ"
**Solution**:
- Check reference sequence orientation
- Verify amplicon primers match reference flanks

### Issue: "Module 4 consensus generation fails silently"
**Solution**:
- This was a DECIPHER bug (fixed in v2.0.1)
- Update to latest version with fixed `lib/consensus.R`
- Bug was caused by `includeNonLetters=FALSE` parameter

### Issue: "Module 6 fails with 'no package called writexl'"
**Solution**:
- This was fixed in v2.0.1 by replacing with `openxlsx`
- Update to latest `06_range_analysis.Rmd`
- Alternative: Install writexl with `install.packages("writexl")`

### Issue: "Cannot find minimap2 on HPC"
**Solution**:
- Compile from source (see HPC Deployment section)
- Add to PATH in ~/.bashrc
- Verify with `minimap2 --version`

---

## Advanced Configuration

### Custom Read Filtering
```r
trim_max_mismatch = 3      # Adapter matching tolerance
trim_with_indels = TRUE     # Allow indels in adapters
```

### Clustering Options
```r
cluster_by = c("repeat", "haplotype")  # Two-step clustering
haplotype_cluster_max = 10              # Max haplotypes
repeat_cluster_max = 20                 # Max repeats
```

### Waterfall Customisation
```r
waterfall_downsample = 1000              # Reads per plot
waterfall_rm_flank_length_outliers = TRUE # Remove outliers
waterfall_y_axis_labels = "auto"         # Label density
waterfall_per_cluster = TRUE             # Per-cluster plots
```

---

## Package Dependencies

Duke requires **37 R packages** (27 CRAN + 5 Bioconductor + 5 base R).

**Note:** Module 6 previously used `writexl` but now uses `openxlsx` for better HPC compatibility.

### CRAN Packages (27)

**Core Data Manipulation:**
- `tidyverse`, `data.table`, `dplyr`, `tidyr`, `stringr`, `tibble`, `plyr`, `readxl`

**Visualisation:**
- `ggplot2`, `ggrepel`, `ggridges`, `ggnewscale`, `scales`, `RColorBrewer`, `viridisLite`, `cowplot`

**Reporting:**
- `rmarkdown`, `knitr`, `kableExtra`, `DT`, `htmltools`, `openxlsx`

**Statistical Analysis:**
- `mclust`, `cluster`, `ineq`, `moments`, `pracma`

**Parallel Processing:**
- `pbapply`, `pbmcapply`

**Utilities:**
- `rstudioapi`

### Bioconductor Packages (5)

- `Biostrings` - DNA/RNA sequence manipulation
- `ShortRead` - FASTQ file processing
- `GenomicAlignments` - Genomic alignment manipulation
- `Rsamtools` - BAM file manipulation
- `DECIPHER` - Multiple sequence alignment

### External Dependencies

**Required system tools:**
- `minimap2` - Read alignment (must be in PATH)
- `samtools` - BAM file processing (v1.11+)

**Installation notes:**
- On HPC systems without minimap2 modules, compile from source (see HPC Deployment section)
- samtools typically available as module on HPC systems

### Verification Script

```r
required_cran <- c("tidyverse", "data.table", "dplyr", "tidyr", "stringr", 
                   "tibble", "readxl", "plyr", "ggplot2", "ggrepel", 
                   "ggridges", "ggnewscale", "scales", "RColorBrewer", 
                   "viridisLite", "cowplot", "rmarkdown", "knitr", 
                   "kableExtra", "DT", "htmltools", "openxlsx", "mclust", 
                   "cluster", "ineq", "moments", "pracma", "pbapply", 
                   "pbmcapply", "rstudioapi")

required_bioc <- c("Biostrings", "ShortRead", "GenomicAlignments", 
                   "Rsamtools", "DECIPHER")

# Check CRAN packages
missing_cran <- required_cran[!sapply(required_cran, requireNamespace, quietly = TRUE)]
if (length(missing_cran) > 0) {
  cat("Missing CRAN packages:\n")
  print(missing_cran)
} else {
  cat("All CRAN packages installed!\n")
}

# Check Bioconductor packages
missing_bioc <- required_bioc[!sapply(required_bioc, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) {
  cat("Missing Bioconductor packages:\n")
  print(missing_bioc)
} else {
  cat("All Bioconductor packages installed!\n")
}
```

---

## Citation

If you use Duke Pipeline in your research, please cite:

```
[Your paper citation here]
```

---

## License

[Your license here]

---

## Changelog

### Version 2.0.1 (2025-01-06)

**Bug Fixes:**
- ✅ **Module 4 (consensus.R)**: Fixed DECIPHER bug where `includeNonLetters=FALSE` caused `ConsensusSequence()` to fail silently. Removed problematic parameter (lines 26-33).
- ✅ **Module 6 (06_range_analysis.Rmd)**: Replaced `writexl` package with `openxlsx` for broader HPC compatibility (line 23, line 2052).

**Enhancements:**
- ✅ Added comprehensive HPC deployment documentation
- ✅ Documented minimap2 installation for systems without module
- ✅ Added R dependency installation guide for HPC
- ✅ Created job submission and monitoring workflows
- ✅ Added resource allocation guidelines by dataset size

**Files Modified:**
- `lib/consensus.R`: Removed `includeNonLetters` parameter
- `06_range_analysis.Rmd`: Changed `writexl::write_xlsx()` to `openxlsx::write.xlsx()`
- `README.md`: Added HPC deployment section
- `duke_myriad.sh`: Template job script for UCL Myriad

### Version 2.0 (2024)
- ✅ Modular architecture (7 modules)
- ✅ Descriptive directory naming
- ✅ Enhanced range analysis (Module 6)
- ✅ Distribution visualisation (Module 7)
- ✅ Flexible metric selection
- ✅ Automatic temp cleanup option
- ✅ Resume capability
- ✅ Renamed results/ to module_data/

### Version 1.0
- Initial monolithic script
