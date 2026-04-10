# Duke pipeline

A modular pipeline for amplicon sequencing analysis with comprehensive repeat length characterisation and instability metrics.

**Latest Updates (v2.2.0):**
- ✨ **New parameters** — `check_duplicate_readnames`, `rm_flank_length_outliers`, `flank_iqr_multiplier`, `plot_dpi`, `plot_per_sample`, `export_read_counts`
- 🐛 **Bug fixes** — partial module runs no longer crash summary; reference N-masking now case-insensitive; clustering connection leak fixed
- 🔧 **Simplified** — `log_dir` and `verbose` removed (log always written to `logs/`; all messages always printed)
- 📄 **Streamlined job scripts** — HPC scripts reduced to essentials; README is now the single reference for resource guidance

---

## Quick start

### Three ways to run Duke

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

## Table of contents

- [Installation](#installation)
- [Input file requirements](#input-file-requirements)
- [Run planning](#run-planning)
- [File structure](#file-structure)
- [Run analysis](#run-analysis)
  - [Command-line interface](#1-command-line-interface-recommended)
  - [Script-based](#2-script-based-duke_runr)
  - [Interactive RStudio](#3-interactive-rstudio)
- [HPC deployment](#hpc-deployment)
  - [Job submission](#job-submission)
  - [Monitoring jobs](#monitoring-jobs)
  - [Myriad resource recommendations](#myriad-resource-recommendations)
  - [Kathleen resource recommendations](#kathleen-resource-recommendations)
  - [Resource summary](#resource-summary)
- [Parameters](#parameters)
- [Module overview](#module-overview)
- [Common workflows](#common-workflows)
- [Output structure](#output-structure)
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

## Input file requirements

Duke requires sequencing data and a reference sequence. Additional files are needed depending on which features you use.

### Required files

#### 1. Sequencing data

**Location:** Specified by `--dir_data` parameter  
**Formats:** FASTQ, FASTQ.gz, FASTA, FA, BAM  
**Organisation:** Files can be in subdirectories (recursive search enabled by default)

```bash
# Example data directory structure:
data/
├── sample_001.bam
├── sample_002.bam
├── batch_01/
│   ├── sample_003.fastq.gz
│   └── sample_004.fastq.gz
└── batch_02/
    └── sample_005.bam
```

#### 2. Reference sequence

**Location:** Specified by `--path_ref` parameter  
**Format:** FASTA file  
**Requirements:**
- Must contain your target amplicon sequence
- Include sufficient flanking sequence around repeat region
- Use `NNNNN` or `nnnnn` (upper or lowercase) as the repeat region separator

**Demo Reference:** `www/HTTset20.fasta` - HTT exon 1 reference for CAG repeat analysis

**Example HTTset20.fasta structure:**
```fasta
>HTTset20
...CGGTGCTGAGCGGCGCCGCGAGTCGGCCCGAGGCCTCCGGGGACTGCCGTGCCGGGCGGGAGACCGCC
ATGGCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTCNNNNNCAACAGCCGCCACC
GCCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGCCGCC...
```

The reference contains:
- Full HTT exon 1 sequence with flanking regions
- `nnnnn` (5+ consecutive N/n, upper or lowercase) replaces the CAG repeat tract
- This masks the repeat region, allowing Duke to detect repeats de novo from sequencing data
- Flanking sequences provide alignment context on both sides of the repeat

---

### Optional files

#### 3. Adapter patterns (required if trimming enabled)

**Location:** Specified by `--path_trim_patterns` parameter  
**Format:** CSV with columns: `adapter_name`, `adapter_sequence`  
**When needed:** Only required if `--trim TRUE` (default behaviour)

**Demo Adapters:** `www/adapters.csv` - Common Illumina, PacBio, and Nextera adapters

**Example adapters.csv:**
```csv
adapter_name,adapter_sequence
illumina_universal_read1,AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
illumina_indexed_read2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
illumina_p5,AATGATACGGCGACCACCGAGATCTACAC
illumina_p7,CAAGCAGAAGACGGCATACGAGAT
nextera_transposase,CTGTCTCTTATACACATCT
truseq_universal,AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
pacbio_smrtbell,ACGAGCGTAGCGTG
ampliconez_f,ACACTCTTTCCCTACACGACGCTCTTCCGATC
ampliconez_r,GACTGGAGTTCAGACGTGTGCTCTTCCGATCT
```

**To disable trimming entirely:**
```bash
./duke --trim FALSE ...
```

---

#### 4. Settings spreadsheet (required for Module 6)

**Location:** Specified by `--path_settings` parameter  
**Format:** Excel (.xlsx) or CSV  
**When needed:** Required for Module 6 (Range Analysis). Can skip Module 6 if not provided.

**Demo Settings:** `www/settings_example.xlsx` - Example with HTT CAG analysis ranges

##### Settings File Format

The settings file controls range-based peak detection and instability analysis. Each row represents one sample.

**Required Columns:**

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `file_name` | Text | Sample filename (must exactly match input files) | `m84277_251024_180357_s3.hifi_reads.bc2006.demux.bc1002--bc1050.bam` |
| `analysis_ranges` | Bracketed | Range definitions for peak detection | `[0-35][36-NA]` |
| `floor` | Bracketed | Minimum frequency threshold per range | `[3][3]` |
| `max_peaks` | Bracketed | Maximum peaks to detect per range | `[2][1]` |

**Optional Columns:**

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `manual_control_repeat_length` | Bracketed | Known repeat lengths for validation | `[14][100]` |
| `group` | Text | Sample grouping for comparisons | `patient_A`, `control`, `treated` |
| `group_control_sample` | Logical | Flag this sample as group control | `TRUE` or `FALSE` |
| `time` | Numeric | Timepoint for longitudinal analysis | `0`, `6`, `12`, `24` |
| `exclude` | Logical | Exclude sample from analysis | `TRUE` or `FALSE` |

##### Understanding Bracketed Parameters

Bracketed parameters use square brackets to specify values for each range:

**Format:** `[value1][value2][value3]`

- Each bracketed value corresponds to one range in `analysis_ranges`
- Number of bracketed values must match across all bracketed columns
- Use `NA` for unbounded upper limit in ranges

**Examples:**

```
# Two ranges: normal (0-35) and expanded (36+)
analysis_ranges: [0-35][36-NA]
floor:           [3][3]          # Min 3 reads in each range
max_peaks:       [2][1]          # Find up to 2 peaks in normal, 1 in expanded

# Three ranges: contracted, normal, expanded
analysis_ranges: [0-13][14-35][36-NA]
floor:           [5][3][3]
max_peaks:       [1][2][1]
manual_control:  [10][17][NA]   # Expected: 10 in contracted, 17 in normal
```

##### Complete Example Settings File

**CSV format:**
```csv
file_name,analysis_ranges,floor,max_peaks,manual_control_repeat_length,group,group_control_sample,time,exclude
sample_001.bam,[0-35][36-NA],[3][3],[2][1],[14][100],patient_A,TRUE,0,FALSE
sample_002.bam,[0-35][36-NA],[3][3],[2][1],[14][100],patient_A,FALSE,6,FALSE
sample_003.bam,[0-35][36-NA],[3][3],[2][1],[14][100],patient_A,FALSE,12,FALSE
sample_004.bam,[0-35][36-NA],[3][3],[2][1],[14][100],control,TRUE,0,FALSE
sample_005.bam,[0-35][36-NA],[3][3],[2][1],[14][100],control,FALSE,6,FALSE
```

**Excel format:**
| file_name | analysis_ranges | floor | max_peaks | manual_control_repeat_length | group | group_control_sample | time | exclude |
|-----------|----------------|-------|-----------|------------------------------|-------|---------------------|------|---------|
| sample_001.bam | [0-35][36-NA] | [3][3] | [2][1] | [14][100] | patient_A | TRUE | 0 | FALSE |
| sample_002.bam | [0-35][36-NA] | [3][3] | [2][1] | [14][100] | patient_A | FALSE | 6 | FALSE |

##### Creating Your Settings File

1. **Copy the template:**
   ```bash
   cp www/settings_example.xlsx my_settings.xlsx
   ```

2. **Edit file_name column:**
   - Must exactly match your input filenames
   - Include file extension (.bam, .fastq.gz, etc.)
   - No path, just the filename

3. **Define analysis ranges:**
   - Based on your expected repeat distribution
   - For HTT: `[0-35][36-NA]` separates normal from expanded
   - For other repeats: adjust based on biology

4. **Set detection thresholds:**
   - `floor`: Minimum reads to consider a peak (3-5 typical)
   - `max_peaks`: How many peaks to detect per range (1-2 typical)

5. **Add grouping (optional but recommended):**
   - `group`: For comparing patient cohorts, treatments, etc.
   - `group_control_sample`: TRUE for baseline/control samples
   - `time`: Numeric timepoint for longitudinal tracking

6. **Provide to Duke:**
   ```bash
   ./duke --path_settings my_settings.xlsx ...
   ```

##### Skipping Module 6

If you don't have a settings file ready:

```bash
# Run without Module 6 (Range Analysis)
./duke --run_modules 1,2,3,4,5,7 ...

# Or run all modules except 6
./duke --run_modules 1,2,3,4,5,7 --path_settings NA ...
```

You can create the settings file later and re-run just Module 6:

```bash
# Re-run only Module 6 with settings
./duke --run_modules 6 --path_settings my_settings.xlsx --resume TRUE ...
```

---

#### 5. Manual read exclusions (optional)

**Location:** Specified by `--path_manual_exclusions` parameter  
**Format:** Excel (.xlsx) or CSV  
**When needed:** Optional. Use to manually exclude problematic reads identified during QC.

**Demo Exclusions:** `www/read_exclusions_example.xlsx` - Example exclusion list format

##### Exclusions File Format

The exclusions file allows you to remove specific reads from analysis after manual inspection.

**Required Columns:**

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `file_name` | Text | Sample filename (without path) | `sample_001` |
| `read_name` | Text | Exact read identifier from BAM/FASTQ | `m64011_190830_220126/4194640/ccs` |
| `reason` | Text | Explanation (for documentation) | `Low quality alignment` |

**Sheet Name:** If using Excel, sheet must be named **"Exclusions"**

##### Complete Example Exclusions File

**CSV format:**
```csv
file_name,read_name,reason
sample_001,m64011_190830_220126/4194640/ccs,Low quality alignment
sample_001,m64011_190830_220126/4194641/ccs,Outlier repeat length
sample_002,m64011_190830_220126/5123456/ccs,Suspected artifact
sample_003,m64011_190830_220126/7891234/ccs,Chimeric read
sample_003,m64011_190830_220126/7891235/ccs,Off-target amplification
```

**Excel format:**
| file_name | read_name | reason |
|-----------|-----------|--------|
| sample_001 | m64011_190830_220126/4194640/ccs | Low quality alignment |
| sample_001 | m64011_190830_220126/4194641/ccs | Outlier repeat length |

##### Common Exclusion Reasons

- Low quality alignment
- Outlier repeat length (suspected sequencing error)
- Suspected artifact (non-biological signal)
- Chimeric read (multiple templates)
- Off-target amplification
- Manual inspection flagged
- Failed visual QC

##### Workflow for Creating Exclusions

1. **Run initial analysis without exclusions:**
   ```bash
   ./duke --path_manual_exclusions NA ...
   ```

2. **Review QC outputs:**
   - Module 2: Alignment plots (`02_alignment/plots/`)
   - Module 5: Waterfall plots (`05_waterfall/`)
   - Look for obvious outliers or problematic reads

3. **Extract read names from BAM files:**
   ```bash
   # View read names in a BAM file
   samtools view sample_001.bam | cut -f1 | head
   
   # Example output:
   # m64011_190830_220126/4194640/ccs
   # m64011_190830_220126/4194641/ccs
   # m64011_190830_220126/4194642/ccs
   ```

4. **Create exclusions file:**
   ```bash
   # Copy template
   cp www/read_exclusions_example.xlsx my_exclusions.xlsx
   
   # Edit with identified reads
   # - file_name: Must match input (without .bam/.fastq extension)
   # - read_name: Copy exactly from samtools output
   # - reason: Document why (for your records)
   ```

5. **Re-run analysis with exclusions:**
   ```bash
   ./duke --path_manual_exclusions my_exclusions.xlsx --resume TRUE ...
   ```

##### Important Notes

- **Read names must match exactly:** Include all slashes, numbers, and suffixes
- **file_name matching:** 
  - For BAM: `sample_001.bam` → use `sample_001` in exclusions file
  - For FASTQ: `sample_001.fastq.gz` → use `sample_001` in exclusions file
- **When exclusions are applied:** After alignment (Module 2), before repeat detection
- **Effect on downstream modules:** Excluded reads are removed from all subsequent analysis
- **Resume behaviour:** Exclusions require re-running from Module 2 onwards

##### Verifying Exclusions

After running with exclusions, check Module 2 output:

```
Manual exclusion complete!
Alignments excluded: 5
```

If no reads are excluded, check:
1. Read names match exactly (use `samtools view` to verify)
2. File names match (without path, with/without extension as needed)
3. Excel sheet is named "Exclusions" (case-sensitive)

---

## Run planning

Empirical data from multiple Revio flow cells can help guide experimental design for HTT amplicon sequencing.

### Expected yield per flow cell

Based on 4 complete flow cells processed through Duke:

| Samples per flow cell | Total reads | Reads per sample |
|-----------------------|-------------|------------------|
| 250 | ~7.2 million | ~28,700 |
| 287 | ~8.7 million | ~30,200 |
| 288 | ~7.5 million | ~26,200 |
| 381 | ~8.5 million | ~22,400 |

### Key observations

- Each Revio flow cell yields approximately **7–9 million aligned reads** for HTT amplicons
- Read depth per sample scales inversely with sample count
- At ~250 samples, expect **~27,000–30,000 reads per sample**
- At ~380 samples, expect **~22,000 reads per sample**
- Even at high multiplexing, coverage remains excellent for repeat length calling

### Recommendations

For most HTT CAG repeat analyses, **>5,000 reads per sample** provides robust allele calling and instability metrics. This means a single Revio flow cell can accommodate **250–400 samples** whilst maintaining coverage.

**Planning guidance:**

| Project size | Flow cells | Expected reads/sample |
|--------------|------------|----------------------|
| ≤300 samples | 1 | ~25,000–30,000 |
| 300–400 samples | 1 | ~20,000–25,000 |
| >400 samples | Consider splitting across 2 | ~25,000+ |

---

## File structure

```
duke/
├── duke                          # CLI wrapper executable
├── README.md                     # This file
├── lib/                          # Function libraries
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
├── www/                          # Demo and example files
│   ├── HTTset20.fasta            # Demo HTT reference
│   ├── adapters.csv              # Demo adapter patterns
│   ├── read_exclusions_example.xlsx
│   └── settings_example.xlsx
├── logs/                         # Job output logs
├── demo/                         # Demo datasets
└── archive/                      # Old versions
```

**Note:** Library files prefixed with module numbers for easy identification (00-07)

---

## Run analysis

### 1. Command-line interface (recommended)

#### Quick start - minimal example

```bash
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/output \
  --path_ref /path/to/reference.fasta
```

**Note:** Module 6 (Range Analysis) requires a settings spreadsheet. Without `--path_settings`, Module 6 will fail with a clear error. Either add `--path_settings /path/to/settings.xlsx` to include it, or explicitly skip it with `--run_modules 1,2,3,4,5,7`. See [Settings Spreadsheet](#4-settings-spreadsheet-required-for-module-6) for format details.

#### Local test run example
```bash
# Navigate to Duke directory
cd ~/Scratch/bin/duke

# Run Duke with demo reference and adapters
./duke \
  --dir_data demo/2025.12.17_pb_test/data \
  --dir_out demo/2025.12.17_pb_test/result_duke \
  --path_ref www/HTTset20.fasta \
  --path_trim_patterns www/adapters.csv \
  --threads 3 \
  --resume TRUE \
  --remove_intermediate TRUE \
  --cleanup_temp FALSE
```

#### With adapter trimming
```bash
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/output \
  --path_ref /path/to/reference.fasta \
  --path_trim_patterns /path/to/adapters.csv
```

#### Without adapter trimming
```bash
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/output \
  --path_ref /path/to/reference.fasta \
  --trim FALSE
```

#### See all options
```bash
./duke --help
```

**Key Features:**
- ✅ 54+ configurable parameters
- ✅ Resume enabled by default (skip completed modules)
- ✅ Comprehensive help with examples
- ✅ No file editing required

---

### 2. Script-based (duke_run.R)

**Local:**
```r
# Edit scripts/duke_run.R
params$dir_data <- "demo/2025.12.17_pb_test/data"
params$dir_out <- "demo/2025.12.17_pb_test/result_duke"
params$path_ref <- "www/HTTset20.fasta"
params$path_trim_patterns <- "www/adapters.csv"
params$trim <- TRUE

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

## HPC deployment

### Cluster selection

| Cluster | Best For | Cores | Scheduler | Script |
|---------|----------|-------|-----------|--------|
| **Myriad** | ≤500 samples | 1-36 | SGE (`-pe smp`) | `duke_myriad.sh` |
| **Old Kathleen** | 1000+ samples | 80-160 | SGE (`-pe mpi`) | `duke_kathleen.sh` |
| **New Kathleen** | 1000+ samples | 80-160 | Slurm | `duke_kathleen_slurm.sh` |

---

### Job submission

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

### Monitoring jobs

#### Check job status

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

#### Monitor log files

```bash
# Watch job output in real-time
tail -f logs/duke_*.out

# Check for errors
grep -i error logs/duke_*.err

# View completed log
less logs/duke_*.out
```

#### Cancel jobs

**SGE:**
```bash
qdel <JOB_ID>
```

**Slurm:**
```bash
scancel <JOB_ID>
```

---

### Myriad resource recommendations

Based on benchmark testing (January 2026) with datasets of 83-310 samples.

**Recommended configuration:**
```bash
#$ -pe smp 12
#$ -l mem=64G        # 64GB per core = 768GB total
#$ -l tmpfs=500G
#$ -l h_rt=48:00:00
```

With Duke command:
```bash
./duke --threads 12 ...
```

#### Expected runtimes (Myriad)

| Samples | Runtime | Example |
|---------|---------|---------|
| ~100 | 8-10h | jasmine (102 samples): 8.1h |
| ~300 | ~10h | pg (298 samples): 9.4h |
| 300+ | 28-35h | fs (310 samples): 28.3h |

#### Myriad benchmark results

Testing compared different core/memory configurations:

| Dataset | Samples | 12c/64G | 6c/128G | 24c/32G |
|---------|---------|---------|---------|---------|
| jasmine | 102 | **8.1h** | 8.9h | 13.0h |
| pg | 298 | **9.4h** | 19.9h | 17.5h |
| fs | 310 | 28-34h | **27.1h** | Failed |
| yas | 83 | 9.6h | 15.1h | **5.7h** |

#### Key findings

1. **12 cores with 64GB/core is optimal for most datasets** - provides the best balance of parallelisation and memory

2. **More RAM per core doesn't always help** - 6c/128G was slower than 12c/64G for the pg dataset (19.9h vs 9.4h) despite having more memory per core

3. **Too many cores with less memory can fail** - 24c/32G failed to complete the fs dataset within 48h

4. **Very large datasets (300+ samples) are slow regardless** - expect 24-48h; consider Kathleen for these

#### Configurations to avoid (Myriad)

```bash
# Not recommended - too few cores despite high memory
#$ -pe smp 6
#$ -l mem=128G       # Slower for most datasets

# Not recommended - can cause failures on large datasets  
#$ -pe smp 24
#$ -l mem=32G        # fs dataset failed with this config
```

---

### Kathleen resource recommendations

For datasets over 500 samples, use Kathleen with 80+ cores.

**Old Kathleen (SGE):**
```bash
#$ -pe mpi 80
#$ -l mem=4G         # 320GB total
#$ -l h_rt=48:00:00
# Use --threads 80 in Duke
```

**New Kathleen (Slurm):**
```bash
#SBATCH --nodes=2
#SBATCH --ntasks=80
#SBATCH --mem-per-cpu=4G
#SBATCH --time=48:00:00
# Use --threads 80 in Duke
```

**Memory workaround for Kathleen** (if mem=4G+ unavailable):
```bash
# Request double cores, use half for processing
#$ -pe mpi 160
#$ -l mem=2G        # 320GB total
#$ -l h_rt=48:00:00
# Use --threads 80 in Duke (not 160!)
```

---

### Resource summary

| Dataset size | Samples | Cluster | Cores | Memory/core | tmpfs | Runtime |
|--------------|---------|---------|-------|-------------|-------|---------|
| Small | <100 | Myriad | 12 | 64GB | 500GB | 8-10h |
| Medium | 100-300 | Myriad | 12 | 64GB | 500GB | 8-10h |
| Large | 300-500 | Myriad | 12 | 64GB | 500GB | 28-48h |
| Very large | 500-1000 | Kathleen | 80 | 4GB | - | 16-24h |
| Massive | 1000+ | Kathleen | 80-160 | 4GB | - | 24-48h |

**Key Notes:**
- Myriad maximum: 36 cores (but 12 cores is optimal based on benchmarks)
- Kathleen: Request cores in multiples of 40
- Match `--threads` parameter to core count (except memory workaround above)

---

### SGE vs Slurm command reference

| Action | SGE (Myriad / Old Kathleen) | Slurm (New Kathleen) |
|--------|----------------------------|----------------------|
| Make executable | `chmod +x script.sh` | Not needed |
| Submit job | `qsub script.sh` | `sbatch script.sh` |
| Check jobs | `qstat -u $USER` | `squeue -u $USER` |
| Job details | `qstat -j <JOB_ID>` | `scontrol show job <JOB_ID>` |
| Cancel job | `qdel <JOB_ID>` | `scancel <JOB_ID>` |
| After completion | `qacct -j <JOB_ID>` | `sacct -j <JOB_ID>` |
| Job ID variable | `$JOB_ID` | `$SLURM_JOB_ID` |

---

## Parameters

### New parameters in v2.2.0

- `check_duplicate_readnames` — scan for reads with the same name and keep only the longest copy (Module 1)
- `rm_flank_length_outliers` — filter reads with outlier combined flank lengths in Module 3; propagates to all downstream modules
- `flank_iqr_multiplier` — IQR multiplier for flank length outlier detection (only active when `rm_flank_length_outliers = TRUE`)
- `plot_dpi` — resolution for all `ggsave` diagnostic output (150 = screen quality, 300 = print quality)
- `plot_per_sample` — gate on/off per-sample plot generation in Modules 2, 3, and 4
- `export_read_counts` — write per-read repeat counts to `.tsv.gz` files in `03_repeat_detection/read_counts/` (default `TRUE`)

**Removed in v2.2.0:** `log_dir` (hard-coded to `logs/`), `verbose` (all messages always printed), `waterfall_per_sample`, `repeat_histogram_per_sample`, `repeat_scatter_per_sample` (unimplemented).

### Essential parameters

```r
# Required
dir_data = "/path/to/data"            # Input directory
dir_out = "/path/to/output"           # Output directory
path_ref = "/path/to/reference.fasta" # Reference FASTA (NNNNN or nnnnn separator)

# Trimming (if trim = TRUE)
path_trim_patterns = "/path/to/adapters.csv"

# Runtime
threads = 2                           # CPU cores
resume = TRUE                         # Skip completed modules (default)
```

### Complete parameter reference

See `./duke --help` for comprehensive documentation of all 54+ parameters, including:

#### File paths
- `dir_data` - Input directory **(required)**
- `dir_out` - Output directory **(required)**
- `path_ref` - Reference FASTA **(required)** - See `www/HTTset20.fasta` for demo
- `path_trim_patterns` - Adapter file (required if `trim=TRUE`) - See `www/adapters.csv` for demo
- `path_settings` - Settings Excel (required for Module 6) - See [Settings File Format](#settings-file-format)
- `path_manual_exclusions` - Read exclusions Excel (optional) - See [Exclusions File Format](#exclusions-file-format)

#### Import options
- `import_patterns` - File extensions to import
- `downsample` - Limit reads per sample (NA = all)
- `select_one_of_pair` - For paired-end: "R1", "R2", or NA
- `check_duplicate_readnames` - Remove duplicate read names, keeping longest copy (default: TRUE; safe to set FALSE for PacBio CCS)

#### Adapter trimming
- `trim` - Enable/disable (default: TRUE)
- `trim_max_mismatch` - Max errors (default: 3)
- `trim_with_indels` - Allow indels (default: TRUE)

#### Alignment
- `minimap2_args` - minimap2 settings (default: "-t 2 -x sr -w 1")
- `visualise_alignment` - Generate plots (default: TRUE)
- `visualise_alignment_downsample` - Max reads to plot (default: 1000)

#### Repeat detection
All parameters have comprehensive documentation with examples in `scripts/duke_run.R` or via `./duke --help`:

- `rpt_pattern` - Repeat motif (default: "CAG")
- `rpt_min_repeats` - Minimum count (default: 2)
- `rpt_max_mismatch` - Error tolerance (default: 0)
- `rpt_start_perfect_repeats` - Perfect at start (default: 2)
- `rpt_end_perfect_repeats` - Perfect at end (default: 2)
- `rpt_max_gap` - Max gap within tract (default: 6)
- `rpt_max_tract_gap` - Max gap between tracts (default: 18)
- `rpt_return_option` - "longest" or "all"

#### Flank length QC (Module 3)
- `rm_flank_length_outliers` - Filter reads with outlier combined flank lengths (default: FALSE — review flank plots first)
- `flank_iqr_multiplier` - IQR multiplier for outlier detection (default: 1.5; only active when `rm_flank_length_outliers = TRUE`)

#### Clustering
- `cluster` - Enable clustering (default: TRUE)
- `cluster_by` - Method: "repeat", "haplotype", or both
- `haplotype_cluster_max` - Max haplotypes (default: 10)
- `repeat_cluster_max` - Max repeat lengths (default: 20)

#### Plot output
- `plot_dpi` - Resolution for all diagnostic plots (default: 150; use 300 for publication figures)
- `plot_per_sample` - Generate per-sample plot files in Modules 2, 3, 4 (default: TRUE; consider FALSE for >50 samples)

#### Data export
- `export_read_counts` - Export per-read repeat counts to `03_repeat_detection/read_counts/{sample}.tsv.gz` (default: TRUE). Each file contains `read_id`, `repeat_count_full`, `repeat_count_match`, `repeat_count_tracts`. Useful for loading data in Python, Excel, or other downstream tools.

#### Module control
- `waterfall` - Generate waterfall plots (default: TRUE)
- `repeat_histogram` - Generate histograms (default: TRUE)
- `repeat_scatter` - Generate scatter plots (default: TRUE)

#### Runtime
- `threads` - CPU cores (default: 80 HPC / 2 local)
- `resume` - Skip completed modules (default: TRUE)
- `remove_intermediate` - Free RAM between modules (default: TRUE)
- `cleanup_temp` - Delete temp files on completion (default: FALSE)
- `run_modules` - Which modules to run (default: c(1:7))

---

## Module overview

### Module 1: Import and QC
- Imports FASTQ/FASTA/BAM files
- **Optionally trims adapters** (controlled by `trim` parameter)
- **Optionally checks and removes duplicate read names**, keeping longest copy (controlled by `check_duplicate_readnames`; default TRUE; safe to disable for PacBio CCS)
- Optional downsampling
- QC metrics and plots

### Module 2: Alignment
- minimap2 alignment to reference
- Strand correction
- Coverage visualisation
- **Optional manual read exclusion** (controlled by `path_manual_exclusions`)

### Module 3: Repeat detection
- Identifies repeat tracts
- Tolerates sequencing errors
- Multiple counting methods
- **Flank length QC** — per-sample distribution plots always generated; optional outlier filtering via `rm_flank_length_outliers` (default FALSE — review plots before enabling). Filtered reads propagate to all downstream modules.

### Module 4: Allele calling
- **⚠️ Most memory-intensive module** - requires 4GB+ per core
- Clusters by repeat/haplotype using parallel processing
- Consensus sequences
- Variant calling (VCF export)

### Module 5: Waterfall plots
- Visual read inspection
- Per-sample and per-cluster plots
- Configurable downsampling

### Module 6: Range analysis

**Requires:** Settings spreadsheet (see [Settings Spreadsheet](#4-settings-spreadsheet-required-for-module-6))

- **Range-based peak detection:** Identifies modal peaks within user-defined repeat length ranges
- **Instability metrics:** Calculates expansion index, contraction index, and somatic instability measures
- **Control comparisons:** Tests differences from control samples or baseline timepoints
- **Statistical analysis:** Per-group and per-timepoint comparisons
- **Comprehensive tables:** Exports detailed Excel summaries with range statistics
- **Flexible range definitions:** Supports multiple ranges per sample (e.g., normal vs expanded)

**Settings file controls:**
- Which ranges to analyse (e.g., `[0-35][36-NA]` for normal and expanded)
- Peak detection thresholds (`floor`, `max_peaks`)
- Sample grouping for comparative analysis (`group`, `time`)
- Control sample selection (`group_control_sample`)

**To skip Module 6:**
```bash
./duke --run_modules 1,2,3,4,5,7 ...
```

**To run Module 6 only (after settings file created):**
```bash
./duke --run_modules 6 --path_settings my_settings.xlsx --resume TRUE ...
```

### Module 7: Repeat visualisation
- Frequency histograms
- Scatter/violin plots
- Publication-ready figures

---

## Module diagnostics

Every module ends with a **`# Module diagnostics`** section in the HTML report, appearing just before Session info.

### Timing table

All seven modules record wall-clock time for each major computation section and display a summary table with four columns:

| Column | Description |
|---|---|
| Section | Name of the computation step |
| Seconds | Raw elapsed seconds for that step |
| Elapsed (mm:ss) | Step duration formatted as mm:ss or hh:mm:ss |
| Cumulative (mm:ss) | Running total from module start |

The **Total** row (bold) gives the complete module wall time. This makes it straightforward to identify which step dominates runtime and whether caching (resume) is working.

### Output manifest

Modules 1, 2, 3, 4, and 6 — those that save an `.RData` file — also display an **output manifest** table listing every object in the module's output list:

| Column | Description |
|---|---|
| Object | R object name (as accessed via `moduleN_output$...`) |
| Description | What the object contains and how it should be used downstream |
| Class | R class (data.frame, tbl_df, list, etc.) |
| In-memory size | Human-readable size (B / KB / MB / GB) |

The table caption shows the `.RData` file size on disk and full path.

### Excel export

Both tables are appended as additional sheets to each module's output Excel file (`timing` and `manifest`), so timing data is available alongside the analysis results without opening the HTML.

Modules 5 and 7 (waterfall and repeat visualisation) save only a `timing` sheet, as they produce no `.RData` output (only plots).

---

## Common workflows

### Workflow 1: Discovery mode (first run)

**Use when:** Exploring new data or validating amplicon sequencing quality

```bash
# Run without settings file to explore data
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/results \
  --path_ref www/HTTset20.fasta \
  --path_trim_patterns www/adapters.csv \
  --run_modules 1,2,3,4,5,7 \
  --threads 12

# Review outputs:
# - Module 1: QC metrics, read lengths
# - Module 2: Alignment quality, coverage
# - Module 3: Repeat detection
# - Module 4: Allele clustering
# - Module 5: Waterfall plots (visual inspection)
# - Module 7: Repeat distribution histograms

# Use results to:
# 1. Identify appropriate analysis ranges for settings file
# 2. Flag problematic reads for exclusion
# 3. Validate reference sequence and adapter trimming
```

### Workflow 2: Full analysis with range metrics

**Use when:** You know your expected ranges and need quantitative instability metrics

```bash
# 1. Create settings file based on discovery run
#    Define ranges, floors, max_peaks, groups
#    Example: [0-35][36-NA] for normal/expanded HTT CAG

# 2. Run complete analysis with settings
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/results \
  --path_ref www/HTTset20.fasta \
  --path_trim_patterns www/adapters.csv \
  --path_settings /path/to/settings.xlsx \
  --threads 12

# 3. Review Module 6 outputs:
#    - Modal peaks per range
#    - Instability metrics (expansion/contraction indices)
#    - Control comparisons
#    - Excel summary tables
```

### Workflow 3: Quality control with read exclusions

**Use when:** Manual QC reveals problematic reads that need exclusion

```bash
# 1. Initial run (discovery mode)
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/results_initial \
  --path_ref www/HTTset20.fasta

# 2. Review Module 2 and 5 plots
#    Identify outlier reads visually

# 3. Extract problematic read names
samtools view sample_001.bam | grep "suspicious_pattern" | cut -f1 > flagged_reads.txt

# 4. Create exclusions file
#    Format: file_name, read_name, reason
#    See: www/read_exclusions_example.xlsx

# 5. Re-run with exclusions
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/results_filtered \
  --path_ref www/HTTset20.fasta \
  --path_manual_exclusions /path/to/exclusions.xlsx \
  --resume FALSE  # Force full re-run

# 6. Compare filtered vs unfiltered results
#    Check Module 2 output: "Alignments excluded: X"
```

### Workflow 4: Longitudinal analysis

**Use when:** Analysing timepoint series or treatment responses

```bash
# 1. Create settings file with grouping
#    Required columns: group, time, group_control_sample
#    Example:
#      file_name,group,time,group_control_sample,analysis_ranges,floor,max_peaks
#      patient_A_t0.bam,patient_A,0,TRUE,[0-35][36-NA],[3][3],[2][1]
#      patient_A_t6.bam,patient_A,6,FALSE,[0-35][36-NA],[3][3],[2][1]
#      patient_A_t12.bam,patient_A,12,FALSE,[0-35][36-NA],[3][3],[2][1]

# 2. Run analysis
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/results \
  --path_ref www/HTTset20.fasta \
  --path_settings /path/to/settings.xlsx \
  --control_sample_selection flagged  # Use group_control_sample = TRUE

# 3. Module 6 outputs will include:
#    - Baseline (t=0) control comparisons
#    - Timepoint progression tracking
#    - Instability metrics over time
#    - Per-group statistical summaries
```

### Workflow 5: Large dataset HPC run

**Use when:** Processing hundreds of samples on computing cluster

```bash
# 1. Prepare settings file with all samples
#    One row per sample (100+ rows typical)

# 2. Edit HPC job script
nano scripts/duke_myriad.sh  # or duke_kathleen_slurm.sh

# Update paths in script:
#   --path_ref ~/Scratch/refs/HTTset20.fasta
#   --dir_data /home/user/Scratch/data/large_batch
#   --dir_out /home/user/Scratch/results/large_batch  
#   --path_settings /home/user/Scratch/settings/batch_settings.xlsx
#   --threads 36  # Match cluster allocation

# 3. Submit job
qsub scripts/duke_myriad.sh  # SGE
# or
sbatch scripts/duke_kathleen_slurm.sh  # Slurm

# 4. Monitor progress
watch -n 5 'qstat -u $USER'  # SGE
watch -n 5 'squeue -u $USER'  # Slurm

# 5. Check logs during run
tail -f logs/duke_*.out

# 6. Resume if job times out
#    Duke automatically resumes from last completed module
qsub scripts/duke_myriad.sh  # Resume parameter = TRUE by default
```

### Workflow 6: Parameter optimisation

**Use when:** Fine-tuning repeat detection or clustering parameters

```bash
# 1. Test run with downsampling
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/test_params \
  --path_ref www/HTTset20.fasta \
  --downsample 500 \
  --rpt_max_mismatch 0 \
  --threads 4

# 2. Review Module 3 repeat detection
#    Check if repeats are detected correctly

# 3. Adjust parameters and re-run
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/test_params_v2 \
  --path_ref www/HTTset20.fasta \
  --downsample 500 \
  --rpt_max_mismatch 1 \  # Allow 1 mismatch
  --rpt_max_gap 9 \        # Increase gap tolerance
  --threads 4

# 4. Compare results across parameter sets
#    Module 7 histograms are useful for comparison

# 5. Run full dataset with optimised parameters
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/results_final \
  --path_ref www/HTTset20.fasta \
  --rpt_max_mismatch 1 \
  --rpt_max_gap 9 \
  --threads 12
```

---

## Output structure

```
result_duke/
├── 01_import_and_qc.html             # HTML reports (one per module)
├── 02_alignment.html
├── 03_repeat_detection.html
├── 04_allele_calling.html
├── 05_waterfall.html
├── 06_range_analysis.html
├── 07_repeat_visualisation.html
│
├── module_data/                      # Cached RData (for resume)
│   ├── 01_import_qc_results.RData
│   ├── 02_alignment_results.RData
│   ├── 03_repeat_detection_results.RData
│   ├── 04_allele_calling_results.RData
│   └── 06_range_analysis_results.RData
│
├── 01_import_qc/
│   ├── qc.xlsx                       # count_summary | trim_summary | timing | manifest
│   └── plots/
├── 02_alignment/
│   ├── alignment.xlsx                # alignment_summary | flank_filtering | strand_summary
│   │                                 # alignment_quality | sequence_segments | timing | manifest
│   └── plots/
│       └── coverage_by_sample/
├── 03_repeat_detection/
│   ├── repeat_detection.xlsx         # count_summary | tract_detection | count_method_comparison
│   │                                 # heatmap_sample_index | flank_summary | timing | manifest
│   ├── read_counts/                  # Per-read TSV files (only if export_read_counts = TRUE)
│   │   ├── {sample_1}.tsv.gz         # read_id | repeat_count_full | repeat_count_match | repeat_count_tracts
│   │   └── {sample_N}.tsv.gz
│   └── plots/
│       └── flank/
├── 04_allele_calling/
│   ├── allele_calling.xlsx           # cluster_summary | timing | manifest
│   ├── consensus/
│   │   ├── consensus.xlsx
│   │   ├── consensus_sequences.fasta
│   │   └── vcf_by_sample/
│   └── plots/
├── 05_waterfall/
│   ├── waterfall.xlsx                # timing
│   └── plots/
│       └── waterfall_by_cluster/
├── 06_range_analysis/
│   ├── range_analysis.xlsx           # range_summary | modal_peaks | distribution_summary
│   │                                 # instability_metrics | group_controls | timing | manifest
│   └── plots/
└── 07_repeat_visualisation/
    ├── repeat_visualisation.xlsx     # timing
    └── plots/
        ├── histograms_by_sample/
        ├── scatter_by_sample/
        └── density_by_sample/
```

---

## Troubleshooting

### Common issues

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

# Or use demo reference:
./duke --path_ref www/HTTset20.fasta ...
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

**"Settings file not found" or "Settings file required for Module 6"**
```bash
# CAUSE: path_settings not provided or file doesn't exist
# SOLUTION 1: Skip Module 6 if you don't need range analysis
./duke --run_modules 1,2,3,4,5,7 ...

# SOLUTION 2: Create settings file
# 1. Copy template
cp www/settings_example.xlsx my_settings.xlsx

# 2. Edit with your samples and analysis ranges
# - file_name: Must match your input files
# - analysis_ranges: Define your ranges, e.g., [0-35][36-NA]
# - floor, max_peaks: Set thresholds

# 3. Provide to Duke
./duke --path_settings my_settings.xlsx ...
```

**"Cannot parse analysis_ranges" or "Bracketed parameter error"**
```bash
# CAUSE: Incorrect bracketed parameter format in settings file
# SOLUTION: Check format - must use [value1][value2] with NO spaces

# ✅ Correct:
#   analysis_ranges: [0-35][36-NA]
#   floor:           [3][3]
#   max_peaks:       [2][1]

# ❌ Incorrect:
#   [0-35] [36-NA]     (space between brackets)
#   [0-35,36-NA]       (comma instead of ][)
#   [0-35][36-100]     (should be [36-NA] for unbounded)
#   0-35][36-NA        (missing opening bracket)

# TIP: Count your brackets - should be equal number of [ and ]
```

**"Mismatch in bracketed parameter lengths"**
```bash
# CAUSE: Different number of ranges in bracketed parameters
# SOLUTION: Ensure all bracketed columns have the same number of values

# ✅ Correct (both have 2 ranges):
#   analysis_ranges: [0-35][36-NA]
#   floor:           [3][3]
#   max_peaks:       [2][1]

# ❌ Incorrect (mismatched counts):
#   analysis_ranges: [0-35][36-NA]      (2 ranges)
#   floor:           [3]                (1 range - WRONG!)
#   max_peaks:       [2][1]             (2 ranges)

# Each row must have consistent bracketed parameter counts
```

**"file_name not found in data" in settings file**
```bash
# CAUSE: Filename in settings doesn't match actual input files
# SOLUTION: Check exact filename match

# List your input files:
ls -1 /path/to/data

# Example output:
#   m84277_251024_180357_s3.hifi_reads.bc2006.demux.bc1002--bc1050.bam
#   m84277_251024_180357_s3.hifi_reads.bc2006.demux.bc1002--bc1051.bam

# Settings file must use EXACT filenames:
# file_name column: m84277_251024_180357_s3.hifi_reads.bc2006.demux.bc1002--bc1050.bam
#                   (not: bc1050.bam or sample_1 or any abbreviation)

# Common issues:
# - Missing file extension (.bam, .fastq.gz)
# - Extra/missing characters in filename
# - Path included when only filename needed
```

**Manual exclusions not working (reads not excluded)**
```bash
# CAUSE 1: Read names don't match BAM file format exactly
# SOLUTION: Extract exact read names from your BAM files

# Get actual read names:
samtools view input.bam | cut -f1 | head -5

# Example PacBio CCS read name:
#   m64011_190830_220126/4194640/ccs
#   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   Copy this EXACTLY into exclusions file

# Example Illumina read name:
#   M00123:456:000000000-ABCDE:1:1101:15440:1234 1:N:0:1
#   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   Copy this EXACTLY (including spaces and colons)

# CAUSE 2: file_name doesn't match
# SOLUTION: Use base filename without extension
# BAM file: sample_001.bam
# Exclusions: file_name should be "sample_001" (not "sample_001.bam")

# CAUSE 3: Excel sheet not named "Exclusions"
# SOLUTION: Rename sheet to exactly "Exclusions" (case-sensitive)

# Verify exclusions worked:
# Check Module 2 output for: "Alignments excluded: X"
```

**Want to update settings/exclusions and re-run**
```bash
# Settings file changes require re-running Module 6:
./duke --run_modules 6 --path_settings updated_settings.xlsx --resume TRUE ...

# Exclusions require re-running from Module 2 onwards:
./duke --run_modules 2,3,4,5,6,7 --path_manual_exclusions updated_exclusions.xlsx ...

# Or delete specific module data to force re-run:
rm result_duke/module_data/02_alignment_results.RData
rm result_duke/module_data/06_range_analysis_results.RData
```

---

## Version history

### v2.2.0 (April 2026) - Current
- ✨ **NEW:** Module diagnostics section in every HTML report — timing table per module (all 7) and output manifest (modules 1–4, 6); both exported to each module's Excel file
- ✨ **NEW:** `export_read_counts` parameter (Module 3) — opt-in export of per-read repeat counts to `{sample}.tsv.gz` (gzipped TSV; columns: read_id, repeat_count_full, repeat_count_match, repeat_count_tracts)
- ✨ **NEW:** `check_duplicate_readnames` parameter (Module 1; default TRUE)
- ✨ **NEW:** Flank length QC filtering — `rm_flank_length_outliers`, `flank_iqr_multiplier` (Module 3)
- ✨ **NEW:** Plot output control — `plot_dpi`, `plot_per_sample` (Modules 2, 3, 4)
- ✨ **NEW:** `format_size()`, `format_elapsed()`, `build_timing_table()` helpers in `00_utils.R`
- 🐛 **FIXED:** Modules 5 and 7 no longer have spurious resume checks — they always run when included in `run_modules` (use `run_modules` to skip them)
- 🐛 **FIXED:** Manifest in-memory sizes now use adaptive units (B/KB/MB/GB) — previously showed 0 for small objects
- 🐛 **FIXED:** RData file size shown in manifest caption only, not repeated per row
- 🐛 **FIXED:** Summary section crash when running partial module sets
- 🐛 **FIXED:** Reference N-masking now case-insensitive (`NNNNN` or `nnnnn`)
- 🐛 **FIXED:** Clustering early-exit return type inconsistency
- 🐛 **FIXED:** GMM connection leak in Module 4 clustering (replaced `sink(tempfile())` with `capture.output()`)
- 🔧 **REMOVED:** `log_dir` (hard-coded to `logs/`), `verbose`, `waterfall_per_sample`, `repeat_histogram_per_sample`, `repeat_scatter_per_sample`
- 🔧 **UPDATED:** `consensus_downsample` default 200 → 50; `waterfall_per_cluster` default TRUE → FALSE
- 📄 **STREAMLINED:** HPC job scripts; README is now the single reference for resource guidance

### v2.1.0 (January 2026)
- 🐛 **FIXED:** `knit_root_dir` handling for Rmd files in subdirectories
- 📚 **ENHANCED:** Comprehensive reorganisation documentation
- 📊 **ADDED:** Job monitoring commands in README
- 📖 **IMPROVED:** Complete input file format documentation
- 🎯 **ADDED:** Demo reference (HTTset20.fasta) and adapters (adapters.csv)
- 📋 **ADDED:** Common workflows section with six example workflows

### v2.0.1 (January 2025)
- ✨ **NEW:** Command-line interface
- ✅ **NEW:** `trim` parameter (optional adapter trimming)
- 🔧 **RENAMED:** `visualise_alignment_downsample` (was: `_n_reads`)
- 🐛 **FIXED:** CLI resume detection
- 🐛 **FIXED:** Resume default now TRUE (consistent)
- 📚 **ENHANCED:** Comprehensive repeat parameter documentation
- 📦 **ORGANISED:** Library files with module prefixes (00-07)
- 📊 **IMPROVED:** Logging format

### v2.0.0
- Modular architecture (7 modules)
- Resume capability
- Cleanup options

---

## Citation

If you use the Duke pipeline in your research, please cite:

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

**Duke - Clean, organised, and production-ready!** 🎯
