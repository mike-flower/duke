# Duke pipeline

A modular pipeline for amplicon sequencing analysis with comprehensive repeat length characterisation and instability metrics.

**v2.2.0** — [What's new](#version-history)

---

## Table of contents

- [What is Duke?](#what-is-duke)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Input files](#input-files)
- [Running Duke](#running-duke)
  - [Command-line interface](#1-command-line-interface-recommended)
  - [Script-based](#2-script-based-duke_runr)
  - [Interactive RStudio](#3-interactive-rstudio)
- [Output files](#output-files)
- [Module overview](#module-overview)
- [Parameters](#parameters)
- [Module diagnostics](#module-diagnostics)
- [HPC deployment](#hpc-deployment)
- [Run planning](#run-planning)
- [Common workflows](#common-workflows)
- [Troubleshooting](#troubleshooting)
- [Version history](#version-history)
- [Citation](#citation)
- [Contact](#contact)

---

## What is Duke?

Duke analyses amplicon sequencing data to measure repeat lengths and somatic instability. It is designed for HTT CAG repeat analysis but works with any tandem repeat locus. Given raw sequencing files (BAM or FASTQ), it produces:

- **Per-read repeat counts** — all three counting methods, exported as `.tsv.gz` for downstream use
- **Allele calls** — clustering of reads into alleles based on repeat length and/or haplotype
- **Instability metrics** — expansion and contraction indices relative to control samples
- **Visualisations** — waterfall plots, repeat histograms, scatter plots, density ridges
- **HTML reports** — one per module, with interactive tables and figures

Duke runs as seven sequential modules:

| Module | What it does |
|--------|-------------|
| 1. Import and QC | Read import, quality metrics, adapter trimming |
| 2. Alignment | minimap2 alignment to reference flanks, strand correction |
| 3. Repeat detection | Repeat tract finding and counting; flank QC |
| 4. Allele calling | Read clustering, consensus sequences, variant calling |
| 5. Waterfall plots | Per-read base composition visualisation |
| 6. Range analysis | Peak detection and instability metrics (requires settings file) |
| 7. Repeat visualisation | Histograms, scatter plots, density plots |

Resume is enabled by default — completed modules are automatically skipped on re-run.

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

**HPC:** Point R to your existing package library:
```bash
export R_LIBS_USER=~/R/library
```

**File structure:**
```
duke/
├── duke                          # CLI wrapper executable
├── README.md
├── lib/                          # Function libraries (00_utils.R through 07_visualisation.R)
├── modules/                      # Module Rmd files (01-07)
├── scripts/
│   ├── duke_run.R                # Script-based runner
│   ├── duke_cli.R                # CLI script
│   ├── duke_myriad.sh            # Myriad (SGE) job script
│   ├── duke_kathleen.sh          # Old Kathleen (SGE) job script
│   └── duke_kathleen_slurm.sh    # New Kathleen (Slurm) job script
└── www/                          # Demo files
    ├── HTTset20.fasta
    ├── adapters.csv
    ├── settings_example.xlsx
    └── read_exclusions_example.xlsx
```

---

## Quick start

**Minimal run** (no adapter trimming, no range analysis):
```bash
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/output \
  --path_ref /path/to/reference.fasta \
  --trim FALSE \
  --run_modules 1,2,3,4,5,7
```

**Full run** with trimming and range analysis:
```bash
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/output \
  --path_ref /path/to/reference.fasta \
  --path_trim_patterns /path/to/adapters.csv \
  --path_settings /path/to/settings.xlsx \
  --threads 12
```

**Demo run** using the included test data:
```bash
cd ~/Scratch/bin/duke

./duke \
  --dir_data demo/data \
  --dir_out demo/result_duke \
  --path_ref www/HTTset20.fasta \
  --path_trim_patterns www/adapters.csv \
  --threads 3
```

**See all options:**
```bash
./duke --help
```

---

## Input files

### Required

#### 1. Sequencing data

**Parameter:** `--dir_data`  
**Formats:** BAM, FASTQ, FASTQ.gz, FASTA  
**Note:** Subdirectories are searched recursively by default

#### 2. Reference sequence

**Parameter:** `--path_ref`  
**Format:** FASTA with `NNNNN` (or `nnnnn`) marking the repeat region

```fasta
>HTTset20
...CGGTGCTGAGCGGCGCCGCGAGTCGGCCCGAGGCCTCCGGGGACTGCCGTGCCGGGCGGGAGACCGCC
ATGGCGACCCTGGAAAAGCTGATGAAGGCCTTCGAGTCCCTCAAGTCCTTCNNNNNCAACAGCCGCCACC
GCCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGCCGCC...
```

The 5+ consecutive N/n characters mark the repeat tract. Duke detects repeats de novo from the sequencing data using the flanking sequence for alignment.

**Demo reference:** `www/HTTset20.fasta` — HTT exon 1 for CAG repeat analysis

---

### Optional

#### 3. Adapter patterns (required if `--trim TRUE`)

**Parameter:** `--path_trim_patterns`  
**Format:** CSV with columns `adapter_name` and `adapter_sequence`

```csv
adapter_name,adapter_sequence
illumina_universal_read1,AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
nextera_transposase,CTGTCTCTTATACACATCT
pacbio_smrtbell,ACGAGCGTAGCGTG
```

**Demo adapters:** `www/adapters.csv`

To disable trimming entirely:
```bash
./duke --trim FALSE ...
```

---

#### 4. Settings spreadsheet (required for Module 6)

**Parameter:** `--path_settings`  
**Format:** Excel (.xlsx) or CSV

Module 6 requires a settings file defining analysis ranges, grouping, and control samples. Without it, either skip Module 6 or provide the settings and run Module 6 alone:

```bash
# Skip Module 6
./duke --run_modules 1,2,3,4,5,7 ...

# Run Module 6 only after creating settings
./duke --run_modules 6 --path_settings my_settings.xlsx --resume TRUE ...
```

**Demo settings:** `www/settings_example.xlsx`

##### Settings file format

Each row represents one sample.

**Required columns:**

| Column | Description | Example |
|--------|-------------|---------|
| `file_name` | Sample filename (must match exactly) | `sample_001.bam` |
| `analysis_ranges` | Range definitions for peak detection | `[0-35][36-NA]` |
| `floor` | Minimum read frequency threshold per range | `[3][3]` |
| `max_peaks` | Maximum peaks to detect per range | `[2][1]` |

**Optional columns:**

| Column | Description | Example |
|--------|-------------|---------|
| `manual_control_repeat_length` | Known repeat lengths for validation | `[14][100]` |
| `group` | Sample grouping for comparisons | `patient_A`, `control` |
| `group_control_sample` | Flag as group control | `TRUE` / `FALSE` |
| `time` | Timepoint for longitudinal analysis | `0`, `6`, `12` |
| `exclude` | Exclude from analysis | `TRUE` / `FALSE` |

##### Bracketed parameter format

Each bracketed value corresponds to one analysis range. The number of values must be consistent across all bracketed columns in a row.

```
# Two ranges: normal (0-35) and expanded (36+)
analysis_ranges: [0-35][36-NA]
floor:           [3][3]
max_peaks:       [2][1]

# Three ranges
analysis_ranges: [0-13][14-35][36-NA]
floor:           [5][3][3]
max_peaks:       [1][2][1]
```

**Rules:** No spaces between brackets. Use `NA` for unbounded upper limits. Equal numbers of `[` and `]` per cell.

##### Creating your settings file

1. Copy the template: `cp www/settings_example.xlsx my_settings.xlsx`
2. Set `file_name` to match your input files exactly (including extension)
3. Define `analysis_ranges` based on your biology (e.g. `[0-35][36-NA]` for normal/expanded HTT)
4. Set `floor` (minimum reads per peak; typically 3–5) and `max_peaks` (typically 1–2 per range)
5. Add `group`, `group_control_sample`, and `time` columns if doing comparative or longitudinal analysis

---

#### 5. Manual read exclusions (optional)

**Parameter:** `--path_manual_exclusions`  
**Format:** Excel (.xlsx) or CSV. Excel sheet must be named **"Exclusions"** (case-sensitive).

**Required columns:**

| Column | Description | Example |
|--------|-------------|---------|
| `file_name` | Sample filename without extension | `sample_001` |
| `read_name` | Exact read identifier | `m64011_190830_220126/4194640/ccs` |
| `reason` | Documentation | `Low quality alignment` |

Get exact read names from BAM files:
```bash
samtools view sample_001.bam | cut -f1 | head
```

Read names must match exactly — copy directly from `samtools` output.

**Demo exclusions:** `www/read_exclusions_example.xlsx`

---

## Running Duke

### 1. Command-line interface (recommended)

No file editing required. Pass all parameters on the command line.

```bash
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/output \
  --path_ref /path/to/reference.fasta \
  --path_trim_patterns /path/to/adapters.csv \
  --path_settings /path/to/settings.xlsx \
  --threads 12
```

```bash
./duke --help    # Show all 54+ options with defaults
```

### 2. Script-based (duke_run.R)

Edit `scripts/duke_run.R` directly, then run:

```bash
# Local
Rscript scripts/duke_run.R

# HPC — edit paths in the script, then submit
nano scripts/duke_myriad.sh
qsub scripts/duke_myriad.sh
```

All parameters and their defaults are documented with examples in `duke_run.R`.

### 3. Interactive RStudio

Open `scripts/duke_run.R` in RStudio, edit parameters, and source or run line-by-line. Useful for development and debugging individual modules.

---

## Output files

Results are written to `--dir_out`. Each module produces an HTML report and an Excel file.

```
result_duke/
├── 01_import_and_qc.html
├── 02_alignment.html
├── 03_repeat_detection.html
├── 04_allele_calling.html
├── 05_waterfall.html
├── 06_range_analysis.html
├── 07_repeat_visualisation.html
│
├── module_data/                      # Cached RData files (for resume)
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
│   ├── read_counts/                  # Per-read repeat counts (one .tsv.gz per sample)
│   │   └── {sample}.tsv.gz          # read_id | repeat_count_full | repeat_count_match | repeat_count_tracts
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

**HTML reports** are the primary output. Each contains interactive tables and a `# Module diagnostics` section at the end with timing and an output manifest (see [Module diagnostics](#module-diagnostics)).

**Excel files** contain the same data as the HTML tables. Each includes `timing` and `manifest` sheets.

**Per-read count files** (`read_counts/{sample}.tsv.gz`) are written by Module 3 when `export_read_counts = TRUE` (default). One file per sample, one row per read, with all three repeat count methods. Load in Python with `pd.read_csv('sample.tsv.gz', sep='\t')` or in R with `readr::read_tsv('sample.tsv.gz')`.

---

## Module overview

### Module 1: Import and QC

Imports sequencing files, runs quality metrics, and optionally trims adapters.

- Accepts BAM, FASTQ, FASTQ.gz, FASTA
- Duplicate read name detection (keeps longest; safe to disable for PacBio CCS with `--check_duplicate_readnames FALSE`)
- Optional downsampling (`--downsample`)
- Adapter trimming from both ends using `Biostrings::trimLRPatterns`
- **Output:** `qc.xlsx` with read counts and trim statistics; QC plots

### Module 2: Alignment

Aligns reads to left and right reference flanks independently using minimap2.

- Short-read mode (`-x sr -w 1`) works for both ONT and PacBio amplicons
- Strand correction — reverse-complements antisense reads; marks corrected reads as `+c`
- Left→right coordinate order validation removes spurious alignments
- Optional manual read exclusion (`--path_manual_exclusions`)
- Optional coverage visualisation plots
- **Output:** `alignment.xlsx` with alignment quality metrics; coverage plots

### Module 3: Repeat detection

Finds and counts repeat tracts in the sequence between the two flanks.

- Tolerates within-tract gaps and interruptions (configurable)
- Three counting methods calculated for all reads:
  - `repeat_count_full` — tract length ÷ pattern length (recommended; robust for ONT)
  - `repeat_count_match` — exact pattern matches only (stringent; best for PacBio HiFi)
  - `repeat_count_tracts` — tract structure as string e.g. `"3,2,2"` (QC only; not numeric)
- Flank length QC — total combined flank length should be stable per sample; optional outlier filtering (`--rm_flank_length_outliers`)
- **Output:** `repeat_detection.xlsx`; per-read `.tsv.gz` count files; flank QC plots

### Module 4: Allele calling

Clusters reads into alleles using repeat length and/or haplotype sequence.

- ⚠️ Most memory-intensive module — requires 4 GB+ per core
- Default: cluster by repeat length only (`--cluster_by repeat`)
- Optional haplotype sub-clustering using Levenshtein distance on flank sequences
- Consensus sequences generated per cluster using `DECIPHER`
- Variant calling versus reference (VCF export)
- **Output:** `allele_calling.xlsx`; `consensus/` directory with FASTA and VCF files; cluster plots

### Module 5: Waterfall plots

Read-level visualisation showing base composition across left flank, repeat region, and right flank.

- Reads stacked vertically, ordered by repeat length — allele structure and somatic instability immediately visible
- Configurable downsampling (`--waterfall_downsample`; default 1,000 reads)
- Optional per-cluster plots (`--waterfall_per_cluster`)
- Plots saved to `05_waterfall/plots/` — not printed inline in HTML
- **Output:** One PNG per sample

### Module 6: Range analysis

Quantitative analysis of repeat distributions within user-defined length ranges. Requires a settings file.

- Range-based peak detection identifies modal repeat lengths per allele per range
- ~30 distribution statistics per sample per range (mean, median, SD, skewness, CV, tail ratios...)
- Instability metrics: expansion index, contraction index, z-scores relative to controls
- Group and timepoint comparisons
- **Output:** `range_analysis.xlsx` with comprehensive tables; instability and distribution plots

### Module 7: Repeat visualisation

Publication-quality repeat length distribution figures.

- Full-range and per-range histograms (one file per sample) — controlled by `repeat_histogram`
- Cohort summary box+jitter and per-sample horizontal scatter/violin plots — controlled by `repeat_scatter`
- Per-sample ggridges density plots, full range and per analysis range — controlled by `repeat_density`
- Metric overlays (modal, mean, median) on scatter and density plots — controlled by `repeat_distribution_metrics`
- All plots saved to disk, not shown inline in HTML
- **Output:** PNG files organised by plot type in `07_repeat_visualisation/plots/`

---

## Parameters

### Essential parameters

```bash
# Required
--dir_data           Input directory (BAM, FASTQ, etc.)
--dir_out            Output directory
--path_ref           Reference FASTA with NNNNN repeat separator

# Required if trimming enabled (trim = TRUE by default)
--path_trim_patterns Adapter CSV file

# Required for Module 6
--path_settings      Settings Excel file

# Common runtime options
--threads            CPU cores (default: 12)
--resume             Skip completed modules (default: TRUE)
--run_modules        Which modules to run, e.g. 1,2,3,4,5,7 (default: 1,2,3,4,5,6,7)
```

### Complete parameter reference

#### File paths
- `--dir_data` — Input directory **(required)**
- `--dir_out` — Output directory **(required)**
- `--path_ref` — Reference FASTA **(required)**
- `--path_trim_patterns` — Adapter CSV (required if `--trim TRUE`)
- `--path_settings` — Settings Excel (required for Module 6)
- `--path_manual_exclusions` — Read exclusions file (optional)

#### Import options
- `--import_patterns` — File extensions to import
- `--downsample` — Limit reads per sample (NA = all)
- `--select_one_of_pair` — For paired-end: `"R1"`, `"R2"`, or NA
- `--check_duplicate_readnames` — Remove duplicate read names, keeping longest (default: TRUE; safe to disable for PacBio CCS)

#### Adapter trimming
- `--trim` — Enable/disable (default: TRUE)
- `--trim_max_mismatch` — Max errors allowed in adapter match (default: 3)
- `--trim_with_indels` — Allow indels in adapter matching (default: TRUE)

#### Alignment
- `--minimap2_args` — minimap2 settings (default: `"-t 2 -x sr -w 1"`)
- `--visualise_alignment` — Generate coverage plots (default: TRUE)
- `--visualise_alignment_downsample` — Max reads to plot (default: 1,000)

#### Repeat detection
- `--rpt_pattern` — Repeat motif (default: `"CAG"`)
- `--rpt_min_repeats` — Minimum tract length in repeat units (default: 2)
- `--rpt_max_mismatch` — Error tolerance per repeat unit (default: 0)
- `--rpt_start_perfect_repeats` — Perfect repeats required at tract start (default: 2)
- `--rpt_end_perfect_repeats` — Perfect repeats required at tract end (default: 2)
- `--rpt_max_gap` — Max inserted bases within a tract (default: 6)
- `--rpt_max_tract_gap` — Max gap between tracts to merge them (default: 18)
- `--rpt_return_option` — `"longest"` or `"all"` (default: `"longest"`)
- `--repeat_count_method` — Primary counting method (default: `"repeat_count_full"`)
- `--na_repeat_handling` — How to handle reads with no tract: `"convert_to_zero"`, `"filter"`, or `"flag_only"` (default: `"convert_to_zero"`)
- `--export_read_counts` — Export per-read counts to `.tsv.gz` (default: TRUE)

#### Flank length QC (Module 3)
- `--rm_flank_length_outliers` — Filter outlier combined flank lengths (default: FALSE — review flank plots first)
- `--flank_iqr_multiplier` — IQR multiplier for outlier detection (default: 1.5; only active when `rm_flank_length_outliers = TRUE`)

#### Clustering
- `--cluster` — Enable clustering (default: TRUE)
- `--cluster_by` — `"repeat"`, `"haplotype"`, or both (default: `"repeat"`)
- `--haplotype_cluster_max` — Max haplotype clusters (default: 10)
- `--repeat_cluster_max` — Max repeat clusters (default: 20)
- `--repeat_cluster_downsample` — Max reads for cluster optimisation (default: 10,000)
- `--haplotype_cluster_downsample` — Max reads for haplotype optimisation (default: 5,000)
- `--haplotype_region` — Which flanks to use: `"left"`, `"right"`, or `"both"` (default: `"both"`)
- `--haplotype_method` — Distance metric: `"levenshtein"` or `"hamming"` (default: `"levenshtein"`)

#### Consensus and variant calling
- `--cluster_consensus` — Generate consensus sequences (default: TRUE)
- `--consensus_threshold` — Minimum proportion for base call (default: 0.5)
- `--consensus_downsample` — Max reads per consensus (default: 50)
- `--call_variants` — Call variants vs reference (default: TRUE)

#### Waterfall plots (Module 5)
- `--waterfall` — Generate waterfall plots (default: TRUE)
- `--waterfall_downsample` — Max reads per sample (default: 1,000)
- `--waterfall_per_cluster` — Generate per-cluster plots (default: FALSE)
- `--waterfall_y_axis_labels` — Label density: `"auto"` or integer (default: `"auto"`)

#### Range analysis (Module 6)
- `--range_peak_span` — Peak smoothing window (default: 3)
- `--group_control` — Enable control comparisons (default: TRUE)
- `--control_sample_selection` — `"flagged"`, `"earliest"`, or `"all"` (default: `"flagged"`)
- `--control_setpoint_metric` — `"modal_length"`, `"mean_length"`, or `"median_length"` (default: `"median_length"`)
- `--control_aggregation_method` — `"mean"`, `"median"`, or `"trimmed_mean"` (default: `"median"`)

#### Repeat visualisation (Module 7)
- `--repeat_histogram` — Generate frequency histograms (full-range and per analysis range) (default: TRUE)
- `--repeat_histogram_binwidth` — Bin width in repeat units (default: 1)
- `--repeat_scatter` — Generate cohort summary box+jitter and per-sample scatter/violin plots (default: TRUE)
- `--repeat_density` — Generate per-sample ggridges density plots, full range and per analysis range (default: TRUE)
- `--repeat_distribution_metrics` — Metric overlays on scatter and density plots: `"modal_length"`, `"mean_length"`, `"median_length"` (default: all three)

#### Plot output
- `--plot_dpi` — Resolution for diagnostic plots (default: 150; use 300 for publication)
- `--plot_per_sample` — Generate per-sample plot files in Modules 2, 3, and 4 (Module 2: coverage/strand/alignment/segment; Module 3: scatter/histogram/violin/density; Module 4: violin/scatter). Consider FALSE for large datasets (>50 samples) (default: TRUE)

#### Runtime
- `--threads` — CPU cores (default: 12)
- `--resume` — Skip completed modules (default: TRUE)
- `--remove_intermediate` — Free RAM between modules (default: TRUE)
- `--remove_temp` — Delete temp files on completion (default: FALSE)
- `--run_modules` — Which modules to run (default: `1,2,3,4,5,6,7`)

---

## Module diagnostics

Every module ends with a `# Module diagnostics` section in the HTML report, just before Session info.

### Timing table

All seven modules record wall-clock time for each major computation step:

| Column | Description |
|---|---|
| Section | Name of the computation step |
| Seconds | Raw elapsed seconds |
| Elapsed (mm:ss) | Step duration formatted as mm:ss or hh:mm:ss |
| Cumulative (mm:ss) | Running total from module start |

The **Total** row (bold) gives the complete module wall time. Useful for identifying bottlenecks and confirming resume is working correctly.

### Output manifest

Modules 1, 2, 3, 4, and 6 also display an output manifest listing every R object saved to the `.RData` file:

| Column | Description |
|---|---|
| Object | R object name (accessible via `moduleN_output$...`) |
| Description | Contents and intended downstream use |
| Class | R class (data.frame, tbl_df, list, etc.) |
| In-memory size | Human-readable size (B / KB / MB / GB) |

The table caption shows the `.RData` file size on disk and full path.

### Excel export

Both tables are appended to each module's Excel file as `timing` and `manifest` sheets. Modules 5 and 7 (plot-only modules) include only a `timing` sheet.

---

## HPC deployment

### Cluster selection

| Cluster | Best for | Cores | Scheduler | Script |
|---------|----------|-------|-----------|--------|
| **Myriad** | ≤500 samples | 1–36 | SGE (`-pe smp`) | `duke_myriad.sh` |
| **Old Kathleen** | 1000+ samples | 80–160 | SGE (`-pe mpi`) | `duke_kathleen.sh` |
| **New Kathleen** | 1000+ samples | 80–160 | Slurm | `duke_kathleen_slurm.sh` |

### Job submission

Edit the data/output paths in the job script, then submit. The `--threads` parameter is set automatically from `$NSLOTS` (SGE) or `$SLURM_NTASKS` (Slurm) — you only need to change the core count in the scheduler header line.

#### Myriad (SGE)
```bash
nano scripts/duke_myriad.sh    # Edit data/output paths
qsub scripts/duke_myriad.sh
```

#### Old Kathleen (SGE)
```bash
nano scripts/duke_kathleen.sh
qsub scripts/duke_kathleen.sh
```

#### New Kathleen (Slurm)
```bash
nano scripts/duke_kathleen_slurm.sh
sbatch scripts/duke_kathleen_slurm.sh
```

### Monitoring jobs

**SGE (Myriad / Old Kathleen):**
```bash
qstat -u $USER                  # Check job status
watch -n 5 'qstat -u $USER'    # Watch queue
qstat -j <JOB_ID>              # Detailed info
qdel <JOB_ID>                  # Cancel
qacct -j <JOB_ID>              # After completion
tail -f logs/duke_*.out         # Watch log in real-time
grep -i error logs/duke_*.err  # Check for errors
```

**Slurm (New Kathleen):**
```bash
squeue -u $USER
watch -n 5 'squeue -u $USER'
scontrol show job <JOB_ID>
scancel <JOB_ID>
sacct -j <JOB_ID>
tail -f logs/duke_*.out
```

### SGE vs Slurm reference

| Action | SGE | Slurm |
|--------|-----|-------|
| Submit | `qsub script.sh` | `sbatch script.sh` |
| Check jobs | `qstat -u $USER` | `squeue -u $USER` |
| Job details | `qstat -j <JOB_ID>` | `scontrol show job <JOB_ID>` |
| Cancel | `qdel <JOB_ID>` | `scancel <JOB_ID>` |
| After completion | `qacct -j <JOB_ID>` | `sacct -j <JOB_ID>` |
| Job ID variable | `$JOB_ID` | `$SLURM_JOB_ID` |

### Resource recommendations

#### Myriad

**Recommended configuration:**
```bash
#$ -pe smp 12
#$ -l mem=64G        # 64 GB per core = 768 GB total
#$ -l tmpfs=500G
#$ -l h_rt=48:00:00
```

With Duke:
```bash
./duke --threads $NSLOTS ...
```

**Run modes**

| Mode | Modules | Flags | When to use |
|------|---------|-------|-------------|
| Full | M1–M7 | `--trim TRUE --waterfall TRUE` (defaults) | First-pass analysis; when per-read waterfall inspection is needed for QC |
| Slim | M1–M4, M6–M7 | `--trim FALSE --waterfall FALSE` | When adapters are absent or already stripped, and waterfall plots are not required; saves 5–30% runtime |

In slim mode, Module 1 still runs in full but skips the adapter-trimming step. Module 5 (waterfall plots) is omitted entirely. All other modules and outputs are identical to a full run.

**Benchmark results** (Myriad, 12 cores / 64 GB, v2.2.0):

All timings are from complete fresh runs (all modules executed from scratch, no cached results re-used). v1.0 = January 2026; v2.2.0 = April 2026. "—" = configuration not run for that dataset. csf timings are approximate, derived from a concurrent multi-job log.

| Dataset | Samples | Total reads | Reads / sample | v1.0 full | v2.2.0 full | v2.2.0 slim | Change (full) |
|---------|--------:|------------:|---------------:|----------:|------------:|------------:|:-------------:|
| jasmine | 102 | ~2.9 M | ~28 k | 8.1 h | — | — | — |
| yas | 77 ¹ | 3.6 M | ~46 k | 9.6 h | 5.3 h | **5.0 h** | −45% |
| jd | 131 | 4.4 M | ~34 k | — | 6.3 h | — | new dataset |
| pg | 285 ¹ | 3.6 M | ~13 k | 9.4 h | 8.3 h | **5.9 h** | −12% |
| csf | ~100 | ~4 M | ~40 k | — | ~17.5 h ² | ~3.9 h ² | new dataset |
| fs | 310 | 10.6 M | ~34 k | 28–34 h | 27.6 h | — | ~same |

¹ Sample count differs slightly from v1.0 due to QC exclusions applied between runs.  
² Timing approximate; two jobs ran concurrently and shared a log file.

**Key findings:**

- **Total reads is the strongest runtime predictor.** Datasets under 5 M reads consistently finish in 5–9 h; the fs dataset at 10.6 M reads takes ~28 h regardless of pipeline version. Sample count is a secondary factor.
- **v2.2.0 is materially faster for most datasets.** The yas dataset improved by 45%; pg by 12%. The fs dataset (highest read depth) showed no meaningful change, consistent with it being I/O- and read-throughput-bound rather than compute-bound.
- **Slim mode offers significant savings for larger sample counts.** The pg dataset (285 samples) went from 8.3 h to 5.9 h (−29%) with `--trim FALSE --waterfall FALSE`. Savings are smaller for low sample counts (yas: 5.3 → 5.0 h, −5%).
- **12 cores / 64 GB remains optimal.** More RAM per core does not improve performance for most datasets; more cores with less memory can cause failures on large datasets.

**Configurations to avoid:**
```bash
# Slower for most datasets
#$ -pe smp 6
#$ -l mem=128G

# Can fail on large datasets
#$ -pe smp 24
#$ -l mem=32G
```

**Runtime estimation:**

Use total reads (R, in millions) and sample count (N) to estimate v2.2.0 full-run time on Myriad at 12c / 64G:

```
Runtime_full (h)  ≈  8.94 − 2.68R + 0.377R² + 0.014N
Runtime_slim (h)  ≈  Runtime_full × max(0.5,  1 − 0.001N)
```

Uncertainty is approximately ±25% within the calibrated range (N ≤ 310, R ≤ 10.6 M). For larger datasets, treat estimates as indicative only.

**Quick reference:**

| Total reads | ~100 samples | ~200 samples | ~300 samples |
|-------------|:------------:|:------------:|:------------:|
| 3–4 M | 5–6 h | 6–8 h | 7–9 h |
| 5–7 M | 7–10 h | 9–12 h | 11–14 h |
| 8–11 M | 15–22 h | 17–25 h | 20–28 h |
| >11 M | 28–40 h | 30–42 h | 33–45 h |

For `h_rt`, add a 20% safety margin: if you estimate 10 h, request 12 h.

#### Kathleen

**Old Kathleen (SGE):**
```bash
#$ -pe mpi 80
#$ -l mem=4G
#$ -l h_rt=48:00:00
# --threads $NSLOTS  (set automatically from scheduler)
```

**New Kathleen (Slurm):**
```bash
#SBATCH --nodes=2
#SBATCH --ntasks=80
#SBATCH --mem-per-cpu=4G
#SBATCH --time=48:00:00
# --threads $NSLOTS  (set automatically from scheduler)
```

**Memory workaround** (if 4 GB/core unavailable):
```bash
#$ -pe mpi 160
#$ -l mem=2G         # 320 GB total
# --threads 80  (not 160 — override required for memory workaround)
```

#### Resource summary

Myriad runtimes are from v2.2.0 benchmarks (full mode, 12c / 64G). Slim mode (`--trim FALSE --waterfall FALSE`) reduces runtimes by ~5–30% depending on sample count. Kathleen runtimes are indicative estimates only — no benchmark data has been collected on Kathleen.

| Dataset size | Samples | Cluster | Cores | Memory/core | tmpfs | Runtime (full) | Runtime (slim) |
|--------------|--------:|---------|------:|------------:|------:|:--------------:|:--------------:|
| Small | <100 | Myriad | 12 | 64 GB | 500 GB | 5–9 h | 4–8 h |
| Medium | 100–300 | Myriad | 12 | 64 GB | 500 GB | 6–10 h | 5–8 h |
| Large | 300–500 | Myriad | 12 | 64 GB | 500 GB | 15–35 h | 10–25 h |
| Very large | 500–1000 | Kathleen | 80 | 4 GB | — | est. 16–24 h | — |
| Massive | 1000+ | Kathleen | 80–160 | 4 GB | — | est. 24–48 h | — |

> Runtime depends strongly on **total reads**, not sample count alone. Two datasets of 300 samples can differ by 3× in runtime if one has 3× the read depth. Use the estimation formula in the Myriad section above.

---

## Run planning

Empirical data from multiple PacBio Revio flow cells for HTT amplicon sequencing:

| Samples per flow cell | Total reads | Reads per sample |
|-----------------------|-------------|------------------|
| 250 | ~7.2 million | ~28,700 |
| 287 | ~8.7 million | ~30,200 |
| 288 | ~7.5 million | ~26,200 |
| 381 | ~8.5 million | ~22,400 |

Each Revio flow cell yields approximately 7–9 million aligned reads. For most analyses, **>5,000 reads per sample** provides robust allele calling and instability metrics.

| Project size | Flow cells | Expected reads/sample |
|--------------|------------|----------------------|
| ≤300 samples | 1 | ~25,000–30,000 |
| 300–400 samples | 1 | ~20,000–25,000 |
| >400 samples | Consider 2 | ~25,000+ |

---

## Common workflows

### Workflow 1: Discovery run (first time with new data)

Run without a settings file to explore the data. Use the repeat histograms in Module 3 and Module 7 to identify the allele distribution and define appropriate analysis ranges for Module 6.

```bash
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/results \
  --path_ref www/HTTset20.fasta \
  --path_trim_patterns www/adapters.csv \
  --run_modules 1,2,3,4,5,7 \
  --threads 12
```

### Workflow 2: Full analysis with instability metrics

```bash
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/results \
  --path_ref www/HTTset20.fasta \
  --path_trim_patterns www/adapters.csv \
  --path_settings /path/to/settings.xlsx \
  --threads 12
```

### Workflow 3: Quality control with read exclusions

```bash
# 1. Initial run
./duke --dir_data /path/to/data --dir_out /path/to/results --path_ref www/HTTset20.fasta

# 2. Review Module 2 alignment plots and Module 5 waterfall plots
#    Identify problematic reads visually

# 3. Extract read names
samtools view sample_001.bam | cut -f1 | head

# 4. Create exclusions file (see Input files section)

# 5. Re-run from Module 2 with exclusions
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/results \
  --path_ref www/HTTset20.fasta \
  --path_manual_exclusions /path/to/exclusions.xlsx \
  --run_modules 2,3,4,5,6,7
```

### Workflow 4: Longitudinal analysis

Settings file with `group`, `time`, and `group_control_sample` columns:

```csv
file_name,group,time,group_control_sample,analysis_ranges,floor,max_peaks
patient_A_t0.bam,patient_A,0,TRUE,[0-35][36-NA],[3][3],[2][1]
patient_A_t6.bam,patient_A,6,FALSE,[0-35][36-NA],[3][3],[2][1]
patient_A_t12.bam,patient_A,12,FALSE,[0-35][36-NA],[3][3],[2][1]
```

```bash
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/results \
  --path_ref www/HTTset20.fasta \
  --path_settings /path/to/settings.xlsx \
  --control_sample_selection flagged
```

Module 6 outputs include timepoint progression, baseline comparisons, and instability metrics over time.

### Workflow 5: Large dataset HPC run

```bash
# Edit paths in job script
nano scripts/duke_myriad.sh

# Submit
qsub scripts/duke_myriad.sh

# Monitor
watch -n 5 'qstat -u $USER'
tail -f logs/duke_*.out

# Resume if job times out (resume = TRUE by default)
qsub scripts/duke_myriad.sh
```

### Workflow 6: Re-run specific modules

```bash
# Update settings and re-run Module 6 only
./duke --run_modules 6 --path_settings updated_settings.xlsx --resume TRUE ...

# Re-run from Module 2 with updated exclusions
./duke --run_modules 2,3,4,5,6,7 --path_manual_exclusions updated_exclusions.xlsx ...

# Force re-run of a specific module by deleting its cached output
rm result_duke/module_data/03_repeat_detection_results.RData
./duke --run_modules 3,4,5,6,7 ...
```

### Workflow 7: Parameter optimisation

```bash
# Test with downsampling to iterate quickly
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/test_params \
  --path_ref www/HTTset20.fasta \
  --downsample 500 \
  --rpt_max_mismatch 0 \
  --threads 4

# Adjust and re-run
./duke \
  --dir_data /path/to/data \
  --dir_out /path/to/test_params_v2 \
  --path_ref www/HTTset20.fasta \
  --downsample 500 \
  --rpt_max_mismatch 1 \
  --rpt_max_gap 9 \
  --threads 4
```

---

## Troubleshooting

**"Error: --path_ref is required"**
```bash
./duke --path_ref /path/to/reference.fasta ...
# Or use the demo: ./duke --path_ref www/HTTset20.fasta ...
```

**"duke: command not found"**
```bash
chmod +x duke
cd ~/Scratch/bin/duke
./duke --help
```

**"Error: there is no package called 'openxlsx'"**
```bash
export R_LIBS_USER=~/R/library
# Or install — see Installation section
```

**"scripts/duke_cli.R not found"**
```bash
head duke    # Should show: Rscript "$SCRIPT_DIR/scripts/duke_cli.R" "$@"
# Run from the duke root directory
```

**"Settings file not found" or Module 6 fails without settings**
```bash
# Skip Module 6
./duke --run_modules 1,2,3,4,5,7 ...

# Or create a settings file and provide it
./duke --path_settings my_settings.xlsx ...
```

**"Cannot parse analysis_ranges" or "Bracketed parameter error"**

Format must be `[value1][value2]` with no spaces between brackets. Use `NA` for unbounded upper limits. All bracketed columns must have the same number of values per row.

```
✅  [0-35][36-NA]
❌  [0-35] [36-NA]    (space between brackets)
❌  [0-35,36-NA]      (comma instead of ][)
❌  [0-35][36-100]    (should be [36-NA] for unbounded)
```

**"Mismatch in bracketed parameter lengths"**

```
✅  analysis_ranges: [0-35][36-NA]   floor: [3][3]   max_peaks: [2][1]
❌  analysis_ranges: [0-35][36-NA]   floor: [3]       max_peaks: [2][1]
```

**"file_name not found" in settings file**
```bash
ls -1 /path/to/data    # Check exact filenames including extension
# file_name must match exactly, e.g. sample_001.bam (not sample_001)
```

**Manual exclusions not working**
```bash
samtools view input.bam | cut -f1 | head -5    # Get exact read names
# file_name: use base name without extension (sample_001, not sample_001.bam)
# Excel sheet must be named "Exclusions" (case-sensitive)
# Check Module 2 output for: "Alignments excluded: X"
```

**Job stuck at Module 4 (very slow)**

Insufficient memory. Use at least 4 GB/core, 8 GB recommended:
```bash
#$ -l mem=8G
```

**"Cannot allocate memory" (Kathleen)**
```bash
# Request double cores, use half for processing
#$ -pe mpi 160
#$ -l mem=2G
# --threads $NSLOTS  (set automatically from scheduler)
```

**"Error: cannot open the connection" when running Rmd**
```bash
# Verify scripts have the knit_root_dir fix (v2.1.0+)
grep "knit_root_dir" scripts/duke_run.R    # Should return a line
```

---

## Version history

### v2.2.0 (April 2026) — Current
- ✨ **NEW:** Module diagnostics in every HTML report — timing table (all 7 modules) and output manifest (modules 1–4, 6); both exported to each module's Excel file
- ✨ **NEW:** `export_read_counts` — per-read repeat counts exported to `{sample}.tsv.gz` (default TRUE)
- ✨ **NEW:** `check_duplicate_readnames` parameter (Module 1; default TRUE)
- ✨ **NEW:** Flank length QC filtering — `rm_flank_length_outliers`, `flank_iqr_multiplier` (Module 3)
- ✨ **NEW:** Plot output control — `plot_dpi`, `plot_per_sample` (Modules 2, 3, 4)
- ✨ **NEW:** `format_size()`, `format_elapsed()`, `build_timing_table()` helpers in `00_utils.R`
- 🐛 **FIXED:** Modules 5 and 7 no longer have spurious resume checks
- 🐛 **FIXED:** Manifest in-memory sizes now use adaptive units (B/KB/MB/GB)
- 🐛 **FIXED:** Summary section crash when running partial module sets
- 🐛 **FIXED:** Reference N-masking now case-insensitive (`NNNNN` or `nnnnn`)
- 🐛 **FIXED:** Clustering early-exit return type inconsistency
- 🐛 **FIXED:** GMM connection leak in Module 4 clustering
- 🔧 **REMOVED:** `log_dir`, `verbose`, `waterfall_per_sample`, `repeat_histogram_per_sample`, `repeat_scatter_per_sample`
- 🔧 **UPDATED:** `consensus_downsample` default 200 → 50; `waterfall_per_cluster` default TRUE → FALSE
- 📄 **STREAMLINED:** HPC job scripts; README restructured for new users

### v2.1.0 (January 2026)
- 🐛 **FIXED:** `knit_root_dir` handling for Rmd files in subdirectories
- 📚 **ENHANCED:** Comprehensive reorganisation documentation
- 📊 **ADDED:** Job monitoring commands in README
- 📖 **IMPROVED:** Complete input file format documentation
- 🎯 **ADDED:** Demo reference (HTTset20.fasta) and adapters (adapters.csv)
- 📋 **ADDED:** Common workflows section

### v2.0.1 (January 2025)
- ✨ **NEW:** Command-line interface
- ✅ **NEW:** `trim` parameter (optional adapter trimming)
- 🔧 **RENAMED:** `visualise_alignment_downsample` (was: `_n_reads`)
- 🐛 **FIXED:** CLI resume detection; resume default now TRUE
- 📦 **ORGANISED:** Library files with module prefixes (00–07)

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
UCL Queen Square Institute of Neurology, London, UK

- Email: michael.flower@ucl.ac.uk
- UCL Profile: https://profiles.ucl.ac.uk/45681-michael-flower
- ORCID: https://orcid.org/0000-0001-5568-6239
- GitHub: https://github.com/mike-flower/duke
- Issues: https://github.com/mike-flower/duke/issues
