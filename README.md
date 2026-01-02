# Duke Pipeline 2.0

A modular pipeline for amplicon sequencing analysis with comprehensive repeat length characterisation and instability metrics.

## Table of Contents

**Quick Start:**
- [Installation](#installation)
- [File Structure](#file-structure)
- [Run Analysis](#run-analysis)
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
| `floor` | Peak detection threshold(s), one per range | `[14]` or `[14][100]` |
| `max_peaks` | Maximum peaks per range | `[3]` or `[3][2]` |

### Optional Columns

| Column | Description | Example |
|--------|-------------|---------|
| `group` | Sample grouping for comparative analysis | `patient1`, `control` |
| `group_control_sample` | Any non-blank value = control sample | `1`, `TRUE`, or blank |
| `time` | Timepoint for temporal analysis | `0`, `6`, `12` |
| `manual_control_repeat_length` | Manual control setpoint(s), one per range | `[24][120]` or blank |
| `exclude` | Any non-blank value = exclude from analysis | `1`, `TRUE`, or blank |

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
| `floor` | `[val1][val2]` | `[14][100]` |
| `max_peaks` | `[val1][val2]` | `[3][2]` |
| `manual_control_repeat_length` | `[val1][val2]` or blank | `[24][120]` |

**Example row:**
```
file_name: sample_001.fastq.gz
analysis_ranges: [0-35][36-200]
floor: [14][100]
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
│   └── qc/
│       └── qc.xlsx
├── 02_alignment.html
├── 02_alignment/
│   ├── plots/
│   └── alignment/
│       └── alignment_qc.xlsx
├── 03_repeat_detection.html
├── 03_repeat_detection/
│   ├── plots/
│   └── repeat_detection/
│       └── repeat_summaries.xlsx
├── 04_allele_calling.html
├── 04_allele_calling/
│   ├── plots/
│   ├── consensus/
│   │   ├── consensus_sequences.fasta
│   │   └── consensus_summary.xlsx
│   └── variants/
│       ├── *.vcf
│       └── variant_summary.xlsx
├── 05_waterfall.html
├── 05_waterfall/
│   └── plots/
├── 06_range_analysis.html
├── 06_range_analysis/
│   ├── plots/
│   │   ├── distribution_by_sample/
│   │   ├── distribution_by_range/
│   │   ├── instability_by_sample/
│   │   └── peaks_by_sample/
│   └── range_analysis/
│       └── range_analysis_results.xlsx
├── 07_repeat_visualisation.html
├── 07_repeat_visualisation/
│   └── plots/
│       ├── scatter_combined.png
│       ├── histograms_by_sample/
│       └── scatter_by_sample/
├── module_data/
│   └── *.RData (for resume functionality)
└── temp/
    └── (intermediate files, removed if cleanup_temp = TRUE)
```

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

### Version 2.0 (Current)
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
