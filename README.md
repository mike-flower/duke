# Duke

A modular pipeline for amplicon sequencing of tandem repeats: per-read repeat counts, optional UMI collapse, allele structure and somatic instability metrics, with an HTML report for every step.

**v3.3.13** – what changed in each version is in [CHANGELOG.md](CHANGELOG.md).

**New to Duke?** Read [What Duke does](#what-duke-does), then follow [Installation](#installation) and [Quick start](#quick-start). The rest of this page is reference: look things up as you need them.

---

## Contents

- [What Duke does](#what-duke-does)
- [Terms used in Duke](#terms-used-in-duke)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Input files](#input-files)
- [Output files](#output-files)
- [Several repeat tracts in one run](#several-repeat-tracts-in-one-run)
- [UMIs](#umis)
- [How repeats are counted](#how-repeats-are-counted)
- [Instability metrics](#instability-metrics)
- [Parameters](#parameters)
- [Resume](#resume)
- [Running on Myriad](#running-on-myriad)
- [Common workflows](#common-workflows)
- [Troubleshooting](#troubleshooting)
- [Citation and contact](#citation-and-contact)

---

## What Duke does

| Module | Report | What it does | Runs when |
|---|---|---|---|
| 1. Import and QC | `01_import_qc.html` | Import BAM / FASTQ / FASTA; read counts and lengths; PacBio passes and quality; optional np filter, de-duplication, downsampling and adapter trimming | always |
| 2. Alignment | `02_alignment.html` | minimap2 to the left and right flanks of each reference; read accounting from input to analysed reads; flank boundary, MAPQ and strand QC | always |
| 3. Repeat detection | `03_repeat_detection.html` | Finds every tract in every read; counts by three methods; modal length by method; sequence between tracts; composite measures; per-read export | always |
| 4. UMI | `04_umi.html` | Extracts UMIs, clusters them into molecules and collapses each molecule's reads to one length | `umi = TRUE` |
| 5. Allele structure | `05_allele_structure.html` | Clustering into alleles; each allele's consensus sequence and repeat notation; phased flank variants with read support; a phased VCF, allele FASTA and allele-tagged BAMs; waterfall plots | `cluster = TRUE` or `waterfall = TRUE` |
| 6. Instability | `06_instability.html` | Distributions and instability metrics per sample, range and measure, for reads and molecules; group controls; timecourses; tract-against-tract plots | a settings file is given |
| 7. Repeat lengths | `07_repeat_lengths.html` | Every repeat-length figure in one place: cohort heatmaps, density and histogram figures by group (all lengths and each range), optional per-sample figures, and one card per sample | always |

At the end of every run Duke writes `00_index.html`, which links every report, table and figure, and `run_summary.xlsx`, one row per sample from every module. Reads are grey and molecules teal in every figure. Every report has a timing table and an output manifest, and every module writes its tables to one Excel file.

### Terms used in Duke

| Term | Meaning |
|---|---|
| sample, file stem | One input file. Its name without the extension (the file stem) is the sample's identifier in every table and in `00_index.html`. Figures use a short label: the file stem without the dot-separated parts that every sample shares |
| tract | One run of a repeat motif in the amplicon, e.g. the CAG tract. One tract by default (`rpt_pattern`); several with a [tracts table](#several-repeat-tracts-in-one-run) |
| composite | A length made from several tracts per read: `sum` (their counts added) or `span` (first tract start to last tract end) |
| measure | Any tract or composite. Every result is given per measure |
| primary measure | The first tract. It is the one used to select reads into ranges, for clustering and for waterfalls |
| UMI, molecule | A UMI is a random barcode added to each template molecule before PCR; reads with the same UMI are copies of one molecule. Module 4 collapses them to one length per molecule, which removes PCR bias and most PCR errors |
| unit | Whether a result counts reads or molecules. A *combined* figure shows both (reads grey, molecules teal) |
| range, range01, range02, … | A band of repeat lengths set in the settings file, e.g. `[0-35][36-NA]` (normal and expanded). Ranges are named `range01`, `range02`, … in the order of their brackets, and both bounds are included |
| setpoint | The length that instability is measured from: the sample's own median (`sample_relative`) or that of its group's control samples (`control_relative`) |
| allele 1, allele 2 | The two largest clusters of a sample in Module 5, numbered by length (1 = shorter); further clusters are 3 upwards |
| np | The number of passes behind a PacBio HiFi read (more passes, more accurate). Reported up to a cap (60 on Revio): a read with more passes is reported as 60. A Kinnex (skera) segment carries the pass count of the whole concatenated read it was cut from |

---

## Installation

**1. Get Duke.**

```bash
git clone https://github.com/mike-flower/duke.git
cd duke
```

**2. Install the R packages** (in R, once):

```r
install.packages(c("data.table", "dplyr", "tidyr", "tibble", "ggplot2", "patchwork", "scales", "RColorBrewer",
                   "DT", "rmarkdown", "knitr", "openxlsx", "readxl", "mclust", "cluster", "pracma",
                   "digest", "pbapply", "pbmcapply", "ragg"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Biostrings", "ShortRead", "Rsamtools", "GenomicAlignments", "DECIPHER", "pwalign"))
```

`ragg` draws figures much faster (without it Duke uses the cairo device). `pwalign` exists from Bioconductor 3.19; on older versions skip it and Duke uses the aligner in `Biostrings`. `rstudioapi` is optional: with it, `scripts/duke_run_local.R` finds the Duke folder when sourced in RStudio.

**3. Install the command-line tools.** minimap2, samtools and pandoc must be on the `PATH`; bgzip and tabix (from htslib) are optional and add a compressed, indexed copy of the VCF. On a Mac, `brew install minimap2 samtools pandoc htslib`; with conda, `conda install -c bioconda minimap2 samtools htslib`. RStudio comes with pandoc, but the command line needs its own copy.

**4. Check.** `./duke --help` lists every option. If it fails, the message names the missing package or tool.

**On Myriad**, install on a login node (compute nodes have no internet access), into your own R library:

```bash
cd ~/Scratch/bin
git clone https://github.com/mike-flower/duke.git
git clone https://github.com/lh3/minimap2.git
cd minimap2
make
cd ..
module load r/recommended
mkdir -p ~/R/library
export R_LIBS_USER=~/R/library
R
```

then run the `install.packages` and `BiocManager::install` lines above in R and quit. The job scripts load R and samtools, set `R_LIBS_USER` and put `~/Scratch/bin/minimap2` on the `PATH`. To update Duke later: `cd ~/Scratch/bin/duke` and `git pull`.

**Layout:**
```
duke/
├── duke                      # command-line wrapper
├── lib/                      # functions (params.R holds every parameter, its default and help)
├── modules/                  # one Rmd per module
├── scripts/
│   ├── duke_run.R            # all parameters with defaults, cluster paths
│   ├── duke_run_local.R      # the same, laptop paths
│   ├── duke_cli.R            # used by ./duke
│   ├── duke_myriad.sh, duke_kathleen.sh, duke_kathleen_slurm.sh   # job scripts
│   ├── update_readme_params.R  # rewrites the parameter tables in this README from lib/params.R
│   └── update_run_scripts.R    # rewrites duke_run.R and duke_run_local.R from lib/params.R
├── logs/                     # one folder per run: run log, resolved parameters and tracts
└── www/                      # example reference, tracts table and settings files
```

---

## Quick start

You need a folder of reads (one BAM, FASTQ or FASTA file per sample) and a reference with the repeat masked by `N`s (`www/HTTset20_cagcaa.fasta` for HTT CAG). Give full paths: Duke runs from its own folder, so a relative path is read from there (`www/...` works for that reason).

**1. A first look, without a settings file.** Write each run to a new folder, named with the Duke version, so earlier results are kept.

```bash
./duke --dir_data /path/to/run1/data --dir_out /path/to/run1/result_duke_v3.3.13 --path_ref www/HTTset20_cagcaa.fasta --threads 8
```

Add `--umi TRUE` for a UMI design and `--cluster TRUE` for alleles, consensus sequences and variants. Without a settings file every module runs except Module 6. Open `00_index.html` in the output folder: Modules 1 and 2 show whether the reads and alignment are good, and the Module 3 and Module 7 length plots show where each sample's alleles are, so where to put the ranges.

**2. Write a settings file.** Copy `www/settings_example.xlsx`: one row per sample with its file name, ranges, and optionally group, time and control status (see [Settings file](#settings-file-path_settings-needed-for-module-6)).

**3. Run again with the settings file**, into the same folder:

```bash
./duke --dir_data /path/to/run1/data --dir_out /path/to/run1/result_duke_v3.3.13 --path_ref www/HTTset20_cagcaa.fasta --path_settings /path/to/run1/settings_duke.xlsx --threads 8
```

[Resume](#resume) keeps Modules 1–5, as nothing that affects them changed, and runs Modules 6 and 7. Edit the settings and run again as often as needed.

**Other ways to run.** Add `--dry_run` to any command to check the input paths and print the parameters it would use, without running anything. `./duke --help` lists every option by module with its default. From R or RStudio, edit the paths in `scripts/duke_run_local.R` and source it; every parameter is listed there with its default and a one-line explanation. On Myriad, see [Running on Myriad](#running-on-myriad).

---

## Input files

### Sequencing data (`dir_data`)

BAM, FASTQ (`.fastq`, `.fq`, gzipped or not) or FASTA, searched recursively. The sample name (`file_stem`) is the file name without its extension.

### Reference (`path_ref`)

A FASTA with one block of `N` (five or more) where the repeat is. Duke aligns each read to the sequence either side of the block, so the block can be any length. For HTT, `www/HTTset20_cagcaa.fasta` is the reference for the CAG tract: its N block covers the CAG tract and the CAA CAG after it, so the end of the tract is found by tract detection rather than fixed by the flank alignment. `www/HTTset20_cag_ccg.fasta` covers the CAG and the CCG tract together, for CAG and CCG profiling. `www/HTTset20.fasta` (N block on the CAG tract only) is kept for comparison with earlier runs; on the 2026 jasmine run it gave the same modal length as `HTTset20_cagcaa` in every sample, by reads and by molecules.

### Settings file (`path_settings`, needed for Module 6)

One row per sample, `.xlsx` (first sheet) or `.csv`:

| Column | Required | Example | Meaning |
|---|---|---|---|
| `file_name` | yes | `sample_01.bam` | The input file (with or without extension) |
| `analysis_ranges` | yes | `[0-35][36-NA]` | Ranges on the primary measure (the first tract in the tracts table, or the single tract from `rpt_pattern`), in units of its motif, one bracket each, named `range01`, `range02`, … in order. Both bounds are included (`[36-NA]` is 36 and above); `NA` = unbounded |
| `floor` | no | `[3][3]` | Per range: leave out lengths with fewer reads (or molecules) than this; a value below 1 is a share of the commonest length's count (`0.1` = 10%, as in Lee et al. 2010). Blank, 0 or 1 keeps every length |
| `max_peaks` | no | `[2][1]` | Modal peaks reported per range |
| `group` | no | `QS5` | Group for controls, timecourses and plots |
| `time` | no | `0`, `49` | Time point |
| `group_control_sample` | no | `TRUE` | The sample is a control for its group |
| `exclude` | no | `TRUE` | Leave the sample out of Modules 6 and 7 |
| `manual_control_repeat_length` | no | `[14][120]` | Control setpoint per range; where given, it replaces the group control for that sample and range. Leave blank to use the group control |
| `analysis_ranges_<measure>` (and `floor_<measure>`, `max_peaks_<measure>`, `manual_control_repeat_length_<measure>`) | no | `[NA-8][9-NA]` | Ranges for another tract or composite on its own scale (see below) |

Rows whose file is not in `dir_data` are ignored, so one platemap-wide settings file can serve every run; if a file appears in more than one row, the first row is used. Group and range names may contain any characters; those that are unsafe in file names are replaced with `-` in output paths only. Examples: `www/settings_example.xlsx`, `www/settings_example_two_tracts.xlsx` (copy one and replace the rows with your samples).

### Other optional files

- `path_tracts` – the tracts table (next section).
- `path_trim_patterns` – adapter CSV (`adapter_name`, `adapter_sequence`) for `trim = TRUE`; `www/adapters.csv`.
- `path_manual_exclusions` – reads to drop (`file_name`, `read_name`, `reason`; Excel sheet `Exclusions`); `www/read_exclusions_example.xlsx`.

### Path parameters

<!-- params:Paths:start -->

| Parameter | Default | Format | What it does |
|---|---|---|---|
| `dir_data` | – | text | Input directory of BAM / FASTQ / FASTA files (required) |
| `dir_out` | – | text | Output directory (required) |
| `path_ref` | – | text | Reference FASTA with NNNNN marking the repeat region (required unless every tract in path_tracts names its own reference) |
| `path_settings` | – | text | Settings .xlsx/.csv (ranges, groups, controls, time). Needed for Module 6 |
| `path_tracts` | – | text | Optional tracts table (.csv or .xlsx sheet `tracts`; optional sheet `composites`). NA = one tract from rpt_pattern and path_ref |
| `path_trim_patterns` | – | text | Adapter CSV (adapter_name, adapter_sequence), used when trim = TRUE |
| `path_manual_exclusions` | – | text | Optional read exclusions file (file_name, read_name, reason) |

<!-- params:Paths:end -->

---

## Output files

Open **`00_index.html`** first: it links every report, table and figure, by sample, range and group. **`run_summary.xlsx`** has one row per sample with the results of every module.

```
result_duke/
├── 01_import_qc.html             01_import_qc/qc.xlsx, plots/
├── 02_alignment.html             02_alignment/alignment.xlsx (read_accounting), plots/ (alignment/), bam/ (keep_bam)
├── 03_repeat_detection.html      03_repeat_detection/repeat_detection.xlsx, read_repeat_counts/{sample}.tsv.gz, plots/
├── 04_umi.html                   04_umi/umi.xlsx, molecules/{sample}.tsv.gz, plots/ (reads_vs_molecules/)
├── 05_allele_structure.html      05_allele_structure/allele_structure.xlsx, consensus/alleles.fasta, variants/alleles.vcf (.vcf.gz), bam/, plots/ (alleles/, clusters/, waterfall/)
├── 06_instability.html           06_instability/instability.xlsx, instability_metrics.tsv, distribution_summary.tsv, plots/ (timecourse/, pair_trend/, pair_correlation/, pair_scatter/)
├── 07_repeat_lengths.html        07_repeat_lengths/repeat_lengths.xlsx, plots/ (heatmap/, group_density/, group_histogram/, sample_ranges/, sample_card/; optional sample_density/, sample_histogram/)
├── 00_index.html                 start here: links every report, table and figure, by sample, range and group
├── run_summary.xlsx              one row per sample from every module, headline metrics, column and file guide (+ .tsv copies)
├── module_data/                  one .RData per module (used by resume)
└── temp/                         per-sample caches (remove_temp = TRUE deletes it)
```

Every module keeps its tables and data files at the top of its folder and its figures under `plots/`.

### Plots

Every saved figure is named **`{plot type}_{sample, group or measure}_{measure}_{range}_{unit}.png`**, with the parts that apply, so a file keeps its meaning wherever it is copied and sorts by plot type, then sample or group. Plot types with one file per sample, range, group or measure are saved in a folder named after the type (e.g. `05_allele_structure/plots/waterfall/waterfall_{sample}.png`); plot types with one file for the run sit in `plots/` itself. A type folder with many files is split into subfolders, so each holds one set to flick through: first by unit (`combined/`, `reads/`, `molecules/`), then, for types with one file per sample or group in each range, by range (`all-lengths/`, `range01/`, …), and for types with one file per range, by measure or tract pair. For example `07_repeat_lengths/plots/sample_density/combined/range02/sample_density_{sample}_CAG_range02_combined.png`, `06_instability/plots/timecourse/combined/CAG/timecourse_CAG_range02_combined.png`, `07_repeat_lengths/plots/sample_card/molecules/sample_card_{sample}_molecules.png`.

**Units.** A plot that can show both reads and molecules is drawn once for each version in `plot_units`: `combined` (reads grey and molecules teal in one figure, the default), `reads` or `molecules`; the version is the last part of the file name. `--plot_units combined,reads,molecules` gives all three; `--plot_units reads,molecules` gives the two separately instead of combined. Plots that cannot overlay the two (heatmaps, per-sample tract scatter) give a reads and a molecules file for `combined`. Without UMIs, only reads are drawn.

**Choosing plots.** `plots` sets which plot types are saved: `all` (the default: every type marked *yes* under *default* below), a list of types (`--plots timecourse,sample_ranges,waterfall`), or `all` with types left out or added (`--plots all,-waterfall,sample_density`). A per-sample, per-range or per-group type that is not selected is not drawn at all, which saves time on large runs; the reports still show their own figures. Changing `plots` or `plot_units` re-runs only the modules whose plots change. The older switches still work: `waterfall = FALSE`, `sample_cards = FALSE`, `visualise_alignment = FALSE` and `plot_per_sample = FALSE` leave the matching types out, and `waterfall_per_cluster = TRUE` and `plot_sample_ranges = TRUE` add theirs.

| Plot type | Module | Files | Reads / molecules versions | Default | What it shows |
|---|---|---|---|---|---|
| `reads_per_sample` | 1 | one for the run | – | yes | Reads per sample, one labelled row each, imported and kept |
| `read_length_distribution` | 1 | one for the run | – | yes | Read length distribution pooled and per sample |
| `quality_vs_length` | 1 | one for the run | – | yes | Mean base quality against read length per sample |
| `pacbio_np_rq` | 1 | one for the run | – | yes | PacBio pass count and read quality per sample |
| `pacbio_np_distribution` | 1 | one for the run | – | yes | PacBio pass count of every read: the share at each pass count, the cap marked, with the reads shorter and longer than the median length compared |
| `pacbio_np_vs_length` | 1 | one for the run | – | yes | PacBio pass count against read length: the median and the middle 50% and 80% of reads in each length band, and the share at the cap |
| `kinnex_segments` | 1 | one for the run | – | yes | Kinnex reads only: segments per concatenated read in this run, and pass count and read quality by where the segment starts in its concatenated read |
| `adapter_trimming` | 1 | one for the run | – | yes | Adapter trimming outcome |
| `read_accounting` | 2 | one for the run | – | yes | Reads kept at each alignment step |
| `flank_boundary_shortfall` | 2 | one for the run | – | yes | Flank alignments that stop short of the repeat |
| `mapping_quality` | 2 | one for the run | – | yes | Mapping quality per flank |
| `strand_by_step` | 2 | one for the run | – | yes | Share of forward-strand reads at each step |
| `alignment` | 2 | one per sample | – | yes | Each read's flank alignments (one per sample and reference) |
| `detection_rate` | 3 | one for the run | – | yes | Reads with each tract found |
| `flank_length` | 3 | one for the run | – | yes | Flank length between primer and tract per sample |
| `strand_by_length` | 3 | one for the run | – | yes | Strand balance by tract length |
| `count_method_agreement` | 3 | one for the run | – | yes | How often the counting methods differ by tract length |
| `count_method_gap` | 3 | one for the run | – | yes | How far the counting methods differ by tract length |
| `break_positions` | 3 | one for the run | – | yes | Where broken tracts are broken |
| `between_tract_sequences` | 3 | one for the run | – | yes | Sequences between adjacent tracts |
| `extraction_outcome` | 4 | one for the run | – | yes | UMI extraction outcome per sample |
| `umi_length` | 4 | one for the run | – | yes | UMI length per sample |
| `umi_neighbours` | 4 | one for the run | – | yes | UMIs one mismatch from another against chance |
| `length_rule_scorecard` | 4 | one for the run | – | yes | Molecule length rules scored |
| `length_rule_comparison` | 4 | one for the run | – | yes | Molecule length rules compared |
| `shortened_flag_by_length` | 4 | one for the run | – | yes | Molecules flagged as shortened by length |
| `shortened_flag_by_depth_band` | 4 | one for the run | – | yes | Molecules flagged as shortened by length and depth |
| `shortened_filter_distribution` | 4 | one for the run | – | yes | Shortened-read share per molecule |
| `shortened_flag_per_sample` | 4 | one for the run | – | yes | Molecules flagged as shortened per sample |
| `reads_per_molecule` | 4 | one for the run | – | yes | Reads per molecule |
| `molecules_per_sample` | 4 | one for the run | – | yes | Molecules found and kept per sample |
| `molecules_and_depth` | 4 | one for the run | – | yes | Molecules against read depth per sample |
| `within_molecule_error` | 4 | one for the run | – | yes | Read length minus molecule length by molecule length |
| `within_molecule_error_hist` | 4 | one for the run | – | yes | Read length minus molecule length (histogram) |
| `reads_far_from_molecule` | 4 | one for the run | – | yes | Reads far from their molecule's length per sample |
| `reads_vs_molecules_cohort` | 4 | one for the run | – | yes | Read and molecule length distributions for the run |
| `reads_vs_molecules` | 4 | one per sample | – | yes | Read and molecule length distributions (one per sample) |
| `allele_shares` | 5 | one for the run | – | yes | Each sample's reads and molecules split between its alleles |
| `clusters_overview` | 5 | one for the run | – | yes | Clusters in every sample |
| `clusters` | 5 | one per sample | – | yes | One sample's reads by cluster |
| `variants_by_sample` | 5 | one for the run | – | yes | Variants in alleles 1 and 2 of every sample |
| `alleles` | 5 | one per sample | – | yes | One sample's alleles and variants |
| `waterfall` | 5 | one per sample | – | yes | Reads stacked by length with each base coloured. With UMIs (clustering on molecules), each row is one molecule, ordered by its molecule length and drawn with its read closest to that length: one read's bases, not a consensus |
| `waterfall_cluster` | 5 | one per sample | – | no | Waterfall of each cluster of each sample |
| `experiment_map` | 6 | one for the run | – | yes | Samples by group and time |
| `counts_per_range` | 6 | one for the run | yes | yes | Reads and molecules in each range |
| `control_setpoints` | 6 | one for the run | yes | yes | Control samples and group setpoints |
| `timecourse` | 6 | one per range | yes | yes | Headline metrics against time (one per measure and range) |
| `instability_index_by_range` | 6 | one for the run | yes | yes | Instability index in every range |
| `pair_trend` | 6 | one per range | yes | yes | Change in one tract against change in another (one per range) |
| `pair_correlation` | 6 | one per range | yes | yes | Per-sample correlation between two tracts (one per range) |
| `pair_scatter` | 6 | one per sample | yes | yes | One sample's lengths of two tracts per read or molecule (every range) |
| `pair_map` | 6 | one per tract pair | yes | yes | Two tracts' lengths against each other for every read or molecule, all lengths, a panel per group |
| `pair_map_group` | 6 | one per group | yes | yes | One group's samples: two tracts' lengths against each other, a panel per sample |
| `heatmap` | 7 | one per measure | yes | yes | Every sample's length distribution (one per measure) |
| `group_density` | 7 | one per group | yes | yes | A group's samples as ridges (all lengths and each range) |
| `group_histogram` | 7 | one per group | yes | yes | A group's samples as histograms (all lengths and each range) |
| `sample_ranges` | 7 | one per sample | yes | yes | One sample's length in every range |
| `kinnex_position` | 7 | one for the run | – | yes | Kinnex reads only: the share of long segments (20+ units above the sample's median) by the segment's order within its concatenated read, and the repeat length by start position, with a test of whether the long segment comes first |
| `sample_card` | 7 | one per sample | yes | yes | One sample's every measure (all lengths and the ranges drawn in the report) |
| `sample_density` | 7 | one per sample | yes | no | One sample's ridge for all lengths and each range (one file each) |
| `sample_histogram` | 7 | one per sample | yes | no | One sample's histogram for all lengths and each range (one file each) |

Key tables:

- `run_summary.xlsx`: one row per sample with the results of every module (sheet `samples`), the headline metrics for every range and unit, and guides to its columns and to every output file.
- `06_instability/instability.xlsx`: `instability_metrics` (one row per sample × unit × measure × range × setpoint type; the column names used by v2 scripts are kept), `distribution_summary`, `modal_peaks`, `group_controls`, `pair_stats`, `range_table`.
- `03_repeat_detection/read_repeat_counts/{sample}.tsv.gz`: one row per read with every tract's counts and structure, the composites, and the v2 columns `repeat_count_full`, `repeat_count_match`, `repeat_count_tracts`.
- `02_alignment/alignment.xlsx`, sheet `read_accounting`: reads per sample at every step from input to analysed.

Logs go to `logs/<date_time>/` in the Duke folder, not the output folder, so runs from one Duke folder share it: the run log, the resolved parameters as R code and the resolved tracts table. The run log's first lines give the output folder. Warnings that would otherwise reach the end of the run (from report text, figure sizes and rendering; those inside the reports' code are silenced) are written to the log where they happen, each different message once with the module it came from (e.g. `Warning (Module 6): ...`), and the end of the run gives their count per module, with the full list in `warnings-<date_time>.tsv` beside the run log.

---

## Several repeat tracts in one run

The default is one tract: `rpt_pattern` (CAG) in `path_ref`. To measure more, give a tracts table (`path_tracts`, `.csv` or `.xlsx` with a sheet `tracts`), one row per tract in 5' to 3' order. Any number of tracts works.

| Column | Required | Meaning |
|---|---|---|
| `tract` | yes | Name used in every column, plot and file (letters, digits, `_`, `-`) |
| `motif` | yes | Repeat unit, e.g. `CAG` |
| `path_ref` | no | Reference whose N block holds this tract; blank = `path_ref` |
| `count_method`, `min_repeats`, `max_mismatch`, `start_perfect_repeats`, `end_perfect_repeats`, `max_gap`, `max_tract_gap`, `return_option`, `flank_pad`, `na_handling` | no | Per-tract overrides of the global `rpt_*` settings; blank = the global value |

**One reference or several.** Tracts that share a reference are found inside its one N block, each searched downstream of the one before; this suits adjacent tracts such as HTT CAG and CCG (`www/HTTset20_cag_ccg.fasta`), and Duke also reports the sequence between them per read (e.g. `CAACAGCCGCCA`), which identifies the interruption variant. Tracts far apart, or in different amplicons, get their own reference each; reads are then aligned to every reference and the tracts are joined by read name. Both can be mixed.

**Composite measures.** A combined length per read, set with `composites` (`--composites "CAG_CCG_sum=CAG+CCG:sum,CAG_CCG_span=CAG+CCG:span"`) or a `composites` sheet in the tracts workbook (`name`, `tracts`, `mode`):

- `sum` – the tracts' counts added, each counted by its `count_method` (`full` by default, so interruptions inside a tract are counted); in bp when the motifs differ in length;
- `span` – from the start of the first tract to the end of the last, divided by the motif length, so any sequence between the tracts is included (in bp when the motifs differ in length). For HTT, span = CAG + CCG + 4 for the standard allele.

**What you get for every tract and composite (each a "measure"):**

- Module 3: detection rate, counts by every method, modal by method, and the per-read export (`CAG_full`, `CAG_match`, `CAG_pure`, `CCG_full`, ..., `CAG_CCG_sum`, `CAG_CCG_span`).
- Module 6: distributions and the full set of instability metrics, per measure. A measure without its own `analysis_ranges_<measure>` column uses the primary tract's ranges: the reads are chosen by their CAG length (e.g. the expanded allele) and the metrics are calculated on that measure. With its own column the ranges are on its own scale. The primary range selects reads only; it does not cut the other measure, so in a CAG range of 40–50 a read with CAG 48 and CCG 10 has a sum of 58 and is kept. All lengths are in units of the measure's own motif, so with a CAG primary and a pentanucleotide tract, `analysis_ranges` is in CAG units and `analysis_ranges_<pentanucleotide tract>` in 5 bp units.
- Module 6, comparing tracts: for each pair (`measure_pairs`, default every pair of tracts), the mean change in tract B against the change in tract A from the sample median, per group, plus a per-read scatter and a Spearman correlation per sample. This shows whether CAG instability carries CCG with it.
- Module 7: a heatmap per measure, and a row per measure on each sample card.

The first tract is the primary measure: it is used for the read selection above, the legacy `repeat_count_*` columns, clustering and waterfalls. Example tables: `www/tracts_example.csv`, `www/tracts_example.xlsx`.

---

## UMIs

`umi = TRUE` turns on Module 4. The read layout assumed is `flank_5 – UMI – flank_3 – amplicon`, as in the custom PacBio UMI design (v4.07):

| Parameter | Default | Design element |
|---|---|---|
| `umi_flank_5` | `GTAAAACGACGGCCAGT` | M13-F |
| `umi_flank_3` | `GGCGACCCTGGAAAAGCTGATGAA` | ciosiF+1 |
| `umi_len_min`, `umi_len_max` | 12, 16 | N14 UMI, allowing for indels |

Set these to match your own oligos. Reads are oriented to the reference first, so the flanks are searched on the forward strand.

**Steps.** Each flank is found by exact match, then with up to `umi_flank_max_mismatch` edits (indels allowed); the UMI is the sequence between them. UMIs are clustered with the umi_tools directional method (a UMI joins a neighbour within `umi_cluster_max_dist` edits that has at least 2n – 1 of its reads). Each UMI family is one molecule. Molecules with fewer than `umi_min_reads` reads are dropped. Each molecule's length is `top_cluster` by default: counting down from the longest read, the first window of lengths 2 × `umi_top_cluster_w` + `umi_top_cluster_wfrac` × length units deep that holds at least `umi_top_cluster_k` reads and at least `umi_top_cluster_min_share` (default 20%) of the molecule's reads, and the median of the reads in it, so a few stray long reads cannot set the length of a large molecule. PCR slippage mostly shortens reads, so the longest well-supported length is the best estimate of the template. A length halfway between two whole units goes to the one more of the molecule's reads sit on. Finally, a molecule where at least `umi_shortened_min_share` (20%) of the reads are at least `umi_shortened_min_units` (20) units longer than the molecule itself is flagged and, with `umi_shortened_filter` (on by default), dropped: such a molecule most likely came from a longer template whose first copies were shortened in PCR, so its length is too short. This is judged from the molecule's own reads, with no alleles involved.

**Why the shortened-molecule filter (`umi_shortened_filter`).** In the 2026 jasmine run (66 samples, a 14-CAG and a ~130-CAG allele), 64% of the molecules called between the two alleles held at least 20% of their reads 20 or more units longer than the molecule, against about 4% of molecules at either allele; 65% of those long reads sat within 15 units of the sample's expanded peak. Almost all of these molecules had only 3–5 reads (79% flagged at 3–5 reads, 1–2% at 11 or more), so they are expanded templates sampled too few times for the full-length reads to be found. The filter cut the share of molecules above 30 CAG that sat between the alleles from a median of 14.7% to 5.7% per sample, close to the 3.9–4.0% seen when `umi_min_reads` was raised to 6 or 11 instead, while keeping the shallow molecules that show no long reads. Its cost is that about 11% of genuine 3–5-read molecules are also flagged, because one stray read is already 20% of them; this is similar at every length, so it does not bias the distribution. Set `umi_shortened_filter = FALSE` to keep these molecules; Module 4 always shows what the filter removes.

**Why repeat length is not used to form molecules.** Molecules are formed from the UMI sequences alone, and nothing about alleles is assumed. PCR slippage and deletion can make reads of one molecule differ by tens of units (a 125-CAG molecule can give reads of 10 CAG or fewer), which is exactly the error collapsing is meant to remove, so length cannot say which reads belong together. Two molecules sharing a UMI by chance is rare (Module 4 reports the expected rate), and PCR chimeras are a minority of a family's reads. Earlier v3 versions split a UMI family whose reads fell on both sides of an allele boundary; this was removed in v3.0.4 because it turned slipped copies of one molecule into a second, artificial molecule and assumed every sample has one short and one long allele.

**Report.** Extraction outcome per sample (found, no 5' flank, no 3' flank, wrong length), UMI length, reads per molecule, UMI neighbours against the number expected by chance (a check that clustering merges errors, not distinct UMIs), within-molecule read error against molecule length for every measure, the share of each sample's reads more than 10 units from their molecule (a measure of PCR damage), molecules holding reads much longer than themselves, a comparison of the six length rules (`top_cluster` by reads alone and by reads and share, `max_supported`, `mode`, `median` and `mean`: tests of each, and where they differ), and reads against molecules.

**Downstream.** Modules 6 and 7 report every distribution and metric for reads and for molecules (`unit` column). Each measure is collapsed separately, so a molecule's `CAG_CCG_sum` is the collapsed sum, not the sum of the collapsed CAG and CCG. `04_umi/molecules/{sample}.tsv.gz` lists every molecule.

### UMI parameters (Module 4)

<!-- params:UMI:start -->

| Parameter | Default | Format | What it does |
|---|---|---|---|
| `umi` | `FALSE` | TRUE / FALSE | Extract UMIs and collapse reads to molecules |
| `umi_flank_5` | `GTAAAACGACGGCCAGT` | text | Sequence immediately 5' of the UMI (default M13-F) |
| `umi_flank_3` | `GGCGACCCTGGAAAAGCTGATGAA` | text | Sequence immediately 3' of the UMI (default ciosiF+1) |
| `umi_flank_max_mismatch` | `3` | whole number | Edits allowed when locating each UMI flank (indel tolerant) |
| `umi_search_window` | `250` | whole number | The 5' flank must start within this many bases of the read start |
| `umi_len_min` | `12` | whole number | Shortest UMI accepted |
| `umi_len_max` | `16` | whole number | Longest UMI accepted |
| `umi_cluster_max_dist` | `1` | whole number | Edit distance for directional UMI clustering (0 = exact match; 1 or 2) |
| `umi_min_reads` | `3` | whole number | Minimum reads for a molecule to be kept |
| `umi_collapse_stat` | `top_cluster` | text | Rule that turns the reads of a molecule into its one length: `top_cluster`, `max_supported`, `mode`, `median` or `mean` (all are compared in the Module 4 report) |
| `umi_top_cluster_k` | `2` | whole number | top_cluster: reads required in the window |
| `umi_top_cluster_w` | `3` | number | top_cluster: the window reaches 2 × w + wfrac × length repeats below its top read; w is the constant part |
| `umi_top_cluster_wfrac` | `0.1` | number | top_cluster: the length-proportional part of the window (see umi_top_cluster_w) |
| `umi_top_cluster_min_share` | `0.2` | number | top_cluster: share of the molecule's reads the window must also hold, so a few stray long reads cannot set the length (0 = only umi_top_cluster_k reads needed) |
| `umi_shortened_filter` | `TRUE` | TRUE / FALSE | Drop molecules where many reads are much longer than the molecule (most likely a longer template whose reads were shortened in PCR): the shortened-molecule filter |
| `umi_shortened_min_units` | `20` | number | Shortened-molecule filter: how many units longer than its molecule a read must be to count as long |
| `umi_shortened_min_share` | `0.2` | number | Shortened-molecule filter: share of the reads of a molecule that must be long for the molecule to be flagged as shortened |

<!-- params:UMI:end -->

---

## How repeats are counted

Duke counts the **repeat tract itself**: the run of CAG units, without the interruptions and sequence that follow it. For a 47Q HTT allele (`(CAG)n CAA CAG CCG CCA (CCG)7`): standard 45 CAG, DOI (duplication of the interruption, CAA CAG CAA CAG) 43 CAG, LOI (loss of interruption, CAA → CAG) 47 CAG.

A tract is a run of motif matches, with gaps of up to `rpt_max_gap` bases inside a run and `rpt_max_tract_gap` between runs merged. It must start and end with `rpt_start_perfect_repeats` / `rpt_end_perfect_repeats` perfect units, which trims e.g. a trailing `CAACAG` off the CAG tract. The search covers the sequence between the two flank alignments plus `rpt_flank_pad` bases (30) either side, because the aligner sometimes places the last unit inside a flank (the cause of the 1-unit under-count of LOI alleles in v2).

Three counts per tract:

| Method | Definition |
|---|---|
| `full` (default) | tract length in bp ÷ motif length |
| `match` | exact motif matches in the tract |
| `pure` | the longest uninterrupted run |

**Why `full` is the default.** A single wrong base inside the tract breaks it into two runs. `full` runs across it (the gap is within `rpt_max_gap`), `match` loses one unit, and `pure` drops to the longer of the two pieces. Sequencing errors inside long repeats are common enough for this to matter: in the 2026 jasmine MSH3ko/MLH1ko UMI run (PacBio HiFi, read quality about Q42 overall), in one sample 6% of reads with a 14-CAG tract but 58% of reads with a ~130-CAG tract had a broken tract, with breaks at random positions, so the median `pure` count of the long allele was 88 against 132 for `full` and 131 for `match`. `pure` is only the right count where the question is a real interruption, and Module 3 shows which situation you are in.

Module 3 compares the methods read by read, so no assumption about alleles or controls is needed: the share of reads where `match` or `pure` differs from `full` at each tract length (per sample and pooled), each count against `full` for every read, and where along the tract the breaks fall (flat = sequencing errors; a sharp peak = a real interruption). The tables are in `repeat_detection.xlsx` (`method_agreement`, `method_by_length`, `modal_by_method`).

### Repeat detection parameters (Module 3)

These settings apply to every tract unless the tracts table overrides them for one tract. The tracts-table column is the parameter name without `rpt_` (`flank_pad`, `max_gap`, ...), except `count_method` for `repeat_count_method` and `na_handling` for `na_repeat_handling`.

<!-- params:Repeats:start -->

| Parameter | Default | Format | What it does |
|---|---|---|---|
| `rpt_pattern` | `CAG` | text | Repeat motif (used when there is no tracts table) |
| `rpt_tract_name` | – | text | Name for the single tract (default: the motif) |
| `rpt_min_repeats` | `2` | whole number | Minimum tract length in repeat units |
| `rpt_max_mismatch` | `0` | whole number | Mismatches allowed per repeat unit |
| `rpt_start_perfect_repeats` | `2` | whole number | Perfect units required at the start of a tract |
| `rpt_end_perfect_repeats` | `2` | whole number | Perfect units required at the end of a tract |
| `rpt_max_gap` | `6` | whole number | Maximum bases between units inside a tract |
| `rpt_max_tract_gap` | `18` | whole number | Maximum gap between sub-tracts before they are merged |
| `rpt_return_option` | `longest` | text | Which tract to keep when several are found: `longest` or `first` |
| `rpt_flank_pad` | `30` | whole number | Bases either side of the flank boundaries searched for the tract, so a tract can extend past where the aligner placed the flank |
| `repeat_count_method` | `full` | text | Primary count: `full` (tract length / motif length), `match` (exact units) or `pure` (longest uninterrupted run) |
| `na_repeat_handling` | `keep` | text | Reads that align but have no tract: `keep` (no length; counted in the detection rate only), `convert_to_zero` (length 0, as in v2; falls into ranges starting at 0 or NA) or `filter` (removed) |
| `composites` | – | text list | Combined measures, e.g. `CAG_CCG_sum=CAG+CCG:sum` (each tract's count added, under its count_method) or `CAG_CCG_span=CAG+CCG:span` (first tract start to last tract end) |
| `export_read_counts` | `TRUE` | TRUE / FALSE | Write per-read lengths to 03_repeat_detection/read_repeat_counts/ |
| `rm_flank_length_outliers` | `FALSE` | TRUE / FALSE | Remove reads with outlying combined flank length |
| `flank_iqr_multiplier` | `1.5` | number | IQR multiplier for the flank length filter |

<!-- params:Repeats:end -->

---

## Instability metrics

For each sample × measure × range × unit, with setpoint s and d = length – s per read:

| Metric | Definition |
|---|---|
| `instability_index` (II) | mean(d) = EI + CI. With a modal setpoint this is the instability index of Lee et al. (2010); their threshold of 10% of the modal peak is `floor = 0.1` in the settings file |
| `expansion_index` (EI) | sum of positive d ÷ n |
| `contraction_index` (CI) | sum of negative d ÷ n |
| `mean_change_above_setpoint` | mean d of the reads above s: how large a typical expansion is (`mean_expansion_size` before v3.1; "expansion index" in v2) |
| `mean_change_below_setpoint` | mean d of the reads below s (`mean_contraction_size` before v3.1; "contraction index" in v2) |
| `share_at_setpoint`, `share_above_setpoint`, `share_below_setpoint` | share within ±0.5 of s, above, below (`prop_at`, `prop_above`, `prop_below` before v3.1); `n_at_setpoint` etc. are the counts |
| `share_kplus_above_setpoint`, `share_kplus_below_setpoint` | share at least k units above / below s, e.g. `share_5plus_above_setpoint` (`expansion_thresholds`, default 5 and 10; `prop_above_k` before v3.1) |
| `expansion_ratio`, `contraction_ratio` | reads above / at s, below / at s (unstable when few reads sit at s) |
| `modal_minus_setpoint`, `median_minus_setpoint`, `mean_minus_setpoint` | the sample's modal, median and mean minus s (`modal_diff` etc. before v3.1) |
| `z_score` | (sample value – control setpoint) ÷ SD across the control samples; blank with manual controls or a single control |

Every column, with its units, is also listed in the `definitions` sheet of `instability.xlsx`.

Two setpoints are reported (`setpoint_type`):

- `sample_relative` – the sample's own `control_setpoint_metric` (median by default; `modal_length` gives Lee's II);
- `control_relative` – that metric aggregated (`control_aggregation_method`) across the group's control samples, so the metrics measure change since the control (e.g. t = 0). Where the settings file gives a `manual_control_repeat_length` for a sample and range, it replaces the group control there (a blank cell uses the group control); `setpoint_source` records which was used, and Module 6 lists the samples where a manual length was used.

Distribution statistics (modal, median, mean, SD, IQR, 5th and 95th percentiles, skewness, kurtosis) are in `distribution_summary`.

### Instability parameters (Module 6)

<!-- params:Instability:start -->

| Parameter | Default | Format | What it does |
|---|---|---|---|
| `range_peak_span` | `3` | whole number | Window for modal peak detection |
| `range_min_n` | `50` | whole number | Fewest reads (or molecules) a sample needs in a range for that range's statistics and metrics to be calculated; below it they are left blank, and a control sample is left out of the setpoint |
| `group_control` | `TRUE` | TRUE / FALSE | Calculate group control setpoints |
| `control_sample_selection` | `flagged` | text | `flagged`, `earliest` or `all` |
| `control_setpoint_metric` | `median_length` | text | `modal_length`, `mean_length` or `median_length` |
| `control_aggregation_method` | `median` | text | `mean`, `median` or `trimmed_mean` across control samples |
| `expansion_thresholds` | `5, 10` | whole-number list | k for proportion of reads at least k units above / below the setpoint |
| `cohort_plot_by` | `group, time` | text list | Dimensions for cohort plots: group, time |
| `time_label` | `Time (days)` | text | Axis label for the time column of the settings file, with its unit |
| `headline_metrics` | `modal_length, median_length, mean_length, instability_index, expansion_index, contraction_index, share_5plus_above_setpoint, share_5plus_below_setpoint` | text list | Metrics drawn as timecourses by group |
| `headline_setpoint_type` | `control_relative` | text | Setpoint type for headline plots: `control_relative` or `sample_relative` |
| `plot_ranges` | – | text list | Ranges drawn inline (default: the first two ranges of the primary measure) |
| `measure_pairs` | – | text list | Pairs of measures to plot against each other per read, e.g. `CAG:CCG` (default: every pair of tracts) |
| `pair_max_points` | `4000` | whole number | Maximum reads per sample in per-read pair plots |

<!-- params:Instability:end -->

---

## Parameters

Every parameter is listed in a table with its default, format and what it does. The tables for repeat detection, UMIs, instability and paths are in those sections; the rest are below.

| Parameters | Table |
|---|---|
| Paths (`dir_data`, `path_ref`, `path_settings`, ...) | [Path parameters](#path-parameters) |
| Module 1: import and QC | [Import parameters](#import-parameters-module-1) |
| Module 2: alignment | [Alignment parameters](#alignment-parameters-module-2) |
| Module 3: repeat detection (`rpt_*`, counting, composites) | [Repeat detection parameters](#repeat-detection-parameters-module-3) |
| Module 4: UMI (`umi_*`) | [UMI parameters](#umi-parameters-module-4) |
| Module 5: allele structure (clustering, waterfalls) | [Allele structure parameters](#allele-structure-parameters-module-5) |
| Module 6: instability (controls, metrics, plots) | [Instability parameters](#instability-parameters-module-6) |
| Module 7: repeat lengths | [Repeat length parameters](#repeat-length-parameters-module-7) |
| Threads, resume, modules, plots | [Runtime parameters](#runtime-parameters) |

**Setting a parameter.** On the command line, `--name value`: `--rpt_flank_pad 40`, `--umi TRUE`. A list is comma-separated, `--expansion_thresholds 5,10,20`, and a whole-number list can use ranges, `--run_modules 1:3,6`. In `scripts/duke_run.R` or `duke_run_local.R`, `name = value` (`rpt_flank_pad = 40L`, `expansion_thresholds = c(5, 10, 20)`). `NA` means not set (for a path) or no limit (for a filter); `–` means empty. `./duke --help` shows the same information, and `--dry_run` prints the values a run would use.

The tables are generated from `lib/params.R`, which defines every parameter; after changing it, run `Rscript scripts/update_readme_params.R` to update them.

### Most often changed

| Parameter | Default | |
|---|---|---|
| `threads` | 4 | CPU cores |
| `run_modules` | 1:7 | e.g. `--run_modules 6,7` |
| `rpt_pattern` | CAG | motif for a single tract |
| `path_tracts`, `composites` | – | several tracts and combined measures |
| `repeat_count_method` | full | full, match or pure |
| `umi` | FALSE | Module 4 |
| `cluster` | FALSE | allele clustering, consensus and variants (Module 5) |
| `waterfall` | TRUE | waterfall plots (Module 5) |
| `trim` | FALSE | adapter trimming; demultiplexed HiFi reads are already trimmed |
| `check_duplicate_readnames` | FALSE | PacBio and ONT read names are unique |
| `np_min` | NA | minimum PacBio passes |
| `na_repeat_handling` | keep | reads with no tract: keep (no length), convert_to_zero (length 0, as in v2) or filter |
| `headline_metrics` | modal, median, mean, II, EI, CI, share_5plus_above_setpoint, share_5plus_below_setpoint | timecourse rows in Module 6 |
| `plot_ranges` | first two ranges | ranges drawn in Module 6 and on sample cards |
| `sample_cards` | TRUE | one figure per sample in Module 7 |

Parameters from v2 that no longer exist are ignored with a message; the v2 count names (`repeat_count_full` etc.) are still accepted.

### Import parameters (Module 1)

<!-- params:Import:start -->

| Parameter | Default | Format | What it does |
|---|---|---|---|
| `import_patterns` | `\.fastq$, \.fastq\.gz$, \.fq$, \.fq\.gz$, \.fasta$, \.fa$, \.bam$` | text list | Regex patterns for files to import |
| `import_recursive` | `TRUE` | TRUE / FALSE | Search dir_data recursively |
| `r1_pattern` | `_R1_` | text | Pattern identifying R1 files in paired-end data |
| `r2_pattern` | `_R2_` | text | Pattern identifying R2 files in paired-end data |
| `select_one_of_pair` | `R1` | text | Paired-end: import `R1`, `R2` or `all` |
| `downsample` | `NA` | whole number | Maximum reads per sample (NA = all) |
| `np_min` | `NA` | whole number | Drop PacBio HiFi reads with fewer passes than this (NA = no filter) |
| `pacbio_scatter_downsample` | `50000` | whole number | Not used since v3.3.2 (the pass-count plot now counts every read); kept so older parameter files still run |
| `check_duplicate_readnames` | `FALSE` | TRUE / FALSE | Remove duplicate read names, keeping the longest. Off by default: PacBio and ONT read names are unique, and the check removed 0 reads in 36 M |
| `trim` | `FALSE` | TRUE / FALSE | Trim adapters from read ends. Off by default: demultiplexed HiFi reads are already trimmed (20 of 36 M reads were trimmed in the msn run) |
| `trim_max_mismatch` | `3` | whole number | Maximum errors in an adapter match |
| `trim_with_indels` | `TRUE` | TRUE / FALSE | Allow indels in adapter matching |

<!-- params:Import:end -->

### Alignment parameters (Module 2)

<!-- params:Alignment:start -->

| Parameter | Default | Format | What it does |
|---|---|---|---|
| `minimap2_args` | `-t 1 -x sr -w 1` | text | minimap2 arguments for aligning reads to each flank |
| `keep_bam` | `FALSE` | TRUE / FALSE | Keep the per-sample flank BAMs in 02_alignment/bam/ for IGV |
| `visualise_alignment` | `TRUE` | TRUE / FALSE | Write per-sample flank alignment plots |
| `visualise_alignment_downsample` | `1000` | whole number | Maximum reads per alignment plot |

<!-- params:Alignment:end -->

### Allele structure parameters (Module 5)

`cluster_by = haplotype` is described in [Common workflows](#common-workflows).

<!-- params:Allele structure:start -->

| Parameter | Default | Format | What it does |
|---|---|---|---|
| `cluster` | `FALSE` | TRUE / FALSE | Cluster reads into alleles, build consensus sequences and call variants |
| `cluster_by` | `repeat` | text list | `repeat`, `haplotype`, or `repeat,haplotype` |
| `cluster_unit` | `auto` | text | What the repeat lengths are clustered on: `molecules` (the length of each kept UMI molecule, and every read takes the cluster of its molecule), `reads`, or `auto` (molecules when umi = TRUE, otherwise reads) |
| `repeat_cluster_max` | `6` | whole number | Maximum repeat-length clusters per sample |
| `repeat_cluster_min_prop` | `0.05` | number | Clusters holding fewer than this proportion of reads (or molecules) are merged into their nearest neighbour, unless far from every other cluster (more than 10 units and a quarter of the length) and holding at least 1% and 20: such a cluster is a minor allele and is kept |
| `repeat_cluster_downsample` | `10000` | whole number | Reads used to choose the number of repeat clusters |
| `haplotype_cluster_max` | `6` | whole number | Maximum haplotype clusters per sample |
| `haplotype_cluster_downsample` | `3000` | whole number | Reads used to define haplotype clusters |
| `haplotype_region` | `both` | text | Flank(s) used for haplotypes: `left`, `right` or `both` |
| `haplotype_trim_length` | `auto` | text | `auto`, a number of bases, or `NA` for no trimming |
| `consensus_threshold` | `0.5` | number | Minimum agreement for a consensus base |
| `consensus_downsample` | `50` | whole number | Maximum reads per consensus alignment |
| `consensus_parts` | `FALSE` | TRUE / FALSE | Also write each allele's left flank, repeat region and right flank as `alleles_left.fasta`, `alleles_mid.fasta` and `alleles_right.fasta` beside `alleles.fasta` |
| `call_variants` | `TRUE` | TRUE / FALSE | Call variants in each cluster consensus against the reference |
| `allele_minor_share` | `0.2` | number | Also report a substitution carried by at least this share of an allele's reads when it is not in the allele's consensus (NA = off) |
| `allele_support_reads` | `5000` | whole number | Maximum reads per allele used to count each variant's read support |
| `allele_bam` | `TRUE` | TRUE / FALSE | Write a BAM per sample of reads aligned to the whole reference, tagged by allele (HP) |
| `allele_bam_reads` | `2000` | whole number | Maximum reads per allele in the allele BAMs |
| `allele_bam_method` | `auto` | text | How reads are placed on the reference in the allele BAMs: `flanks` aligns each read's left and right flank to the reference flanks and writes the repeat region as one insertion at the N block (for short amplicons, which minimap2 leaves unmapped); `minimap2` aligns whole reads; `auto` uses flanks when the reads' flanks total under 1,000 bp |
| `allele_bam_minimap2` | `-x map-hifi` | text | minimap2 options for the allele BAMs when minimap2 is used (e.g. -x map-ont for nanopore reads) |
| `waterfall` | `TRUE` | TRUE / FALSE | Write per-sample waterfall plots |
| `waterfall_downsample` | `1000` | whole number | Maximum reads per waterfall |
| `waterfall_flank_bp` | `60` | whole number | Bases of flank drawn either side of the tract in waterfalls |
| `waterfall_per_cluster` | `FALSE` | TRUE / FALSE | Also write one waterfall per cluster |

<!-- params:Allele structure:end -->

### Repeat length parameters (Module 7)

<!-- params:Repeat lengths:start -->

| Parameter | Default | Format | What it does |
|---|---|---|---|
| `sample_cards` | `TRUE` | TRUE / FALSE | Write one summary figure per sample |
| `repeat_histogram_binwidth` | `1` | whole number | Histogram bin width in repeat units |
| `plot_every_range` | `TRUE` | TRUE / FALSE | Save the by-group distribution figures for every analysis range of the primary measure (FALSE: only plot_ranges, the ones shown in the report) |
| `plot_sample_ranges` | `FALSE` | TRUE / FALSE | Also draw a density and a histogram figure for every sample, for all lengths and for each range drawn (many figures on a large run) |
| `plot_dots_max` | `2000` | whole number | Reads (or molecules) drawn as dots per sample in the density (ridge) figures |
| `plot_dots_unit` | `molecules` | text | In density figures that show reads and molecules together (plot_units = combined): whether the dots are the `molecules` (default) or the `reads`. Figures of reads alone or molecules alone always use their own unit |

<!-- params:Repeat lengths:end -->

### Runtime parameters

<!-- params:Runtime:start -->

| Parameter | Default | Format | What it does |
|---|---|---|---|
| `reuse_from` | – | text | An earlier run's output folder to take module results from, so those modules are not re-run. A module is reused only when its results were made with the same settings as this run (its own and every upstream module's); it keeps that run's report and figures |
| `reuse_modules` | – | whole-number list | With reuse_from: the modules to reuse, e.g. 1:4 (default: every module made by this Duke version whose results still match; a module made by an earlier version is reused only when named here). A module named here that cannot be reused stops the run before it starts |
| `redraw` | – | whole-number list | Modules whose report, tables and figures are drawn again from their saved results (module_data/), without importing or recomputing, e.g. 1,4 after a Duke update that changed only their figures, or with new plots / plot_units. Modules 1, 4 and 7 can be redrawn (Module 7 only draws, so it simply re-runs). The saved results must have been made with the current settings; they can come from reuse_from |
| `plots` | `all` | text list | Plot types to save: `all` (default; every type except sample_density, sample_histogram and waterfall_cluster), a list of types, or `all` with `-type` to leave types out, e.g. `all,-waterfall,sample_density`. The README lists every type |
| `plot_units` | `combined` | text list | For plots that can show reads and molecules: `combined` (both in one figure, default), `reads`, `molecules`, or a list, e.g. `combined,reads,molecules` for all three |
| `plot_dpi` | `150` | whole number | Resolution of saved plots |
| `plot_per_sample` | `TRUE` | TRUE / FALSE | Write the per-sample plots (one figure per sample) in Modules 2, 4, 5, 6 and 7 |
| `threads` | `4` | whole number | CPU cores |
| `resume` | `TRUE` | TRUE / FALSE | Reuse a module's results when the parameters that affect it are unchanged |
| `free_memory` | `TRUE` | TRUE / FALSE | Drop large objects from memory once a module has finished with them (at present the imported read sequences, after Module 2 has aligned them). Changes no results; called remove_intermediate before v3.3.2 |
| `remove_temp` | `FALSE` | TRUE / FALSE | Delete temp/ at the end of the run. temp/ holds per-sample working caches (UMI collapsing, clustering) that let a re-run after a crash pick up part-way through a module; a re-run with changed settings adds new caches beside the old ones, so delete temp/ (or set this TRUE) once a run is final |
| `run_modules` | `1:7` | whole-number list | Modules to run |

<!-- params:Runtime:end -->

---

## Resume

With `resume = TRUE` (default) a module is reused when its stored parameter hash matches. Each hash covers the parameters that affect that module and the hashes of the modules it depends on, plus the reference and settings files and the list of input files (names, sizes, dates), so changing `rpt_flank_pad` re-runs Modules 3 to 7 but changing a plot setting in Module 7 re-runs Module 7 only, and adding or editing the settings file re-runs Modules 6 and 7 only (Module 4 uses it to group one figure, which keeps the groups it had when Module 4 last ran). Changes to Duke's own code are not detected within one folder: after updating Duke, write to a new output folder (recommended; use `reuse_from` to take unchanged modules from the old one) or use `--resume FALSE`. If you run only some modules after changing a parameter that affects an earlier one (e.g. `--run_modules 6` after changing a tract setting), Duke stops and names the module to include, rather than building on out-of-date results.

### Reusing modules from an earlier run (`reuse_from`, `reuse_modules`)

To write a new run to a new folder while taking some modules from an earlier run, rather than re-running them, give the earlier run's output folder:

```
--dir_out .../result_duke_v3.3.2_new --reuse_from .../result_duke_v3.3.1_old --reuse_modules 1:4
```

Before any module runs, Duke copies each module to reuse from the earlier folder into the new one: its saved data (`module_data/`), its output folder, its report and its working caches (`temp/`). The saved data, the largest part, are hard-linked rather than copied when both folders are on the same file system (as on Myriad's Scratch): the new folder gets its own entry for the same file, which takes seconds instead of many minutes and no extra space, and the file stays in the new folder if the earlier folder is deleted (deleting the earlier folder therefore frees that space only once every folder linking to it is deleted). It then skips those modules as it would on resume, and runs the rest on top of them. The earlier folder is only read, never changed: a module that runs again deletes its linked file first and writes a new one.

- A module is reused only when its results were made with the same settings as this run: its own and every upstream module's, the same input files, reference and settings file. So a reused module holds the same data a re-run would give. The plot choice (`plots`, `plot_units`) is not part of that check: a reused module keeps the earlier run's figures.
- **`reuse_modules` not given:** every module whose results still match **and** that was made by this Duke version is reused. A module made by an earlier version is not, because the newer version may have changed its report or figures; the log says so for each module.
- **`reuse_modules` given** (e.g. `1:4`): those modules are reused even if an earlier Duke version made them (their reports and figures are then that version's, and the log says which version). If a named module cannot be reused (it is not in the earlier folder, or its settings differ), Duke stops before copying or running anything and says why, so a long job does not start on a wrong assumption.
- The log lists each module as reused, or not reused and why. Each module's `module_data/*.hash` file records the Duke version that made it, the version that last redrew its report and figures (after a redraw), and, for a reused module, the folder it came from.

Typical use: a new Duke version changes only the later modules (e.g. clustering in Module 5), and Modules 1–4 (import, alignment, repeat detection, UMI collapsing, the slow part) can be taken from the previous run.

### Redrawing a module's report and figures (`redraw`)

A reused or resumed module keeps the report and figures it was made with. When a newer Duke version changes only a module's figures or report, or you want other `plots` / `plot_units`, `redraw` draws them again from the module's saved results, without importing or recomputing:

```
--dir_out .../result_new --reuse_from .../result_old --reuse_modules 1:4 --redraw 1,4
```

- **Module 1** loads its saved reads and tables instead of importing the files again (and skips the np, duplicate-name, downsampling and trimming steps, whose results are saved). **Module 4** loads its saved molecules instead of extracting UMIs and forming molecules again; its length-rule tests, error and distribution figures are recalculated from them, which takes a few minutes on a large run. **Module 7** only draws, so redrawing it simply re-runs it. Modules 2, 3, 5 and 6 compute their results inside the report and cannot be redrawn; leave them out of `redraw` to re-run them.
- The saved results must have been made with the current settings (the same check as resume). They can be this folder's own or come from `reuse_from` in the same run. If they do not match, Duke stops before drawing and says so.
- The saved results (`module_data/`) are not changed, so the modules after a redrawn one stay current; Module 4's per-sample `molecules/` tables are kept as they are. The report opens with a note that it was redrawn, and the module's `.hash` file records the version that made the results and the version that redrew them.
- The tables are the same as the run that made the results. Figures with randomly jittered points can place the points slightly differently.
- A redrawn module's code is the current version's, so a redraw is for changes to figures and text: if a Duke update changes how a module calculates its results, re-run it instead (the CHANGELOG says which).

---

## Running on Myriad

Duke runs on one node, with `--threads` set to the cores the job asks for. Install it first ([Installation](#installation), *On Myriad*).

**1. Make a job script** for each run, from `scripts/duke_myriad.sh`, and keep it with your other job scripts. Edit:

- `-N`: the job name, e.g. `duke-cagccg-337` (change it for each new version or run, so the queue and logs say which is which);
- `-wd`: your Duke folder (`~/Scratch/bin/duke`), so the scheduler's logs go to its `logs/`;
- `-M`: your email address, for the start, end and abort emails;
- the resources: `-pe smp` (cores), `-l mem` (memory **per core** on Myriad) and `-l h_rt` (time). The template asks for 12 cores × 64 GB; a run of tens of samples needs far less (6 cores and `mem=32G`), and a run that only rebuilds reports from reused modules needs 1 core and 8 GB;
- the `export PATH=...minimap2` line, if minimap2 is not in `~/Scratch/bin/minimap2`;
- the paths after `./duke`: `--dir_data`, a new `--dir_out` for each run, `--path_ref`, `--path_settings`, and any other options (the script lists common ones as comments).

**2. Submit and follow it.**

```bash
qsub /path/to/your_job_script.sh
qstat -u $USER
tail -f ~/Scratch/bin/duke/logs/duke_<job id>.out
```

Duke's progress goes to the scheduler's `logs/duke_<job id>.out` and to the run log, `logs/<date_time>/duke_run-<date_time>.log`; R's messages go to the run log only; warnings go to both, labelled with their module, and are counted at the end of the run. Loading R on a compute node can take a few minutes before the first Duke line appears, and new lines can take a minute to show on the login node. Module 1 writes nothing to `dir_out` until it has read every file. `qstat -j <id> | grep usage` shows CPU time and memory so far, `qdel <id>` stops a job and `qacct -j <id>` gives the resources a finished job used. A waiting job's resources can be changed with `qalter`, e.g. `qalter -pe smp 1 -l mem=8G <id>`.

**3. Copy the results back**, e.g. from your own computer: `rsync -av <user>@myriad.rc.ucl.ac.uk:/path/to/result_duke_v3.3.13 ~/results/`. Add `--exclude temp/ --exclude module_data/` to leave out the caches, and `--exclude '*.bam*'` to leave out the BAMs. `00_index.html` works on the copy.

**Other clusters.** `scripts/duke_kathleen.sh` (SGE, `-pe mpi`) and `scripts/duke_kathleen_slurm.sh` (Slurm, `sbatch`; monitor with `squeue -u $USER`, `scancel <id>`, `sacct -j <id>`) ask for one whole node, as Duke cannot spread one run across nodes. A run that fits on Myriad is simplest there.

**Run time.** On 2 cores, 14 PacBio samples of 2,500 reads run Modules 1 to 7 with clustering in about 3 minutes. On a full run (Myriad, 6 cores, 66 PacBio samples, 21 M reads, 2026 jasmine UMI run), Module 1 took 28 min, Module 2 35 min to 1 h 12 min (depending on what else shared the node) and Module 3 about 1 h, against about 1.5 h, 6.5 h and 7.5 h for the same data in v2.3.1. Run time grows roughly in proportion to the number of reads; trimming and clustering, when switched on, add to it. Module 1 holds the reads in memory: about 12.5 GB for 21 M reads.

---

## Common workflows

**First look at new data** – run without a settings file; every module except 6 runs (4 only with `umi = TRUE`, 5 with `cluster = TRUE` or `waterfall = TRUE`), and the Module 3 and 7 plots show where to put the ranges (see [Quick start](#quick-start)).

**Add or change the settings** – the same command with `--path_settings new.xlsx`; Modules 1–5 are reused and 6 and 7 re-run. `--run_modules 6,7` does the same without checking the earlier modules.

**New figures after a Duke update** – a new output folder with `--reuse_from` the old one and `--reuse_modules` the modules whose results are unchanged; add `--redraw 1,4` to draw Modules 1 and 4 again from those results with the new version's figures, rather than importing and collapsing again ([Redrawing](#redrawing-a-modules-report-and-figures-redraw)).

**CAG and CCG in one run** – `--path_ref www/HTTset20_cag_ccg.fasta --path_tracts www/tracts_example.csv --composites "CAG_CCG_sum=CAG+CCG:sum,CAG_CCG_span=CAG+CCG:span"`; add `analysis_ranges_CCG` to the settings for CCG ranges on the CCG scale.

**UMI run** – `--umi TRUE`, with `--umi_flank_5 / --umi_flank_3 / --umi_len_min / --umi_len_max` if your design differs from the default.

**Allele structure** – `--cluster TRUE` for clusters, consensus sequences and VCFs. Add `--cluster_by repeat,haplotype` to split each length cluster by flank haplotype: Duke finds the 15-mers that mark real SNPs and indels (present in 10–90% of reads, and in blocks that always occur together), scores each read on them and clusters the scores. It finds haplotypes carried by more than about 10–20% of reads, and reads that combine the two alleles' flank variants (PCR chimeras) come out as their own cluster.

**Combine runs done separately in v2** – no longer needed: one run with a tracts table joins the tracts by read name.

---

## Troubleshooting

- **"Reference must contain exactly one N block"** – mask the repeat with one run of N; for separate tracts use separate references in the tracts table.
- **A tract is detected in few reads** (Module 3 detection plot) – check the motif and that the N block covers the tract; lower `rpt_start_perfect_repeats` / `rpt_end_perfect_repeats` for very interrupted tracts.
- **Few UMIs extracted** – the Module 4 outcome plot says which flank is missing; check the flank sequences and orientation, and `umi_search_window`.
- **"Settings file not found"** – `path_settings` was given but the file is not there; Duke stops before running anything. `--dry_run` checks every input path.
- **Module 6 skipped** – no settings file was given (`path_settings` unset).
- **"None of the ... file names in the settings file match"** – `file_name` must be the input file name, with or without its extension; the message shows examples from both sides.
- **"No analysis ranges could be read from the settings file"** – the matching rows have no `analysis_ranges`.
- **Results did not change after updating Duke** – resume reused modules made by the old version; write to a new output folder, or use `--resume FALSE`.
- **"Module N needs Module M results made with the current parameters"** – a parameter of an earlier module changed; add that module to `run_modules` (or run every module: resume reuses what is still current).
- **"... must be ..." at the start** – a setting with fixed choices is misspelt; the message lists the choices.
- **"umi_min_reads must exceed umi_top_cluster_k"** – with exactly k reads `top_cluster` returns their median; raise `umi_min_reads` or change the statistic.

---

## Citation and contact

Hölbling, B.V. et al. A multimodal screening platform for endogenous dipeptide repeat proteins in C9orf72 patient iPSC-neurons. *Cell Reports* (2025). https://www.sciencedirect.com/science/article/pii/S2211124725004668

Michael Flower, UCL Queen Square Institute of Neurology – michael.flower@ucl.ac.uk – https://github.com/mike-flower/duke
