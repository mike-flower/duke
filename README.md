# Duke

A modular pipeline for amplicon sequencing of tandem repeats: per-read repeat counts, optional UMI collapse, allele structure and somatic instability metrics, with an HTML report for every step.

**v3.0.7** – see [Changes from v2.3.1](#changes-from-v231)

---

## Contents

- [What Duke does](#what-duke-does)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Input files](#input-files)
- [Several repeat tracts in one run](#several-repeat-tracts-in-one-run)
- [UMIs](#umis)
- [How repeats are counted](#how-repeats-are-counted)
- [Instability metrics](#instability-metrics)
- [Output files](#output-files)
- [Parameters](#parameters)
- [Resume](#resume)
- [HPC deployment](#hpc-deployment)
- [Common workflows](#common-workflows)
- [Troubleshooting](#troubleshooting)
- [Changes from v2.3.1](#changes-from-v231)
- [Citation and contact](#citation-and-contact)

---

## What Duke does

| Module | Report | What it does | Runs when |
|---|---|---|---|
| 1. Import and QC | `01_import_qc.html` | Import BAM / FASTQ / FASTA; read counts and lengths; PacBio passes and quality; optional np filter, de-duplication, downsampling and adapter trimming | always |
| 2. Alignment | `02_alignment.html` | minimap2 to the left and right flanks of each reference; read accounting from input to analysed reads; flank boundary, MAPQ and strand QC | always |
| 3. Repeat detection | `03_repeat_detection.html` | Finds every tract in every read; counts by three methods; modal length by method; sequence between tracts; composite measures; per-read export | always |
| 4. UMI | `04_umi.html` | Extracts UMIs, clusters them into molecules and collapses each molecule's reads to one length | `umi = TRUE` |
| 5. Allele structure | `05_allele_structure.html` | Clustering into alleles, consensus sequences, flank variants (VCF), waterfall plots | `cluster = TRUE` or `waterfall = TRUE` |
| 6. Instability | `06_instability.html` | Distributions and instability metrics per sample, range and measure, for reads and molecules; group controls; timecourses; tract-against-tract plots | a settings file is given |
| 7. Repeat lengths | `07_repeat_lengths.html` | Every repeat-length figure in one place: cohort heatmaps, ridges by group and one card per sample | always |

Reads are grey and molecules teal in every figure. Every report has a timing table and an output manifest, and every module writes its tables to one Excel file.

---

## Installation

```r
install.packages(c("data.table", "dplyr", "tidyr", "tibble", "ggplot2", "ggridges", "patchwork", "scales",
                   "DT", "rmarkdown", "knitr", "openxlsx", "readxl", "mclust", "cluster", "pracma",
                   "digest", "pbapply", "pbmcapply"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Biostrings", "ShortRead", "Rsamtools", "GenomicAlignments", "DECIPHER", "pwalign"))
```

`pwalign` exists from Bioconductor 3.19; on older versions Duke uses the aligner in `Biostrings`. minimap2, samtools and pandoc must be on the `PATH`.

**HPC:** `export R_LIBS_USER=~/R/library` and load the R and samtools modules (the job scripts do this).

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
│   ├── update_readme_params.R  # rewrites the parameter tables in this README from lib/params.R
│   ├── duke_myriad.sh, duke_kathleen.sh, duke_kathleen_slurm.sh
└── www/                      # example reference, tracts table and settings files
```

---

## Quick start

```bash
./duke --dir_data ~/data/run1 --dir_out ~/data/run1/result_duke \
       --path_ref www/HTTset20_cagcaa.fasta --path_settings ~/data/run1/settings_duke.xlsx --threads 8
```

`./duke --help` lists every option by module with its default. Add `--dry_run` to any command to print the resolved parameters and tracts without running anything.

From R or RStudio, edit the paths in `scripts/duke_run_local.R` and source it. Every parameter is listed there with its default and a one-line explanation.

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
| `analysis_ranges` | yes | `[0-35][36-NA]` | Ranges on the primary measure (the first tract in the tracts table, or the single tract from `rpt_pattern`), in units of its motif, one bracket each; `NA` = unbounded |
| `floor` | no | `[3][3]` | Remove lengths with fewer reads than this (a value below 1 is a fraction of the modal bin); 1 keeps all |
| `max_peaks` | no | `[2][1]` | Modal peaks reported per range |
| `group` | no | `QS5` | Group for controls, timecourses and plots |
| `time` | no | `0`, `49` | Time point |
| `group_control_sample` | no | `TRUE` | The sample is a control for its group |
| `exclude` | no | `TRUE` | Leave the sample out of Module 6 |
| `manual_control_repeat_length` | no | `[14][120]` | Control setpoint per range; used only with `use_manual_controls = TRUE` |
| `analysis_ranges_<measure>` (and `floor_<measure>`, `max_peaks_<measure>`, `manual_control_repeat_length_<measure>`) | no | `[NA-8][9-NA]` | Ranges for another tract or composite on its own scale (see below) |

Rows whose file is not in `dir_data` are ignored, so one platemap-wide settings file can serve every run; if a file appears in more than one row, the first row is used. Group and range names may contain any characters; those that are unsafe in file names are replaced with `-` in output paths only. Examples: `www/settings_example.xlsx`, `www/settings_example_two_tracts.xlsx`.

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

- `sum` – the tracts' counts added (pure repeat content only; in bp when the motifs differ in length);
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

**Steps.** Each flank is found by exact match, then with up to `umi_flank_max_mismatch` edits (indels allowed); the UMI is the sequence between them. UMIs are clustered with the umi_tools directional method (a UMI joins a neighbour within `umi_cluster_max_dist` edits that has at least 2n – 1 of its reads). Each UMI family is one molecule. Molecules with fewer than `umi_min_reads` reads are dropped. Each molecule's length is `top_cluster` by default: the median of the longest window of reads of half-width `w` + `wfrac` × length that holds at least `umi_top_cluster_k` reads and at least `umi_top_cluster_min_share` (default 20%) of the molecule's reads, so a few stray long reads cannot set the length of a large molecule. PCR slippage mostly shortens reads, so the longest well-supported length is the best estimate of the template. A length halfway between two whole units goes to the one more of the molecule's reads sit on. Finally, a molecule where at least `umi_long_read_share` (20%) of the reads are at least `umi_long_read_units` (20) units longer than the molecule itself is flagged and, with `umi_long_read_filter` (on by default), dropped: such a molecule most likely came from a longer template whose first copies were shortened in PCR, so its length is too short. This is judged from the molecule's own reads, with no alleles involved.

**Why the long-read filter.** In the 2026 jasmine run (66 samples, a 14-CAG and a ~130-CAG allele), 64% of the molecules called between the two alleles held at least 20% of their reads 20 or more units longer than the molecule, against about 4% of molecules at either allele; 65% of those long reads sat within 15 units of the sample's expanded peak. Almost all of these molecules had only 3–5 reads (79% flagged at 3–5 reads, 1–2% at 11 or more), so they are expanded templates sampled too few times for the full-length reads to be found. The filter cut the share of molecules above 30 CAG that sat between the alleles from a median of 14.7% to 5.7% per sample, close to the 3.9–4.0% seen when `umi_min_reads` was raised to 6 or 11 instead, while keeping the shallow molecules that show no long reads. Its cost is that about 11% of genuine 3–5-read molecules are also flagged, because one stray read is already 20% of them; this is similar at every length, so it does not bias the distribution. Set `umi_long_read_filter = FALSE` to keep these molecules; Module 4 always shows what the filter removes.

**Why repeat length is not used to form molecules.** Molecules are formed from the UMI sequences alone, and nothing about alleles is assumed. PCR slippage and deletion can make reads of one molecule differ by tens of units (a 125-CAG molecule can give reads of 10 CAG or fewer), which is exactly the error collapsing is meant to remove, so length cannot say which reads belong together. Two molecules sharing a UMI by chance is rare (Module 4 reports the expected rate), and PCR chimeras are a minority of a family's reads. Earlier v3 versions split a UMI family whose reads fell on both sides of an allele boundary; this was removed in v3.0.4 because it turned slipped copies of one molecule into a second, artificial molecule and assumed every sample has one short and one long allele.

**Report.** Extraction outcome per sample (found, no 5' flank, no 3' flank, wrong length), UMI length, reads per molecule, UMI neighbours against the number expected by chance (a check that clustering merges errors, not distinct UMIs), within-molecule read error against molecule length for every measure, the share of each sample's reads more than 10 units from their molecule (a measure of PCR damage), molecules holding reads much longer than themselves, a comparison of the six length rules (tests of each, and where they differ), and reads against molecules.

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
| `umi_top_cluster_w` | `3` | number | top_cluster: constant part of the window half-width (repeats) |
| `umi_top_cluster_wfrac` | `0.1` | number | top_cluster: length-proportional part of the window half-width |
| `umi_top_cluster_min_share` | `0.2` | number | top_cluster: share of the molecule's reads the window must also hold, so a few stray long reads cannot set the length (0 = only umi_top_cluster_k reads needed) |
| `umi_long_read_filter` | `TRUE` | TRUE / FALSE | Drop molecules where many reads are much longer than the molecule (most likely a longer template whose reads were shortened in PCR) |
| `umi_long_read_units` | `20` | number | Long-read filter: how many units longer than its molecule a read must be to count as long |
| `umi_long_read_share` | `0.2` | number | Long-read filter: share of the reads of a molecule that must be long for it to be flagged |

<!-- params:UMI:end -->

---

## How repeats are counted

Duke counts the **repeat tract itself**: the run of CAG units, without the interruptions and sequence that follow it. For a 47Q HTT allele (`(CAG)n CAA CAG CCG CCA (CCG)7`): standard 45 CAG, DOI (loss of CAA CAG duplicated) 43 CAG, LOI (loss of interruption, CAA → CAG) 47 CAG.

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
| `composites` | – | text list | Combined measures, e.g. `CAG_CCG_sum=CAG+CCG:sum` (pure units added) or `CAG_CCG_span=CAG+CCG:span` (first tract start to last tract end) |
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
| `mean_expansion_size` | mean d of the reads above s (was "expansion index" in v2) |
| `mean_contraction_size` | mean d of the reads below s (was "contraction index" in v2) |
| `prop_at`, `prop_above`, `prop_below` | proportion within ±0.5 of s, above, below |
| `prop_above_k`, `prop_below_k` | proportion at least k units above / below s (`expansion_thresholds`, default 5 and 10) |
| `expansion_ratio`, `contraction_ratio` | reads above / at s, below / at s (unstable when few reads sit at s) |
| `modal_diff`, `median_diff`, `mean_diff` | the sample's modal, median and mean minus the group control's |
| `z_score` | (sample value – control setpoint) ÷ SD across the control samples; blank with manual controls or a single control |

Two setpoints are reported (`setpoint_type`):

- `sample_relative` – the sample's own `control_setpoint_metric` (median by default; `modal_length` gives Lee's II);
- `control_relative` – that metric aggregated (`control_aggregation_method`) across the group's control samples, so the metrics measure change since the control (e.g. t = 0). With `use_manual_controls = TRUE` a `manual_control_repeat_length` replaces it; `setpoint_source` records which was used.

Distribution statistics (modal, median, mean, SD, IQR, 5th and 95th percentiles, skewness, kurtosis) are in `distribution_summary`.

### Instability parameters (Module 6)

<!-- params:Instability:start -->

| Parameter | Default | Format | What it does |
|---|---|---|---|
| `range_peak_span` | `3` | whole number | Window for modal peak detection |
| `group_control` | `TRUE` | TRUE / FALSE | Calculate group control setpoints |
| `control_sample_selection` | `flagged` | text | `flagged`, `earliest` or `all` |
| `control_setpoint_metric` | `median_length` | text | `modal_length`, `mean_length` or `median_length` |
| `control_aggregation_method` | `median` | text | `mean`, `median` or `trimmed_mean` across control samples |
| `use_manual_controls` | `FALSE` | TRUE / FALSE | Use manual_control_repeat_length from the settings file as the control setpoint (replacing the group control) for those samples and ranges |
| `expansion_thresholds` | `5, 10` | whole-number list | k for proportion of reads at least k units above / below the setpoint |
| `cohort_plot_by` | `group, time` | text list | Dimensions for cohort plots: group, time |
| `headline_metrics` | `modal_length, median_length, mean_length, instability_index, expansion_index, prop_above_5` | text list | Metrics drawn as timecourses by group |
| `headline_setpoint_type` | `control_relative` | text | Setpoint type for headline plots: `control_relative` or `sample_relative` |
| `plot_ranges` | – | text list | Ranges drawn inline (default: the first two ranges of the primary measure) |
| `measure_pairs` | – | text list | Pairs of measures to plot against each other per read, e.g. `CAG:CCG` (default: every pair of tracts) |
| `pair_max_points` | `4000` | whole number | Maximum reads per sample in per-read pair plots |

<!-- params:Instability:end -->

---

## Output files

```
result_duke/
├── 01_import_qc.html             01_import_qc/qc.xlsx, plots/
├── 02_alignment.html             02_alignment/alignment.xlsx (read_accounting), plots/, bam/ (keep_bam)
├── 03_repeat_detection.html      03_repeat_detection/repeat_detection.xlsx, read_repeat_counts/{sample}.tsv.gz, plots/
├── 04_umi.html                   04_umi/umi.xlsx, molecules/{sample}.tsv.gz, plots/
├── 05_allele_structure.html      05_allele_structure/allele_structure.xlsx, consensus/, variants/*.vcf, plots/waterfall/
├── 06_instability.html           06_instability/instability.xlsx, plots/ (timecourse/, pairs/)
├── 07_repeat_lengths.html        07_repeat_lengths/repeat_lengths.xlsx, plots/ (heatmaps, ridges/, sample_cards/)
├── module_data/                  one .RData per module (used by resume)
└── temp/                         per-sample caches (remove_temp = TRUE deletes it)
```

Key tables:

- `06_instability/instability.xlsx`: `instability_metrics` (one row per sample × unit × measure × range × setpoint type; the column names used by v2 scripts are kept), `distribution_summary`, `modal_peaks`, `group_controls`, `pair_stats`, `range_table`.
- `03_repeat_detection/read_repeat_counts/{sample}.tsv.gz`: one row per read with every tract's counts and structure, the composites, and the v2 columns `repeat_count_full`, `repeat_count_match`, `repeat_count_tracts`.
- `02_alignment/alignment.xlsx`, sheet `read_accounting`: reads per sample at every step from input to analysed.

Logs go to `logs/<date_time>/` in the Duke folder: the run log, the resolved parameters as R code and the resolved tracts table.

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
| `use_manual_controls` | FALSE | settings-file control lengths replace the group control |
| `headline_metrics` | modal, median, mean, II, EI, prop_above_5 | timecourse rows in Module 6 |
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
| `pacbio_scatter_downsample` | `50000` | whole number | Maximum points in the PacBio passes-vs-length plot |
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
| `repeat_cluster_max` | `6` | whole number | Maximum repeat-length clusters per sample |
| `repeat_cluster_min_prop` | `0.05` | number | Clusters holding fewer than this proportion of reads are merged into their nearest neighbour |
| `repeat_cluster_downsample` | `10000` | whole number | Reads used to choose the number of repeat clusters |
| `haplotype_cluster_max` | `6` | whole number | Maximum haplotype clusters per sample |
| `haplotype_cluster_downsample` | `3000` | whole number | Reads used to define haplotype clusters |
| `haplotype_region` | `both` | text | Flank(s) used for haplotypes: `left`, `right` or `both` |
| `haplotype_trim_length` | `auto` | text | `auto`, a number of bases, or `NA` for no trimming |
| `consensus_threshold` | `0.5` | number | Minimum agreement for a consensus base |
| `consensus_downsample` | `50` | whole number | Maximum reads per consensus alignment |
| `call_variants` | `TRUE` | TRUE / FALSE | Call variants in each cluster consensus against the reference |
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

<!-- params:Repeat lengths:end -->

### Runtime parameters

<!-- params:Runtime:start -->

| Parameter | Default | Format | What it does |
|---|---|---|---|
| `plot_dpi` | `150` | whole number | Resolution of saved plots |
| `plot_per_sample` | `TRUE` | TRUE / FALSE | Write per-sample plots in Modules 2 and 3 |
| `threads` | `4` | whole number | CPU cores |
| `resume` | `TRUE` | TRUE / FALSE | Reuse a module's results when the parameters that affect it are unchanged |
| `remove_intermediate` | `TRUE` | TRUE / FALSE | Free memory between steps |
| `remove_temp` | `FALSE` | TRUE / FALSE | Delete temp/ (per-sample caches) at the end |
| `run_modules` | `1:7` | whole-number list | Modules to run |

<!-- params:Runtime:end -->

---

## Resume

With `resume = TRUE` (default) a module is reused when its stored parameter hash matches. Each hash covers the parameters that affect that module and the hashes of the modules it depends on, plus the reference and settings files and the list of input files (names, sizes, dates), so changing `rpt_flank_pad` re-runs Modules 3 to 7 but changing a plot setting in Module 7 re-runs Module 7 only. Changes to Duke's own code are not detected: after updating Duke, delete `module_data/` (or use `--resume FALSE`) to re-run everything. If you run only some modules after changing a parameter that affects an earlier one (e.g. `--run_modules 6` after changing a tract setting), Duke stops and names the module to include, rather than building on out-of-date results.

---

## HPC deployment

| Cluster | Best for | Scheduler | Script |
|---|---|---|---|
| Myriad | up to ~500 samples | SGE (`-pe smp`) | `scripts/duke_myriad.sh` |
| Old Kathleen | 1000+ samples | SGE (`-pe mpi`) | `scripts/duke_kathleen.sh` |
| New Kathleen | 1000+ samples | Slurm | `scripts/duke_kathleen_slurm.sh` |

Edit the paths in the script, then `qsub scripts/duke_myriad.sh` (or `sbatch` for Slurm). `--threads` is set from `$NSLOTS` / `$SLURM_NTASKS`. Each script lists the optional additions (tracts, UMIs, clustering, trimming) as comments.

Monitoring: `qstat -u $USER`, `qstat -j <id> | grep usage` (CPU time and data read so far), `qdel <id>`, `qacct -j <id>` (SGE); `squeue -u $USER`, `scancel <id>`, `sacct -j <id>` (Slurm). Duke's progress goes to the scheduler's `logs/duke_<id>.out` and to the run log, `logs/<date_time>/duke_run-<date_time>.log`; R's warnings and messages go to the run log only. Loading R on a compute node can take a few minutes before the first Duke line appears, and on shared filesystems new lines can take a minute to show on the login node. Module 1 writes nothing to `dir_out` until it has read every file.

**Run time.** On 2 cores, 14 PacBio samples of 2,500 reads run Modules 1 to 7 with clustering in about 3 minutes. On a full run (Myriad, 6 cores, 66 PacBio samples, 21 M reads, 2026 jasmine UMI run), Module 1 took 28 min, Module 2 35 min to 1 h 12 min (depending on what else shared the node) and Module 3 about 1 h, against about 1.5 h, 6.5 h and 7.5 h for the same data in v2.3.1. Run time grows roughly in proportion to the number of reads; trimming and clustering, when switched on, add to it. Module 1 holds the reads in memory: about 12.5 GB for 21 M reads.

---

## Common workflows

**First look at new data** – run without a settings file; every module except 6 runs (4 only with `umi = TRUE`), and the Module 3 and 7 plots show where to put the ranges.

**Add or change the settings** – `./duke ... --path_settings new.xlsx --run_modules 6,7` (Modules 1–3 are reused).

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
- **Results did not change after updating Duke** – resume reused old modules; delete `module_data/`.
- **"umi_min_reads must exceed umi_top_cluster_k"** – with exactly k reads `top_cluster` returns their median; raise `umi_min_reads` or change the statistic.

---

## Changes from v2.3.1

**New**
- Several tracts in one run (any number), on one shared reference or separate references, with the sequence between adjacent tracts reported per read, `sum` and `span` composite measures, per-tract metrics, and tract-against-tract plots and correlations.
- Optional UMI module (Module 4): flank-based extraction with user-set flanks (default M13-F / ciosiF+1), directional clustering, `top_cluster` collapse; every Module 6 and 7 result for reads and for molecules.
- Modal, median and mean length by every counting method (Module 3).
- Read accounting from input to analysed reads (Module 2), carried onto every sample card.
- Sample cards: one figure per sample with every measure; cohort heatmaps; timecourses by group for the headline metrics.
- `--dry_run`; `./duke --help` grouped by module; logs per run with the resolved parameters and tracts.

**Changed**
- Counts are pure units with a 30 bp search pad beyond the flank boundaries: LOI alleles are no longer under-counted by 1. Expect +1 in most LOI reads and in about 1% of standard-allele reads (reads with a sequencing error that turns CAA into CAG).
- `expansion_index` and `contraction_index` are now the sums of positive and negative changes divided by n, so they add up to II; the v2 quantities are now `mean_expansion_size` and `mean_contraction_size`. `prop_above_k` / `prop_below_k` added.
- Manual control lengths are opt-in (`use_manual_controls`); in v2 they were applied whenever present.
- Defaults: `trim`, `check_duplicate_readnames` and `cluster` now FALSE. `na_repeat_handling` is now `keep`: reads that span both flanks but have no tract count against the detection rate but get no length, instead of a length of 0 that fell into ranges starting at 0 or NA (`convert_to_zero` restores v2).
- Resume refuses to run a module on upstream results made with other parameters (it names the module to add to `run_modules`), and adding or replacing input files re-runs Module 1.
- Modules reorganised: 4 (allele calling) and 5 (waterfalls) are now Module 5 (allele structure); Module 6 is instability; all repeat-length figures are in Module 7.
- Module 3 is vectorised (15 s against several minutes on the test set); Module 7 writes about one figure per sample instead of ~18,000 PNGs.
- Resume uses parameter hashes per module, so a changed parameter re-runs exactly the modules it affects.

**Fixed after the first v3 build**
- Haplotype clustering compared the full ~2.5 kb flanks of every pair of reads by edit distance (as in v2), which takes hours per sample. It now uses variant k-mers: 14 PacBio samples in 1 min 20 s, and 20,000 reads in under 10 s. Tested with planted variants: haplotypes in 20–50% of reads are found with 98–100% accuracy, and on real data the haplotype split matches the CAG allele in 87–99% of reads.
- BAMs are read with samtools instead of Rsamtools, so the "[samopen] no @SQ lines in the header" notice no longer appears for every PacBio BAM. The reads, qualities and PacBio tags are identical.
- The reports no longer load MathJax, which removes pandoc's `--mathjax` warning.

**v3.0.1**
- Module 1 reads each BAM once: the PacBio tags (np, rq, ec, zm, bq) are kept from the same pass as the sequences instead of being read in a second pass over the whole file. Results are identical; the BAM reading step takes about a third less time.

**v3.0.2**
- UMI clustering (Module 4) assigned UMIs to their groups by name, which slows down with the square of the number of distinct UMIs: 63,000 UMIs took 100 s and a sample with 300,000 about 40 minutes. It now works by position: 315,000 UMIs in 24 s, with identical groups.

**v3.0.3**
- Module 3 compares the counting methods read by read: the share of reads where `match` or `pure` differs from `full` at each tract length (per sample and pooled), each count against `full` for every read, and where along the tract the breaks fall, with `method_by_length` added to `repeat_detection.xlsx`. See [Why `full` is the default](#how-repeats-are-counted).
- UMI allele boundary (Module 4): where the two alleles are far apart, the gap between them is nearly empty and the old rule placed the boundary anywhere in it; in the 2026 jasmine run it fell between 20 and 106 CAG depending on the sample, and where it fell inside the lower tail of the expanded allele (about 10 of 66 samples) it split up to 7,900 UMI families per sample that belonged to one molecule. The boundary is now the first point above the shorter allele where the density reaches the floor of the dip (about 20 CAG for a 14/135 sample), found with a smoothing bandwidth of at least one unit so that the gaps between whole-number lengths are not mistaken for it. This changes the Module 4 molecules for those samples. The read-against-molecule plots show the boundary.
- Module 4 report: a new section, *Families with reads on both alleles*, shows how many UMI families are mixed (per sample, as families and as reads), how large the minority side is, the lengths on each side of every mixed family, and where the minority reads fall, so you can see whether they come from template switching (minority reads at the other allele's peak) or from large slippage or deletion (spread between the alleles), and whether splitting (`umi_split_mixed`, on by default) suits the data. Per-sample counts are in `umi.xlsx` (`mixed_families`).
- Module 4 report: UMI length as one labelled bar per length (share and read count, accepted and rejected); reads per molecule as one bar per count from 1 to 30, for molecules and for the reads they hold, with the median, commonest and middle 50% of kept molecules, plus molecules against depth per sample; UMI neighbours as the observed against the expected share per sample, with how the expected share is calculated; within-molecule read error against molecule length for every measure (tracts and combined measures) instead of a short/long split; the read-against-molecule distribution as one whole-run summary with one example sample (every sample is still saved).
- Module 5: the clustering and consensus sections say "Not run" when `cluster = FALSE` instead of describing results that are absent; the waterfall y-axis marks distinct lengths instead of read quantiles (which repeated the dominant allele's length).
- Module 7: heatmaps on a log colour scale from 0.01% with the top 0.1% of lengths cut, so minor alleles and expansion tails are visible; ridge plots add one panel per plotted range and stop at the 99.9th percentile; samples are labelled by the part of the file name that differs between them.
- Each module clears its own `plots/` folder before writing, so plots from an earlier run cannot be mistaken for new ones.
- A `path_settings` that does not exist stops the run at the start (before, Module 6 was skipped without an error), and a settings file whose names match none of the samples says so. `--dry_run` checks that every input path exists.
- Reports: means are written beside each box (as percentages where the axis is a percentage) instead of over the points; read quality is shown as QV; the read-length histogram has a log count axis; the detection plot shows reads without a tract; the UMI length plot labels only bars with at least 0.1%; Module 1 plots mean rather than median read length against quality; the flank-boundary columns in `alignment.xlsx` are named `prop_left_short` and `prop_right_short`.

**v3.0.4**
- UMI molecules (Module 4) no longer use repeat length: each UMI family is one molecule. The per-sample allele boundary, the labelling of reads as short or long and the splitting of families with reads on both sides are removed, with `umi_split_mixed` and `umi_allele_min_sep` (both are ignored with a message if set). Splitting turned slipped or deleted copies of one molecule into a second, artificial molecule, and assumed one short and one long allele per sample, which fails for controls, multi-clone samples and multi-tract runs. See [Why repeat length is not used to form molecules](#umis).
- Module 4 report: the *Families with reads on both alleles* section is replaced by the share of each sample's reads more than 10 units shorter or longer than their own molecule (`reads_far_from_molecule` in `umi.xlsx`), and a comparison of molecule length under four length rules (`length_rule_comparison`). New parameter `umi_top_cluster_min_share` (default 0, the previous behaviour): the `top_cluster` window must also hold this share of the molecule's reads, so a few stray long reads cannot set the length of a large molecule.
- A molecule length halfway between two whole units (the median of an even number of reads) now goes to the unit more of the molecule's reads sit on (at random if level). Before, plots rounded halves to the even number, which gave molecule histograms a saw-tooth.
- Module 4 report: reads per molecule shows 1–10 reads singly, then bands up to more than 500; the UMI neighbour plot's lines now join the points they belong to; the UMI length axis stops at 100%.
- Module 3: break positions are placed at random within their unit, so short tracts no longer give a comb of spikes that looks like interruptions; the reads-without-a-tract plot gives enough decimal places for very small rates and states the count; the disagreement plot draws each length bin at its centre; the *Sequence between tracts* section failed with more than one tract and is fixed.
- Module 6: timecourse titles give the range as written (`<=35`) rather than `NA-35`. Module 7: sample-card panels stop at the 99.9th percentile of the reads; molecule heatmaps say "molecules", not "reads", for the lengths not drawn.

**v3.0.5**
- `umi_top_cluster_min_share` is now 0.2 by default: the `top_cluster` window must hold at least 2 reads **and** 20% of the molecule's reads. In the 2026 jasmine run this changed 0.5% of molecules (3,502 of 675,633): 1,329 molecules with a median of 40 reads, most of them at 14 CAG, had been called at the expanded allele because of a few long reads (PCR chimeras) and are now called short; 572 molecules called above the expanded allele were brought down to it; and 837 molecules called at the expanded allele by 2–3 reads, whose other reads sat between the alleles, are now called between. Set it to 0 for the v3.0.4 behaviour.
- `umi_collapse_stat` can be `top_cluster`, `max_supported` (the longest exact length held by `umi_top_cluster_k` reads), `mode`, `median` or `mean`. Whatever is chosen, Module 4 calculates every rule for every molecule and reports them side by side in *Choosing a molecule's length*: three tests of each rule made on the molecules' own reads with no assumption about alleles (agreement between random halves of a molecule's reads, dependence on read depth, and how many molecules move when stray reads from another molecule are added, as a PCR chimera would), and, at each molecule length, how many molecules each other rule would make shorter or longer (`length_rule_scorecard` and `length_rule_comparison` in `umi.xlsx`). The tests follow the scorecard of the July 2026 development script, without its allele labels.
- Module 4 report: a new section, *Molecules holding much longer reads*, gives the share of molecules at each length that hold reads at least 20 units longer than the molecule itself (`molecules_with_long_reads` in `umi.xlsx`), as a check for early-cycle PCR deletions creating artificial short molecules. The within-molecule error distribution shows the full range on a log axis with the shares of reads shorter than, equal to and longer than their molecule. Reads per molecule gives each bar's share.
- Module 3 report: the counting-method disagreement plot adds the curve expected from a single error chance per base, fitted to the data, and gives that chance; the flank-length plot splits reads into exact, 1, 2–5, 6–10, 11–50 and more than 50 bases off the commonest flank length.
- Module 6 report: the metric definitions are a full-width table; *Group controls* explains how the setpoints are calculated (separately for reads and molecules, from the chosen control samples); *Length by counting method* is removed, as Module 3 covers it read by read (and the `count_methods` sheet with it); *Measures against each other* is renamed *Comparing tracts*, with a description.
- Every Module 3 and 4 plot now opens with *How to read this* (what the axes are, in plain terms) and *What it shows*.

**v3.0.6**
- A worker process that is killed during a parallel step (usually for running out of memory) now stops the run with a message. Before, the samples that worker held were dropped from that module without any error.
- Report text works for any data: no wording assumes two alleles, a normal and an expanded allele, or a particular length. Module 3's *Counts against each other* explains the **pure** wedge with a worked example, draws the half-way line on the **pure** panel, and adds a sentence calculated from the run; *Where the tract is broken* explains the dips at the tract ends from `rpt_start_perfect_repeats` / `rpt_end_perfect_repeats`; the modal-by-method note no longer points to the Module 6 section removed in v3.0.5.
- Module 4: the second **top_cluster** rule in *Choosing a molecule's length* uses the run's `umi_top_cluster_min_share` (20% when that is 0) instead of a fixed 20%, and is marked "(used)" for any share.
- Module 4 report: *Reads per molecule* uses a log axis with even bins and marks the median and mean of kept molecules (with the middle 50% and 80% in the subtitle); new per-sample box plots of molecules kept, median reads per kept molecule and share of reads used (`molecules_per_sample.png`). *Error distribution* now shows, for each distance *k*, the share of reads at least *k* units shorter or longer than their molecule, which reads more easily than the histogram. *How well each rule behaves* says for each test whether lower or closer to 0 is better, draws the rule used in orange, and marks the best rule in each test; its table names each column's direction and has a "used" column. *Where the rules differ* names the rule used as the reference (it has no line of its own), gives each rule a fixed, distinct colour, and its table marks the reference row.

**v3.0.7**
- HTT reference: `www/HTTset20_cagcaa.fasta` (N block over the CAG tract and the CAA CAG after it) is added and is now the recommended CAG reference; `www/HTTset20.fasta` is kept for comparison with earlier runs.
- A `sum` composite of tracts whose motifs differ in length is now in bp, as `span` already was; before, units of different sizes were added.
- Module 6 explains that a primary range selects reads but does not cut other measures, and that lengths are in each measure's own motif units.
- New long-read filter for UMI molecules (`umi_long_read_filter`, on by default; `umi_long_read_units` = 20, `umi_long_read_share` = 0.2): a molecule where at least 20% of the reads are at least 20 units longer than the molecule is dropped, as it most likely came from a longer template whose reads were shortened in PCR. Judged from the molecule's own reads; no alleles involved. See [Why the long-read filter](#umis) for the jasmine evidence.
- Module 4 report: *Molecules built from a longer template* replaces *Molecules holding much longer reads*: the share of molecules flagged at each length, by reads per molecule; the molecule length distribution with and without the filter, and with `umi_min_reads` raised to 6 and 11 for comparison; the share flagged per sample. *Reads per molecule* colours molecules kept, dropped for too few reads and dropped by the filter. `umi.xlsx` sheets `long_read_flag_by_length` and `long_read_flag_per_sample`; the molecule table carries each molecule's primary length, long reads, their share and the flags.
- Module 5 fixes, found on a run with CAG and CCG tracts: the tract consensus is now built from reads at the cluster's modal length, and terminal gaps count in every consensus, so a consensus stops where most reads stop (before, it ran to the end of the longest sampled read, e.g. ~590 bp for a 136-unit allele). Flank variants are called with free end gaps on both sequences, so sequence between the primary tract and the reference flank (the rest of an N block that covers several tracts) is no longer forced onto the flank; on that run this removed ~2,350 false calls (~19 per cluster), leaving 20. Reads with no primary tract are left out of the clusters (they formed 1–2 read clusters) and counted in the text. *Clusters per sample* is a bar chart.

**Removed**
- Per-sample × per-range histogram, scatter and density PNGs (replaced by the sample cards and heatmaps). Parameters for them are ignored with a message.

---

## Citation and contact

Hölbling, B.V. et al. A multimodal screening platform for endogenous dipeptide repeat proteins in C9orf72 patient iPSC-neurons. *Cell Reports* (2025). https://www.sciencedirect.com/science/article/pii/S2211124725004668

Michael Flower, UCL Queen Square Institute of Neurology – michael.flower@ucl.ac.uk – https://github.com/mike-flower/duke
