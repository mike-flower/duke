# Duke

A modular pipeline for amplicon sequencing of tandem repeats: per-read repeat counts, optional UMI collapse, allele structure and somatic instability metrics, with an HTML report for every step.

**v3.0.0** – see [Changes from v2.3.1](#changes-from-v231)

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
| 4. UMI | `04_umi.html` | Extracts UMIs, clusters them, splits mixed families and collapses reads to molecules | `umi = TRUE` |
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
│   ├── duke_myriad.sh, duke_kathleen.sh, duke_kathleen_slurm.sh
└── www/                      # example reference, tracts table and settings files
```

---

## Quick start

```bash
./duke --dir_data ~/data/run1 --dir_out ~/data/run1/result_duke \
       --path_ref www/HTTset20.fasta --path_settings ~/data/run1/settings_duke.xlsx --threads 8
```

`./duke --help` lists every option by module with its default. Add `--dry_run` to any command to print the resolved parameters and tracts without running anything.

From R or RStudio, edit the paths in `scripts/duke_run_local.R` and source it. Every parameter is listed there with its default and a one-line explanation.

---

## Input files

### Sequencing data (`dir_data`)

BAM, FASTQ (`.fastq`, `.fq`, gzipped or not) or FASTA, searched recursively. The sample name (`file_stem`) is the file name without its extension.

### Reference (`path_ref`)

A FASTA with one block of `N` (five or more) where the repeat is. Duke aligns each read to the sequence either side of the block, so the block can be any length. `www/HTTset20.fasta` covers the HTT CAG tract; `www/HTTset20_cag_ccg.fasta` covers the CAG and the CCG tract together.

### Settings file (`path_settings`, needed for Module 6)

One row per sample, `.xlsx` (first sheet) or `.csv`:

| Column | Required | Example | Meaning |
|---|---|---|---|
| `file_name` | yes | `sample_01.bam` | The input file (with or without extension) |
| `analysis_ranges` | yes | `[0-35][36-NA]` | Ranges on the primary measure, one bracket each; `NA` = unbounded |
| `floor` | no | `[3][3]` | Remove lengths with fewer reads than this (a value below 1 is a fraction of the modal bin); 1 keeps all |
| `max_peaks` | no | `[2][1]` | Modal peaks reported per range |
| `group` | no | `QS5` | Group for controls, timecourses and plots |
| `time` | no | `0`, `49` | Time point |
| `group_control_sample` | no | `TRUE` | The sample is a control for its group |
| `exclude` | no | `TRUE` | Leave the sample out of Module 6 |
| `manual_control_repeat_length` | no | `[14][120]` | Control setpoint per range; used only with `use_manual_controls = TRUE` |
| `analysis_ranges_<measure>` (and `floor_<measure>`, `max_peaks_<measure>`, `manual_control_repeat_length_<measure>`) | no | `[NA-8][9-NA]` | Ranges for another tract or composite on its own scale (see below) |

Group and range names may contain any characters; those that are unsafe in file names are replaced with `-` in output paths only. Examples: `www/settings_example.xlsx`, `www/settings_example_two_tracts.xlsx`.

### Other optional files

- `path_tracts` – the tracts table (next section).
- `path_trim_patterns` – adapter CSV (`adapter_name`, `adapter_sequence`) for `trim = TRUE`; `www/adapters.csv`.
- `path_manual_exclusions` – reads to drop (`file_name`, `read_name`, `reason`; Excel sheet `Exclusions`); `www/read_exclusions_example.xlsx`.

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

- `sum` – the tracts' counts added (pure repeat content only);
- `span` – from the start of the first tract to the end of the last, divided by the motif length, so any sequence between the tracts is included (in bp when the motifs differ in length). For HTT, span = CAG + CCG + 4 for the standard allele.

**What you get for every tract and composite (each a "measure"):**

- Module 3: detection rate, counts by every method, modal by method, and the per-read export (`CAG_full`, `CAG_match`, `CAG_pure`, `CCG_full`, ..., `CAG_CCG_sum`, `CAG_CCG_span`).
- Module 6: distributions and the full set of instability metrics, per measure. A measure without its own `analysis_ranges_<measure>` column uses the primary tract's ranges: the reads are chosen by their CAG length (e.g. the expanded allele) and the metrics are calculated on that measure. With its own column the ranges are on its own scale.
- Module 6, tract against tract: for each pair (`measure_pairs`, default every pair of tracts), the mean change in tract B against the change in tract A from the sample median, per group, plus a per-read scatter and a Spearman correlation per sample. This shows whether CAG instability carries CCG with it.
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

**Steps.** Each flank is found by exact match, then with up to `umi_flank_max_mismatch` edits (indels allowed); the UMI is the sequence between them. UMIs are clustered with the umi_tools directional method (a UMI joins a neighbour within `umi_cluster_max_dist` edits that has at least 2n – 1 of its reads). Each UMI family is then split by allele: the allele boundary is found per sample from the valley in the repeat-length distribution, and a family with at least two reads on each allele is split into two molecules (`umi_split_mixed`). Molecules with fewer than `umi_min_reads` reads are dropped. Each molecule's length is `top_cluster` by default: the densest window of reads (`umi_top_cluster_k` reads within ±(`w` + `wfrac` × length)), which ignores PCR stutter and chimeric outliers better than the median.

**Report.** Extraction outcome per sample (found, no 5' flank, no 3' flank, wrong length), UMI length, reads per molecule, UMI neighbours against the number expected by chance (a check that clustering merges errors, not distinct UMIs), within-molecule length error by allele, and reads against molecules.

**Downstream.** Modules 6 and 7 report every distribution and metric for reads and for molecules (`unit` column). Each measure is collapsed separately, so a molecule's `CAG_CCG_sum` is the collapsed sum, not the sum of the collapsed CAG and CCG. `04_umi/molecules/{sample}.tsv.gz` lists every molecule.

---

## How repeats are counted

Duke counts **pure repeat units**. For a 47Q HTT allele (`(CAG)n CAA CAG CCG CCA (CCG)7`): standard 45 CAG, DOI (loss of CAA CAG duplicated) 43 CAG, LOI (loss of interruption, CAA → CAG) 47 CAG.

A tract is a run of motif matches, with gaps of up to `rpt_max_gap` bases inside a run and `rpt_max_tract_gap` between runs merged. It must start and end with `rpt_start_perfect_repeats` / `rpt_end_perfect_repeats` perfect units, which trims e.g. a trailing `CAACAG` off the CAG tract. The search covers the sequence between the two flank alignments plus `rpt_flank_pad` bases (30) either side, because the aligner sometimes places the last unit inside a flank (the cause of the 1-unit under-count of LOI alleles in v2).

Three counts per tract:

| Method | Definition |
|---|---|
| `full` (default) | tract length in bp ÷ motif length |
| `match` | exact motif matches in the tract |
| `pure` | the longest uninterrupted run |

Module 3 reports the modal length of every sample by every method (`modal_by_method`), and how often the methods agree.

---

## Instability metrics

For each sample × measure × range × unit, with setpoint s and d = length – s per read:

| Metric | Definition |
|---|---|
| `instability_index` (II) | mean(d) = EI + CI (Lee et al. 2010) |
| `expansion_index` (EI) | sum of positive d ÷ n |
| `contraction_index` (CI) | sum of negative d ÷ n |
| `mean_expansion_size` | mean d of the reads above s (was "expansion index" in v2) |
| `mean_contraction_size` | mean d of the reads below s (was "contraction index" in v2) |
| `prop_at`, `prop_above`, `prop_below` | proportion within ±0.5 of s, above, below |
| `prop_above_k`, `prop_below_k` | proportion at least k units above / below s (`expansion_thresholds`, default 5 and 10) |
| `expansion_ratio`, `contraction_ratio` | reads above / at s, below / at s (unstable when few reads sit at s) |
| `modal_diff`, `median_diff`, `mean_diff`, `z_score` | difference from the group control |

Two setpoints are reported (`setpoint_type`):

- `sample_relative` – the sample's own `control_setpoint_metric` (median by default; `modal_length` gives Lee's II);
- `control_relative` – that metric aggregated (`control_aggregation_method`) across the group's control samples, so the metrics measure change since the control (e.g. t = 0). With `use_manual_controls = TRUE` a `manual_control_repeat_length` replaces it; `setpoint_source` records which was used.

Distribution statistics (modal, median, mean, SD, IQR, 5th and 95th percentiles, skewness, kurtosis) are in `distribution_summary`.

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

- `06_instability/instability.xlsx`: `instability_metrics` (one row per sample × unit × measure × range × setpoint type; the column names used by v2 scripts are kept), `distribution_summary`, `modal_peaks`, `group_controls`, `count_methods`, `pair_stats`, `range_table`.
- `03_repeat_detection/read_repeat_counts/{sample}.tsv.gz`: one row per read with every tract's counts and structure, the composites, and the v2 columns `repeat_count_full`, `repeat_count_match`, `repeat_count_tracts`.
- `02_alignment/alignment.xlsx`, sheet `read_accounting`: reads per sample at every step from input to analysed.

Logs go to `logs/<date_time>/` in the Duke folder: the run log, the resolved parameters as R code and the resolved tracts table.

---

## Parameters

Every parameter, its default and a one-line explanation is in `scripts/duke_run.R` and `./duke --help`; both are generated from `lib/params.R`. The ones most often changed:

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

Monitoring: `qstat -u $USER`, `qdel <id>`, `qacct -j <id>` (SGE); `squeue -u $USER`, `scancel <id>`, `sacct -j <id>` (Slurm); `tail -f logs/duke_*.out`.

**Run time.** Module 3 is vectorised (seconds per sample) and Module 7 writes one card per sample instead of one plot per sample × range × plot type, so Modules 3 to 7 take minutes rather than hours. Module 2 (minimap2) and, when enabled, trimming and clustering dominate. On 2 cores, 14 PacBio samples of 2,500 reads run Modules 1 to 7 with clustering in about 3 minutes.

---

## Common workflows

**First look at new data** – run without a settings file; Modules 1–5 and 7 run, and the Module 3 and 7 plots show where to put the ranges.

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
- **Module 6 skipped** – no settings file, or no file names in it match the data (`missing` samples are listed in the report).
- **Results did not change after updating Duke** – resume reused old modules; delete `module_data/`.
- **"umi_min_reads must exceed umi_top_cluster_k"** – with exactly k reads `top_cluster` returns their median; raise `umi_min_reads` or change the statistic.

---

## Changes from v2.3.1

**New**
- Several tracts in one run (any number), on one shared reference or separate references, with the sequence between adjacent tracts reported per read, `sum` and `span` composite measures, per-tract metrics, and tract-against-tract plots and correlations.
- Optional UMI module (Module 4): flank-based extraction with user-set flanks (default M13-F / ciosiF+1), directional clustering, allele-aware splitting of mixed families, `top_cluster` collapse; every Module 6 and 7 result for reads and for molecules.
- Modal, median and mean length by every counting method (Module 3 and `count_methods` in Module 6).
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

**Removed**
- Per-sample × per-range histogram, scatter and density PNGs (replaced by the sample cards and heatmaps). Parameters for them are ignored with a message.

---

## Citation and contact

Hölbling, B.V. et al. A multimodal screening platform for endogenous dipeptide repeat proteins in C9orf72 patient iPSC-neurons. *Cell Reports* (2025). https://www.sciencedirect.com/science/article/pii/S2211124725004668

Michael Flower, UCL Queen Square Institute of Neurology – michael.flower@ucl.ac.uk – https://github.com/mike-flower/duke
