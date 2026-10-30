#!/usr/bin/env Rscript
# ==============================================================================
# Duke - run script
# ==============================================================================
# Cluster version (Myriad / Kathleen paths). duke_run_local.R is the laptop version.
# Edit the parameters below, then run from RStudio (Source) or with
#   Rscript scripts/<this file>
# Every parameter is listed with its default; see README.md for detail. The
# same parameters are available on the command line: ./duke --help
# ==============================================================================

# Find the Duke root (the folder holding lib/ and modules/)
duke_root <- local({
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
  } else {
    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    if (length(f)) dirname(dirname(normalizePath(f))) else getwd()
  }
})
setwd(duke_root)
source("lib/load_all.R")

params <- list(
  # --- Paths -------------------------------------------------------------------
  # Input directory of BAM / FASTQ / FASTA files (required)
  dir_data = "/home/skgtmdf/Scratch/data/2025.12.17_pb_test/data",
  # Output directory (required)
  dir_out = "/home/skgtmdf/Scratch/data/2025.12.17_pb_test/result_duke",
  # Reference FASTA with NNNNN marking the repeat region (required unless every
  # tract in path_tracts names its own reference)
  path_ref = "/home/skgtmdf/Scratch/refs/HTTset20/HTTset20.fasta",
  # Settings .xlsx/.csv (ranges, groups, controls, time). Needed for Module 6
  path_settings = "/home/skgtmdf/Scratch/data/2025.12.17_pb_test/settings/settings_duke.xlsx",
  # Optional tracts table (.csv or .xlsx sheet 'tracts'; optional sheet
  # 'composites'). NA = one tract from rpt_pattern and path_ref
  path_tracts = NA_character_,
  # Adapter CSV (adapter_name, adapter_sequence), used when trim = TRUE
  path_trim_patterns = "/home/skgtmdf/Scratch/refs/adapters/adapters.csv",
  # Optional read exclusions file (file_name, read_name, reason)
  path_manual_exclusions = NA_character_,

  # --- Import ------------------------------------------------------------------
  # Regex patterns for files to import
  import_patterns = c("\\.fastq$", "\\.fastq\\.gz$", "\\.fq$", "\\.fq\\.gz$", "\\.fasta$", "\\.fa$", "\\.bam$"),
  # Search dir_data recursively
  import_recursive = TRUE,
  # Pattern identifying R1 files in paired-end data
  r1_pattern = "_R1_",
  # Pattern identifying R2 files in paired-end data
  r2_pattern = "_R2_",
  # Paired-end: import 'R1', 'R2' or 'all'
  select_one_of_pair = "R1",
  # Maximum reads per sample (NA = all)
  downsample = NA_integer_,
  # Drop PacBio HiFi reads with fewer passes than this (NA = no filter)
  np_min = NA_integer_,
  # Maximum points in the PacBio passes-vs-length plot
  pacbio_scatter_downsample = 50000L,
  # Remove duplicate read names, keeping the longest. Off by default: PacBio
  # and ONT read names are unique, and the check removed 0 reads in 36 M
  check_duplicate_readnames = FALSE,
  # Trim adapters from read ends. Off by default: demultiplexed HiFi reads are
  # already trimmed (20 of 36 M reads were trimmed in the msn run)
  trim = FALSE,
  # Maximum errors in an adapter match
  trim_max_mismatch = 3L,
  # Allow indels in adapter matching
  trim_with_indels = TRUE,

  # --- Alignment ---------------------------------------------------------------
  # minimap2 arguments for aligning reads to each flank
  minimap2_args = "-t 1 -x sr -w 1",
  # Keep the per-sample flank BAMs in 02_alignment/bam/ for IGV
  keep_bam = FALSE,
  # Write per-sample flank alignment plots
  visualise_alignment = TRUE,
  # Maximum reads per alignment plot
  visualise_alignment_downsample = 1000L,

  # --- Repeats -----------------------------------------------------------------
  # Repeat motif (used when there is no tracts table)
  rpt_pattern = "CAG",
  # Name for the single tract (default: the motif)
  rpt_tract_name = NA_character_,
  # Minimum tract length in repeat units
  rpt_min_repeats = 2L,
  # Mismatches allowed per repeat unit
  rpt_max_mismatch = 0L,
  # Perfect units required at the start of a tract
  rpt_start_perfect_repeats = 2L,
  # Perfect units required at the end of a tract
  rpt_end_perfect_repeats = 2L,
  # Maximum bases between units inside a tract
  rpt_max_gap = 6L,
  # Maximum gap between sub-tracts before they are merged
  rpt_max_tract_gap = 18L,
  # Which tract to keep when several are found: 'longest' or 'first'
  rpt_return_option = "longest",
  # Bases either side of the flank boundaries searched for the tract, so a
  # tract can extend past where the aligner placed the flank
  rpt_flank_pad = 30L,
  # Primary count: 'full' (tract length / motif length), 'match' (exact units)
  # or 'pure' (longest uninterrupted run)
  repeat_count_method = "full",
  # Reads with no tract: 'convert_to_zero', 'filter' or 'keep'
  na_repeat_handling = "convert_to_zero",
  # Combined measures, e.g. 'CAG_CCG_sum=CAG+CCG:sum' (pure units added) or
  # 'CAG_CCG_span=CAG+CCG:span' (first tract start to last tract end)
  composites = character(0),
  # Write per-read lengths to 03_repeat_detection/read_repeat_counts/
  export_read_counts = TRUE,
  # Remove reads with outlying combined flank length
  rm_flank_length_outliers = FALSE,
  # IQR multiplier for the flank length filter
  flank_iqr_multiplier = 1.5,

  # --- UMI ---------------------------------------------------------------------
  # Extract UMIs and collapse reads to molecules
  umi = FALSE,
  # Sequence immediately 5' of the UMI (default M13-F)
  umi_flank_5 = "GTAAAACGACGGCCAGT",
  # Sequence immediately 3' of the UMI (default ciosiF+1)
  umi_flank_3 = "GGCGACCCTGGAAAAGCTGATGAA",
  # Edits allowed when locating each UMI flank (indel tolerant)
  umi_flank_max_mismatch = 3L,
  # The 5' flank must start within this many bases of the read start
  umi_search_window = 250L,
  # Shortest UMI accepted
  umi_len_min = 12L,
  # Longest UMI accepted
  umi_len_max = 16L,
  # Edit distance for directional UMI clustering (0 = exact match; 1 or 2)
  umi_cluster_max_dist = 1L,
  # Minimum reads for a molecule to be kept
  umi_min_reads = 3L,
  # Per-molecule length: 'top_cluster', 'median', 'mode' or 'mean'
  umi_collapse_stat = "top_cluster",
  # top_cluster: reads required in the window
  umi_top_cluster_k = 2L,
  # top_cluster: constant part of the window half-width (repeats)
  umi_top_cluster_w = 3,
  # top_cluster: length-proportional part of the window half-width
  umi_top_cluster_wfrac = 0.1,
  # top_cluster: share of the molecule's reads the window must also hold (0 = only
  # umi_top_cluster_k reads needed)
  umi_top_cluster_min_share = 0.2,

  # --- Allele structure --------------------------------------------------------
  # Cluster reads into alleles, build consensus sequences and call variants
  cluster = FALSE,
  # 'repeat', 'haplotype', or 'repeat,haplotype'
  cluster_by = "repeat",
  # Maximum repeat-length clusters per sample
  repeat_cluster_max = 6L,
  # Clusters holding fewer than this proportion of reads are merged into their
  # nearest neighbour
  repeat_cluster_min_prop = 0.05,
  # Reads used to choose the number of repeat clusters
  repeat_cluster_downsample = 10000L,
  # Maximum haplotype clusters per sample
  haplotype_cluster_max = 6L,
  # Reads used to define haplotype clusters
  haplotype_cluster_downsample = 3000L,
  # Flank(s) used for haplotypes: 'left', 'right' or 'both'
  haplotype_region = "both",
  # 'auto', a number of bases, or 'NA' for no trimming
  haplotype_trim_length = "auto",
  # Minimum agreement for a consensus base
  consensus_threshold = 0.5,
  # Maximum reads per consensus alignment
  consensus_downsample = 50L,
  # Call variants in each cluster consensus against the reference
  call_variants = TRUE,
  # Write per-sample waterfall plots
  waterfall = TRUE,
  # Maximum reads per waterfall
  waterfall_downsample = 1000L,
  # Bases of flank drawn either side of the tract in waterfalls
  waterfall_flank_bp = 60L,
  # Also write one waterfall per cluster
  waterfall_per_cluster = FALSE,

  # --- Instability -------------------------------------------------------------
  # Window for modal peak detection
  range_peak_span = 3L,
  # Calculate group control setpoints
  group_control = TRUE,
  # 'flagged', 'earliest' or 'all'
  control_sample_selection = "flagged",
  # 'modal_length', 'mean_length' or 'median_length'
  control_setpoint_metric = "median_length",
  # 'mean', 'median' or 'trimmed_mean' across control samples
  control_aggregation_method = "median",
  # Use manual_control_repeat_length from the settings file as the control
  # setpoint (replacing the group control) for those samples and ranges
  use_manual_controls = FALSE,
  # k for proportion of reads at least k units above / below the setpoint
  expansion_thresholds = c(5L, 10L),
  # Dimensions for cohort plots: group, time
  cohort_plot_by = c("group", "time"),
  # Metrics drawn as timecourses by group
  headline_metrics = c("modal_length", "median_length", "mean_length", "instability_index", "expansion_index", "prop_above_5"),
  # Setpoint type for headline plots: 'control_relative' or 'sample_relative'
  headline_setpoint_type = "control_relative",
  # Ranges drawn inline (default: the first two ranges of the primary measure)
  plot_ranges = character(0),
  # Pairs of measures to plot against each other per read, e.g. 'CAG:CCG'
  # (default: every pair of tracts)
  measure_pairs = character(0),
  # Maximum reads per sample in per-read pair plots
  pair_max_points = 4000L,

  # --- Repeat lengths ----------------------------------------------------------
  # Write one summary figure per sample
  sample_cards = TRUE,
  # Histogram bin width in repeat units
  repeat_histogram_binwidth = 1L,

  # --- Runtime -----------------------------------------------------------------
  # Resolution of saved plots
  plot_dpi = 150L,
  # Write per-sample plots in Modules 2 and 3
  plot_per_sample = TRUE,
  # CPU cores
  threads = 80L,
  # Reuse a module's results when the parameters that affect it are unchanged
  resume = TRUE,
  # Free memory between steps
  remove_intermediate = TRUE,
  # Delete temp/ (per-sample caches) at the end
  remove_temp = FALSE,
  # Modules to run
  run_modules = 1:7
)

run_duke(params, duke_root)
