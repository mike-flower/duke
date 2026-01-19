#!/usr/bin/env Rscript
#
# Duke Pipeline - Command-Line Interface
# 
# This script provides command-line access to ALL Duke Pipeline parameters
# with the EXACT same structure, names, and comments as duke_run.R
#
# Usage:
#   Rscript duke_cli.R --dir_data /path/to/data --dir_out /path/to/output [options]
#
# Examples:
#   Rscript duke_cli.R --dir_data ~/data --dir_out ~/results
#   Rscript duke_cli.R --dir_data ~/data --dir_out ~/results --threads 80 --resume TRUE
#   Rscript duke_cli.R --help
#

# ==============================================================================
# Set Working Directory to Duke Installation Root
# ==============================================================================

# This script lives in scripts/ subdirectory, so we need to go up one level
# to reach the duke installation root where lib/ and modules/ are located

args_temp <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_temp, value = TRUE)

if (length(file_arg) > 0) {
  # Extract script path and get duke root (parent of scripts/)
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(normalizePath(script_path))
  duke_root <- dirname(script_dir)  # Go up one level from scripts/
  setwd(duke_root)
} else {
  # Fallback: assume we're already in duke root
  cat("Warning: Could not detect script location. Assuming current directory is duke root.\n")
}

cat("Duke installation directory:", getwd(), "\n\n")

# ==============================================================================
# Parse Command-Line Arguments
# ==============================================================================

library(optparse)

option_list <- list(
  
  # ===========================================================================
  # File paths and directories
  # ===========================================================================
  
  make_option("--dir_data", type="character", default=NULL,
              help="Input data directory (FASTQ files) [REQUIRED]", metavar="PATH"),
  
  make_option("--dir_out", type="character", default=NULL,
              help="Output results directory [REQUIRED]", metavar="PATH"),
  
  make_option("--path_ref", type="character", default=NULL,
              help="Reference FASTA file with NNNNN separator [REQUIRED]", 
              metavar="PATH"),
  
  # Optional paths
  make_option("--path_trim_patterns", type="character", default=NULL,
              help="Path to adapter trimming patterns file", metavar="PATH"),
  
  make_option("--path_manual_exclusions", type="character", default=NULL,
              help="Path to file with reads to exclude", metavar="PATH"),
  
  # ===========================================================================
  # File import options
  # ===========================================================================
  
  make_option("--import_patterns", type="character",
              default="\\.fastq$,\\.fastq.gz$,\\.fasta$,\\.fa$,\\.bam$",
              help="File patterns to import (comma-separated) [default: %default]",
              metavar="PATTERNS"),
  
  make_option("--import_recursive", type="logical", default=TRUE,
              help="Recursively search directories [default: %default]"),
  
  # Paired-end read handling
  make_option("--r1_pattern", type="character", default="_R1_",
              help="R1 read identifier pattern [default: %default]", metavar="PATTERN"),
  
  make_option("--r2_pattern", type="character", default="_R2_",
              help="R2 read identifier pattern [default: %default]", metavar="PATTERN"),
  
  make_option("--select_one_of_pair", type="character", default="R1",
              help="Keep only R1, R2, or 'all' to keep all [default: %default]", 
              metavar="OPTION"),
  
  # Downsampling (set to NA to disable)
  make_option("--downsample", type="integer", default=NULL,
              help="Downsample to N reads per sample (e.g., 1000) [default: NA = disabled]",
              metavar="N"),
  
  # ===========================================================================
  # Adapter trimming
  # ===========================================================================
  
  make_option("--trim", type="logical", default=TRUE,
              help="Enable adapter trimming (requires --path_trim_patterns if TRUE) [default: %default]"),
  
  make_option("--trim_max_mismatch", type="integer", default=3,
              help="Max mismatches in adapter matching (0-5) [default: %default]", metavar="N"),
  
  make_option("--trim_with_indels", type="logical", default=TRUE,
              help="Allow indels in adapter matching [default: %default]"),
  
  # ===========================================================================
  # Alignment with minimap2
  # ===========================================================================
  
  make_option("--minimap2_args", type="character", default="-t 2 -x sr -w 1",
              help="minimap2 arguments [default: %default]", metavar="ARGS"),
  
  make_option("--visualise_alignment", type="logical", default=TRUE,
              help="Generate alignment coverage plots [default: %default]"),
  
  make_option("--visualise_alignment_downsample", type="integer", default=1000,
              help="Max reads to plot per sample (NA = plot all, may be slow) [default: %default]"),
  
  # ===========================================================================
  # Repeat detection and counting
  # ===========================================================================
  # These parameters control how CAG/CTG repeats are identified in sequences
  #
  # OVERVIEW:
  # Duke finds repeat tracts by scanning for imperfect runs of a motif (e.g., CAG)
  # It tolerates small interruptions (mismatches, gaps) but requires perfect repeats
  # at tract boundaries to confidently identify start/end positions.
  #
  # KEY CONCEPTS:
  # - "Perfect repeat" = exact match to motif (e.g., "CAG")
  # - "Mismatch" = substitution in motif (e.g., "CAT" instead of "CAG")
  # - "Gap" = inserted bases between repeats (e.g., "CAG-TA-CAG" has 2bp gap)
  # - "Tract" = continuous run of repeats (may contain small gaps/mismatches)
  #
  # PARAMETER GUIDE:
  
  make_option("--rpt_pattern", type="character", default="CAG",
              help="Repeat motif to search for (e.g., 'CAG' for HD, 'CTG' for DM1) [default: %default]", 
              metavar="MOTIF"),
  
  make_option("--rpt_min_repeats", type="integer", default=2,
              help="Minimum repeat count to consider valid tract. Lower = more sensitive but more false positives [default: %default]", 
              metavar="N"),
  
  make_option("--rpt_max_mismatch", type="integer", default=0,
              help="Max substitution errors tolerated in each repeat unit. 0=strict perfect matching, 1+=tolerant to sequencing errors [default: %default]", 
              metavar="N"),
  
  make_option("--rpt_start_perfect_repeats", type="integer", default=2,
              help="Number of perfect (error-free) repeats required at tract START. Higher = more confident boundaries [default: %default]", 
              metavar="N"),
  
  make_option("--rpt_end_perfect_repeats", type="integer", default=2,
              help="Number of perfect (error-free) repeats required at tract END. Higher = more confident boundaries [default: %default]", 
              metavar="N"),
  
  make_option("--rpt_max_gap", type="integer", default=6,
              help="Max inserted bases allowed WITHIN a single repeat tract (e.g., 'CAG-TA-CAG' = 2bp gap). Higher = more tolerant [default: %default]", 
              metavar="N"),
  
  make_option("--rpt_max_tract_gap", type="integer", default=18,
              help="Max gap between separate repeat tracts to merge them (e.g., tract1--18bp--tract2). Higher = merges distant tracts [default: %default]", 
              metavar="N"),
  
  make_option("--rpt_return_option", type="character", default="longest",
              help="Which tract to return if multiple found: 'longest' (use biggest) or 'all' (keep all) [default: %default]", 
              metavar="OPTION"),
  
  make_option("--repeat_count_method", type="character", default="repeat_count_full",
              help="Repeat counting method: repeat_count_full/repeat_count_match [default: %default]",
              metavar="METHOD"),
  
  make_option("--na_repeat_handling", type="character", default="convert_to_zero",
              help="Handle NA repeats: convert_to_zero/filter/flag_only [default: %default]",
              metavar="METHOD"),
  
  # ===========================================================================
  # Clustering
  # ===========================================================================
  
  make_option("--cluster", type="logical", default=TRUE,
              help="Perform clustering [default: %default]"),
  
  make_option("--cluster_by", type="character", default="repeat",
              help="Clustering strategy: repeat/haplotype/repeat,haplotype/none [default: %default]",
              metavar="STRATEGY"),
  
  make_option("--haplotype_cluster_max", type="integer", default=10,
              help="Max haplotypes to consider [default: %default]", metavar="N"),
  
  make_option("--repeat_cluster_max", type="integer", default=20,
              help="Max repeats to consider [default: %default]", metavar="N"),
  
  make_option("--repeat_cluster_downsample", type="integer", default=10000,
              help="Max reads for repeat clustering optimisation (speeds up >10k reads) [default: %default]", 
              metavar="N"),
  
  make_option("--haplotype_cluster_downsample", type="integer", default=5000,
              help="Max reads for haplotype clustering optimisation (speeds up >5k reads) [default: %default]", 
              metavar="N"),
  
  make_option("--haplotype_region", type="character", default="both",
              help="Haplotype region: left/right/both [default: %default]", metavar="REGION"),
  
  make_option("--haplotype_method", type="character", default="levenshtein",
              help="Haplotype distance method: hamming/levenshtein [default: %default]",
              metavar="METHOD"),
  
  make_option("--haplotype_trim_length", type="character", default="auto",
              help="Haplotype trim length: auto/NA/numeric [default: %default]",
              metavar="LENGTH"),
  
  # ===========================================================================
  # Consensus generation and variant calling
  # ===========================================================================
  
  make_option("--cluster_consensus", type="logical", default=TRUE,
              help="Generate consensus sequences [default: %default]"),
  
  make_option("--consensus_threshold", type="numeric", default=0.5,
              help="Minimum proportion for consensus base [default: %default]", metavar="FRAC"),
  
  make_option("--consensus_downsample", type="integer", default=50,
              help="Max reads per consensus [default: %default]", metavar="N"),
  
  make_option("--call_variants", type="logical", default=TRUE,
              help="Call variants in consensus sequences [default: %default]"),
  
  # ===========================================================================
  # Waterfall plots (Module 5)
  # ===========================================================================
  
  make_option("--waterfall", type="logical", default=TRUE,
              help="Generate waterfall plots [default: %default]"),
  
  make_option("--waterfall_downsample", type="integer", default=1000,
              help="Max reads per waterfall plot [default: %default]", metavar="N"),
  
  make_option("--waterfall_rm_flank_length_outliers", type="logical", default=TRUE,
              help="Remove flank length outliers from waterfall [default: %default]"),
  
  make_option("--waterfall_y_axis_labels", type="character", default="auto",
              help="Waterfall y-axis labels: auto/all/integer [default: %default]",
              metavar="OPTION"),
  
  make_option("--waterfall_per_cluster", type="logical", default=TRUE,
              help="Generate per-cluster waterfall plots [default: %default]"),
  
  make_option("--dna_colours", type="character",
              default="A=darkgreen,C=blue,T=red,G=black,-=lightgrey,N=grey",
              help="DNA base colors (base=colour,...) [default: %default]", metavar="COLORS"),
  
  # ===========================================================================
  # Range analysis & instability metrics (Module 6)
  # ===========================================================================
  
  make_option("--path_settings", type="character", default=NULL,
              help="Settings Excel file for Module 6", metavar="PATH"),
  
  make_option("--range_peak_span", type="integer", default=3,
              help="Peak detection smoothing window [default: %default]", metavar="N"),
  
  make_option("--group_control", type="logical", default=TRUE,
              help="Enable group control analysis [default: %default]"),
  
  make_option("--control_sample_selection", type="character", default="flagged",
              help="Control sample selection: flagged/earliest/all [default: %default]",
              metavar="METHOD"),
  
  make_option("--control_setpoint_metric", type="character", default="median_length",
              help="Control setpoint metric: modal_length/mean_length/median_length [default: %default]",
              metavar="METRIC"),
  
  make_option("--control_aggregation_method", type="character", default="median",
              help="Control aggregation method: mean/median/trimmed_mean [default: %default]",
              metavar="METHOD"),
  
  # ===========================================================================
  # Module 7: Repeat Distribution Visualisation
  # ===========================================================================
  
  make_option("--repeat_histogram", type="logical", default=TRUE,
              help="Generate repeat frequency histograms [default: %default]"),
  
  make_option("--repeat_histogram_binwidth", type="integer", default=1,
              help="Histogram bin width [default: %default]", metavar="N"),
  
  make_option("--repeat_scatter", type="logical", default=TRUE,
              help="Generate scatter/violin plots [default: %default]"),
  
  make_option("--repeat_distribution_metrics", type="character",
              default="modal_length,mean_length,median_length",
              help="Distribution metrics (comma-separated) [default: %default]",
              metavar="METRICS"),
  
  # ===========================================================================
  # Runtime settings
  # ===========================================================================
  
  make_option("--threads", type="integer", default=12,
              help="Number of CPU threads [default: %default]", metavar="N"),
  
  make_option("--resume", type="logical", default=TRUE,
              help="Resume from previous run [default: %default]"),
  
  # ===========================================================================
  # Memory and disk management
  # ===========================================================================
  
  make_option("--remove_intermediate", type="logical", default=TRUE,
              help="Free up RAM during pipeline execution [default: %default]"),
  
  make_option("--cleanup_temp", type="logical", default=TRUE,
              help="Free up disk space after completion [default: %default]"),
  
  # ===========================================================================
  # CLI-specific options (not in duke_run.R)
  # ===========================================================================
  
  make_option("--run_modules", type="character", default="1,2,3,4,5,6,7",
              help="Modules to run (comma-separated, e.g., '1,2,3' or '5,6,7') [default: %default]",
              metavar="MODULES"),
  
  make_option("--log_dir", type="character", default="logs",
              help="Log directory [default: %default]", metavar="PATH"),
  
  make_option("--dry_run", type="logical", default=FALSE,
              help="Show configuration without running [default: %default]"),
  
  make_option("--verbose", type="logical", default=TRUE,
              help="Verbose output [default: %default]")
)

parser <- OptionParser(
  usage = "%prog --dir_data DATA --dir_out OUTPUT [options]",
  option_list = option_list,
  description = "\nDuke Pipeline - Command-Line Interface\n\nAll parameters match duke_run.R exactly.",
  epilogue = "For documentation, see README.md"
)

args <- parse_args(parser)

# ==============================================================================
# Validate Required Arguments
# ==============================================================================

if (is.null(args$dir_data)) {
  stop("\nError: --dir_data is required\n",
       "Usage: Rscript duke_cli.R --dir_data /path/to/data --dir_out /path/to/output\n")
}

if (is.null(args$dir_out)) {
  stop("\nError: --dir_out is required\n",
       "Usage: Rscript duke_cli.R --dir_data /path/to/data --dir_out /path/to/output\n")
}

# ==============================================================================
# Validate Required Arguments
# ==============================================================================

if (is.null(args$dir_data)) {
  stop("Error: --dir_data is required\n")
}

if (is.null(args$dir_out)) {
  stop("Error: --dir_out is required\n")
}

if (is.null(args$path_ref)) {
  stop("Error: --path_ref is required\n  Specify reference FASTA file with: --path_ref /path/to/reference.fasta\n")
}

# ==============================================================================
# Parse Complex Arguments
# ==============================================================================

import_patterns <- strsplit(args$import_patterns, ",")[[1]]
cluster_by_parsed <- strsplit(args$cluster_by, ",")[[1]]
if (length(cluster_by_parsed) == 1 && cluster_by_parsed[1] == "none") {
  cluster_by_parsed <- "none"
}
repeat_distribution_metrics <- strsplit(args$repeat_distribution_metrics, ",")[[1]]

# Parse run_modules
run_modules <- as.integer(strsplit(args$run_modules, ",")[[1]])
if (any(run_modules < 1 | run_modules > 7)) {
  stop("Error: --run_modules must contain values between 1-7\n")
}
if (args$verbose) {
  cat("Modules to run:", paste(run_modules, collapse=", "), "\n")
}

dna_colours_parsed <- c()
if (!is.null(args$dna_colours)) {
  pairs <- strsplit(args$dna_colours, ",")[[1]]
  for (pair in pairs) {
    parts <- strsplit(pair, "=")[[1]]
    if (length(parts) == 2) {
      dna_colours_parsed[parts[1]] <- parts[2]
    }
  }
}

if (!is.null(args$select_one_of_pair) && args$select_one_of_pair == "all") {
  args$select_one_of_pair <- NA
}

if (!is.null(args$haplotype_trim_length) && args$haplotype_trim_length == "NA") {
  args$haplotype_trim_length <- NA
} else if (!is.null(args$haplotype_trim_length) && args$haplotype_trim_length != "auto") {
  args$haplotype_trim_length <- as.numeric(args$haplotype_trim_length)
}

# Auto-detect settings
if (is.null(args$path_settings)) {
  possible <- c(
    file.path(dirname(args$dir_data), "settings", "settings_duke.xlsx"),
    file.path(args$dir_data, "settings_duke.xlsx")
  )
  for (p in possible) {
    if (file.exists(p)) {
      args$path_settings <- p
      if (args$verbose) cat("Auto-detected settings:", p, "\n")
      break
    }
  }
}

# ==============================================================================
# Build params List (EXACTLY matches duke_run.R)
# ==============================================================================

params <- list(
  dir_data = normalizePath(args$dir_data, mustWork = TRUE),
  dir_out = args$dir_out,
  path_ref = normalizePath(args$path_ref, mustWork = TRUE),
  path_trim_patterns = if (!is.null(args$path_trim_patterns)) 
    normalizePath(args$path_trim_patterns, mustWork = TRUE) else NA,
  path_manual_exclusions = if (!is.null(args$path_manual_exclusions))
    normalizePath(args$path_manual_exclusions, mustWork = TRUE) else NA,
  import_patterns = import_patterns,
  import_recursive = args$import_recursive,
  r1_pattern = args$r1_pattern,
  r2_pattern = args$r2_pattern,
  select_one_of_pair = args$select_one_of_pair,
  downsample = if (!is.null(args$downsample)) args$downsample else NA,
  trim = args$trim,
  trim_max_mismatch = args$trim_max_mismatch,
  trim_with_indels = args$trim_with_indels,
  minimap2_args = args$minimap2_args,
  visualise_alignment = args$visualise_alignment,
  visualise_alignment_downsample = args$visualise_alignment_downsample,
  rpt_pattern = args$rpt_pattern,
  rpt_min_repeats = args$rpt_min_repeats,
  rpt_max_mismatch = args$rpt_max_mismatch,
  rpt_start_perfect_repeats = args$rpt_start_perfect_repeats,
  rpt_end_perfect_repeats = args$rpt_end_perfect_repeats,
  rpt_max_gap = args$rpt_max_gap,
  rpt_max_tract_gap = args$rpt_max_tract_gap,
  rpt_return_option = args$rpt_return_option,
  repeat_count_method = args$repeat_count_method,
  na_repeat_handling = args$na_repeat_handling,
  cluster = args$cluster,
  cluster_by = cluster_by_parsed,
  haplotype_cluster_max = args$haplotype_cluster_max,
  repeat_cluster_max = args$repeat_cluster_max,
  repeat_cluster_downsample = args$repeat_cluster_downsample,
  haplotype_cluster_downsample = args$haplotype_cluster_downsample,
  haplotype_region = args$haplotype_region,
  haplotype_method = args$haplotype_method,
  haplotype_trim_length = args$haplotype_trim_length,
  cluster_consensus = args$cluster_consensus,
  consensus_threshold = args$consensus_threshold,
  consensus_downsample = args$consensus_downsample,
  call_variants = args$call_variants,
  waterfall = args$waterfall,
  waterfall_downsample = args$waterfall_downsample,
  waterfall_rm_flank_length_outliers = args$waterfall_rm_flank_length_outliers,
  waterfall_y_axis_labels = args$waterfall_y_axis_labels,
  waterfall_per_cluster = args$waterfall_per_cluster,
  dna_colours = dna_colours_parsed,
  path_settings = if (!is.null(args$path_settings))
    normalizePath(args$path_settings, mustWork = TRUE) else NULL,
  range_peak_span = args$range_peak_span,
  group_control = args$group_control,
  control_sample_selection = args$control_sample_selection,
  control_setpoint_metric = args$control_setpoint_metric,
  control_aggregation_method = args$control_aggregation_method,
  repeat_histogram = args$repeat_histogram,
  repeat_histogram_binwidth = args$repeat_histogram_binwidth,
  repeat_scatter = args$repeat_scatter,
  repeat_distribution_metrics = repeat_distribution_metrics,
  threads = args$threads,
  resume = args$resume,
  remove_intermediate = args$remove_intermediate,
  cleanup_temp = args$cleanup_temp
)

# ==============================================================================
# Setup Logging (EXACTLY matches duke_run.R)
# ==============================================================================

log_dir <- file.path(getwd(), args$log_dir)
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

log_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_log_dir <- file.path(log_dir, log_timestamp)
if (!dir.exists(run_log_dir)) dir.create(run_log_dir, recursive = TRUE)

log_file <- file.path(run_log_dir, paste0("duke_run-", log_timestamp, ".log"))
params_file <- file.path(run_log_dir, paste0("params-", log_timestamp, ".R"))

cat("# Duke Pipeline Run Parameters\n", file = params_file)
cat("# Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", 
    file = params_file, append = TRUE)
cat("params <- ", file = params_file, append = TRUE)
capture.output(dput(params), file = params_file, append = TRUE)

log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
start_time <- Sys.time()

cat("\n=================================================================\n")
cat("                    DUKE PIPELINE                                \n")
cat("=================================================================\n\n")
cat("Run started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Log file:", log_file, "\n")
cat("Parameters saved:", params_file, "\n\n")

# ==============================================================================
# Dry Run - Print All Parameters
# ==============================================================================

if (args$dry_run) {
  cat("\n=================================================================\n")
  cat("                    DRY RUN - CONFIGURATION                      \n")
  cat("=================================================================\n\n")
  
  cat("*** This is a DRY RUN - pipeline will NOT execute ***\n\n")
  
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("FILE PATHS AND DIRECTORIES\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  cat(sprintf("%-35s: %s\n", "dir_data", params$dir_data))
  cat(sprintf("%-35s: %s\n", "dir_out", params$dir_out))
  cat(sprintf("%-35s: %s\n", "path_ref", params$path_ref))
  cat(sprintf("%-35s: %s\n", "path_trim_patterns", 
              if(is.null(params$path_trim_patterns)) "NULL" else params$path_trim_patterns))
  cat(sprintf("%-35s: %s\n", "path_manual_exclusions", 
              if(is.null(params$path_manual_exclusions)) "NULL" else params$path_manual_exclusions))
  cat(sprintf("%-35s: %s\n", "path_settings", 
              if(is.null(params$path_settings)) "NULL" else params$path_settings))
  
  cat("\n═══════════════════════════════════════════════════════════════\n")
  cat("FILE IMPORT OPTIONS\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  cat(sprintf("%-35s: %s\n", "import_patterns", paste(params$import_patterns, collapse = ", ")))
  cat(sprintf("%-35s: %s\n", "import_recursive", params$import_recursive))
  cat(sprintf("%-35s: %s\n", "r1_pattern", params$r1_pattern))
  cat(sprintf("%-35s: %s\n", "r2_pattern", params$r2_pattern))
  cat(sprintf("%-35s: %s\n", "select_one_of_pair", params$select_one_of_pair))
  cat(sprintf("%-35s: %s\n", "downsample", if(is.na(params$downsample)) "NA (disabled)" else params$downsample))
  
  cat("\n═══════════════════════════════════════════════════════════════\n")
  cat("ADAPTER TRIMMING\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  cat(sprintf("%-35s: %s\n", "trim", params$trim))
  cat(sprintf("%-35s: %s\n", "trim_max_mismatch", params$trim_max_mismatch))
  cat(sprintf("%-35s: %s\n", "trim_with_indels", params$trim_with_indels))
  
  cat("\n═══════════════════════════════════════════════════════════════\n")
  cat("ALIGNMENT\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  cat(sprintf("%-35s: %s\n", "minimap2_args", params$minimap2_args))
  cat(sprintf("%-35s: %s\n", "visualise_alignment", params$visualise_alignment))
  cat(sprintf("%-35s: %s\n", "visualise_alignment_downsample", 
              if(is.na(params$visualise_alignment_downsample)) "NA (plot all)" else params$visualise_alignment_downsample))
  
  cat("\n═══════════════════════════════════════════════════════════════\n")
  cat("REPEAT DETECTION\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  cat(sprintf("%-35s: %s\n", "rpt_pattern", params$rpt_pattern))
  cat(sprintf("%-35s: %s\n", "rpt_min_repeats", params$rpt_min_repeats))
  cat(sprintf("%-35s: %s\n", "rpt_max_mismatch", params$rpt_max_mismatch))
  cat(sprintf("%-35s: %s\n", "rpt_start_perfect_repeats", params$rpt_start_perfect_repeats))
  cat(sprintf("%-35s: %s\n", "rpt_end_perfect_repeats", params$rpt_end_perfect_repeats))
  cat(sprintf("%-35s: %s\n", "rpt_max_gap", params$rpt_max_gap))
  cat(sprintf("%-35s: %s\n", "rpt_max_tract_gap", params$rpt_max_tract_gap))
  cat(sprintf("%-35s: %s\n", "rpt_return_option", params$rpt_return_option))
  cat(sprintf("%-35s: %s\n", "repeat_count_method", params$repeat_count_method))
  cat(sprintf("%-35s: %s\n", "na_repeat_handling", params$na_repeat_handling))
  
  cat("\n═══════════════════════════════════════════════════════════════\n")
  cat("CLUSTERING\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  cat(sprintf("%-35s: %s\n", "cluster", params$cluster))
  cat(sprintf("%-35s: %s\n", "cluster_by", paste(params$cluster_by, collapse = ", ")))
  cat(sprintf("%-35s: %s\n", "haplotype_cluster_max", params$haplotype_cluster_max))
  cat(sprintf("%-35s: %s\n", "repeat_cluster_max", params$repeat_cluster_max))
  cat(sprintf("%-35s: %s\n", "repeat_cluster_downsample", params$repeat_cluster_downsample))
  cat(sprintf("%-35s: %s\n", "haplotype_cluster_downsample", params$haplotype_cluster_downsample))
  cat(sprintf("%-35s: %s\n", "haplotype_region", params$haplotype_region))
  cat(sprintf("%-35s: %s\n", "haplotype_method", params$haplotype_method))
  cat(sprintf("%-35s: %s\n", "haplotype_trim_length", params$haplotype_trim_length))
  
  cat("\n═══════════════════════════════════════════════════════════════\n")
  cat("CONSENSUS AND VARIANT CALLING\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  cat(sprintf("%-35s: %s\n", "cluster_consensus", params$cluster_consensus))
  cat(sprintf("%-35s: %s\n", "consensus_threshold", params$consensus_threshold))
  cat(sprintf("%-35s: %s\n", "consensus_downsample", 
              if(is.na(params$consensus_downsample)) "NA (use all)" else params$consensus_downsample))
  cat(sprintf("%-35s: %s\n", "call_variants", params$call_variants))
  
  cat("\n═══════════════════════════════════════════════════════════════\n")
  cat("WATERFALL PLOTS (MODULE 5)\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  cat(sprintf("%-35s: %s\n", "waterfall", params$waterfall))
  cat(sprintf("%-35s: %s\n", "waterfall_downsample", 
              if(is.na(params$waterfall_downsample)) "NA (plot all)" else params$waterfall_downsample))
  cat(sprintf("%-35s: %s\n", "waterfall_rm_flank_length_outliers", params$waterfall_rm_flank_length_outliers))
  cat(sprintf("%-35s: %s\n", "waterfall_y_axis_labels", params$waterfall_y_axis_labels))
  cat(sprintf("%-35s: %s\n", "waterfall_per_cluster", params$waterfall_per_cluster))
  
  cat("\n═══════════════════════════════════════════════════════════════\n")
  cat("RANGE ANALYSIS (MODULE 6)\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  cat(sprintf("%-35s: %s\n", "range_peak_span", params$range_peak_span))
  cat(sprintf("%-35s: %s\n", "group_control", params$group_control))
  cat(sprintf("%-35s: %s\n", "control_sample_selection", params$control_sample_selection))
  cat(sprintf("%-35s: %s\n", "control_setpoint_metric", params$control_setpoint_metric))
  cat(sprintf("%-35s: %s\n", "control_aggregation_method", params$control_aggregation_method))
  
  cat("\n═══════════════════════════════════════════════════════════════\n")
  cat("REPEAT VISUALISATION (MODULE 7)\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  cat(sprintf("%-35s: %s\n", "repeat_histogram", params$repeat_histogram))
  cat(sprintf("%-35s: %s\n", "repeat_histogram_binwidth", params$repeat_histogram_binwidth))
  cat(sprintf("%-35s: %s\n", "repeat_scatter", params$repeat_scatter))
  cat(sprintf("%-35s: %s\n", "repeat_distribution_metrics", 
              paste(params$repeat_distribution_metrics, collapse = ", ")))
  
  cat("\n═══════════════════════════════════════════════════════════════\n")
  cat("RUNTIME AND EXECUTION\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  cat(sprintf("%-35s: %s\n", "threads", params$threads))
  cat(sprintf("%-35s: %s\n", "resume", params$resume))
  cat(sprintf("%-35s: %s\n", "remove_intermediate", params$remove_intermediate))
  cat(sprintf("%-35s: %s\n", "cleanup_temp", params$cleanup_temp))
  cat(sprintf("%-35s: %s\n", "run_modules", paste(run_modules, collapse = ", ")))
  cat(sprintf("%-35s: %s\n", "log_dir", log_dir))
  cat(sprintf("%-35s: %s\n", "verbose", args$verbose))
  
  cat("\n═══════════════════════════════════════════════════════════════\n")
  cat("EXECUTION PLAN\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  module_names <- c("Module 1: Import and QC", 
                    "Module 2: Alignment", 
                    "Module 3: Repeat Detection",
                    "Module 4: Allele Calling", 
                    "Module 5: Waterfall Plots", 
                    "Module 6: Range Analysis",
                    "Module 7: Repeat Visualisation")
  
  cat("Modules to run:\n")
  for (mod in run_modules) {
    status <- ""
    if (mod == 5 && !params$waterfall) {
      status <- " [SKIP - waterfall disabled]"
    } else if (params$resume) {
      output <- file.path(params$dir_out, "module_data", 
                         paste0("0", mod, "_", c("import_qc", "alignment", "repeat_detection",
                                                  "allele_calling", "waterfall", "range_analysis",
                                                  "repeat_visualisation")[mod], "_results.RData"))
      if (file.exists(output)) {
        status <- " [SKIP - results found]"
      } else {
        status <- " [RUN]"
      }
    } else {
      status <- " [RUN]"
    }
    cat(sprintf("  %s%s\n", module_names[mod], status))
  }
  
  cat("\n═══════════════════════════════════════════════════════════════\n")
  cat("DRY RUN COMPLETE\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  cat("To execute this configuration, remove --dry_run flag:\n")
  cat("  ./duke [same arguments without --dry_run]\n\n")
  
  quit(status = 0)
}

cat("Configuration:\n")
cat("  Data directory:", params$dir_data, "\n")
cat("  Output directory:", params$dir_out, "\n")
cat("  Reference:", params$path_ref, "\n")
cat("  Threads:", params$threads, "\n")
cat("  Resume:", params$resume, "\n\n")

library(rmarkdown)
lib_files <- list.files("lib", pattern = "\\.R$", full.names = TRUE)
for (lf in lib_files) source(lf)

if (!dir.exists(params$dir_out)) dir.create(params$dir_out, recursive = TRUE)

# ==============================================================================
# Run Modules
# ==============================================================================

tryCatch({
  # Module names and files
  module_names <- c("Import and QC", "Alignment", "Repeat Detection",
                    "Allele Calling", "Waterfall Plots", "Range Analysis",
                    "Repeat Visualisation")
  
  rmd_files <- c("modules/01_import_and_qc.Rmd", "modules/02_alignment.Rmd", "modules/03_repeat_detection.Rmd",
                 "modules/04_allele_calling.Rmd", "modules/05_waterfall.Rmd", "modules/06_range_analysis.Rmd",
                 "modules/07_repeat_visualisation.Rmd")
  
  output_files <- c("01_import_qc_results.RData", "02_alignment_results.RData",
                    "03_repeat_detection_results.RData", "04_allele_calling_results.RData",
                    "05_waterfall_results.RData", "06_range_analysis_results.RData",
                    "07_repeat_visualisation_results.RData")
  
  # Only run selected modules
  for (mod in run_modules) {
    cat("\n-----------------------------------------------------------------\n")
    cat("Module", mod, ":", module_names[mod], "\n")
    cat("-----------------------------------------------------------------\n\n")
    
    output <- file.path(params$dir_out, "module_data", output_files[mod])
    
    # Skip Module 5 if waterfall disabled
    if (mod == 5 && !params$waterfall) {
      cat("Skipping Module 5 (waterfall = FALSE)...\n")
      next
    }
    
    # Skip if resume enabled and output exists
    if (params$resume && file.exists(output)) {
      cat("Skipping Module", mod, "(results found)...\n")
    } else {
      rmarkdown::render(rmd_files[mod], output_dir = params$dir_out, 
                        params = params, envir = new.env(),
                        knit_root_dir = getwd())
      cat("\nModule", mod, "complete!\n")
    }
  }
  
  # Cleanup
  if (params$cleanup_temp) {
    temp_dir <- file.path(params$dir_out, "temp")
    if (dir.exists(temp_dir)) {
      unlink(temp_dir, recursive = TRUE)
      cat("\nRemoved temp directory\n")
    }
  }
  
  end_time <- Sys.time()
  duration_secs <- as.numeric(difftime(end_time, start_time, units = "secs"))
  duration_hours <- floor(duration_secs / 3600)
  duration_mins <- floor((duration_secs %% 3600) / 60)
  duration_secs_remainder <- round(duration_secs %% 60)
  
  cat("\n=================================================================\n")
  cat("                    PIPELINE COMPLETED                           \n")
  cat("=================================================================\n\n")
  cat("Run finished:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
  cat(sprintf("Duration: %02d:%02d:%02d (HH:MM:SS)\n\n", duration_hours, duration_mins, duration_secs_remainder))
  
}, error = function(e) {
  cat("\n=================================================================\n")
  cat("                    PIPELINE FAILED                              \n")
  cat("=================================================================\n\n")
  cat("Error:", conditionMessage(e), "\n\n")
  sink(type = "message"); sink(type = "output"); close(log_con)
  quit(status = 1)
})

sink(type = "message"); sink(type = "output"); close(log_con)
quit(status = 0)
