#!/usr/bin/env Rscript
# ==============================================================================
# Duke Pipeline - Enhanced Runner Script
# ==============================================================================
# Runs the modular Duke pipeline for amplicon sequencing analysis
# ==============================================================================

# Set the working directory where this script is saved
# Works from both RStudio and command line
if (interactive()) {
  # Running in RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  } else {
    # Running interactively but not in RStudio
    cat("Note: Running interactively. Make sure you're in the pipeline directory.\n")
  }
} else {
  # Running from command line (Rscript)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg)
    setwd(dirname(normalizePath(script_path)))
  } else {
    cat("Note: Could not detect script location. Make sure you're in the pipeline directory.\n")
  }
}

# Clear objects
rm(list = setdiff(ls(), "env_duke"))


# Load required libraries
library(rmarkdown)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

params <- list(
  
  # -----------------------------------------------------------------------------
  # File paths and directories
  # -----------------------------------------------------------------------------
  dir_data = "/home/skgtmdf/Scratch/data/2025.12.17_pb_test/data",
  dir_out = "/home/skgtmdf/Scratch/data/2025.12.17_pb_test/result_duke",
  
  path_ref = "/home/skgtmdf/Scratch/refs/HTTset20/HTTset20.fasta",  # Reference with NNNNN separator
  
  # Optional paths
  path_trim_patterns = "/home/skgtmdf/Scratch/refs/adapters/adapters.csv",
  path_manual_exclusions = NA,  # Path to file with reads to exclude
  
  # -----------------------------------------------------------------------------
  # File import options
  # -----------------------------------------------------------------------------
  import_patterns = c("\\.fastq$", "\\.fastq.gz$", "\\.fasta$", "\\.fa$", "\\.bam$"),
  import_recursive = TRUE,
  
  # Paired-end read handling
  r1_pattern = "_R1_",
  r2_pattern = "_R2_",
  select_one_of_pair = "R1",  # Keep only R1, or set to NA to keep all
  
  # Downsampling (set to NA to disable)
  downsample = NA,  # e.g., 1000 to keep 1000 reads per sample
  
  # -----------------------------------------------------------------------------
  # Adapter trimming (from both 5' and 3' ends using Biostrings::trimLRPatterns)
  # -----------------------------------------------------------------------------
  
  # trim_max_mismatch: Maximum mismatches allowed when matching adapters
  #   - 0 = Perfect match only
  #   - 3 = Good for ONT (~10% error for 30bp adapter)
  #   - 5 = Very tolerant (risk of false matches)
  #   Example: "AGATCGGAAGAGC" with max_mismatch=3 matches "AGATCGGTAGAGC"
  
  # trim_with_indels: Allow insertions/deletions (gaps) in adapter matching
  #   - TRUE = More sensitive, better for ONT (finds adapters with indel errors)
  #   - FALSE = Faster, better for Illumina/PacBio HiFi (no indel tolerance)
  #   Example: "AGATCGG--AAGAGC" (2bp deletion) matches if TRUE
  
  trim_max_mismatch = 3,
  trim_with_indels = TRUE,
  
  # -----------------------------------------------------------------------------
  # Alignment with minimap2
  # -----------------------------------------------------------------------------
  minimap2_args = "-t 2 -x sr -w 1",
  # Using short-read mode (-x sr) with maximum sensitivity (-w 1) for amplicon sequencing
  # This setting works optimally for variable flank lengths (short to long amplicons)
  # 
  # Why -x sr for ONT amplicons?
  #   - Amplicon targets are short sequences (flanks typically 100-1500bp)
  #   - Short-read mode provides better sensitivity than long-read mode for these targets
  #   - Works regardless of sequencing platform (ONT, PacBio, or Illumina)
  # 
  # Why -w 1?
  #   - Minimizer window size = 1 gives maximum sensitivity
  #   - Ensures alignment detection for both short (<300bp) and long (>800bp) flanks
  #   - Critical for detecting all reads in amplicon data
  # 
  # Threading: -t <int> sets minimap2 internal threads (separate from params$threads)
  # 
  # Required flags (added automatically): -a (SAM output), --MD (mismatches), --cs (CIGAR)
  # 
  # Alternative presets for other applications:
  #   "-t 2 -x map-ont -k 15"  # For long genomic reads (>1kb), not amplicons
  #   "-t 2 -x map-hifi"       # For PacBio HiFi genomic data
  #   "-t 2 -x map-pb"         # For PacBio CLR genomic data
  
  # Visualisation options (new in enhanced Module 2)
  visualise_alignment = TRUE,  # Plot raw alignments (can be slow for large datasets)
  visualise_alignment_corrected = TRUE,  # Plot strand-corrected alignments
  
  # -----------------------------------------------------------------------------
  # Repeat detection and counting
  # -----------------------------------------------------------------------------
  rpt_pattern = "CAG",
  rpt_min_repeats = 2,
  rpt_max_mismatch = 0,
  rpt_start_perfect_repeats = 2,
  rpt_end_perfect_repeats = 2,
  rpt_max_gap = 6,
  rpt_max_tract_gap = 18,
  rpt_return_option = "longest",
  
  # Repeat counting method (used for clustering and distribution analysis):
  # 
  # Three methods are calculated, but only one is used as the primary count:
  # 
  # "repeat_count_full" (RECOMMENDED for ONT data):
  #   - Formula: Tract length ÷ pattern length
  #   - Example: 51 bases ÷ 3bp = 17 CAG repeats
  #   - Includes: All bases in the tract (tolerates sequencing errors)
  #   - Best for: Oxford Nanopore data, robust counting after tract isolation
  # 
  # "repeat_count_match":
  #   - Formula: Count of exact pattern matches only
  #   - Example: Counts only perfect "CAG" sequences, ignoring errors
  #   - Excludes: Sequencing errors and interruptions
  #   - Best for: High-accuracy data (PacBio HiFi, Illumina), stringent counting
  # 
  # "repeat_count_tracts":
  #   - Formula: Reports structure as "3,2,2" (tract sizes separated by gaps)
  #   - Example: "3,2,2" means three tracts containing 3, 2, and 2 perfect repeats
  #   - Note: Returns a STRING, not a number - cannot be used as repeat_count_method
  #   - Available in data for QC visualization only (any gap >0bp breaks a tract)
  #
  repeat_count_method = "repeat_count_full",
  
  # How to handle reads where no repeat tract was found (repeat_count = NA):
  # 
  # These could be:
  #   - Sequencing artifacts (e.g., duplicated flanks, chimeric reads)
  #   - Genuine deletions (e.g., CRISPR-edited samples with repeat removed)
  #   - Off-target amplification (non-HTT sequences)
  # 
  # "convert_to_zero":
  #   - Treat as 0 repeats (no tract present)
  #   - Best for: CRISPR-edited samples where deletions are expected
  #   - Preserves all reads but includes artifacts in the 0-count bin
  #   - Default (matches original duke.Rmd behavior)
  # 
  # "filter":
  #   - Remove these reads entirely
  #   - Best for: Non-edited samples, prioritising data quality over completeness
  #   - Conservative approach that removes artifacts but also genuine deletions
  # 
  # "flag_only":
  #   - Keep NA as-is, add 'no_tract_found' flag column
  #   - Best for: Exploratory analysis, want to inspect these reads manually
  #   - Allows filtering later on a per-sample basis
  #   - Note: Downstream analyses must handle NAs
  #
  na_repeat_handling = "convert_to_zero",
  
  # -----------------------------------------------------------------------------
  # Clustering
  # -----------------------------------------------------------------------------
  cluster = TRUE,
  
  # Clustering strategy (character vector):
  #   c("repeat", "haplotype") - Cluster by repeat first, then haplotype within each repeat group (RECOMMENDED)
  #   c("haplotype", "repeat") - Cluster by haplotype first, then repeat within each haplotype group
  #   "repeat"                 - Cluster by repeat count only (fast)
  #   "haplotype"              - Cluster by flanking sequences only
  #   "none"                   - No clustering (all reads in single group)
  #
  # Vector order specifies clustering order:
  #   - c("repeat", "haplotype") = repeat-first strategy
  #   - c("haplotype", "repeat") = haplotype-first strategy
  #
  # Why repeat-first is recommended for long reads (ONT/PacBio):
  #   - Repeat counting is more robust to sequencing errors
  #   - Prevents error amplification in haplotype clustering
  #   - Enforces biological constraint: one haplotype per repeat length (diploid)
  #
  # Examples:
  #   cluster_by = c("repeat", "haplotype")  # Two-step: repeat then haplotype
  #   cluster_by = "repeat"                  # Single: repeat only (fastest)
  #   cluster_by = "haplotype"               # Single: haplotype only
  cluster_by = "repeat",
  
  # Separate cluster limits for repeat and haplotype:
  haplotype_cluster_max = 10,   # Max haplotypes to consider (auto-detect)
  repeat_cluster_max = 20,      # Max repeats to consider (allow more for instability)
  
  # Haplotype clustering options (used when cluster_by includes "haplotype")
  haplotype_region = "both",        # "left", "right", or "both" flanking regions
  haplotype_method = "levenshtein", # "hamming" or "levenshtein" distance
  haplotype_trim_length = "auto",   # Trimming options:
                                     #   NA = no trim (uses full flanks)
                                     #   "auto" = trim to modal (commonest) length per flank
                                     #   numeric (e.g., 200) = trim to N bp total (split between flanks)
                                     # Trimming keeps bases adjacent to repeat (most informative)
  haplotype_subsample = 250,         # Subsampling for speed (large datasets):
                                     #   NA = no subsampling (compute full distance matrix)
                                     #   5000 = cluster 5000 representatives, assign remainder
                                     # Speeds up clustering for >5000 sequences
  
  # -----------------------------------------------------------------------------
  # Consensus generation and variant calling
  # -----------------------------------------------------------------------------
  cluster_consensus = TRUE,
  consensus_threshold = 0.5,   # Minimum proportion for consensus base (0.5 = majority)
  consensus_downsample = 50,  # Max reads per consensus (NA = use all)
  
  # Variant calling (compare consensus to reference)
  call_variants = TRUE,        # Call variants in consensus sequences
  
  # -----------------------------------------------------------------------------
  # Waterfall plots (Module 5)
  # -----------------------------------------------------------------------------
  waterfall = TRUE,  # Generate waterfall plots (visual read inspection)
  
  # waterfall_downsample: Maximum reads to plot per sample
  #   - 5000 = plot up to 5000 reads (recommended for performance)
  #   - NA = plot all reads (may be slow for large datasets)
  #   Note: Downsampling is random but reproducible (seed = 123)
  waterfall_downsample = 1000,
  
  # waterfall_rm_flank_length_outliers: Remove reads with unusual flank lengths
  #   - TRUE = remove outliers using 1.5 × IQR method per cluster (recommended)
  #   - FALSE = plot all reads regardless of flank length
  #   Helps remove chimeric reads and sequencing artifacts
  waterfall_rm_flank_length_outliers = TRUE,
  
  # waterfall_y_axis_labels: Controls y-axis labeling density
  #   - "auto" = dynamic based on read count (default, from original Duke)
  #              <30 reads: every read, 30-50: every 3rd, 50-100: every 5th,
  #              100-200: every 10th, 200-500: every 20th, 500-1000: every 50th,
  #              1000-5000: every 100th, >5000: every 500th
  #   - "all" = label every single read (only for small datasets)
  #   - Integer (e.g. 10, 20, 50) = target number of labels
  waterfall_y_axis_labels = "auto",
  
  # waterfall_per_cluster: Generate separate waterfall plot for each cluster
  #   - TRUE = create per-cluster plots (can create many files)
  #   - FALSE = only per-sample plots (default, faster)
  waterfall_per_cluster = TRUE,
  
  # DNA base colors for waterfall plots
  dna_colours = c("A" = "darkgreen", "C" = "blue", "T" = "red", 
                  "G" = "black", "-" = "lightgrey", "N" = "grey"),
  
  # -----------------------------------------------------------------------------
  # Range analysis & instability metrics (Module 6)
  # -----------------------------------------------------------------------------
  # path_settings: Path to settings Excel file
  #   Required for Module 6
  #   Should contain: analysis_ranges, floor, max_peaks, group, group_control_sample
  path_settings = "/home/skgtmdf/Scratch/data/2025.12.17_pb_test/settings/settings_duke.xlsx",
  
  # range_peak_span: Peak detection smoothing window (odd number)
  #   - 3 = sharp peaks (default)
  #   - 5-7 = smoother peaks, less sensitive to noise
  range_peak_span = 3,
  
  # group_control: Enable group control analysis
  #   - TRUE = calculate control setpoints per group
  #   - FALSE = skip instability metrics
  group_control = TRUE,
  
  # control_sample_selection: Method for selecting which samples to use as controls
  #   - "flagged" = use samples with group_control_sample = 1 (default)
  #   - "earliest" = use first timepoint (time == min)
  #   - "all" = use all samples in group
  control_sample_selection = "flagged",
  
  # control_setpoint_metric: Which summary statistic to use as the control setpoint
  #   - "modal_length" = most frequent repeat length (default, robust)
  #   - "mean_length" = arithmetic mean (sensitive to outliers)
  #   - "median_length" = median repeat length (robust middle value)
  control_setpoint_metric = "median_length",
  
  # control_aggregation_method: How to aggregate control samples' values into a single setpoint
  #   - "mean" = arithmetic mean of control values (default)
  #   - "median" = median of control values (robust to outliers)
  #   - "trimmed_mean" = 10% trimmed mean (excludes extreme 10%)
  control_aggregation_method = "median",
  
  # -----------------------------------------------------------------------------
  # Module 7: Repeat Distribution Visualisation
  # -----------------------------------------------------------------------------
  # repeat_histogram: Generate frequency histogram plots
  repeat_histogram = TRUE,
  
  # repeat_histogram_binwidth: Bin width for histogram (in repeats)
  repeat_histogram_binwidth = 1,
  
  # repeat_scatter: Generate scatter/violin plots showing individual reads
  repeat_scatter = TRUE,
  
  # repeat_distribution_metrics: Which summary metrics to show as points/markers
  #   Options: "modal_length", "mean_length", "median_length"
  #   Can specify multiple, e.g., c("modal_length", "mean_length")
  repeat_distribution_metrics = c("modal_length", "mean_length", "median_length"),
  
  # -----------------------------------------------------------------------------
  # Runtime settings
  # -----------------------------------------------------------------------------
  threads = 6,
  resume = TRUE,
  
  # -----------------------------------------------------------------------------
  # Memory and disk management
  # -----------------------------------------------------------------------------
  # These parameters control cleanup during and after the pipeline run
  
  # remove_intermediate: Free up RAM during pipeline execution
  #   - Deletes large R objects (sequences, alignments) from memory after use
  #   - Applied in Modules 1-4 where memory usage is highest
  #   - Useful for: large datasets (many samples/deep coverage) or low-RAM systems
  #   - Safe to enable: data is saved to disk before deletion
  remove_intermediate = TRUE,
  
  # cleanup_temp: Free up disk space after pipeline completion  
  #   - Deletes temp/ directory containing intermediate RData files
  #   - Applied at end of pipeline after all modules successfully complete
  #   - Useful for: repeated runs, limited disk space, production workflows
  #   - Keep FALSE for: debugging, inspecting intermediate results
  cleanup_temp = FALSE
)

# ==============================================================================
# LOGGING SETUP
# ==============================================================================

# Create logs directory at root level (same location as duke_run.R)
log_dir <- file.path(getwd(), "logs")
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}

# Generate timestamp and create subdirectory for this run
log_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_log_dir <- file.path(log_dir, log_timestamp)
if (!dir.exists(run_log_dir)) {
  dir.create(run_log_dir, recursive = TRUE)
}

# Log files in subdirectory WITH timestamp in filename
log_file <- file.path(run_log_dir, paste0(log_timestamp, "_duke_run.log"))
params_file <- file.path(run_log_dir, paste0(log_timestamp, "_params.R"))

# Save parameters to file for reproducibility
cat("# Duke Pipeline Run Parameters\n", file = params_file)
cat("# Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", 
    file = params_file, append = TRUE)
cat("params <- ", file = params_file, append = TRUE)
capture.output(dput(params), file = params_file, append = TRUE)

# Start capturing output to log file
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)  # split=TRUE keeps console output
sink(log_con, type = "message")

# Record start time
start_time <- Sys.time()

# ==============================================================================
# RUN PIPELINE
# ==============================================================================

cat("\n")
cat("=================================================================\n")
cat("                    DUKE PIPELINE 2.0                            \n")
cat("=================================================================\n")
cat("\n")
cat("Run started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Log file:", log_file, "\n")
cat("Parameters saved:", params_file, "\n")
cat("\n")
cat("Configuration:\n")
cat("  Data directory:", params$dir_data, "\n")
cat("  Output directory:", params$dir_out, "\n")
cat("  Reference:", params$path_ref, "\n")
cat("  Threads:", params$threads, "\n")
cat("  Resume:", params$resume, "\n")
cat("  Cleanup temp:", params$cleanup_temp, "\n")
cat("  Remove intermediate:", params$remove_intermediate, "\n")
cat("\n")

# Create output directory
if (!dir.exists(params$dir_out)) {
  dir.create(params$dir_out, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# Module 1: Import and QC
# -----------------------------------------------------------------------------
cat("\n")
cat("-----------------------------------------------------------------\n")
cat("Module 1: Import and QC\n")
cat("-----------------------------------------------------------------\n")

module1_output <- file.path(params$dir_out, "module_data", "01_import_qc_results.RData")

if (params$resume && file.exists(module1_output)) {
  cat("Skipping Module 1 (results found)...\n")
} else {
  rmarkdown::render(
    input = "01_import_and_qc.Rmd",
    output_dir = params$dir_out,
    params = params,
    envir = new.env()
  )
  cat("\nModule 1 complete!\n")
}

# -----------------------------------------------------------------------------
# Module 2: Alignment and Processing
# -----------------------------------------------------------------------------
cat("\n")
cat("-----------------------------------------------------------------\n")
cat("Module 2: Alignment and Processing\n")
cat("-----------------------------------------------------------------\n")

module2_output <- file.path(params$dir_out, "module_data", "02_alignment_results.RData")

if (params$resume && file.exists(module2_output)) {
  cat("Skipping Module 2 (results found)...\n")
} else {
  rmarkdown::render(
    input = "02_alignment.Rmd",
    output_dir = params$dir_out,
    params = params,
    envir = new.env()
  )
  cat("\nModule 2 complete!\n")
}

# -----------------------------------------------------------------------------
# Module 3: Repeat Detection
# -----------------------------------------------------------------------------
cat("\n")
cat("-----------------------------------------------------------------\n")
cat("Module 3: Repeat Detection\n")
cat("-----------------------------------------------------------------\n")

module3_output <- file.path(params$dir_out, "module_data", "03_repeat_detection_results.RData")

if (params$resume && file.exists(module3_output)) {
  cat("Skipping Module 3 (results found)...\n")
} else {
  rmarkdown::render(
    input = "03_repeat_detection.Rmd",
    output_dir = params$dir_out,
    params = params,
    envir = new.env()
  )
  cat("\nModule 3 complete!\n")
}

# -----------------------------------------------------------------------------
# Module 4: Allele Calling
# -----------------------------------------------------------------------------
cat("\n")
cat("-----------------------------------------------------------------\n")
cat("Module 4: Allele Calling\n")
cat("-----------------------------------------------------------------\n")

module4_output <- file.path(params$dir_out, "module_data", "04_allele_calling_results.RData")

if (params$resume && file.exists(module4_output)) {
  cat("Skipping Module 4 (results found)...\n")
} else {
  rmarkdown::render(
    input = "04_allele_calling.Rmd",
    output_dir = params$dir_out,
    params = params,
    envir = new.env()
  )
  cat("\nModule 4 complete!\n")
}

# -----------------------------------------------------------------------------
# Module 5: Waterfall Plots (Optional)
# -----------------------------------------------------------------------------
if (params$waterfall) {
  cat("\n")
  cat("-----------------------------------------------------------------\n")
  cat("Module 5: Waterfall Plots\n")
  cat("-----------------------------------------------------------------\n")
  
  module5_output <- file.path(params$dir_out, "module_data", "05_waterfall_results.RData")
  
  if (params$resume && file.exists(module5_output)) {
    cat("Skipping Module 5 (results found)...\n")
  } else {
    rmarkdown::render(
      input = "05_waterfall.Rmd",
      output_dir = params$dir_out,
      params = params,
      envir = new.env()
    )
    cat("\nModule 5 complete!\n")
  }
} else {
  cat("\n")
  cat("Skipping Module 5 (waterfall = FALSE)...\n")
}

# -----------------------------------------------------------------------------
# Module 6: Range Analysis & Instability Metrics
# -----------------------------------------------------------------------------
cat("\n")
cat("-----------------------------------------------------------------\n")
cat("Module 6: Range Analysis & Instability Metrics\n")
cat("-----------------------------------------------------------------\n")

module6_output <- file.path(params$dir_out, "module_data", "06_range_analysis_results.RData")

if (params$resume && file.exists(module6_output)) {
  cat("Skipping Module 6 (results found)...\n")
} else {
  rmarkdown::render(
    input = "06_range_analysis.Rmd",
    output_dir = params$dir_out,
    params = params,
    envir = new.env()
  )
  cat("\nModule 6 complete!\n")
}

# -----------------------------------------------------------------------------
# Module 7: Repeat Distribution Visualisation
# -----------------------------------------------------------------------------
cat("\n")
cat("-----------------------------------------------------------------\n")
cat("Module 7: Repeat Distribution Visualisation\n")
cat("-----------------------------------------------------------------\n")

module7_output <- file.path(params$dir_out, "module_data", "07_repeat_visualisation_results.RData")

if (params$resume && file.exists(module7_output)) {
  cat("Skipping Module 7 (results found)...\n")
} else {
  rmarkdown::render(
    input = "07_repeat_visualisation.Rmd",
    output_dir = params$dir_out,
    params = params,
    envir = new.env()
  )
  cat("\nModule 7 complete!\n")
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n")
cat("=================================================================\n")
cat("                      PIPELINE COMPLETE                          \n")
cat("=================================================================\n")
cat("\n")
cat("Results saved to:", params$dir_out, "\n")
cat("\n")
cat("HTML reports:\n")
cat("  - ", file.path(params$dir_out, "01_import_and_qc.html"), "\n")
cat("  - ", file.path(params$dir_out, "02_alignment.html"), "\n")
cat("  - ", file.path(params$dir_out, "03_repeat_detection.html"), "\n")
cat("  - ", file.path(params$dir_out, "04_allele_calling.html"), "\n")
if (params$waterfall) {
  cat("  - ", file.path(params$dir_out, "05_waterfall.html"), "\n")
}
cat("  - ", file.path(params$dir_out, "06_range_analysis.html"), "\n")
cat("  - ", file.path(params$dir_out, "07_repeat_visualisation.html"), "\n")
cat("\n")
cat("RData files:\n")
cat("  - ", module1_output, "\n")
cat("  - ", module2_output, "\n")
cat("  - ", module3_output, "\n")
cat("  - ", module4_output, "\n")
if (params$waterfall && exists("module5_output")) {
  cat("  - ", module5_output, "\n")
}
cat("  - ", module6_output, "\n")
cat("  - ", module7_output, "\n")
cat("\n")
cat("Excel exports:\n")
cat("  - ", file.path(params$dir_out, "06_range_analysis", "range_analysis_results.xlsx"), "\n")
cat("\n")

# ==============================================================================
# CLEANUP
# ==============================================================================

if (params$cleanup_temp) {
  cat("\n")
  cat("=================================================================\n")
  cat("                      CLEANUP                                     \n")
  cat("=================================================================\n")
  cat("\n")
  
  temp_dir <- file.path(params$dir_out, "temp")
  
  if (dir.exists(temp_dir)) {
    cat("Removing temporary files...\n")
    
    # Get size before deletion
    temp_files <- list.files(temp_dir, recursive = TRUE, full.names = TRUE)
    temp_size <- sum(file.size(temp_files), na.rm = TRUE) / 1024^2  # MB
    
    # Remove temp directory
    unlink(temp_dir, recursive = TRUE)
    
    cat(sprintf("  Removed %.1f MB from temp directory\n", temp_size))
    cat("  Temp directory deleted\n")
  } else {
    cat("No temp directory found\n")
  }
  
  cat("\n")
}

cat("\n")
cat("=================================================================\n")
cat("                      ALL DONE!                                  \n")
cat("=================================================================\n")
cat("\n")

# ==============================================================================
# FINISH LOGGING
# ==============================================================================

# Calculate and log duration (before closing log file)
end_time <- Sys.time()
duration <- end_time - start_time

cat("Run finished:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total duration:", format(duration), "\n")
cat("\n")

# Stop capturing output
sink(type = "message")
sink(type = "output")
close(log_con)

# Final summary to console only
cat("\n")
cat("Pipeline completed successfully!\n")
cat("Log file saved:", log_file, "\n")
cat("Parameters saved:", params_file, "\n")
cat("\n")
