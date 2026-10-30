# ==============================================================================
# Duke - parameter specification
# ==============================================================================
# One table defines every parameter: its default, type, the module it belongs
# to and its help text. duke_run.R overrides the defaults, the CLI builds its
# options from this table, and run_duke() validates against it - so the three
# can never drift apart.
#
# Types: chr, int, num, lgl (scalars); chrvec, intvec, numvec (vectors;
# comma-separated on the command line).
# ==============================================================================

.p <- function(name, default, type, group, help) {
  list(name = name, default = default, type = type, group = group, help = help)
}

duke_param_spec <- function() {
  list(
    # --- Paths -----------------------------------------------------------------
    .p("dir_data", NA_character_, "chr", "Paths", "Input directory of BAM / FASTQ / FASTA files (required)"),
    .p("dir_out", NA_character_, "chr", "Paths", "Output directory (required)"),
    .p("path_ref", NA_character_, "chr", "Paths",
       "Reference FASTA with NNNNN marking the repeat region (required unless every tract in path_tracts names its own reference)"),
    .p("path_settings", NA_character_, "chr", "Paths", "Settings .xlsx/.csv (ranges, groups, controls, time). Needed for Module 6"),
    .p("path_tracts", NA_character_, "chr", "Paths",
       "Optional tracts table (.csv or .xlsx sheet 'tracts'; optional sheet 'composites'). NA = one tract from rpt_pattern and path_ref"),
    .p("path_trim_patterns", NA_character_, "chr", "Paths", "Adapter CSV (adapter_name, adapter_sequence), used when trim = TRUE"),
    .p("path_manual_exclusions", NA_character_, "chr", "Paths", "Optional read exclusions file (file_name, read_name, reason)"),

    # --- Module 1: import and QC ---------------------------------------------------
    .p("import_patterns", c("\\.fastq$", "\\.fastq\\.gz$", "\\.fq$", "\\.fq\\.gz$", "\\.fasta$", "\\.fa$", "\\.bam$"),
       "chrvec", "Import", "Regex patterns for files to import"),
    .p("import_recursive", TRUE, "lgl", "Import", "Search dir_data recursively"),
    .p("r1_pattern", "_R1_", "chr", "Import", "Pattern identifying R1 files in paired-end data"),
    .p("r2_pattern", "_R2_", "chr", "Import", "Pattern identifying R2 files in paired-end data"),
    .p("select_one_of_pair", "R1", "chr", "Import", "Paired-end: import 'R1', 'R2' or 'all'"),
    .p("downsample", NA_integer_, "int", "Import", "Maximum reads per sample (NA = all)"),
    .p("np_min", NA_integer_, "int", "Import", "Drop PacBio HiFi reads with fewer passes than this (NA = no filter)"),
    .p("pacbio_scatter_downsample", 50000L, "int", "Import", "Maximum points in the PacBio passes-vs-length plot"),
    .p("check_duplicate_readnames", FALSE, "lgl", "Import",
       "Remove duplicate read names, keeping the longest. Off by default: PacBio and ONT read names are unique, and the check removed 0 reads in 36 M"),
    .p("trim", FALSE, "lgl", "Import",
       "Trim adapters from read ends. Off by default: demultiplexed HiFi reads are already trimmed (20 of 36 M reads were trimmed in the msn run)"),
    .p("trim_max_mismatch", 3L, "int", "Import", "Maximum errors in an adapter match"),
    .p("trim_with_indels", TRUE, "lgl", "Import", "Allow indels in adapter matching"),

    # --- Module 2: alignment ---------------------------------------------------------
    .p("minimap2_args", "-t 1 -x sr -w 1", "chr", "Alignment", "minimap2 arguments for aligning reads to each flank"),
    .p("keep_bam", FALSE, "lgl", "Alignment", "Keep the per-sample flank BAMs in 02_alignment/bam/ for IGV"),
    .p("visualise_alignment", TRUE, "lgl", "Alignment", "Write per-sample flank alignment plots"),
    .p("visualise_alignment_downsample", 1000L, "int", "Alignment", "Maximum reads per alignment plot"),

    # --- Module 3: repeat detection --------------------------------------------------
    .p("rpt_pattern", "CAG", "chr", "Repeats", "Repeat motif (used when there is no tracts table)"),
    .p("rpt_tract_name", NA_character_, "chr", "Repeats", "Name for the single tract (default: the motif)"),
    .p("rpt_min_repeats", 2L, "int", "Repeats", "Minimum tract length in repeat units"),
    .p("rpt_max_mismatch", 0L, "int", "Repeats", "Mismatches allowed per repeat unit"),
    .p("rpt_start_perfect_repeats", 2L, "int", "Repeats", "Perfect units required at the start of a tract"),
    .p("rpt_end_perfect_repeats", 2L, "int", "Repeats", "Perfect units required at the end of a tract"),
    .p("rpt_max_gap", 6L, "int", "Repeats", "Maximum bases between units inside a tract"),
    .p("rpt_max_tract_gap", 18L, "int", "Repeats", "Maximum gap between sub-tracts before they are merged"),
    .p("rpt_return_option", "longest", "chr", "Repeats", "Which tract to keep when several are found: 'longest' or 'first'"),
    .p("rpt_flank_pad", 30L, "int", "Repeats",
       "Bases either side of the flank boundaries searched for the tract, so a tract can extend past where the aligner placed the flank"),
    .p("repeat_count_method", "full", "chr", "Repeats",
       "Primary count: 'full' (tract length / motif length), 'match' (exact units) or 'pure' (longest uninterrupted run)"),
    .p("na_repeat_handling", "keep", "chr", "Repeats",
       "Reads that align but have no tract: 'keep' (no length; counted in the detection rate only), 'convert_to_zero' (length 0, as in v2; falls into ranges starting at 0 or NA) or 'filter' (removed)"),
    .p("composites", character(0), "chrvec", "Repeats",
       "Combined measures, e.g. 'CAG_CCG_sum=CAG+CCG:sum' (pure units added) or 'CAG_CCG_span=CAG+CCG:span' (first tract start to last tract end)"),
    .p("export_read_counts", TRUE, "lgl", "Repeats", "Write per-read lengths to 03_repeat_detection/read_repeat_counts/"),
    .p("rm_flank_length_outliers", FALSE, "lgl", "Repeats", "Remove reads with outlying combined flank length"),
    .p("flank_iqr_multiplier", 1.5, "num", "Repeats", "IQR multiplier for the flank length filter"),

    # --- Module 4: UMI ---------------------------------------------------------------------
    .p("umi", FALSE, "lgl", "UMI", "Extract UMIs and collapse reads to molecules"),
    .p("umi_flank_5", "GTAAAACGACGGCCAGT", "chr", "UMI", "Sequence immediately 5' of the UMI (default M13-F)"),
    .p("umi_flank_3", "GGCGACCCTGGAAAAGCTGATGAA", "chr", "UMI", "Sequence immediately 3' of the UMI (default ciosiF+1)"),
    .p("umi_flank_max_mismatch", 3L, "int", "UMI", "Edits allowed when locating each UMI flank (indel tolerant)"),
    .p("umi_search_window", 250L, "int", "UMI", "The 5' flank must start within this many bases of the read start"),
    .p("umi_len_min", 12L, "int", "UMI", "Shortest UMI accepted"),
    .p("umi_len_max", 16L, "int", "UMI", "Longest UMI accepted"),
    .p("umi_cluster_max_dist", 1L, "int", "UMI", "Edit distance for directional UMI clustering (0 = exact match; 1 or 2)"),
    .p("umi_min_reads", 3L, "int", "UMI", "Minimum reads for a molecule to be kept"),
    .p("umi_collapse_stat", "top_cluster", "chr", "UMI", "Rule that turns the reads of a molecule into its one length: 'top_cluster', 'max_supported', 'mode', 'median' or 'mean' (all are compared in the Module 4 report)"),
    .p("umi_top_cluster_k", 2L, "int", "UMI", "top_cluster: reads required in the window"),
    .p("umi_top_cluster_w", 3, "num", "UMI", "top_cluster: constant part of the window half-width (repeats)"),
    .p("umi_top_cluster_wfrac", 0.10, "num", "UMI", "top_cluster: length-proportional part of the window half-width"),
    .p("umi_top_cluster_min_share", 0.2, "num", "UMI", "top_cluster: share of the molecule's reads the window must also hold, so a few stray long reads cannot set the length (0 = only umi_top_cluster_k reads needed)"),
    .p("umi_long_read_filter", TRUE, "lgl", "UMI", "Drop molecules where many reads are much longer than the molecule (most likely a longer template whose reads were shortened in PCR)"),
    .p("umi_long_read_units", 20, "num", "UMI", "Long-read filter: how many units longer than its molecule a read must be to count as long"),
    .p("umi_long_read_share", 0.2, "num", "UMI", "Long-read filter: share of the reads of a molecule that must be long for it to be flagged"),

    # --- Module 5: allele structure ---------------------------------------------------------
    .p("cluster", FALSE, "lgl", "Allele structure", "Cluster reads into alleles, build consensus sequences and call variants"),
    .p("cluster_by", "repeat", "chrvec", "Allele structure", "'repeat', 'haplotype', or 'repeat,haplotype'"),
    .p("repeat_cluster_max", 6L, "int", "Allele structure", "Maximum repeat-length clusters per sample"),
    .p("repeat_cluster_min_prop", 0.05, "num", "Allele structure", "Clusters holding fewer than this proportion of reads are merged into their nearest neighbour"),
    .p("repeat_cluster_downsample", 10000L, "int", "Allele structure", "Reads used to choose the number of repeat clusters"),
    .p("haplotype_cluster_max", 6L, "int", "Allele structure", "Maximum haplotype clusters per sample"),
    .p("haplotype_cluster_downsample", 3000L, "int", "Allele structure", "Reads used to define haplotype clusters"),
    .p("haplotype_region", "both", "chr", "Allele structure", "Flank(s) used for haplotypes: 'left', 'right' or 'both'"),
    .p("haplotype_trim_length", "auto", "chr", "Allele structure", "'auto', a number of bases, or 'NA' for no trimming"),
    .p("consensus_threshold", 0.5, "num", "Allele structure", "Minimum agreement for a consensus base"),
    .p("consensus_downsample", 50L, "int", "Allele structure", "Maximum reads per consensus alignment"),
    .p("call_variants", TRUE, "lgl", "Allele structure", "Call variants in each cluster consensus against the reference"),
    .p("waterfall", TRUE, "lgl", "Allele structure", "Write per-sample waterfall plots"),
    .p("waterfall_downsample", 1000L, "int", "Allele structure", "Maximum reads per waterfall"),
    .p("waterfall_flank_bp", 60L, "int", "Allele structure", "Bases of flank drawn either side of the tract in waterfalls"),
    .p("waterfall_per_cluster", FALSE, "lgl", "Allele structure", "Also write one waterfall per cluster"),

    # --- Module 6: instability --------------------------------------------------------------
    .p("range_peak_span", 3L, "int", "Instability", "Window for modal peak detection"),
    .p("group_control", TRUE, "lgl", "Instability", "Calculate group control setpoints"),
    .p("control_sample_selection", "flagged", "chr", "Instability", "'flagged', 'earliest' or 'all'"),
    .p("control_setpoint_metric", "median_length", "chr", "Instability", "'modal_length', 'mean_length' or 'median_length'"),
    .p("control_aggregation_method", "median", "chr", "Instability", "'mean', 'median' or 'trimmed_mean' across control samples"),
    .p("use_manual_controls", FALSE, "lgl", "Instability",
       "Use manual_control_repeat_length from the settings file as the control setpoint (replacing the group control) for those samples and ranges"),
    .p("expansion_thresholds", c(5L, 10L), "intvec", "Instability", "k for proportion of reads at least k units above / below the setpoint"),
    .p("cohort_plot_by", c("group", "time"), "chrvec", "Instability", "Dimensions for cohort plots: group, time"),
    .p("headline_metrics", c("modal_length", "median_length", "mean_length", "instability_index",
                             "expansion_index", "prop_above_5"),
       "chrvec", "Instability", "Metrics drawn as timecourses by group"),
    .p("headline_setpoint_type", "control_relative", "chr", "Instability", "Setpoint type for headline plots: 'control_relative' or 'sample_relative'"),
    .p("plot_ranges", character(0), "chrvec", "Instability", "Ranges drawn inline (default: the first two ranges of the primary measure)"),
    .p("measure_pairs", character(0), "chrvec", "Instability", "Pairs of measures to plot against each other per read, e.g. 'CAG:CCG' (default: every pair of tracts)"),
    .p("pair_max_points", 4000L, "int", "Instability", "Maximum reads per sample in per-read pair plots"),

    # --- Module 7: repeat lengths -----------------------------------------------------------
    .p("sample_cards", TRUE, "lgl", "Repeat lengths", "Write one summary figure per sample"),
    .p("repeat_histogram_binwidth", 1L, "int", "Repeat lengths", "Histogram bin width in repeat units"),

    # --- Output and runtime ---------------------------------------------------------------
    .p("plot_dpi", 150L, "int", "Runtime", "Resolution of saved plots"),
    .p("plot_per_sample", TRUE, "lgl", "Runtime", "Write per-sample plots in Modules 2 and 3"),
    .p("threads", 4L, "int", "Runtime", "CPU cores"),
    .p("resume", TRUE, "lgl", "Runtime", "Reuse a module's results when the parameters that affect it are unchanged"),
    .p("remove_intermediate", TRUE, "lgl", "Runtime", "Free memory between steps"),
    .p("remove_temp", FALSE, "lgl", "Runtime", "Delete temp/ (per-sample caches) at the end"),
    .p("run_modules", 1:7, "intvec", "Runtime", "Modules to run")
  )
}

duke_defaults <- function() {
  spec <- duke_param_spec()
  setNames(lapply(spec, `[[`, "default"), vapply(spec, `[[`, "", "name"))
}

# Parameters renamed or removed since v2.3.1, mapped to their replacements
.legacy_params <- list(
  repeat_histogram = NULL, repeat_scatter = NULL, repeat_density = NULL,
  repeat_distribution_metrics = NULL, cluster_consensus = NULL, dna_colours = NULL,
  waterfall_y_axis_labels = NULL, haplotype_method = NULL, visualise_alignment_n_reads = NULL,
  umi_allele_min_sep = NULL, umi_split_mixed = NULL
)

# Merge user params into defaults, coerce types and validate
resolve_params <- function(user = list()) {
  spec <- duke_param_spec()
  types <- setNames(vapply(spec, `[[`, "", "type"), vapply(spec, `[[`, "", "name"))
  p <- duke_defaults()
  unknown <- setdiff(names(user), names(p))
  legacy  <- intersect(unknown, names(.legacy_params))
  if (length(legacy)) message("Ignoring parameters no longer used: ", paste(legacy, collapse = ", "))
  unknown <- setdiff(unknown, legacy)
  if (length(unknown)) warning("Unknown parameters ignored: ", paste(unknown, collapse = ", "))
  for (n in intersect(names(user), names(p))) p[[n]] <- user[[n]]

  coerce <- function(x, type) {
    if (is.null(x)) return(x)
    switch(type,
      chr = if (length(x) == 0 || (length(x) == 1 && (is.na(x) || identical(x, "") || identical(toupper(x), "NA")))) NA_character_ else as.character(x)[1],
      int = if (length(x) == 0 || is.na(x[1]) || identical(toupper(as.character(x[1])), "NA")) NA_integer_ else as.integer(x[1]),
      num = if (length(x) == 0 || is.na(x[1]) || identical(toupper(as.character(x[1])), "NA")) NA_real_ else as.numeric(x[1]),
      lgl = { v <- toupper(trimws(as.character(x[1])))
              if (!v %in% c("TRUE", "FALSE", "T", "F", "1", "0", "YES", "NO")) stop("A TRUE/FALSE parameter was given '", x[1], "'")
              v %in% c("TRUE", "T", "1", "YES") },
      chrvec = { x <- as.character(x); if (length(x) == 1 && grepl("[,;]", x)) x <- trimws(strsplit(x, "[,;]")[[1]]); x[nzchar(x) & !is.na(x)] },
      intvec = { if (is.character(x)) x <- trimws(strsplit(paste(x, collapse = ","), ",")[[1]]); as.integer(x) },
      numvec = { if (is.character(x)) x <- trimws(strsplit(paste(x, collapse = ","), ",")[[1]]); as.numeric(x) },
      x)
  }
  for (n in names(p)) p[[n]] <- coerce(p[[n]], types[[n]])

  # Back-compatible count method names
  p$repeat_count_method <- sub("^repeat_count_", "", p$repeat_count_method)
  if (!p$repeat_count_method %in% c("full", "match", "pure"))
    stop("repeat_count_method must be 'full', 'match' or 'pure'")
  if (identical(p$na_repeat_handling, NA_character_)) p$na_repeat_handling <- "keep"
  p$cluster_by <- unlist(strsplit(p$cluster_by, ","))
  if (is.na(p$dir_data) || is.na(p$dir_out)) stop("dir_data and dir_out are required")
  if (isTRUE(p$umi) && p$umi_min_reads <= p$umi_top_cluster_k && p$umi_collapse_stat == "top_cluster")
    stop("umi_min_reads (", p$umi_min_reads, ") must exceed umi_top_cluster_k (", p$umi_top_cluster_k,
         "): with exactly k reads top_cluster returns their median, which at k = 2 is their mean")
  if (!p$umi_cluster_max_dist %in% 0:2) stop("umi_cluster_max_dist must be 0, 1 or 2")
  if (!p$umi_collapse_stat %in% c("top_cluster", "max_supported", "mode", "median", "mean"))
    stop("umi_collapse_stat must be 'top_cluster', 'max_supported', 'mode', 'median' or 'mean'")
  if (!is.finite(p$umi_long_read_units) || p$umi_long_read_units <= 0) stop("umi_long_read_units must be a positive number")
  if (!is.finite(p$umi_long_read_share) || p$umi_long_read_share <= 0 || p$umi_long_read_share > 1) stop("umi_long_read_share must be above 0 and at most 1")
  if (!is.finite(p$umi_top_cluster_min_share) || p$umi_top_cluster_min_share < 0 || p$umi_top_cluster_min_share > 1)
    stop("umi_top_cluster_min_share must be between 0 and 1")
  p$dir_out <- normalizePath(p$dir_out, mustWork = FALSE)
  p
}

# ------------------------------------------------------------------------------
# Tracts and composite measures
# ------------------------------------------------------------------------------
# Returns a data.frame with one row per tract, in 5'->3' order, and every
# repeat-finding setting resolved (tract value, else the global rpt_* value).
resolve_tracts <- function(p) {
  cols <- c(min_repeats = "rpt_min_repeats", max_mismatch = "rpt_max_mismatch",
            start_perfect_repeats = "rpt_start_perfect_repeats",
            end_perfect_repeats = "rpt_end_perfect_repeats", max_gap = "rpt_max_gap",
            max_tract_gap = "rpt_max_tract_gap", return_option = "rpt_return_option",
            count_method = "repeat_count_method", na_handling = "na_repeat_handling",
            flank_pad = "rpt_flank_pad")
  if (is.na(p$path_tracts)) {
    tr <- data.frame(tract = p$rpt_tract_name %||% p$rpt_pattern, motif = p$rpt_pattern,
                     path_ref = p$path_ref, stringsAsFactors = FALSE)
  } else {
    if (!file.exists(p$path_tracts)) stop("Tracts file not found: ", p$path_tracts)
    ext <- tolower(tools::file_ext(p$path_tracts))
    tr <- if (ext == "csv") utils::read.csv(p$path_tracts, stringsAsFactors = FALSE, check.names = FALSE)
          else as.data.frame(readxl::read_excel(p$path_tracts, sheet = "tracts"))
    names(tr) <- tolower(trimws(names(tr)))
    if (!all(c("tract", "motif") %in% names(tr))) stop("Tracts table needs columns 'tract' and 'motif'")
    if (!"path_ref" %in% names(tr)) tr$path_ref <- p$path_ref
    tr$path_ref[is.na(tr$path_ref) | tr$path_ref == ""] <- p$path_ref
  }
  for (k in names(cols)) {
    if (!k %in% names(tr)) tr[[k]] <- NA
    blank <- is.na(tr[[k]]) | tr[[k]] == ""
    tr[[k]][blank] <- p[[cols[[k]]]]
  }
  int_cols <- c("min_repeats", "max_mismatch", "start_perfect_repeats", "end_perfect_repeats",
                "max_gap", "max_tract_gap", "flank_pad")
  for (k in int_cols) tr[[k]] <- as.integer(tr[[k]])
  tr$count_method <- sub("^repeat_count_", "", tr$count_method)
  tr$motif <- toupper(trimws(tr$motif))
  tr$tract <- trimws(as.character(tr$tract))
  if (any(duplicated(tr$tract))) stop("Tract names must be unique")
  if (any(safe_path(tr$tract) != tr$tract)) stop("Tract names may contain only letters, digits, '_' and '-'")
  if (any(is.na(tr$path_ref))) stop("Every tract needs a reference: set path_ref or the tracts table path_ref column")
  tr$path_ref <- normalizePath(tr$path_ref, mustWork = TRUE)
  tr$ref_group <- match(tr$path_ref, unique(tr$path_ref))
  tr$order_in_ref <- stats::ave(seq_len(nrow(tr)), tr$ref_group, FUN = seq_along)
  tr
}

# Composite specification: "name=A+B:sum" (also read from a 'composites' sheet)
resolve_composites <- function(p, tracts) {
  comp <- data.frame(name = character(), tracts = character(), mode = character(), stringsAsFactors = FALSE)
  if (length(p$composites)) {
    for (s in p$composites) {
      m <- regmatches(s, regexec("^\\s*([^=]+)=\\s*([^:]+):\\s*(sum|span)\\s*$", s))[[1]]
      if (length(m) != 4) stop("Cannot parse composite '", s, "'. Use name=TRACT1+TRACT2:sum or :span")
      comp <- rbind(comp, data.frame(name = trimws(m[2]), tracts = gsub("\\s", "", m[3]), mode = m[4]))
    }
  }
  if (!is.na(p$path_tracts) && tolower(tools::file_ext(p$path_tracts)) %in% c("xlsx", "xls") &&
      "composites" %in% readxl::excel_sheets(p$path_tracts)) {
    cs <- as.data.frame(readxl::read_excel(p$path_tracts, sheet = "composites"))
    names(cs) <- tolower(names(cs))
    if (nrow(cs)) comp <- rbind(comp, data.frame(name = cs$name, tracts = gsub("\\s", "", cs$tracts), mode = tolower(cs$mode)))
  }
  if (nrow(comp) == 0) return(comp)
  comp$components <- strsplit(comp$tracts, "+", fixed = TRUE)
  bad <- setdiff(unlist(comp$components), tracts$tract)
  if (length(bad)) stop("Composite refers to unknown tract(s): ", paste(bad, collapse = ", "))
  if (any(lengths(comp$components) < 2)) stop("A composite needs at least two tracts")
  if (any(comp$name %in% tracts$tract) || any(duplicated(comp$name))) stop("Composite names must be unique and differ from tract names")
  comp
}
