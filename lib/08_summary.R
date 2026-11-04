# ==============================================================================
# Duke - run summary
# ==============================================================================
# One workbook at the top of the output folder, run_summary.xlsx (with plain-text
# copies of its tables), that brings together what every module found for each
# sample, so a run can be read, filtered and taken into other analyses without
# opening each module's workbook. It is built from the module workbooks after the
# modules have run, so it is always up to date even when modules were resumed.
# ==============================================================================

write_run_summary <- function(p) {
  suppressPackageStartupMessages(library(data.table))   # attached in a pipeline run already; needed when called on its own
  dir_out <- p$dir_out
  rd <- function(mod, file, sheet) {
    f <- file.path(dir_out, mod, file)
    if (!file.exists(f) || !sheet %in% readxl::excel_sheets(f)) return(NULL)
    x <- tryCatch(as.data.table(readxl::read_excel(f, sheet)), error = function(e) NULL)
    if (is.null(x) || !nrow(x)) NULL else x
  }
  primary <- p$tracts_table$tract[1]
  stage <- rd("01_import_qc", "qc.xlsx", "stage_counts")
  qc <- rd("01_import_qc", "qc.xlsx", "qc_summary")
  acc <- rd("02_alignment", "alignment.xlsx", "read_accounting")
  det <- rd("03_repeat_detection", "repeat_detection.xlsx", "detection")
  umi <- rd("04_umi", "umi.xlsx", "summary")
  al <- rd("05_allele_structure", "allele_structure.xlsx", "alleles")
  im <- rd("06_instability", "instability.xlsx", "instability_metrics")
  so <- rd("07_repeat_lengths", "repeat_lengths.xlsx", "sample_order")
  stems <- unique(c(stage$file_stem, acc$file_stem, det$file_stem))
  if (!length(stems)) return(invisible(NULL))

  s <- data.table(file_stem = stems)
  lab <- short_labels(stems)
  s[, short_label := unname(lab[file_stem])]
  if (!is.null(im)) s <- unique(im[, .(file_stem, group, time, is_control)], by = "file_stem")[s, on = "file_stem"]
  if (!is.null(stage)) s <- stage[, .(file_stem, reads_imported = imported, reads_after_qc = final)][s, on = "file_stem"]
  if (!is.null(qc)) s <- qc[, intersect(c("file_stem", "median_length", "median_np", "median_qv"), names(qc)), with = FALSE][s, on = "file_stem"]
  if (!is.null(qc)) setnames(s, "median_length", "median_read_length_bp", skip_absent = TRUE)
  if (!is.null(acc)) s <- acc[, .(reads_aligned = sum(final)), by = file_stem][s, on = "file_stem"]
  if (!is.null(det)) {
    dw <- dcast(det, file_stem ~ tract, value.var = c("with_tract", "detection_rate"))
    setnames(dw, names(dw), sub("^with_tract_(.*)$", "reads_with_\\1", sub("^detection_rate_(.*)$", "\\1_detection_rate", names(dw))))
    s <- dw[s, on = "file_stem"]
  }
  if (!is.null(umi)) {
    u <- umi[, intersect(c("file_stem", "n_reads_with_umi", "molecules", "kept", "reads_in_kept", "median_reads_per_kept"), names(umi)), with = FALSE]
    if (all(c("reads_in_kept", "n_reads_with_umi") %in% names(u))) u[, share_of_umi_reads_in_kept_molecules := reads_in_kept / n_reads_with_umi]
    setnames(u, c("n_reads_with_umi", "molecules", "kept", "reads_in_kept", "median_reads_per_kept"),
             c("reads_with_umi", "umi_molecules", "molecules_kept", "reads_in_kept_molecules", "median_reads_per_kept_molecule"), skip_absent = TRUE)
    s <- u[s, on = "file_stem"]
  }
  mc <- paste0("modal_", primary)
  if (!is.null(al) && mc %in% names(al)) {
    a12 <- al[allele %in% 1:2]
    vv <- intersect(c(mc, "proportion", "molecule_proportion", "repeat_notation"), names(a12))
    aw <- dcast(a12, file_stem ~ allele, value.var = vv)
    if (length(vv) == 1) setnames(aw, c("1", "2"), paste0(vv, "_", 1:2), skip_absent = TRUE)
    setnames(aw, names(aw), sub(paste0("^", mc, "_(\\d)$"), paste0("allele\\1_", primary),
                                sub("^molecule_proportion_(\\d)$", "allele\\1_molecule_share",
                                    sub("^proportion_(\\d)$", "allele\\1_read_share", sub("^repeat_notation_(\\d)$", "allele\\1_notation", names(aw))))))
    s <- aw[s, on = "file_stem"]
    s <- al[allele > 2, .(other_clusters = .N), by = file_stem][s, on = "file_stem"]
    s[is.na(other_clusters) & file_stem %in% al$file_stem, other_clusters := 0L]
  }
  # Columns in pipeline order, so the table reads left to right as the run went
  ord <- c("file_stem", "short_label", "group", "time", "is_control", "reads_imported", "reads_after_qc",
           "median_read_length_bp", "median_np", "median_qv", "reads_aligned",
           grep("^reads_with_", names(s), value = TRUE), grep("_detection_rate$", names(s), value = TRUE),
           "reads_with_umi", "umi_molecules", "molecules_kept", "reads_in_kept_molecules", "share_of_umi_reads_in_kept_molecules", "median_reads_per_kept_molecule",
           grep("^allele1_", names(s), value = TRUE), grep("^allele2_", names(s), value = TRUE), "other_clusters")
  ord <- intersect(ord, names(s))
  setcolorder(s, c(ord, setdiff(names(s), ord)))
  if (!is.null(so)) { o <- match(s$file_stem, so$file_stem); s <- s[order(o, s$file_stem)] }

  # Headline metrics for the primary measure: one row per sample, unit and range, with
  # the columns most people want; every metric is in 06_instability/instability_metrics.tsv
  hl <- NULL
  if (!is.null(im)) {
    keep <- intersect(c("file_stem", "group", "time", "is_control", "unit", "measure", "analysis_range_name", "setpoint_type", "n_used",
                        unique(c("modal_length", "median_length", "mean_length", "setpoint_length", p$headline_metrics))), names(im))
    # The setpoint type of the headline plots (control_relative by default); every type is
    # in 06_instability/instability_metrics.tsv
    hl <- im[measure == primary, keep, with = FALSE]
    if ("setpoint_type" %in% names(hl) && isTRUE(p$headline_setpoint_type %in% hl$setpoint_type)) hl <- hl[setpoint_type == p$headline_setpoint_type]
    rt <- rd("06_instability", "instability.xlsx", "range_table")
    # each sample's own bounds for each range (a range can be set differently per sample)
    if (!is.null(rt) && "range_text" %in% names(rt))
      hl <- unique(rt[measure == primary, .(file_stem, analysis_range_name = range_name, range = sub(paste0("^", primary, " "), "", range_text))],
                   by = c("file_stem", "analysis_range_name"))[hl, on = .(file_stem, analysis_range_name)]
    hl[, short_label := unname(lab[file_stem])]
    setcolorder(hl, c("file_stem", "short_label", intersect(c("group", "time", "is_control", "unit", "measure", "analysis_range_name", "range", "setpoint_type"), names(hl))))
    setorder(hl, unit, analysis_range_name, group, time, file_stem, na.last = TRUE)
  }

  files <- data.table(
    where = c("00_index.html", "run_summary.xlsx / run_summary_samples.tsv", "run_summary_headline_metrics.tsv",
              "0N_*.html", "01_import_qc/", "02_alignment/", "03_repeat_detection/", "03_repeat_detection/read_repeat_counts/",
              "04_umi/", "04_umi/molecules/", "05_allele_structure/", "06_instability/", "06_instability/instability_metrics.tsv",
              "07_repeat_lengths/", "*/plots/"),
    what = c("Open this first: links every report, table and figure, by sample, range and group",
             "This summary: one row per sample, from every module",
             sprintf("The headline metrics for %s (setpoint type %s), one row per sample, unit (reads or molecules) and range", primary, p$headline_setpoint_type),
             "One report per module: what was done, every figure with how to read it, and the run's numbers",
             "Read counts at each import step, read length and quality (qc.xlsx)",
             "Flank alignment and the reads kept (alignment.xlsx; per-sample alignment figures)",
             "Tract detection and counting methods (repeat_detection.xlsx)",
             "Every read's length for every tract and measure, one file per sample (for your own analyses)",
             "UMI extraction and molecule collapsing (umi.xlsx)",
             "Every molecule's length for every measure, one file per sample",
             "Alleles: clusters, consensus sequences (FASTA), phased variants (VCF), allele BAMs, waterfalls (allele_structure.xlsx)",
             "Instability metrics for every sample, measure, range and unit, with definitions (instability.xlsx)",
             "The same metrics as plain text, for R, Python or Prism",
             "Length distribution figures for every group, range and sample (plots/sample_card/ = one page per sample)",
             "Every figure in a report is saved as a PNG under that module's plots/ folder"))
  cols <- data.table(column = c("file_stem", "short_label", "group, time, is_control", "reads_imported, reads_after_qc",
                                "median_read_length_bp, median_np, median_qv", "reads_aligned",
                                "reads_with_<tract>, <tract>_detection_rate", "reads_with_umi, umi_molecules, molecules_kept",
                                "reads_in_kept_molecules, share_of_umi_reads_in_kept_molecules, median_reads_per_kept_molecule",
                                sprintf("allele1_%s, allele2_%s", primary, primary), "allele1_read_share, allele2_read_share",
                                "allele1_molecule_share, allele2_molecule_share",
                                "allele1_notation, allele2_notation", "other_clusters"),
                     meaning = c("The sample's file name without extension: its identifier in every table, and in the plot file names",
                                 "Short label used in figures, where the full name will not fit: the parts of the file name that differ between this run's samples (always unique within the run, but it depends on which samples are in the run, so use file_stem to match samples between runs)",
                                 "From the settings file", "Reads read in, and left after the import filters (Module 1)",
                                 "Median read length, pass count and read quality (QV; PacBio only)",
                                 "Reads with both flanks aligned in the right order (Module 2)",
                                 "Reads in which each tract was found, and as a share of aligned reads (Module 3)",
                                 "Reads holding a UMI, the UMI molecules (families) found, and those kept (enough reads and not flagged as shortened; Module 4)",
                                 "Reads in kept molecules, as a count and as a share of the reads with a UMI, and the kept molecules' median depth",
                                 sprintf("Modal %s length of the shorter (1) and longer (2) allele (Module 5)", primary),
                                 "Share of the sample's clustered reads in each allele",
                                 "Share of the sample's clustered molecules in each allele (UMI runs clustered on molecules; the better measure of allele balance, as each molecule counts once)",
                                 "Each allele's consensus repeat region in repeat notation",
                                 "Clusters beyond the two alleles (a subclone, PCR chimeras or contamination; see Module 5)"))
  out <- list(samples = s, headline_metrics = hl, columns = cols, files = files)
  xl <- file.path(dir_out, "run_summary.xlsx"); if (file.exists(xl)) file.remove(xl)
  write_excel_function(out, xl)
  fwrite(s, file.path(dir_out, "run_summary_samples.tsv"), sep = "\t", na = "NA")
  if (!is.null(hl)) fwrite(hl, file.path(dir_out, "run_summary_headline_metrics.tsv"), sep = "\t", na = "NA")
  invisible(out)
}

# ------------------------------------------------------------------------------
# 00_index.html: one page at the top of the output folder that links everything
# ------------------------------------------------------------------------------
# The reports, the run summary, then every figure arranged so it can be found
# without knowing the folder layout: cohort figures by module (as thumbnails),
# a table with one row per sample linking each of its figures, one row per range
# linking each range's figures, and one row per group. Figures are found by
# their file names, so new figures are picked up without changing this code.
# Every link is relative, so the page works wherever the folder is copied.
write_run_index <- function(p, summary = NULL) {
  suppressPackageStartupMessages(library(data.table))
  dir_out <- p$dir_out
  esc <- function(x) { x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE); gsub(">", "&gt;", x, fixed = TRUE) }
  href <- function(x) gsub("#", "%23", gsub("%", "%25", x, fixed = TRUE), fixed = TRUE)
  a <- function(path, text) sprintf('<a href="%s">%s</a>', href(path), esc(text))
  pngs <- list.files(dir_out, pattern = "\\.png$", recursive = TRUE)
  pngs <- pngs[!grepl("^(temp|module_data)/|_files/", pngs)]
  s <- summary$samples
  stems <- if (!is.null(s)) s$file_stem else character(0)
  # Every figure is named {type}_{sample or group}_{measure}_{range}_{unit}.png (see
  # plot_catalogue). The parts are read back from the name: the type (longest catalogue
  # match), the unit (last part, when combined/reads/molecules), the sample (the longest
  # file stem the rest starts with, so sample_1 is not taken for sample_10) and the group.
  cat_ <- plot_catalogue()
  fig <- data.table(path = pngs, base = sub("\\.png$", "", basename(pngs)))
  fig[, type := plot_type_of(path)]
  fig[, module := cat_$module[match(type, cat_$type)]]
  fig[is.na(module), module := suppressWarnings(as.integer(substr(path, 2, 2)))]
  fig[, rest := ifelse(is.na(type), base, substring(base, nchar(type) + 2))]
  fig[, unit := ifelse(grepl("_(combined|reads|molecules)$", rest), sub(".*_", "", rest), NA_character_)]
  fig[, rest := ifelse(is.na(unit), rest, sub("_(combined|reads|molecules)$", "", rest))]
  so <- stems[order(-nchar(stems))]
  fig[, sample := vapply(rest, function(r) { m <- so[r == so | startsWith(r, paste0(so, "_"))]; if (length(m)) m[1] else NA_character_ }, "", USE.NAMES = FALSE)]
  groups <- if (!is.null(s) && "group" %in% names(s)) unique(stats::na.omit(s$group)) else character(0)
  gs <- safe_path(groups); gs <- gs[order(-nchar(gs))]
  fig[, group := ifelse(type %in% cat_$type[cat_$scope == "group"],
                        vapply(rest, function(r) { m <- gs[r == gs | startsWith(r, paste0(gs, "_"))]; if (length(m)) m[1] else NA_character_ }, "", USE.NAMES = FALSE), NA_character_)]
  fig[, range := ifelse(grepl("(^|_)range[0-9]+(_|$)", rest), sub(".*(range[0-9]+).*", "\\1", rest), NA_character_)]
  unit_word <- c(combined = "reads + molecules", reads = "reads", molecules = "molecules (UMI-collapsed)")
  uw <- function(u) ifelse(is.na(u), "", unit_word[u])
  # A short label for a figure within a row: what is left of its name after the row's own part
  strip_part <- function(r, part) {
    out <- mapply(function(x, y) if (is.na(y)) x else sub(y, "", x, fixed = TRUE), r, part, USE.NAMES = FALSE)
    gsub("^_+|_+$", "", gsub("_{2,}", "_", out))
  }
  ord_unit <- function(u) match(u, c("combined", "reads", "molecules"), nomatch = 0)

  mod_names <- c(`01_import_qc` = "Module 1: Import and QC", `02_alignment` = "Module 2: Alignment", `03_repeat_detection` = "Module 3: Repeat detection",
                 `04_umi` = "Module 4: UMI", `05_allele_structure` = "Module 5: Allele structure", `06_instability` = "Module 6: Instability",
                 `07_repeat_lengths` = "Module 7: Repeat lengths")
  mod_desc <- c(`01_import_qc` = "reads in, filters, read length and quality", `02_alignment` = "flank alignment and the reads kept",
                `03_repeat_detection` = "tract detection, counting methods, flanks", `04_umi` = "UMI extraction and collapsing reads into molecules",
                `05_allele_structure` = "allele clusters, consensus sequences, variants, waterfalls", `06_instability` = "instability metrics by sample, range and group; timecourses; tract comparison",
                `07_repeat_lengths` = "length distributions: heatmaps, groups, samples, every range")
  h <- c('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">',
         '<title>Duke run index</title><style>',
         'body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:24px auto;max-width:1400px;padding:0 16px;color:#222;background:#fff}',
         'h1{margin-bottom:4px}h2{margin-top:36px;border-bottom:1px solid #ccc;padding-bottom:4px}h3{margin-top:24px}',
         '.meta{color:#666}table{border-collapse:collapse;font-size:13px;margin:8px 0}th,td{border:1px solid #ddd;padding:4px 8px;text-align:left;vertical-align:top}',
         'th{background:#f3f3f3;position:sticky;top:0}.wrap{overflow-x:auto;max-height:80vh}.thumbs{display:flex;flex-wrap:wrap;gap:12px}',
         '.thumbs figure{margin:0;width:240px}.thumbs img{width:240px;height:160px;object-fit:contain;border:1px solid #ddd;background:#fff}',
         '.thumbs figcaption{font-size:12px;word-break:break-all}.note{font-size:13px;color:#555;max-width:900px}a{color:#0b5cad}',
         '.r{color:#555}.m{color:#00838f}code{font-size:12px;word-break:break-all}',
         # Contents bar: stays at the top of the window while scrolling
         'nav.toc{position:sticky;top:0;z-index:10;background:#fff;border-bottom:1px solid #ccc;padding:6px 0;margin:8px 0 0;font-size:14px}nav.toc a{margin-right:16px;text-decoration:none}',
         'h2{scroll-margin-top:44px}',
         # Sample by sample: the sample, group and time columns stay in view when the table is scrolled sideways
         '.fix1{position:sticky;left:0;background:#fff;z-index:1;width:300px;min-width:300px;max-width:300px;word-break:break-all;box-shadow:2px 0 0 #bbb}.fix1 small{color:#666;word-break:normal}',
         'th.fix1{background:#f3f3f3;z-index:3}th{z-index:2}',
         '</style></head><body>',
         '<h1 id="top">Duke run index</h1>',
         sprintf('<p class="meta">Duke %s &middot; %s &middot; %d samples<br>Output folder: <code>%s</code></p>', esc(p$duke_version %||% ""), format(Sys.time(), "%Y-%m-%d %H:%M"),
                 length(stems), esc(dir_out)),
         '<nav class="toc"><a href="#reports">Reports</a><a href="#summary">Run summary</a><a href="#samples">Sample by sample</a><a href="#ranges">Range by range</a><a href="#groups">Group by group</a><a href="#run">Whole-run figures</a><a href="#top">Top</a></nav>',
         '<p class="note">Everything this run produced, linked from one page. Start with the reports and the run summary. Figures are grouped below by what they show: the whole run, one sample, one range, one group. Wherever both are drawn, <span class="r">reads</span> are grey and <span class="m">molecules (UMI-collapsed: one consensus length per UMI)</span> teal, and file and folder names say which.</p>')

  h <- c(h, '<h2 id="reports">Reports</h2><table><tr><th>Report</th><th>What it covers</th><th>Tables</th></tr>')
  for (m in names(mod_names)) {
    rep <- paste0(m, ".html"); if (!file.exists(file.path(dir_out, rep))) next
    # the module's tables, plus the allele FASTA and VCF one folder down in Module 5
    xl <- c(list.files(file.path(dir_out, m), pattern = "\\.(xlsx|tsv)$"),
            grep("^(consensus|variants)/[^/]+\\.(fasta|vcf)$", list.files(file.path(dir_out, m), recursive = TRUE), value = TRUE))
    h <- c(h, sprintf("<tr><td>%s</td><td>%s</td><td>%s</td></tr>", a(rep, mod_names[[m]]), esc(mod_desc[[m]]),
                      paste(vapply(xl, function(x) a(file.path(m, x), x), ""), collapse = "<br>")))
  }
  h <- c(h, "</table>")
  rs <- c("run_summary.xlsx", "run_summary_samples.tsv", "run_summary_headline_metrics.tsv")
  rs <- rs[file.exists(file.path(dir_out, rs))]
  if (length(rs)) h <- c(h, '<h2 id="summary">Run summary</h2><p class="note">One row per sample with what every module found (reads at each step, molecules, alleles), the headline metrics for every range, and a guide to the columns. The .tsv files are the same tables as plain text.</p><p>',
                         paste(vapply(rs, function(x) a(x, x), ""), collapse = " &middot; "), "</p>")

  # Sample by sample: one column per kind of per-sample figure
  sp <- fig[!is.na(sample)]
  if (nrow(sp)) {
    sp[, kind := trimws(paste(type, strip_part(rest, sample), uw(unit)))]
    kinds <- unique(sp[order(module, type, ord_unit(unit), kind), .(module, kind)])
    h <- c(h, '<h2 id="samples">Sample by sample</h2><p class="note">One row per sample, one column per kind of per-sample figure: the plot type, then the measure, then which units it shows (<i>reads + molecules</i>, <i>reads</i>, or <i>molecules</i>). File names follow the same order: <code>{plot type}_{sample}_{measure}_{unit}.png</code>.</p><div class="wrap"><table><tr><th class="fix1">Sample (file name)<br><small>group, time</small></th>',
           paste0("<th>M", kinds$module, " ", esc(kinds$kind), "</th>", collapse = ""), "</tr>")
    for (i in seq_along(stems)) {
      st <- stems[i]
      cells <- vapply(kinds$kind, function(k) { f <- sp[sample == st & kind == k, path]; if (length(f)) a(f[1], "view") else "" }, "")
      # The sample column (with its group and time) stays in view when the table is scrolled sideways
      h <- c(h, sprintf("<tr><td class=\"fix1\" title=\"%s\"><b>%s</b><br><small>%s%s</small></td>%s</tr>", esc(s$short_label[i]), esc(st),
                        esc(display_group(s$group[i] %||% "")), if (!is.na(s$time[i] %||% NA)) paste0(", t = ", esc(as.character(s$time[i]))) else "", paste0("<td>", cells, "</td>", collapse = "")))
    }
    h <- c(h, "</table></div>")
  }

  # Range by range: every figure of a range that is not one sample's
  rp <- fig[!is.na(range) & is.na(sample)]
  if (nrow(rp)) {
    rt <- tryCatch(as.data.table(readxl::read_excel(file.path(dir_out, "06_instability", "instability.xlsx"), "range_table")), error = function(e) NULL)
    # A range can be set differently for different samples in the settings file: list each
    # version, with the number of samples that have it
    rtxt <- if (!is.null(rt) && "range_text" %in% names(rt)) {
      rt[measure == p$tracts_table$tract[1], .(n = uniqueN(file_stem)), by = .(range_name, range_text)][
        , .(range_text = if (.N == 1) range_text else paste0("varies by sample: ", paste(sprintf("%s (%d)", range_text, n), collapse = "; "))), by = range_name]
    } else NULL
    rp[, lab0 := strip_part(rest, range)][, label := trimws(paste(lab0, uw(unit)))]
    h <- c(h, '<h2 id="ranges">Range by range</h2><p class="note">One row per analysis range: every figure drawn for that range (Module 6 timecourses and tract comparison, Module 7 group figures), by plot type. Each sample in every range is under <i>Sample by sample</i> (sample_ranges).</p><div class="wrap"><table><tr><th>Range</th><th>Lengths</th><th>Figures</th></tr>')
    for (r in sort(unique(rp$range))) {
      x <- rp[range == r][order(module, type, lab0, ord_unit(unit))]
      cell <- paste(vapply(unique(x$type), function(t) paste0("<b>", esc(t), "</b>: ", paste(vapply(which(x$type == t), function(k) a(x$path[k], x$label[k]), ""), collapse = " &middot; ")), ""), collapse = "<br>")
      txt <- if (!is.null(rtxt)) rtxt$range_text[match(r, rtxt$range_name)] else ""
      h <- c(h, sprintf("<tr><td>%s</td><td>%s</td><td>%s</td></tr>", r, esc(ifelse(is.na(txt), "", txt)), cell))
    }
    h <- c(h, "</table></div>")
  }

  # Group by group
  gp <- fig[!is.na(group)]
  if (nrow(gp)) {
    gp[, lab0 := strip_part(rest, group)][, label := trimws(paste(lab0, uw(unit)))]
    h <- c(h, '<h2 id="groups">Group by group</h2><p class="note">Module 7: each group&rsquo;s samples stacked in time order, for all lengths and for each range. <i>group_density</i> = ridges with every read (or molecule) as a dot; <i>group_histogram</i> = one row per sample.</p><div class="wrap"><table><tr><th>Group</th><th>Figures</th></tr>')
    for (g in sort(unique(gp$group))) {
      x <- gp[group == g][order(type, lab0, ord_unit(unit))]
      cell <- paste(vapply(unique(x$type), function(t) paste0("<b>", esc(t), "</b>: ", paste(vapply(which(x$type == t), function(k) a(x$path[k], x$label[k]), ""), collapse = " &middot; ")), ""), collapse = "<br>")
      h <- c(h, sprintf("<tr><td>%s</td><td>%s</td></tr>", esc(g), cell))
    }
    h <- c(h, "</table></div>")
  }

  cp <- fig[is.na(sample) & is.na(range) & is.na(group)][order(module, type, ord_unit(unit), base)]
  if (nrow(cp)) {
    h <- c(h, '<h2 id="run">Whole-run figures</h2><p class="note">Figures that summarise every sample at once, by module. Click a thumbnail for the full figure.</p>')
    for (m in names(mod_names)) {
      f <- cp[startsWith(path, paste0(m, "/")), path]; if (!length(f)) next
      h <- c(h, sprintf("<h3>%s</h3><div class=\"thumbs\">", esc(mod_names[[m]])))
      h <- c(h, sprintf('<figure><a href="%s"><img loading="lazy" src="%s" alt=""></a><figcaption>%s</figcaption></figure>', href(f), href(f), esc(basename(f))))
      h <- c(h, "</div>")
    }
  }
  h <- c(h, "</body></html>")
  # 00_ so that it sorts above the module reports; an index.html from before v3.3.3 is replaced
  if (file.exists(file.path(dir_out, "index.html"))) file.remove(file.path(dir_out, "index.html"))
  con <- file(file.path(dir_out, "00_index.html"), open = "w", encoding = "UTF-8"); on.exit(close(con))
  writeLines(enc2utf8(h), con)
  invisible(file.path(dir_out, "00_index.html"))
}
