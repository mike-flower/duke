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
  # Only modules whose results belong to this run (p$current_modules, set by run_duke);
  # a workbook left by an earlier run with other settings is ignored
  cur <- p$current_modules %||% as.character(1:7)
  rd <- function(mod, file, sheet) {
    if (!substr(mod, 2, 2) %in% cur) return(NULL)
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
  # Readable numbers: shares and rates to 4 decimal places, QV to 1
  rnd <- function(d) { for (cn in names(d)) if (is.double(d[[cn]])) data.table::set(d, j = cn, value = round(d[[cn]], if (cn == "median_qv") 1 else 4)); d }
  s <- rnd(s)
  # molecule columns are empty without UMIs
  for (cn in grep("molecule_share$", names(s), value = TRUE)) if (all(is.na(s[[cn]]))) s[, (cn) := NULL]

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
    # reads first, as in every figure
    hl[, unit_order := match(unit, c("reads", "molecules"))]
    setorder(hl, unit_order, analysis_range_name, group, time, file_stem, na.last = TRUE)
    hl[, unit_order := NULL]
    hl <- rnd(hl)
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
                                 "Reads with a UMI and a primary-tract length (the reads Module 4 collapses), the UMI molecules (families) found, and those kept (enough reads and not flagged as shortened; Module 4)",
                                 "Reads in kept molecules, as a count and as a share of the reads with a UMI, and the kept molecules' median depth",
                                 sprintf("Modal %s length of the shorter (1) and longer (2) allele (Module 5)", primary),
                                 "Share of the sample's clustered reads in each allele",
                                 "Share of the sample's clustered molecules in each allele (UMI runs clustered on molecules; the better measure of allele balance, as each molecule counts once)",
                                 "Each allele's consensus repeat region in repeat notation",
                                 "Clusters beyond the two alleles (a subclone, PCR chimeras or contamination; see Module 5)"))
  # List only the columns and outputs this run has (e.g. no UMI rows in a run without UMIs)
  cols <- cols[vapply(strsplit(column, ",\\s*"), function(cc) any(grepl("<", cc)) || any(cc %in% names(s)), TRUE)]
  files <- files[vapply(sub(" .*$", "", where), function(w) grepl("\\*|^run_summary|^00_index", w) || file.exists(file.path(dir_out, w)), TRUE)]
  if (is.null(hl)) files <- files[!grepl("headline", where)]
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
  pngs <- pngs[substr(pngs, 2, 2) %in% (p$current_modules %||% as.character(1:7))]   # this run's modules only
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
  # Only types drawn per unit end in _combined/_reads/_molecules (a file stem can end in
  # e.g. "hifi_reads" too, and must not lose it)
  has_unit <- fig$type %in% cat_$type[cat_$units]
  fig[, unit := ifelse(has_unit & grepl("_(combined|reads|molecules)$", rest), sub(".*_", "", rest), NA_character_)]
  fig[, rest := ifelse(is.na(unit), rest, sub("_(combined|reads|molecules)$", "", rest))]
  # File names carry each stem as safe_path() writes it (spaces etc. replaced), so match on that
  stem_fs <- stats::setNames(stems, safe_path(stems))
  so <- names(stem_fs)[order(-nchar(names(stem_fs)))]
  fig[, sample := vapply(rest, function(r) { m <- so[r == so | startsWith(r, paste0(so, "_"))]; if (length(m)) unname(stem_fs[m[1]]) else NA_character_ }, "", USE.NAMES = FALSE)]
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
  # ---- The page ----------------------------------------------------------------
  # Figures are laid out as grids: what varies (sample, group, measure or range) down the
  # side, and the plot type and unit across the top, so a figure is found by its row and
  # column rather than by reading link text. Each range, group and module folds away, and
  # a search box keeps only the figures whose names contain every word typed.
  desc_of <- function(t) { d <- cat_$description[match(t, cat_$type)]; ifelse(is.na(d), "", d) }
  type_ord <- function(t) { o <- match(t, cat_$type); ifelse(is.na(o), 999L, o) }
  # Without UMIs there are no molecules, and the combined figures show reads only
  umi_run <- isTRUE(p$umi) || any(fig$unit %in% "molecules", na.rm = TRUE)
  unit_lab <- c(combined = if (umi_run) "reads + molecules" else "reads", reads = "reads", molecules = "molecules")
  unit_cls <- function(u) ifelse(is.na(u), "u-n", paste0("u-", substr(u, 1, 1)))
  # Group names as given in the settings file, by the form they take in file names
  grp_name <- stats::setNames(display_group(groups), safe_path(groups))
  gname <- function(g) ifelse(is.na(g), NA_character_, ifelse(g %in% names(grp_name), grp_name[g], g))
  # What the search box matches: the file path with its separators as spaces, the unit in
  # words, the module, and the sample's short label and group as given
  stem_short <- if (!is.null(s) && "short_label" %in% names(s)) stats::setNames(as.character(s$short_label), s$file_stem) else character(0)
  stem_group <- if (!is.null(s) && "group" %in% names(s)) stats::setNames(display_group(s$group), s$file_stem) else character(0)
  fig[, key := tolower(paste(path, gsub("[_/.-]+", " ", sub("\\.png$", "", path)), paste0("m", module),
                             ifelse(is.na(unit), "", unit_lab[unit]),
                             ifelse(is.na(sample), "", paste(stem_short[sample], stem_group[sample])),
                             ifelse(is.na(group), "", gname(group))))]
  fig[is.na(key), key := tolower(path)]
  link <- function(path, unit, key) sprintf('<a class="v %s" href="%s" title="%s" data-k="%s">view</a>',
                                            unit_cls(unit), href(path), esc(basename(path)), esc(key))
  # One grid: x holds the figures, already in row order; rk and rl give each figure's row
  # (key and label); sb splits a plot type's columns (e.g. by measure) above the unit
  grid <- function(x, rk, rl, sb = rep("", nrow(x)), corner = "", sticky = FALSE) {
    x <- copy(x)[, `:=`(rk = rk, rl = rl, sb = sb)]
    x[, ck := paste(type, sb, ifelse(is.na(unit), "", unit), sep = "\r")]
    cols <- unique(x[order(type_ord(type), module, sb, ord_unit(unit)), .(type, module, sb, unit, ck)])
    has_sb <- any(nzchar(cols$sb)); has_u <- any(!is.na(cols$unit))
    nh <- 1L + has_sb + has_u
    runs <- function(v) { r <- rle(v); list(v = r$values, n = r$lengths, first = cumsum(c(1L, head(r$lengths, -1)))) }
    rt <- runs(cols$type)
    hdr <- c('<tr>', sprintf('<th class="rl" rowspan="%d">%s</th>', nh, corner),
             sprintf('<th class="t" colspan="%d" title="%s">%s <small>M%s</small></th>', rt$n, esc(desc_of(rt$v)), esc(rt$v), cols$module[rt$first]), '</tr>')
    if (has_sb) { rs <- runs(paste(cols$type, cols$sb, sep = "\r"))
      hdr <- c(hdr, '<tr>', sprintf('<th class="s" colspan="%d">%s</th>', rs$n, esc(cols$sb[rs$first])), '</tr>') }
    if (has_u) hdr <- c(hdr, '<tr>', sprintf('<th class="u %s">%s</th>', unit_cls(cols$unit), ifelse(is.na(cols$unit), "", unit_lab[cols$unit])), '</tr>')
    rows <- unique(x$rk); by_row <- split(x, factor(x$rk, levels = rows))
    body <- vapply(rows, function(r) {
      y <- by_row[[r]]
      cells <- vapply(cols$ck, function(k) { i <- which(y$ck == k); if (length(i)) link(y$path[i[1]], y$unit[i[1]], y$key[i[1]]) else "" }, "")
      paste0('<tr class="fr"><td class="rl">', y$rl[1], '</td>', paste0('<td>', cells, '</td>', collapse = ""), '</tr>')
    }, "")
    paste0('<div class="gw', if (sticky) ' tall', '"><table class="grid">', paste(hdr, collapse = ""), paste(body, collapse = ""), '</table></div>')
  }
  # Plot types that share the same rows go in one grid; others get their own
  grids_by_rows <- function(x, rk, rl, corner) {
    sets <- vapply(split(rk, x$type), function(v) paste(sort(unique(v)), collapse = "\r"), "")
    fams <- split(names(sets), unname(sets))
    fams <- fams[order(vapply(fams, function(t) min(type_ord(t)), 0))]
    vapply(fams, function(t) { i <- which(x$type %in% t); grid(x[i], rk[i], rl[i], corner = corner) }, "")
  }
  card <- function(id, title, extra, n, body, open = FALSE)
    paste0('<details class="card" id="', id, '"', if (open) ' open data-open', '><summary><span class="ct">', title, '</span>',
           if (nzchar(extra)) paste0('<span class="cx">', extra, '</span>'), '<span class="cn">', n, ' figure', if (n != 1) 's', '</span></summary>',
           '<div class="cb">', paste(body, collapse = ""), '</div></details>')
  fold_buttons <- function(sec) sprintf('<span class="fold"><button data-x="open" data-s="%s">Open all</button><button data-x="close" data-s="%s">Close all</button></span>', sec, sec)

  h <- c('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">',
         '<title>Duke run index</title><style>',
         ':root{--ink:#1f2328;--mute:#5f6b7a;--line:#d8dee4;--soft:#f6f8fa;--link:#0b5cad;--r:#5f6b7a;--rbg:#eef0f3;--m:#00838f;--mbg:#e0f3f4;--c:#3949ab;--cbg:#e8eaf6}',
         '*{box-sizing:border-box}body{font-family:-apple-system,"Segoe UI",Helvetica,Arial,sans-serif;margin:0;color:var(--ink);background:#fff;font-size:14px;line-height:1.45}',
         '.page{max-width:1500px;margin:0 auto;padding:0 20px 60px}',
         'header.top{padding:22px 0 6px}h1{margin:0;font-size:26px}.meta{color:var(--mute);margin:4px 0 0}code{font-size:12px;word-break:break-all}',
         'nav.bar{position:sticky;top:0;z-index:20;background:#fff;border-bottom:1px solid var(--line);display:flex;flex-wrap:wrap;align-items:center;gap:6px 18px;padding:8px 0;margin-bottom:6px}',
         'nav.bar a{color:var(--link);text-decoration:none;font-weight:500}nav.bar a:hover{text-decoration:underline}',
         '.find{margin-left:auto;display:flex;align-items:center;gap:8px}.find input{width:340px;max-width:60vw;padding:6px 10px;border:1px solid #b8c2cc;border-radius:6px;font-size:14px}',
         '.find input:focus{outline:2px solid #9cc3f5;border-color:#6aa3e8}#qn{color:var(--mute);font-size:12px;white-space:nowrap}',
         'h2{font-size:20px;margin:34px 0 4px;padding-bottom:4px;border-bottom:2px solid var(--line);scroll-margin-top:56px;display:flex;align-items:baseline;gap:12px}',
         '.note{color:var(--mute);max-width:980px;margin:6px 0 12px}a{color:var(--link)}',
         '.key{display:flex;gap:8px;flex-wrap:wrap;margin:8px 0 0;align-items:center;color:var(--mute);font-size:12px}',
         '.fold{margin-left:auto;display:flex;gap:6px}.fold button{font-size:12px;padding:3px 10px;border:1px solid var(--line);background:var(--soft);border-radius:5px;cursor:pointer;color:var(--ink)}.fold button:hover{background:#eaeef2}',
         # Reports as cards
         '.reports{display:grid;grid-template-columns:repeat(auto-fill,minmax(310px,1fr));gap:12px}',
         '.rep{border:1px solid var(--line);border-radius:8px;padding:12px 14px;background:#fff}.rep h3{margin:0 0 4px;font-size:15px}.rep p{margin:0 0 8px;color:var(--mute);font-size:13px}',
         '.rep .tabs{display:flex;flex-wrap:wrap;gap:4px 10px;font-size:12px}',
         '.pills{display:flex;gap:8px;flex-wrap:wrap}.pills a{border:1px solid var(--line);border-radius:6px;padding:5px 12px;text-decoration:none;background:var(--soft)}',
         # Fold-away cards
         'details.card{border:1px solid var(--line);border-radius:8px;margin:8px 0;background:#fff}',
         'details.card>summary{cursor:pointer;list-style:none;display:flex;flex-wrap:wrap;gap:4px 14px;align-items:baseline;padding:9px 14px;background:var(--soft);border-radius:8px}',
         'details.card[open]>summary{border-bottom:1px solid var(--line);border-radius:8px 8px 0 0}',
         'details.card>summary::-webkit-details-marker{display:none}details.card>summary:before{content:"\\25B8";color:var(--mute);width:12px}details.card[open]>summary:before{content:"\\25BE"}',
         '.ct{font-weight:600}.cx{color:var(--mute);font-size:13px}.cn{margin-left:auto;color:var(--mute);font-size:12px}.cb{padding:10px 14px 14px}',
         # Grids
         '.gw{overflow:auto;margin:6px 0 14px;border:1px solid var(--line);border-radius:6px;max-width:100%}.gw.tall{max-height:78vh}',
         'table.grid{border-collapse:separate;border-spacing:0;font-size:13px}',
         '.grid th,.grid td{border-bottom:1px solid #eaeef2;border-right:1px solid #eaeef2;padding:4px 8px;text-align:center;vertical-align:middle;white-space:nowrap}',
         '.grid th{background:var(--soft);position:sticky;top:0;z-index:2;font-weight:600}.grid th.t small{color:var(--mute);font-weight:400;margin-left:4px}',
         '.grid th.s{font-weight:500;color:var(--ink)}.grid th.u{font-weight:500;font-size:12px}',
         '.grid tr:nth-child(2) th{top:29px}.grid tr:nth-child(3) th{top:57px}',
         '.grid td.rl,.grid th.rl{position:sticky;left:0;text-align:left;background:#fff;z-index:1;white-space:normal;min-width:220px;max-width:420px;box-shadow:1px 0 0 var(--line)}',
         '.grid th.rl{background:var(--soft);z-index:3;color:var(--mute);font-weight:500}',
         '.grid tr.fr:hover td{background:#f9fbfd}.grid tr.fr:hover td.rl{background:#f3f7fb}',
         '.rl .sub{color:var(--mute);font-size:12px;margin-left:6px}.rl small{color:var(--mute)}.rl .len{color:var(--mute);font-size:12px;margin-left:6px}',
         # Unit colours: column headings and links
         '.u-c{color:var(--c)}.u-r{color:var(--r)}.u-m{color:var(--m)}',
         'a.v{display:inline-block;padding:1px 9px;border-radius:10px;font-size:12px;text-decoration:none;border:1px solid transparent}',
         'a.v.u-c{background:var(--cbg)}a.v.u-r{background:var(--rbg)}a.v.u-m{background:var(--mbg)}a.v.u-n{background:#eef3fb;color:var(--link)}a.v:hover{border-color:currentColor}',
         '.chip{display:inline-block;padding:1px 9px;border-radius:10px;font-size:12px}.chip.u-c{background:var(--cbg)}.chip.u-r{background:var(--rbg)}.chip.u-m{background:var(--mbg)}',
         # Summary-figure thumbnails
         '.thumbs{display:grid;grid-template-columns:repeat(auto-fill,minmax(230px,1fr));gap:14px}',
         '.thumbs figure{margin:0}.thumbs img{width:100%;height:160px;object-fit:contain;border:1px solid var(--line);border-radius:6px;background:#fff}',
         '.thumbs figcaption{font-size:12px;word-break:break-all;margin-top:4px}.thumbs figcaption span{display:block;color:var(--mute);word-break:normal}',
         # Search: what does not match is hidden (kept in place in grids, so columns stay aligned)
         'a.v.off{visibility:hidden}figure.off{display:none}section.fs.empty{display:none}.none{color:var(--mute);font-style:italic}',
         '</style></head><body><div class="page">',
         '<header class="top" id="top"><h1>Duke run index</h1>',
         sprintf('<p class="meta">Duke %s &middot; %s &middot; %d samples &middot; output folder <a href="./" title="Open the output folder"><code>%s</code></a></p></header>',
                 esc(p$duke_version %||% ""), format(Sys.time(), "%Y-%m-%d %H:%M"), length(stems), esc(dir_out)),
         '<nav class="bar"><a href="#reports">Reports</a><a href="#summary">Run summary</a><a href="#samples">Samples</a><a href="#ranges">Ranges</a><a href="#groups">Groups</a><a href="#run">Summary figures</a><a href="#top">Top</a>',
         '<span class="find"><input id="q" type="search" placeholder="Find figures, e.g. range03 molecules" autocomplete="off"><span id="qn"></span></span></nav>',
         '<p class="note">Everything this run produced, linked from one page. Start with the reports and the run summary. Figures are arranged below by what they show: one sample (<i>Samples</i>), one range (<i>Ranges</i>), one group (<i>Groups</i>), or every sample at once (<i>Summary figures</i>). In each grid, find the row (the sample, group, measure or range) and the column (the plot type, then the unit). Type in the search box to keep only the figures whose file names contain every word you type (a sample, group, range, measure, plot type or unit).</p>',
         if (umi_run) '<div class="key">Units: <span class="chip u-c">reads + molecules</span><span class="chip u-r">reads</span><span class="chip u-m">molecules (UMI-collapsed: one consensus length per UMI)</span></div>')

  # Reports, one card per module
  h <- c(h, '<h2 id="reports">Reports</h2><div class="reports">')
  for (m in names(mod_names)) {
    rep <- paste0(m, ".html"); if (!file.exists(file.path(dir_out, rep)) || !substr(m, 2, 2) %in% (p$current_modules %||% as.character(1:7))) next
    # the module's tables, plus the allele FASTA and VCF one folder down in Module 5
    xl <- c(list.files(file.path(dir_out, m), pattern = "\\.(xlsx|tsv)$"),
            grep("^(consensus|variants)/[^/]+\\.(fasta|vcf)$", list.files(file.path(dir_out, m), recursive = TRUE), value = TRUE))
    h <- c(h, sprintf('<div class="rep"><h3>%s</h3><p>%s</p><div class="tabs">%s</div></div>', a(rep, mod_names[[m]]), esc(mod_desc[[m]]),
                      paste(vapply(xl, function(x) a(file.path(m, x), x), ""), collapse = "")))
  }
  h <- c(h, "</div>")
  rs <- c("run_summary.xlsx", "run_summary_samples.tsv", "run_summary_headline_metrics.tsv")
  rs <- rs[file.exists(file.path(dir_out, rs))]
  if (length(rs)) h <- c(h, '<h2 id="summary">Run summary</h2><p class="note">One row per sample with what every module found (reads at each step, molecules, alleles), the headline metrics for every range, and a guide to the columns. The .tsv files are the same tables as plain text.</p><div class="pills">',
                         paste(vapply(rs, function(x) a(x, x), ""), collapse = ""), "</div>")

  # Samples: one row per sample; columns by plot type, then measure, then unit
  sp <- fig[!is.na(sample)]
  if (nrow(sp)) {
    sp[, sb := strip_part(rest, safe_path(sample))]
    sp[, so := match(sample, stems)]
    sp <- sp[order(so, type_ord(type), sb, ord_unit(unit))]
    rl <- vapply(seq_len(nrow(sp)), function(i) {
      j <- sp$so[i]
      sprintf('<b title="%s">%s</b><br><small>%s%s</small>', esc(s$short_label[j] %||% ""), esc(sp$sample[i]), esc(display_group(s$group[j] %||% "")),
              if (!is.na(s$time[j] %||% NA)) paste0(", t = ", esc(as.character(s$time[j]))) else "")
    }, "")
    h <- c(h, '<section class="fs" id="sec-samples"><h2 id="samples">Samples</h2><p class="note">One row per sample (with its group and time); one column per per-sample figure, headed by the plot type (hover for what it shows), then the measure, then the unit. File names follow the same order: <code>{plot type}_{sample}_{measure}_{unit}.png</code>.</p>',
           grid(sp, sp$sample, rl, sp$sb, corner = "Sample <small>(group, time)</small>", sticky = TRUE), '</section>')
  }

  # Ranges: one card per range; inside, the range's figures by group or measure
  rp <- fig[!is.na(range) & is.na(sample)]
  rtxt <- NULL; rt <- NULL
  if (nrow(rp) || any(!is.na(fig$range))) {
    rt <- tryCatch(as.data.table(readxl::read_excel(file.path(dir_out, "06_instability", "instability.xlsx"), "range_table")), error = function(e) NULL)
    # A range can be set differently for different samples in the settings file: list each
    # version, with the number of samples that have it
    rtxt <- if (!is.null(rt) && "range_text" %in% names(rt)) {
      rt[measure == p$tracts_table$tract[1], .(n = uniqueN(file_stem)), by = .(range_name, range_text)][
        , .(range_text = if (.N == 1) range_text else paste0("varies by sample: ", paste(sprintf("%s (%d)", range_text, n), collapse = "; "))), by = range_name]
    } else NULL
  }
  range_len <- function(r) { t <- if (!is.null(rtxt)) rtxt$range_text[match(r, rtxt$range_name)] else NA; ifelse(is.na(t), "", t) }
  # The same for the samples of one group, as a group usually shares one setting
  range_len_in <- function(r, st) {
    if (is.null(rt) || !"range_text" %in% names(rt) || !length(st)) return(range_len(r))
    x <- rt[measure == p$tracts_table$tract[1] & file_stem %in% st, .(n = uniqueN(file_stem)), by = .(range_name, range_text)]
    vapply(r, function(r1) { y <- x[range_name == r1]
      if (!nrow(y)) range_len(r1) else if (nrow(y) == 1) y$range_text else paste0("varies by sample: ", paste(sprintf("%s (%d)", y$range_text, y$n), collapse = "; ")) }, "")
  }
  if (nrow(rp)) {
    rp[, rem := strip_part(strip_part(rest, range), group)]
    rp[, gn := gname(group)]
    rp <- rp[order(range, type_ord(type), is.na(gn), gn, rem, ord_unit(unit))]
    h <- c(h, '<section class="fs" id="sec-ranges"><h2 id="ranges">Ranges', fold_buttons("sec-ranges"), '</h2><p class="note">One card per analysis range (click to open): the range&rsquo;s timecourses and tract comparisons (Module 6) by measure, and its group figures (Module 7) by group. Each sample in every range is under <i>Samples</i> (sample_ranges).</p>')
    for (r in unique(rp$range)) {
      x <- rp[range == r]
      rk <- paste(ifelse(is.na(x$gn), "", x$gn), x$rem, sep = "\r")
      rl <- ifelse(is.na(x$gn), sprintf("<b>%s</b>", esc(x$rem)),
                   paste0("<b>", esc(x$gn), "</b>", ifelse(nzchar(x$rem), sprintf('<span class="sub">%s</span>', esc(x$rem)), "")))
      body <- c(grids_by_rows(x[is.na(gn)], rk[is.na(x$gn)], rl[is.na(x$gn)], "Measure"),
                grids_by_rows(x[!is.na(gn)], rk[!is.na(x$gn)], rl[!is.na(x$gn)], "Group"))
      h <- c(h, card(paste0("range-", r), esc(r), esc(range_len(r)), nrow(x), body))
    }
    h <- c(h, '</section>')
  }

  # Groups: one card per group; inside, all lengths and each range
  gp <- fig[!is.na(group)]
  if (nrow(gp)) {
    gp[, rem := strip_part(rest, group)]
    gp[, span := ifelse(!is.na(range), range, ifelse(grepl("(^|_)all-lengths(_|$)", rem), "all lengths", ""))]
    gp[, meas := strip_part(strip_part(rem, range), "all-lengths")]
    gp <- gp[order(group, type_ord(type), meas, span != "all lengths", span, ord_unit(unit))]
    h <- c(h, '<section class="fs" id="sec-groups"><h2 id="groups">Groups', fold_buttons("sec-groups"), '</h2><p class="note">One card per group (click to open): Module 7&rsquo;s figures of the group&rsquo;s samples stacked in time order, for all lengths and for each range. <i>group_density</i> = ridges with every read (or molecule) as a dot; <i>group_histogram</i> = one row per sample.</p>')
    for (g in unique(gp$group)) {
      x <- gp[group == g]
      g_stems <- if (!is.null(s) && "group" %in% names(s)) s$file_stem[safe_path(s$group) == g] else character(0)
      rk <- paste(x$meas, x$span, sep = "\r")
      rl <- paste0(ifelse(nzchar(x$span), sprintf("<b>%s</b>", esc(x$span)), ""),
                   ifelse(!is.na(x$range), sprintf('<span class="len">%s</span>', esc(range_len_in(x$range, g_stems))), ""),
                   ifelse(nzchar(x$meas), sprintf(ifelse(nzchar(x$span), '<span class="sub">%s</span>', '<b>%s</b>'), esc(x$meas)), ""))
      gs_ <- if (!is.null(s) && "group" %in% names(s)) s[display_group(group) == gname(g)] else NULL
      extra <- if (!is.null(gs_) && nrow(gs_)) paste0(nrow(gs_), " sample", if (nrow(gs_) != 1) "s",
                                                      if ("time" %in% names(gs_) && any(!is.na(gs_$time))) paste0(", t = ", paste(sort(unique(stats::na.omit(gs_$time))), collapse = ", ")) else "") else ""
      sp_ <- nzchar(x$span)   # figures of all lengths or a range, then those of neither (e.g. a tract pair)
      h <- c(h, card(paste0("group-", g), esc(gname(g)), esc(extra), nrow(x),
                     c(grids_by_rows(x[sp_], rk[sp_], rl[sp_], "Lengths"), grids_by_rows(x[!sp_], rk[!sp_], rl[!sp_], "Measure"))))
    }
    h <- c(h, '</section>')
  }

  cp <- fig[is.na(sample) & is.na(range) & is.na(group)][order(module, type_ord(type), ord_unit(unit), base)]
  if (nrow(cp)) {
    h <- c(h, paste0('<section class="fs" id="sec-run"><h2 id="run">Summary figures', fold_buttons("sec-run"), '</h2><p class="note">Figures that show every sample in the run together, rather than one sample, range or group: for example read counts and quality per sample (Module 1)', if (dir.exists(file.path(dir_out, "04_umi"))) ", the length distributions of reads against molecules (Module 4)", if (dir.exists(file.path(dir_out, "05_allele_structure"))) ", the clusters in every sample (Module 5)", if (dir.exists(file.path(dir_out, "06_instability"))) ", the instability index in every range (Module 6)", '. By module; click a thumbnail for the full figure.</p>'))
    for (m in names(mod_names)) {
      x <- cp[startsWith(path, paste0(m, "/"))]; if (!nrow(x)) next
      th <- sprintf('<figure data-k="%s"><a href="%s"><img loading="lazy" src="%s" alt=""></a><figcaption>%s<span>%s</span></figcaption></figure>',
                    esc(x$key), href(x$path), href(x$path), esc(basename(x$path)), esc(desc_of(x$type)))
      h <- c(h, card(paste0("run-", m), esc(mod_names[[m]]), "", nrow(x), c('<div class="thumbs">', th, '</div>'), open = TRUE))
    }
    h <- c(h, '</section>')
  }
  h <- c(h, '<p class="none" id="nomatch" hidden>No figure&rsquo;s file name contains all of those words.</p></div>',
         '<script>(function(){',
         'var q=document.getElementById("q"),qn=document.getElementById("qn"),none=document.getElementById("nomatch");',
         'var figs=[].slice.call(document.querySelectorAll("[data-k]"));',
         'function src(e){var a=e.tagName==="A"?e:e.querySelector("a");return a?a.getAttribute("href"):"";}',
         'function uniq(l){var o={};l.forEach(function(e){o[src(e)]=1;});return Object.keys(o).length;}',
         'var total=uniq(figs);',
         'function each(sel,f){[].forEach.call(document.querySelectorAll(sel),f);}',
         'function live(e){return e.querySelector("[data-k]:not(.off)");}',
         'function run(){var t=q.value.toLowerCase().trim().split(/\\s+/).filter(Boolean),c=0;',
         ' var hit=[];figs.forEach(function(e){var k=e.getAttribute("data-k"),ok=t.every(function(w){return k.indexOf(w)>=0;});e.classList.toggle("off",!ok);if(ok)hit.push(e);});c=uniq(hit);',
         ' var on=t.length>0;',
         ' each("tr.fr",function(r){r.hidden=on&&!live(r);});',
         ' each(".gw",function(g){g.hidden=on&&!live(g);});',
         ' each("details.card",function(d){var any=live(d);d.hidden=on&&!any;d.open=on?!!any:d.hasAttribute("data-open");});',
         ' each("section.fs",function(s){s.classList.toggle("empty",on&&!live(s));});',
         ' qn.textContent=on?(c+" of "+total+" figures"):(total+" figures");none.hidden=!(on&&c===0);}',
         'q.addEventListener("input",run);',
         'each("button[data-x]",function(b){b.addEventListener("click",function(){var o=b.getAttribute("data-x")==="open";each("#"+b.getAttribute("data-s")+" details.card",function(d){if(!d.hidden)d.open=o;});});});',
         'run();})();</script>')
  h <- c(h, "</body></html>")
  # 00_ so that it sorts above the module reports; an index.html from before v3.3.3 is replaced
  if (file.exists(file.path(dir_out, "index.html"))) file.remove(file.path(dir_out, "index.html"))
  con <- file(file.path(dir_out, "00_index.html"), open = "w", encoding = "UTF-8"); on.exit(close(con))
  writeLines(enc2utf8(h), con)
  invisible(file.path(dir_out, "00_index.html"))
}
