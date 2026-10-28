# ==============================================================================
# Duke - repeat tract detection and counting
# ==============================================================================
# For each read, each tract is searched for in the region between the two flank
# alignments, extended by `flank_pad` bases into each flank. The extension lets
# a tract run on past the point where the aligner happened to place the flank
# boundary (for example, an LOI allele whose last pure CAG minimap2 assigns to
# the right flank). A tract must overlap the between-flanks region, so repeats
# that lie wholly inside a flank are never picked up.
#
# Tracts sharing a reference are found in table order (5' to 3'): each is
# searched for downstream of the previous one. The sequence between consecutive
# tracts is kept per read.
#
# Finding a tract (all reads at once, vectorised):
#   1. every occurrence of the motif (allowing max_mismatch per unit)
#   2. occurrences no more than max_gap bases apart form a run; runs of at least
#      min_repeats units are kept
#   3. runs no more than max_tract_gap apart are merged into tracts
#   4. each end is trimmed back to the first / last point where
#      start_perfect_repeats / end_perfect_repeats exact units occur in a row
#   5. the longest (or first) tract overlapping the between-flanks region is kept
#
# Three counts per tract:
#   full  = tract length in bases / motif length (rounded). Includes bases of any
#           interruption inside the tract.
#   match = number of exact motif units inside the tract
#   pure  = longest run of exact, uninterrupted units
# ==============================================================================

.self_overlaps <- function(motif) {
  L <- nchar(motif)
  if (L < 2) return(TRUE)
  any(vapply(1:(L - 1), function(k) substr(motif, 1, k) == substr(motif, L - k + 1, L), logical(1)))
}

# Long table of motif matches: read index (i) and start/end in window coordinates
.motif_hits <- function(windows, motif, max_mismatch) {
  dss <- Biostrings::DNAStringSet(windows)
  m <- Biostrings::vmatchPattern(motif, dss, max.mismatch = max_mismatch, fixed = TRUE)
  st <- Biostrings::startIndex(m)
  n <- lengths(st)
  if (sum(n) == 0) return(data.table::data.table(i = integer(), s = integer(), e = integer()))
  L <- nchar(motif)
  h <- data.table::data.table(i = rep.int(seq_along(st), n), s = unlist(st, use.names = FALSE))
  h <- h[s >= 1]
  h[, e := s + L - 1L]
  h <- h[e <= nchar(windows)[i]]
  data.table::setorder(h, i, s)
  if (.self_overlaps(motif) || max_mismatch > 0) {
    # Greedy non-overlapping selection (needed for self-overlapping motifs and
    # for inexact matching, where e.g. CAG and GAG hits can overlap in CAGAG)
    keep <- logical(nrow(h)); last_i <- -1L; last_e <- -1L
    for (k in seq_len(nrow(h))) {
      if (h$i[k] != last_i || h$s[k] > last_e) { keep[k] <- TRUE; last_i <- h$i[k]; last_e <- h$e[k] }
    }
    h <- h[keep]
  }
  h
}

# Find one tract per window. Returns data.table(i, start, end) in window coords
# for windows where a tract overlapping [core_lo, core_hi] was found.
find_tracts_vectorised <- function(windows, core_lo, core_hi, motif, min_repeats = 2L,
                                   max_mismatch = 0L, start_perfect = 2L, end_perfect = 2L,
                                   max_gap = 6L, max_tract_gap = 18L, return_option = "longest") {
  L <- nchar(motif)
  empty <- data.table::data.table(i = integer(), start = integer(), end = integer())
  hits <- .motif_hits(windows, motif, max_mismatch)
  if (nrow(hits) == 0) return(empty)

  # 2. runs
  hits[, gap := s - data.table::shift(e) - 1L, by = i]
  hits[, run := cumsum(is.na(gap) | gap > max_gap), by = i]
  runs <- hits[, .(s = min(s), e = max(e), n = .N), by = .(i, run)][n >= min_repeats]
  # A run may extend into the flank padding, but a run lying wholly inside the
  # padding belongs to the flank and is never merged into the tract
  runs <- runs[s <= core_hi[i] + 1L & e >= core_lo[i] - 1L]
  if (nrow(runs) == 0) return(empty)

  # 3. merge runs into tracts
  data.table::setorder(runs, i, s)
  runs[, gap := s - data.table::shift(e) - 1L, by = i]
  runs[, tr := cumsum(is.na(gap) | gap > max_tract_gap), by = i]
  tracts <- runs[, .(start = min(s), end = max(e)), by = .(i, tr)]

  # 4. perfect-unit trimming from exact matches
  exact <- if (max_mismatch == 0) hits[, .(i, s)] else .motif_hits(windows, motif, 0L)[, .(i, s)]
  if (nrow(exact) == 0 && (start_perfect > 0 || end_perfect > 0)) return(empty)
  if (nrow(exact) > 0 && (start_perfect > 0 || end_perfect > 0)) {
    data.table::setorder(exact, i, s)
    exact[, chain := cumsum(is.na(data.table::shift(s)) | (s - data.table::shift(s)) != L), by = i]
    exact[, `:=`(pos = seq_len(.N), size = .N), by = .(i, chain)]
    exact[, `:=`(fwd = size - pos + 1L, bwd = pos)]
    tracts[, tid := .I]
    # candidate starts: exact unit p inside the tract with >= start_perfect units ahead
    if (start_perfect > 0) {
      cs <- exact[fwd >= start_perfect, .(i, p = s)]
      j <- cs[tracts, on = .(i), allow.cartesian = TRUE, nomatch = NULL][
        p >= start & p + start_perfect * L - 1L <= end, .(new_start = min(p)), by = tid]
      tracts[j, start := i.new_start, on = .(tid)]
      tracts <- tracts[tid %in% j$tid]   # no qualifying perfect start: not a tract
    }
    if (end_perfect > 0) {
      ce <- exact[bwd >= end_perfect, .(i, q = s + L - 1L)]
      j <- ce[tracts, on = .(i), allow.cartesian = TRUE, nomatch = NULL][
        q <= end & q - end_perfect * L + 1L >= start, .(new_end = max(q)), by = tid]
      tracts[j, end := i.new_end, on = .(tid)]
      tracts <- tracts[tid %in% j$tid]
    }
    tracts <- tracts[end >= start]
  }

  # 5. overlap with the between-flanks region, then choose one tract per read
  tracts[, `:=`(clo = core_lo[i], chi = core_hi[i])]
  tracts <- tracts[start <= chi + 1L & end >= clo - 1L]
  if (nrow(tracts) == 0) return(empty)
  if (return_option == "first") {
    data.table::setorder(tracts, i, start)
  } else {
    tracts[, len := end - start + 1L]
    data.table::setorder(tracts, i, -len, start)
  }
  tracts[!duplicated(i), .(i, start, end)]
}

# Counts for a vector of tract sequences
count_tract_units <- function(seqs, motif) {
  L <- nchar(motif)
  n <- length(seqs)
  out <- data.table::data.table(full = rep(NA_real_, n), match = NA_integer_, pure = NA_integer_,
                                structure = NA_character_)
  ok <- !is.na(seqs) & nchar(seqs) > 0
  if (!any(ok)) return(out)
  out$full[ok] <- floor(nchar(seqs[ok]) / L + 0.5)
  h <- .motif_hits(seqs[ok], motif, 0L)
  if (nrow(h) == 0) { out$match[ok] <- 0L; out$pure[ok] <- 0L; out$structure[ok] <- "0"; return(out) }
  # contiguous chains of exact units (any gap breaks the chain)
  h[, chain := cumsum(is.na(data.table::shift(e)) | s != data.table::shift(e) + 1L), by = i]
  ch <- h[, .(n = .N), by = .(i, chain)]
  agg <- ch[, .(match = sum(n), pure = max(n), structure = paste(n, collapse = ",")), by = i]
  idx <- which(ok)
  out$match[idx] <- 0L; out$pure[idx] <- 0L; out$structure[idx] <- "0"
  out$match[idx[agg$i]] <- agg$match
  out$pure[idx[agg$i]] <- agg$pure
  out$structure[idx[agg$i]] <- agg$structure
  out
}

# ------------------------------------------------------------------------------
# All tracts for one sample and one reference group
# ------------------------------------------------------------------------------
# aln: Module 2 table for the sample (one row per read)
# tr:  tracts rows for this reference, in 5'->3' order
process_sample_repeats <- function(aln, tr) {
  q <- aln$query_sequence
  qlen <- nchar(q)
  core_lo_read <- aln$left_qend + 1L
  core_hi_read <- pmax(aln$right_qstart - 1L, aln$left_qend)
  out <- data.table::data.table(file_stem = aln$file_stem, qname = aln$qname, ref_id = aln$ref_id,
                                strand = aln$strand)
  prev_end <- rep(NA_integer_, nrow(aln))   # end of the last tract found so far (search window)
  prev_name <- NULL; prev_en <- NULL        # the immediately preceding tract
  for (k in seq_len(nrow(tr))) {
    t <- tr[k, ]
    pad <- t$flank_pad
    win_lo <- pmax(1L, aln$left_qend + 1L - pad)
    win_hi <- pmin(qlen, aln$right_qstart - 1L + pad)
    # later tracts are searched downstream of the previous tract
    has_prev <- !is.na(prev_end)
    win_lo[has_prev] <- pmax(win_lo[has_prev], prev_end[has_prev] + 1L)
    bad <- win_hi < win_lo
    win_hi[bad] <- win_lo[bad]
    windows <- substring(q, win_lo, win_hi)
    found <- find_tracts_vectorised(
      windows, core_lo = core_lo_read - win_lo + 1L, core_hi = core_hi_read - win_lo + 1L,
      motif = t$motif, min_repeats = t$min_repeats, max_mismatch = t$max_mismatch,
      start_perfect = t$start_perfect_repeats, end_perfect = t$end_perfect_repeats,
      max_gap = t$max_gap, max_tract_gap = t$max_tract_gap, return_option = t$return_option)
    st <- rep(NA_integer_, nrow(aln)); en <- st
    st[found$i] <- found$start + win_lo[found$i] - 1L
    en[found$i] <- found$end + win_lo[found$i] - 1L
    seqs <- ifelse(is.na(st), NA_character_, substring(q, st, en))
    cnt <- count_tract_units(seqs, t$motif)
    nm <- t$tract
    out[, paste0(nm, "_start") := st]
    out[, paste0(nm, "_end") := en]
    out[, paste0(nm, "_full") := cnt$full]
    out[, paste0(nm, "_match") := cnt$match]
    out[, paste0(nm, "_pure") := cnt$pure]
    out[, paste0(nm, "_structure") := cnt$structure]
    out[, paste0(nm, "_seq") := seqs]
    if (!is.null(prev_name)) {
      # only when both this tract and the one immediately before it were found
      bt <- ifelse(!is.na(prev_en) & !is.na(st), substring(q, prev_en + 1L, st - 1L), NA_character_)
      out[, paste0("between_", prev_name, "_", nm) := bt]
    }
    prev_end <- ifelse(is.na(en), prev_end, en)
    prev_name <- nm; prev_en <- en
  }
  out
}

# Primary-tract segments in the layout older Duke modules used (left / mid / right),
# so existing downstream scripts that read alignment_counts_slim keep working.
primary_segments <- function(aln, rt, primary) {
  q <- aln$query_sequence
  st <- rt[[paste0(primary, "_start")]]; en <- rt[[paste0(primary, "_end")]]
  found <- !is.na(st)
  left <- ifelse(found, substring(q, aln$left_qstart, st - 1L), substring(q, aln$left_qstart, aln$left_qend))
  mid <- ifelse(found, substring(q, st, en), NA_character_)
  right <- ifelse(found, substring(q, en + 1L, aln$right_qend),
                  substring(q, aln$left_qend + 1L, aln$right_qend))
  data.table::data.table(left = left, mid = mid, right = right)
}

# ------------------------------------------------------------------------------
# Composite measures from the per-read tract table
# ------------------------------------------------------------------------------
compute_composites <- function(rt, composites, tracts) {
  if (is.null(composites) || nrow(composites) == 0) return(rt)
  for (k in seq_len(nrow(composites))) {
    comps <- composites$components[[k]]
    nm <- composites$name[k]
    if (composites$mode[k] == "sum") {
      vals <- lapply(comps, function(t) {
        cm <- tracts$count_method[tracts$tract == t]
        rt[[paste0(t, "_", cm)]]
      })
      rt[, (nm) := Reduce(`+`, vals)]
    } else {
      starts <- do.call(pmin, lapply(comps, function(t) rt[[paste0(t, "_start")]]))
      ends <- do.call(pmax, lapply(comps, function(t) rt[[paste0(t, "_end")]]))
      bp <- ends - starts + 1L
      L <- unique(nchar(tracts$motif[tracts$tract %in% comps]))
      rt[, paste0(nm, "_bp") := bp]
      rt[, (nm) := if (length(L) == 1) floor(bp / L + 0.5) else bp]
    }
  }
  rt
}
