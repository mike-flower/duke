# ==============================================================================
# Duke - UMI extraction and collapse to molecules
# ==============================================================================
# A UMI family descends from one original DNA molecule, which has one true
# repeat length. The spread of repeat calls within a family is therefore PCR
# slippage plus sequencing error, and collapsing a family to one value gives
# one measurement per molecule.
#
#   1. The UMI is the sequence between two user-set flanks (default M13-F and
#      ciosiF+1), searched on the strand-corrected read. An exact flank match is
#      used where there is one; otherwise an indel-tolerant match.
#   2. UMIs are clustered per sample with the umi_tools "directional" rule:
#      UMI b joins UMI a when they are within umi_cluster_max_dist edits and
#      count(a) >= 2 * count(b) - 1, followed transitively from the most
#      abundant UMI. Neighbours are found with a deletion-key index, so the
#      cost grows with the number of UMIs rather than its square.
#   3. Each UMI family is one molecule. Repeat length plays no part in forming
#      molecules, so nothing is assumed about how many alleles a sample has.
#   4. Molecules with fewer than umi_min_reads reads are dropped, and each
#      measure is collapsed to one value per molecule (default top_cluster,
#      needing at least 2 reads and 20% of the molecule's reads).
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Extraction
# ------------------------------------------------------------------------------
.flank_locate <- function(x, pat, mm, want = c("end", "start")) {
  want <- match.arg(want)
  n <- length(x)
  st <- rep(NA_integer_, n); len <- rep(NA_integer_, n)
  ex <- regexpr(pat, x, fixed = TRUE)
  hit <- ex > 0
  st[hit] <- ex[hit]; len[hit] <- nchar(pat)
  todo <- which(!hit & nchar(x) >= nchar(pat) - mm)
  if (length(todo) && mm > 0) {
    m <- aregexec(pat, x[todo], max.distance = list(all = mm), fixed = TRUE)
    s2 <- vapply(m, function(z) as.integer(z[1]), integer(1))
    l2 <- vapply(m, function(z) as.integer(attr(z, "match.length")[1]), integer(1))
    ok <- s2 > 0
    st[todo[ok]] <- s2[ok]; len[todo[ok]] <- l2[ok]
  }
  if (want == "end") st + len - 1L else st
}

extract_umis <- function(seqs, flank5, flank3, mm = 3L, window = 250L, len_min = 12L, len_max = 16L) {
  seqs[is.na(seqs)] <- ""
  n <- length(seqs)
  head5 <- substr(seqs, 1, window + nchar(flank5))
  e5 <- .flank_locate(head5, flank5, mm, "end")
  found5 <- !is.na(e5)
  # The 3' flank is searched only just downstream of the 5' flank
  reach <- len_max + nchar(flank3) + mm + 5L
  s3 <- rep(NA_integer_, n)
  if (any(found5)) {
    seg <- substring(seqs[found5], e5[found5] + 1L, e5[found5] + reach)
    s3_local <- .flank_locate(seg, flank3, mm, "start")
    s3[found5] <- ifelse(is.na(s3_local), NA_integer_, e5[found5] + s3_local)
  }
  found3 <- !is.na(s3)
  umi <- rep(NA_character_, n)
  ok <- found5 & found3 & s3 > e5 + 1L
  umi[ok] <- substring(seqs[ok], e5[ok] + 1L, s3[ok] - 1L)
  outcome <- rep("no 5' flank", n)
  outcome[found5 & !found3] <- "no 3' flank"
  outcome[found5 & found3 & !ok] <- "flanks overlap"
  outcome[ok] <- "UMI extracted"
  ulen <- nchar(umi)
  rej <- ok & (ulen < len_min | ulen > len_max | grepl("N", umi, fixed = TRUE))
  outcome[rej] <- "UMI wrong length or contains N"
  umi[rej] <- NA_character_
  data.table::data.table(umi = umi, umi_outcome = outcome, umi_length = ifelse(ok, ulen, NA_integer_))
}

# ------------------------------------------------------------------------------
# 2. Directional clustering via a deletion-key neighbour index
# ------------------------------------------------------------------------------
.deletion_keys <- function(u, d) {
  cur <- data.table::data.table(id = seq_along(u), k = u)
  keys <- list(cur)
  for (step in seq_len(d)) {
    L <- nchar(cur$k)
    nxt <- cur[rep(seq_len(nrow(cur)), L)]
    pos <- unlist(lapply(L, seq_len))
    nxt[, k := paste0(substr(k, 1, pos - 1L), substr(k, pos + 1L, nchar(k)))]
    nxt <- unique(nxt)
    keys[[step + 1]] <- nxt
    cur <- nxt
  }
  unique(data.table::rbindlist(keys))
}

umi_neighbour_pairs <- function(u, d) {
  if (d <= 0 || length(u) < 2) return(data.table::data.table(a = integer(), b = integer(), dist = integer()))
  keys <- .deletion_keys(u, d)
  j <- keys[keys, on = "k", allow.cartesian = TRUE, nomatch = NULL][id < i.id, .(a = id, b = i.id)]
  j <- unique(j)
  if (nrow(j) == 0) return(data.table::data.table(a = integer(), b = integer(), dist = integer()))
  j[, dist := mapply(function(x, y) utils::adist(x, y)[1, 1], u[a], u[b])]
  j[dist <= d & dist > 0]
}

cluster_umis_directional <- function(umi_vec, d = 1L) {
  tb <- sort(table(umi_vec), decreasing = TRUE)
  u <- names(tb); cnt <- as.integer(tb)
  map <- setNames(u, u)
  if (d <= 0 || length(u) < 2) return(list(map = map, pairs = NULL, u = u, cnt = cnt))
  pr <- umi_neighbour_pairs(u, d)
  if (nrow(pr) == 0) return(list(map = map, pairs = pr, u = u, cnt = cnt))
  # directed edges hub -> node when count(hub) >= 2 * count(node) - 1
  e1 <- pr[cnt[a] >= 2L * cnt[b] - 1L, .(from = a, to = b)]
  e2 <- pr[cnt[b] >= 2L * cnt[a] - 1L, .(from = b, to = a)]
  edges <- rbind(e1, e2)
  adj <- split(edges$to, factor(edges$from, levels = seq_along(u)))
  # Work on integer positions (u[k] <-> position k): assigning into a named
  # character vector by name costs O(n) per assignment, O(n^2) in all
  hub <- seq_along(u)
  assigned <- rep(FALSE, length(u))
  for (i in which(lengths(adj) > 0L)) {   # u is sorted by abundance, so hubs come first
    if (assigned[i]) next
    queue <- i; assigned[i] <- TRUE
    while (length(queue)) {
      v <- queue[1]; queue <- queue[-1]
      nb <- adj[[v]]
      nb <- nb[!assigned[nb]]
      if (length(nb)) { assigned[nb] <- TRUE; hub[nb] <- i; queue <- c(queue, nb) }
    }
  }
  list(map = setNames(u[hub], u), pairs = pr, u = u, cnt = cnt)
}

# Expected fraction of UMIs with at least one other UMI within d substitutions by chance
expected_chance_neighbours <- function(n, L, d = 1L) {
  p1 <- sum(vapply(0:d, function(k) choose(L, k) * 3^k, numeric(1))) / 4^L
  1 - (1 - p1)^(n - 1)
}

# ------------------------------------------------------------------------------
# 3. Collapse statistics
# ------------------------------------------------------------------------------
.mode_num <- function(v) { rv <- round(v); ux <- unique(rv); as.numeric(ux[which.max(tabulate(match(rv, ux)))]) }

# Median of the highest window of width (2*w + wfrac*x) holding at least k reads and at
# least min_share of the molecule's reads. Stutter grows with length, so the window widens
# with it. PCR slippage mostly shortens reads, so the longest well-supported window is the
# best estimate of the template; min_share stops a few stray long reads (e.g. PCR
# chimeras) from setting the length of a large family.
.top_cluster <- function(v, k = 2L, w = 3, wfrac = 0.10, min_share = 0) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  need <- max(k, ceiling(min_share * length(v)))
  for (x in sort(v, decreasing = TRUE)) {
    win <- v[v >= x - (2 * w + wfrac * x) & v <= x]
    if (length(win) >= need) return(as.numeric(stats::median(win)))
  }
  as.numeric(stats::median(v))
}

# Read lengths are whole units, but a median of an even number of reads can fall halfway
# between two (e.g. 133.5). Such a molecule length is moved to the neighbouring whole unit
# that more of its reads sit on, or, if they are level, to one of the two at random. Always
# rounding up would bias molecules upwards, and R's round() (halves to the even number)
# piles molecules onto even lengths, which shows as a saw-tooth in length histograms.
.whole_unit <- function(v, reads) {
  if (!is.finite(v) || v - floor(v) != 0.5 || any(reads != floor(reads), na.rm = TRUE)) return(v)
  lo <- floor(v); hi <- lo + 1
  n_lo <- sum(reads == lo, na.rm = TRUE); n_hi <- sum(reads == hi, na.rm = TRUE)
  if (n_lo > n_hi) lo else if (n_hi > n_lo) hi else if (stats::runif(1) < 0.5) lo else hi
}

# The longest whole-unit length held by at least k reads
.max_supported <- function(v, k = 2L) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  tb <- table(floor(v + 0.5)); s <- as.numeric(names(tb))[tb >= k]
  if (length(s)) max(s) else as.numeric(stats::median(v))
}

collapse_fun <- function(stat, k = 2L, w = 3, wfrac = 0.10, min_share = 0) {
  switch(stat,
    top_cluster = function(v) .top_cluster(v, k, w, wfrac, min_share),
    max_supported = function(v) .max_supported(v, k),
    median = function(v) as.numeric(stats::median(v, na.rm = TRUE)),
    mean = function(v) mean(v, na.rm = TRUE),
    mode = function(v) .mode_num(v[is.finite(v)]),
    stop("Unknown umi_collapse_stat: ", stat, " (use top_cluster, max_supported, mode, median or mean)"))
}

# Every length rule the report compares, whatever umi_collapse_stat is set to
collapse_rules <- function(p) {
  k <- p$umi_top_cluster_k; w <- p$umi_top_cluster_w; wf <- p$umi_top_cluster_wfrac
  list(
    top_cluster_k = list(label = sprintf("top_cluster, %d reads", k), f = function(v) .top_cluster(v, k, w, wf, 0)),
    top_cluster_share20 = list(label = sprintf("top_cluster, %d reads and 20%%", k), f = function(v) .top_cluster(v, k, w, wf, 0.2)),
    max_supported = list(label = sprintf("longest length on %d reads", k), f = function(v) .max_supported(v, k)),
    mode = list(label = "mode", f = function(v) .mode_num(v[is.finite(v)])),
    median = list(label = "median", f = function(v) as.numeric(stats::median(v, na.rm = TRUE))),
    mean = list(label = "mean", f = function(v) mean(v, na.rm = TRUE)))
}

# How well each length rule behaves, judged from the molecules' own reads and needing no
# assumption about alleles. On a random subset of molecules (lower is better for all three):
#   split_half: give a molecule's reads at random to two halves and collapse each; the mean
#     difference between the halves (molecules with >= 6 reads). A rule that a few reads can
#     move gives two different answers.
#   depth_drift: collapse 20 reads and 3 reads of the same molecule; the mean difference
#     (molecules with >= 20 reads), signed: above 0 means more reads give a longer length. A rule whose answer depends on depth makes deep and shallow
#     samples incomparable.
#   chimera_shift: add stray reads (10% of the molecule's reads, at least 1) taken from a random
#     other molecule of the same sample, as a PCR chimera would, and collapse again; the share
#     of molecules that move by more than 2 units (molecules with >= 6 reads).
collapse_scorecard <- function(rd, rules, n_max = 3000L, seed = 1L) {
  set.seed(seed)
  rd <- rd[is.finite(length)]
  by_mol <- split(rd$length, rd$molecule_id)
  smp <- rd[, .(file_stem = file_stem[1], n = .N), by = molecule_id]
  pick <- function(min_n) { ids <- smp[n >= min_n, molecule_id]; if (length(ids) > n_max) sample(ids, n_max) else ids }
  ids6 <- pick(6L); ids20 <- pick(20L)
  # a random donor molecule from the same sample for each molecule tested for chimeras
  donor <- smp[, .(molecule_id, donor = if (.N > 1) molecule_id[(seq_len(.N) + sample.int(.N - 1, 1) - 1) %% .N + 1] else NA_character_), by = file_stem]
  donor <- setNames(donor$donor, donor$molecule_id)
  rbindlist(lapply(names(rules), function(r) {
    f <- rules[[r]]$f
    sh <- if (length(ids6)) mean(vapply(ids6, function(m) { v <- by_mol[[m]]; i <- sample.int(length(v)); h <- length(v) %/% 2
      abs(f(v[i[seq_len(h)]]) - f(v[i[(h + 1):(2 * h)]])) }, numeric(1)), na.rm = TRUE) else NA_real_
    dd <- if (length(ids20)) mean(vapply(ids20, function(m) { v <- by_mol[[m]]; f(sample(v, 20L)) - f(sample(v, 3L)) }, numeric(1)), na.rm = TRUE) else NA_real_
    cs <- if (length(ids6)) mean(vapply(ids6, function(m) { v <- by_mol[[m]]; d <- donor[[m]]
      if (is.na(d)) return(NA_real_); dv <- by_mol[[d]]; extra <- dv[sample.int(length(dv), max(1L, ceiling(0.1 * length(v))), replace = TRUE)]
      abs(f(c(v, extra)) - f(v)) > 2 }, numeric(1)), na.rm = TRUE) else NA_real_
    data.table::data.table(rule = r, label = rules[[r]]$label, split_half = sh, depth_drift = dd, chimera_shift = cs,
                           n_split_half = length(ids6), n_depth_drift = length(ids20))
  }))
}

# ------------------------------------------------------------------------------
# One sample: reads with UMIs -> molecules
# ------------------------------------------------------------------------------
# reads: data.table(qname, umi)   lengths: long data.table(qname, measure, length)
# primary: the measure whose length is compared across collapse rules
form_molecules <- function(sample, reads, lengths, p, primary) {
  r <- reads[!is.na(umi)]
  empty <- list(reads = NULL, molecules = NULL, molecule_lengths = NULL, clustering = NULL, collapse_compare = NULL)
  if (nrow(r) == 0) return(empty)
  cl <- cluster_umis_directional(r$umi, p$umi_cluster_max_dist)
  r[, umi_group := unname(cl$map[umi])]
  # One UMI family is one molecule. Repeat length plays no part in forming molecules:
  # PCR slippage and deletion can make reads of one molecule differ by tens of units, which
  # is the error that collapsing is meant to remove, so length cannot say which reads belong
  # together.
  r[, molecule_id := paste(sample, umi_group, sep = "|")]
  mol <- r[, .(umi = umi_group[1], n_reads = .N), by = molecule_id]
  mol[, kept := n_reads >= p$umi_min_reads]
  mol[, file_stem := sample]

  f <- collapse_fun(p$umi_collapse_stat, p$umi_top_cluster_k, p$umi_top_cluster_w, p$umi_top_cluster_wfrac,
                    p$umi_top_cluster_min_share)
  keep_ids <- mol$molecule_id[mol$kept]
  lk <- lengths[r[molecule_id %in% keep_ids, .(qname, molecule_id)], on = "qname", nomatch = NULL]
  set.seed(sum(utf8ToInt(sample)))   # reproducible tie-breaks in .whole_unit()
  ml <- lk[, .(length = .whole_unit(f(length), length), n_reads = .N, read_sd = stats::sd(length)), by = .(molecule_id, measure)]
  ml[, file_stem := sample]

  # The primary measure's molecule length by other collapse rules, for the report's comparison
  prim <- lk[measure == primary]
  rules <- collapse_rules(p)
  cmp <- prim[, c(list(n_reads = .N), lapply(rules, function(g) .whole_unit(g$f(length), length))), by = molecule_id]
  cmp[, file_stem := sample]

  pairs <- cl$pairs
  clustering <- data.frame(file_stem = sample, n_reads_with_umi = nrow(r),
                           distinct_umis = length(cl$u), umi_groups = data.table::uniqueN(r$umi_group),
                           umis_with_neighbour = if (is.null(pairs) || nrow(pairs) == 0) 0L else data.table::uniqueN(c(pairs$a, pairs$b)),
                           expected_by_chance = expected_chance_neighbours(length(cl$u), stats::median(nchar(cl$u)), max(1, p$umi_cluster_max_dist)) * length(cl$u),
                           molecules = nrow(mol), molecules_kept = sum(mol$kept), stringsAsFactors = FALSE)
  list(reads = r[, .(file_stem = sample, qname, umi, umi_group, molecule_id)],
       molecules = mol, molecule_lengths = ml, clustering = clustering, collapse_compare = cmp)
}
