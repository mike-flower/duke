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
#   3. Each read is labelled with an allele by splitting the primary measure's
#      read distribution at the valley between the two alleles. A family with
#      at least 2 reads on each allele is two molecules sharing a barcode and is
#      split (umi_split_mixed); otherwise its allele is its reads' majority.
#   4. Molecules with fewer than umi_min_reads reads are dropped, and each
#      measure is collapsed to one value per molecule (default top_cluster).
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
# 3. Alleles from the read distribution: centres and the valley between them
# ------------------------------------------------------------------------------
detect_alleles <- function(x, min_sep = 20, min_frac = 0.002, min_concentration = 0.55, min_dip_ratio = 5) {
  x <- x[is.finite(x)]
  single <- function(c1) list(centers = c1, valley = NA_real_)
  if (length(x) < 50) return(single(stats::median(x)))
  dmode <- function(z) { d <- stats::density(z, n = 1024); d$x[which.max(d$y)] }
  c1 <- dmode(x)
  hi <- x[x > c1 + min_sep]; lo <- x[x < c1 - min_sep]
  side <- if (length(hi) >= length(lo)) hi else lo
  if (length(side) < max(20, min_frac * length(x))) return(single(c1))
  c2 <- dmode(side)
  if (mean(abs(side - c2) <= max(5, 0.15 * abs(c2))) < min_concentration) return(single(c1))
  centers <- sort(c(c1, c2))
  # Repeat lengths are whole numbers, so a narrow bandwidth gives a dip between every
  # integer; a bandwidth of at least one unit gives the dip between the alleles only
  d <- stats::density(x, n = 2048, bw = max(stats::bw.nrd0(x), 1))
  idx <- which(d$x > centers[1] & d$x < centers[2])
  if (!length(idx)) return(single(c1))
  # The valley is the first point above the lower allele where the density reaches the
  # floor of the dip: its minimum, or 1% of the smaller peak where the gap between the
  # alleles is empty. Taking the first such point (rather than the middle of the flat
  # region) keeps the boundary just above the lower allele when the alleles are far apart,
  # so reads between them (e.g. contracted copies of the long allele) go with the long allele
  # and the boundary never cuts into the long allele's own stutter.
  h_small <- min(d$y[which.min(abs(d$x - centers[1]))], d$y[which.min(abs(d$x - centers[2]))])
  floor_y <- max(min(d$y[idx]) * 1.05 + 1e-12, 0.01 * h_small)
  valley <- d$x[idx[which(d$y[idx] <= floor_y)[1]]]
  h2 <- d$y[which.min(abs(d$x - c2))]; hv <- min(d$y[idx])
  if (!(h2 > min_dip_ratio * hv)) return(single(c1))
  list(centers = centers, valley = valley)
}

# ------------------------------------------------------------------------------
# 4. Collapse statistics
# ------------------------------------------------------------------------------
.mode_num <- function(v) { rv <- round(v); ux <- unique(rv); as.numeric(ux[which.max(tabulate(match(rv, ux)))]) }

# Median of the highest window of width (2*w + wfrac*x) holding at least k reads.
# Stutter grows with length, so the window widens with it.
.top_cluster <- function(v, k = 2L, w = 3, wfrac = 0.10) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  for (x in sort(v, decreasing = TRUE)) {
    win <- v[v >= x - (2 * w + wfrac * x) & v <= x]
    if (length(win) >= k) return(as.numeric(stats::median(win)))
  }
  as.numeric(stats::median(v))
}

collapse_fun <- function(stat, k = 2L, w = 3, wfrac = 0.10) {
  switch(stat,
    top_cluster = function(v) .top_cluster(v, k, w, wfrac),
    median = function(v) as.numeric(stats::median(v, na.rm = TRUE)),
    mean = function(v) mean(v, na.rm = TRUE),
    mode = function(v) .mode_num(v[is.finite(v)]),
    stop("Unknown umi_collapse_stat: ", stat))
}

# ------------------------------------------------------------------------------
# One sample: reads with UMIs -> molecules
# ------------------------------------------------------------------------------
# reads: data.table(qname, umi, primary)   lengths: long data.table(qname, measure, length)
form_molecules <- function(sample, reads, lengths, p) {
  r <- reads[!is.na(umi)]
  empty <- list(reads = NULL, molecules = NULL, molecule_lengths = NULL, clustering = NULL)
  if (nrow(r) == 0) return(empty)
  cl <- cluster_umis_directional(r$umi, p$umi_cluster_max_dist)
  r[, umi_group := unname(cl$map[umi])]

  al <- detect_alleles(r$primary, min_sep = p$umi_allele_min_sep)
  r[, read_allele := if (is.na(al$valley)) "single" else ifelse(primary < al$valley, "short", "long")]

  fam <- r[, .(n_short = sum(read_allele == "short"), n_long = sum(read_allele == "long"), n = .N), by = umi_group]
  fam[, mixed := pmin(n_short, n_long) >= 2L]
  fam[, family_allele := data.table::fifelse(n_short > n_long, "short",
                          data.table::fifelse(n_long > n_short, "long",
                          data.table::fifelse(n_long == 0 & n_short == 0, "single", "long")))]
  r <- fam[r, on = "umi_group"]
  split_it <- isTRUE(p$umi_split_mixed) & r$mixed
  r[, molecule_allele := data.table::fifelse(split_it, read_allele, family_allele)]
  r[, molecule_id := paste(sample, umi_group, molecule_allele, sep = "|")]

  mol <- r[, .(umi = umi_group[1], allele = molecule_allele[1], n_reads = .N,
               mixed_umi = mixed[1],
               allele_agreement = mean(read_allele == molecule_allele[1] | read_allele == "single")),
           by = molecule_id]
  mol[, kept := n_reads >= p$umi_min_reads]
  mol[, file_stem := sample]

  f <- collapse_fun(p$umi_collapse_stat, p$umi_top_cluster_k, p$umi_top_cluster_w, p$umi_top_cluster_wfrac)
  keep_ids <- mol$molecule_id[mol$kept]
  lk <- lengths[r[molecule_id %in% keep_ids, .(qname, molecule_id)], on = "qname", nomatch = NULL]
  ml <- lk[, .(length = f(length), n_reads = .N, read_sd = stats::sd(length)), by = .(molecule_id, measure)]
  ml[, file_stem := sample]

  pairs <- cl$pairs
  clustering <- data.frame(file_stem = sample, n_reads_with_umi = nrow(r),
                           distinct_umis = length(cl$u), umi_groups = data.table::uniqueN(r$umi_group),
                           umis_with_neighbour = if (is.null(pairs) || nrow(pairs) == 0) 0L else data.table::uniqueN(c(pairs$a, pairs$b)),
                           expected_by_chance = expected_chance_neighbours(length(cl$u), stats::median(nchar(cl$u)), max(1, p$umi_cluster_max_dist)) * length(cl$u),
                           valley = al$valley, allele_centres = paste(round(al$centers, 1), collapse = ", "),
                           molecules = nrow(mol), molecules_kept = sum(mol$kept),
                           mixed_families = sum(fam$mixed), stringsAsFactors = FALSE)
  list(reads = r[, .(file_stem = sample, qname, umi, umi_group, read_allele, molecule_id)],
       molecules = mol, molecule_lengths = ml, clustering = clustering)
}
