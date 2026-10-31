# ==============================================================================
# Duke - allele structure: clustering, consensus, variants, waterfalls
# ==============================================================================

# ------------------------------------------------------------------------------
# Repeat-length clustering
# ------------------------------------------------------------------------------
# The number of clusters is chosen by a Gaussian mixture (BIC) on a subsample,
# the full data are clustered with k-means at that k, and any cluster holding
# fewer than min_prop of the reads is merged into its nearest neighbour. The
# merge step stops the expansion smear of an unstable allele being cut into
# many small clusters (the old version gave 3-20 clusters for 15% of samples).
cluster_by_repeat_length <- function(values, cluster_max = 6L, min_prop = 0.05, downsample = 10000L) {
  out <- rep(NA_integer_, length(values))
  ok <- which(is.finite(values))
  v <- values[ok]
  if (length(unique(v)) < 2) { out[ok] <- 1L; return(out) }
  set.seed(123)
  opt_v <- if (length(v) > downsample) sample(v, downsample) else v
  k <- tryCatch({
    utils::capture.output(g <- mclust::Mclust(opt_v, G = 1:cluster_max, verbose = FALSE))
    g$G
  }, error = function(e) 2L)
  k <- min(k, length(unique(v)))
  cl <- if (k <= 1) rep(1L, length(v)) else { set.seed(123); stats::kmeans(v, centers = k, nstart = 10)$cluster }
  repeat {
    tab <- table(cl); props <- tab / sum(tab)
    if (length(tab) <= 1 || min(props) >= min_prop) break
    small <- as.integer(names(props)[which.min(props)])
    centres <- tapply(v, cl, mean)
    others <- centres[names(centres) != small]
    target <- as.integer(names(others)[which.min(abs(others - centres[as.character(small)]))])
    cl[cl == small] <- target
  }
  centres <- tapply(v, cl, mean)
  cl <- match(cl, as.integer(names(sort(centres))))   # 1 = shortest
  out[ok] <- cl
  out
}

# ------------------------------------------------------------------------------
# Haplotype clustering on flank sequence (bases next to the repeat are kept)
# ------------------------------------------------------------------------------
haplotype_sequences <- function(left, right, region = "both", trim = "auto") {
  # "auto": the length 95% of reads reach, less 20 bp, so every kept base is
  # covered by nearly all reads (alignment ends jitter by a few bases, and the
  # k-mers at a ragged end would otherwise look like a haplotype difference)
  modal <- function(x) { n <- nchar(x[!is.na(x) & nzchar(x)]); if (!length(n)) NA_integer_ else max(0L, as.integer(stats::quantile(n, 0.05, type = 1)) - 20L) }
  tl <- tr <- NA_integer_
  if (!is.na(trim) && tolower(trim) == "auto") { tl <- modal(left); tr <- modal(right) }
  else if (!is.na(suppressWarnings(as.numeric(trim)))) { tl <- tr <- ceiling(as.numeric(trim) / 2) }
  keep_end <- function(x, n) ifelse(!is.na(n) & nchar(x) > n, substr(x, nchar(x) - n + 1, nchar(x)), x)
  keep_start <- function(x, n) ifelse(!is.na(n) & nchar(x) > n, substr(x, 1, n), x)
  left <- keep_end(ifelse(is.na(left), "", left), tl)
  right <- keep_start(ifelse(is.na(right), "", right), tr)
  # The two flanks are joined with an N so that no k-mer spans the junction. A k-mer
  # across the junction would hold the last bases of one flank and the first of the
  # other, so any jitter where the tract starts or ends would change ~k k-mers in the
  # same reads at once – the signature of a real variant – and split one allele into
  # spurious haplotypes (seen as a ~10% minor "haplotype" of the expanded allele that
  # differed only in the bases next to the tract).
  switch(region, left = left, right = right, both = paste0(left, "N", right))
}

# Haplotypes are found from k-mers rather than all-pairs edit distance, which
# on ~2.5 kb flanks costs ~1e13 operations per sample. Reads from different
# haplotypes differ in the k-mers that span each SNP or indel, so:
#   1. k-mers (k = 15) are counted in a subsample of reads; those present in
#      between min_frac and 1 - min_frac of reads are informative (k-mers
#      everyone shares, and rare ones made by sequencing errors, drop out)
#   2. each read becomes a presence/absence profile over the informative
#      k-mers (Aho-Corasick matching, Biostrings::vwhichPDict)
#   3. up to `downsample` reads are clustered (Manhattan distance on the
#      profiles, complete linkage); a split is accepted only when the mean
#      silhouette exceeds 0.5 and every cluster holds at least min_prop of reads
#   4. every other read joins the cluster whose medoid profile is closest
.kmer_table <- function(s, k) {
  w <- pmax(nchar(s) - k + 1L, 0L)
  r <- rep.int(seq_along(s), w)
  st <- sequence(w)
  kt <- unique(data.table::data.table(r = r, kmer = substring(s[r], st, st + k - 1L)))
  kt[!grepl("N", kmer, fixed = TRUE)]  # k-mers across the flank junction or an N base
}

.kmer_profile <- function(s, kmers, pd) {
  hits <- Biostrings::vwhichPDict(pd, Biostrings::DNAStringSet(s))
  M <- matrix(0L, length(s), length(kmers))
  n <- lengths(hits)
  if (sum(n)) M[cbind(rep.int(seq_along(s), n), unlist(hits))] <- 1L
  M
}

cluster_by_haplotype <- function(seqs, cluster_max = 6L, downsample = 3000L, k = 15L,
                                 min_frac = 0.1, min_prop = 0.05, discover_n = 1000L, max_kmers = 2000L) {
  out <- rep(NA_integer_, length(seqs))
  # N marks the junction between the flanks (and any uncalled base); k-mers containing
  # it are not used, so it does not exclude the read
  ok <- which(!is.na(seqs) & nchar(seqs) >= k & !grepl("[^ACGTN]", seqs))
  if (length(ok) < 10) { out[ok] <- 1L; return(out) }
  s <- seqs[ok]
  set.seed(123)
  disc <- if (length(s) > discover_n) sample(length(s), discover_n) else seq_along(s)
  kt <- .kmer_table(s[disc], k)
  f <- kt[, .(frac = .N / length(disc)), by = kmer][frac >= min_frac & frac <= 1 - min_frac]
  if (!nrow(f)) { out[ok] <- 1L; return(out) }
  if (nrow(f) > max_kmers) f <- f[order(abs(frac - 0.5))][seq_len(max_kmers)]
  kmers <- f$kmer
  pd <- Biostrings::PDict(Biostrings::DNAStringSet(kmers))
  idx <- if (length(s) > downsample) sort(sample(length(s), downsample)) else seq_along(s)
  M0 <- .kmer_profile(s[idx], kmers, pd)
  # A real variant gives ~k k-mers per allele present in exactly the same reads
  # (and the two alleles' k-mers are mirror images); sequencing-error k-mers do
  # not. Keep k-mers with at least k/2 strongly correlated partners, group them
  # into variant blocks, and score each read on each block.
  keep <- which(apply(M0, 2, stats::var) > 0)
  if (length(keep) < 2) { out[ok] <- 1L; return(out) }
  C <- suppressWarnings(stats::cor(M0[, keep, drop = FALSE])); C[is.na(C)] <- 0
  strong <- abs(C) >= 0.8
  keep2 <- which(rowSums(strong) - 1 >= k %/% 2)
  if (length(keep2) < 2) { out[ok] <- 1L; return(out) }
  C2 <- C[keep2, keep2, drop = FALSE]
  blk <- stats::cutree(stats::hclust(stats::as.dist(1 - abs(C2)), method = "average"), h = 0.2)
  kept <- keep[keep2]
  # orient each k-mer to its block's first member, so ref- and alt-allele k-mers agree
  sgn <- vapply(seq_along(kept), function(j) { f <- which(blk == blk[j])[1]; if (C2[j, f] < 0) -1 else 1 }, numeric(1))
  # Score = share of the block's k-mers present in the read that belong to the
  # first allele. Only k-mers actually present count: a sequencing error near a
  # variant removes k-mers of both alleles, and absence is not evidence of the
  # other allele. A read with none of the block's k-mers scores 0.5.
  block_scores <- function(P) {
    P <- P[, kept, drop = FALSE]
    vapply(sort(unique(blk)), function(b) {
      j <- blk == b
      a1 <- rowSums(P[, j & sgn > 0, drop = FALSE]); a2 <- rowSums(P[, j & sgn < 0, drop = FALSE])
      ifelse(a1 + a2 > 0, a1 / (a1 + a2), 0.5)
    }, numeric(nrow(P)))
  }
  M <- block_scores(M0); if (!is.matrix(M)) M <- matrix(M, ncol = 1)
  d <- stats::dist(M, method = "manhattan")
  if (max(d) == 0) { out[ok] <- 1L; return(out) }
  hc <- stats::hclust(d, method = "complete")
  k_best <- 1L; best <- 0.5
  for (kk in seq.int(2L, max(2L, min(cluster_max, length(idx) - 1)))) {
    if (kk > cluster_max) break
    cl <- stats::cutree(hc, k = kk)
    if (min(table(cl)) / length(cl) < min_prop) next
    sil <- mean(cluster::silhouette(cl, d)[, "sil_width"])
    if (sil > best) { best <- sil; k_best <- kk }
  }
  sub_cl <- stats::cutree(hc, k = k_best)
  if (k_best == 1L) { out[ok] <- 1L; return(out) }
  dm <- as.matrix(d)
  medoids <- M[vapply(seq_len(k_best), function(c) { m <- which(sub_cl == c); m[which.min(rowSums(dm[m, m, drop = FALSE]))] }, integer(1)), , drop = FALSE]
  full <- integer(length(s)); full[idx] <- sub_cl
  rest <- setdiff(seq_along(s), idx)
  for (ch in split(rest, ceiling(seq_along(rest) / 5000))) {
    R <- block_scores(.kmer_profile(s[ch], kmers, pd)); if (!is.matrix(R)) R <- matrix(R, nrow = length(ch))
    dd <- vapply(seq_len(k_best), function(c) rowSums(abs(sweep(R, 2, medoids[c, ]))), numeric(length(ch)))
    full[ch] <- if (is.matrix(dd)) max.col(-dd, ties.method = "first") else which.min(dd)
  }
  out[ok] <- full
  out
}

# Cluster one sample's reads. Returns cluster_number, repeat_cluster, haplotype_cluster.
cluster_sample <- function(d, p) {
  by <- p$cluster_by
  n <- nrow(d)
  rc <- if ("repeat" %in% by) cluster_by_repeat_length(d$length, p$repeat_cluster_max, p$repeat_cluster_min_prop, p$repeat_cluster_downsample) else rep(1L, n)
  hc <- rep(1L, n)
  if ("haplotype" %in% by) {
    hs <- haplotype_sequences(d$left, d$right, p$haplotype_region, p$haplotype_trim_length)
    if (identical(by[1], "repeat")) {
      for (r in unique(rc[!is.na(rc)])) { i <- which(rc == r); hc[i] <- cluster_by_haplotype(hs[i], p$haplotype_cluster_max, p$haplotype_cluster_downsample) }
    } else {
      hc <- cluster_by_haplotype(hs, p$haplotype_cluster_max, p$haplotype_cluster_downsample)
      if ("repeat" %in% by) for (h in unique(hc[!is.na(hc)])) { i <- which(hc == h); rc[i] <- cluster_by_repeat_length(d$length[i], p$repeat_cluster_max, p$repeat_cluster_min_prop, p$repeat_cluster_downsample) }
    }
  }
  rc[is.na(rc)] <- 0L; hc[is.na(hc)] <- 0L   # 0 = no repeat length / no usable flank sequence
  key <- paste0("R", rc, "-H", hc)
  # cluster 1 = most reads
  ord <- names(sort(table(key), decreasing = TRUE))
  data.frame(cluster_number = match(key, ord), repeat_cluster = paste0("R", rc), haplotype_cluster = paste0("H", hc))
}

# ------------------------------------------------------------------------------
# Consensus (DECIPHER)
# ------------------------------------------------------------------------------
# Terminal gaps count as gaps, so the consensus stops where most reads stop
# rather than running on to the end of the longest read in the sample.
consensus_function <- function(seqs, threshold = 0.5, sample_size = 50L, minor_share = NA) {
  seqs <- seqs[!is.na(seqs) & nchar(seqs) > 0]
  if (length(seqs) == 0) return(list(consensus = NA_character_, n = 0L, minor = NULL))
  if (!is.na(sample_size) && length(seqs) > sample_size) { set.seed(123); seqs <- sample(seqs, sample_size) }
  if (length(seqs) == 1) return(list(consensus = seqs, n = 1L, minor = NULL))
  dss <- Biostrings::DNAStringSet(seqs)
  aln <- suppressMessages(DECIPHER::AlignSeqs(dss, processors = 1, verbose = FALSE))
  cons_g <- as.character(suppressMessages(DECIPHER::ConsensusSequence(aln, threshold = threshold, ambiguity = FALSE,
                                                                      noConsensusChar = "N", includeTerminalGaps = TRUE)))
  cons_g <- gsub("[^ACGT-]", "N", cons_g)
  # Substitutions carried by a minority of the reads (at least minor_share): returned as
  # the position in the (gapless) consensus and the base, at most the 20 strongest. Only
  # bases are considered, not gaps: a gap column near a read end is mostly where reads
  # stop, not a deletion.
  minor <- NULL
  if (!is.na(minor_share) && length(seqs) >= 5) {
    cm <- Biostrings::consensusMatrix(aln, as.prob = FALSE)
    acgt <- cm[intersect(c("A", "C", "G", "T"), rownames(cm)), , drop = FALSE]
    n_col <- colSums(cm)
    cg <- strsplit(cons_g, "")[[1]]
    cons_pos <- cumsum(cg != "-")
    alts <- list()
    for (b in rownames(acgt)) {
      cols <- which(acgt[b, ] / n_col >= minor_share & acgt[b, ] >= 3 & cg != b & cg != "-" & cg != "N")
      if (length(cols)) alts[[b]] <- data.frame(cons_pos = cons_pos[cols], base = b, sample_share = acgt[b, cols] / n_col[cols], stringsAsFactors = FALSE)
    }
    if (length(alts)) { minor <- do.call(rbind, alts); minor <- utils::head(minor[order(-minor$sample_share), , drop = FALSE], 20) }
  }
  list(consensus = gsub("-", "", cons_g), n = length(seqs), minor = minor)
}

# ------------------------------------------------------------------------------
# Variants in flank consensus sequences, in coordinates of the full reference
# ------------------------------------------------------------------------------
# Each flank consensus is aligned to its own reference flank with free end gaps
# on both (overlap), so the repeat never enters the alignment, and any sequence
# between the primary tract and the flank (the rest of an N block that spans
# several tracts, e.g. CAA CAG CCG ... CCT after CAG) is left unaligned rather
# than forced onto the flank as a run of false variants.
# Positions map directly onto the reference. Runs of adjacent inserted or
# deleted bases are merged into single events. Positions are 1-based on the
# full reference including its N block, ready for IGV alongside that FASTA.
call_flank_variants <- function(consensus, ref_flank, offset, edge = 3L) {
  if (is.na(consensus) || nchar(consensus) < 10) return(NULL)
  ns <- if (requireNamespace("pwalign", quietly = TRUE)) "pwalign" else "Biostrings"
  pw <- lapply(setNames(nm = c("pairwiseAlignment", "alignedPattern", "alignedSubject", "subject")),
               function(f) getExportedValue(ns, f))
  aln <- pw$pairwiseAlignment(Biostrings::DNAString(consensus), Biostrings::DNAString(ref_flank),
                              type = "overlap", gapOpening = 5, gapExtension = 2)
  q <- strsplit(as.character(pw$alignedPattern(aln)), "")[[1]]
  r <- strsplit(as.character(pw$alignedSubject(aln)), "")[[1]]
  rpos <- BiocGenerics::start(pw$subject(aln)) - 1L + cumsum(r != "-")
  # Read contexts for each event: `flank_bp` bases of consensus either side, joined to the
  # alternative bases (alt_context) or to the reference bases (ref_context). Counting reads
  # that contain each gives the variant's read support.
  qn <- which(q != "-")
  ctx <- function(i, j, flank_bp = 8L) {
    before <- q[qn[qn < i]]; after <- q[qn[qn > j]]
    left <- paste(utils::tail(before, flank_bp), collapse = ""); right <- paste(utils::head(after, flank_bp), collapse = "")
    data.frame(alt_context = paste0(left, paste(q[i:j][q[i:j] != "-"], collapse = ""), right),
               ref_context = paste0(left, paste(r[i:j][r[i:j] != "-"], collapse = ""), right), stringsAsFactors = FALSE)
  }
  ev <- list(); k <- 1L; i <- 1L; n <- length(q)
  while (i <= n) {
    if (r[i] != "-" && q[i] != "-" && r[i] != q[i] && q[i] != "N" && r[i] != "N") {
      ev[[k]] <- data.frame(pos = rpos[i], ref = r[i], alt = q[i], type = "SNV", ctx(i, i)); k <- k + 1L; i <- i + 1L
    } else if (q[i] == "-") {
      j <- i; while (j < n && q[j + 1] == "-") j <- j + 1L
      anchor_pos <- rpos[i] - 1L
      anchor <- if (anchor_pos >= 1) substr(ref_flank, anchor_pos, anchor_pos) else ""
      ev[[k]] <- data.frame(pos = anchor_pos, ref = paste0(anchor, paste(r[i:j], collapse = "")), alt = anchor, type = "DEL", ctx(i, j)); k <- k + 1L
      i <- j + 1L
    } else if (r[i] == "-") {
      j <- i; while (j < n && r[j + 1] == "-") j <- j + 1L
      anchor_pos <- if (i > 1) rpos[i - 1] else 0L
      anchor <- if (anchor_pos >= 1) substr(ref_flank, anchor_pos, anchor_pos) else ""
      ev[[k]] <- data.frame(pos = anchor_pos, ref = anchor, alt = paste0(anchor, paste(q[i:j], collapse = "")), type = "INS", ctx(i, j)); k <- k + 1L
      i <- j + 1L
    } else i <- i + 1L
  }
  if (!length(ev)) return(NULL)
  v <- do.call(rbind, ev)
  # Drop calls within `edge` bases of either end of the flank: read ends and
  # the tract boundary are where consensus alignment is least reliable
  v <- v[v$pos > edge & v$pos < nchar(ref_flank) - edge & !grepl("N", v$alt), , drop = FALSE]
  v$pos <- v$pos + offset
  v
}

# Reference position (on the full reference, i.e. plus `offset`) of every base of a
# flank consensus, from the same overlap alignment the variants use; NA for bases
# outside the flank or inserted relative to it.
consensus_ref_map <- function(consensus, ref_flank, offset) {
  if (is.na(consensus) || nchar(consensus) < 10) return(integer(0))
  ns <- if (requireNamespace("pwalign", quietly = TRUE)) "pwalign" else "Biostrings"
  pa <- getExportedValue(ns, "pairwiseAlignment")
  a <- pa(Biostrings::DNAString(consensus), Biostrings::DNAString(ref_flank), type = "overlap", gapOpening = 5, gapExtension = 2)
  q <- strsplit(as.character(getExportedValue(ns, "alignedPattern")(a)), "")[[1]]
  r <- strsplit(as.character(getExportedValue(ns, "alignedSubject")(a)), "")[[1]]
  qs <- BiocGenerics::start(getExportedValue(ns, "pattern")(a)); rs <- BiocGenerics::start(getExportedValue(ns, "subject")(a))
  out <- rep(NA_integer_, nchar(consensus))
  qi <- qs - 1L + cumsum(q != "-"); ri <- rs - 1L + cumsum(r != "-")
  ok <- q != "-" & r != "-"
  out[qi[ok]] <- ri[ok] + offset
  out
}

# Minority substitutions of one allele flank as variant rows (pos, ref, alt, type and
# read contexts), from the positions consensus_function found.
minor_variants <- function(consensus, minor, ref_flank, offset, edge = 3L, flank_bp = 8L) {
  if (is.null(minor) || !nrow(minor) || is.na(consensus)) return(NULL)
  map <- consensus_ref_map(consensus, ref_flank, offset)
  if (!length(map)) return(NULL)
  out <- lapply(seq_len(nrow(minor)), function(j) {
    cp <- minor$cons_pos[j]; rp <- map[cp]
    if (is.na(rp) || rp - offset <= edge || rp - offset >= nchar(ref_flank) - edge) return(NULL)
    rb <- substr(ref_flank, rp - offset, rp - offset)
    if (rb == minor$base[j]) return(NULL)   # the minority base is the reference: not a variant
    left <- substr(consensus, max(1, cp - flank_bp), cp - 1); right <- substr(consensus, cp + 1, cp + flank_bp)
    data.frame(pos = rp, ref = rb, alt = minor$base[j], type = "SNV",
               alt_context = paste0(left, minor$base[j], right), ref_context = paste0(left, substr(consensus, cp, cp), right), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  if (is.null(out) || !nrow(out)) NULL else unique(out)
}

# ------------------------------------------------------------------------------
# Alleles: the sequence between the reference flanks, repeat notation, read support
# ------------------------------------------------------------------------------
# The part of an allele between the two reference flanks: the end of the left consensus
# after the left flank, the tract, and the start of the right consensus before the right
# flank (for a reference whose N block spans several tracts, e.g. CAG CAA CAG CCG ... CCT,
# that is everything in the block). Positions come from the same overlap alignments used
# for the variants.
between_flanks <- function(cons_left, cons_mid, cons_right, left_ref, right_ref, n_end) {
  ns <- if (requireNamespace("pwalign", quietly = TRUE)) "pwalign" else "Biostrings"
  pa <- getExportedValue(ns, "pairwiseAlignment"); pat <- getExportedValue(ns, "pattern"); sbj <- getExportedValue(ns, "subject")
  tail_l <- ""; head_r <- ""; left_from <- NA_integer_; right_to <- NA_integer_
  if (!is.na(cons_left) && nchar(cons_left) >= 10) {
    a <- pa(Biostrings::DNAString(cons_left), Biostrings::DNAString(left_ref), type = "overlap", gapOpening = 5, gapExtension = 2)
    e <- BiocGenerics::end(pat(a)); if (e < nchar(cons_left)) tail_l <- substring(cons_left, e + 1)
    left_from <- BiocGenerics::start(sbj(a))
  }
  if (!is.na(cons_right) && nchar(cons_right) >= 10) {
    a <- pa(Biostrings::DNAString(cons_right), Biostrings::DNAString(right_ref), type = "overlap", gapOpening = 5, gapExtension = 2)
    b <- BiocGenerics::start(pat(a)); if (b > 1) head_r <- substr(cons_right, 1, b - 1)
    right_to <- n_end + BiocGenerics::end(sbj(a))
  }
  list(sequence = paste0(tail_l, ifelse(is.na(cons_mid), "", cons_mid), head_r), left_from = left_from, right_to = right_to)
}

# Repeat notation, e.g. "(CAG)17 CAACAG CCGCCA (CCG)7 CCTCCT": runs of at least two copies
# of any tract motif are written as (motif)n, everything else as plain sequence.
repeat_notation <- function(x, motifs) {
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  motifs <- unique(motifs[!is.na(motifs) & nzchar(motifs)])
  hits <- do.call(rbind, lapply(motifs, function(m) {
    g <- gregexpr(paste0("(", m, "){2,}"), x)[[1]]
    if (g[1] == -1) return(NULL)
    data.frame(start = as.integer(g), len = attr(g, "match.length"), motif = m, stringsAsFactors = FALSE)
  }))
  if (is.null(hits)) return(x)
  hits <- hits[order(hits$start, -hits$len), ]
  out <- character(0); pos <- 1L
  for (k in seq_len(nrow(hits))) {
    h <- hits[k, ]; if (h$start < pos) next
    if (h$start > pos) out <- c(out, substr(x, pos, h$start - 1))
    out <- c(out, sprintf("(%s)%d", h$motif, h$len %/% nchar(h$motif)))
    pos <- h$start + h$len
  }
  if (pos <= nchar(x)) out <- c(out, substring(x, pos))
  paste(out, collapse = " ")
}

# Read support for a variant in one allele's reads: reads containing the alternative
# context against reads containing the reference context (reads with neither, e.g. an
# error nearby or a read that stops short, are not counted).
variant_support <- function(reads, alt_context, ref_context) {
  reads <- reads[!is.na(reads)]
  n_alt <- sum(grepl(alt_context, reads, fixed = TRUE)); n_ref <- sum(grepl(ref_context, reads, fixed = TRUE))
  c(alt_reads = n_alt, ref_reads = n_ref, support = if (n_alt + n_ref > 0) n_alt / (n_alt + n_ref) else NA_real_)
}

# Phased VCF, one column per sample. Each sample's alleles 1 and 2 are the two haplotypes
# (GT a|b); the region between the reference flanks is one record whose ALT alleles are
# the allele sequences found in the cohort, with MC (motif count of each tract) and AL
# (length in bp) per haplotype. Flank variants carry HAF (share of each haplotype's
# reads with the variant) and HDP (reads counted per haplotype).
write_allele_vcf <- function(alleles, variants, ref, motifs_txt, path) {
  samples <- unique(alleles$file_stem)
  hdr <- c("##fileformat=VCFv4.2", paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")), "##source=Duke",
           paste0("##reference=", basename(ref$path)), paste0("##contig=<ID=", ref$ref_id, ",length=", nchar(ref$full), ">"),
           "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"SNV, INS, DEL or REPEAT (the region between the reference flanks)\">",
           "##INFO=<ID=MINOR,Number=0,Type=Flag,Description=\"Not in any allele's consensus; carried by a minority of an allele's reads only\">",
           "##FILTER=<ID=LowSupport,Description=\"No allele's reads confirm the variant (consensus variant carried by under half of the allele's counted reads)\">",
           "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotype: allele 1 (shorter) | allele 2 (longer); 1 = the allele's consensus carries the ALT\">",
           "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">",
           paste0("##FORMAT=<ID=MC,Number=.,Type=String,Description=\"Repeat record: motif counts of each tract (", motifs_txt, ") per haplotype, tracts joined by _\">"),
           "##FORMAT=<ID=AL,Number=.,Type=Integer,Description=\"Repeat record: length in bp of each haplotype's sequence between the flanks\">",
           "##FORMAT=<ID=HDP,Number=.,Type=Integer,Description=\"Reads per haplotype (repeat record) or reads counted at the variant per haplotype\">",
           "##FORMAT=<ID=HAF,Number=.,Type=Float,Description=\"Share of each haplotype's counted reads that carry the ALT\">",
           paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples), collapse = "\t"))
  a12 <- alleles[alleles$allele %in% 1:2, ]
  body <- character(0)
  # the repeat record
  anchor_pos <- ref$n_start - 1L; anchor <- substr(ref$full, anchor_pos, anchor_pos)
  ref_block <- substr(ref$full, ref$n_start, ref$n_end)
  alt_seqs <- unique(a12$repeat_sequence[!is.na(a12$repeat_sequence)])
  if (length(alt_seqs)) {
    cells <- vapply(samples, function(sm) {
      x <- a12[a12$file_stem == sm, ]; x <- x[order(x$allele), ]
      if (!nrow(x)) return(".|.:.:.:.:.")
      if (nrow(x) == 1) x <- x[c(1, 1), ]
      gt <- paste(match(x$repeat_sequence, alt_seqs), collapse = "|")
      sprintf("%s:%d:%s:%s:%s", gt, anchor_pos, paste(x$motif_counts, collapse = ","), paste(nchar(x$repeat_sequence), collapse = ","), paste(x$reads, collapse = ","))
    }, "")
    body <- c(body, paste(c(ref$ref_id, anchor_pos, ".", paste0(anchor, ref_block), paste(paste0(anchor, alt_seqs), collapse = ","), ".", "PASS",
                            "TYPE=REPEAT", "GT:PS:MC:AL:HDP", cells), collapse = "\t"))
  }
  # flank variants (alleles 1 and 2; `variants` has a row for each sample, allele 1 and 2 and
  # variant found in either, with the read support measured in both)
  if (!is.null(variants) && nrow(variants)) variants <- variants[variants$allele %in% 1:2, , drop = FALSE]
  if (!is.null(variants) && nrow(variants)) {
    ev <- unique(variants[, c("pos", "ref", "alt", "type")]); ev <- ev[order(ev$pos), ]
    for (i in seq_len(nrow(ev))) {
      v <- variants[variants$pos == ev$pos[i] & variants$ref == ev$ref[i] & variants$alt == ev$alt[i], ]
      v$in_consensus[is.na(v$in_consensus)] <- FALSE
      minor_only <- !any(v$in_consensus)
      v$carried[is.na(v$carried)] <- FALSE; v$confirmed[is.na(v$confirmed)] <- FALSE
      cells <- vapply(samples, function(sm) {
        x <- a12[a12$file_stem == sm, ]; x <- x[order(x$allele), ]
        if (!nrow(x)) return(".|.:.:.:.")
        if (nrow(x) == 1) x <- x[c(1, 1), ]
        vv <- v[v$file_stem == sm, ]
        g <- vapply(x$allele, function(k) { r <- vv[vv$allele == k, ]; if (nrow(r) && r$in_consensus[1]) "1" else "0" }, "")
        dp <- vapply(x$allele, function(k) { r <- vv[vv$allele == k, ]; if (nrow(r)) as.character(r$alt_reads[1] + r$ref_reads[1]) else "." }, "")
        af <- vapply(x$allele, function(k) { r <- vv[vv$allele == k, ]; if (nrow(r) && !is.na(r$support[1])) sprintf("%.3f", r$support[1]) else "." }, "")
        sprintf("%s:%d:%s:%s", paste(g, collapse = "|"), anchor_pos, paste(dp, collapse = ","), paste(af, collapse = ","))
      }, "")
      filt <- if (any(v$carried & v$confirmed, na.rm = TRUE)) "PASS" else "LowSupport"
      body <- c(body, paste(c(ref$ref_id, ev$pos[i], ".", ev$ref[i], ev$alt[i], ".", filt,
                              paste0("TYPE=", ev$type[i], if (minor_only) ";MINOR" else ""), "GT:PS:HDP:HAF", cells), collapse = "\t"))
    }
  }
  # sort by position (the repeat record sits between the two flanks)
  if (length(body)) body <- body[order(as.integer(vapply(strsplit(body, "\t"), `[`, "", 2)))]
  writeLines(c(hdr, body), path)
  if (nzchar(Sys.which("bgzip")) && nzchar(Sys.which("tabix"))) {
    system2("bgzip", c("-f", "-k", path), stdout = FALSE, stderr = FALSE)
    system2("tabix", c("-f", "-p", "vcf", paste0(path, ".gz")), stdout = FALSE, stderr = FALSE)
  }
  invisible(path)
}

# Reads of one sample aligned to the whole reference (minimap2), with HP = allele
# (1 shorter, 2 longer, 3+ other clusters), XR = primary repeat length and, with UMIs,
# MI = molecule. In IGV, group or colour reads by the HP tag to see each allele's reads.
write_allele_bam <- function(d, ref, path, preset = "-x map-hifi", threads = 1L, work_dir = tempdir()) {
  if (!nrow(d)) return(invisible(NULL))
  stem <- file.path(work_dir, basename(sub("\\.bam$", "", path)))
  ref_fa <- paste0(stem, "-ref.fa"); reads_fa <- paste0(stem, "-reads.fa"); sam <- paste0(stem, ".sam")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(stats::setNames(ref$full, ref$ref_id)), ref_fa)
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(stats::setNames(paste0(d$left, d$mid, d$right), d$qname)), reads_fa)
  args <- c("-a", strsplit(preset, "\\s+")[[1]], "-t", threads, ref_fa, reads_fa)
  system2("minimap2", args, stdout = sam, stderr = FALSE)
  x <- readLines(sam); is_h <- startsWith(x, "@")
  rec <- x[!is_h]; q <- sub("\t.*$", "", rec)
  tag <- paste0("\tHP:i:", d$hp, "\tXR:i:", round(d$length), ifelse(is.na(d$molecule_id), "", paste0("\tMI:Z:", d$molecule_id)))
  rec <- paste0(rec, tag[match(q, d$qname)])
  writeLines(c(x[is_h], rec), sam)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  system2("samtools", c("sort", "-o", path, sam), stdout = FALSE, stderr = FALSE)
  system2("samtools", c("index", path), stdout = FALSE, stderr = FALSE)
  unlink(c(ref_fa, reads_fa, sam))
  invisible(path)
}

# One sample's alleles for review: one row per allele. Left and right panels: the stretch
# of each flank the allele's consensus covers (reference coordinates), with each variant
# marked (filled = in the consensus, open = carried by a minority of the reads) and
# labelled with the share of the allele's reads that carry it. Middle panel: the allele's
# repeat region written out in repeat notation.
allele_map_plot <- function(al, vars, ref, title) {
  al <- data.table::as.data.table(al)[order(allele)]
  al[, row := factor(sprintf("Allele %s\n%s reads (%s)", allele, format(reads, big.mark = ","), scales::percent(proportion, 1)),
                     levels = rev(sprintf("Allele %s\n%s reads (%s)", allele, format(reads, big.mark = ","), scales::percent(proportion, 1))))]
  v <- data.table::as.data.table(vars)
  if (nrow(v)) v <- al[, .(file_stem, allele, row)][v, on = .(file_stem, allele), nomatch = NULL]
  side_plot <- function(side) {
    ns1 <- ref$n_start; ne1 <- ref$n_end
    seg <- if (side == "left") al[, .(row, from = left_from, to = ns1 - 1)] else al[, .(row, from = ne1 + 1, to = right_to)]
    seg <- seg[!is.na(from) & !is.na(to) & to >= from]
    ns0 <- ref$n_start; ne0 <- ref$n_end   # `ref` is also a column of the variants table
    vv <- if (nrow(v)) v[if (side == "left") pos < ns0 else pos > ne0] else v
    p <- ggplot2::ggplot(seg, ggplot2::aes(y = row)) + ggplot2::geom_segment(ggplot2::aes(x = from, xend = to, yend = row), linewidth = 3, colour = "grey80")
    # labels alternate above and below the line, so neighbouring variants do not overlap
    if (nrow(vv)) { vv <- vv[order(row, pos)][, lab_v := ifelse(seq_len(.N) %% 2 == 1, -0.7, 1.7), by = row]
                    vv[, alt_lab := ifelse(nchar(alt) > 8, paste0(substr(alt, 1, 6), "..."), alt)][, ref_lab := ifelse(nchar(ref) > 8, paste0(substr(ref, 1, 6), "..."), ref)] }
    if (nrow(vv)) p <- p + ggplot2::geom_point(data = vv, ggplot2::aes(pos, row, shape = in_consensus), size = 3, colour = "#d95f02", fill = "#d95f02") +
      ggplot2::geom_text(data = vv, ggplot2::aes(pos, row, vjust = lab_v, label = sprintf("%d %s>%s\n%s", pos, ref_lab, alt_lab, ifelse(is.na(support), "support n/a", scales::percent(support, 1)))),
                         size = 2.4, lineheight = 0.85) +
      ggplot2::scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1), labels = c(`TRUE` = "in the consensus", `FALSE` = "minority of reads"), name = NULL, drop = FALSE)
    p + ggplot2::scale_x_continuous(labels = scales::comma, breaks = scales::pretty_breaks(4)) +
      ggplot2::labs(x = sprintf("%s flank (reference position)", if (side == "left") "Left" else "Right"), y = NULL) +
      ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
  }
  mid <- ggplot2::ggplot(al, ggplot2::aes(0, row)) +
    ggplot2::geom_label(ggplot2::aes(label = stringr_wrap(repeat_notation, 34)), size = 2.8, label.size = 0, fill = "grey95", lineheight = 0.9) +
    ggplot2::scale_x_continuous(limits = c(-1, 1)) + ggplot2::labs(x = "Between the flanks", y = NULL) +
    ggplot2::theme_bw() + ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
  patchwork::wrap_plots(side_plot("left"), mid + ggplot2::theme(axis.text.y = ggplot2::element_blank()), side_plot("right") + ggplot2::theme(axis.text.y = ggplot2::element_blank()),
                        nrow = 1, widths = c(2, 1.6, 2), guides = "collect") +
    patchwork::plot_annotation(title = title, subtitle = "Grey = flank sequence each allele's consensus covers; orange = variants against the reference, labelled with the share of the allele's reads carrying them") &
    ggplot2::theme(legend.position = "bottom")
}
stringr_wrap <- function(x, width) vapply(x, function(s) if (is.na(s)) "" else paste(strwrap(s, width), collapse = "\n"), "")

write_vcf <- function(variants, columns, ref_name, ref_length, ref_file, path) {
  hdr <- c("##fileformat=VCFv4.2", paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
           "##source=Duke", paste0("##reference=", basename(ref_file)),
           paste0("##contig=<ID=", ref_name, ",length=", ref_length, ">"),
           "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"SNV, INS or DEL\">",
           "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Haploid genotype: 1 = this cluster's consensus carries the variant\">",
           paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", columns), collapse = "\t"))
  if (is.null(variants) || nrow(variants) == 0) { writeLines(hdr, path); return(invisible(path)) }
  ev <- unique(variants[, c("pos", "ref", "alt", "type")])
  ev <- ev[order(ev$pos, ev$ref, ev$alt), ]
  body <- vapply(seq_len(nrow(ev)), function(i) {
    has <- variants$column[variants$pos == ev$pos[i] & variants$ref == ev$ref[i] & variants$alt == ev$alt[i]]
    gt <- ifelse(columns %in% has, "1", "0")
    paste(c(ref_name, ev$pos[i], ".", ev$ref[i], ev$alt[i], ".", "PASS", paste0("TYPE=", ev$type[i]), "GT", gt), collapse = "\t")
  }, "")
  writeLines(c(hdr, body), path)
  invisible(path)
}

# Y-axis labels for a waterfall. Rows are reads, so a common length takes up many
# rows and a rare one few. Labels are placed at the first read reaching a length,
# kept at least 1/30 of the plot apart, choosing multiples of 10 first, then of 5,
# then any whole length, so the axis is labelled evenly up the plot (including
# inside a large allele) rather than only at a few lengths.
waterfall_breaks <- function(len, n_target = 30) {
  n <- length(len); if (!n) return(integer(0))
  gap <- n / n_target
  cand <- sort(unique(floor(len)))
  row_of <- function(v) which(len >= v)[1]
  chosen <- c(1L, n)
  for (step in c(10, 5, 1)) {
    for (v in cand[cand %% step == 0]) {
      r <- row_of(v); if (is.na(r)) next
      if (all(abs(chosen - r) >= gap)) chosen <- c(chosen, r)
    }
  }
  sort(unique(chosen))
}

# One sample's clusters: histogram of the primary length, stacked by cluster, with
# each cluster's modal length and share of reads written above it.
cluster_length_plot <- function(d, cs, title, xlab = "Repeat length") {
  d <- data.table::as.data.table(d)[is.finite(length)]
  d[, cluster := factor(cluster_number, levels = sort(unique(cluster_number)))]
  cs <- data.table::as.data.table(cs)[, cluster := factor(cluster_number, levels = levels(d$cluster))]
  h <- d[, .N, by = .(cluster, length = floor(length + 0.5))]
  top <- h[, .(top = sum(N)), by = length][, max(top)]
  cs[, label := sprintf("C%s: mode %s, %s of reads%s", cluster, round(modal_length), scales::percent(proportion, 0.1),
                        ifelse(!is.na(haplotype_cluster) & haplotype_cluster != "H1", paste0(", ", haplotype_cluster), ""))]
  xr <- c(min(d$length), stats::quantile(d$length, 0.999)) + c(-3, 3)
  ggplot2::ggplot(h, ggplot2::aes(length, N, fill = cluster)) + ggplot2::geom_col(width = 1) +
    ggplot2::geom_vline(data = cs, ggplot2::aes(xintercept = modal_length, colour = cluster), linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
    ggplot2::scale_fill_brewer(palette = "Set1", name = "Cluster", labels = stats::setNames(cs$label, cs$cluster)) +
    ggplot2::scale_colour_brewer(palette = "Set1") +
    ggplot2::scale_y_sqrt(labels = scales::comma, expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::coord_cartesian(xlim = xr) +
    ggplot2::labs(x = xlab, y = "Reads (square-root scale)", title = title,
                  subtitle = "Bars stacked by cluster; dashed line = each cluster's modal length") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "right")
}

# ------------------------------------------------------------------------------
# Waterfall: reads stacked by repeat length, bases coloured. Only the
# flank_bp bases either side of the tract are drawn - the informative region,
# and far quicker than drawing whole 3 kb reads.
# ------------------------------------------------------------------------------
create_waterfall_plot <- function(d, title, flank_bp = 60L, max_reads = 1000L, length_label = "Repeat length") {
  d <- d[!is.na(d$mid), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  n_total <- nrow(d)
  if (!is.na(max_reads) && nrow(d) > max_reads) { set.seed(123); d <- d[sample(nrow(d), max_reads), , drop = FALSE] }
  d <- d[order(d$length, d$mid), , drop = FALSE]
  d$read_order <- seq_len(nrow(d))
  lft <- ifelse(nchar(d$left) > flank_bp, substr(d$left, nchar(d$left) - flank_bp + 1, nchar(d$left)), d$left)
  rgt <- substr(d$right, 1, flank_bp)
  mk <- function(s, offset_fun, seg) {
    ch <- strsplit(s, "")
    n <- lengths(ch)
    data.frame(read_order = rep(d$read_order, n), pos = unlist(Map(offset_fun, n, seq_along(n))),
               base = unlist(ch), segment = seg)
  }
  pd <- rbind(mk(lft, function(n, i) seq_len(n) - n - 1L, "left"),
              mk(d$mid, function(n, i) seq_len(n) - 1L, "tract"),
              mk(rgt, function(n, i) nchar(d$mid[i]) + seq_len(n) - 1L, "right"))
  cols <- c(A = "#1b9e77", C = "#377eb8", G = "grey15", T = "#e41a1c", N = "grey70", `-` = "grey85")
  # Mark the first read at or above each of a set of evenly spaced lengths, so the
  # labels are distinct lengths rather than read quantiles (which repeat the
  # length of a dominant allele many times)
  brk <- waterfall_breaks(d$length)
  ggplot2::ggplot(pd, ggplot2::aes(pos, read_order, fill = base)) +
    ggplot2::geom_raster() +
    ggplot2::geom_vline(xintercept = -0.5, linetype = "dashed", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = cols, name = "Base") +
    ggplot2::scale_y_continuous(breaks = brk, labels = round(d$length[brk]), expand = c(0, 0)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(title = title,
                  subtitle = sprintf("%s of %s reads, ordered by repeat length; %d bp of flank shown either side; dashed line = tract start",
                                     format(nrow(d), big.mark = ","), format(n_total, big.mark = ","), flank_bp),
                  x = "Position relative to tract start (bp)", y = length_label) +
    ggplot2::theme_bw()
}
