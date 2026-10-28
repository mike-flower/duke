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
  switch(region, left = left, right = right, both = paste0(left, right))
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
  unique(data.table::data.table(r = r, kmer = substring(s[r], st, st + k - 1L)))
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
  ok <- which(!is.na(seqs) & nchar(seqs) >= k & !grepl("[^ACGT]", seqs))
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
consensus_function <- function(seqs, threshold = 0.5, sample_size = 50L) {
  seqs <- seqs[!is.na(seqs) & nchar(seqs) > 0]
  if (length(seqs) == 0) return(list(consensus = NA_character_, n = 0L))
  if (!is.na(sample_size) && length(seqs) > sample_size) { set.seed(123); seqs <- sample(seqs, sample_size) }
  if (length(seqs) == 1) return(list(consensus = seqs, n = 1L))
  dss <- Biostrings::DNAStringSet(seqs)
  aln <- suppressMessages(DECIPHER::AlignSeqs(dss, processors = 1, verbose = FALSE))
  cons <- suppressMessages(DECIPHER::ConsensusSequence(aln, threshold = threshold, ambiguity = FALSE,
                                                       noConsensusChar = "N", includeTerminalGaps = FALSE))
  cons <- gsub("-", "", as.character(cons))
  list(consensus = gsub("[^ACGT]", "N", cons), n = length(seqs))
}

# ------------------------------------------------------------------------------
# Variants in flank consensus sequences, in coordinates of the full reference
# ------------------------------------------------------------------------------
# Each flank consensus is aligned to its own reference flank (the consensus
# globally, the flank locally), so the repeat never enters the alignment and
# positions map directly onto the reference. Runs of adjacent inserted or
# deleted bases are merged into single events. Positions are 1-based on the
# full reference including its N block, ready for IGV alongside that FASTA.
call_flank_variants <- function(consensus, ref_flank, offset, edge = 3L) {
  if (is.na(consensus) || nchar(consensus) < 10) return(NULL)
  ns <- if (requireNamespace("pwalign", quietly = TRUE)) "pwalign" else "Biostrings"
  pw <- lapply(setNames(nm = c("pairwiseAlignment", "alignedPattern", "alignedSubject", "subject")),
               function(f) getExportedValue(ns, f))
  aln <- pw$pairwiseAlignment(Biostrings::DNAString(consensus), Biostrings::DNAString(ref_flank),
                              type = "global-local", gapOpening = 5, gapExtension = 2)
  q <- strsplit(as.character(pw$alignedPattern(aln)), "")[[1]]
  r <- strsplit(as.character(pw$alignedSubject(aln)), "")[[1]]
  rpos <- BiocGenerics::start(pw$subject(aln)) - 1L + cumsum(r != "-")
  ev <- list(); k <- 1L; i <- 1L; n <- length(q)
  while (i <= n) {
    if (r[i] != "-" && q[i] != "-" && r[i] != q[i] && q[i] != "N" && r[i] != "N") {
      ev[[k]] <- data.frame(pos = rpos[i], ref = r[i], alt = q[i], type = "SNV"); k <- k + 1L; i <- i + 1L
    } else if (q[i] == "-") {
      j <- i; while (j < n && q[j + 1] == "-") j <- j + 1L
      anchor_pos <- rpos[i] - 1L
      anchor <- if (anchor_pos >= 1) substr(ref_flank, anchor_pos, anchor_pos) else ""
      ev[[k]] <- data.frame(pos = anchor_pos, ref = paste0(anchor, paste(r[i:j], collapse = "")), alt = anchor, type = "DEL"); k <- k + 1L
      i <- j + 1L
    } else if (r[i] == "-") {
      j <- i; while (j < n && r[j + 1] == "-") j <- j + 1L
      anchor_pos <- if (i > 1) rpos[i - 1] else 0L
      anchor <- if (anchor_pos >= 1) substr(ref_flank, anchor_pos, anchor_pos) else ""
      ev[[k]] <- data.frame(pos = anchor_pos, ref = anchor, alt = paste0(anchor, paste(q[i:j], collapse = "")), type = "INS"); k <- k + 1L
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
  brk <- unique(round(seq(1, nrow(d), length.out = min(12, nrow(d)))))
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
