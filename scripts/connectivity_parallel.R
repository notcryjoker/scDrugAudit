# Parallel, NUMERICALLY EQUIVALENT KS calibration for the full EH3226 L1000 ref.
#
# This is calibrate_connectivity_fast() split across columns:
#   * the permutation stream V is generated ONCE on the master with exactly the
#     same set.seed(seed) + sample() sequence as both the formal runner and
#     connectivity_fast.R;
#   * each worker receives ONLY its own block of reference columns;
#   * the master assembles p = (count+1)/(n_perm+1), BH q, and the same
#     data.frame schema as calibrate_connectivity().
#
# MEMORY / SPEED: .kc_worker is a TOP-LEVEL function (environment = globalenv),
# so R does not serialize the surrounding frame that holds the ~1.1 GB ref_db.
# All column blocks are handed to parLapply IN ONE CALL so every worker stays
# busy; the master temporarily holds one extra copy of the matrix as blocks.
#
# Validation: scripts/_verify_parallel.R compares this against
# 11_lincs_connectivity.R --fast on the 100-signature reference for all three
# signature specs (primary / mean / min-ct 2) -> 6/6 identical.
#
# Enable in 11_lincs_connectivity.R with --parallel N [--block-cols M].

.kc_worker <- function(M, V, g, gq, S) {
  suppressPackageStartupMessages(library(matrixStats))
  np1 <- nrow(V); k <- ncol(M)
  sc <- matrix(NA_real_, np1, k)
  for (j in seq_len(k)) {
    idx <- match(g[order(M[g, j], decreasing = TRUE)], gq)
    Q  <- V[, idx, drop = FALSE]
    C  <- matrixStats::rowCumsums(Q) / S
    Cu <- C; Cu[Q <= 0] <- NA_real_
    Cd <- C; Cd[Q >= 0] <- NA_real_
    # Algebraic identity, bit-for-bit safe: for the down-gene positions,
    #   max(-C) - min(-C) == max(C) - min(C)
    # because negation is exact in IEEE-754 and b-a == (-a)-(-b).
    # rowRanges() does one pass instead of separate rowMaxs + rowMins.
    ru <- matrixStats::rowRanges(Cu, na.rm = TRUE)
    rd <- matrixStats::rowRanges(Cd, na.rm = TRUE)
    sc[, j] <- (ru[, 2L] - ru[, 1L]) - (rd[, 2L] - rd[, 1L])
  }
  obs <- sc[1L, ]
  cnt <- vapply(seq_len(k),
                function(j) sum(abs(sc[-1L, j]) >= abs(obs[j])), numeric(1L))
  list(score = obs, count = cnt)
}

calibrate_connectivity_parallel <- function(query, ref_db, n_perm = 2000,
                                            method = "ks", seed = 1,
                                            n_workers = 4L, block_cols = 4000L) {
  if (method != "ks")
    return(.calibrate_connectivity_ref(query, ref_db, n_perm, method, seed))
  if (!requireNamespace("parallel", quietly = TRUE) ||
      !requireNamespace("matrixStats", quietly = TRUE))
    stop("parallel calibration needs parallel + matrixStats")

  ref_db <- as.matrix(ref_db)
  if (is.null(colnames(ref_db))) colnames(ref_db) <- paste0("C", seq_len(ncol(ref_db)))
  gq0 <- names(query); qv <- as.numeric(query)
  gq <- gq0
  g  <- intersect(gq, rownames(ref_db))
  if (length(g) < 5L) stop("fewer than 5 query genes align with the reference")
  S <- sum(abs(qv[match(g, gq)]))
  nc <- ncol(ref_db)

  # EXACT same permutation stream as calibrate_connectivity_fast()
  set.seed(seed)
  V <- matrix(0, n_perm + 1L, length(qv))
  V[1L, ] <- qv
  for (b in seq_len(n_perm)) V[b + 1L, ] <- qv[sample(length(qv))]

  # Permutation columns outside the ranked gene set g are never used. Drop
  # them AFTER drawing the permutation stream on the FULL query: this is a
  # pure projection, Vfull[, match(g[ord], gq)] == Vsmall[, ord].
  V <- V[, match(g, gq), drop = FALSE]
  gq <- g

  blocks <- split(seq_len(nc), ceiling(seq_len(nc) / block_cols))
  nb <- length(blocks)
  nw <- min(as.integer(n_workers), nb)

  # Restartable checkpoints: one small rds per column block, keyed by the
  # query vector, permutation settings and reference column names. A killed
  # run resumes without recomputing finished blocks. Purely additive: with
  # no --ckpt-dir the numerical path is unchanged.
  ckpt <- getOption("scda.ckpt", NULL)
  ckfiles <- NULL
  if (!is.null(ckpt)) {
    dir.create(ckpt, showWarnings = FALSE, recursive = TRUE)
    key <- substr(digest::digest(list(ver = "v2", qv = qv, qnames = gq0,
                                      n_perm = n_perm, seed = seed,
                                      cn = colnames(ref_db),
                                      block_cols = block_cols), algo = "md5"),
                  1L, 12L)
    ckfiles <- file.path(ckpt, sprintf("%s_b%03d.rds", key, seq_len(nb)))
    cat(sprintf("[par] checkpoint key %s | %d blocks in %s\n", key, nb, ckpt))
  }

  cl <- parallel::makeCluster(nw)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # Two blocks per worker per batch: all workers stay busy, while the master
  # holds only 2*nw column blocks at a time instead of a full matrix copy.
  # .kc_worker is top-level so no big closure is serialized to workers.
  step <- max(nw, 1L) * 2L
  out <- vector("list", nb)
  for (s in seq.int(1L, nb, by = step)) {
    ids <- s:min(s + step - 1L, nb)
    todo <- ids
    if (!is.null(ckpt)) {
      cached <- vapply(ids, function(i) {
        if (!file.exists(ckfiles[i])) return(FALSE)
        z <- tryCatch(readRDS(ckfiles[i]), error = function(e) NULL)
        !is.null(z) && length(z$score) == length(blocks[[i]])
      }, logical(1L))
      for (i in ids[cached]) out[[i]] <- readRDS(ckfiles[i])
      todo <- ids[!cached]
    }
    if (length(todo)) {
      mats <- lapply(blocks[todo], function(cols) ref_db[, cols, drop = FALSE])
      z <- parallel::parLapply(cl, mats, .kc_worker,
                               V = V, g = g, gq = gq, S = S)
      out[todo] <- z
      if (!is.null(ckpt))
        for (i in todo) saveRDS(out[[i]], ckfiles[i])
      rm(mats, z); invisible(gc())
    }
    cat(sprintf("[par] %d/%d column blocks done%s\n", max(ids), nb,
                if (length(todo) < length(ids))
                  sprintf(" (%d from checkpoint)", length(ids) - length(todo))
                else ""))
  }
  rm(blocks)

  scores <- unlist(lapply(out, `[[`, "score"), use.names = FALSE)
  counts <- unlist(lapply(out, `[[`, "count"), use.names = FALSE)
  names(scores) <- colnames(ref_db)
  p <- (counts + 1) / (n_perm + 1)
  q <- p.adjust(p, method = "BH")
  data.frame(compound = colnames(ref_db), score = round(scores, 3),
             p_value = round(p, 4), q_value = round(q, 4),
             significant = q < 0.05, stringsAsFactors = FALSE)
}
