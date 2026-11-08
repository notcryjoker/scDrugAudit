# Numerically EQUIVALENT fast KS calibration for audit_connectivity().
#
# Why this is equivalent to connectivity_lincs.R::calibrate_connectivity():
#   connectivity_score(query, ref, "ks") does
#     g   <- intersect(names(query), names(ref))          # fixed gene set
#     rnk <- g[order(ref[g], decreasing = TRUE)]          # depends ONLY on ref
#     q   <- query[rnk]                                   # permutation reorders q
#     w   <- abs(q) / sum(abs(q))
#     cum <- cumsum(w * sign(q))
#   Since w * sign(q) == q / sum(abs(q)), and sum(abs(q)) is invariant under
#   permutation of the query vector, the whole score reduces to
#     cum <- cumsum(q) / S ,  S = sum(abs(query[g]))
#   with es_up / es_down computed on q>0 / q<0 positions.
#   The per-column ranking index is therefore computed ONCE and reused for all
#   permutations, cutting ~7.1M order() calls down to ncol(ref_db) calls.
#
# Validation: scripts/_verify_fast.R compares every column of the returned
# data.frame against the original implementation on a small matrix.
#
# Usage in 11_lincs_connectivity.R: pass --fast to source this file and shadow
# calibrate_connectivity; audit_connectivity() then picks up the fast version.

.calibrate_connectivity_ref <- calibrate_connectivity

calibrate_connectivity_fast <- function(query, ref_db, n_perm = 200,
                                        method = "ks", seed = 1) {
  if (method != "ks")
    return(.calibrate_connectivity_ref(query, ref_db, n_perm, method, seed))

  ref_db <- as.matrix(ref_db)
  if (is.null(colnames(ref_db)))
    colnames(ref_db) <- paste0("C", seq_len(ncol(ref_db)))

  gq <- names(query)
  qv <- as.numeric(query)
  g  <- intersect(gq, rownames(ref_db))
  if (length(g) < 5L) stop("fewer than 5 query genes align with the reference")

  S <- sum(abs(qv[match(g, gq)]))
  nc <- ncol(ref_db)
  if (!requireNamespace("matrixStats", quietly = TRUE))
    stop("fast KS calibration needs matrixStats (rowCumsums/rowMaxs/rowMins)")

  # Reproduce the ORIGINAL permutation stream exactly: set.seed(seed), then
  # draw sample(length(qv)) once per permutation. Row 1 = observed signature.
  set.seed(seed)
  V <- matrix(0, n_perm + 1L, length(qv))
  V[1L, ] <- qv
  for (b in seq_len(n_perm)) V[b + 1L, ] <- qv[sample(length(qv))]

  # For each compound, rank genes ONCE and evaluate all permutations together
  # with C-level row operations. ~nc R-level iterations instead of ~nc*n_perm.
  sc <- matrix(NA_real_, n_perm + 1L, nc)
  for (j in seq_len(nc)) {
    idx <- match(g[order(ref_db[g, j], decreasing = TRUE)], gq)
    Q  <- V[, idx, drop = FALSE]
    C  <- matrixStats::rowCumsums(Q) / S
    Cu <- C;  Cu[Q <= 0] <- NA_real_      # query-up positions only
    Cd <- -C; Cd[Q >= 0] <- NA_real_     # -(cum) at query-down positions
    eu <- matrixStats::rowMaxs(Cu, na.rm = TRUE) -
          matrixStats::rowMins(Cu, na.rm = TRUE)
    ed <- matrixStats::rowMaxs(Cd, na.rm = TRUE) -
          matrixStats::rowMins(Cd, na.rm = TRUE)
    sc[, j] <- eu - ed
  }

  scores <- sc[1L, ]
  names(scores) <- colnames(ref_db)
  null <- sc[-1L, , drop = FALSE]

  p <- vapply(seq_len(nc), function(j)
    (sum(abs(null[, j]) >= abs(scores[j])) + 1) / (n_perm + 1), numeric(1L))
  q <- p.adjust(p, method = "BH")
  data.frame(compound = colnames(ref_db), score = round(scores, 3),
             p_value = round(p, 4), q_value = round(q, 4),
             significant = q < 0.05, stringsAsFactors = FALSE)
}
