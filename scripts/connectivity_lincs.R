# connectivity_lincs.R
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------



#'
#' Connectivity Score
#'
#' Function implementing an audit / analysis step of the scDrugAudit framework.
#'
#' @param query Function argument (see Details).
#' @param ref Function argument (see Details).
#' @param method Function argument (see Details).
#' @return audit / result object.
#' @export

connectivity_score <- function(query, ref, method = c("ks", "cosine")) {
  method <- match.arg(method)
  g <- intersect(names(query), names(ref))
  if (length(g) < 5) return(NA)
  if (method == "cosine") {
    q <- query[g]; r <- ref[g]
    sum(q * r) / (sqrt(sum(q^2)) * sqrt(sum(r^2)))
  } else {
    ord <- order(ref[g], decreasing = TRUE)
    rnk <- g[ord]; q <- query[rnk]
    w <- abs(q) / sum(abs(q))
    cum <- cumsum(w * sign(q))
    up <- which(q > 0); down <- which(q < 0)
    es_up <- if (length(up)) max(cum[up]) - min(cum[up]) else 0
    es_down <- if (length(down)) max(-cum[down]) - min(-cum[down]) else 0
    es_up - es_down
  }
}



#'
#' Calibrate Connectivity
#'
#' Function implementing an audit / analysis step of the scDrugAudit framework.
#'
#' @param query Function argument (see Details).
#' @param ref_db Function argument (see Details).
#' @param n_perm Function argument (see Details).
#' @param method Function argument (see Details).
#' @param seed Function argument (see Details).
#' @return audit / result object.
#' @export

calibrate_connectivity <- function(query, ref_db, n_perm = 200, method = "ks", seed = 1) {
  ref_db <- as.matrix(ref_db)
  if (is.null(colnames(ref_db))) colnames(ref_db) <- paste0("C", seq_len(ncol(ref_db)))
  scores <- sapply(seq_len(ncol(ref_db)), function(j)
    connectivity_score(query, ref_db[, j], method))
  names(scores) <- colnames(ref_db)
  gq <- names(query)
  qvals <- as.numeric(query)
  set.seed(seed)
  null <- matrix(NA, n_perm, ncol(ref_db))
  for (b in seq_len(n_perm)) {
    pq <- qvals[sample(length(qvals))]; names(pq) <- gq
    null[b, ] <- sapply(seq_len(ncol(ref_db)), function(j)
      connectivity_score(pq, ref_db[, j], method))
  }
  p <- sapply(seq_len(ncol(ref_db)), function(j)
    (sum(abs(null[, j]) >= abs(scores[j])) + 1) / (n_perm + 1))
  q <- p.adjust(p, method = "BH")
  data.frame(compound = colnames(ref_db), score = round(scores, 3),
             p_value = round(p, 4), q_value = round(q, 4),
             significant = q < 0.05, stringsAsFactors = FALSE)
}



#'
#' Audit Connectivity
#'
#' Function implementing an audit / analysis step of the scDrugAudit framework.
#'
#' @param query_naive Function argument (see Details).
#' @param query_rigorous Function argument (see Details).
#' @param ref_db Function argument (see Details).
#' @param n_perm Function argument (see Details).
#' @param method Function argument (see Details).
#' @param seed Function argument (see Details).
#' @return audit / result object.
#' @export

audit_connectivity <- function(query_naive, query_rigorous, ref_db,
                               n_perm = 200, method = "ks", seed = 1) {
  cal_n <- calibrate_connectivity(query_naive, ref_db, n_perm, method, seed)
  cal_r <- calibrate_connectivity(query_rigorous, ref_db, n_perm, method, seed + 1)
  top_n <- cal_n$compound[which.max(abs(cal_n$score))]
  top_r <- cal_r$compound[which.max(abs(cal_r$score))]
  list(
    naive = data.frame(
      n_compounds_tested = ncol(ref_db),
      n_significant = sum(cal_n$significant),
      top_compound = top_n,
      top_score = round(max(abs(cal_n$score)), 3),
      top_significant = cal_n$significant[which.max(abs(cal_n$score))],
      stringsAsFactors = FALSE),
    rigorous = data.frame(
      n_compounds_tested = ncol(ref_db),
      n_significant = sum(cal_r$significant),
      top_compound = top_r,
      top_score = round(max(abs(cal_r$score)), 3),
      top_significant = cal_r$significant[which.max(abs(cal_r$score))],
      stringsAsFactors = FALSE),
    calibration_naive = cal_n,
    calibration_rigorous = cal_r
  )
}



#'
#' Load Lincs L1000 Ref
#'
#' Function implementing an audit / analysis step of the scDrugAudit framework.
#'
#' @param path Function argument (see Details).
#' @param url Function argument (see Details).
#' @param gene_id Function argument (see Details).
#' @param query_genes Function argument (see Details).
#' @return audit / result object.
#' @export

load_lincs_l1000_ref <- function(path = NULL, url = NULL,
                                 gene_id = c("entrez", "symbol"),
                                 query_genes = NULL) {
  gene_id <- match.arg(gene_id)
  if (is.null(path) && is.null(url))
    stop("provide a local path or remote url for the LINCS L1000 reference; ",
         "with no data, use make_synthetic_lincs() as a fallback.")
  mat <- if (!is.null(path) && file.exists(path)) {
    readRDS(path)
  } else if (!is.null(url)) {
    if (!requireNamespace("httr", quietly = TRUE) &&
        !requireNamespace("curl", quietly = TRUE))
      stop("remote loading needs httr or curl; or pre-download locally first.")
    tmp <- tempfile(fileext = ".rds")
    if (requireNamespace("httr", quietly = TRUE))
      httr::GET(url, httr::write_disk(tmp, overwrite = TRUE)) else
      curl::curl_download(url, tmp)
    readRDS(tmp)
  } else stop("path does not exist and no usable url was provided.")
  if (!is.matrix(mat) && !inherits(mat, "dgCMatrix"))
    stop("reference must be a genes x compounds matrix.")
  if (!is.null(query_genes) && gene_id == "entrez") {
    if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
        all(grepl("^[0-9]+$", rownames(mat)))) {
      map <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
        keys = rownames(mat), column = "SYMBOL", keytype = "ENTREZID")
      rownames(mat) <- as.character(map[rownames(mat)])
      # drop unmapped / empty and de-duplicate symbols (multi-Entrez -> one SYMBOL)
      keep <- !is.na(rownames(mat)) & nzchar(rownames(mat))
      mat <- mat[keep, , drop = FALSE]
      dup <- duplicated(rownames(mat))
      if (any(dup)) mat <- mat[!dup, , drop = FALSE]
    }
  }
  mat
}



#'
#' Make Synthetic Lincs
#'
#' Function implementing an audit / analysis step of the scDrugAudit framework.
#'
#' @param n_genes Function argument (see Details).
#' @param n_compounds Function argument (see Details).
#' @param seed Function argument (see Details).
#' @return audit / result object.
#' @export

make_synthetic_lincs <- function(n_genes = 200, n_compounds = 60, seed = 20261024) {
  set.seed(seed)
  genes <- paste0("G", seq_len(n_genes))
  mat <- matrix(rnorm(n_genes * n_compounds), nrow = n_genes, ncol = n_compounds,
                dimnames = list(genes, paste0("C", seq_len(n_compounds))))
  list(mat = mat, genes = genes, compounds = colnames(mat))
}



#'
#' Print Connectivity
#'
#' Function implementing an audit / analysis step of the scDrugAudit framework.
#'
#' @param aud Function argument (see Details).
#' @return audit / result object.
#' @export

print_connectivity <- function(aud) {
  cat("\n===== Audit #5: connectivity calibration (FDR-inflated hits) =====\n")
  cat("naive signature: testing", aud$naive$n_compounds_tested, " compounds, ",
      aud$naive$n_significant, " significant at FDR<0.05; ",
      "top =", aud$naive$top_compound, "(score", aud$naive$top_score,
      ", significant:", aud$naive$top_significant, ")\n")
  cat("rigorous signature: testing", aud$rigorous$n_compounds_tested, " compounds, ",
      aud$rigorous$n_significant, " significant at FDR<0.05; ",
      "top =", aud$rigorous$top_compound, "(score", aud$rigorous$top_score,
      ", significant:", aud$rigorous$top_significant, ")\n")
  cat("if naive top hit significant=FALSE, this is the typical paper trap of a high raw score without FDR support (false hit).\n\n")
  invisible(aud)
}
