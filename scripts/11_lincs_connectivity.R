#!/usr/bin/env Rscript
# 11_lincs_connectivity.R
# --------------------------------------------------------------------------
# Battery #5 (LINCS L1000 connectivity calibration) -- local kit runner.
#
# Aligned to the FINALIZED formal spec (library doc
# "scDrugAudit · 电池#5 连接性校准 口径定稿.md") and to the formal runner
# scripts/11_lincs_connectivity.R in sc_drug_bias_benchmark/.
#
# ===== FINALIZED SIGNATURE SPEC (see PROJECT_STATE.md) =====
#   * --query shared (DEFAULT, formal semantics):
#       BOTH queries use the SAME gene backbone (audit$candidates, 5557 genes),
#       but each query takes its signed logFC FROM ITS OWN pipeline
#       ($naive / $rigorous). Controlled comparison: identical gene set,
#       different effect sizes. Requires $naive and $rigorous top-level slots
#       (per-celltype DE table LISTs) in the audit cache.
#       --query identical  = both queries are the SAME vector (diagnostic only,
#                            triggers a BY CONSTRUCTION warning).
#       --query per-pipeline = each pipeline's own significant-gene signature
#                            (naive 606 / rigorous 5557); legacy sensitivity.
#       --query raw = local-kit fallback for a MERGED cache that has only
#                            $candidates; NOT formal semantics (warned loudly).
#   * --agg absmax (DEFAULT, PRIMARY): keep the logFC with largest |logFC| and
#     its sign. --agg mean is the aggregation SENSITIVITY; the difference
#     between them IS the fragility finding for the manuscript and does not
#     depend on any exact cross-CT conflict percentage (86.5% naive / 85.8%
#     rigorous under definition A in the formal cache).
#     EXCLUDED PERMANENTLY: n_donors weighting -- the naive pipeline's
#     data.frame has NO n_donors COLUMN (structural: naive builds no
#     donor-level model); rigorous has one (4-9). Never reintroduce.
#   * --min-ct DEFAULT OFF (0): full 5557-gene backbone. --min-ct 2 is the
#     sensitivity run (921 genes with n_ct>=2). n_ct>=2 is the L1-L4 evidence
#     grading rule, NOT a query filter in the primary run.
#   * --nperm DEFAULT 2000. This does NOT lift the resolution floor: for real
#     L1000 m~12328 the BH rank-1 threshold is ~4e-6, while min attainable p
#     is 5e-3 (200) or 5e-4 (2000). Both are RESOLUTION-LIMITED, so "0
#     significant" is a constructive necessity, NOT a biological finding.
#     2000 only refines the p step, smooths the FDR curve, stabilises ranking.
#
# Usage (formal repo root):
#   Rscript scripts/11_lincs_connectivity.R --ref $REF --query shared --agg absmax \
#     --nperm 2000 --gene-id entrez --outdir results/tables/lincs_primary
#
# Extra switches in this kit:
#   --transpose        ref is compounds x genes (formal runner needs a value;
#                      here it is a plain flag)
#   --fast             numerically-equivalent vectorised KS calibration,
#                      STRONGLY recommended when ref has >1000 compounds
#   --audit PATH       default results/tables/full_audit_real.rds
#   --min-celltypes N  alias of --min-ct ; --n-perm N alias of --nperm
# --------------------------------------------------------------------------

args0 <- commandArgs(FALSE)
fidx <- grep("^--file=", args0)
script_dir <- if (length(fidx)) dirname(sub("^--file=", "", args0[fidx[1]])) else "."
# formal layout has R/connectivity_lincs.R; this kit keeps it in scripts/
cand_paths <- c(file.path(dirname(script_dir), "R", "connectivity_lincs.R"),
                file.path(script_dir, "connectivity_lincs.R"))
src <- cand_paths[file.exists(cand_paths)][1]
if (is.na(src)) stop("cannot locate connectivity_lincs.R (R/ or scripts/)")
source(src)

args <- commandArgs(TRUE)
get_arg <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (!is.na(i) && length(args) >= i + 1) args[i + 1] else default
}
has <- function(name) paste0("--", name) %in% args

ref_path   <- get_arg("ref")
audit_path <- get_arg("audit", "results/tables/full_audit_real.rds")
outdir     <- get_arg("outdir", ".")
gene_id    <- get_arg("gene-id", "symbol")
query_mode <- get_arg("query", "shared")
agg_mode   <- get_arg("agg", "absmax")
min_ct     <- as.integer(get_arg("min-ct", get_arg("min-celltypes", "0")))
n_perm     <- as.integer(get_arg("nperm", get_arg("n-perm", "2000")))
method     <- get_arg("method", "ks")
transpose  <- has("transpose")
fast_mode  <- has("fast")

query_mode <- match.arg(tolower(query_mode),
                        c("shared", "identical", "per-pipeline", "raw"))
agg_mode   <- match.arg(tolower(agg_mode), c("absmax", "mean"))
gene_id    <- match.arg(tolower(gene_id), c("symbol", "entrez"))
if (is.na(min_ct) || min_ct < 0) stop("--min-ct must be a non-negative integer (0 = OFF)")
if (is.na(n_perm) || n_perm < 1)  stop("--nperm must be a positive integer")

if (fast_mode) {
  source(file.path(script_dir, "connectivity_fast.R"))
  calibrate_connectivity <- calibrate_connectivity_fast
  cat("[11] fast KS calibration enabled (numerically equivalent; scripts/_verify_fast.R)\n")
}

# --parallel N: column-block parallel KS calibration. Numerically identical to
# --fast (same permutation stream, same running sums); needed for the full
# 45956-signature EH3226 reference. Validated by scripts/_verify_parallel.R.
parallel_n <- as.integer(get_arg("parallel", "1"))
block_cols <- as.integer(get_arg("block-cols", "2000"))
# Optional restartable block-level checkpoints for the 45956-signature run.
ckpt_dir <- get_arg("ckpt-dir", NULL)
if (!is.null(ckpt_dir)) options(scda.ckpt = ckpt_dir)
if (is.na(parallel_n) || parallel_n < 1) stop("--parallel must be >= 1")
if (is.na(block_cols) || block_cols < 1) stop("--block-cols must be >= 1")
if (parallel_n > 1L) {
  if (!file.exists(file.path(script_dir, "connectivity_fast.R")))
    stop("--parallel needs connectivity_fast.R in the script directory")
  source(file.path(script_dir, "connectivity_fast.R"))
  source(file.path(script_dir, "connectivity_parallel.R"))
  calibrate_connectivity <- function(query, ref_db, n_perm, method, seed)
    calibrate_connectivity_parallel(query, ref_db, n_perm, method, seed,
                                    n_workers = parallel_n, block_cols = block_cols)
  cat(sprintf("[11] parallel KS calibration enabled (%d workers, %d columns/block)\n",
              parallel_n, block_cols))
}

stopifnot("missing --ref"   = !is.null(ref_path)   && file.exists(ref_path))
stopifnot("missing --audit" = !is.null(audit_path) && file.exists(audit_path))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[11] ref=%s gene_id=%s audit=%s outdir=%s query=%s agg=%s min_ct=%d nperm=%d method=%s transpose=%s fast=%s\n",
            ref_path, gene_id, audit_path, outdir, query_mode, agg_mode,
            min_ct, n_perm, method, transpose, fast_mode))

# ---- reference matrix ----
mat <- readRDS(ref_path)
if (!is.matrix(mat) && !inherits(mat, "dgCMatrix"))
  stop("reference must be a genes x compounds matrix")
if (transpose) mat <- t(mat)
cat("[11] ref dim:", paste(dim(mat), collapse = " x "),
    " | rownames sample:", paste(head(rownames(mat)), collapse = ","), "\n")
if (gene_id == "entrez") {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("gene_id=entrez requires org.Hs.eg.db (BiocManager::install('org.Hs.eg.db'))")
  map <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                               keys = rownames(mat), column = "SYMBOL",
                               keytype = "ENTREZID")
  rownames(mat) <- unname(as.character(map[rownames(mat)]))
  # sync 2026-11-06 with formal library connectivity_lincs.R:
  # drop NA / empty symbols, then de-duplicate (multi-Entrez -> one SYMBOL)
  keep <- !is.na(rownames(mat)) & nzchar(rownames(mat))
  mat <- mat[keep, , drop = FALSE]
  dup <- duplicated(rownames(mat))
  if (any(dup)) mat <- mat[!dup, , drop = FALSE]
  cat("[11] Entrez -> Symbol mapped; dim:", paste(dim(mat), collapse = " x "),
      "| unmapped dropped:", sum(!keep), "| duplicate symbols dropped:", sum(dup), "\n")
}
if (ncol(mat) > 1000 && !fast_mode)
  cat("[11] WARNING: ref has", ncol(mat), "compounds and --fast is OFF; add --fast.\n")

# ---- audit cache ----
aud <- readRDS(audit_path)
if (!is.list(aud) || is.null(aud$candidates))
  stop("audit cache must be a list with $candidates")
has_pipelines <- all(c("naive", "rigorous") %in% names(aud))
cat(sprintf("[11] audit cache: %s | $naive/$rigorous present: %s\n",
            if (has_pipelines) "FORMAL (two pipelines)" else "MERGED (local kit)",
            has_pipelines))

# ---- collapse per-celltype DE tables into ONE named signed-logFC vector ----
# Mirrors the formal runner's collapse(): accepts a LIST of data.frames OR a
# single data.frame; explicit de-duplication (never name-subset a named vector,
# which silently keeps only the first row per gene in R).
collapse <- function(lst, agg, min_ct, filter_sig, label) {
  df <- if (is.data.frame(lst)) lst else do.call(rbind, lst)
  df$gene <- as.character(df$gene)
  if (filter_sig && "padj" %in% names(df))
    df <- df[!is.na(df$padj) & df$padj < 0.05 &
             !is.na(df$logFC) & abs(df$logFC) > 0.25, , drop = FALSE]
  n_rows <- nrow(df); n_dup <- sum(duplicated(df$gene))
  if (min_ct > 1 && "celltype" %in% names(df)) {
    nct <- tapply(as.character(df$celltype), df$gene,
                  function(z) length(unique(z)))
    df <- df[df$gene %in% names(nct)[nct >= min_ct], , drop = FALSE]
  }
  g <- df$gene; fc <- as.numeric(df$logFC)
  ok <- !is.na(g) & nzchar(g) & is.finite(fc)
  g <- g[ok]; fc <- fc[ok]
  sp <- split(fc, factor(g, levels = sort(unique(g))))
  v <- switch(agg,
    absmax = vapply(sp, function(z) z[which.max(abs(z))], numeric(1)),
    mean   = vapply(sp, mean, numeric(1)))
  names(v) <- names(sp)
  stopifnot(!any(duplicated(names(v))))
  cat(sprintf("[11]   [%s] rows=%d | dup gene rows=%d | collapsed by %s | unique genes=%d\n",
              label, n_rows, n_dup, agg, length(v)))
  v
}

cand <- aud$candidates
backbone <- unique(as.character(cand$gene))
if (min_ct > 1) {
  if ("celltype" %in% names(cand)) {
    nct <- tapply(as.character(cand$celltype), as.character(cand$gene),
                  function(z) length(unique(z)))
    keep <- names(nct)[nct >= min_ct]
  } else if (!is.null(aud$evidence) && "n_celltypes" %in% names(aud$evidence)) {
    e <- setNames(aud$evidence$n_celltypes, as.character(aud$evidence$gene))
    keep <- names(e)[!is.na(e) & e >= min_ct]
  } else stop("--min-ct requires candidates$celltype or evidence$n_celltypes")
  # sync 2026-11-06 with the FORMAL library runner: its backbone comes from
  # tapply() and is therefore in sorted-gene order. Keeping $candidates
  # appearance order here makes the KS running sum depend on how ties in the
  # reference are broken (see results/tables/_minct2_order_diagnosis.txt).
  backbone <- unique(keep)
}
cat(sprintf("[11] candidates backbone: %d unique genes%s\n",
            length(backbone), if (min_ct > 1) sprintf(" (min_ct>=%d)", min_ct) else ""))

# ---- build the TWO queries ----
if (has_pipelines) {
  pool_n <- collapse(aud$naive,    agg_mode, 0, FALSE, "naive pool")
  pool_r <- collapse(aud$rigorous, agg_mode, 0, FALSE, "rigorous pool")
  if (query_mode == "shared") {
    genes <- intersect(backbone, intersect(names(pool_n), names(pool_r)))
    q_naive <- pool_n[genes]; q_rigorous <- pool_r[genes]
    cat(sprintf("[11] shared backbone: %d genes present in BOTH pipelines\n", length(genes)))
  } else if (query_mode == "identical") {
    qc <- collapse(cand, agg_mode, min_ct, FALSE, "candidates")
    q_naive <- q_rigorous <- qc
    cat("[11] identical mode: qn and qr are the SAME vector (diagnostic only)\n")
  } else if (query_mode == "per-pipeline") {
    q_naive    <- collapse(aud$naive,    agg_mode, min_ct, TRUE, "naive signature")
    q_rigorous <- collapse(aud$rigorous, agg_mode, min_ct, TRUE, "rigorous signature")
  } else {
    stop("--query raw is only for the merged local-kit cache")
  }
  formal_semantics <- TRUE
} else {
  if (query_mode == "per-pipeline")
    stop("--query per-pipeline requires $naive/$rigorous in the cache (formal cache).")
  cat("[11] !!! LOCAL-KIT APPROXIMATION: cache has no $naive/$rigorous, so both\n")
  cat("[11] !!! queries are built from the single merged $candidates logFC.\n")
  cat("[11] !!! This is NOT the formal 'shared backbone / own logFC' semantics.\n")
  qc <- collapse(cand, agg_mode, min_ct, FALSE, "candidates (merged)")
  q_naive <- q_rigorous <- qc
  formal_semantics <- FALSE
}

# align to reference
q_naive    <- q_naive[intersect(names(q_naive), rownames(mat))]
q_rigorous <- q_rigorous[intersect(names(q_rigorous), rownames(mat))]
cat(sprintf("[11] after ref overlap -> naive=%d rigorous=%d genes\n",
            length(q_naive), length(q_rigorous)))
if (length(q_naive) < 5 || length(q_rigorous) < 5)
  stop("fewer than 5 query genes align with ref; check --gene-id / --transpose")

q_identical <- identical(unname(q_naive), unname(q_rigorous))
if (q_identical) {
  cat("[11] ############################################################\n")
  cat("[11] # WARNING: qn and qr are IDENTICAL. Every compound gets the SAME\n")
  cat("[11] # raw score; the only difference is the permutation seed.\n")
  cat("[11] # The naive-vs-rigorous contrast is BY CONSTRUCTION, not evidence.\n")
  cat("[11] ############################################################\n")
  if (!formal_semantics) {
    # Same signature by construction, not by choice: mirror the calibration
    # instead of letting seed/seed+1 fabricate a spurious difference.
    mirror_results <- TRUE
    cat("[11] merged-cache fallback: rigorous calibration mirrored from naive.\n")
  } else mirror_results <- FALSE
} else {
  common <- intersect(names(q_naive), names(q_rigorous))
  d <- abs(q_naive[common] - q_rigorous[common])
  cat(sprintf("[11] query difference: %d/%d shared genes differ | sign flips=%d | max|d|=%.3f\n",
              sum(d > 1e-9), length(common),
              sum(sign(q_naive[common]) != sign(q_rigorous[common])), max(d)))
  mirror_results <- FALSE
}

# ---- permutation-resolution diagnostic (formal spec §3.4) ----
m_cmp <- ncol(mat)
p_min <- 1 / (n_perm + 1)
bh1   <- 0.05 / m_cmp
res_limited <- p_min > bh1
k_needed <- ceiling(p_min * m_cmp / 0.05)
res_block <- function() {
  cat("\n----- Permutation-resolution diagnostic -----\n")
  cat(sprintf("compounds tested (m) = %d | n_perm = %d\n", m_cmp, n_perm))
  cat(sprintf("smallest attainable permutation p = 1/(n_perm+1) = %.5g\n", p_min))
  cat(sprintf("BH threshold for top-ranked compound = 0.05/m = %.5g\n", bh1))
  if (res_limited) {
    cat(sprintf("=> RESOLUTION-LIMITED: need >= %d compounds tied at the minimal p\n", k_needed))
    cat("   for even the first FDR<0.05 hit. A '0 significant' result is forced\n")
    cat("   BY CONSTRUCTION, NOT evidence of absent drug connectivity.\n")
  } else cat("=> resolution adequate for single-compound significance.\n")
  cat("---------------------------------------------\n")
}
res_block()

# ---- run battery #5 ----
res <- audit_connectivity(q_naive, q_rigorous, ref_db = mat,
                          n_perm = n_perm, method = method)
if (mirror_results) {
  res$rigorous <- res$naive
  res$calibration_rigorous <- res$calibration_naive
}
print_connectivity(res)

# ---- outputs: formal filenames + local-kit extras ----
write.csv(res$calibration_naive,
          file.path(outdir, "lincs_calibration_naive.csv"), row.names = FALSE)
write.csv(res$calibration_rigorous,
          file.path(outdir, "lincs_calibration_rigorous.csv"), row.names = FALSE)
# backwards-compatible names used by earlier kit runs
write.csv(res$calibration_naive,
          file.path(outdir, "connectivity_calibration_naive.csv"), row.names = FALSE)
write.csv(res$calibration_rigorous,
          file.path(outdir, "connectivity_calibration_rigorous.csv"), row.names = FALSE)
aud$connectivity <- res
saveRDS(aud, file.path(outdir, "full_audit_real_connectivity.rds"))

sink(file.path(outdir, "lincs_summary.txt"))
cat("Battery #5 (LINCS L1000 connectivity) - run summary\n")
cat(sprintf("reference: %s | dim %s\n", basename(ref_path), paste(dim(mat), collapse = " x ")))
cat(sprintf("audit cache: %s\n", if (formal_semantics) "FORMAL ($naive/$rigorous)" else
                                  "MERGED local kit ($candidates only) - APPROXIMATION"))
cat(sprintf("query mode: %s | aggregation: %s | min_celltypes: %s\n",
            query_mode, agg_mode, if (min_ct > 1) as.character(min_ct) else "off"))
cat(sprintf("queries identical: %s | formal semantics: %s\n", q_identical, formal_semantics))
cat(sprintf("naive query genes overlapping ref: %d | rigorous: %d\n",
            length(q_naive), length(q_rigorous)))
cat(sprintf("method: %s | n_perm: %d | gene-id: %s\n", method, n_perm, gene_id))
cat(sprintf("naive: %d significant | top %s (score %.3f, q %.4g)\n",
            sum(res$calibration_naive$significant),
            res$naive$top_compound, res$naive$top_score,
            min(res$calibration_naive$q_value)))
cat(sprintf("rigorous: %d significant | top %s (score %.3f, q %.4g)\n",
            sum(res$calibration_rigorous$significant),
            res$rigorous$top_compound, res$rigorous$top_score,
            min(res$calibration_rigorous$q_value)))
res_block()
sink()
cat(sprintf("[11] DONE -> %s (lincs_calibration_*.csv, lincs_summary.txt, full_audit_real_connectivity.rds)\n",
            outdir))
