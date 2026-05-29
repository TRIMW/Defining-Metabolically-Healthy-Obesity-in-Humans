# Step 2: Trajectory inference (slingshot) and differential expression
# along pseudotime (tradeSeq + condiments).
#
# Required packages:
#   BiocManager::install(c("zellkonverter","SingleCellExperiment","slingshot",
#                          "tradeSeq","condiments","BiocParallel"))
#   install.packages(c("ggplot2","dplyr","viridis","patchwork"))
#
# Usage:  Rscript pseudotime_02_tradeseq.R
#   or open in RStudio and run interactively.
#
# Checkpoints (saved automatically, loaded on re-run):
#   ckpt_01_pregam.rds   – after slingshot + subsampling (skip to Stage 1)
#   ckpt_02_stage1.rds   – after Stage 1 fitGAM         (skip to Stage 2)
#   ckpt_03_stage2.rds   – after Stage 2 fitGAM         (skip to tests)
#
# To rerun from scratch:  set FORCE_RERUN <- TRUE  (or delete the .rds files)

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(slingshot)
  library(tradeSeq)
  library(condiments)
  library(BiocParallel)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(viridis)
  library(patchwork)
  library(ggrepel)
})

set.seed(42)
N_CORES     <- 4
register(SnowParam(N_CORES, type = "SOCK"))   # reduce to 1 on Windows
BASE_OUT    <- "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq/outputs/pseudotime/no_myo/nonrefitGAM"
FORCE_RERUN <- FALSE   # TRUE = ignore all checkpoints and recompute everything

CP1 <- file.path(BASE_OUT, "ckpt_01_pregam.rds")   # pre-GAM state
CP2 <- file.path(BASE_OUT, "ckpt_02_stage1.rds")   # Stage 1 done
CP3 <- file.path(BASE_OUT, "ckpt_03_stage2.rds")   # Stage 2 done

resume_from <- if (!FORCE_RERUN && file.exists(CP3)) 3L else
               if (!FORCE_RERUN && file.exists(CP2)) 2L else
               if (!FORCE_RERUN && file.exists(CP1)) 1L else 0L
cat(sprintf("Checkpoint level: %d  (0 = from scratch, 3 = skip to tests)\n\n",
            resume_from))


# ── 1. Load the object written by pseudotime_01_preprocessing.py ─────────────
cat("Loading processed h5ad...\n")
sce <- readH5AD(file.path(BASE_OUT, "fap_preadip_processed.h5ad"),
                reader = "R",
                verbose = FALSE)

assayNames(sce)[assayNames(sce) == "X"] <- "logcounts"
counts_mat  <- assay(sce, "counts")            # keep sparse (dgCMatrix)
pca_harmony <- reducedDim(sce, "X_pca_harmony")

cat(sprintf("Cells: %d    Genes: %d\n", ncol(sce), nrow(sce)))
print(table(sce$cell_type_, sce$depot))


# ══════════════════════════════════════════════════════════════════════════════
# CHECKPOINT 1 BLOCK: slingshot + UMAP + condiments + evaluateK + subsampling
# ══════════════════════════════════════════════════════════════════════════════
if (resume_from >= 1) {
  cat("\nLoading checkpoint 1 (pre-GAM state)...\n")
  cp1               <- readRDS(CP1)
  pseudotime_mat    <- cp1$pseudotime_mat
  cell_weights      <- cp1$cell_weights
  umap_df           <- cp1$umap_df
  curve_coords      <- cp1$curve_coords
  counts_mat        <- cp1$counts_mat        # HVG-filtered counts
  nknots            <- cp1$nknots
  keep_idx          <- cp1$keep_idx
  counts_sub        <- cp1$counts_sub
  pseudotime_sub    <- cp1$pseudotime_sub
  weights_sub       <- cp1$weights_sub
  condition_sub     <- cp1$condition_sub
  U_sex_sub         <- cp1$U_sex_sub
  condition_sub_raw <- cp1$condition_sub_raw
  om_idx            <- cp1$om_idx
  sq_idx            <- cp1$sq_idx
  meta_full         <- cp1$meta_full
  cat(sprintf("  %d subsampled cells, k=%d knots\n",
              length(keep_idx), nknots))
  cat("Condiments results (from prior run):\n")
  cat("Progression test — Cohort:\n"); print(cp1$prog_cohort)
  cat("Progression test — Depot:\n");  print(cp1$prog_depot)
  cat("Differentiation test — Cohort:\n"); print(cp1$diff_cohort)

} else {

  # ── 2. Slingshot trajectory  (FAPs → Pre_adipocytes) ───────────────────────
  cat("\nFitting slingshot trajectory...\n")
  sce <- slingshot(
    sce,
    clusterLabels  = sce$cell_type_,
    reducedDim     = "X_pca_harmony",
    start.clus     = "FAPs",
    end.clus       = "Pre_adipocytes", #c(, "Myofibroblasts"),
    approx_points  = 300,
    extend         = "n",
  )
  pseudotime_mat <- slingPseudotime(sce, na = FALSE)
  cell_weights   <- slingCurveWeights(sce)
  cat(sprintf("Lineages: %d\n", ncol(pseudotime_mat)))

  # ── 3. UMAP plots: slingshot curve overlaid ─────────────────────────────────
  umap_df <- as.data.frame(reducedDim(sce, "X_umap_harmony"))
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$cell_type  <- sce$cell_type_
  umap_df$cohort     <- sce$cohort
  umap_df$depot      <- sce$depot
  umap_df$pseudotime <- pseudotime_mat[, 1]

  umap_mat <- reducedDim(sce, "X_umap_harmony")
  curve_coords <- lapply(seq_along(slingCurves(sce)), function(i) {
    pt       <- pseudotime_mat[, i]
    finite   <- is.finite(pt)
    ord      <- order(pt[finite])
    umap_ord <- umap_mat[finite, , drop = FALSE][ord, ]
    t_seq    <- seq_len(nrow(umap_ord))
    sp1 <- smooth.spline(t_seq, umap_ord[, 1], df = 15)
    sp2 <- smooth.spline(t_seq, umap_ord[, 2], df = 15)
    t_new <- seq(1, nrow(umap_ord), length.out = 300)
    data.frame(UMAP1 = predict(sp1, t_new)$y,
               UMAP2 = predict(sp2, t_new)$y)
  })

  make_umap_plot <- function(color_var, palette = NULL, title = color_var) {
    p <- ggplot(umap_df, aes(UMAP1, UMAP2)) +
      geom_point(aes(color = .data[[color_var]]), size = 0.6, alpha = 0.7) +
      theme_classic(base_size = 11) +
      ggtitle(title) +
      theme(legend.key.size = unit(0.4, "cm"))
    if (!is.null(palette)) p <- p + scale_color_manual(values = palette)
    # for (crv in curve_coords) {
    #   p <- p + geom_path(data = crv, aes(UMAP1, UMAP2),
    #                      color = "black", linewidth = 0.8, inherit.aes = FALSE)
    # }
    p
  }

  p_ct    <- make_umap_plot("cell_type",
                             c(FAPs="#FF9800", Pre_adipocytes="#4CAF50", "Myofibroblasts"="#6e0465"),
                             "Cell type")
  p_coh   <- make_umap_plot("cohort",
                             c(Healthy="#2196F3", Unhealthy="#F44336"))
  p_depot <- make_umap_plot("depot",
                             c(Omentum="#FF9800", Subcutaneous="#9C27B0"))
  p_pt    <- ggplot(umap_df, aes(UMAP1, UMAP2, color = pseudotime)) +
    geom_point(size = 0.6, alpha = 0.7) +
    scale_color_viridis(name = "Pseudotime") +
    theme_classic(base_size = 11) + ggtitle("Pseudotime")
  # for (crv in curve_coords) {
  #   p_pt <- p_pt + geom_path(data = crv, aes(UMAP1, UMAP2),
  #                             color = "white", linewidth = 0.8, inherit.aes = FALSE)
  # }
  umap_panel <- (p_ct | p_coh | p_depot | p_pt)
  ggsave(file.path(BASE_OUT, "slingshot_umap.pdf"),
         umap_panel, width = 25, height = 5)
  cat("  Saved: slingshot_umap.pdf\n")

  # ── 4. condiments: differential trajectory topology ──────────────────────────
  cat("\nRunning condiments tests...\n")
  prog_cohort <- progressionTest(
    pseudotime_mat, cellWeights = cell_weights,
    conditions = sce$cohort,
  )
  cat("Progression test — Cohort:\n"); print(prog_cohort)

  prog_depot <- progressionTest(
    pseudotime_mat, cellWeights = cell_weights,
    conditions = sce$depot,
  )
  cat("Progression test — Depot:\n"); print(prog_depot)

  # diff_cohort <- differentiationTest(
  #   pseudotime_mat, cellWeights = cell_weights,
  #   conditions = sce$cohort
  # )
  # cat("Differentiation test — Cohort:\n"); print(diff_cohort)

  # ── 5. tradeSeq: HVG filter + evaluateK + subsampling + metadata ─────────────
  if ("highly_variable" %in% colnames(rowData(sce))) {
    hvg_mask <- rowData(sce)$highly_variable
  } else {
    message("'highly_variable' not in rowData — using top 3000 variance genes")
    rv       <- sparseMatrixStats::rowVars(assay(sce, "logcounts"))
    hvg_mask <- rank(-rv, ties.method = "first") <= 3000
  }
  counts_mat <- counts_mat[hvg_mask, ]
  cat(sprintf("GAM fitting on %d HVGs (of %d total) — ~%.0fx speedup\n",
              sum(hvg_mask), nrow(sce), nrow(sce) / sum(hvg_mask)))

  # cat("\nEvaluating number of knots (this may take a few minutes)...\n")
  # k_eval <- evaluateK(
  #   counts      = counts_mat,
  #   pseudotime  = pseudotime_mat,
  #   cellWeights = cell_weights,
  #   k           = 3:7,
  #   nGenes      = 300,
  #   verbose     = FALSE,
  #   BPPARAM     = SnowParam(N_CORES, type = "SOCK"),
  # )
  # pdf(file.path(BASE_OUT, "tradeseq_knot_selection.pdf"), width = 6, height = 4)
  # plot(k_eval)
  # dev.off()
  # cat("  Saved: tradeseq_knot_selection.pdf  — inspect to pick k\n")
  nknots <- 5   # adjust after inspecting the plot, then delete CP1 to rerun

  U_sex         <- model.matrix(~ sex, data = data.frame(sex = sce$sex))[, -1, drop = FALSE]
  condition_var <- factor(paste(sce$cohort, sce$depot, sep = "_"))

  MAX_CELLS <- 5000
  if (ncol(sce) > MAX_CELLS) {
    set.seed(42)
    strata   <- paste(sce$cell_type_, sce$cohort, sce$depot, sep = "_")
    keep_idx <- unlist(tapply(seq_len(ncol(sce)), strata, function(idx) {
      sample(idx, max(1L, round(length(idx) * MAX_CELLS / ncol(sce))))
    }))
    keep_idx <- sort(keep_idx)
    cat(sprintf("Subsampled %d → %d cells (stratified)\n",
                ncol(sce), length(keep_idx)))
  } else {
    keep_idx <- seq_len(ncol(sce))
  }
  counts_sub     <- counts_mat[, keep_idx]
  pseudotime_sub <- pseudotime_mat[keep_idx, , drop = FALSE]
  weights_sub    <- cell_weights[keep_idx,  , drop = FALSE]
  condition_sub  <- condition_var[keep_idx]
  U_sex_sub      <- U_sex[keep_idx, , drop = FALSE]

  # Capture metadata NOW — Stage 2 fitGAM will overwrite `sce` with a tradeSeq
  # object that has no cohort/depot/sex columns.
  condition_sub_raw <- data.frame(
    cohort = sce$cohort[keep_idx],
    depot  = sce$depot[keep_idx],
    sex    = sce$sex[keep_idx]
  )
  om_idx <- which(condition_sub_raw$depot == "Omentum")
  sq_idx <- which(condition_sub_raw$depot == "Subcutaneous")

  meta_full <- data.frame(
    cell_type   = sce$cell_type_,
    cohort      = sce$cohort,
    depot       = sce$depot,
    sex         = sce$sex,
    participant = sce$common_participant
  )

  saveRDS(list(
    pseudotime_mat    = pseudotime_mat,
    cell_weights      = cell_weights,
    umap_df           = umap_df,
    curve_coords      = curve_coords,
    counts_mat        = counts_mat,
    nknots            = nknots,
    keep_idx          = keep_idx,
    counts_sub        = counts_sub,
    pseudotime_sub    = pseudotime_sub,
    weights_sub       = weights_sub,
    condition_sub     = condition_sub,
    U_sex_sub         = U_sex_sub,
    condition_sub_raw = condition_sub_raw,
    om_idx            = om_idx,
    sq_idx            = sq_idx,
    meta_full         = meta_full,
    prog_cohort       = prog_cohort,
    prog_depot        = prog_depot,
    diff_cohort       = NULL #diff_cohort
  ), CP1)
  cat(sprintf("  Checkpoint 1 saved → %s\n", basename(CP1)))
}


# ══════════════════════════════════════════════════════════════════════════════
# CHECKPOINT 2 BLOCK: Stage 1 fitGAM
# ══════════════════════════════════════════════════════════════════════════════
if (resume_from >= 2) {
  cat("\nLoading checkpoint 2 (Stage 1 results)...\n")
  cp2           <- readRDS(CP2)
  sce_s1        <- cp2$sce_s1
  assoc_stage1  <- cp2$assoc_stage1
  sig_genes     <- cp2$sig_genes
  counts_stage2 <- cp2$counts_stage2
  assoc_test    <- cp2$assoc_test
  cat(sprintf("  %d trajectory-associated genes\n", length(sig_genes)))

} else {
  cat(sprintf("\n[Stage 1] Trajectory GAMs on %d HVGs, %d cells (no conditions)...\n",
              nrow(counts_sub), length(keep_idx)))
  sce_s1 <- fitGAM(
    counts      = counts_sub,
    pseudotime  = pseudotime_sub,
    cellWeights = weights_sub,
    nknots      = nknots,
    BPPARAM     = SnowParam(N_CORES, type = "SOCK"),
    verbose     = TRUE,
  )

  assoc_stage1 <- associationTest(sce_s1, l2fc = log2(1.5))
  sig_genes <- rownames(assoc_stage1)[
    !is.na(assoc_stage1$pvalue) &
    p.adjust(assoc_stage1$pvalue, "BH") < 0.05
  ]
  if (length(sig_genes) < 100) {
    sig_genes <- rownames(assoc_stage1)[order(assoc_stage1$pvalue)][
      seq_len(min(500L, nrow(assoc_stage1)))
    ]
    cat("  FDR < 0.05 gave < 100 genes; using top 500 by p-value.\n")
  }
  cat(sprintf("  %d / %d HVGs are trajectory-associated — proceeding to Stage 2.\n",
              length(sig_genes), nrow(counts_sub)))

  counts_stage2 <- counts_sub[sig_genes, ]
  assoc_test    <- assoc_stage1[order(assoc_stage1$pvalue), ]

  saveRDS(list(
    sce_s1        = sce_s1,
    assoc_stage1  = assoc_stage1,
    sig_genes     = sig_genes,
    counts_stage2 = counts_stage2,
    assoc_test    = assoc_test
  ), CP2)
  cat(sprintf("  Checkpoint 2 saved → %s\n", basename(CP2)))
}


# ══════════════════════════════════════════════════════════════════════════════
# CHECKPOINT 3 BLOCK: Stage 2 fitGAM
# ══════════════════════════════════════════════════════════════════════════════
if (resume_from >= 3) {
  cat("\nLoading checkpoint 3 (Stage 2 results)...\n")
  cp3           <- readRDS(CP3)
  sce           <- cp3$sce
  counts_stage2 <- cp3$counts_stage2
  cat("  Stage 2 GAMs loaded.\n")

} else {
  cat(sprintf("\n[Stage 2] Condition GAMs on %d genes (k=%d, %d conditions)...\n",
              length(sig_genes), nknots, nlevels(condition_sub)))
  sce <- fitGAM(
    counts      = counts_stage2,
    pseudotime  = pseudotime_sub,
    cellWeights = weights_sub,
    nknots      = nknots,
    conditions  = condition_sub,
    U           = U_sex_sub,
    BPPARAM     = SnowParam(N_CORES, type = "SOCK"),
    verbose     = TRUE,
  )
  cat("Stage 2 GAM fitting complete.\n")

  saveRDS(list(sce = sce, counts_stage2 = counts_stage2), CP3)
  cat(sprintf("  Checkpoint 3 saved → %s\n", basename(CP3)))
}


# ── 6. Statistical tests ──────────────────────────────────────────────────────
cat("\nRunning tradeSeq tests...\n")

# A. Association: reuse Stage 1 result (covers all HVGs, not just sig_genes)
assoc_test <- assoc_stage1[order(assoc_stage1$pvalue), ]

# B. Condition tests — single call on the Stage 2 model; no refitting needed.
# global=TRUE  → omnibus test (any condition differs across all 4 groups)
# pairwise=TRUE → all pairwise columns added: waldStat_ivsj / df_ivsj / pvalue_ivsj
#
# The 4-level condition factor levels (alphabetical order):
#   1 = Healthy_Omentum   2 = Healthy_Subcutaneous
#   3 = Unhealthy_Omentum 4 = Unhealthy_Subcutaneous
# Verify with the printout below; adjust extract_pair() calls if order differs.
cat("  Condition factor levels:\n")
for (i in seq_along(levels(condition_sub)))
  cat(sprintf("    %d: %s\n", i, levels(condition_sub)[i]))

cat("  Running conditionTest (global + pairwise)...\n")
cond_all <- conditionTest(sce, global = TRUE, pairwise = TRUE, l2fc = log2(1.5))

# Helper: pull a tidy data.frame for one pairwise contrast
extract_pair <- function(i, j) {
  sfx  <- sprintf("conds%dvs%d", i, j)
  cols <- paste0(c("waldStat_", "df_", "pvalue_"), sfx)
  res  <- cond_all[, cols, drop = FALSE]
  colnames(res) <- c("waldStat", "df", "pvalue")
  res[order(res$pvalue), ]
}

# Omnibus: any of the 4 conditions differs
cond_test_overall   <- cond_all[, c("waldStat", "df", "pvalue"), drop = FALSE]
cond_test_overall   <- cond_test_overall[order(cond_test_overall$pvalue), ]

# Cohort effect (Healthy vs Unhealthy) within each depot
cond_test_cohort_om <- extract_pair(1, 3)   # Healthy_Om  vs Unhealthy_Om
cond_test_cohort_sq <- extract_pair(2, 4)   # Healthy_SQ  vs Unhealthy_SQ

# Depot effect (Omentum vs Subcutaneous) within each cohort
cond_test_depot_h   <- extract_pair(1, 2)   # Healthy_Om  vs Healthy_SQ
cond_test_depot_u   <- extract_pair(3, 4)   # Unhealthy_Om vs Unhealthy_SQ

# C. Start-vs-end test (single lineage; diffEndTest requires ≥2 lineages)
end_test <- startVsEndTest(sce, pseudotimeValues = c(0, 1), l2fc = log2(1.5))
end_test <- end_test[order(end_test$pvalue), ]

# # D. Pattern test: genes with different expression shapes across conditions
# pat_test <- patternTest(sce, l2fc = log2(1.5))
# pat_test <- pat_test[order(pat_test$pvalue), ]

write.csv(assoc_test,           file.path(BASE_OUT, "tradeseq_association_test.csv"))
write.csv(cond_test_overall,    file.path(BASE_OUT, "tradeseq_condition_overall.csv"))
write.csv(cond_test_cohort_om,  file.path(BASE_OUT, "tradeseq_condition_cohort_omentum.csv"))
write.csv(cond_test_cohort_sq,  file.path(BASE_OUT, "tradeseq_condition_cohort_subcutaneous.csv"))
write.csv(cond_test_depot_h,    file.path(BASE_OUT, "tradeseq_condition_depot_healthy.csv"))
write.csv(cond_test_depot_u,    file.path(BASE_OUT, "tradeseq_condition_depot_unhealthy.csv"))
write.csv(end_test,             file.path(BASE_OUT, "tradeseq_endpoint_test.csv"))
cat("  Saved all test results.\n")

fdr <- 0.05
summarise_test <- function(df, label) {
  n <- sum(p.adjust(df$pvalue, "BH") < fdr, na.rm = TRUE)
  cat(sprintf("  %-50s %d genes\n", label, n))
}
cat(sprintf("\nFDR < %.2f summary:\n", fdr))
summarise_test(assoc_test,          "Association (trajectory):")
summarise_test(cond_test_overall,   "Condition — omnibus (4-level):")
summarise_test(cond_test_cohort_om, "Condition — cohort within Omentum:")
summarise_test(cond_test_cohort_sq, "Condition — cohort within Subcutaneous:")
summarise_test(cond_test_depot_h,   "Condition — depot within Healthy:")
summarise_test(cond_test_depot_u,   "Condition — depot within Unhealthy:")
summarise_test(end_test,            "Start vs end:")


# ── 7. Visualise top hits ─────────────────────────────────────────────────────
top_assoc          <- rownames(assoc_test)[1:min(9, nrow(assoc_test))]
top_cond_cohort_om <- rownames(cond_test_cohort_om)[1:min(9, nrow(cond_test_cohort_om))]
top_cond_depot_h   <- rownames(cond_test_depot_h)[1:min(9, nrow(cond_test_depot_h))]

# plotSmoothers with a condition-fitted model only accepts one gene at a time;
# loop and tile with patchwork.
plot_smoother_grid <- function(genes, file, ncol = 3) {
  genes <- genes[genes %in% rownames(counts_stage2)]
  plots <- lapply(genes, function(g)
    plotSmoothers(sce, counts_stage2, gene = g, alpha = 0.4, border = TRUE) +
      ggtitle(g) + theme(plot.title = element_text(size = 9))
  )
  panel <- wrap_plots(plots, ncol = ncol)
  pdf(file, width = 18, height = ceiling(length(genes) / ncol) * 3)
  print(panel)
  dev.off()
}

plot_smoother_grid(top_assoc,          file.path(BASE_OUT, "tradeseq_top_associated_genes.pdf"))
plot_smoother_grid(top_cond_cohort_om, file.path(BASE_OUT, "tradeseq_top_cohort_genes.pdf"))
plot_smoother_grid(top_cond_depot_h,   file.path(BASE_OUT, "tradeseq_top_depot_genes.pdf"))
cat("  Saved: tradeseq_top_*_genes.pdf\n")

pt_df <- data.frame(
  pseudotime  = pseudotime_mat[, 1],
  cell_type   = meta_full$cell_type,
  cohort      = meta_full$cohort,
  depot       = meta_full$depot,
  sex         = meta_full$sex,
  participant = meta_full$participant
)

p_pt_box <- ggplot(pt_df, aes(x = cohort, y = pseudotime, fill = cohort)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.4) +
  facet_grid(cell_type ~ depot) +
  scale_fill_manual(values = c(Healthy="#2196F3", Unhealthy="#F44336")) +
  theme_classic(base_size = 11) +
  ggtitle("Pseudotime by cohort, depot and cell type") +
  ylab("Slingshot pseudotime")

ggsave(file.path(BASE_OUT, "pseudotime_violin_by_group.pdf"),
       p_pt_box, width = 8, height = 6)
cat("  Saved: pseudotime_violin_by_group.pdf\n")


# ── 7b–c. DE scatter plots: log2FC and signed −log10(p) ──────────────────────
cat("\nComputing per-condition predicted expression for log2FC...\n")

# predictSmooth returns fitted count-scale predictions (exp of the NB log-link)
# for each condition × pseudotime point. Average across pseudotime → one value
# per gene × condition; log2 ratios give signed fold changes.
genes_fc <- Reduce(intersect, list(
  rownames(cond_test_cohort_om), rownames(cond_test_cohort_sq),
  rownames(cond_test_depot_h),   rownames(cond_test_depot_u)
))

pred_tidy <- predictSmooth(sce, gene = genes_fc, nPoints = 50, tidy = TRUE)

pred_mean <- pred_tidy %>%
  group_by(gene, condition) %>%
  summarise(mean_expr = mean(yhat, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = mean_expr)

# Identify condition columns by pattern (robust to factor level order)
lv    <- levels(condition_sub)
lv_ho <- grep("Healthy.*Om|Om.*Healthy",      lv, value = TRUE)  # Healthy_Omentum
lv_hs <- grep("Healthy.*Sub|Sub.*Healthy",    lv, value = TRUE)  # Healthy_Subcutaneous
lv_uo <- grep("Unhealthy.*Om|Om.*Unhealthy",  lv, value = TRUE)  # Unhealthy_Omentum
lv_us <- grep("Unhealthy.*Sub|Sub.*Unhealthy",lv, value = TRUE)  # Unhealthy_Subcutaneous

# Shared scatter helper — takes a prepared data.frame with columns
#   xvar, yvar, sig_cat, and genes to label
make_de_scatter <- function(df, xvar, yvar, xlab, ylab, title, pal, top_genes) {
  ggplot(df, aes(.data[[xvar]], .data[[yvar]], color = sig_cat)) +
    geom_point(data = \(d) filter(d, sig_cat == "Neither"),
               alpha = 0.25, size = 0.7) +
    geom_point(data = \(d) filter(d, sig_cat != "Neither"),
               alpha = 0.75, size = 1.1) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                linewidth = 0.3, colour = "grey40") +
    geom_text_repel(
      data          = \(d) filter(d, gene %in% top_genes),
      aes(label     = gene),
      size          = 2.4,
      max.overlaps  = Inf,
      min.segment.length = 0,
      segment.size  = 0.3,
      segment.color = "grey50",
      box.padding   = 0.35,
      point.padding = 0.2,
      show.legend   = FALSE
    ) +
    scale_color_manual(values = pal, name = NULL) +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 10)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 10)) +
    coord_cartesian(clip = "off") +
    labs(x = xlab, y = ylab, title = title) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "top",
      legend.text     = element_text(size = 8),
      plot.margin     = margin(t = 8, r = 35, b = 8, l = 8, unit = "mm")
    )
}

save_scatter_pair <- function(df, fc_x, fc_y, sp_x, sp_y,
                               xlab_fc, ylab_fc, xlab_sp, ylab_sp,
                               title_fc, title_sp, pal, top_genes, stem) {
  p_fc <- make_de_scatter(df, fc_x, fc_y, xlab_fc, ylab_fc, title_fc, pal, top_genes)
  p_sp <- make_de_scatter(df, sp_x, sp_y, xlab_sp, ylab_sp, title_sp, pal, top_genes)
  ggsave(file.path(BASE_OUT, paste0(stem, "_log2fc.pdf")),   p_fc, width = 6, height = 5.5)
  ggsave(file.path(BASE_OUT, paste0(stem, "_signed_p.pdf")), p_sp, width = 6, height = 5.5)
  ggsave(file.path(BASE_OUT, paste0(stem, "_combined.pdf")), p_fc | p_sp, width = 12, height = 5.5)
  cat(sprintf("  Saved: %s_*.pdf\n", stem))
}

top_genes_from <- function(df, pv1, pv2, n = 35) {
  df %>%
    filter(sig_cat != "Neither") %>%
    mutate(rs = rank(.data[[pv1]], na.last = TRUE) +
               rank(.data[[pv2]], na.last = TRUE)) %>%
    arrange(rs) %>% slice_head(n = n) %>% pull(gene)
}

# ── 7b. Cohort effect (Unhealthy vs Healthy) within each depot ───────────────
cohort_df <- pred_mean %>%
  mutate(
    log2fc_om = log2(.data[[lv_uo]] / .data[[lv_ho]]),
    log2fc_sq = log2(.data[[lv_us]] / .data[[lv_hs]])
  ) %>%
  left_join(data.frame(gene     = rownames(cond_test_cohort_om),
                       pvalue_om = cond_test_cohort_om$pvalue,
                       padj_om   = p.adjust(cond_test_cohort_om$pvalue, "BH"),
                       stringsAsFactors = FALSE), by = "gene") %>%
  left_join(data.frame(gene     = rownames(cond_test_cohort_sq),
                       pvalue_sq = cond_test_cohort_sq$pvalue,
                       padj_sq   = p.adjust(cond_test_cohort_sq$pvalue, "BH"),
                       stringsAsFactors = FALSE), by = "gene") %>%
  mutate(
    sig_cat = factor(case_when(
      padj_om < 0.05 & padj_sq < 0.05 ~ "Both",
      padj_om < 0.05                   ~ "Omentum only",
      padj_sq < 0.05                   ~ "Subcutaneous only",
      TRUE                             ~ "Neither"
    ), levels = c("Both", "Omentum only", "Subcutaneous only", "Neither")),
    signed_p_om = -log10(pmax(pvalue_om, 1e-300)) * sign(log2fc_om),
    signed_p_sq = -log10(pmax(pvalue_sq, 1e-300)) * sign(log2fc_sq)
  )

cohort_pal <- c("Both"="#E91E63", "Omentum only"="#FF9800",
                "Subcutaneous only"="#9C27B0", "Neither"="grey75")

save_scatter_pair(
  cohort_df,
  "log2fc_om", "log2fc_sq", "signed_p_om", "signed_p_sq",
  "log2FC (Unhealthy/Healthy) — Omentum",
  "log2FC (Unhealthy/Healthy) — Subcutaneous",
  expression(-log[10](p) %.% sign(FC) ~ "— Omentum"),
  expression(-log[10](p) %.% sign(FC) ~ "— Subcutaneous"),
  "Cohort DE: within Omentum vs within Subcutaneous",
  "Cohort DE: signed significance",
  cohort_pal,
  top_genes_from(cohort_df, "pvalue_om", "pvalue_sq"),
  "cohort_de"
)

write.csv(cohort_df[, c("gene", "log2fc_om", "log2fc_sq",
                        "pvalue_om", "padj_om", "pvalue_sq", "padj_sq", "sig_cat")],
          file.path(BASE_OUT, "cohort_de_fc_summary.csv"), row.names = FALSE)
cat("  Saved: cohort_de_fc_summary.csv\n")

# ── 7c. Depot effect (Subcutaneous vs Omentum) within each cohort ────────────
depot_df <- pred_mean %>%
  mutate(
    log2fc_h = log2(.data[[lv_hs]] / .data[[lv_ho]]),   # SQ/Om within Healthy
    log2fc_u = log2(.data[[lv_us]] / .data[[lv_uo]])    # SQ/Om within Unhealthy
  ) %>%
  left_join(data.frame(gene    = rownames(cond_test_depot_h),
                       pvalue_h = cond_test_depot_h$pvalue,
                       padj_h   = p.adjust(cond_test_depot_h$pvalue, "BH"),
                       stringsAsFactors = FALSE), by = "gene") %>%
  left_join(data.frame(gene    = rownames(cond_test_depot_u),
                       pvalue_u = cond_test_depot_u$pvalue,
                       padj_u   = p.adjust(cond_test_depot_u$pvalue, "BH"),
                       stringsAsFactors = FALSE), by = "gene") %>%
  mutate(
    sig_cat = factor(case_when(
      padj_h < 0.05 & padj_u < 0.05 ~ "Both",
      padj_h < 0.05                  ~ "Healthy only",
      padj_u < 0.05                  ~ "Unhealthy only",
      TRUE                           ~ "Neither"
    ), levels = c("Both", "Healthy only", "Unhealthy only", "Neither")),
    signed_p_h = -log10(pmax(pvalue_h, 1e-300)) * sign(log2fc_h),
    signed_p_u = -log10(pmax(pvalue_u, 1e-300)) * sign(log2fc_u)
  )

depot_pal <- c("Both"="#9C27B0", "Healthy only"="#2196F3",
               "Unhealthy only"="#F44336", "Neither"="grey75")
save_scatter_pair(
  depot_df,
  "log2fc_h", "log2fc_u", "signed_p_h", "signed_p_u",
  "log2FC (Subcutaneous/Omentum) — Healthy",
  "log2FC (Subcutaneous/Omentum) — Unhealthy",
  expression(-log[10](p) %.% sign(FC) ~ "— Healthy"),
  expression(-log[10](p) %.% sign(FC) ~ "— Unhealthy"),
  "Depot DE: within Healthy vs within Unhealthy",
  "Depot DE: signed significance",
  depot_pal,
  top_genes_from(depot_df, "pvalue_h", "pvalue_u"),
  "depot_de"
)

write.csv(depot_df[, c("gene", "log2fc_h", "log2fc_u",
                       "pvalue_h", "padj_h", "pvalue_u", "padj_u", "sig_cat")],
          file.path(BASE_OUT, "depot_de_fc_summary.csv"), row.names = FALSE)
cat("  Saved: depot_de_fc_summary.csv\n")


# ── 8. Save the full SCE object ───────────────────────────────────────────────
saveRDS(sce, file.path(BASE_OUT, "sce_tradeseq_fitted.rds"))
cat("\nSaved: sce_tradeseq_fitted.rds\n")
cat("\nDone.\n")


### Some genes trajectories
### 


# laod cp3:
CP3 <- file.path(BASE_OUT, "ckpt_03_stage2.rds")   # Stage 2 done
cp3 <- readRDS(CP3)
sce <- cp3$sce
counts_stage2 <- cp3$counts_stage2

genes <- strsplit("AADAC, CCL2, IL6, PI16, DEPP1, SPARCL1, CLU", ", ")[[1]]

plot_smoother_grid(genes, file.path(BASE_OUT, "tradeseq_top_interesting_genes.pdf"))
