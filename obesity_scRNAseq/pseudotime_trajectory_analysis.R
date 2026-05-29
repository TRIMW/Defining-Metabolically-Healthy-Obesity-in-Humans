# pseudotime_trajectory_analysis.R
#
# Three trajectory analyses using Slingshot + condiments + tradeSeq:
#   1. Omentum : adipocyte progenitor differentiation  (FAPs → Pre_adipocytes)
#   2. SubQ    : adipocyte progenitor differentiation  (FAPs → Pre-adipocytes)
#   3. Omentum : mesothelial transition                (Mesothelial_cells → MMT_cells)
#
# Differential trajectory (T2D vs Healthy):
#   condiments::progressionTest()     — do cells progress differently along trajectory?
#   condiments::differentiationTest() — do cells split to different fates (branching only)?
#   tradeSeq::conditionTest()         — which genes differ in their pseudotime kinetics?
#   Sex is included as a covariate in tradeSeq GAM where available (omentum only).

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(scater)
  library(scuttle)
  library(scran)
  library(limma)          # removeBatchEffect for sex regression before PCA
  library(harmony)        # cohort integration before trajectory fitting
  library(slingshot)
  library(condiments)
  library(tradeSeq)
  library(tidyverse)
  library(BiocParallel)
  library(patchwork)
  library(pheatmap)
  library(ggrepel)
})

# ── Parameters ──────────────────────────────────────────────────────────────────
BASE_DIR   <- "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq"
INPUT_DIR  <- file.path(BASE_DIR, "inputs/anndatas")
OUTPUT_DIR <- file.path(BASE_DIR, "results/pseudotime")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_CORES      <- 4
N_HVG        <- 2000
N_PC         <- 30
N_KNOTS      <- 6        # starting point; tune with evaluateK() below
L2FC_THRESH  <- log2(1.5) # |log2FC| > 1 filter in conditionTest
set.seed(42)
register(SnowParam(N_CORES, type = "SOCK"))  # MulticoreParam (fork) is unreliable on macOS

# ── Load datasets ────────────────────────────────────────────────────────────────
cat("Loading omentum...\n")
sce_om <- readH5AD(file.path(INPUT_DIR, "omentum_reannotated_finely_truncated.h5ad"),
                   use_hdf5 = FALSE)
# zellkonverter reads 'counts' layer into assay("counts"); counts() accessor works.

cat("Loading subq...\n")
sce_sq <- readH5AD(file.path(INPUT_DIR, "sq_reannotated_finely_truncated.h5ad"),
                   use_hdf5 = FALSE)

# ── Helper functions ─────────────────────────────────────────────────────────────

# Subset → normalize → (optionally regress sex) → recompute PCA + UMAP.
# We always recompute rather than reuse the global PCA: the global embedding
# reflects all cell types and compresses the local progenitor/transition structure.
prep_subset <- function(sce,
                        cell_types,
                        cell_type_col   = "cell_type_",
                        regress_sex     = FALSE,
                        integrate_cohort = FALSE,
                        n_hvg           = N_HVG,
                        n_pc            = N_PC) {

  sub <- sce[, sce[[cell_type_col]] %in% cell_types]
  cat("  Cells:", ncol(sub), "\n")
  cat("  Cell type × cohort breakdown:\n")
  print(table(sub[[cell_type_col]], sub$cohort))

  sub <- logNormCounts(sub, assay.type = "counts")
  dec  <- modelGeneVar(sub)
  hvgs <- getTopHVGs(dec, n = n_hvg)

  if (regress_sex && !all(is.na(sub$sex))) {
    # When integrating cohort, only protect cell-type structure (not cohort) during
    # sex regression — Harmony will handle cohort correction afterwards.
    protect_mat <- if (integrate_cohort) model.matrix(~ sub[[cell_type_col]])
                   else model.matrix(~ sub$cohort)
    sex_mat  <- model.matrix(~ sub$sex)[, -1, drop = FALSE]
    corrected <- removeBatchEffect(logcounts(sub)[hvgs, ],
                                   covariates = sex_mat,
                                   design     = protect_mat)
    assay(sub, "logcounts_nosex", withDimnames = FALSE) <-
      Matrix::Matrix(0, nrow = nrow(sub), ncol = ncol(sub))
    assay(sub, "logcounts_nosex")[hvgs, ] <- corrected
    sub <- runPCA(sub, exprs_values = "logcounts_nosex", subset_row = hvgs,
                  ncomponents = n_pc, name = "PCA")
  } else {
    sub <- runPCA(sub, subset_row = hvgs, ncomponents = n_pc, name = "PCA")
  }

  # Cohort integration: correct the PCA embedding so the trajectory reflects
  # cell state (progenitor → committed) rather than Healthy/Unhealthy separation.
  # Gene expression (used by tradeSeq) is untouched — only the embedding changes.
  if (integrate_cohort) {
    harmony_vars <- "cohort"
    if (!all(is.na(sub$sex))) harmony_vars <- c(harmony_vars, "sex")

    # Try Harmony; fall back to linear PCA correction if LAPACK is unavailable.
    # The fallback applies removeBatchEffect to the PC coordinates — equivalent
    # to Harmony's first iteration and sufficient for most datasets.
    sub <- tryCatch({
      cat("  Running Harmony", as.character(packageVersion("harmony")),
          "(", paste(harmony_vars, collapse = " + "), ")...\n")
      RunHarmony(sub, group.by.vars = harmony_vars,
                 reduction.save = "HARMONY", verbose = FALSE)
    }, error = function(e) {
      cat("  Harmony", as.character(packageVersion("harmony")),
          "failed (", conditionMessage(e), ")\n")
      cat("  Falling back to linear PCA correction via removeBatchEffect\n")
      pca_mat   <- reducedDim(sub, "PCA")
      batch2    <- if (length(harmony_vars) > 1) sub$sex else NULL
      corrected <- t(removeBatchEffect(t(pca_mat),
                                       batch  = sub$cohort,
                                       batch2 = batch2))
      reducedDim(sub, "HARMONY") <- corrected
      sub
    })

    sub <- runUMAP(sub, dimred = "HARMONY", n_dimred = min(n_pc, 20), name = "UMAP")
    metadata(sub)$traj_dimred <- "HARMONY"
    cat("  Trajectory will use HARMONY embedding\n")
  } else {
    sub <- runUMAP(sub, dimred = "PCA", n_dimred = min(n_pc, 20), name = "UMAP")
    metadata(sub)$traj_dimred <- "PCA"
  }

  metadata(sub)$hvgs <- hvgs
  sub
}

# Fit Slingshot trajectories on a prepped SCE.
# start_clus must match a level in sub[[cell_type_col]].
# end_clus is optional (NULL = infer from data); specify to force a terminal state.
fit_slingshot <- function(sub,
                          start_clus,
                          end_clus      = NULL,
                          cell_type_col = "cell_type_") {
  # Use the dimred stored by prep_subset (PCA or HARMONY).
  dimred <- if (!is.null(metadata(sub)$traj_dimred)) metadata(sub)$traj_dimred else "PCA"
  cat("  Slingshot fitting on:", dimred, "\n")
  slingshot(
    data          = sub,
    clusterLabels = cell_type_col,
    reducedDim    = dimred,
    start.clus    = start_clus,
    end.clus      = end_clus,
    approx_points = 200
  )
}

# condiments differential trajectory tests.
# conditions: factor with levels c("Healthy","Unhealthy")
run_condiments <- function(sds, conditions, nrep = 500, out_prefix) {
  cat("  condiments: imbalance score...\n")
  imb <- imbalance_score(sds, dimred = "UMAP", conditions = conditions, k = 10)

  cat("  condiments: progressionTest...\n")
  prog <- progressionTest(sds, conditions = conditions,
                          global = TRUE, lineages = TRUE, rep = nrep)

  diff <- NULL
  if (length(slingLineages(sds)) > 1) {
    cat("  condiments: differentiationTest (", length(slingLineages(sds)), "lineages)...\n")
    diff <- differentiationTest(sds, conditions = conditions) #, nrep = nrep)
  }

  saveRDS(list(imbalance = imb, progression = prog, differentiation = diff),
          paste0(out_prefix, "_condiments.rds"))
  cat("  progressionTest:\n"); print(prog)
  if (!is.null(diff)) { cat("  differentiationTest:\n"); print(diff) }

  list(imbalance = imb, progression = prog, differentiation = diff)
}

# Tune number of knots then fit tradeSeq GAM.
# U: covariate matrix (sex); pass NULL when sex unavailable.
# Returns list(sce_gam, condition_res, pattern_res).
run_tradeseq <- function(sds, conditions, U = NULL,
                         n_knots = N_KNOTS, out_prefix) {
  sds_obj <- as.SlingshotDataSet(sds)

  # Filter to HVGs that also pass a minimum expression threshold.
  # fitGAM on all ~26k genes would take days; ~2k HVGs takes ~30-60 min.
  hvgs    <- metadata(sds)$hvgs
  cnt_all <- counts(sds)
  expressed <- rowSums(cnt_all > 0) >= 0.05 * ncol(cnt_all)
  keep    <- intersect(hvgs, rownames(cnt_all)[expressed])
  cnt     <- cnt_all[keep, ]
  cat("  Genes for GAM:", nrow(cnt), "(HVGs ∩ expressed in ≥5% of cells)\n")

  # ── Optional: uncomment to tune knots (slow but informative) ──────────────────
  # k_eval <- evaluateK(counts = cnt, sds = sds_obj, conditions = conditions,
  #                     nGenes = 300, verbose = FALSE, parallel = TRUE)
  # pdf(paste0(out_prefix, "_evaluateK.pdf"))
  # plot_evalutateK_results(k_eval)   # typo is in the tradeSeq package itself
  # dev.off()
  # n_knots <- as.integer(readline("Enter chosen n_knots: "))

  cat("  Fitting GAM (nknots =", n_knots, ")...\n")
  sce_gam <- fitGAM(
    counts     = cnt,
    sds        = sds_obj,
    conditions = conditions,
    U          = U,
    nknots     = n_knots,
    parallel   = TRUE,
    BPPARAM    = SnowParam(N_CORES, type = "SOCK"),
    verbose    = FALSE
  )
  saveRDS(sce_gam, paste0(out_prefix, "_sce_gam.rds"))

  cat("  conditionTest...\n")
  cond_res <- conditionTest(sce_gam, l2fc = L2FC_THRESH) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    arrange(pvalue)

  # patternTest compares patterns across lineages — only valid with ≥2 lineages.
  n_lineages <- length(slingLineages(sds))
  patt_res   <- NULL
  if (n_lineages >= 2) {
    cat("  patternTest (", n_lineages, "lineages)...\n")
    patt_res <- patternTest(sce_gam) %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      arrange(pvalue)
    write_csv(patt_res, paste0(out_prefix, "_patternTest.csv"))
  } else {
    cat("  patternTest skipped: single lineage\n")
  }

  write_csv(cond_res, paste0(out_prefix, "_conditionTest.csv"))

  # Genes associated with pseudotime (any non-flat pattern along trajectory).
  cat("  associationTest...\n")
  assoc_res <- associationTest(sce_gam, global = TRUE, lineages = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    mutate(padj = p.adjust(pvalue, method = "BH")) %>%
    arrange(pvalue)
  write_csv(assoc_res, paste0(out_prefix, "_associationTest.csv"))

  # Genes that differ between trajectory start and end (with direction).
  # logFC is derived from predictSmooth: log2(mean fitted value at end / start).
  cat("  startVsEndTest...\n")
  sve_res <- startVsEndTest(sce_gam, global = TRUE, lineages = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    mutate(padj = p.adjust(pvalue, method = "BH"))

  # Compute direction: predictSmooth at 2 points (start, end) → mean log2FC
  pred2 <- predictSmooth(sce_gam, gene = sve_res$gene, nPoints = 2, tidy = FALSE)
  # pred2 is genes × (2 * n_lineages); odd columns = start, even = end
  start_cols <- seq(1, ncol(pred2), by = 2)
  end_cols   <- seq(2, ncol(pred2), by = 2)
  logfc      <- rowMeans(log2(pred2[, end_cols, drop = FALSE] + 1) -
                         log2(pred2[, start_cols, drop = FALSE] + 1))
  sve_res$logFC_start_to_end <- logfc[sve_res$gene]
  sve_res <- arrange(sve_res, pvalue)
  write_csv(sve_res, paste0(out_prefix, "_startVsEndTest.csv"))

  cat("  Top start→end genes (up along trajectory):\n")
  print(head(sve_res %>% filter(!is.na(pvalue), logFC_start_to_end > 0), 10))
  cat("  Top start→end genes (down along trajectory):\n")
  print(head(sve_res %>% filter(!is.na(pvalue), logFC_start_to_end < 0), 10))

  cat("  Top condition-differential genes:\n")
  print(head(cond_res %>% filter(!is.na(pvalue)), 20))

  list(sce_gam    = sce_gam,    condition_res = cond_res,
       pattern_res = patt_res,  assoc_res     = assoc_res,
       sve_res     = sve_res,   cnt           = cnt)
}

# Diagnostic plots: UMAP coloured by pseudotime / cohort / cell type,
# plus pseudotime density split by cohort.
plot_trajectory <- function(sds, cell_type_col, title, out_prefix) {
  umap_mat <- reducedDim(sds, "UMAP")
  pt_mat   <- slingPseudotime(sds)

  df <- data.frame(
    UMAP1     = umap_mat[, 1],
    UMAP2     = umap_mat[, 2],
    pseudotime = rowMeans(pt_mat, na.rm = TRUE),
    cohort    = factor(sds$cohort, levels = c("Healthy", "Unhealthy")),
    cell_type = sds[[cell_type_col]]
  )

  p_pt <- ggplot(df, aes(UMAP1, UMAP2, colour = pseudotime)) +
    geom_point(size = 0.4, alpha = 0.6) +
    scale_colour_viridis_c(option = "C", name = "Pseudotime") +
    labs(subtitle = "Pseudotime") + theme_classic(base_size = 9)

  p_coh <- ggplot(df, aes(UMAP1, UMAP2, colour = cohort)) +
    geom_point(size = 0.4, alpha = 0.5) +
    scale_colour_manual(values = c(Healthy = "#3182bd", Unhealthy = "#de2d26")) +
    labs(subtitle = "Cohort") + theme_classic(base_size = 9)

  p_ct <- ggplot(df, aes(UMAP1, UMAP2, colour = cell_type)) +
    geom_point(size = 0.4, alpha = 0.6) +
    labs(subtitle = "Cell type") + theme_classic(base_size = 9)

  p_dens <- ggplot(df, aes(x = pseudotime, fill = cohort)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c(Healthy = "#3182bd", Unhealthy = "#de2d26")) +
    labs(subtitle = "Pseudotime density by cohort") +
    theme_classic(base_size = 9)

  # Multi-lineage: one density panel per lineage
  if (ncol(pt_mat) > 1) {
    pt_long <- pt_mat %>%
      as.data.frame() %>%
      mutate(cohort = df$cohort) %>%
      pivot_longer(-cohort, names_to = "lineage", values_to = "pseudotime") %>%
      drop_na(pseudotime)

    p_dens <- ggplot(pt_long, aes(x = pseudotime, fill = cohort)) +
      geom_density(alpha = 0.5) +
      scale_fill_manual(values = c(Healthy = "#3182bd", Unhealthy = "#de2d26")) +
      facet_wrap(~lineage) +
      labs(subtitle = "Pseudotime density by cohort × lineage") +
      theme_classic(base_size = 9)
  }

  p_combined <- (p_pt | p_coh | p_ct) / p_dens +
    plot_annotation(title = title)

  ggsave(paste0(out_prefix, "_trajectory_umap.pdf"), p_combined,
         width = 14, height = 9)
  invisible(p_combined)
}

# Visualise tradeSeq conditionTest results:
#   1. Fitted smoother curves per condition (plotSmoothers) — WHERE along the
#      trajectory does T2D expression diverge from Healthy?
#   2. UMAP coloured by expression, faceted by cohort — spatial context.
#   3. Pseudotime-ordered heatmap — global view of gene dynamics across trajectory.
#   4. Significance overview — ranked genes with FDR labels.
#   5. Start-vs-end waterfall — genes ranked by log2FC start→end, colour = direction.
plot_tradeseq_results <- function(sce_gam, cond_res, sds, cnt,
                                   cell_type_col, title, out_prefix, ct_colors,
                                   sve_res = NULL, n_top = 20,
                                   fontsize_row = 16,
                                   fontsize = 14,
                                   annotation_height = 20) {
  cohort_labels <- c("Healthy" = "Ob", "Unhealthy" = "Ob+T2D")
  coh_cols <- c("#A1A1A1", "#FDC707")
  names(coh_cols) <- cohort_labels[c("Healthy", "Unhealthy")]

  top_genes <- cond_res %>%
    filter(!is.na(pvalue)) %>%
    mutate(padj = p.adjust(pvalue, method = "BH")) %>%
    arrange(pvalue) %>%
    head(n_top) %>%
    pull(gene)

  if (length(top_genes) == 0) {
    cat("  No genes to plot\n"); return(invisible(NULL))
  }
  cat("  Plotting", length(top_genes), "top condition-DE genes\n")

  # # ── 1. Smoother curves (fitted GAM per condition) ─────────────────────────────
  # # Blue = Healthy, Red = Unhealthy. One panel per gene, 12 per page.
  # # curvesCols length must match n_conditions × n_lineages.
  # n_lin    <- length(slingLineages(sds))
  # crv_cols <- rep(c("#3182bd", "#de2d26"), times = n_lin)  # Healthy then Unhealthy per lineage
  # 
  # smoother_plots <- lapply(top_genes, function(g) {
  #   tryCatch(
  #     plotSmoothers(sce_gam, counts = cnt, gene = g, curvesCols = crv_cols) +
  #       ggtitle(g) + theme_classic(base_size = 8) +
  #       theme(legend.position = "none", plot.title = element_text(size = 9, face = "bold")),
  #     error = function(e) NULL
  #   )
  # })
  # smoother_plots <- Filter(Negate(is.null), smoother_plots)
  # 
  # n_per_page <- 12
  # for (i in seq(1, length(smoother_plots), by = n_per_page)) {
  #   chunk <- smoother_plots[i:min(i + n_per_page - 1, length(smoother_plots))]
  #   p <- wrap_plots(chunk, ncol = 4) +
  #     plot_annotation(title    = paste(title, "— fitted expression by pseudotime"),
  #                     subtitle = "Blue = Healthy  |  Red = Unhealthy")
  #   ggsave(paste0(out_prefix, "_smoothers_p", ceiling(i / n_per_page), ".pdf"),
  #          p, width = 16, height = ceiling(length(chunk) / 4) * 4 + 1)
  # }
  # 
  # # ── 2. UMAP coloured by expression, split by cohort ───────────────────────────
  # umap_mat  <- reducedDim(sds, "UMAP")
  # umap_base <- data.frame(UMAP1 = umap_mat[, 1], UMAP2 = umap_mat[, 2],
  #                          cohort = sds$cohort)
  # 
  # umap_plots <- lapply(head(top_genes, 9), function(g) {
  #   if (!g %in% rownames(logcounts(sds))) return(NULL)
  #   df      <- umap_base
  #   df$expr <- as.numeric(logcounts(sds)[g, ])
  #   df$expr <- pmin(df$expr, quantile(df$expr, 0.99))  # cap extreme outliers
  #   ggplot(df, aes(UMAP1, UMAP2, colour = expr)) +
  #     geom_point(size = 0.3, alpha = 0.7) +
  #     scale_colour_gradient(low = "grey92", high = "#b2182b", name = "log\ncounts") +
  #     facet_wrap(~cohort, nrow = 1) +
  #     ggtitle(g) + theme_classic(base_size = 8) +
  #     theme(strip.text = element_text(size = 7),
  #           legend.key.height = unit(0.4, "cm"))
  # })
  # umap_plots <- Filter(Negate(is.null), umap_plots)
  # if (length(umap_plots) > 0) {
  #   p_umap <- wrap_plots(umap_plots, ncol = 3) +
  #     plot_annotation(title = paste(title, "— gene expression on UMAP"))
  #   ggsave(paste0(out_prefix, "_umap_expression.pdf"), p_umap,
  #          width = 18, height = ceiling(length(umap_plots) / 3) * 5)
  # }

  # ── 3. Pseudotime-ordered heatmap ─────────────────────────────────────────────
  # Cells ordered by mean pseudotime; colour = row-scaled logcounts (capped ±3 SD).
  pt_mean  <- rowMeans(slingPseudotime(sds), na.rm = TRUE)
  cell_ord <- order(pt_mean)
  hm_genes <- intersect(top_genes, rownames(logcounts(sds)))

  if (length(hm_genes) > 0) {
    expr_mat    <- as.matrix(logcounts(sds)[hm_genes, cell_ord])
    expr_scaled <- t(scale(t(expr_mat)))
    expr_scaled  <- pmax(pmin(expr_scaled, 3), -3)

    ann_col <- data.frame(
      Cohort     = cohort_labels[sds$cohort[cell_ord]],
      CellType   = sds[[cell_type_col]][cell_ord],
      Pseudotime = pt_mean[cell_ord],
      row.names  = colnames(sds)[cell_ord]
    )
    pdf(paste0(out_prefix, "_pseudotime_heatmap.pdf"), width = 14, height = 6)
    pheatmap(
      expr_scaled,
      cluster_cols       = FALSE,
      cluster_rows       = TRUE,
      show_colnames      = FALSE,
      annotation_col     = ann_col,
      annotation_colors  = list(Cohort = coh_cols, CellType = ct_colors),
      color              = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
      fontsize           = fontsize,
      fontsize_row       = fontsize_row,
      annotation_height  = annotation_height,
      main               = paste(title, "— top condition-DE genes × pseudotime")
    )
    dev.off()
  }

  # # ── 4. Significance overview ───────────────────────────────────────────────────
  # # Ranked by p-value; top 10 genes labelled; FDR < 0.05 highlighted in red.
  # p_sig <- cond_res %>%
  #   filter(!is.na(pvalue)) %>%
  #   mutate(padj = p.adjust(pvalue, method = "BH"),
  #          rank  = row_number(),
  #          label = ifelse(gene %in% head(top_genes, 10), gene, NA_character_)) %>%
  #   ggplot(aes(rank, -log10(pvalue))) +
  #   geom_point(aes(colour = padj < 0.05), size = 0.8, alpha = 0.7) +
  #   geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20) +
  #   scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#b2182b"),
  #                       name = "FDR < 0.05") +
  #   labs(x = "Gene rank", y = expression(-log[10](p)),
  #        title = paste(title, "— conditionTest overview")) +
  #   theme_classic(base_size = 9)
  # ggsave(paste0(out_prefix, "_conditionTest_overview.pdf"), p_sig,
  #        width = 8, height = 5)
  # 
  # # ── 5. Start-vs-end waterfall ──────────────────────────────────────────────────
  # # Top significant genes ranked by |logFC|; bars coloured by direction.
  # # Genes going UP toward the terminal state (e.g. mature adipocytes) are red;
  # # genes going DOWN (enriched in progenitors) are blue.
  # if (!is.null(sve_res)) {
  #   n_bar <- min(40, nrow(sve_res))
  #   p_sve <- sve_res %>%
  #     filter(!is.na(pvalue)) %>%
  #     mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  #     filter(padj < 0.05) %>%
  #     slice_max(abs(logFC_start_to_end), n = n_bar) %>%
  #     mutate(gene = fct_reorder(gene, logFC_start_to_end),
  #            direction = ifelse(logFC_start_to_end > 0, "Up at end", "Up at start")) %>%
  #     ggplot(aes(logFC_start_to_end, gene, fill = direction)) +
  #     geom_col() +
  #     scale_fill_manual(values = c("Up at end" = "#b2182b", "Up at start" = "#2166ac")) +
  #     labs(x = "log2FC (end vs start of trajectory)", y = NULL, fill = NULL,
  #          title = paste(title, "— start vs end of trajectory"),
  #          subtitle = paste("FDR < 0.05 genes, ranked by |log2FC|; n =",
  #                           sum(sve_res$padj < 0.05, na.rm = TRUE), "significant")) +
  #     theme_classic(base_size = 9) +
  #     theme(axis.text.y = element_text(size = 7))
  #   ggsave(paste0(out_prefix, "_startVsEnd_waterfall.pdf"), p_sve,
  #          width = 7, height = max(4, n_bar * 0.22 + 2))
  # }
  # 
  # invisible(NULL)
}

# ── Analysis 1: Omentum — Adipocyte progenitor trajectory ────────────────────────
# Cell types: FAPs (root) → Pre_adipocytes
# Sex regressed out from PCA to reduce its confounding on trajectory shape.
cat("\n=== Analysis 1: Omentum adipocyte progenitor trajectory ===\n")

out1 <- file.path(OUTPUT_DIR, "om_adipocyte")
dir.create(out1, recursive = TRUE, showWarnings = FALSE)

sub_om_adip <- prep_subset(
  sce_om,
  cell_types       = c("FAPs", "Pre_adipocytes"),
  cell_type_col    = "cell_type_",
  regress_sex      = TRUE,
  integrate_cohort = TRUE
)
saveRDS(sub_om_adip, file.path(out1, "sub_sce.rds"))

sds_om_adip <- fit_slingshot(
  sub_om_adip,
  start_clus    = "FAPs",
  end_clus = "Pre_adipocytes",
  cell_type_col = "cell_type_"
)
saveRDS(sds_om_adip, file.path(out1, "slingshot.rds"))

cat("  Slingshot lineages found:", length(slingLineages(sds_om_adip)), "\n")
cat("  Lineage structure:\n"); print(slingLineages(sds_om_adip))

plot_trajectory(sds_om_adip, "cell_type_", "Omentum: adipocyte progenitors",
                file.path(out1, "om_adip"))

cond_om_adip <- factor(sds_om_adip$cohort, levels = c("Healthy", "Unhealthy"))

cm_om_adip <- run_condiments(sds_om_adip, cond_om_adip, nrep = 500,
                              out_prefix = file.path(out1, "om_adip"))

# Sex design matrix for tradeSeq covariate (cells × 1)
sex_U_om_adip <- model.matrix(~ sds_om_adip$sex)[, -1, drop = FALSE]

ts_om_adip <- run_tradeseq(sds_om_adip, cond_om_adip, U = sex_U_om_adip,
                            n_knots = N_KNOTS,
                            out_prefix = file.path(out1, "om_adip"))

plot_tradeseq_results(ts_om_adip$sce_gam, ts_om_adip$condition_res,
                      sds_om_adip, ts_om_adip$cnt,
                      cell_type_col = "cell_type_",
                      title = "Omentum: adipocyte progenitors",
                      out_prefix = file.path(out1, "om_adip"),
                      sve_res = ts_om_adip$sve_res)

# ── Analysis 2: SubQ — Adipocyte progenitor trajectory ──────────────────────────
# NOTE: subq participants are very unequal in cell number
# (SQ_Healthy_3: 80, SQ_Unhealthy_4: 62 total cells).
# Interpret condiments results cautiously; nrep is kept lower for speed.
cat("\n=== Analysis 2: SubQ adipocyte progenitor trajectory ===\n")

out2 <- file.path(OUTPUT_DIR, "sq_adipocyte")
dir.create(out2, recursive = TRUE, showWarnings = FALSE)

# SubQ cell type labels use "Pre-adipocytes" (hyphen), not "Pre_adipocytes"
sub_sq_adip <- prep_subset(
  sce_sq,
  cell_types       = c("FAPs", "Pre_adipocytes"),
  cell_type_col    = "cell_type_",
  regress_sex      = FALSE,   # sex not available in subq
  integrate_cohort = TRUE
)
saveRDS(sub_sq_adip, file.path(out2, "sub_sce.rds"))

sds_sq_adip <- fit_slingshot(
  sub_sq_adip,
  start_clus    = "FAPs",
  end_clus = "Pre_adipocytes",
  cell_type_col = "cell_type_"
)
saveRDS(sds_sq_adip, file.path(out2, "slingshot.rds"))

cat("  Slingshot lineages found:", length(slingLineages(sds_sq_adip)), "\n")

plot_trajectory(sds_sq_adip, "cell_type_", "SubQ: adipocyte progenitors",
                file.path(out2, "sq_adip"))

cond_sq_adip <- factor(sds_sq_adip$cohort, levels = c("Healthy", "Unhealthy"))

cm_sq_adip <- run_condiments(sds_sq_adip, cond_sq_adip, nrep = 500,
                              out_prefix = file.path(out2, "sq_adip"))


# Sex design matrix for tradeSeq covariate (cells × 1)
sex_U_om_adip <- model.matrix(~ sds_sq_adip$sex)[, -1, drop = FALSE]

ts_sq_adip <- run_tradeseq(sds_sq_adip, cond_sq_adip, U = sex_U_om_adip,
                            n_knots = N_KNOTS,
                            out_prefix = file.path(out2, "sq_adip"))

plot_tradeseq_results(ts_sq_adip$sce_gam, ts_sq_adip$condition_res,
                      sds_sq_adip, ts_sq_adip$cnt,
                      cell_type_col = "cell_type_",
                      title = "SubQ: adipocyte progenitors",
                      out_prefix = file.path(out2, "sq_adip"),
                      sve_res = ts_sq_adip$sve_res)

# ── Analysis 3: Omentum — Mesothelial transition ─────────────────────────────────
# Mesothelial_cells (root) → MMT_cells (mesothelial-to-mesenchymal transition)
#   + Myofibroblasts as a second potential endpoint.
#
# IMPORTANT CAVEAT: MMT_cells = 37 cells. Slingshot will fit a trajectory but
# power for condiments/tradeSeq tests is very low. Treat results as exploratory.
# Consider also running on broad cell_type_ which includes Mesothelial (5009 cells)
# together with Myofibroblasts (697) as a more powered alternative.
cat("\n=== Analysis 3: Omentum mesothelial transition ===\n")

out3 <- file.path(OUTPUT_DIR, "om_mesothelial")
dir.create(out3, recursive = TRUE, showWarnings = FALSE)

sub_om_meso <- prep_subset(
  sce_om,
  cell_types       = c("Mesothelial", "IGFBP2_cells"),
  cell_type_col    = "cell_type_",
  regress_sex      = TRUE,
  integrate_cohort = TRUE
)
saveRDS(sub_om_meso, file.path(out3, "sub_sce.rds"))

sds_om_meso <- fit_slingshot(
  sub_om_meso,
  start_clus    = "Mesothelial",
  end_clus      = "IGFBP2_cells",
  cell_type_col = "cell_type_"
)
saveRDS(sds_om_meso, file.path(out3, "slingshot.rds"))

cat("  Slingshot lineages found:", length(slingLineages(sds_om_meso)), "\n")
cat("  Lineage structure:\n"); print(slingLineages(sds_om_meso))

plot_trajectory(sds_om_meso, "cell_type_", "Omentum: mesothelial transition",
                file.path(out3, "om_meso"))

cond_om_meso <- factor(sds_om_meso$cohort, levels = c("Healthy", "Unhealthy"))

cm_om_meso <- run_condiments(sds_om_meso, cond_om_meso, nrep = 500,
                              out_prefix = file.path(out3, "om_meso"))

sex_U_om_meso <- model.matrix(~ sds_om_meso$sex)[, -1, drop = FALSE]

ts_om_meso <- run_tradeseq(sds_om_meso, cond_om_meso, U = sex_U_om_meso,
                            n_knots = N_KNOTS,
                            out_prefix = file.path(out3, "om_meso"))

plot_tradeseq_results(ts_om_meso$sce_gam, ts_om_meso$condition_res,
                      sds_om_meso, ts_om_meso$cnt,
                      cell_type_col = "cell_type2_",
                      title = "Omentum: mesothelial transition",
                      out_prefix = file.path(out3, "om_meso"),
                      sve_res = ts_om_meso$sve_res
                      )

# ── Analysis 4: Omentum — Extended mesothelial trajectory (4 cell types) ─────────
#
# WHY this is different from Analysis 3:
#   Analysis 3 used only Mesothelial + IGFBP2_cells. The problem: mesothelial
#   cells differ hugely between Healthy and Unhealthy, creating alternating
#   "waves" of condition-enriched pseudotime regions. All tradeSeq conditionTest
#   hits were dominated by mesothelial-cell differences rather than trajectory
#   dynamics downstream.
#
# TWO FIXES implemented here:
#
#   Fix A — downstream-only tradeSeq:
#     Slingshot is fit on all 4 cell types (mesothelial defines the root and
#     sets the pseudotime scale). tradeSeq is then fit ONLY on the non-mesothelial
#     cells (IGFBP2_cells, FAPs, Pre_adipocytes). conditionTest now tests for
#     condition differences in the downstream trajectory only, completely bypassing
#     the mesothelial "wave" problem.
#
#   Fix B — separate trajectory per condition + smoother comparison:
#     Slingshot is fit independently on Healthy and Unhealthy cells. Each
#     pseudotime is rank-normalised to [0, 1]. tradeSeq GAMs are fit separately
#     per condition (downstream cells only). predictSmooth is used to get fitted
#     expression curves on a common pseudotime grid; the max absolute difference
#     between curves identifies genes with the most divergent dynamics.
#     associationTest per condition identifies genes that are trajectory-associated
#     in one condition but not the other.
#
# Trajectory direction: Mesothelial (root) → IGFBP2_cells → FAPs → Pre_adipocytes
#
cat("\n=== Analysis 4: Omentum extended mesothelial trajectory (4 cell types) ===\n")

CELL_TYPES_4CT <- c("Mesothelial", "IGFBP2_cells", "FAPs", "Pre_adipocytes")

out4 <- file.path(OUTPUT_DIR, "om_meso_extended")
dir.create(out4, recursive = TRUE, showWarnings = FALSE)

# ── 4.1  Prep — Harmony on participant + cohort (aligns conditions for trajectory)
sub_4ct <- prep_subset(
  sce_om,
  cell_types       = CELL_TYPES_4CT,
  cell_type_col    = "cell_type_",
  regress_sex      = TRUE,
  integrate_cohort = TRUE    # Harmony on cohort + sex → shared trajectory embedding
)
saveRDS(sub_4ct, file.path(out4, "sub_sce.rds"))

# sub_4ct <- readRDS(file.path(out4, "sub_sce.rds"))

cat("  Cell counts:\n")
print(table(sub_4ct$cell_type_, sub_4ct$cohort))

# ── 4.2  Slingshot on combined data (all 4 cell types)
sds_4ct <- fit_slingshot(
  sub_4ct,
  start_clus    = "Mesothelial",
  end_clus      = "Pre_adipocytes",
  cell_type_col = "cell_type_"
)
saveRDS(sds_4ct, file.path(out4, "slingshot_4ct.rds"))
sds_4ct <- readRDS(file.path(out4, "slingshot_4ct.rds"))

cat("  Lineages:\n"); print(slingLineages(sds_4ct))
plot_trajectory(sds_4ct, "cell_type_", "Omentum: extended mesothelial (4ct)",
                file.path(out4, "om_meso_ext"))

cond_4ct <- factor(sds_4ct$cohort, levels = c("Healthy", "Unhealthy"))

# condiments: does progression differ between conditions along the full trajectory?
cm_4ct <- run_condiments(sds_4ct, cond_4ct, nrep = 500,
                          out_prefix = file.path(out4, "om_meso_ext"))

# running tradeSeq on all 4 cell-types

sex_U_4ct <- model.matrix(~ sds_4ct$sex)[, -1, drop = FALSE]

ts_om_4ct <- run_tradeseq(sds_4ct, cond_4ct, U = sex_U_4ct,
                           n_knots = N_KNOTS,
                           out_prefix = file.path(out4, "om_meso_ext"))

ct_colors <- c(
  "IGFBP2_cells" = "#FD1486", # pink
  "Mesothelial" = "#00546B", # dark green
  "Pre_adipocytes" = "#941DEC", # violet
  "FAPs" = "#6EB8E7" # light blue
)

plot_tradeseq_results(ts_om_4ct$sce_gam, ts_om_4ct$condition_res,
                      sds_4ct, ts_om_4ct$cnt,
                      cell_type_col = "cell_type_",
                      title = "Omentum: mesothelial-FAPs transition",
                      out_prefix = file.path(out4, "om_meso_4ct"),
                      sve_res = ts_om_4ct$sve_res,
                      ct_colors = ct_colors
                      )

# ── 4.3  Fix A: downstream-only tradeSeq ──────────────────────────────────────────
# Key idea: Slingshot pseudotime comes from the 4-ct trajectory (mesothelial sets
# the root), but the GAM is fit ONLY on non-mesothelial cells. This prevents the
# large mesothelial H/U difference from dominating conditionTest.

cat("\n--- Fix A: downstream-only tradeSeq (IGFBP2_cells + FAPs + Pre_adipocytes) ---\n")

not_meso <- sds_4ct$cell_type_ != "Mesothelial"
cat("  Downstream cells:", sum(not_meso), "\n")
cat(table(sds_4ct$cell_type_[not_meso], sds_4ct$cohort[not_meso]))

# Pseudotime matrix for downstream cells only
pt_downstream      <- slingPseudotime(sds_4ct)[not_meso, , drop = FALSE]
weights_downstream <- slingCurveWeights(sds_4ct)[not_meso, , drop = FALSE]
cond_downstream    <- factor(sds_4ct$cohort[not_meso], levels = c("Healthy", "Unhealthy"))

# Build count matrix for downstream cells (HVGs ∩ expressed)
hvgs_4ct   <- metadata(sds_4ct)$hvgs
cnt_all_4ct <- counts(sds_4ct)
expressed_4ct <- rowSums(cnt_all_4ct > 0) >= 0.05 * ncol(cnt_all_4ct)
keep_genes  <- intersect(hvgs_4ct, rownames(cnt_all_4ct)[expressed_4ct])
cnt_downstream <- cnt_all_4ct[keep_genes, not_meso]

sex_U_downstream <- model.matrix(~ sds_4ct$sex[not_meso])[, -1, drop = FALSE]

cat("  Genes for downstream GAM:", nrow(cnt_downstream), "\n")
cat("  Fitting GAM on downstream cells only...\n")

sce_gam_ds <- fitGAM(
  counts      = cnt_downstream,
  pseudotime  = pt_downstream,
  cellWeights = weights_downstream,
  conditions  = cond_downstream,
  U           = sex_U_downstream,
  nknots      = N_KNOTS,
  parallel    = TRUE,
  BPPARAM     = SnowParam(N_CORES, type = "SOCK"),
  verbose     = FALSE
)
saveRDS(sce_gam_ds, file.path(out4, "sce_gam_downstream.rds"))

cond_res_ds <- conditionTest(sce_gam_ds, l2fc = L2FC_THRESH) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  arrange(pvalue)

assoc_res_ds <- associationTest(sce_gam_ds, global = TRUE, lineages = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  arrange(pvalue)

sve_res_ds <- startVsEndTest(sce_gam_ds, global = TRUE, lineages = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(padj = p.adjust(pvalue, method = "BH"))

pred2_ds <- predictSmooth(sce_gam_ds, gene = sve_res_ds$gene, nPoints = 2, tidy = FALSE)
start_cols_ds <- seq(1, ncol(pred2_ds), by = 2)
end_cols_ds   <- seq(2, ncol(pred2_ds), by = 2)
logfc_ds      <- rowMeans(log2(pred2_ds[, end_cols_ds, drop = FALSE] + 1) -
                           log2(pred2_ds[, start_cols_ds, drop = FALSE] + 1))
sve_res_ds$logFC_start_to_end <- logfc_ds[sve_res_ds$gene]
sve_res_ds <- arrange(sve_res_ds, pvalue)

write_csv(cond_res_ds,  file.path(out4, "downstream_conditionTest.csv"))
write_csv(assoc_res_ds, file.path(out4, "downstream_associationTest.csv"))
write_csv(sve_res_ds,   file.path(out4, "downstream_startVsEndTest.csv"))

cat("  Top condition-differential genes (downstream only):\n")
print(head(cond_res_ds %>% filter(!is.na(pvalue)) %>%
             mutate(padj = p.adjust(pvalue, "BH")), 20))

# Build a temporary SCE that carries the downstream cells for plotting helpers
sds_ds_plot <- sds_4ct[, not_meso]
plot_tradeseq_results(
  sce_gam_ds, cond_res_ds, sds_ds_plot, cnt_downstream,
  cell_type_col = "cell_type_",
  title = "Omentum meso extended: downstream tradeSeq",
  out_prefix = file.path(out4, "downstream"),
  sve_res = sve_res_ds
)

# ── 4.4  Fix B: separate trajectory per condition + smoother comparison ───────────
# Fit Slingshot independently on H and U cells, then compare fitted expression
# curves after rank-normalising pseudotime to [0,1].
#
# This is NOT a formal statistical test; it is an exploratory ranking of genes
# whose dynamics differ most between conditions. Use as a discovery tool and
# validate hits with Fix A conditionTest or external data.

cat("\n--- Fix B: separate trajectory per condition ---\n")

# Helper: [0, 1] min-max normalisation of pseudotime (per lineage column)
normalise_pt <- function(pt_mat) {
  apply(pt_mat, 2, function(x) {
    valid <- !is.na(x)
    out   <- rep(NA_real_, length(x))
    rng   <- range(x[valid])
    if (diff(rng) == 0) return(out)
    out[valid] <- (x[valid] - rng[1]) / diff(rng)
    out
  })
}

fit_separate <- function(sce_full, cohort_label, cell_types,
                          cell_type_col = "cell_type_",
                          start_clus = "Mesothelial",
                          n_hvg = N_HVG, n_pc = N_PC) {
  cat("  Fitting separate trajectory for", cohort_label, "...\n")
  sub <- sce_full[, sce_full$cohort == cohort_label &
                    sce_full[[cell_type_col]] %in% cell_types]
  cat("  Cells:", ncol(sub), "\n")

  sub <- logNormCounts(sub, assay.type = "counts")
  dec  <- modelGeneVar(sub)
  hvgs <- getTopHVGs(dec, n = n_hvg)
  sub  <- runPCA(sub, subset_row = hvgs, ncomponents = n_pc, name = "PCA")

  # Harmony on participant only (single condition, no cohort correction needed)
  n_participants <- length(unique(sub$common_participant))
  if (n_participants > 1) {
    sub <- tryCatch(
      RunHarmony(sub, group.by.vars = "common_participant",
                 reduction.save = "HARMONY", verbose = FALSE),
      error = function(e) {
        cat("  Harmony failed, using PCA\n")
        reducedDim(sub, "HARMONY") <- reducedDim(sub, "PCA")
        sub
      }
    )
  } else {
    reducedDim(sub, "HARMONY") <- reducedDim(sub, "PCA")
  }

  sub <- runUMAP(sub, dimred = "HARMONY", n_dimred = min(n_pc, 20), name = "UMAP")
  metadata(sub)$hvgs     <- hvgs
  metadata(sub)$traj_dimred <- "HARMONY"

  sds <- slingshot(sub,
                   clusterLabels = cell_type_col,
                   reducedDim    = "HARMONY",
                   start.clus    = start_clus,
                   approx_points = 200)
  list(sds = sds, hvgs = hvgs)
}

fit_H <- fit_separate(sce_om, "Healthy",   CELL_TYPES_4CT)
fit_U <- fit_separate(sce_om, "Unhealthy", CELL_TYPES_4CT)

sds_H <- fit_H$sds
sds_U <- fit_U$sds

saveRDS(sds_H, file.path(out4, "slingshot_Healthy.rds"))
saveRDS(sds_U, file.path(out4, "slingshot_Unhealthy.rds"))

# Rank-normalise pseudotime to [0,1] for comparability
pt_H_norm <- normalise_pt(slingPseudotime(sds_H))
pt_U_norm <- normalise_pt(slingPseudotime(sds_U))

cat("  Healthy pseudotime range after normalisation:",
    range(pt_H_norm, na.rm = TRUE), "\n")
cat("  Unhealthy pseudotime range after normalisation:",
    range(pt_U_norm, na.rm = TRUE), "\n")

# Downstream cells only (exclude mesothelial) for each condition
not_meso_H <- sds_H$cell_type_ != "Mesothelial"
not_meso_U <- sds_U$cell_type_ != "Mesothelial"

# Shared HVGs expressed in ≥5% of cells in EACH condition
make_keep <- function(sds, hvgs, not_meso_mask) {
  cnt <- counts(sds)[, not_meso_mask]
  expressed <- rowSums(cnt > 0) >= 0.05 * ncol(cnt)
  intersect(hvgs, rownames(cnt)[expressed])
}
keep_H <- make_keep(sds_H, fit_H$hvgs, not_meso_H)
keep_U <- make_keep(sds_U, fit_U$hvgs, not_meso_U)
keep_sep <- intersect(keep_H, keep_U)
cat("  Genes shared between conditions for separate GAMs:", length(keep_sep), "\n")

cat("  Fitting Healthy downstream GAM...\n")
sce_gam_H <- fitGAM(
  counts      = counts(sds_H)[keep_sep, not_meso_H],
  pseudotime  = pt_H_norm[not_meso_H, , drop = FALSE],
  cellWeights = slingCurveWeights(sds_H)[not_meso_H, , drop = FALSE],
  nknots      = N_KNOTS,
  parallel    = TRUE,
  BPPARAM     = SnowParam(N_CORES, type = "SOCK"),
  verbose     = FALSE
)

cat("  Fitting Unhealthy downstream GAM...\n")
sce_gam_U <- fitGAM(
  counts      = counts(sds_U)[keep_sep, not_meso_U],
  pseudotime  = pt_U_norm[not_meso_U, , drop = FALSE],
  cellWeights = slingCurveWeights(sds_U)[not_meso_U, , drop = FALSE],
  nknots      = N_KNOTS,
  parallel    = TRUE,
  BPPARAM     = SnowParam(N_CORES, type = "SOCK"),
  verbose     = FALSE
)

saveRDS(sce_gam_H, file.path(out4, "sce_gam_Healthy_downstream.rds"))
saveRDS(sce_gam_U, file.path(out4, "sce_gam_Unhealthy_downstream.rds"))

# Compare fitted smoothers on a common [0,1] pseudotime grid
N_SMOOTH_PTS <- 100
genes_sep <- intersect(rownames(sce_gam_H), rownames(sce_gam_U))
cat("  Comparing smoothers for", length(genes_sep), "genes ...\n")

smooth_H_mat <- predictSmooth(sce_gam_H, gene = genes_sep,
                               nPoints = N_SMOOTH_PTS, tidy = FALSE)
smooth_U_mat <- predictSmooth(sce_gam_U, gene = genes_sep,
                               nPoints = N_SMOOTH_PTS, tidy = FALSE)

# If multiple lineages, smooth_H_mat has ncol = n_lin * N_SMOOTH_PTS; average lineages
n_lin_H <- length(slingLineages(sds_H))
n_lin_U <- length(slingLineages(sds_U))

collapse_lineages <- function(mat, n_lin, n_pts) {
  if (n_lin == 1) return(mat)
  # Stack lineages: columns are [lin1_t1...lin1_tN, lin2_t1...lin2_tN, ...]
  lin_mats <- lapply(seq_len(n_lin), function(i) {
    mat[, ((i - 1) * n_pts + 1):(i * n_pts)]
  })
  Reduce("+", lin_mats) / n_lin  # gene × n_pts
}

sm_H <- collapse_lineages(smooth_H_mat, n_lin_H, N_SMOOTH_PTS)
sm_U <- collapse_lineages(smooth_U_mat, n_lin_U, N_SMOOTH_PTS)

# Gene-level summary: max |difference| and AUC of |difference| between curves
diff_abs   <- abs(sm_H - sm_U)
max_diff   <- apply(diff_abs, 1, max)
auc_diff   <- rowSums(diff_abs) / N_SMOOTH_PTS
cor_curves <- diag(cor(t(sm_H), t(sm_U)))  # per-gene correlation between curves

smoother_comparison <- data.frame(
  gene       = genes_sep,
  max_abs_diff = max_diff,
  auc_abs_diff = auc_diff,
  curve_cor    = cor_curves
) %>% arrange(desc(max_abs_diff))

write_csv(smoother_comparison, file.path(out4, "separate_smoother_comparison.csv"))
cat("  Top genes with most divergent dynamics (Healthy vs Unhealthy):\n")
print(head(smoother_comparison, 20))

# associationTest per condition: trajectory-associated genes in H but not U
assoc_H <- associationTest(sce_gam_H, global = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(padj_H = p.adjust(pvalue, "BH")) %>%
  arrange(pvalue)

assoc_U <- associationTest(sce_gam_U, global = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(padj_U = p.adjust(pvalue, "BH")) %>%
  arrange(pvalue)

write_csv(assoc_H, file.path(out4, "separate_associationTest_Healthy.csv"))
write_csv(assoc_U, file.path(out4, "separate_associationTest_Unhealthy.csv"))

# Merge: genes significant in one condition but not the other
assoc_merged <- assoc_H %>%
  select(gene, padj_H) %>%
  full_join(assoc_U %>% select(gene, padj_U), by = "gene") %>%
  left_join(smoother_comparison %>% select(gene, max_abs_diff, curve_cor), by = "gene") %>%
  mutate(
    sig_H = !is.na(padj_H) & padj_H < 0.05,
    sig_U = !is.na(padj_U) & padj_U < 0.05,
    category = case_when(
      sig_H & !sig_U ~ "Healthy-specific",
      sig_U & !sig_H ~ "Unhealthy-specific",
      sig_H &  sig_U ~ "Both",
      TRUE           ~ "Neither"
    )
  ) %>%
  arrange(desc(max_abs_diff))

write_csv(assoc_merged, file.path(out4, "separate_association_comparison.csv"))
cat("  Condition-specific trajectory genes:\n")
print(table(assoc_merged$category))

# Plot: top divergent genes — smoother curves for H and U overlaid
top_div_genes <- head(smoother_comparison$gene, 16)
pt_grid <- seq(0, 1, length.out = N_SMOOTH_PTS)

plot_list <- lapply(top_div_genes, function(g) {
  df <- data.frame(
    pseudotime = rep(pt_grid, 2),
    expr       = c(sm_H[g, ], sm_U[g, ]),
    cohort     = rep(c("Healthy", "Unhealthy"), each = N_SMOOTH_PTS)
  )
  ggplot(df, aes(pseudotime, expr, colour = cohort)) +
    geom_line(lw = 1.2) +
    scale_colour_manual(values = c(Healthy = "#3182bd", Unhealthy = "#de2d26")) +
    labs(title = g, x = "Normalised pseudotime [0,1]", y = "Fitted expression") +
    theme_classic(base_size = 8) +
    theme(legend.position = "none",
          plot.title = element_text(size = 9, face = "bold"))
})

p_div <- wrap_plots(plot_list, ncol = 4) +
  plot_annotation(
    title    = "Top divergent genes: separate H vs U smoothers",
    subtitle = "Blue = Healthy | Red = Unhealthy | x-axis = rank-normalised pseudotime [0,1]"
  )
ggsave(file.path(out4, "separate_top_divergent_smoothers.pdf"), p_div,
       width = 16, height = ceiling(length(plot_list) / 4) * 4 + 1)

# ── Analysis 4.5: Milo — Differential Neighbourhood Abundance ────────────────────
#
# Milo does not use pseudotime: it builds a KNN graph, samples representative
# neighbourhoods, counts cells per sample per neighbourhood, then fits an edgeR
# NB GLM per neighbourhood to test whether it is over/under-represented in T2D.
#
# EMBEDDING CHOICE — protect-cohort PCA:
#   Participants are nested within cohort, so Harmony on participant also pulls
#   conditions together (removing cohort signal as a side effect). Instead we use
#   limma::removeBatchEffect(batch = participant, design = ~cohort) applied to
#   the raw PCA coordinates. This removes participant batch effects while explicitly
#   preserving the cohort signal — the correct embedding for a DA test.
#
# INSTALL (if needed):
#   BiocManager::install(c("miloR", "edgeR"))

suppressPackageStartupMessages({
  library(miloR)
  library(edgeR)
})

cat("\n=== Analysis 4.5: Milo DA (omentum 4-cell-type mesothelial) ===\n")

out_milo <- file.path(OUTPUT_DIR, "om_meso_milo")
dir.create(out_milo, recursive = TRUE, showWarnings = FALSE)

# ── 4.5.1  Protect-cohort PCA ──────────────────────────────────────────────────
# sub_4ct was produced by prep_subset with integrate_cohort=TRUE, so it carries
# reducedDim "PCA" (raw) and "HARMONY" (harmony-corrected).
# We apply removeBatchEffect to the raw PCA, protecting cohort (and sex if available).

pca_raw <- reducedDim(sub_4ct, "PCA")
n_pcs   <- ncol(pca_raw)

protect_design <- if (!all(is.na(sub_4ct$sex))) {
  model.matrix(~ cohort + sex, data = as.data.frame(colData(sub_4ct)))
} else {
  model.matrix(~ cohort, data = as.data.frame(colData(sub_4ct)))
}

pca_protected <- t(removeBatchEffect(
  t(pca_raw),
  batch  = sub_4ct$common_participant,
  design = protect_design
))
reducedDim(sub_4ct, "PCA_PROTECTED") <- pca_protected
cat("  Protect-cohort PCA computed:", nrow(pca_protected), "cells ×", ncol(pca_protected), "PCs\n")

# ── UMAP sanity check: does protect-cohort embedding look right? ───────────────
# Goal: participants mixed within each cell type (batch removed),
#       Healthy and Unhealthy still visibly separated (cohort preserved).
# Compare against the HARMONY embedding (which removes both participant AND cohort).

sub_4ct <- runUMAP(sub_4ct, dimred = "PCA_PROTECTED",
                   n_dimred = min(n_pcs, 20), name = "UMAP_PROTECTED")

ct_cols4 <- c(Mesothelial    = "#4393c3",
              IGFBP2_cells   = "#92c5de",
              FAPs           = "#f4a582",
              Pre_adipocytes = "#d6604d")
coh_cols <- c(Healthy = "#3182bd", Unhealthy = "#de2d26")

make_umap_df <- function(sce, dimred) {
  umap <- reducedDim(sce, dimred)
  data.frame(
    UMAP1       = umap[, 1],
    UMAP2       = umap[, 2],
    cell_type   = sce$cell_type_,
    cohort      = sce$cohort,
    participant = sce$common_participant
  )
}

df_prot <- make_umap_df(sub_4ct, "UMAP_PROTECTED") %>%
  mutate(embedding = "Protect-cohort\n(participant removed, cohort kept)")
df_harm <- make_umap_df(sub_4ct, "UMAP") %>%
  mutate(embedding = "HARMONY\n(participant + cohort removed)")

umap_compare <- bind_rows(df_prot, df_harm) %>%
  mutate(embedding = factor(embedding, levels = c(
    "Protect-cohort\n(participant removed, cohort kept)",
    "HARMONY\n(participant + cohort removed)"
  )))

p_ct <- ggplot(umap_compare, aes(UMAP1, UMAP2, colour = cell_type)) +
  geom_point(size = 0.3, alpha = 0.4) +
  scale_colour_manual(values = ct_cols4, name = "Cell type") +
  facet_wrap(~embedding, scales = "free") +
  labs(title = "Cell type") +
  theme_classic(base_size = 9) +
  theme(strip.text = element_text(face = "bold"),
        legend.key.size = unit(0.4, "cm"))

p_coh <- ggplot(umap_compare, aes(UMAP1, UMAP2, colour = cohort)) +
  geom_point(size = 0.3, alpha = 0.4) +
  scale_colour_manual(values = coh_cols, name = "Cohort") +
  facet_wrap(~embedding, scales = "free") +
  labs(title = "Cohort") +
  theme_classic(base_size = 9) +
  theme(strip.text = element_text(face = "bold"),
        legend.key.size = unit(0.4, "cm"))

p_part <- ggplot(umap_compare, aes(UMAP1, UMAP2, colour = participant)) +
  geom_point(size = 0.3, alpha = 0.4) +
  facet_wrap(~embedding, scales = "free") +
  labs(title = "Participant") +
  theme_classic(base_size = 9) +
  theme(strip.text = element_text(face = "bold"),
        legend.key.size = unit(0.4, "cm"))

p_embedding_check <- (p_ct / p_coh / p_part) +
  plot_annotation(
    title    = "Embedding comparison: protect-cohort vs HARMONY",
    subtitle = paste(
      "Protect-cohort: participants should be mixed within cell types;",
      "H/U should remain separated.",
      "\nHARMONY: both participant AND cohort effects removed (conditions aligned).",
      "\nUse protect-cohort for Milo; use HARMONY for trajectory/DPT."
    )
  )

ggsave(file.path(out_milo, "embedding_comparison_umap.pdf"),
       p_embedding_check, width = 12, height = 14)
cat("  Embedding comparison UMAP saved.\n")

# ── 4.5.2  Build Milo object and KNN graph ─────────────────────────────────────
milo_obj <- Milo(sub_4ct)

K_MILO <- 15
D_MILO <- min(n_pcs, 20)

milo_obj <- buildGraph(milo_obj, k = K_MILO, d = D_MILO,
                        reduced.dim = "PCA_PROTECTED")

# ── 4.5.3  Sample representative neighbourhoods ────────────────────────────────
milo_obj <- makeNhoods(milo_obj, prop = 0.2, k = K_MILO, d = D_MILO,
                        refined = TRUE, reduced_dims = "PCA_PROTECTED")
cat("  Neighbourhoods sampled:", ncol(nhoods(milo_obj)), "\n")

# ── 4.5.4  Count cells per neighbourhood per sample ───────────────────────────
milo_obj <- countCells(milo_obj,
                        meta.data = as.data.frame(colData(milo_obj)),
                        sample    = "common_participant")

# ── 4.5.5  Sample-level design table ──────────────────────────────────────────
sample_meta <- colData(milo_obj) %>%
  as.data.frame() %>%
  select(common_participant, cohort, sex) %>%
  distinct() %>%
  mutate(cohort = factor(cohort, levels = c("Healthy", "Unhealthy"))) %>%
  arrange(common_participant)
rownames(sample_meta) <- sample_meta$common_participant

cat("  Sample metadata:\n"); print(sample_meta)

# ── 4.5.6  DA test (edgeR NB GLM per neighbourhood) ──────────────────────────
# fdr.weighting="graph-overlap" applies the spatial FDR correction that accounts
# for neighbourhood overlap (Dann et al. 2022).
milo_design <- if ("sex" %in% colnames(sample_meta) && !all(is.na(sample_meta$sex))) {
  ~ cohort + sex
} else {
  ~ cohort
}
cat("  Milo design:", deparse(milo_design), "\n")

milo_res <- testNhoods(
  milo_obj,
  design          = milo_design,
  design.df       = sample_meta,
  model.contrasts = "cohortUnhealthy",   # Unhealthy vs Healthy
  fdr.weighting   = "graph-overlap"
)
spatialFDR_threshold <- 0.2
cat("  Significant DA neighbourhoods (SpatialFDR <spatialFDR_threshold):",
    sum(milo_res$SpatialFDR < spatialFDR_threshold, na.rm = TRUE), "\n")
cat("  Unhealthy-enriched:", sum(milo_res$SpatialFDR < spatialFDR_threshold & milo_res$logFC > 0, na.rm=TRUE), "\n")
cat("  Healthy-enriched:",   sum(milo_res$SpatialFDR < spatialFDR_threshold & milo_res$logFC < 0, na.rm=TRUE), "\n")

# ── 4.5.7  Annotate neighbourhoods with dominant cell type ─────────────────────
milo_obj  <- buildNhoodGraph(milo_obj)
milo_res  <- annotateNhoods(milo_obj, milo_res, coldata_col = "cell_type_")

# Mixed neighbourhoods: fraction < 0.7 of dominant cell type → label "Mixed"
milo_res$cell_type_clean <- ifelse(
  milo_res$cell_type__fraction < 0.7, "Mixed", milo_res$cell_type_
)

cat("  DA neighbourhood summary by cell type:\n")
print(
  milo_res %>%
    group_by(cell_type_clean) %>%
    dplyr::summarise(
      n_nhoods      = dplyr::n(),
      n_sig         = sum(SpatialFDR < spatialFDR_threshold, na.rm = TRUE),
      n_enriched_U  = sum(SpatialFDR < spatialFDR_threshold & logFC > 0, na.rm = TRUE),
      n_enriched_H  = sum(SpatialFDR < spatialFDR_threshold & logFC < 0, na.rm = TRUE),
      .groups = "drop"
    )
)

write_csv(milo_res, file.path(out_milo, "milo_nhood_results.csv"))

# ── 4.5.8  Plots ───────────────────────────────────────────────────────────────
coh_cols <- c(Healthy = "#3182bd", Unhealthy = "#de2d26")

# Beeswarm: logFC per neighbourhood, grouped by cell type
p_bee <- plotDAbeeswarm(milo_res, group.by = "cell_type_clean", alpha = spatialFDR_threshold) +
  scale_colour_gradient2(low = "#3182bd", mid = "white", high = "#de2d26",
                         midpoint = 0, name = "logFC\n(U vs H)") +
  labs(title = "Milo: DA beeswarm by cell type",
       subtitle = "Each point = one neighbourhood; positive = Unhealthy-enriched") +
  theme_classic(base_size = 10)

# Neighbourhood graph on UMAP: nodes = neighbourhoods, coloured by logFC
p_graph <- plotNhoodGraphDA(milo_obj, milo_res, alpha = spatialFDR_threshold) +
  scale_fill_gradient2(low = "#3182bd", mid = "white", high = "#de2d26",
                       midpoint = 0, name = "logFC\n(U vs H)") +
  labs(title = "Milo: DA on UMAP (protect-cohort embedding)") +
  theme_classic(base_size = 10)

# Volcano: logFC vs -log10(SpatialFDR), coloured by cell type
ct_cols <- c(Mesothelial    = "#4393c3",
             IGFBP2_cells   = "#92c5de",
             FAPs           = "#f4a582",
             Pre_adipocytes = "#d6604d",
             Mixed          = "#999999")

p_vol <- milo_res %>%
  mutate(sig = SpatialFDR < spatialFDR_threshold,
         label = ifelse(sig & abs(logFC) > quantile(abs(logFC[sig]), 0.75, na.rm=TRUE),
                        cell_type_clean, NA_character_)) %>%
  ggplot(aes(logFC, -log10(SpatialFDR + 1e-10),
             colour = cell_type_clean, shape = sig)) +
  geom_point(size = 1.2, alpha = 0.7) +
  geom_hline(yintercept = -log10(spatialFDR_threshold), lty = 2, lw = 0.7) +
  geom_vline(xintercept = 0, lty = 3, colour = "grey50") +
  scale_colour_manual(values = ct_cols, name = "Cell type") +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), guide = "none") +
  labs(x = "logFC (Unhealthy vs Healthy)",
       y = expression(-log[10](SpatialFDR)),
       title = "Milo: neighbourhood DA volcano") +
  theme_classic(base_size = 10)

p_milo <- (p_bee | p_graph) / p_vol +
  plot_annotation(title = "Milo — Differential Neighbourhood Abundance",
                  subtitle = "4-cell-type mesothelial trajectory, omentum")

ggsave(file.path(out_milo, "milo_results.pdf"), p_milo,
       width = 14, height = 12)

# Is the mean neighbourhood logFC different from 0?
t.test(milo_res$logFC, mu = 0)

# Per cell type: is there a directional bias?
milo_res %>%
  group_by(cell_type_clean) %>%
  summarise(
    mean_logFC = mean(logFC, na.rm = TRUE),
    t_p        = t.test(logFC, mu = 0)$p.value,
    n_nhoods   = dplyr::n()
  )

cat("  Milo outputs saved to:", out_milo, "\n")

# ── Analysis 5: Monocle3 (install instructions included) ──────────────────────────
# Monocle3 provides an independent trajectory algorithm (learn_graph / DDRTree)
# and fit_models for GAM-based condition × pseudotime interaction tests.
#
# INSTALL (if not already):
#   if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("BiocGenerics", "DelayedArray", "DelayedMatrixStats",
#                          "limma", "lme4", "S4Vectors", "SingleCellExperiment",
#                          "SummarizedExperiment", "batchelor", "HDF5Array",
#                          "terra", "ggrastr"))
#   if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")
#   devtools::install_github("cole-trapnell-lab/monocle3")

if (requireNamespace("monocle3", quietly = TRUE)) {
  suppressPackageStartupMessages(library(monocle3))

  cat("\n=== Analysis 5: Monocle3 (omentum 4-cell-type mesothelial trajectory) ===\n")

  out5 <- file.path(OUTPUT_DIR, "om_meso_monocle3")
  dir.create(out5, recursive = TRUE, showWarnings = FALSE)

  # Build CDS from the prepped SCE (sub_4ct already has Harmony embedding)
  gene_meta <- as.data.frame(rowData(sub_4ct))
  if (!"gene_short_name" %in% colnames(gene_meta))
    gene_meta$gene_short_name <- rownames(gene_meta)

  cds <- new_cell_data_set(
    expression_data = counts(sub_4ct),
    cell_metadata   = as.data.frame(colData(sub_4ct)),
    gene_metadata   = gene_meta
  )

  # Inject existing Harmony UMAP (avoids re-running Monocle3's own UMAP)
  reducedDims(cds)[["UMAP"]] <- reducedDim(sub_4ct, "UMAP")

  cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = 1e-3)
  cds <- learn_graph(cds, use_partition = FALSE)  # single connected graph

  # Root cells: Mesothelial cells closest to cluster centroid in UMAP space
  meso_cells  <- which(cds$cell_type_ == "Mesothelial")
  umap_meso   <- reducedDims(cds)[["UMAP"]][meso_cells, ]
  centroid_m  <- colMeans(umap_meso)
  root_cell   <- meso_cells[which.min(rowSums((umap_meso - centroid_m)^2))]

  # In Monocle3 you typically call order_cells() interactively.
  # Here we pass the closest principal graph node to the root cell as root_pr_nodes.
  # The principal graph nodes are stored in principal_graph_aux.
  pr_nodes <- t(principal_graph_aux(cds)[["UMAP"]]$dp_mst)
  umap_root <- reducedDims(cds)[["UMAP"]][root_cell, , drop = FALSE]
  dist_to_pr <- sqrt(rowSums((pr_nodes - as.vector(umap_root))^2))
  root_pr    <- rownames(pr_nodes)[which.min(dist_to_pr)]
  cds        <- order_cells(cds, root_pr_nodes = root_pr)

  saveRDS(cds, file.path(out5, "cds.rds"))
  save_monocle_objects(
    cds,
    out5, #file.path(out5, "cds.rds")
  )

  # Graph-test: genes associated with the trajectory
  pr_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = N_CORES)
  pr_test_res <- pr_test_res %>%
    rownames_to_column("gene") %>%
    mutate(padj = p.adjust(q_value, "BH")) %>%    # q_value already BH; redundant but explicit
    arrange(q_value)

  write_csv(pr_test_res, file.path(out5, "graph_test_results.csv"))
  cat("  Top trajectory-associated genes (Monocle3 graph test):\n")
  print(head(pr_test_res, 20))

  # Condition × pseudotime interaction: fit GAM per gene
  # model: expr ~ splines::ns(pseudotime, df=3) * cohort
  # (Monocle3's fit_models handles this natively)
  cat("  Fitting condition × pseudotime interaction models...\n")

  # Subset to HVGs for speed
  cds_hvg <- cds[keep_genes, ]  # keep_genes defined in Analysis 4 above

  gene_fits <- fit_models(
    cds_hvg,
    model_formula_str = "~splines::ns(pseudotime, df=3) * cohort",
    cores             = N_CORES,
    verbose           = FALSE
  )

  fit_coefs <- coefficient_table(gene_fits)
  interaction_res <- fit_coefs %>%
    filter(grepl(":cohort", term)) %>%      # interaction terms
    select(gene_short_name, term, estimate, std_err, test_val, p_value, q_value) %>%
    arrange(q_value)

  write_csv(interaction_res, file.path(out5, "condition_pseudotime_interaction.csv"))
  cat("  Top condition × pseudotime interaction genes:\n")
  print(head(interaction_res, 20))

  # Visualise top trajectory-associated genes on UMAP
  top_traj_genes <- head(pr_test_res$gene, 9)
  pdf(file.path(out5, "monocle3_top_genes_umap.pdf"), width = 14, height = 14)
  print(plot_cells(cds,
    genes             = top_traj_genes,
    show_trajectory_graph = TRUE,
    label_cell_groups     = FALSE
  ))
  dev.off()

  pdf(file.path(out5, "monocle3_pseudotime_umap.pdf"), width = 10, height = 8)
  print(plot_cells(cds,
    color_cells_by    = "pseudotime",
    show_trajectory_graph = TRUE,
    label_cell_groups = FALSE,
    trajectory_graph_color = "grey40"
  ))
  dev.off()

  cat("  Monocle3 outputs in:", out5, "\n")

} else {
  cat("\n  [Skipping Analysis 5 — monocle3 not installed]\n")
  cat("  Install instructions in script header above.\n")
}

# ── Quick cross-analysis summary ─────────────────────────────────────────────────
cat("\n=== Summary: progressionTest p-values ===\n")
for (nm in c("om_adip", "sq_adip", "om_meso")) {
  rds <- file.path(OUTPUT_DIR,
                   c(om_adip = "om_adipocyte", sq_adip = "sq_adipocyte",
                     om_meso = "om_mesothelial")[nm],
                   paste0(nm, "_condiments.rds"))
  if (file.exists(rds)) {
    res <- readRDS(rds)
    cat(nm, "— global progressionTest p =",
        res$progression$p.value[res$progression$lineage == "global"], "\n")
  }
}

# Extended analysis summary
ext_cond_file <- file.path(out4, "downstream_conditionTest.csv")
if (file.exists(ext_cond_file)) {
  ext_res <- read_csv(ext_cond_file, show_col_types = FALSE)
  n_sig   <- sum(!is.na(ext_res$pvalue) & p.adjust(ext_res$pvalue, "BH") < 0.05)
  cat("om_meso_ext downstream conditionTest — FDR<0.05 genes:", n_sig, "\n")
}

cat("\nAll outputs in:", OUTPUT_DIR, "\n")

# ── Analysis 6: CytoTRACE — Differentiation State Scoring ────────────────────────
#
# CytoTRACE estimates differentiation potential from transcriptional diversity:
# cells expressing more genes are less differentiated (higher score = more stem-like).
# Score range: 0 (most differentiated) → 1 (least differentiated / most progenitor-like).
#
# WHY: answers the trajectory direction question without RNA velocity.
# If Mesothelial cells score highest → they are the least differentiated → root model
# supported. If FAPs or IGFBP2_cells score highest → reconsider the root assumption.
# Also tests whether T2D shifts differentiation state within each cell type.
#
# INSTALL:
#   devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r")
#
# NOTE: CytoTRACE2 is run on the 4-cell-type subset (sub_4ct) so scores are
# relative to these populations. For scores relative to all omentum cell types,
# re-run on sce_om with the full dataset (slower but more globally calibrated).

if (requireNamespace("CytoTRACE2", quietly = TRUE)) {
  library(CytoTRACE2)

  cat("\n=== Analysis 6: CytoTRACE2 (omentum 4-cell-type mesothelial) ===\n")

  out6 <- file.path(OUTPUT_DIR, "om_meso_cytotrace")
  dir.create(out6, recursive = TRUE, showWarnings = FALSE)

  # ── 6.1  Run CytoTRACE2 ───────────────────────────────────────────────────────
  # Input: genes × cells raw count matrix (genes as rows, cells as columns).
  # cytotrace2() expects a data.frame with genes as rows; pass raw counts.
  cat("  Running CytoTRACE2 on", ncol(sub_4ct), "cells ×", nrow(sub_4ct), "genes...\n")

  ct_result <- cytotrace2(as.data.frame(as.matrix(counts(sub_4ct))),
    ncores   = N_CORES, species = "human"
  )

  # Add scores to the SCE object
  sub_4ct$cytotrace        <- ct_result$CytoTRACE2_Score         # 0=differentiated, 1=progenitor
  sub_4ct$cytotrace_potency <- ct_result$CytoTRACE2_Potency # ordered potency category

  # Save per-cell scores
  cytotrace_df <- data.frame(
    cell_id        = colnames(sub_4ct),
    cytotrace      = sub_4ct$cytotrace,
    potency_class  = as.character(sub_4ct$cytotrace_potency),
    cell_type      = sub_4ct$cell_type_,
    cohort         = sub_4ct$cohort,
    participant    = sub_4ct$common_participant,
    pseudotime     = rowMeans(slingPseudotime(sds_4ct), na.rm = TRUE)
  )
  write_csv(cytotrace_df, file.path(out6, "cytotrace_scores.csv"))

  cat("  Mean CytoTRACE2 score per cell type (higher = more progenitor-like):\n")
  print(
    cytotrace_df %>%
      group_by(cell_type) %>%
      dplyr::summarise(
        mean_score   = mean(cytotrace),
        median_score = median(cytotrace),
        n            = dplyr::n(),
        .groups      = "drop"
      ) %>%
      arrange(dplyr::desc(mean_score))
  )

  # ── 6.2  Plots ────────────────────────────────────────────────────────────────
  CT_ORDER_4 <- c("Mesothelial", "IGFBP2_cells", "FAPs", "Pre_adipocytes")
  ct_cols4   <- c(Mesothelial    = "#4393c3",
                  IGFBP2_cells   = "#92c5de",
                  FAPs           = "#f4a582",
                  Pre_adipocytes = "#d6604d")
  coh_cols   <- c(Healthy = "#3182bd", Unhealthy = "#de2d26")

  cytotrace_df <- cytotrace_df %>%
    mutate(
      cell_type = factor(cell_type, levels = CT_ORDER_4),
      cohort    = factor(cohort,    levels = c("Healthy", "Unhealthy"))
    )

  # Panel A: CytoTRACE score by cell type — answers the root question
  p_ct_viol <- ggplot(cytotrace_df,
                      aes(cell_type, cytotrace, fill = cell_type)) +
    geom_violin(scale = "width", alpha = 0.8, colour = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA,
                 fill = "white", colour = "grey30", linewidth = 0.5) +
    scale_fill_manual(values = ct_cols4, guide = "none") +
    scale_x_discrete(name = NULL) +
    scale_y_continuous(name = "CytoTRACE score\n(1 = most progenitor-like)") +
    labs(title    = "CytoTRACE2: differentiation state by cell type",
         subtitle = "Cell type with highest score is the most progenitor-like → expected root") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  # Panel B: CytoTRACE by cell type × cohort — does T2D alter differentiation state?
  p_ct_cond <- ggplot(cytotrace_df,
                      aes(cohort, cytotrace, fill = cohort)) +
    geom_violin(scale = "width", alpha = 0.8, colour = NA) +
    geom_boxplot(width = 0.15, outlier.shape = NA,
                 fill = "white", colour = "grey30", linewidth = 0.5) +
    scale_fill_manual(values = coh_cols, guide = "none") +
    scale_y_continuous(name = "CytoTRACE score") +
    facet_wrap(~cell_type, nrow = 1, scales = "free_x") +
    labs(title    = "CytoTRACE2: T2D effect on differentiation state per cell type",
         subtitle = "Shift toward 1 in Unhealthy = more progenitor-like (less committed)") +
    theme_classic(base_size = 9) +
    theme(strip.text = element_text(face = "bold"))

  # Wilcoxon test per cell type: Healthy vs Unhealthy CytoTRACE
  wilcox_res <- cytotrace_df %>%
    group_by(cell_type) %>%
    dplyr::summarise(
      mean_H   = mean(cytotrace[cohort == "Healthy"]),
      mean_U   = mean(cytotrace[cohort == "Unhealthy"]),
      delta    = mean_U - mean_H,
      W_pval   = wilcox.test(
        cytotrace[cohort == "Healthy"],
        cytotrace[cohort == "Unhealthy"]
      )$p.value,
      .groups  = "drop"
    ) %>%
    mutate(padj = p.adjust(W_pval, "BH"))
  cat("  CytoTRACE2: Healthy vs Unhealthy per cell type (Wilcoxon):\n")
  print(wilcox_res)
  write_csv(wilcox_res, file.path(out6, "cytotrace_HvU_wilcoxon.csv"))

  # Panel C: CytoTRACE vs pseudotime — correlation validates trajectory direction
  # Expected: negative correlation (higher CytoTRACE = earlier pseudotime / closer to root)
  cor_res <- cytotrace_df %>%
    filter(!is.na(pseudotime)) %>%
    group_by(cell_type) %>%
    dplyr::summarise(
      r    = cor(cytotrace, pseudotime, method = "spearman"),
      n    = dplyr::n(),
      .groups = "drop"
    )
  cat("  Spearman correlation CytoTRACE2 ~ pseudotime per cell type:\n")
  print(cor_res)

  p_ct_pt <- ggplot(
    cytotrace_df %>% filter(!is.na(pseudotime)) %>%
      left_join(cor_res, by = "cell_type"),
    aes(pseudotime, cytotrace, colour = cell_type)
  ) +
    geom_point(size = 0.3, alpha = 0.15) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    geom_text(
      data = cor_res,
      aes(x = Inf, y = Inf,
          label = sprintf("ρ = %.2f", r),
          colour = cell_type),
      hjust = 1.1, vjust = 1.5, size = 3, fontface = "bold", inherit.aes = FALSE
    ) +
    scale_colour_manual(values = ct_cols4, guide = "none") +
    facet_wrap(~cell_type, nrow = 1) +
    labs(x = "DPT pseudotime", y = "CytoTRACE score",
         title    = "CytoTRACE2 vs pseudotime",
         subtitle = "Negative ρ = higher CytoTRACE2 at early pseudotime → trajectory direction confirmed") +
    theme_classic(base_size = 9) +
    theme(strip.text = element_text(face = "bold"))

  # Panel D: CytoTRACE on UMAP, coloured continuously
  umap_df <- reducedDim(sub_4ct, "UMAP") %>%
    as.data.frame() %>%
    setNames(c("UMAP1", "UMAP2")) %>%
    bind_cols(cytotrace_df %>% select(cytotrace, cell_type, cohort))

  p_ct_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = cytotrace)) +
    geom_point(size = 0.4, alpha = 0.6) +
    scale_colour_viridis_c(option = "magma", direction = -1,
                           name = "CytoTRACE2\n(1=progenitor)") +
    facet_wrap(~cohort) +
    labs(title = "CytoTRACE2 on UMAP by condition") +
    theme_classic(base_size = 9)

  # Assemble and save
  p_cytotrace <- (p_ct_viol | p_ct_umap) /
                  p_ct_cond /
                  p_ct_pt +
    plot_annotation(title    = "CytoTRACE2 — Differentiation State Analysis",
                    subtitle = "4-cell-type mesothelial trajectory, omentum")

  ggsave(file.path(out6, "cytotrace_results.pdf"), p_cytotrace,
         width = 14, height = 14)

  # ── 6.3  Potency category distribution ────────────────────────────────────────
  # CytoTRACE2 assigns each cell an ordered potency class
  # (Differentiated → Unipotent → Oligopotent → Multipotent → Pluripotent → Totipotent).
  potency_tbl <- cytotrace_df %>%
    dplyr::count(cell_type, potency_class) %>%
    group_by(cell_type) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()
  cat("  CytoTRACE2 potency class distribution per cell type:\n")
  print(potency_tbl)
  write_csv(potency_tbl, file.path(out6, "cytotrace2_potency_distribution.csv"))

  p_potency <- ggplot(potency_tbl,
                      aes(x = factor(cell_type, levels = CT_ORDER_4),
                          y = pct, fill = potency_class)) +
    geom_col(position = "stack", colour = "white", linewidth = 0.2) +
    scale_fill_viridis_d(option = "plasma", direction = -1,
                         name = "Potency class") +
    scale_y_continuous(name = "% of cells") +
    scale_x_discrete(name = NULL) +
    labs(title    = "CytoTRACE2 potency class by cell type",
         subtitle = "Higher potency (purple) = more progenitor-like") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))

  ggsave(file.path(out6, "cytotrace2_potency_barplot.pdf"), p_potency,
         width = 7, height = 5)

  cat("  CytoTRACE2 outputs saved to:", out6, "\n")

} else {
  cat("\n  [Skipping Analysis 6 — CytoTRACE2 not installed]\n")
  cat("  Install: devtools::install_github(\"digitalcytometry/cytotrace2\", subdir = \"cytotrace2_r\")\n")
}

# ── Analysis 7: Slingshot infers the MST freely ────────────────────────
#
# Let Slingshot infer the MST freely
# The most direct test. Run Slingshot without specifying start.clus or end.clus — 
# it fits a minimum spanning tree on the cluster centroids 
# in reduced space and picks the root geometrically (the most "extreme" cluster). 
# Run it three ways and compare:

# A: free topology — no constraints at all
sds_free <- slingshot(sub_4ct, clusterLabels='cell_type_', reducedDim='HARMONY')

# B: IGFBP2_cells as root
sds_igfbp2 <- slingshot(sub_4ct, clusterLabels='cell_type_', reducedDim='HARMONY',
                        start.clus='IGFBP2_cells')

# C: FAPs as root with two endpoints
sds_fap <- slingshot(sub_4ct, clusterLabels='cell_type_', reducedDim='HARMONY',
                     start.clus='FAPs')

# Inspect the inferred lineages in each case
slingLineages(sds_free)    # what does the data prefer?
slingLineages(sds_igfbp2)
slingLineages(sds_fap)

saveRDS(sds_igfbp2, file.path(out4, "slingshot_igfbp2.rds"))
sds_igfbp2 <- readRDS(file.path(out4, "slingshot_igfbp2.rds"))

# If sds_free produces a Y-shape with FAPs or IGFBP2_cells as the internal node, 
# that's evidence for the branching model. 
# Check slingBranchID() to see which cells are assigned to which arm.


plot_trajectory(sds_igfbp2, "cell_type_", "Omentum: IGFBP2 start",
                file.path(out4, "om_igfbp2start"))

cond_igfbp2 <- factor(sds_igfbp2$cohort, levels = c("Healthy", "Unhealthy"))

# condiments: does progression differ between conditions along the full trajectory?
cm_igfbp2 <- run_condiments(sds_igfbp2, cond_igfbp2, nrep = 500,
                         out_prefix = file.path(out4, "om_igfbp2start"))

# ── condiments visualisation for sds_igfbp2 ───────────────────────────────────
#
# Four panels:
#   A. Imbalance score on UMAP     — spatial separation of H vs U cells
#   B. Progression test            — do conditions advance differently along each lineage?
#   C. Topology test               — does the branching structure itself differ?
#   D. Differential differentiation — do conditions choose different fates?

out_cm <- file.path(out4, "condiments_igfbp2")
dir.create(out_cm, recursive = TRUE, showWarnings = FALSE)

ct_cols_cm <- c(Mesothelial    = "#4393c3",
                IGFBP2_cells   = "#92c5de",
                FAPs           = "#f4a582",
                Pre_adipocytes = "#d6604d")
coh_cols_cm <- c(Healthy = "#3182bd", Unhealthy = "#de2d26")

# helper: safely pull global progression p
.prog_global_p <- function(prog) {
  idx <- which(prog$lineage %in% c("All", "global") |
               suppressWarnings(is.na(as.integer(prog$lineage))))
  if (length(idx)) prog$p.value[idx[1]] else prog$p.value[1]
}

# ── A. Imbalance score on UMAP ────────────────────────────────────────────────
# Positive → overrepresented in Unhealthy; negative → overrepresented in Healthy.
# imbalance_score() returns $scores and $scaled_scores; use scaled for plotting.
imb_scores <- if (!is.null(cm_igfbp2$imbalance$scaled_scores))
                cm_igfbp2$imbalance$scaled_scores
              else
                cm_igfbp2$imbalance$scores

umap_cm <- data.frame(
  UMAP1     = reducedDim(sds_igfbp2, "UMAP")[, 1],
  UMAP2     = reducedDim(sds_igfbp2, "UMAP")[, 2],
  imbalance = imb_scores$scaled_scores,
  cell_type = sds_igfbp2$cell_type_,
  cohort    = sds_igfbp2$cohort
)

p_imb_umap <- ggplot(umap_cm, aes(UMAP1, UMAP2, colour = imbalance)) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_colour_gradient2(low  = "#3182bd", mid = "grey90", high = "#de2d26",
                         midpoint = 0,
                         name = "Imbalance\n(+= Unhealthy\n−= Healthy)") +
  labs(title    = "Imbalance score",
       subtitle = "Spatial separation of conditions on UMAP") +
  theme_classic(base_size = 10)

p_imb_ct <- ggplot(umap_cm, aes(UMAP1, UMAP2, colour = imbalance)) +
  geom_point(size = 0.4, alpha = 0.6) +
  scale_colour_gradient2(low = "#3182bd", mid = "grey90", high = "#de2d26",
                         midpoint = 0, name = "Imbalance") +
  facet_wrap(~cell_type, nrow = 1) +
  labs(title = "Imbalance score by cell type") +
  theme_classic(base_size = 9) +
  theme(strip.text = element_text(face = "bold"))

# ── B. Progression test ───────────────────────────────────────────────────────
pt_mat  <- slingPseudotime(sds_igfbp2)
n_lin   <- ncol(pt_mat)
lin_nms <- paste0("Lineage ", seq_len(n_lin))

pt_long <- do.call(rbind, lapply(seq_len(n_lin), function(l) {
  data.frame(
    pseudotime = pt_mat[, l],
    cohort     = sds_igfbp2$cohort,
    cell_type  = sds_igfbp2$cell_type_,
    lineage    = lin_nms[l]
  )
})) %>% filter(!is.na(pseudotime))

global_p_prog <- .prog_global_p(cm_igfbp2$progression)

p_prog_dens <- ggplot(pt_long, aes(pseudotime, colour = cohort, fill = cohort)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  scale_colour_manual(values = coh_cols_cm) +
  scale_fill_manual(values = coh_cols_cm) +
  facet_wrap(~lineage, scales = "free_x", nrow = 1) +
  labs(x        = "Pseudotime",
       y        = "Density",
       title    = "progressionTest: pseudotime distribution by condition",
       subtitle = sprintf("Global p = %.3g", global_p_prog)) +
  theme_classic(base_size = 9) +
  theme(legend.position = "top")

# -log10 p-value bar per lineage (including global)
prog_bar_df <- cm_igfbp2$progression %>%
  mutate(
    label        = ifelse(lineage %in% c("All", "global"), "Global", paste0("Lineage ", lineage)),
    neg_log10_p  = -log10(pmax(p.value, 1e-4)),
    sig          = ifelse(p.value < 0.05, "p < 0.05", "n.s.")
  )

p_prog_pval <- ggplot(prog_bar_df, aes(x = label, y = neg_log10_p, fill = sig)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("p < 0.05" = "#de2d26", "n.s." = "grey70"),
                    name = NULL) +
  labs(x = NULL, y = "-log10(p-value)",
       title = "progressionTest p-values") +
  theme_classic(base_size = 10) +
  theme(legend.position = "top")

# ── C. Topology test ──────────────────────────────────────────────────────────
# Tests whether the branching structure (MST topology) differs between conditions.
# Set RUN_TOPOLOGY_TEST <- TRUE to enable (slow: ~200 permutations of full UMAP).
RUN_TOPOLOGY_TEST <- FALSE

if (RUN_TOPOLOGY_TEST) {
  cat("  condiments: topologyTest (IGFBP2-start)...\n")
  topo_res <- topologyTest(sds_igfbp2, conditions = cond_igfbp2, rep = 200)
  cat("  topologyTest:\n")
  print(topo_res)
  write_csv(as.data.frame(topo_res), file.path(out_cm, "topology_test.csv"))
  topo_p <- if ("p.value" %in% names(topo_res)) topo_res$p.value[1] else NA_real_
} else {
  cat("  [topologyTest skipped — set RUN_TOPOLOGY_TEST <- TRUE to enable]\n")
  topo_res <- NULL
  topo_p   <- NA_real_
}

p_topo <- ggplot(
  data.frame(test = "Topology", p = topo_p,
             neg_log10_p = -log10(pmax(topo_p, 1e-4)),
             sig = ifelse(!is.na(topo_p) & topo_p < 0.05, "p < 0.05", "n.s.")),
  aes(x = test, y = neg_log10_p, fill = sig)
) +
  geom_col(width = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("p < 0.05" = "#de2d26", "n.s." = "grey70"),
                    name = NULL) +
  annotate("text", x = 1, y = -log10(pmax(topo_p, 1e-4)) + 0.1,
           label = if (is.na(topo_p)) "skipped" else sprintf("p = %.3g", topo_p),
           vjust = 0, size = 3.5) +
  labs(x = NULL, y = "-log10(p-value)",
       title    = "topologyTest",
       subtitle = if (RUN_TOPOLOGY_TEST) "Does branching structure differ by condition?"
                  else "Set RUN_TOPOLOGY_TEST <- TRUE to run") +
  theme_classic(base_size = 10) +
  theme(legend.position = "top")

# ── D. Differential differentiation: fate probability ─────────────────────────
# slingCurveWeights (normalised to sum = 1 per cell) = fate probability per lineage.
if (!is.null(cm_igfbp2$differentiation) && n_lin > 1) {

  wts_raw <- slingCurveWeights(sds_igfbp2)
  wts     <- sweep(wts_raw, 1, rowSums(wts_raw), "/")   # row-normalise

  wts_long <- do.call(rbind, lapply(seq_len(n_lin), function(l) {
    data.frame(
      fate_prob = wts[, l],
      lineage   = lin_nms[l],
      cohort    = sds_igfbp2$cohort,
      cell_type = sds_igfbp2$cell_type_
    )
  }))

  diff_p <- if ("p.value" %in% names(cm_igfbp2$differentiation))
               cm_igfbp2$differentiation$p.value[1] else NA_real_

  p_fate_viol <- ggplot(wts_long, aes(cohort, fate_prob, fill = cohort)) +
    geom_violin(scale = "width", alpha = 0.8, colour = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA,
                 fill = "white", colour = "grey30", linewidth = 0.4) +
    scale_fill_manual(values = coh_cols_cm, guide = "none") +
    facet_wrap(~lineage, nrow = 1) +
    labs(y        = "Fate probability (curve weight)",
         x        = NULL,
         title    = "differentiationTest: fate probabilities by condition",
         subtitle = sprintf("differentiationTest p = %.3g", diff_p)) +
    theme_classic(base_size = 9) +
    theme(strip.text = element_text(face = "bold"))

  fate_umap_plots <- lapply(seq_len(n_lin), function(l) {
    df <- umap_cm
    df$fate_prob <- wts[, l]
    ggplot(df, aes(UMAP1, UMAP2, colour = fate_prob)) +
      geom_point(size = 0.4, alpha = 0.6) +
      scale_colour_viridis_c(option = "inferno", name = "Fate\nprob.") +
      facet_wrap(~cohort) +
      labs(title = paste0(lin_nms[l], " fate probability")) +
      theme_classic(base_size = 8) +
      theme(strip.text = element_text(face = "bold"))
  })

  # Combine into a single wrap_plots call to avoid patchwork flattening
  # a nested patchwork/p_fate_viol structure before plot_layout is applied.
  p_diffn <- wrap_plots(
    c(fate_umap_plots, list(p_fate_viol)),
    ncol    = 1,
    heights = c(rep(3, n_lin), 2)
  )

  ggsave(file.path(out_cm, "differentiation_test.pdf"),
         p_diffn, width = 12, height = n_lin * 4 + 4)

} else {
  p_fate_viol <- NULL
}

# ── Assemble condiments summary figure ────────────────────────────────────────
diff_p_label <- if (!is.null(cm_igfbp2$differentiation) && "p.value" %in% names(cm_igfbp2$differentiation))
                  sprintf(" | differentiationTest p = %.3g", cm_igfbp2$differentiation$p.value[1])
                else ""

# Wrap layout in explicit parens before + plot_annotation:
# without parens, R's operator precedence makes + bind to p_prog_dens only.
p_condiments_main <- ((p_imb_umap | p_prog_pval ) /  # p_topo
                       p_imb_ct /
                       p_prog_dens) +
  plot_annotation(
    title    = "Condiments — differential trajectory analysis (IGFBP2-start)",
    subtitle = paste0(
      "progressionTest global p = ", formatC(global_p_prog, format = "g", digits = 3),
      #" | topologyTest p = ",        formatC(topo_p,       format = "g", digits = 3),
      diff_p_label
    )
  )

ggsave(file.path(out_cm, "condiments_summary.pdf"),
       p_condiments_main, width = 16, height = 14)

cat("  condiments visualisations saved to:", out_cm, "\n")

# running tradeSeq on all 4 cell-types

sex_U_igfbp2 <- model.matrix(~ sds_igfbp2$sex)[, -1, drop = FALSE]

ts_igfbp2 <- run_tradeseq(sds_igfbp2, cond_igfbp2, U = sex_U_igfbp2,
                          n_knots = N_KNOTS,
                          out_prefix = file.path(out4, "om_igfbp2start"))

plot_tradeseq_results(ts_igfbp2$sce_gam, ts_igfbp2$condition_res,
                      sds_igfbp2, ts_igfbp2$cnt,
                      cell_type_col = "cell_type_",
                      title = "Omentum: mesothelial-FAPs transition",
                      out_prefix = file.path(out4, "om_igfbp2start"),
                      sve_res = ts_igfbp2$sve_res)
