#!/usr/bin/env Rscript
# MOFA analysis: omentum vs. subcutaneous adipose pseudo-bulk
# Two runs: fine cell types (cell_type2_) and broad cell types (cell_type_)
#
# Key design: one MOFA model per cell type (2 views: Omentum, Subcutaneous).
# Keeps cell-type-composition variation out of the factors so that each model
# asks specifically: "which gene programs are co-regulated between depots for
# this cell type?"  A second "all-cell-types combined" model is also fit for
# a global overview.
#
# Input:  pseudo-bulk h5ad files (participant × cell_type aggregated)
# Output: MOFA model (.rds + .hdf5) and diagnostic PDFs

suppressPackageStartupMessages({
  library(zellkonverter)         # readH5AD → SingleCellExperiment
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(DESeq2)
  library(MOFA2)
  library(tidyverse)
})

# ── Parameters ─────────────────────────────────────────────────────────────────
INPUT_DIR  <- "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq/inputs/anndatas"
OUTPUT_DIR <- "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq/results/MOFA"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

MIN_COUNT       <- 10    # raw counts: gene must reach this in >= MIN_FRAC of samples
MIN_FRAC        <- 0.5   # fraction of within-cell-type samples that must meet MIN_COUNT
N_TOP_FEATURES  <- 500   # top-variance genes per cell type (n=10 → keep signal:noise high)
N_FACTORS       <- 5     # ≤ n_samples/2; 8 overfits with n=10
CONVERGENCE     <- "slow"
SEED            <- 42

# ── Input files ────────────────────────────────────────────────────────────────
FILES <- list(
  fine = list(
    omentum = file.path(INPUT_DIR, "omentum_reannotated_finely_pb_cell_type2_.h5ad"),
    sq      = file.path(INPUT_DIR, "sq_reannotated_finely_pb_cell_type2_.h5ad")
  ),
  broad = list(
    omentum = file.path(INPUT_DIR, "omentum_reannotated_finely_pb_cell_type_.h5ad"),
    sq      = file.path(INPUT_DIR, "sq_reannotated_finely_pb_cell_type_.h5ad")
  )
)

# cell type annotation column per resolution
CT_COL <- list(fine = "cell_type2_", broad = "cell_type_")

# ── I/O helper ─────────────────────────────────────────────────────────────────
read_pseudobulk <- function(path) {
  sce <- readH5AD(path, use_hdf5 = FALSE)

  # Build Entrez-ID → gene-symbol lookup from var metadata
  var_df     <- as.data.frame(rowData(sce))
  id_to_name <- if ("entrez_name" %in% colnames(var_df)) {
    setNames(as.character(var_df$entrez_name), rownames(sce))
  } else {
    setNames(rownames(sce), rownames(sce))   # fall back: ID = name
  }

  list(
    counts     = assay(sce, "X"),             # genes × obs (participant_celltype)
    obs        = as.data.frame(colData(sce)), # obs metadata
    id_to_name = id_to_name                  # named vector: entrez_id → symbol
  )
}

# ── Per-cell-type normalisation ─────────────────────────────────────────────────
filter_and_rlog <- function(count_mat) {
  # count_mat: genes × n_participants integer matrix for one cell type
  # Returns rlog matrix (same dims) or NULL when too few genes pass

  keep <- rowSums(count_mat >= MIN_COUNT) >= ceiling(ncol(count_mat) * MIN_FRAC)
  if (sum(keep) < 10) {
    message("    Only ", sum(keep), " genes pass filter — skipping")
    return(NULL)
  }
  count_mat <- count_mat[keep, , drop = FALSE]
  message("    ", sum(keep), " genes retained")

  storage.mode(count_mat) <- "integer"
  coldata <- data.frame(row.names = colnames(count_mat),
                        sample    = colnames(count_mat))

  dds <- DESeqDataSetFromMatrix(count_mat, colData = coldata, design = ~1)
  rld <- rlog(dds, blind = TRUE)
  assay(rld)
}

select_top_variable <- function(mat, n) {
  rv  <- apply(mat, 1, var)
  top <- order(rv, decreasing = TRUE)[seq_len(min(n, nrow(mat)))]
  mat[top, , drop = FALSE]
}

# ── Build per-cell-type rlog matrices for one tissue ──────────────────────────
# Returns list(ct_matrices = named list of (genes × participants) matrices,
#              sample_meta = data.frame,
#              all_parts   = character vector of all participant IDs)
build_tissue_data <- function(path, ct_col,
                              participant_col = "common_participant") {
  message("  Reading: ", basename(path))
  pb         <- read_pseudobulk(path)
  obs        <- pb$obs
  cnts       <- pb$counts
  id_to_name <- pb$id_to_name

  participants <- obs[[participant_col]]
  cell_types   <- obs[[ct_col]]
  all_parts    <- sort(unique(participants))
  unique_cts   <- sort(unique(cell_types))

  message("  ", length(all_parts), " participants | ",
          length(unique_cts), " cell types | ", nrow(cnts), " genes")

  sample_meta <- obs %>%
    dplyr::select(all_of(c(participant_col, "cohort", "sex"))) %>%
    dplyr::distinct() %>%
    dplyr::arrange(.data[[participant_col]]) %>%
    dplyr::rename(sample = all_of(participant_col)) %>%
    dplyr::mutate(
      cohort_num = as.integer(cohort == "Unhealthy"),
      sex_num    = as.integer(sex == "M")
    )

  ct_matrices <- list()

  for (ct in unique_cts) {
    message("  [", ct, "]")
    idx      <- which(cell_types == ct)
    ct_parts <- participants[idx]
    ct_cnts  <- cnts[, idx, drop = FALSE]
    colnames(ct_cnts) <- ct_parts

    if (ncol(ct_cnts) < 3) { message("    < 3 participants — skipping"); next }

    rlog_mat <- tryCatch(
      filter_and_rlog(ct_cnts),
      error = function(e) { message("    Error: ", e$message); NULL }
    )
    if (is.null(rlog_mat)) next

    if (nrow(rlog_mat) > N_TOP_FEATURES)
      rlog_mat <- select_top_variable(rlog_mat, N_TOP_FEATURES)

    # Entrez ID → gene symbol (fall back to ID on NA or duplicate)
    ids        <- rownames(rlog_mat)
    symbols    <- id_to_name[ids]
    use_sym    <- !is.na(symbols) & nzchar(symbols)
    gene_label <- make.unique(ifelse(use_sym, symbols, ids), sep = "_dup")
    rownames(rlog_mat) <- gene_label          # clean symbol names; no __CT suffix here

    ct_matrices[[ct]] <- rlog_mat
    message("    Retained ", nrow(rlog_mat), " features for MOFA")
  }

  if (length(ct_matrices) == 0)
    stop("No cell types passed QC for: ", basename(path))

  list(ct_matrices = ct_matrices, sample_meta = sample_meta, all_parts = all_parts)
}

# ── Stack all cell-type matrices into one view (for the combined model) ────────
# Appends __CT to rownames so genes from different cell types don't collide.
build_combined_matrix <- function(tissue_data) {
  all_parts <- tissue_data$all_parts
  blocks    <- lapply(names(tissue_data$ct_matrices), function(ct) {
    m        <- tissue_data$ct_matrices[[ct]]
    ct_label <- gsub("[^A-Za-z0-9]", "_", ct)
    rownames(m) <- paste0(rownames(m), "__", ct_label)
    m
  })
  all_features <- unlist(lapply(blocks, rownames), use.names = FALSE)
  full_mat <- matrix(NA_real_, nrow = length(all_features), ncol = length(all_parts),
                     dimnames = list(all_features, all_parts))
  for (m in blocks) {
    samps <- intersect(colnames(m), all_parts)
    full_mat[rownames(m), samps] <- m[, samps]
  }
  pct_obs <- 100 * mean(!is.na(full_mat))
  message("  Combined matrix: ", nrow(full_mat), " features × ", ncol(full_mat),
          " participants  (", round(pct_obs, 1), "% observed)")
  full_mat
}

# ── Run MOFA and save ──────────────────────────────────────────────────────────
# mofa_data   : named list of (features × participants) matrices, one per view
# sample_meta : data.frame with columns sample / cohort / sex / cohort_num / sex_num
run_mofa <- function(mofa_data, sample_meta, output_prefix) {

  common_parts <- Reduce(intersect, lapply(mofa_data, colnames))
  stopifnot(length(common_parts) >= 4)
  message("  Common participants (n=", length(common_parts), "): ",
          paste(common_parts, collapse = ", "))

  mofa_data <- lapply(mofa_data, function(m) m[, common_parts, drop = FALSE])
  for (nm in names(mofa_data))
    message("  ", nm, ": ", nrow(mofa_data[[nm]]), " features")

  mofa <- create_mofa(mofa_data)

  data_opts              <- get_default_data_options(mofa)
  data_opts$scale_views  <- TRUE    # normalise variance across views

  model_opts             <- get_default_model_options(mofa)
  model_opts$num_factors <- N_FACTORS

  train_opts                  <- get_default_training_options(mofa)
  train_opts$convergence_mode <- CONVERGENCE
  train_opts$seed             <- SEED
  train_opts$maxiter          <- 2000

  mofa <- prepare_mofa(mofa,
    data_options     = data_opts,
    model_options    = model_opts,
    training_options = train_opts
  )

  hdf5_path <- paste0(output_prefix, "_mofa.hdf5")
  mofa <- MOFA2::run_mofa(mofa, outfile = hdf5_path, use_basilisk = TRUE)

  meta <- sample_meta %>% dplyr::filter(sample %in% common_parts)
  samples_metadata(mofa) <- as.data.frame(meta)

  saveRDS(mofa, paste0(output_prefix, "_mofa.rds"))
  message("  Model saved → ", basename(paste0(output_prefix, "_mofa.rds")))

  mofa
}

# ── Diagnostic plots ───────────────────────────────────────────────────────────
plot_mofa_results <- function(mofa, output_prefix) {

  n_factors      <- get_dimensions(mofa)[["K"]]
  cohort_colours <- c("Healthy" = "#4393C3", "Unhealthy" = "#D6604D")
  view_names     <- views_names(mofa)

  # ── 1. Variance explained ──
  pdf(paste0(output_prefix, "_variance_explained.pdf"), width = 10, height = 5)
  var_plots <- plot_variance_explained(mofa, plot_total = TRUE)
  if (is.list(var_plots)) lapply(var_plots, print) else print(var_plots)
  dev.off()

  # ── 2. Factor–covariate associations ──
  # Pearson r + Wilcoxon p-value (n=5 vs 5 makes Pearson unreliable alone).
  # Wilcoxon minimum possible two-tailed p at n=5 vs 5 is 0.008, so
  # interpret effect size (r, Cohen's d) rather than p-values.
  factor_scores <- get_factors(mofa, factors = "all")[[1]]   # samples × factors
  meta          <- samples_metadata(mofa)
  meta          <- meta[match(rownames(factor_scores), meta$sample), ]

  cohort_assoc <- do.call(rbind, lapply(seq_len(n_factors), function(f) {
    scores  <- factor_scores[, f]
    healthy <- scores[meta$cohort == "Healthy"]
    unheal  <- scores[meta$cohort == "Unhealthy"]
    wt      <- suppressWarnings(wilcox.test(healthy, unheal, exact = FALSE))
    d       <- (mean(unheal) - mean(healthy)) /
               sqrt((var(unheal) * (length(unheal) - 1) +
                     var(healthy) * (length(healthy) - 1)) /
                    (length(unheal) + length(healthy) - 2))
    r_cor   <- suppressWarnings(cor(scores, meta$cohort_num, method = "pearson"))
    data.frame(factor = paste0("Factor", f), wilcox_p = wt$p.value,
               cohens_d = d, pearson_r = r_cor)
  }))
  write.csv(cohort_assoc,
            paste0(output_prefix, "_cohort_association.csv"), row.names = FALSE)

  pdf(paste0(output_prefix, "_factor_covariate_cor.pdf"), width = 8, height = 5)
  tryCatch({
    p <- correlate_factors_with_covariates(
      mofa, covariates = c("cohort_num", "sex_num"), plot = "log_pval")
    print(p)
  }, error = function(e) message("Covariate correlation failed: ", e$message))
  dev.off()

  # ── 3. Individual factor violin plots coloured by cohort ──
  pdf(paste0(output_prefix, "_factor_values.pdf"), width = 5, height = 4)
  for (f in seq_len(n_factors)) {
    pval_lab <- sprintf("Wilcoxon p=%.3f  d=%.2f",
                        cohort_assoc$wilcox_p[f], cohort_assoc$cohens_d[f])
    p <- plot_factor(mofa, factor = f, color_by = "cohort", shape_by = "sex",
                     dot_size = 3, dodge = TRUE,
                     add_violin = TRUE, violin_alpha = 0.25) +
      scale_color_manual(values = cohort_colours) +
      ggtitle(paste0("Factor ", f, "  |  ", pval_lab))
    print(p)
  }
  dev.off()

  # ── 4. Factor scatter coloured by cohort ──
  pdf(paste0(output_prefix, "_factor_scatter.pdf"), width = 6, height = 5)
  for (f in seq_len(max(1, n_factors - 1))) {
    p <- plot_factors(mofa, factors = c(f, f + 1),
                      color_by = "cohort", shape_by = "sex", dot_size = 4) +
      scale_color_manual(values = cohort_colours) +
      ggtitle(paste("Factor", f, "vs Factor", f + 1))
    print(p)
  }
  dev.off()

  # ── 5. Top feature weights per factor and view ──
  pdf(paste0(output_prefix, "_top_weights.pdf"), width = 9, height = 5)
  for (f in seq_len(n_factors)) {
    for (view in view_names) {
      tryCatch({
        p <- plot_top_weights(mofa, view = view, factor = f,
                              nfeatures = 20, scale = TRUE, abs = FALSE) +
          ggtitle(paste("Factor", f, "–", view, "top weights"))
        print(p)
      }, error = function(e) NULL)
    }
  }
  dev.off()

  message("  Plots written.")
}

# ── Cross-depot Spearman correlation (per cell type, per gene) ─────────────────
# For each shared cell type, correlates each gene's pseudobulk expression in
# omentum vs subcutaneous across participants.
# Genes with high |rho| are candidates for systemic regulation.
run_cross_depot_correlation <- function(om_data, sq_data, output_prefix) {

  shared_cts   <- intersect(names(om_data$ct_matrices), names(sq_data$ct_matrices))
  shared_parts <- intersect(om_data$all_parts, sq_data$all_parts)
  message("  Cross-depot correlation: ", length(shared_cts), " shared cell types, n=",
          length(shared_parts))

  all_results <- list()

  for (ct in shared_cts) {
    om_mat <- om_data$ct_matrices[[ct]][, shared_parts, drop = FALSE]
    sq_mat <- sq_data$ct_matrices[[ct]]
    shared_genes <- intersect(rownames(om_mat), rownames(sq_mat))
    if (length(shared_genes) < 10) next
    om_mat <- om_mat[shared_genes, ]
    sq_mat <- sq_mat[shared_genes, shared_parts, drop = FALSE]

    # Spearman rho per gene across participants
    rho_vec <- sapply(shared_genes, function(g) {
      suppressWarnings(cor(om_mat[g, ], sq_mat[g, ], method = "spearman"))
    })

    # Permutation-based p-value (1000 shuffles); important given n=10
    set.seed(SEED)
    n_perm  <- 1000
    perm_rho <- replicate(n_perm, {
      perm_idx <- sample(length(shared_parts))
      sapply(shared_genes, function(g)
        suppressWarnings(cor(om_mat[g, ], sq_mat[g, perm_idx], method = "spearman"))
      )
    })
    # perm_rho is genes × permutations
    perm_p <- rowMeans(abs(perm_rho) >= abs(rho_vec))

    res <- data.frame(
      gene      = shared_genes,
      cell_type = ct,
      spearman_r = rho_vec,
      perm_pval  = perm_p,
      row.names  = NULL
    )
    res <- res[order(-abs(res$spearman_r)), ]
    all_results[[ct]] <- res
    message("    ", ct, ": ", sum(perm_p < 0.05, na.rm = TRUE),
            " genes perm-p < 0.05 (top: ", head(shared_genes[order(-abs(rho_vec))], 3), ")")
  }

  result_df <- do.call(rbind, all_results)
  write.csv(result_df, paste0(output_prefix, "_crossdepot_spearman.csv"), row.names = FALSE)

  # Plot: rho distribution per cell type + top-correlated genes
  pdf(paste0(output_prefix, "_crossdepot_correlation.pdf"), width = 10, height = 5)

  p_dist <- ggplot(result_df, aes(x = spearman_r, fill = cell_type)) +
    geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
    facet_wrap(~cell_type, scales = "free_y") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_classic(base_size = 10) +
    xlab("Cross-depot Spearman ρ") + ylab("Genes") +
    theme(legend.position = "none") +
    ggtitle("Cross-depot gene correlation distribution per cell type")
  print(p_dist)

  # Dot plot: top 15 positively and negatively correlated genes per cell type
  top_genes <- result_df %>%
    dplyr::group_by(cell_type) %>%
    dplyr::slice_max(order_by = abs(spearman_r), n = 15) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(gene = reorder(gene, spearman_r))

  p_top <- ggplot(top_genes, aes(x = spearman_r, y = gene, color = spearman_r)) +
    geom_point(size = 2.5) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_gradient2(low = "#D6604D", mid = "grey80", high = "#4393C3",
                          midpoint = 0, name = "ρ") +
    facet_wrap(~cell_type, scales = "free_y") +
    theme_classic(base_size = 9) +
    xlab("Cross-depot Spearman ρ") + ylab(NULL) +
    ggtitle("Top cross-depot correlated genes per cell type")
  print(p_top)

  dev.off()
  message("  Saved: ", basename(paste0(output_prefix, "_crossdepot_spearman.csv")))
  message("  Saved: ", basename(paste0(output_prefix, "_crossdepot_correlation.pdf")))

  invisible(result_df)
}

# ── Main loop ──────────────────────────────────────────────────────────────────
for (resolution in c("broad", "fine")) {

  message("\n", strrep("═", 65))
  message("  Resolution: ", toupper(resolution))
  message(strrep("═", 65))

  ct_col <- CT_COL[[resolution]]

  message("\n[Loading Omentum]")
  om <- build_tissue_data(FILES[[resolution]]$omentum, ct_col)

  message("\n[Loading Subcutaneous]")
  sq <- build_tissue_data(FILES[[resolution]]$sq, ct_col)

  # ── A. Combined model: all cell types in one big view per depot ─────────────
  # Good for a global overview of inter-depot co-variation.
  # Caveat: factors are dominated by cell-type composition differences across
  # participants, so cohort signal will be hard to see here.
  message("\n[A] Combined model — all cell types")
  mofa_data_combined <- list(
    Omentum      = build_combined_matrix(om),
    Subcutaneous = build_combined_matrix(sq)
  )
  prefix_combined <- file.path(OUTPUT_DIR, paste0(resolution, "_combined"))
  mofa_combined   <- run_mofa(mofa_data_combined, om$sample_meta, prefix_combined)
  plot_mofa_results(mofa_combined, prefix_combined)

  # ── B. Per-cell-type models: one MOFA per shared cell type ─────────────────
  # Each model has exactly 2 views (Omentum, Subcutaneous) for one cell type.
  # Cell-type composition variation is removed entirely; factors capture only
  # within-cell-type cross-depot co-regulation.
  shared_cts <- intersect(names(om$ct_matrices), names(sq$ct_matrices))
  message("\n[B] Per-cell-type models (", length(shared_cts), " shared cell types)")

  for (ct in shared_cts) {
    ct_safe  <- gsub("[^A-Za-z0-9]", "_", ct)
    prefix_ct <- file.path(OUTPUT_DIR, paste0(resolution, "_", ct_safe))
    message("\n  --- ", ct, " ---")

    mofa_data_ct <- list(
      Omentum      = om$ct_matrices[[ct]],
      Subcutaneous = sq$ct_matrices[[ct]]
    )
    tryCatch({
      mofa_ct <- run_mofa(mofa_data_ct, om$sample_meta, prefix_ct)
      plot_mofa_results(mofa_ct, prefix_ct)
    }, error = function(e) message("  Skipping ", ct, ": ", e$message))
  }

  # ── C. Cross-depot Spearman correlation (no MOFA, direct gene-level) ───────
  message("\n[C] Cross-depot Spearman correlation")
  prefix_cor <- file.path(OUTPUT_DIR, paste0(resolution, "_crossdepot"))
  run_cross_depot_correlation(om, sq, prefix_cor)

  message("\nFinished resolution: ", resolution)
}

message("\nAll done. Results in: ", OUTPUT_DIR)
