# Sys.unsetenv("GITHUB_PAT")
# devtools::install_github(repo = "https://github.com/livnatje/DIALOGUE")

suppressPackageStartupMessages({
  library(zellkonverter)         # readH5AD → SingleCellExperiment
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(DESeq2)
  library(tidyverse)
  library(devtools)
  library(nnls)
  library(DIALOGUE)
})

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("sparseMatrixStats")

# __________________________________________________________________________
# Run a toy example
# __________________________________________________________________________
rA <- readRDS("/Library/Frameworks/R.framework/Versions/4.5-x86_64/Resources/library/DIALOGUE/data/test.example.rds")

summary(rA)


setwd("~/Desktop/") # use the path to where the DIALOGUE package is

param <- DLG.get.param(
  k = 2,
  results.dir = "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq/results/DIALOGUE/toy",
  n.genes = 100,
  conf = "cellQ",
  pheno = "pathology"
  )

R <- DIALOGUE.run(rA = rA,main = "toyExample", param = param, plot.flag = T)

summary(R)

# The gene sets in MCP1 and MCP2
summary(unlist(R$MCPs,recursive = F))

# The genes up-regulated in cell type A in MCP1
R$MCPs$MCP1$A.down

# The MCPs' association with a specific phenotype of interest, given as phenoZ in DIALOGUE.run.
# The magnitude of the value is -log10(p-value)
# and its sign denotes whether the association is positive or negative.
print(round(R$phenoZ,2))

# __________________________________________________________________________
# Pairwise DIALOGUE: all tissue × annotation-column × cell-type-pair combos
# Output: results/DIALOGUE/{tissue}/{ct_col}/{ct1}.vs.{ct2}/
# __________________________________________________________________________
suppressPackageStartupMessages({
  library(irlba)
  library(sparseMatrixStats)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(patchwork)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Global parameters ─────────────────────────────────────────────────────────
MIN_CELLS_PER_DONOR <- 10   # min cells per donor to include that donor
N_HVG               <- 3000
N_PCS               <- 30
MIN_DONORS          <- 8
N_MCP               <- 6   # MCPs to find per pair
# NULL shows all gene labels; set to a character vector to show only those genes
GENES_TO_LABEL      <- NULL

INPUT_ROOT  <- "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq/inputs/anndatas/"
OUTPUT_ROOT <- "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq/results/DIALOGUE/"

# Per-tissue config: h5ad path, participant column, sex availability, and
# which annotation columns to iterate over.
# SubQ obs has 'cell_type' / 'cell_type1' (no trailing underscore); add
# 'cell_type_' / 'cell_type2_' here if those columns exist in your h5ad.
TISSUE_CONFIGS <- list(
  omentum = list(
    h5ad        = paste0(INPUT_ROOT, "omentum_reannotated_finely_truncated.h5ad"),
    participant = "common_participant",
    has_sex     = TRUE,
    ct_cols     = c("cell_type_"),
    # Per-column remapping applied before building DIALOGUE objects.
    # Each entry is a named list of target_label -> grepl pattern.
    # Unmatched cell types are kept as-is.
    ct_remap    = list(
      cell_type2_ = c(
        "CD4_T_cells" = "CD4",
        "CD8_T_cells" = "CD8",
        "NK_cells"    = "NK"
      )
    )
  )#,
  # subq = list(
  #   h5ad        = paste0(INPUT_ROOT, "sq_reannotated_finely_truncated.h5ad"),
  #   participant = "common_participant",
  #   has_sex     = TRUE,
  #   ct_cols     = c("cell_type_") #, "cell_type2_"
  # )
)

# ── Helper functions ──────────────────────────────────────────────────────────

# Seurat-style LogNormalize: log1p(counts / colSum * 10 000).
log_normalize <- function(counts_mat) {
  lib_size <- Matrix::colSums(counts_mat)
  lib_size[lib_size == 0] <- 1
  normed <- Matrix::t(Matrix::t(counts_mat) / lib_size) * 1e4
  as.matrix(log1p(normed))
}

# Build one DIALOGUE cell.type object.
# check_sex: if TRUE and a 'sex' column exists, skip cell types with only 1 sex.
make_dialogue_ct <- function(sce, ct_name,
                             cell_type_col   = "cell_type_",
                             participant_col = "participant",
                             check_sex       = TRUE,
                             min_cells       = MIN_CELLS_PER_DONOR,
                             n_hvg           = N_HVG,
                             n_pcs           = N_PCS,
                             min_donors      = MIN_DONORS) {

  mask        <- colData(sce)[[cell_type_col]] == ct_name
  sce_ct      <- sce[, mask]
  meta        <- as.data.frame(colData(sce_ct))
  samples_vec <- meta[[participant_col]]

  donor_n     <- table(samples_vec)
  keep_donors <- names(donor_n)[donor_n >= min_cells]
  if (length(keep_donors) < min_donors) {
    message("  Skipping '", ct_name, "': only ", length(keep_donors),
            " donors >= ", min_cells, " cells")
    return(NULL)
  }

  keep_cells  <- samples_vec %in% keep_donors
  sce_ct      <- sce_ct[, keep_cells]
  meta        <- meta[keep_cells, , drop = FALSE]
  samples_vec <- meta[[participant_col]]

  if (length(unique(meta[["cohort"]])) < 2) {
    message("  Skipping '", ct_name, "': only 1 cohort passes cell threshold")
    return(NULL)
  }
  if (check_sex && "sex" %in% names(meta) && length(unique(meta[["sex"]])) < 2) {
    message("  Skipping '", ct_name, "': only 1 sex passes cell threshold")
    return(NULL)
  }

  raw     <- assay(sce_ct, "counts")
  tpm_mat <- log_normalize(raw)
  colnames(tpm_mat) <- colnames(sce_ct)

  donor_avg  <- sapply(keep_donors, function(s)
    rowMeans(tpm_mat[, samples_vec == s, drop = FALSE]))
  gene_sd    <- matrixStats::rowSds(donor_avg)
  keep_genes <- !is.na(gene_sd) & gene_sd > 0
  tpm_mat    <- tpm_mat[keep_genes, ]
  raw_filt   <- raw[keep_genes, ]

  gene_var  <- sparseMatrixStats::rowVars(raw_filt)
  hvg_idx   <- order(gene_var, decreasing = TRUE)[seq_len(min(n_hvg, nrow(raw_filt)))]
  hvg_sc    <- t(scale(t(tpm_mat[hvg_idx, ])))
  hvg_sc[hvg_sc >  10] <-  10
  hvg_sc[hvg_sc < -10] <- -10

  n_pcs_use <- min(n_pcs, nrow(hvg_sc) - 1, ncol(hvg_sc) - 1)
  pca_res   <- irlba::prcomp_irlba(t(hvg_sc), n = n_pcs_use,
                                   center = FALSE, scale. = FALSE)
  X_mat     <- pca_res$x
  rownames(X_mat) <- colnames(sce_ct)
  colnames(X_mat) <- paste0("PC", seq_len(n_pcs_use))

  cellQ_vec <- log10(Matrix::colSums(raw_filt) + 1)

  meta_out <- data.frame(
    row.names      = colnames(sce_ct),
    donor          = samples_vec,
    cohort_numeric = as.integer(meta[["cohort"]] == "Unhealthy"),
    cellQ          = cellQ_vec
  )
  if (check_sex && "sex" %in% names(meta))
    meta_out$sex <- meta[["sex"]]

  message("  '", ct_name, "': ", ncol(sce_ct), " cells, ",
          length(keep_donors), " donors, ", sum(keep_genes), " genes")

  make.cell.type(name = ct_name, tpm = tpm_mat, samples = samples_vec,
                 X = X_mat, metadata = meta_out, cellQ = cellQ_vec)
}

# DIALOGUE strips non-alphanumeric characters when building gene-list keys.
sanitize_ct <- function(s) gsub("[^A-Za-z0-9]", "", s)

# Show which donors have cells in each cell type and return whether enough
# donors are shared across ALL cell types to run DIALOGUE.
# verbose = TRUE prints the full donor × cell-type presence matrix.
check_donor_coverage <- function(ct_list, min_donors = 5L, verbose = TRUE) {
  all_donors <- sort(unique(unlist(lapply(ct_list, function(ct) unique(ct@samples)))))
  ct_names   <- names(ct_list)

  mat <- outer(all_donors, ct_names,
               FUN = Vectorize(function(d, ct) d %in% ct_list[[ct]]@samples))
  dimnames(mat) <- list(donor = all_donors, cell_type = ct_names)

  n_per_donor <- rowSums(mat)
  n_shared    <- sum(n_per_donor == length(ct_names))

  if (verbose) {
    cat("\n  Donor × cell-type coverage (1 = present):\n")
    print(mat * 1L)
    cat("\n  Cell types per donor:\n")
    print(sort(n_per_donor))
    cat("\n  Donors in ALL", length(ct_names), "cell types:", n_shared, "\n")
    cat(  "  Donors in >= half:          ",
          sum(n_per_donor >= length(ct_names) / 2), "\n\n")
  } else {
    message("    Shared donors (all cell types): ", n_shared)
  }

  n_shared >= min_donors
}

# Per-donor MCP score: mean(up tpm) - mean(down tpm) averaged over cells.
compute_donor_scores <- function(ct_obj, up_genes, down_genes) {
  tpm     <- ct_obj@tpm
  samples <- ct_obj@samples
  meta    <- ct_obj@metadata

  up_use   <- intersect(up_genes,   rownames(tpm))
  down_use <- intersect(down_genes, rownames(tpm))

  score_up   <- if (length(up_use)   > 0) colMeans(tpm[up_use,   , drop = FALSE]) else rep(0, ncol(tpm))
  score_down <- if (length(down_use) > 0) colMeans(tpm[down_use, , drop = FALSE]) else rep(0, ncol(tpm))
  cell_score <- score_up - score_down

  donors <- sort(unique(samples))
  data.frame(
    donor     = donors,
    score     = sapply(donors, function(d) mean(cell_score[samples == d])),
    cohort    = sapply(donors, function(d) {
      num <- meta$cohort_numeric[which(samples == d)[1]]
      if (is.na(num)) "Unknown" else if (num == 1L) "Ob+T2D" else "Ob"
    }),
    cell_type = ct_obj@name,
    stringsAsFactors = FALSE
  )
}

# Generate all three plots for one DIALOGUE result.
# R         : DIALOGUE.run() output
# ct_subset : named list of cell.type objects for this pair
# pair_dir  : root output directory for this pair (already created)
plot_mcp_results <- function(R, ct_subset, pair_dir) {
  PLOT_DIR    <- file.path(pair_dir, "plots")
  dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
  ct_names    <- names(ct_subset)
  cohort_cols <- c(Ob = "#A1A1A1", "Ob+T2D" = "#FDC707")

  for (mcp in names(R$MCPs)) {
    mcp_genes <- R$MCPs[[mcp]]

    # ── donor scores ──────────────────────────────────────────────────────────
    score_list <- lapply(ct_names, function(ct) {
      key    <- sanitize_ct(ct)
      up_g   <- mcp_genes[[paste0(key, ".up")]]   %||% character(0)
      down_g <- mcp_genes[[paste0(key, ".down")]] %||% character(0)
      if (length(up_g) == 0 && length(down_g) == 0) return(NULL)
      compute_donor_scores(ct_subset[[ct]], up_g, down_g)
    })
    score_list <- Filter(Negate(is.null), score_list)
    if (length(score_list) == 0) next

    scores_df <- bind_rows(score_list) %>%
      mutate(cohort = recode(cohort, Healthy = "Ob", Unhealthy = "Ob+T2D"),
             cohort = factor(cohort, levels = c("Ob", "Ob+T2D")))

    # ── Plot 1: boxplot + jitter by cohort, faceted by cell type ─────────────
    p_cohort <- ggplot(scores_df, aes(x = cohort, y = score, fill = cohort)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.45) +
      geom_jitter(aes(color = cohort), width = 0.12, size = 2.2, alpha = 0.85) +
      facet_wrap(~ cell_type, scales = "free_y") +
      scale_fill_manual(values  = cohort_cols) +
      scale_color_manual(values = c(Ob = "#A1A1A1", "Ob+T2D" = "#FDC707")) +
      labs(title = paste(mcp, "— donor scores by cohort"),
           x = NULL, y = "MCP score (up − down gene mean)") +
      theme_classic(base_size = 16) +
      theme(legend.position  = "none",
            strip.background = element_blank(),
            strip.text       = element_text(face = "bold", size = 16),
            axis.text        = element_text(size = 14),
            axis.title       = element_text(size = 15),
            plot.title       = element_text(size = 16, face = "bold"))

    ggsave(file.path(PLOT_DIR, paste0(mcp, "_cohort_by_celltype.pdf")),
           p_cohort, width = 1.2 * length(score_list), height = 6.5)

    # ── Plot 2: cross-cell-type score correlation ─────────────────────────────
    if (length(score_list) >= 2) {
      ct_present <- unique(scores_df$cell_type)
      corr_plots <- lapply(combn(ct_present, 2, simplify = FALSE), function(pair) {
        wide <- scores_df %>%
          filter(cell_type %in% pair) %>%
          dplyr::select(donor, score, cohort, cell_type) %>%
          pivot_wider(names_from = cell_type, values_from = score) %>%
          drop_na()
        ct1 <- pair[1]; ct2 <- pair[2]
        if (nrow(wide) < 3) return(NULL)
        r_val <- round(cor(wide[[ct1]], wide[[ct2]], use = "complete.obs"), 2)
        ggplot(wide, aes(x = .data[[ct1]], y = .data[[ct2]], color = cohort)) +
          geom_point(size = 2.5, alpha = 0.85) +
          geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                      color = "grey35", linewidth = 0.7, alpha = 0.15) +
          scale_color_manual(values = cohort_cols) +
          annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
                   label = paste0("r = ", r_val), size = 3.5) +
          labs(title = paste(mcp, ":", ct1, "vs", ct2),
               x = paste(ct1, "score"), y = paste(ct2, "score"), color = NULL) +
          theme_classic(base_size = 11) +
          theme(legend.position = "bottom")
      })
      corr_plots <- Filter(Negate(is.null), corr_plots)
      if (length(corr_plots) > 0) {
        n_col <- min(3L, length(corr_plots))
        ggsave(
          file.path(PLOT_DIR, paste0(mcp, "_cross_celltype_correlation.pdf")),
          wrap_plots(corr_plots, ncol = n_col),
          width = 5 * n_col, height = 5 * ceiling(length(corr_plots) / n_col)
        )
      }
    }

    # ── Plot 3: per-donor gene z-score heatmap ────────────────────────────────
    heat_blocks <- list()
    for (ct in ct_names) {
      key       <- sanitize_ct(ct)
      ct_obj    <- ct_subset[[ct]]
      tpm_ct    <- ct_obj@tpm
      samp      <- ct_obj@samples
      donors_ct <- sort(unique(samp))

      donor_means_ct <- sapply(donors_ct, function(d)
        rowMeans(tpm_ct[, samp == d, drop = FALSE]))
      rownames(donor_means_ct) <- rownames(tpm_ct)
      colnames(donor_means_ct) <- donors_ct

      for (direction in c("up", "down")) {
        genes     <- mcp_genes[[paste0(key, ".", direction)]]
        if (is.null(genes) || length(genes) == 0) next
        genes_use <- intersect(genes, rownames(tpm_ct))
        if (length(genes_use) == 0) next
        block_z <- t(scale(t(donor_means_ct[genes_use, , drop = FALSE])))
        heat_blocks[[paste0(ct, " (", direction, ")")]] <- list(mat = block_z)
      }
    }
    GENES_TO_LABEL <- unlist(strsplit(
      "GPX4, PRDX5, TXNIP, SELENOP, IL6, TNFAIP6, IRF1, ITLN2, TXNIP, COL3A1, NR2F2, ADAMTS4, ADAMTS9, BMP2, COL4A1, TXNIP, PTGES, PLA2G5, SFRP4, PTX3, SERPINE1, XBP1",
      ", "))
    
    if (length(heat_blocks) > 0) {
      cohort_map <- local({
        m <- character(0)
        for (ct in ct_names) {
          meta <- ct_subset[[ct]]@metadata
          new  <- setNames(ifelse(meta$cohort_numeric == 1L, "Ob+T2D", "Ob"),
                           meta$donor)
          m <- c(m, new[setdiff(names(new), names(m))])
        }
        m
      })
      
      all_donors   <- sort(unique(unlist(lapply(heat_blocks, function(b) colnames(b$mat)))))
      donor_cohort <- cohort_map[all_donors]
      d_ord        <- order(donor_cohort)
      all_donors   <- all_donors[d_ord]
      donor_cohort <- donor_cohort[d_ord]

      sex_map <- local({
        m <- character(0)
        for (ct in ct_names) {
          meta <- ct_subset[[ct]]@metadata
          if (!"sex" %in% names(meta)) next
          new <- tapply(meta$sex, meta$donor, function(x) as.character(x[1]))
          m   <- c(m, new[setdiff(names(new), names(m))])
        }
        if (length(m) == 0) NULL else m
      })
      donor_sex <- if (!is.null(sex_map)) sex_map[all_donors] else NULL

      mats <- lapply(heat_blocks, function(blk) {
        m       <- matrix(NA_real_, nrow = nrow(blk$mat), ncol = length(all_donors),
                          dimnames = list(rownames(blk$mat), all_donors))
        dst_idx <- match(colnames(blk$mat), all_donors)
        valid   <- !is.na(dst_idx)
        if (any(valid)) m[, dst_idx[valid]] <- blk$mat[, which(valid), drop = FALSE]
        m
      })
      
      full_mat  <- do.call(rbind, mats)
      row_split <- factor(rep(names(heat_blocks), sapply(mats, nrow)),
                          levels = names(heat_blocks))
      
      leg_gp <- list(title_gp = gpar(fontsize = 18, fontface = "bold"),
                     labels_gp = gpar(fontsize = 20))
      if (!is.null(donor_sex)) {
        row_ann <- rowAnnotation(
          Cohort = donor_cohort,
          Sex    = donor_sex,
          col    = list(
            Cohort = c(Ob = "#A1A1A1", "Ob+T2D" = "#FDC707"),
            Sex    = c(F = "#E78AC3", M = "#66C2A5")
          ),
          show_annotation_name    = FALSE,
          annotation_legend_param = list(Cohort = leg_gp, Sex = leg_gp)
        )
      } else {
        row_ann <- rowAnnotation(
          Cohort = donor_cohort,
          col    = list(Cohort = c(Ob = "#A1A1A1", "Ob+T2D" = "#FDC707")),
          show_annotation_name    = FALSE,
          annotation_legend_param = list(Cohort = leg_gp)
        )
      }
      col_fun    <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
      genes_all <- rownames(full_mat)
      if (is.null(GENES_TO_LABEL)) {
        mark_ann  <- NULL
        show_cols <- TRUE
      } else {
        at_idx   <- which(genes_all %in% GENES_TO_LABEL)
        mark_ann <- HeatmapAnnotation(
          mark = anno_mark(
            at        = at_idx,
            labels    = genes_all[at_idx],
            which     = "column",
            side      = "bottom",
            labels_gp = gpar(fontsize = 20),
            lines_gp  = gpar(col = "black", lwd = 0.8),
            padding   = unit(2, "mm"),
            link_width = unit(8, "mm")
          )
        )
        show_cols <- FALSE
      }
      
      ht <- Heatmap(t(full_mat),           # donors × genes
                    name                  = "z-score",
                    col                   = col_fun,
                    left_annotation       = row_ann,
                    bottom_annotation     = mark_ann,
                    column_split          = row_split,
                    column_title_rot      = 90,
                    column_title_gp       = gpar(fontsize = 22),
                    cluster_rows          = TRUE,
                    cluster_columns       = FALSE,     # no dendrogram within gene blocks
                    cluster_column_slices = FALSE,     # preserve block order
                    show_row_names        = FALSE,
                    show_column_names     = show_cols,
                    row_names_gp          = gpar(fontsize = 26),
                    column_names_gp       = gpar(fontsize = 11),
                    na_col                = "grey90",
                    border                = TRUE,
                    border_gp             = gpar(col = "grey30", lwd = 0.8),
                    heatmap_legend_param  = list(
                      title          = "z-score",
                      title_gp       = gpar(fontsize = 16, fontface = "bold"),
                      labels_gp      = gpar(fontsize = 18),
                      legend_height  = unit(8, "cm"),
                      grid_width     = unit(0.6, "cm")
                    ),
                    clustering_distance_rows = "pearson",
                    clustering_method_rows   = "ward.D2",
      )
      
      pdf(file.path(PLOT_DIR, paste0(mcp, "_gene_zscore_heatmap.pdf")),
          width  = 22,
          height = 10)
      draw(ht, column_title    = paste(mcp, "— per-donor gene z-score"),
           column_title_gp = gpar(fontsize = 11, fontface = "bold"))
      dev.off()
    }

    message("    Saved plots for ", mcp)
  }
  message("  Plots written to: ", PLOT_DIR)
}

# ── Main loop ─────────────────────────────────────────────────────────────────
summary_records <- list()   # accumulates phenoZ rows for pairwise summary
all_results     <- list()   # accumulates all-celltypes results per tissue/ct_col

for (tissue_name in names(TISSUE_CONFIGS)) {
  cfg <- TISSUE_CONFIGS[[tissue_name]]

  if (!file.exists(cfg$h5ad)) {
    message("Skipping ", tissue_name, ": h5ad not found at\n  ", cfg$h5ad)
    next
  }

  message("\n", strrep("=", 60))
  message("Loading ", tissue_name)
  sce <- readH5AD(cfg$h5ad, use_hdf5 = FALSE)

  for (ct_col in cfg$ct_cols) {
    if (!ct_col %in% names(colData(sce))) {
      message("Skipping ", tissue_name, " / ", ct_col, ": column not present")
      next
    }

    # Apply optional subtype-merging remap for this column.
    # A remapped column is written into colData so make_dialogue_ct can use it;
    # the original column is never modified.
    effective_col <- ct_col
    remap <- cfg$ct_remap[[ct_col]]
    if (!is.null(remap)) {
      raw_vals    <- as.character(colData(sce)[[ct_col]])
      merged_vals <- raw_vals
      for (target in names(remap)) {
        pattern             <- remap[[target]]
        merged_vals[grepl(pattern, raw_vals, ignore.case = TRUE)] <- target
      }
      effective_col                    <- paste0(ct_col, "merged_")
      colData(sce)[[effective_col]]    <- merged_vals
      n_merged <- sum(merged_vals != raw_vals)
      message("  Remap applied to ", ct_col, ": ", n_merged,
              " cells reassigned (merged subtypes: ",
              paste(names(remap), collapse = ", "), ")")
    }

    cell_types <- sort(unique(na.omit(as.character(colData(sce)[[effective_col]]))))
    message("\n", strrep("-", 50))
    message(tissue_name, " / ", ct_col, ": ", length(cell_types),
            " cell types — building DIALOGUE objects")

    # Build all cell-type objects once; reuse across every pair.
    ct_objects <- lapply(setNames(cell_types, cell_types), function(ct) {
      make_dialogue_ct(sce, ct,
                       cell_type_col   = effective_col,
                       participant_col = cfg$participant,
                       check_sex       = cfg$has_sex)
    })
    ct_objects <- Filter(Negate(is.null), ct_objects)

    if (length(ct_objects) < 2) {
      message("  Fewer than 2 cell types passed QC — skipping")
      next
    }

    # ── All-cell-types run ───────────────────────────────────────────────────
    {
      run_dir      <- file.path(OUTPUT_ROOT, tissue_name, ct_col, "all_celltypes")
      dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
      # Use our own cache file — DIALOGUE's internal DIALOGUE2_*.rds is saved
      # before MCPs are computed and will be missing R$MCPs on reload.
      our_rds_path <- file.path(run_dir, "R_complete.rds")

      message("\n  [", tissue_name, " / ", ct_col, "] All-celltypes run (",
              length(ct_objects), " cell types)")

      if (!check_donor_coverage(ct_objects, min_donors = 5L, verbose = TRUE)) {
        message("    Skipping: fewer than 5 donors shared across all cell types")
        next
      }

      if (file.exists(our_rds_path)) {
        message("    Cached result found — loading")
        R_all <- readRDS(our_rds_path)
      } else {
        param_all <- DLG.get.param(
          k           = N_MCP,
          results.dir = run_dir,
          n.genes     = 200,
          conf        = "cellQ",
          pheno       = "cohort_numeric",
          bypass.emp  = TRUE
        )
        R_all <- tryCatch(
          DIALOGUE.run(rA = ct_objects, main = "all_celltypes",
                       param = param_all, plot.flag = FALSE),
          error = function(e) {
            message("    DIALOGUE.run failed: ", conditionMessage(e))
            NULL
          }
        )
        # Save our own copy immediately after DIALOGUE.run returns, while the
        # complete object (including MCPs) is still in memory.
        if (!is.null(R_all)) saveRDS(R_all, our_rds_path)
      }

      if (!is.null(R_all) && length(R_all$MCPs) > 0) {
        message("    phenoZ: ",
                paste(names(R_all$phenoZ), round(R_all$phenoZ, 2),
                      sep = "=", collapse = " | "))
        plot_mcp_results(R_all, ct_objects, run_dir)
        grp_key <- paste(tissue_name, ct_col, sep = "///")
        all_results[[grp_key]] <- list(
          tissue      = tissue_name,
          ct_col      = ct_col,
          n_celltypes = length(ct_objects),
          R           = R_all
        )
      } else {
        message("    No MCPs found")
      }
    }

    # # ── Pairwise runs ────────────────────────────────────────────────────────
    # pairs <- combn(names(ct_objects), 2, simplify = FALSE)
    # message("\n  [", tissue_name, " / ", ct_col, "] Pairwise: ",
    #         length(ct_objects), " cell types → ", length(pairs), " pairs")
    # 
    # for (pair in pairs) {
    #   ct1 <- pair[1]; ct2 <- pair[2]
    #   pair_label <- paste0(ct1, "_vs_", ct2)           # used as DIALOGUE 'main'
    #   pair_dir   <- file.path(OUTPUT_ROOT, tissue_name,
    #                           ct_col, paste0(ct1, ".vs.", ct2))
    #   dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
    # 
    #   message("\n  [", tissue_name, " / ", ct_col, "] ", ct1, " vs ", ct2)
    # 
    #   ct_subset <- ct_objects[c(ct1, ct2)]
    #
    #   if (!check_donor_coverage(ct_subset, min_donors = 5L, verbose = FALSE)) {
    #     message("    Skipping: fewer than 5 donors shared between ", ct1, " and ", ct2)
    #     next
    #   }
    #
    #   # Use our own cache — DIALOGUE's DIALOGUE2_*.rds lacks R$MCPs (saved
    #   # before MCPs are computed internally).
    #   our_rds_path <- file.path(pair_dir, "R_complete.rds")
    #   if (file.exists(our_rds_path)) {
    #     message("    Cached result found — loading")
    #     R <- readRDS(our_rds_path)
    #   } else {
    #     param <- DLG.get.param(
    #       k           = N_MCP,
    #       results.dir = pair_dir,
    #       n.genes     = 200,
    #       conf        = "cellQ",
    #       pheno       = "cohort_numeric",
    #       bypass.emp  = TRUE
    #     )
    #     R <- tryCatch(
    #       DIALOGUE.run(rA = ct_subset, main = pair_label,
    #                    param = param, plot.flag = FALSE),
    #       error = function(e) {
    #         message("    DIALOGUE.run failed: ", conditionMessage(e))
    #         NULL
    #       }
    #     )
    #     if (!is.null(R)) saveRDS(R, our_rds_path)
    #   }
    # 
    #   if (is.null(R) || length(R$MCPs) == 0) {
    #     message("    No MCPs found — skipping plots")
    #     next
    #   }
    # 
    #   message("    phenoZ: ",
    #           paste(names(R$phenoZ), round(R$phenoZ, 2), sep = "=", collapse = " | "))
    # 
    #   plot_mcp_results(R, ct_subset, pair_dir)
    # 
    #   # Accumulate phenoZ for cross-pair summary
    #   pz <- R$phenoZ
    #   summary_records <- c(summary_records, lapply(names(pz), function(mcp_name) {
    #     data.frame(tissue = tissue_name, ct_col = ct_col,
    #                ct1 = ct1, ct2 = ct2,
    #                pair   = paste0(ct1, " vs ", ct2),
    #                mcp    = mcp_name,
    #                phenoZ = pz[[mcp_name]],
    #                stringsAsFactors = FALSE)
    #   }))
    # }
  }
}

# ── Cross-pair summary: which pairs associate most strongly with cohort ────────
if (length(summary_records) > 0) {
  summary_df <- bind_rows(summary_records)
  TOP_N      <- 20L

  for (grp_key in unique(paste(summary_df$tissue, summary_df$ct_col, sep = "///"))) {
    tis     <- strsplit(grp_key, "///")[[1L]][1L]
    col     <- strsplit(grp_key, "///")[[1L]][2L]
    sub     <- filter(summary_df, tissue == tis, ct_col == col)
    out_dir <- file.path(OUTPUT_ROOT, tis, col)

    # ── Flat CSV (all pairs, all MCPs, sorted by |phenoZ|) ───────────────────
    write.csv(arrange(sub, desc(abs(phenoZ))),
              file.path(out_dir, "phenoZ_summary.csv"), row.names = FALSE)

    # ── Console table ─────────────────────────────────────────────────────────
    cat("\n", strrep("=", 60), "\n", sep = "")
    cat("phenoZ summary: ", tis, " / ", col, "\n\n", sep = "")
    sub %>%
      group_by(pair) %>%
      summarise(max_abs_phenoZ = max(abs(phenoZ)),
                best_mcp       = mcp[which.max(abs(phenoZ))],
                best_phenoZ    = phenoZ[which.max(abs(phenoZ))],
                .groups = "drop") %>%
      arrange(desc(max_abs_phenoZ)) %>%
      print(n = Inf)

    # ── Bar chart: top N pairs by max |phenoZ| ────────────────────────────────
    # Pairs ordered by their maximum |phenoZ| across MCPs; one facet per MCP.
    top_pairs <- sub %>%
      group_by(pair) %>%
      summarise(max_abs = max(abs(phenoZ)), .groups = "drop") %>%
      slice_max(max_abs, n = TOP_N, with_ties = FALSE) %>%
      arrange(max_abs) %>%   # ascending so highest is at top of horizontal bars
      pull(pair)

    p_bar <- sub %>%
      filter(pair %in% top_pairs) %>%
      mutate(pair = factor(pair, levels = top_pairs)) %>%
      ggplot(aes(x = phenoZ, y = pair, fill = phenoZ > 0)) +
        geom_col() +
        geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey30") +
        facet_wrap(~ mcp, ncol = length(unique(sub$mcp))) +
        scale_fill_manual(
          values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
          labels = c(`TRUE` = "Unhealthy↑", `FALSE` = "Healthy↑"),
          name   = NULL
        ) +
        labs(
          title    = paste0(tis, " / ", col,
                            "  —  phenoZ  (top ", TOP_N, " pairs)"),
          subtitle = "Sorted by max |phenoZ| across MCPs. phenoZ = −log10(p) × sign",
          x = "phenoZ", y = NULL
        ) +
        theme_classic(base_size = 10) +
        theme(strip.background = element_blank(),
              strip.text       = element_text(face = "bold"),
              legend.position  = "bottom")

    ggsave(file.path(out_dir, "phenoZ_top_pairs_bar.pdf"), p_bar,
           width  = 4 * length(unique(sub$mcp)),
           height = max(4, length(top_pairs) * 0.35 + 2))

    # ── Heatmap: all pairs × MCPs ─────────────────────────────────────────────
    # Rows sorted by max |phenoZ| descending so the strongest pairs are at top.
    wide <- sub %>%
      pivot_wider(id_cols = pair, names_from = mcp, values_from = phenoZ) %>%
      column_to_rownames("pair")
    row_order <- order(apply(abs(wide), 1L, max, na.rm = TRUE), decreasing = TRUE)
    wide      <- wide[row_order, , drop = FALSE]

    lim        <- max(abs(unlist(wide)), na.rm = TRUE)
    col_fun_ph <- colorRamp2(c(-lim, 0, lim), c("#2166AC", "white", "#B2182B"))

    ht_sum <- Heatmap(
      as.matrix(wide),
      name                 = "phenoZ",
      col                  = col_fun_ph,
      cluster_rows         = FALSE,
      cluster_columns      = FALSE,
      show_row_names       = TRUE,
      show_column_names    = TRUE,
      row_names_gp         = gpar(fontsize = max(5, min(9, 200 / nrow(wide)))),
      column_names_gp      = gpar(fontsize = 9, fontface = "bold"),
      heatmap_legend_param = list(title = "phenoZ"),
      column_title         = paste0(tis, " / ", col,
                                    "  —  phenoZ: all pairs × MCPs",
                                    "\n(sorted by max |phenoZ|)"),
      column_title_gp      = gpar(fontsize = 11, fontface = "bold"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        val <- as.matrix(wide)[i, j]
        if (!is.na(val))
          grid.text(sprintf("%.1f", val), x, y, gp = gpar(fontsize = 6))
      }
    )

    pdf(file.path(out_dir, "phenoZ_all_pairs_heatmap.pdf"),
        width  = max(5, ncol(wide) * 1.5 + 3),
        height = max(4, nrow(wide) * 0.28 + 2))
    draw(ht_sum)
    dev.off()

    message("Summary saved for ", tis, " / ", col)
  }
}

# ── All-cell-types phenoZ summary ─────────────────────────────────────────────
if (length(all_results) > 0) {
  pz_df <- bind_rows(lapply(names(all_results), function(key) {
    r  <- all_results[[key]]
    pz <- r$R$phenoZ
    if (length(pz) == 0) return(NULL)
    mcp_names <- if (is.matrix(pz) && !is.null(rownames(pz)) && !is.null(colnames(pz)))
      paste(rep(colnames(pz), each = nrow(pz)), rep(rownames(pz), times = ncol(pz)), sep = "__")
    else if (!is.null(names(pz)))
      names(pz)
    else
      paste0("MCP", seq_along(pz))
    data.frame(tissue      = if (length(r$tissue)      == 0) NA_character_ else r$tissue,
               ct_col      = if (length(r$ct_col)      == 0) NA_character_ else r$ct_col,
               n_celltypes = if (length(r$n_celltypes) == 0) NA_integer_   else r$n_celltypes,
               mcp         = mcp_names,
               phenoZ      = as.numeric(pz),
               stringsAsFactors = FALSE)
  }))

  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("All-cell-types DIALOGUE — phenoZ summary\n\n")
  print(arrange(pz_df, tissue, ct_col, desc(abs(phenoZ)))) #, n = Inf)

  write.csv(pz_df,
            file.path(OUTPUT_ROOT, "all_celltypes_phenoZ_summary.csv"),
            row.names = FALSE)

  # Bar chart: one facet per tissue/ct_col combination, bars = MCPs
  p_all <- pz_df %>%
    mutate(facet_label = paste0(tissue, "\n", ct_col,
                                "\n(", n_celltypes, " cell types)")) %>%
    ggplot(aes(x = phenoZ, y = mcp, fill = phenoZ > 0)) +
      geom_col() +
      geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey30") +
      facet_wrap(~ facet_label, ncol = length(all_results)) +
      scale_fill_manual(
        values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"),
        labels = c(`TRUE` = "Unhealthy↑", `FALSE` = "Healthy↑"),
        name   = NULL
      ) +
      labs(
        title    = "All-cell-types DIALOGUE — MCP association with cohort",
        subtitle = "phenoZ = −log10(p) × sign;  positive = enriched in Unhealthy",
        x = "phenoZ", y = NULL
      ) +
      theme_classic(base_size = 11) +
      theme(strip.background = element_blank(),
            strip.text       = element_text(face = "bold"),
            legend.position  = "bottom",
            panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.3))

  ggsave(
    file.path(OUTPUT_ROOT, "all_celltypes_phenoZ_comparison.pdf"),
    p_all,
    width  = max(4, 3.5 * length(all_results)),
    height = max(3, N_MCP * 1.7 + 2)
  )

  message("All-cell-types summary saved to ", OUTPUT_ROOT)
}

message("\nAll done.")



df <- data.frame(lapply(R_all$MCPs[["MCP3"]], function(x) {
  x <- unlist(x)
  length(x) <- max(lengths(R_all$MCPs[["MCP3"]]))
  return(x)
}))
