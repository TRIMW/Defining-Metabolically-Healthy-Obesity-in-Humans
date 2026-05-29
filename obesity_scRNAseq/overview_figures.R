## Overview figures: Healthy vs Unhealthy DE across cell types and tissues
## Panels per DE method:
##   A  – DEG count barplot
##   B  – top-gene bubble plot
##   C  – faceted volcano plots
##   D  – top-DEG logFC heatmap (ggplot)
##   E  – Hallmark NES heatmap  (ComplexHeatmap, with dendrogram)
##   F  – KEGG NES heatmap      (ComplexHeatmap, with dendrogram)

library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

setwd("/Users/willtrim/Documents/projs/Wills_analysis")

# ── Paths ─────────────────────────────────────────────────────────────────────
DE_ROOT  <- "obesity_scRNAseq/results/DE"
OUT_DIR  <- "obesity_scRNAseq/results/overview_figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TISSUES   <- c("omentum", "subq")
CT_LEVEL  <- "cell_type_"   # "cell_type_" = broad  |  "cell_type2_" = fine
FDR_CUT   <- 0.05
L2FC_CUT  <- 0.5
TOP_GENES <- 5              # genes per cell type in bubble / logFC heatmap
TOP_PATHS <- 30             # max pathways in NES heatmaps

# ── Cell-type colour palette (edit to match your UMAPs) ───────────────────────
CT_COLOURS <- c(
  B_cells          = "#E41A1C",
  ECs              = "#377EB8",
  FAPs             = "#4DAF4A",
  IGFBP2_cells     = "#984EA3",
  Mesothelial      = "#FF7F00",
  Myeloid_cells    = "#A65628",
  Myofibroblasts   = "#F781BF",
  NK_cells         = "#999999",
  `Pre-adipocytes` = "#66C2A5",
  Pre_adipocytes   = "#66C2A5",
  Plasma_cells     = "#FC8D62",
  T_cells          = "#8DA0CB"
)

ct_colour <- function(ct_names) {
  cols <- CT_COLOURS[ct_names]
  cols[is.na(cols)] <- "grey70"
  unname(cols)
}

# ── Colour constants ───────────────────────────────────────────────────────────
dir_colours <- c(Up = "#B2182B", Down = "#2166AC")
nes_col_fun <- colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B"))


# Ribosomal genes heavily overlap viral pathway gene sets → artifactual hits
viral_re <- regex(
  "(sars.cov|type 1 diabetes|autoimmune thyroid disease|toxoplasmosis|covid|coronavirus|influenza infection|Parkinson Disease|Amyotrohic lateral sclerosis|viral myocarditis|Gastric acid secretion|cardiomyopathy)", ignore_case = TRUE)


# ── Gene name lookup (Ensembl → symbol) ───────────────────────────────────────
gene_lookup <- read_csv(
  "obesity_scRNAseq/inputs/gProfiler_hsapiens_2026-04-29_19-37-55.csv",
  show_col_types = FALSE
) %>%
  filter(name != "" & !is.na(name)) %>%
  select(ensembl_id = initial_alias, gene_name = name) %>%
  distinct(ensembl_id, .keep_all = TRUE)

ensembl2gene <- setNames(gene_lookup$gene_name, gene_lookup$ensembl_id)

resolve_gene <- function(ensembl_ids,
                         lut = get("ensembl2gene", envir = globalenv())) {
  result <- lut[ensembl_ids]
  ifelse(is.na(result), ensembl_ids, result)
}

# ── Data loaders ──────────────────────────────────────────────────────────────
load_de <- function(method, tissue, ct_level = CT_LEVEL) {
  suffix <- switch(ct_level,
    "cell_type_"  = "_per_cell_type_.csv",
    "cell_type2_" = "_per_cell_type2_.csv"
  )
  fname <- switch(method,
    deseq2 = paste0("deseq2_results", suffix),
    edgeR  = paste0("edgeR_results",  suffix),
    MAST   = paste0("MAST_results",   suffix)
  )
  raw <- read_csv(file.path(DE_ROOT, tissue, fname), show_col_types = FALSE)

  if (method == "deseq2") {
    raw %>% transmute(gene = variable, logFC = log_fc, padj = adj_p_value,
                      cell_type)
  } else if (method == "edgeR") {
    raw %>% transmute(gene = resolve_gene(ensembl_id), logFC, padj = FDR,
                      cell_type)
  } else {
    raw %>% transmute(gene = primerid, logFC = coef, padj = FDR, cell_type)
  }
}

load_fgsea <- function(method, tissue, db, ct_level = CT_LEVEL) {
  method_dir <- switch(method,
    deseq2 = "deseq2", edgeR = "edger", MAST = "MAST"
  )
  dir_path <- file.path(DE_ROOT, tissue, "fgsea", ct_level, method_dir)
  files    <- list.files(dir_path, pattern = paste0("\\.", db, "\\.csv$"),
                         full.names = TRUE)
  map_dfr(files, function(f) {
    ct <- sub(paste0("\\.", db, "\\.csv$"), "", basename(f))
    read_csv(f, show_col_types = FALSE) %>% mutate(cell_type = ct)
  })
}

# ── Panel A: DEG count barplot ─────────────────────────────────────────────────
make_deg_count_plot <- function(method) {
  de_all <- map_dfr(TISSUES, function(tis) {
    load_de(method, tis) %>% mutate(tissue = tis)
  }) %>%
    filter(!is.na(padj)) %>%
    mutate(sig = padj < FDR_CUT,
           direction = if_else(logFC > 0, "Up", "Down"))

  counts <- de_all %>%
    filter(sig) %>%
    dplyr::count(tissue, cell_type, direction) %>%
    mutate(n = if_else(direction == "Down", -n, n),
           tissue = factor(tissue, levels = TISSUES))

  order_ct <- de_all %>%
    filter(sig) %>%
    dplyr::count(cell_type, name = "total") %>%
    arrange(dplyr::desc(total)) %>%
    pull("cell_type")

  counts <- counts %>%
    mutate(cell_type = factor(cell_type, levels = order_ct))

  ggplot(counts, aes(x = cell_type, y = n, fill = direction)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
    facet_wrap(~ tissue, ncol = 1, scales = "free_x",
               labeller = as_labeller(c(omentum = "Omentum",
                                        subq    = "Subcutaneous"))) +
    scale_fill_manual(values = dir_colours, name = NULL) +
    scale_y_continuous(
      breaks = c(-2000, -500, -100, -20, 0, 20, 100, 500, 2000),
      labels = \(x) abs(x),
      trans  = scales::pseudo_log_trans(base = 10, sigma = 1)
    ) +
    scale_x_discrete(labels = \(x) gsub("_", " ", x)) +
    labs(title    = paste0(method, " – Significant DEGs (FDR < ", FDR_CUT, ")"),
         subtitle = "Unhealthy vs Healthy; positive = up in Unhealthy",
         x = NULL, y = "# DEGs") +
    theme_classic(base_size = 11) +
    theme(axis.text.x       = element_text(angle = 40, hjust = 1),
          strip.background   = element_blank(),
          strip.text         = element_text(face = "bold", colour = "grey20"),
          legend.position    = "top",
          panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
          axis.text         = element_text(size = 20),
          axis.title.y      = element_text(size = 20, margin = margin(r = 12)),
          legend.title      = element_text(size = 18),
          legend.text       = element_text(size = 16)
    )
}

# ── Panel B: Top-gene bubble plot ─────────────────────────────────────────────
make_bubble_plot <- function(method) {
  de_all <- map_dfr(TISSUES, function(tis) {
    load_de(method, tis) %>% mutate(tissue = tis)
  }) %>% filter(!is.na(padj), padj < FDR_CUT)

  top_genes <- de_all %>%
    group_by(tissue, cell_type) %>%
    slice_max(abs(logFC), n = TOP_GENES, with_ties = FALSE) %>%
    ungroup()

  gene_order <- top_genes %>%
    arrange(dplyr::desc(abs(logFC))) %>%
    pull("gene") %>%
    unique()

  top_genes <- top_genes %>%
    mutate(gene      = factor(gene, levels = rev(gene_order)),
           neg_log_p = -log10(padj),
           direction = if_else(logFC > 0, "Up", "Down"),
           tissue    = factor(tissue, levels = TISSUES))

  ggplot(top_genes, aes(x = cell_type, y = gene,
                        size = neg_log_p, colour = direction)) +
    geom_point(alpha = 0.8) +
    facet_wrap(~ tissue, scales = "free_x",
               labeller = as_labeller(c(omentum = "Omentum",
                                        subq    = "Subcutaneous"))) +
    scale_colour_manual(values = dir_colours, name = "Direction") +
    scale_size_continuous(name = "−log10(FDR)", range = c(1.5, 6)) +
    labs(title    = paste0(method, " – Top ", TOP_GENES,
                           " DEGs per cell type (FDR < ", FDR_CUT, ")"),
         subtitle = "Size = significance; colour = direction in Unhealthy",
         x = NULL, y = NULL) +
    theme_classic(base_size = 10) +
    theme(axis.text.x     = element_text(angle = 40, hjust = 1),
          axis.text.y     = element_text(size = 7),
          strip.background = element_blank(),
          strip.text       = element_text(face = "bold"),
          legend.position  = "right",
          panel.grid.major = element_line(colour = "grey90", linewidth = 0.25))
}

# ── Panel C: Faceted volcano plots ────────────────────────────────────────────
make_volcano_plot <- function(method) {
  de_all <- map_dfr(TISSUES, function(tis) {
    load_de(method, tis) %>% mutate(tissue = tis)
  }) %>%
    filter(!is.na(padj), !is.na(logFC)) %>%
    mutate(
      neg_log_p = pmin(-log10(padj), 15),   # cap for display
      sig       = padj < FDR_CUT,
      direction = case_when(
        sig & logFC >  0 ~ "Up",
        sig & logFC <  0 ~ "Down",
        TRUE             ~ "NS"
      ),
      tissue = factor(tissue, levels = TISSUES,
                      labels = c("Omentum", "Subcutaneous"))
    )

  # label top genes per panel
  labels <- de_all %>%
    dplyr::filter(sig) %>%
    dplyr::group_by(tissue, cell_type) %>%
    dplyr::slice_max(abs(logFC), n = 3, with_ties = FALSE) %>%
    dplyr::ungroup()

  ggplot(de_all, aes(x = logFC, y = neg_log_p, colour = direction)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_hline(yintercept = -log10(FDR_CUT), linetype = "dashed",
               linewidth = 0.3, colour = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
    ggrepel::geom_text_repel(
      data = labels, aes(label = gene),
      size = 2, max.overlaps = 10, show.legend = FALSE,
      segment.size = 0.2
    ) +
    facet_grid(tissue ~ cell_type) +
    scale_colour_manual(
      values = c(Up = "#B2182B", Down = "#2166AC", NS = "grey70"),
      name = NULL
    ) +
    labs(title    = paste0(method, " – Volcano plots by cell type"),
         subtitle = paste0("Dashed line: FDR = ", FDR_CUT,
                           ";  top 3 genes per panel labelled"),
         x = "log2 fold-change (Unhealthy / Healthy)", y = "−log10(FDR)") +
    theme_classic(base_size = 8) +
    theme(strip.background  = element_blank(),
          strip.text.x      = element_text(face = "bold", size = 7),
          strip.text.y      = element_text(face = "bold", size = 8),
          legend.position   = "top",
          panel.grid.major  = element_line(colour = "grey92", linewidth = 0.2),
          axis.text         = element_text(size = 6))
}

# ── Panel C2: Per-cell-type volcano plots (MAST only) ─────────────────────────
# One ggplot per cell type, faceted by tissue, top n_top genes labelled.
# Returns a named list of ggplots; caller writes a multi-page PDF.
make_volcano_per_celltype_mast <- function(n_top = 10) {
  de_all <- map_dfr(TISSUES, function(tis) {
    load_de("MAST", tis) %>% 
      dplyr::mutate(tissue = tis)
  }) %>%
    dplyr::filter(!is.na(padj), !is.na(logFC)) %>%
    dplyr::mutate(
      neg_log_p = pmin(-log10(padj), 1000),
      sig       = padj < FDR_CUT,
      direction = case_when(
        sig & logFC >  L2FC_CUT ~ "Up",
        sig & logFC <  -L2FC_CUT ~ "Down",
        TRUE             ~ "NS"
      ),
      tissue = factor(tissue, levels = TISSUES,
                      labels = c("Omentum", "Subcutaneous"))
    )

  cell_types <- sort(unique(de_all$cell_type))

  plots <- lapply(cell_types, function(ct) {
    df_ct <- de_all %>% filter(cell_type == ct)

    labels <- df_ct %>%
      dplyr::filter(sig) %>%
      dplyr::group_by(tissue) %>%
      dplyr::slice_max(abs(logFC), n = n_top, with_ties = FALSE) %>%
      dplyr::ungroup()

    n_up   <- sum(df_ct$direction == "Up",   na.rm = TRUE)
    n_down <- sum(df_ct$direction == "Down",  na.rm = TRUE)

    ggplot(df_ct, aes(x = logFC, y = neg_log_p, colour = direction)) +
      geom_point(size = 0.8, alpha = 0.45) +
      geom_hline(yintercept = -log10(FDR_CUT), linetype = "dashed",
                 linewidth  = 0.3, colour = "grey50") +
      geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
      ggrepel::geom_text_repel(
        data          = labels,
        aes(label     = gene),
        size          = 5.8,
        max.overlaps  = 20,
        show.legend   = FALSE,
        segment.size  = 0.2,
        segment.alpha = 0.5
      ) +
      facet_wrap(~tissue, nrow = 1, scales = "fixed") +
      scale_colour_manual(
        values = c(Up = "#B2182B", Down = "#2166AC", NS = "grey70"),
        name   = NULL
      ) +
      labs(
        title    = ct,
        subtitle = sprintf(
          "Top %d sig. genes per tissue labelled  |  FDR < %.2f  |  Up: %d  Down: %d",
          n_top, FDR_CUT, n_up, n_down
        ),
        x = "log2 fold-change (Ob+T2D / Ob)",
        y = "−log10(FDR)"
      ) +
      theme_classic(base_size = 10) +
      theme(
        strip.background  = element_blank(),
        strip.text        = element_text(face = "bold", size = 20),
        legend.position   = "top",
        plot.title        = element_text(face = "bold", size = 12,
                                         colour = ct_colour(ct)),
        plot.subtitle     = element_text(size = 8, colour = "grey40"),
        panel.grid.major  = element_line(colour = "grey93", linewidth = 0.2),
        axis.text         = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1, "cm")
      ) +
     guides(colour = guide_legend(override.aes = list(size = 10)))
  })

  names(plots) <- cell_types
  plots
}

# ── Panel D: Top-DEG logFC heatmap (ggplot) ───────────────────────────────────
make_logfc_heatmap <- function(method) {
  de_all <- map_dfr(TISSUES, function(tis) {
    load_de(method, tis) %>% mutate(tissue = tis)
  }) %>% filter(!is.na(padj), padj < FDR_CUT)

  top_genes <- de_all %>%
    group_by(tissue, cell_type) %>%
    slice_max(abs(logFC), n = TOP_GENES, with_ties = FALSE) %>%
    ungroup() %>%
    pull("gene") %>%
    unique()

  plot_data <- de_all %>%
    filter(gene %in% top_genes) %>%
    mutate(
      logFC_plot = pmax(pmin(logFC, 4), -4),   # clamp for colour scale
      tissue     = factor(tissue, levels = TISSUES,
                          labels = c("Omentum", "Subcutaneous"))
    )

  gene_order <- plot_data %>%
    arrange(dplyr::desc(abs(logFC_plot))) %>%
    pull("gene") %>%
    unique()

  plot_data <- plot_data %>%
    mutate(gene = factor(gene, levels = rev(gene_order)))

  ggplot(plot_data, aes(x = cell_type, y = gene, fill = logFC_plot)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    facet_wrap(~ tissue, scales = "free_x") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "log2FC",
                         limits = c(-4, 4), oob = scales::squish) +
    labs(title    = paste0(method, " – Top DEG logFC heatmap (FDR < ", FDR_CUT, ")"),
         subtitle = paste0("Top ", TOP_GENES, " genes per cell type by |logFC|"),
         x = NULL, y = NULL) +
    theme_classic(base_size = 9) +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y      = element_text(size = 6),
          strip.background  = element_blank(),
          strip.text        = element_text(face = "bold"),
          legend.position   = "right",
          panel.grid        = element_blank())
}

# ── Panels E & F: NES heatmaps via ComplexHeatmap ─────────────────────────────

# Build an NES matrix (pathways × cell_types) for one tissue
.nes_mat <- function(gsea_df, tis, paths) {
  df <- gsea_df %>%
    filter(tissue == tis, pathway %in% paths) %>%
    mutate(NES = ifelse(is.na(NES), 0, NES)) %>%
    dplyr::select(pathway, cell_type, NES)
  wide <- tidyr::pivot_wider(df, names_from = "cell_type", values_from = "NES",
                              values_fill = 0, values_fn = mean)
  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- wide$pathway
  mat
}

# Build a significance-label matrix aligned to a reference NES matrix
.sig_mat <- function(gsea_df, tis, ref_mat) {
  df <- gsea_df %>%
    filter(tissue == tis, pathway %in% rownames(ref_mat)) %>%
    mutate(star = dplyr::case_when(
      !is.na(padj) & padj < 0.01    ~ "**",
      !is.na(padj) & padj < FDR_CUT ~ "*",
      TRUE ~ ""
    )) %>%
    dplyr::select(pathway, cell_type, star)
  wide <- tidyr::pivot_wider(df, names_from = "cell_type", values_from = "star",
                              values_fill = "", values_fn = dplyr::first)
  raw <- as.matrix(wide[, -1, drop = FALSE])
  rownames(raw) <- wide$pathway

  # align rows and cols to ref_mat
  miss_r <- setdiff(rownames(ref_mat), rownames(raw))
  if (length(miss_r))
    raw <- rbind(raw, matrix("", length(miss_r), ncol(raw),
                             dimnames = list(miss_r, colnames(raw))))
  miss_c <- setdiff(colnames(ref_mat), colnames(raw))
  if (length(miss_c))
    raw <- cbind(raw, matrix("", nrow(raw), length(miss_c),
                             dimnames = list(rownames(raw), miss_c)))
  raw[rownames(ref_mat), colnames(ref_mat), drop = FALSE]
}

# Pad a matrix to include all rows in `all_paths` (fill = 0)
.pad_rows <- function(mat, all_paths) {
  miss <- setdiff(all_paths, rownames(mat))
  if (length(miss))
    mat <- rbind(mat, matrix(0, length(miss), ncol(mat),
                             dimnames = list(miss, colnames(mat))))
  mat[all_paths, , drop = FALSE]
}

make_nes_heatmap <- function(method, db, keep_paths = NULL, show_ct_annot = TRUE) {
  db_label <- switch(db, hallmark = "Hallmark", kegg = "KEGG", db)

  gsea_all <- map_dfr(TISSUES, function(tis) {
    load_fgsea(method, tis, db) %>% mutate(tissue = tis)
  })

  if (db == "hallmark")
    gsea_all <- gsea_all %>%
      mutate(pathway = str_remove(pathway, "^HALLMARK_") %>%
                       str_replace_all("_", " ") %>%
                       str_to_title())

  # Select pathways
  if (!is.null(keep_paths)) {
    pattern  <- paste(keep_paths, collapse = "|")
    sig_paths <- gsea_all %>%
      filter(str_detect(pathway, regex(pattern, ignore_case = TRUE))) %>%
      pull("pathway") %>%
      unique()
  } else {
    sig_paths <- gsea_all %>%
      filter(!is.na(padj), padj < FDR_CUT, !str_detect(pathway, viral_re)) %>%
      pull("pathway") %>%
      unique()

    if (length(sig_paths) == 0) {
      message(method, "/", db, ": no sig pathways — showing top ", TOP_PATHS,
              " by p-value")
      sig_paths <- gsea_all %>%
        arrange(pval) %>%
        pull("pathway") %>%
        unique() %>%
        head(TOP_PATHS)
    } else if (length(sig_paths) > TOP_PATHS) {
      sig_paths <- gsea_all %>%
        filter(!is.na(padj), padj < FDR_CUT) %>%
        dplyr::count(pathway, sort = TRUE) %>%
        head(TOP_PATHS) %>%
        pull("pathway")
    }
  }

  mat_om <- .nes_mat(gsea_all, "omentum", sig_paths)
  mat_sq <- .nes_mat(gsea_all, "subq",    sig_paths)
  all_paths <- union(rownames(mat_om), rownames(mat_sq))
  mat_om    <- .pad_rows(mat_om, all_paths)
  mat_sq    <- .pad_rows(mat_sq, all_paths)

  # Cluster rows on combined matrix
  row_clust <- hclust(dist(cbind(mat_om, mat_sq)), method = "ward.D2")

  sig_om <- .sig_mat(gsea_all, "omentum", mat_om)
  sig_sq <- .sig_mat(gsea_all, "subq",    mat_sq)

  # Column annotation helper
  make_col_ann <- function(mat) {
    cts    <- colnames(mat)
    colors <- setNames(ct_colour(cts), cts)
    HeatmapAnnotation(
      `Cell type` = cts,
      col         = list(`Cell type` = colors),
      show_legend = FALSE,
      show_annotation_name = FALSE,
      annotation_height = unit(4, "mm")
    )
  }

  common_args <- list(
    col              = nes_col_fun,
    cluster_rows     = row_clust,
    cluster_columns  = FALSE,
    row_names_gp     = gpar(fontsize = 16),
    column_names_gp  = gpar(fontsize = 18),
    column_names_rot = 45,
    heatmap_legend_param = list(
      title_gp      = gpar(fontsize = 15, fontface = "bold"),
      labels_gp     = gpar(fontsize = 13),
      legend_height = unit(6, "cm")
    )
  )

  ht_om <- do.call(Heatmap, c(list(
    mat_om,
    name            = "NES",
    show_row_dend   = TRUE,
    row_dend_width  = unit(2, "cm"),
    row_names_side  = "left",
    row_labels      = str_wrap(rownames(mat_om), width = 20),
    column_labels   = gsub("_", " ", colnames(mat_om)),
    top_annotation  = if (show_ct_annot) make_col_ann(mat_om) else NULL,
    column_title    = "Omentum",
    column_title_gp = gpar(fontface = "bold"),
    cell_fun = function(j, i, x, y, ...) {
      lbl <- sig_om[i, j]
      if (nzchar(lbl)) grid.text(lbl, x, y, gp = gpar(fontsize = 7))
    }
  ), common_args))

  ht_sq <- do.call(Heatmap, c(list(
    mat_sq,
    name            = "NES ",   # trailing space → separate legend entry avoided
    show_row_names  = FALSE,
    column_labels   = gsub("_", " ", colnames(mat_sq)),
    top_annotation  = if (show_ct_annot) make_col_ann(mat_sq) else NULL,
    column_title    = "Subcutaneous",
    column_title_gp = gpar(fontface = "bold"),
    cell_fun = function(j, i, x, y, ...) {
      lbl <- sig_sq[i, j]
      if (nzchar(lbl)) grid.text(lbl, x, y, gp = gpar(fontsize = 7))
    }
  ), common_args))

  # Legend for cell type colours
  all_cts   <- union(colnames(mat_om), colnames(mat_sq))
  ct_legend <- Legend(
    labels    = all_cts,
    legend_gp = gpar(fill = ct_colour(all_cts)),
    title     = "Cell type",
    nrow      = ceiling(length(all_cts) / 3)
  )

  list(
    ht       = ht_om + ht_sq,
    legends  = list(ct_legend),
    title    = paste0(method, " – ", db_label,
                      " pathway enrichment (NES)\n",
                      "* FDR < 0.05   ** FDR < 0.01   (OB+T2D vs Ob)")
  )
}

save_nes_heatmap <- function(obj, path, width = 14, height = 10) {
  pdf(path, width = width, height = height)
  draw(obj$ht,
       column_title          = obj$title,
       column_title_side     = "top",
       column_title_gp       = gpar(fontsize = 11, fontface = "bold"),
       annotation_legend_list = obj$legends,
       merge_legend          = TRUE,
       heatmap_legend_side   = "right",
       annotation_legend_side = "bottom")
  dev.off()
  invisible(path)
}

pathways_list <- unlist(strsplit("Oxidative phosphorylation
Thermogenesis
Aldosterone synthesis and secretion
Focal adhesion
Cytoskeleton in muscle cells
Regulation of actin cytoskeleton
IgSF CAM signaling
Phagosome
Complement and coagulation cascades
Antigen processing and Presentation
Cell adhesion molecule (CAM) interaction", "\n"))

# ── Render and save ───────────────────────────────────────────────────────────
METHODS <- c("MAST") #"deseq2", "edgeR", 

for (method in METHODS) {
  message("Rendering figures for: ", method)

  fig_A <- make_deg_count_plot(method)
  fig_B <- make_bubble_plot(method)
  fig_C <- make_volcano_plot(method)
  fig_D <- make_logfc_heatmap(method)
  fig_E <- make_nes_heatmap(method, "hallmark")
  fig_F <- make_nes_heatmap(method, "kegg", keep_paths=pathways_list, show_ct_annot=FALSE)

  # ggplot panels saved with ggsave
  ggsave(file.path(OUT_DIR, paste0(method, "_A_deg_counts.pdf")),
         fig_A, width = 5, height = 10, device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0(method, "_B_top_genes.pdf")),
         fig_B, width = 14, height = 9, device = cairo_pdf)
  ggsave(file.path(OUT_DIR, paste0(method, "_C_volcano.pdf")),
         fig_C, width = 22, height = 7, device = cairo_pdf)

  # C2: per-cell-type volcano (MAST only) — one page per cell type
  if (method == "MAST") {
    volcano_ct_list <- make_volcano_per_celltype_mast(n_top = 10)
    out_c2 <- file.path(OUT_DIR, "MAST_C2_volcano_per_celltype.pdf")
    pdf(out_c2, width = 5, height = 6)
    invisible(lapply(volcano_ct_list, print))
    dev.off()
    message("  Saved: MAST_C2_volcano_per_celltype.pdf  (",
            length(volcano_ct_list), " pages)")
  }

  ggsave(file.path(OUT_DIR, paste0(method, "_D_logfc_heatmap.pdf")),
         fig_D, width = 14, height = 10, device = cairo_pdf)

  # ComplexHeatmap panels saved with pdf + draw
  save_nes_heatmap(fig_E,
    file.path(OUT_DIR, paste0(method, "_E_hallmark.pdf")),
    width = 14, height = 10)
  save_nes_heatmap(fig_F,
    file.path(OUT_DIR, paste0(method, "_F_kegg.pdf")),
    width = 16, height = 8)

  # ggplot-only combined (A + B top row, C + D bottom rows)
  combined_gg <- (fig_A | fig_B) / fig_C / fig_D +
    plot_annotation(
      title = paste0("Overview: Healthy vs Unhealthy (", method, ")"),
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    ) +
    plot_layout(heights = c(1.2, 1.0, 1.2))

  ggsave(file.path(OUT_DIR, paste0(method, "_overview_ggplot.pdf")),
         combined_gg, width = 22, height = 26, device = cairo_pdf)

  message("  Saved: ", OUT_DIR, "/", method, "_*.pdf")
}

message("Done.")
