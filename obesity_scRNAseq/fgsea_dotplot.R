library(tidyverse)
library(ggplot2)
library(ggdendro)
library(aplot)


load_hallmark <- function(ct, dir, db = "hallmark") {
  f <- file.path(dir, paste0(ct, ".", db, ".csv"))
  if (!file.exists(f)) { warning("Missing: ", f); return(NULL) }
  read_csv(f, show_col_types = FALSE) |>
    mutate(cell_type = ct,
           pathway = str_remove(pathway, "^HALLMARK_") %>%
            str_remove("^GOBP_") %>%
             str_replace_all("_", " ") %>%
             str_to_title())
}

lymphoid_cts <- c(
  "mNK2_cell",            "mNK1_cells",           "Effector_CD8_T_cells", "mNK3_cells",           "CD8_TRM",              "CD4_T_cells",         
  "CM_CD4_T_cells",       "GD_T_cells",           "CD8_T_cells",          "Mast_cells_",          "B_cells",              "EM_CD4_T_cells",      
  "NA_CD4_T_cells",       "Plasma_cells", "Tregs", "PLGC2_cells"
)
myeloid_cts <- c( "Monocytes",   "Neutrophils", "Macrophages", "cDCs",        "pDCs"  )

others_cts <- c("LECs",           "Myofibroblasts", "VECs",           "FAPs",           "Pre_adipocytes")


cts_list <- list(
  "others" = others_cts,
  "myeloid" = myeloid_cts,
  "lymphoid" = lymphoid_cts
)

cts_name <- "lymphoid"

# ── Configuration ─────────────────────────────────────────────────────────────
DATA_DIR    <- "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq/results/DE/omentum/fgsea/cell_type2_/MAST"
TOP_N       <- 10        # top N enriched pathways per cell type (by |NES|, padj < 0.05)
PADJ_CUTOFF <- 0.1
viral_re <- regex(
  "(sars.cov|type 1 diabetes|autoimmune thyroid disease|toxoplasmosis|covid|coronavirus|influenza infection|Parkinson Disease|Amyotrohic lateral sclerosis|viral myocarditis|Gastric acid secretion|cardiomyopathy)", ignore_case = TRUE)


for(db_name in c("hallmark", "bp", "kegg", "reactome")) {
  for (cts_name in names(cts_list)) {
  
    OUT_FILE    <- file.path(DATA_DIR, paste0(db_name, "_dotplot_", cts_name, "_top", TOP_N, ".pdf"))
    
    CELL_TYPES  <- cts_list[[cts_name]]

    # ── Load & tidy ───────────────────────────────────────────────────────────────
    dat <- map_dfr(CELL_TYPES, load_hallmark, dir = DATA_DIR, db = db_name)
    
    # ── Select top pathways per cell type ────────────────────────────────────────
    top_paths <- dat |>
      filter(padj < PADJ_CUTOFF, !str_detect(pathway, viral_re)) |>
      group_by(cell_type) |>
      slice_max(order_by = abs(NES), n = TOP_N, with_ties = FALSE) |>
      ungroup() |>
      pull(pathway) |>
      unique()
    
    plot_dat <- dat |>
      filter(pathway %in% top_paths) |>
      mutate(
        color_val = sign(NES) * (-log10(padj + 1e-300)),
        sig       = padj < PADJ_CUTOFF
      )
    
    # ── Cluster rows (pathways) and columns (cell types) ─────────────────────────
    mat_nes <- plot_dat |>
      dplyr::select(pathway, cell_type, NES) |>
      pivot_wider(names_from = cell_type, values_from = NES, values_fill = 0) |>
      column_to_rownames("pathway") |>
      as.matrix()
    
    # keep only columns for selected cell types that are present
    mat_nes <- mat_nes[, intersect(CELL_TYPES, colnames(mat_nes)), drop = FALSE]
    
    row_hc  <- hclust(dist(mat_nes),      method = "ward.D2")
    col_hc  <- hclust(dist(t(mat_nes)),   method = "ward.D2")
    
    row_order <- rownames(mat_nes)[row_hc$order]
    col_order <- colnames(mat_nes)[col_hc$order]
    
    n_col <- length(col_order)
    n_row <- length(row_order)
    
    # Integer positions for x/y — same coordinate space for main plot and dendrograms.
    # Limits mirror ggplot2's discrete-axis default (expand by 0.5 on each side).
    x_scale <- scale_x_continuous(
      breaks = seq_len(n_col), labels = col_order,
      limits = c(0.5, n_col + 0.5), expand = c(0, 0)
    )
    y_scale <- scale_y_continuous(
      breaks = seq_len(n_row), labels = row_order,
      limits = c(0.5, n_row + 0.5), expand = c(0, 0)
    )
    
    plot_dat <- plot_dat |>
      mutate(
        x_int = match(cell_type, col_order),
        y_int = match(pathway,   row_order)
      )
    
    # ── Build dendrograms ─────────────────────────────────────────────────────────
    # Reuse the same scale objects so all three plots share identical numeric limits.
    # aplot::insert_top / insert_left then aligns panels without a scale-type clash.
    
    col_ddata  <- dendro_data(as.dendrogram(col_hc), type = "rectangle")
    col_dendro <- ggplot(segment(col_ddata)) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      x_scale +
      theme_void()
    
    # Swap x<->y in aes (instead of coord_flip) so the shared side is the y-axis.
    # Negate the depth axis so the root faces left.
    row_ddata  <- dendro_data(as.dendrogram(row_hc), type = "rectangle")
    row_dendro <- ggplot(segment(row_ddata)) +
      geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend)) +
      y_scale +
      theme_void()
    
    # ── Color scale limits (symmetric around 0) ───────────────────────────────────
    clim <- max(abs(plot_dat$color_val), na.rm = TRUE)
    
    # ── Main dotplot ──────────────────────────────────────────────────────────────
    main_plot <- ggplot(plot_dat, aes(x = x_int, y = y_int)) +
      geom_point(aes(
        size  = abs(NES),
        color = color_val,
        shape = ifelse(sig, 16, 1)
      )) +
      scale_shape_identity() +
      scale_size_continuous(
        name   = "|NES|",
        range  = c(1, 8),
        limits = c(0, max(abs(plot_dat$NES), na.rm = TRUE))
      ) +
      scale_color_gradient2(
        name     = "sign(NES) ×\n–log₁₀(padj)",
        low      = "#2166AC",
        mid      = "white",
        high     = "#D6604D",
        midpoint = 0,
        limits   = c(-clim, clim)
      ) +
      x_scale + y_scale +
      labs(x = NULL, y = NULL) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y       = element_text(size = 8),
        panel.grid.major  = element_line(color = "grey90", linewidth = 0.3),
        panel.grid.minor  = element_blank(),
        legend.position   = "right"
      )
    
    # ── Assemble: aplot aligns adjacent plots by their shared axis labels ─────────
    final <- main_plot |>
      insert_top(col_dendro,  height = 0.12) |>
      insert_left(row_dendro, width  = 0.15)
    
    # ── Save ─────────────────────────────────────────────────────────────────────
    n_rows <- length(row_order)
    n_cols <- length(col_order)
    pdf(OUT_FILE,
        width  = max(6,  n_cols * 1.0 + 4),
        height = max(5,  n_rows * 0.35 + 3))
    print(final)
    dev.off()
    message("Saved: ", OUT_FILE)
  }
}
