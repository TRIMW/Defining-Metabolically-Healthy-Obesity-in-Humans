library(tidyverse)
library(ggplot2)
library(ggdendro)
library(aplot)
library(tidytable)

# ── Configuration ─────────────────────────────────────────────────────────────
TISSUES <- c(
  Omentum = "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq/results/DE/omentum/fgsea/cell_type2_/MAST",
  Subq    = "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq/results/DE/subq/fgsea/cell_type2_/MAST"
)
OUT_DIR     <- "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq/results/DE"
TOP_N       <- 15
PADJ_CUTOFF <- 0.3
CT_GAP      <- 0.7   # extra space (in column-width units) between cell-type groups

TISSUE_COLORS <- c(Omentum = "#dba065", Subq = "#65db82")

viral_re <- regex(
  "(sars.cov|Infection|Influenza|Tuberculosis|type I diabetes|type 1 diabetes|autoimmune thyroid disease|toxoplasmosis|covid|coronavirus|influenza infection|Parkinson Disease|Amyotrophic lateral sclerosis|viral myocarditis|Gastric acid secretion|cardiomyopathy)",
  ignore_case = TRUE
)

GO_LEVEL_MIN <- 3

go_bp_allowed <- read_csv(
  "/Users/willtrim/Documents/reference/GO/go_bp_id_to_name.csv",
  show_col_types = FALSE
) |>
  filter(obsolete == "FALSE", level >= GO_LEVEL_MIN) |>
  pull(msigdb_name)   # raw GOBP_XXX names, matching the CSV pathway column before transformation

# Cell types shared across tissues (tissue-specific ones will appear as a single column)
lymphoid_cts <- c(
  "B_cells", "Plasma_cells", "mNK1_cells",           "mNK2_cell",            "mNK3_cells",        "GD_T_cells",          "CD8_TRM",              "Effector_CD8_T_cells",  "CD8_T_cells",          "EM_CD4_T_cells",      
  "NA_CD4_T_cells",          "CD4_T_cells",         "CM_CD4_T_cells",       "Tregs", "PLGC2_cells"
)
myeloid_cts <- c( "cDCs",  "cDC1s", "cDC2Bs",  "pDCs", "Macrophages",  "Monocytes",   "Neutrophils")

others_cts <- c("LECs",           "Myofibroblasts", "VECs",           "FAPs",           "Pre_adipocytes")


cts_list <- list(
  lymphoid = lymphoid_cts,
  myeloid  = myeloid_cts,
  others   = others_cts
)

# ── Helpers ───────────────────────────────────────────────────────────────────
load_hallmark <- function(ct, dir, tissue, db = "hallmark", keep = NULL) {
  f <- file.path(dir, paste0(ct, ".", db, ".csv"))
  if (!file.exists(f)) {
    # print(paste0("file does not exists: ", f))
    return(NULL)
  }
  # else {
    # print(paste0("file exists: ", f))
  # }
  df <- read_csv(f, show_col_types = FALSE) |>
    mutate(across(c(pval, padj, log2err, ES, NES, size), as.numeric))
  if (db == "bp") df <- filter(df, pathway %in% go_bp_allowed)
  df <- df |>
    mutate(
      cell_type = ct,
      tissue    = tissue,
      pathway   = str_remove(pathway, "^HALLMARK_") |>
                  str_remove("^GOBP_")              |>
                  str_remove("^REACTOME_")              |>
                  str_replace_all("_", " ")         |>
                  str_to_title()
    )
  
  if (!is.null(keep)) df <- filter(df, pathway %in% keep)
  return(df)
}

make_dotplot <- function(dat, top_paths, col_info, tissue_colors, padj_cutoff, title,
                         clim = NULL, nes_max = NULL, pathway_order = NULL) {

  plot_dat <- dat |>
    filter(pathway %in% top_paths) |>
    mutate(
      color_val = sign(NES) * (-log10(padj + 1e-300)),
      sig       = padj < padj_cutoff
    ) |>
    left_join(col_info |> select(cell_type, tissue, x_pos), by = c("cell_type", "tissue"))

  if (nrow(plot_dat) == 0) return(NULL)

  if (!is.null(pathway_order)) {
    # ── Use caller-supplied order; retain only pathways present in the data ────
    row_order <- pathway_order[pathway_order %in% plot_dat$pathway]
    row_hc    <- NULL
  } else {
    # ── Cluster pathways on the full NES matrix ────────────────────────────────
    mat_nes <- plot_dat |>
      select(pathway, x_pos, NES) |>
      pivot_wider(names_from = x_pos, values_from = NES, values_fill = 0) |>
      column_to_rownames("pathway") |>
      as.matrix()

    row_hc    <- hclust(dist(mat_nes), method = "ward.D2")
    row_order <- rownames(mat_nes)[row_hc$order]
  }
  n_row <- length(row_order)

  plot_dat <- plot_dat |> mutate(y_int = match(pathway, row_order))

  # ── Shared scales ────────────────────────────────────────────────────────────
  x_lim <- c(min(col_info$x_pos) - 1.25, max(col_info$x_pos) + 0.5)
  x_scale <- scale_x_continuous(limits = x_lim, expand = c(0, 0))
  y_scale <- scale_y_continuous(
    breaks = seq_len(n_row), labels = str_wrap(row_order, width = 50),
    limits = c(0.5, n_row + 0.5), expand = c(0, 0)
  )

  if (is.null(clim))    clim    <- max(abs(plot_dat$color_val),   na.rm = TRUE)
  if (is.null(nes_max)) nes_max <- max(abs(plot_dat$NES),         na.rm = TRUE)

  # Dashed separators at midpoints between cell-type groups
  ct_sep <- col_info |>
    summarize(sep = max(x_pos) + 0.5 + CT_GAP / 2, .by = ct_idx) |>
    filter(ct_idx < max(ct_idx)) |>
    pull(sep)

  # ── Main dotplot ─────────────────────────────────────────────────────────────
  main_plot <- ggplot(plot_dat, aes(x = x_pos, y = y_int)) +
    geom_vline(xintercept = ct_sep, color = "grey70", linewidth = 0.4, linetype = "dashed") +
    geom_point(aes(size = abs(NES), color = color_val, shape = ifelse(sig, 16, 1))) +
    scale_shape_identity() +
    scale_size_continuous(
      name   = "|NES|",
      range  = c(1, 8),
      limits = c(0, nes_max)
    ) +
    scale_color_gradient2(
      name     = expression(atop("sign(NES) ×", -log[10](padj))),
      low      = "#2166AC", mid = "white", high = "#D6604D",
      midpoint = 0, limits = c(-clim, clim),
      trans    = scales::pseudo_log_trans(sigma = 0.5),
      guide    = guide_colorbar(barheight = unit(5, "cm"))
    ) +
    x_scale + y_scale +
    labs(x = NULL, y = NULL, title = title) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.text.y        = element_text(size = 18),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      legend.position    = "right",
      legend.title       = element_text(size = 19),
      legend.text        = element_text(size = 16),
      plot.title         = element_text(size = 14, face = "bold")
    )

  # ── Tissue annotation strip (top) ────────────────────────────────────────────
  tissue_annot <- ggplot(col_info, aes(x = x_pos, y = 1, fill = tissue)) +
    geom_tile(width = 1, height = 1, color = "white", linewidth = 0.5) +
    scale_fill_manual(values = tissue_colors, name = "Tissue") +
    x_scale +
    theme_void() +
    theme(
      legend.position  = "right",
      legend.key.size  = unit(0.5, "cm"),
      legend.text      = element_text(size = 18),
      legend.title     = element_text(size = 19)
    )

  # ── Cell-type label strip (bottom) ───────────────────────────────────────────
  ct_label_dat <- col_info |>
    summarize(x_center = mean(x_pos), .by = c(cell_type, ct_idx)) |>
    mutate(label = gsub("_", " ", cell_type))

  ct_label_strip <- ggplot(ct_label_dat, aes(x = x_center, y = 0, label = label)) +
    geom_text(angle = 40, hjust = 1, vjust = 1, size = 6.0) +
    x_scale +
    expand_limits(y = c(-2, 0)) +
    theme_void()

  # ── Assemble ─────────────────────────────────────────────────────────────────
  fig <- main_plot |>
    insert_top(tissue_annot,      height = 0.04) |>
    insert_bottom(ct_label_strip, height = 0.44)

  if (!is.null(row_hc)) {
    row_ddata  <- dendro_data(as.dendrogram(row_hc), type = "rectangle")
    row_dendro <- ggplot(segment(row_ddata)) +
      geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend)) +
      y_scale +
      theme_void()
    fig <- fig |> insert_left(row_dendro, width = 0.10)
  }

  fig
}

# ── Main loop ─────────────────────────────────────────────────────────────────
for (db_name in c("hallmark", "bp", "kegg", "reactome")) {
  for (cts_name in names(cts_list)) {

    CELL_TYPES <- cts_list[[cts_name]]

    # Load from both tissues; NULL entries (missing files) are dropped by map_dfr
    dat <- imap_dfr(TISSUES, \(dir, tissue)
      map_dfr(CELL_TYPES, load_hallmark, dir = dir, tissue = tissue, db = db_name)
    )

    if (nrow(dat) == 0) next

    # Top pathways: select per (cell_type, tissue), then union
    top_paths <- dat |>
      filter(padj < PADJ_CUTOFF, !str_detect(pathway, viral_re)) |>
      group_by(cell_type, tissue) |>
      dplyr::slice_max(abs(NES), n = TOP_N, with_ties = FALSE) |>
      ungroup() |>
      pull(pathway) |>
      unique()

    if (length(top_paths) == 0) next

    # Column layout: paired by cell type (omentum | subq), gap between groups.
    # semi_join drops combos where no data was loaded (tissue-specific cell types
    # end up with a single column).
    col_info <- crossing(cell_type = CELL_TYPES, tissue = names(TISSUES)) |>
      semi_join(dat, by = c("cell_type", "tissue")) |>
      mutate(
        ct_idx = match(cell_type, CELL_TYPES),
        ti_idx = match(tissue, names(TISSUES))
      ) |>
      arrange(ct_idx, ti_idx) |>
      mutate(x_pos = row_number() + (ct_idx - 1) * CT_GAP)

    fig <- make_dotplot(
      dat, top_paths, col_info, TISSUE_COLORS, PADJ_CUTOFF,
      title = paste0(db_name, " — ", cts_name, " cell types")
    )
    if (is.null(fig)) next

    n_cols <- nrow(col_info)
    n_rows <- length(unique(top_paths))
    out <- file.path(OUT_DIR,
                     paste0(db_name, "_dotplot_", cts_name, "_two_tissues_top", TOP_N, ".pdf"))
    pdf(out,
        width  = max(8,  n_cols * 0.85 + 5),
        height = max(5,  n_rows * 0.55 + 3))
    print(fig)
    dev.off()
    message("Saved: ", out)
  }
}


# ── Curated multi-DB dotplots (two tissues) ───────────────────────────────────
# Hand-picked pathways from GOBP / HALLMARK / KEGG combined into one dotplot
# per lineage, using the same paired omentum | subq column layout above.
# Pathway labels get a [H/G/K] suffix only when the same name appears from
# more than one database.

curated_paths <- list(
  myeloid = list(
    bp = c(
      "Antigen Processing And Presentation Of Peptide Or Polysaccharide Antigen Via Mhc Class Ii",
      "Antigen Processing And Presentation"
    ),
    hallmark = c(
      "Interferon Alpha Response",
      "Hypoxia",
      "Tnfa Signaling Via Nfkb",
      "Reactive Oxygen Species Pathway",
      "Interferon Gamma Response"
    ),
    kegg = c(
      "Cell Adhesion Molecule (Cam) Interaction"
    )
  ),
  lymphoid = list(
    hallmark = c(
      "Tnfa Signaling Via Nfkb",
      "Tgf Beta Signaling"
    ),
    bp = c(
      "Leukocyte Chemotaxis",
      "Leukocyte Migration",
      "Leukocyte Mediated Cytotoxicity",
      "Antigen Processing And Presentation"
    ),
    kegg = c(
      "Oxidative Phosphorylation",
      "Chemical Carcinogenesis − Reactive Oxygen Species"
    )
  )
)

# curated_paths <- list(
#   others = list(
#     bp = c(
#       "GOBP_RESPONSE_TO_LIPID",
#       "GOBP_EPIGENETIC_REGULATION_OF_GENE_EXPRESSION",
#       "GOBP_AEROBIC_RESPIRATION",
#       "GOBP_CELL_ADHESION",
#       "GOBP_CELL_CELL_ADHESION",
#       "GOBP_PHOSPHORYLATION",
#       "GOBP_FAT_CELL_DIFFERENTIATION",
#       "GOBP_POSITIVE_REGULATION_OF_CELL_ADHESION",
#       "GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_CELL_CELL_ADHESION",
#       "GOBP_LEUKOCYTE_CELL_CELL_ADHESION",
#       "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
#       "GOBP_REGULATION_OF_T_CELL_ACTIVATION",
#       "GOBP_REGULATION_OF_CELL_DEVELOPMENT",
#       "GOBP_RESPONSE_TO_GROWTH_FACTOR",
#       "GOBP_TYPE_II_INTERFERON_PRODUCTION",
#       "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
#       "GOBP_COLLAGEN_FIBRIL_ORGANIZATION"
#     ),
#     hallmark = strsplit(
#       "HALLMARK_TNFA_SIGNALING_VIA_NFKB
# HALLMARK_INTERFERON_GAMMA_RESPONSE
# HALLMARK_INFLAMMATORY_RESPONSE
# HALLMARK_INTERFERON_ALPHA_RESPONSE
# HALLMARK_CHOLESTEROL_HOMEOSTASIS
# HALLMARK_TNFA_SIGNALING_VIA_NFKB
# HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
# HALLMARK_INFLAMMATORY_RESPONSE
# HALLMARK_TGF_BETA_SIGNALING
# HALLMARK_INTERFERON_GAMMA_RESPONSE
# HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "\n"),
#     kegg = strsplit(
#       "IL-17 signaling pathway
# PI3K-Akt signaling pathway
# Cadherin signaling
# Antigen processing and presentation
# Focal adhesion
# ECM-receptor interaction
# Integrin signaling
# Insulin resistance
# Leukocyte transendothelial migration
# Antigen processing and presentation
# Th17 cell differentiation
# Cell adhesion molecule (CAM) interaction
# Cadherin signaling
# Lipid and atherosclerosis
# Vascular smooth muscle contraction
# PPAR signaling pathway
# Cell adhesion molecule (CAM) interaction
# cGMP-PKG signaling pathway
# Thermogenesis", "\n"),
#     reactome = c(
#       "REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING"
#     )
#   )
# )

# pathways_order <- unlist(strsplit("Aerobic Respiration
# Thermogenesis
# Ppar Signaling Pathway
# Lipid And Atherosclerosis
# Response To Lipid
# Tnfa Signaling Via Nfkb
# Il−17 Signaling Pathway
# Interferon Alpha Response
# Interferon Gamma Response
# Antigen Processing And Presentation
# Inflammatory Response
# Tgf Beta Signaling
# Vascular Smooth Muscle Contraction
# Fat Cell Differentiation
# Epithelial Mesenchymal Transition
# Cell Cell Adhesion
# Positive Regulation Of Cell Adhesion
# Leukocyte Cell Cell Adhesion
# Leukocyte Transendothelial Migration", "\n"))

db_abbrev <- c(hallmark = "H", bp = "G", kegg = "K", reactome = "R")

# ── Global color / size scales (computed once across both lineages) ────────────
all_curated_raw <- map_dfr(names(curated_paths), \(cts_name) {
  keep_paths <- curated_paths[[cts_name]]
  imap_dfr(TISSUES, \(dir, tissue)
    map_dfr(names(keep_paths), \(db)
      map_dfr(cts_list[[cts_name]], load_hallmark, dir = dir, tissue = tissue, db = db, keep = unlist(keep_paths[[db]]))
    )
  )
})

global_nes_max <- max(abs(all_curated_raw$NES), na.rm = TRUE)
global_clim    <- max(
  abs(sign(all_curated_raw$NES) * (-log10(all_curated_raw$padj + 1e-300))),
  na.rm = TRUE
)

for (cts_name in c("myeloid", "lymphoid")) {
  CELL_TYPES <- cts_list[[cts_name]]
  keep_paths <- curated_paths[[cts_name]]

  # Load curated pathways from both tissues across all specified databases
  combined_dat <- imap_dfr(TISSUES, \(dir, tissue)
    map_dfr(names(keep_paths), \(db)
      map_dfr(CELL_TYPES, load_hallmark, dir = dir, tissue = tissue, db = db, keep = unlist(keep_paths[[db]])) |>
        mutate(db = db)
    )
  )
  
  # combined_dat <- combined_dat %>%
  #   filter(pathway %in% pathways_order)

  if (nrow(combined_dat) == 0) { warning("No curated data for: ", cts_name); next }

  # Append [H/G/K] only for names that appear from more than one database
  path_db_n <- combined_dat |>
    distinct(pathway, db) |>
    count(pathway, name = "n_db")

  combined_dat <- combined_dat |>
    left_join(path_db_n, by = "pathway") |>
    mutate(pathway = ifelse(
      n_db > 1,
      paste0(pathway, " [", db_abbrev[db], "]"),
      pathway
    ))

  top_paths <- unique(combined_dat$pathway)

  # Column layout: same paired (omentum | subq) structure as the main loop
  col_info <- crossing(cell_type = CELL_TYPES, tissue = names(TISSUES)) |>
    semi_join(combined_dat, by = c("cell_type", "tissue")) |>
    dplyr::mutate(
      ct_idx = match(cell_type, CELL_TYPES),
      ti_idx = match(tissue, names(TISSUES))
    ) |>
    dplyr::arrange(ct_idx, ti_idx) |>
    dplyr::mutate(x_pos = dplyr::row_number() + (ct_idx - 1) * CT_GAP)

  fig <- make_dotplot(
    combined_dat, top_paths, col_info, TISSUE_COLORS, PADJ_CUTOFF,
    title   = NULL, #paste0("Curated multi-DB — ", cts_name, " cell types"),
    clim    = global_clim,
    nes_max = global_nes_max,
    pathway_order = NULL#pathways_order
  )
  if (is.null(fig)) next

  n_cols <- nrow(col_info)
  n_rows <- length(top_paths)
  out <- file.path(OUT_DIR, paste0("curated_multidb_dotplot_", cts_name, "_two_tissues.pdf"))
  pdf(out,
      width  = max(5,  n_cols * 0.4 + 8),
      height = max(5,  n_rows * 0.5 + 3))
  print(fig)
  dev.off()
  message("Saved: ", out)
}

