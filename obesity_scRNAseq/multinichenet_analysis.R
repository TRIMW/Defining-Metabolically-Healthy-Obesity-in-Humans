# library(devtools)
# devtools::install_github("saeyslab/multinichenetr")

suppressPackageStartupMessages({
  library(multinichenetr)
  library(nichenetr)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
})

# Local override of multinichenetr::make_sample_lr_prod_activity_plots.
# Extra args:
#   row_spacing           numeric in "lines" units; gap *between* sender-receiver
#                         blocks (default 0.05 — package default was 0.25)
#   inner_row_spacing     numeric in discrete-axis units; padding above/below each
#                         row *within* a block (default 0.4 — ggplot2 default is
#                         0.6; use 0 to pack rows edge-to-edge)
#   sender_receiver_order optional character vector of "Sender --> Receiver"
#                         strings; controls the order of the row-facet blocks.
#                         Any values not present in the data are silently dropped.
#   group_labels          optional named character vector mapping original group
#                         names to display labels, e.g.
#                         c("Healthy" = "Ob", "Unhealthy" = "Ob+T2D").
#                         Also renames sample columns whose prefix matches a key
#                         (e.g. "Healthy_0" -> "Ob_0").
#   expr_limits           optional length-2 numeric vector fixing the color scale
#                         limits for the expression dotplot (p1), e.g. c(-3, 3).
#                         Defaults to the data range of the current figure.
#   cs_limits             optional length-2 numeric vector fixing the color scale
#                         limits for the celltype-specificity dotplot (p_cs),
#                         e.g. c(0, 1). Defaults to the data range.
make_sample_lr_prod_activity_plots <- function(
    prioritization_tables, prioritized_tbl_oi,
    widths = NULL,
    row_spacing = 0.05,
    inner_row_spacing = 0.4,
    sender_receiver_order = NULL,
    group_labels = NULL,
    expr_limits = NULL,
    cs_limits = NULL) {

  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  sample_data <- prioritization_tables$sample_prioritization_tbl %>%
    dplyr::filter(id %in% prioritized_tbl_oi$id) %>%
    dplyr::mutate(
      sender_receiver = paste(sender, receiver, sep = " --> "),
      lr_interaction  = paste(ligand, receptor, sep = " - ")
    ) %>%
    dplyr::arrange(receiver) %>%
    dplyr::group_by(receiver) %>%
    dplyr::arrange(sender, .by_group = TRUE)

  sr_levels <- if (!is.null(sender_receiver_order)) {
    intersect(sender_receiver_order, unique(sample_data$sender_receiver))
  } else {
    unique(sample_data$sender_receiver)
  }
  sample_data <- sample_data %>%
    dplyr::filter(sender_receiver %in% sr_levels) %>%
    dplyr::mutate(sender_receiver = factor(sender_receiver, levels = sr_levels))

  group_data <- prioritization_tables$group_prioritization_table_source %>%
    dplyr::mutate(
      sender_receiver = paste(sender, receiver, sep = " --> "),
      lr_interaction  = paste(ligand, receptor, sep = " - ")
    ) %>%
    dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction,
                    group, activity, activity_scaled, direction_regulation,
                    prioritization_score) %>%
    dplyr::filter(id %in% sample_data$id) %>%
    dplyr::arrange(receiver) %>%
    dplyr::group_by(receiver) %>%
    dplyr::arrange(sender, .by_group = TRUE) %>%
    dplyr::mutate(sender_receiver = factor(sender_receiver, levels = sr_levels))

  group_data_celltype_specificity <- prioritization_tables$group_prioritization_tbl %>%
    dplyr::mutate(
      sender_receiver = paste(sender, receiver, sep = " --> "),
      lr_interaction  = paste(ligand, receptor, sep = " - ")
    ) %>%
    dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction,
                    group, scaled_pb_ligand, scaled_pb_receptor) %>%
    dplyr::filter(id %in% sample_data$id) %>%
    dplyr::arrange(receiver) %>%
    dplyr::group_by(receiver) %>%
    dplyr::arrange(sender, .by_group = TRUE) %>%
    dplyr::mutate(sender_receiver = factor(sender_receiver, levels = sr_levels))

  group_data_frac_expression <- prioritization_tables$group_prioritization_table_source %>%
    dplyr::mutate(
      sender_receiver = paste(sender, receiver, sep = " --> "),
      lr_interaction  = paste(ligand, receptor, sep = " - ")
    ) %>%
    dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction,
                    group, fraction_ligand_group, fraction_receptor_group) %>%
    dplyr::filter(id %in% sample_data$id) %>%
    dplyr::arrange(receiver) %>%
    dplyr::group_by(receiver) %>%
    dplyr::arrange(sender, .by_group = TRUE) %>%
    dplyr::mutate(sender_receiver = factor(sender_receiver, levels = sr_levels))

  group_data <- group_data %>%
    inner_join(group_data_celltype_specificity) %>%
    inner_join(group_data_frac_expression) %>%
    dplyr::mutate(sender_receiver = factor(sender_receiver, levels = sr_levels))

  rm(group_data_celltype_specificity, group_data_frac_expression)

  if (!is.null(group_labels)) {
    recode_group <- function(x) {
      for (old in names(group_labels)) x <- gsub(paste0("^", old, "$"), group_labels[old], x)
      x
    }
    recode_sample <- function(x) {
      for (old in names(group_labels)) x <- gsub(paste0("^", old, "_"), paste0(group_labels[old], "_"), x)
      x
    }
    sample_data <- sample_data %>%
      dplyr::mutate(group  = recode_group(group),
                    sample = recode_sample(sample))
    group_data <- group_data %>%
      dplyr::mutate(group = recode_group(group))
  }

  keep_sender_receiver_values        <- c(0.25, 0.9, 1.75, 4)
  names(keep_sender_receiver_values) <- levels(sample_data$keep_sender_receiver)

  sp_y    <- unit(row_spacing, "lines")
  y_expnd <- expansion(add = inner_row_spacing)

  p1 <- sample_data %>%
    ggplot(aes(sample, lr_interaction,
               color = scaled_LR_pb_prod, size = keep_sender_receiver)) +
    geom_point() +
    facet_grid(sender_receiver ~ group, scales = "free", space = "free",
               switch = "y") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(expand = y_expnd) +
    theme_light() +
    theme(
      axis.ticks = element_blank(), axis.title = element_blank(),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9, angle = 90, hjust = 0),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.spacing.x  = unit(0.4, "lines"), panel.spacing.y = sp_y,
      strip.text.x.top  = element_text(size = 10, color = "black",
                                        face = "bold", angle = 0),
      strip.text.y.left = element_text(size = 9,  color = "black",
                                        face = "bold", angle = 0),
      strip.background  = element_rect(color = "darkgrey", fill = "whitesmoke",
                                        size = 1.5, linetype = "solid")
    ) +
    labs(color = "Scaled L-R\npseudobulk exprs product",
         size  = "Sufficient presence\nof sender & receiver") +
    scale_size_manual(values = keep_sender_receiver_values)

  max_lfc         <- abs(sample_data$scaled_LR_pb_prod) %>% max()
  custom_scale_p1 <- scale_color_gradientn(
    colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),
    values  = c(0, 0.35, 0.485, 0.5, 0.515, 0.65, 1),
    limits  = if (!is.null(expr_limits)) expr_limits
              else c(-1 * max_lfc, max_lfc)
  )
  p1 <- p1 + custom_scale_p1

  p2 <- group_data %>%
    ggplot(aes(direction_regulation, lr_interaction, fill = activity_scaled)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(sender_receiver ~ group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(expand = y_expnd) +
    theme_light() +
    theme(
      axis.ticks = element_blank(), axis.title = element_blank(),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9, angle = 90, hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.spacing.x  = unit(0.2, "lines"), panel.spacing.y = sp_y,
      strip.text.x     = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y     = element_blank(),
      strip.background = element_rect(color = "darkgrey", fill = "whitesmoke",
                                       size = 1.5, linetype = "solid")
    ) +
    labs(fill = "Scaled Ligand\nActivity in Receiver")

  max_activity    <- abs(group_data$activity_scaled) %>% max(na.rm = TRUE)
  custom_scale_p2 <- scale_fill_gradientn(
    colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]),
    values  = c(0, 0.51, 0.575, 0.625, 0.675, 0.725, 1),
    limits  = c(-1 * max_activity, max_activity)
  )
  p2 <- p2 + custom_scale_p2

  p3 <- group_data %>%
    ggplot(aes(direction_regulation, lr_interaction, fill = activity)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(sender_receiver ~ group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(expand = y_expnd) +
    theme_light() +
    theme(
      axis.ticks = element_blank(), axis.title = element_blank(),
      axis.title.y = element_blank(), axis.text.y = element_blank(),
      axis.text.x  = element_text(size = 9, angle = 90, hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.spacing.x  = unit(0.2, "lines"), panel.spacing.y = sp_y,
      strip.text.x     = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y     = element_blank(),
      strip.background = element_rect(color = "darkgrey", fill = "whitesmoke",
                                       size = 1.5, linetype = "solid")
    ) +
    labs(fill = "Ligand\nActivity in Receiver")
  p3 <- p3 + scale_fill_gradient2(low = "white", mid = "white",
                                   high = "darkorange", midpoint = 0)

  cs_data <- group_data %>%
    distinct(sender_receiver, lr_interaction, group,
             scaled_pb_ligand, scaled_pb_receptor) %>%
    tidyr::gather(LR, celltype_specificity, scaled_pb_ligand:scaled_pb_receptor)
  cs_data$LR[cs_data$LR == "scaled_pb_ligand"]   <- "ligand"
  cs_data$LR[cs_data$LR == "scaled_pb_receptor"] <- "receptor"

  frac_data <- group_data %>%
    distinct(sender_receiver, lr_interaction, group,
             fraction_ligand_group, fraction_receptor_group) %>%
    tidyr::gather(LR, fraction_expression,
                  fraction_ligand_group:fraction_receptor_group)
  frac_data$LR[frac_data$LR == "fraction_ligand_group"]   <- "ligand"
  frac_data$LR[frac_data$LR == "fraction_receptor_group"] <- "receptor"

  cs_data <- cs_data %>% inner_join(frac_data)

  p_cs <- cs_data %>%
    ggplot(aes(LR, lr_interaction,
               color = celltype_specificity, size = fraction_expression)) +
    geom_point() +
    facet_grid(sender_receiver ~ group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(expand = y_expnd) +
    theme_light() +
    viridis::scale_color_viridis(limits = cs_limits) +
    theme(
      axis.ticks = element_blank(), axis.title = element_blank(),
      axis.title.y = element_blank(), axis.text.y = element_blank(),
      axis.text.x  = element_text(size = 9, angle = 90, hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.spacing.x  = unit(0.2, "lines"), panel.spacing.y = sp_y,
      strip.text.x     = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y     = element_blank(),
      strip.background = element_rect(color = "darkgrey", fill = "whitesmoke",
                                       size = 1.5, linetype = "solid")
    ) +
    labs(color = "Scaled celltype specificity",
         size  = "Fraction of expression")

  if (!is.null(widths)) {
    p <- patchwork::wrap_plots(p1, p_cs, nrow = 1,
                               guides = "collect", widths = widths)
  } else {
    p <- patchwork::wrap_plots(p1, p_cs, nrow = 1,
                               guides = "collect",
                               widths = c(
                                 sample_data$sample %>% unique() %>% length(),
                                 2 * (sample_data$group %>% unique() %>% length())
                               ))
  }
  return(p)
}

# ---- Parameters ---------------------------------------------------------------
# MultiNicheNet runs its own pseudobulk DE (edgeR LRT) internally, per cell type.
# It then scores each L-R pair by combining:
#   (1) ligand activity  — NicheNet prior-network AUPR against DE gene targets
#   (2) ligand expression in sender     — continuous, from pseudobulk
#   (3) receptor expression in receiver — continuous, from pseudobulk
#   (4) DE log-fold change of ligand / receptor between conditions
# All four components feed a single prioritization score (0–1).


ORGANISM <- "human"

TISSUES <- list(
  omentum = list(
    h5ad     = "inputs/anndatas/omentum_reannotated_finely_truncated.h5ad",
    part_col = "common_participant"
  ),
  subq = list(
    h5ad     = "inputs/anndatas/sq_reannotated_finely_truncated.h5ad",
    part_col = "common_participant"
  )
)

CELL_TYPE_COLS <- c("cell_type_", "cell_type2_") # 

COHORT_COL <- "cohort"

# Contrasts: muscat/edgeR syntax.
# "Unhealthy-Healthy": positive logFC = up in Unhealthy
# "Healthy-Unhealthy": positive logFC = up in Healthy
CONTRASTS_OI <- c("'Unhealthy-Healthy','Healthy-Unhealthy'")
CONTRAST_TBL <- tibble(
  contrast = c("Unhealthy-Healthy", "Healthy-Unhealthy"),
  group    = c("Unhealthy",         "Healthy")
)

# Filtering / DE parameters
MIN_CELLS         <- 10     # min cells per cell type per sample
FRACTION_CUTOFF   <- 0.05   # min fraction of samples expressing LR pair
MIN_SAMPLE_PROP   <- 0.50   # min proportion of samples a cell type must appear in
LOGFC_THRESHOLD   <- 0.50   # |log2FC| threshold for DE genes fed to ligand activity
P_VAL_THRESHOLD   <- 0.05   # FDR threshold for DE genes
TOP_N_TARGETS     <- 250    # top target genes per ligand for activity scoring
TOP_N_LR          <- 50     # top L-R pairs to visualise per contrast
N_CORES           <- 4      # parallel cores for NicheNet activity scoring
SCENARIO <- "regular"

BASE_DIR     <- "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq"
NICHENET_DIR <- file.path(BASE_DIR, "inputs/nichenet_v2")    # LT matrix (Zenodo 7074291)
NICHENET_MNN <- file.path(BASE_DIR, "inputs/nichenet_mnn")   # LR network (Zenodo 10229222)
OUTPUT_DIR   <- file.path(BASE_DIR, "results/multinichenet", SCENARIO)


# ---- Load NicheNet networks ---------------------------------------------------
# Per the multinichenetr vignette (saeyslab/multinichenetr, pairwise_analysis_MISC):
#   - ligand_target_matrix : Zenodo 7074291  (NicheNet v2)
#   - lr_network           : Zenodo 10229222 (human allInfo, with alias resolution)
# Both have target genes as ROWS and ligands as COLUMNS after processing.

dir.create(NICHENET_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(NICHENET_MNN, recursive = TRUE, showWarnings = FALSE)

# 1. Ligand-target matrix — Zenodo 7074291
lt_file <- file.path(NICHENET_DIR, "ligand_target_matrix_nsga2r_final.rds")
if (!file.exists(lt_file)) {
  cat("Downloading ligand_target_matrix from Zenodo 7074291...\n")
  download.file(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds",
    lt_file, mode = "wb"
  )
}

# 2. LR network — Zenodo 10229222
lr_file <- file.path(NICHENET_MNN, "lr_network_human_allInfo_30112033.rds")
if (!file.exists(lr_file)) {
  cat("Downloading lr_network from Zenodo 10229222...\n")
  download.file(
    "https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds",
    lr_file, mode = "wb"
  )
}

cat("Loading and preprocessing networks...\n")

# Load LR network and resolve gene aliases → make syntactically valid names
lr_network_all <- readRDS(lr_file) %>%
  mutate(
    ligand   = convert_alias_to_symbols(ligand,   organism = ORGANISM),
    receptor = convert_alias_to_symbols(receptor, organism = ORGANISM)
  ) %>%
  mutate(
    ligand   = make.names(ligand),
    receptor = make.names(receptor)
  )

lr_network <- lr_network_all %>% distinct(ligand, receptor)

# Load LT matrix and apply same name harmonisation
ligand_target_matrix <- readRDS(lt_file)
colnames(ligand_target_matrix) <- colnames(ligand_target_matrix) %>%
  convert_alias_to_symbols(organism = ORGANISM) %>% make.names()
rownames(ligand_target_matrix) <- rownames(ligand_target_matrix) %>%
  convert_alias_to_symbols(organism = ORGANISM) %>% make.names()

# Keep only LR pairs whose ligand appears in the LT matrix, then trim the matrix
lr_network <- lr_network %>% filter(ligand %in% colnames(ligand_target_matrix))
ligand_target_matrix <- ligand_target_matrix[, lr_network$ligand %>% unique()]

cat("  LT matrix:", nrow(ligand_target_matrix), "targets x",
    ncol(ligand_target_matrix), "ligands\n")
cat("  LR network:", nrow(lr_network), "pairs |",
    n_distinct(lr_network$ligand), "ligands\n")

# ---- Network–SCE compatibility check -----------------------------------------
check_network_sce <- function(sce, lt_matrix, lr_net) {
  # After harmonisation: colnames = ligands, rownames = target genes
  sce_genes   <- rownames(sce)
  lt_ligands  <- colnames(lt_matrix)
  lt_targets  <- rownames(lt_matrix)

  n_lig <- length(intersect(sce_genes, lt_ligands))
  n_tgt <- length(intersect(sce_genes, lt_targets))

  cat("  [Network check after SCE harmonisation]\n")
  cat("    SCE genes:          ", length(sce_genes), "\n")
  cat("    LT ligands (cols):  ", length(lt_ligands), "  (in SCE:", n_lig, ")\n")
  cat("    LT targets (rows):  ", length(lt_targets), "  (in SCE:", n_tgt, ")\n")
  cat("    LR ligand examples: ", paste(head(lr_net$ligand, 5), collapse = ", "), "\n")

  if (n_tgt < 1000)
    warning("Very few target genes match SCE (", n_tgt,
            "). Check gene name harmonisation — was makenames_SCE() applied?")
  invisible(list(n_lig = n_lig, n_tgt = n_tgt))
}

colors_lookup <-  c(
  "Plasma_cells" = "#6e5f5b", 
  "B_cells" = "#fc7447", # orange
  "IGFBP2_cells" = "#FD1486", # pink
  "IGFBP2" = "#FD1486", # pink
  "Mesothelial" = "#00546B", # dark green
  "Mesothelial_cells" = "#00546B", # dark green
  "Pre-adipocytes" = "#941DEC", # violet
  "Pre_adipocytes" = "#941DEC", # violet
  "T_cells" = "#FED000", # yellow
  "FAPs" = "#6EB8E7", # light blue
  "NK_cells" = "#408EF7", # dark blue
  "Myeloid_cells" = "#a12e2a", # dark red
  "ECs" = "#7bde18", # light green
  "Myofibroblasts" = "#de1818",
  'Monocytes' = '#E63946',
  'Neutrophils' = '#2196F3',
  'pDCs' = '#2DC653',
  'cDCs' = '#FF9800',
  'pDCS' = '#9C27B0',
  'LAMs' = '#00BCD4',
  'Macrophages' = '#00BCD4',
  'PVMs' = '#F06292',
  'CD9_ATM' = '#FFEB3B',
  'cDC1s' = '#795548',
  'cDC2Bs' = '#607D8B',
  'CD4_T_cells' = '#6D4C41',
  'Effector_CD8_T_cells' = '#FF5722',
  'CD8_T_cells' = '#3F51B5',
  'CD8_TRM' = '#8BC34A',
  'mNK1_cells' = '#F44336',
  'mNK2_cells' = '#673AB7',
  'mNK2_cell' = '#408EF7',
  'mNK3_cells' = '#F06292',
  'NA_CD4_T_cells' = '#CDDC39',
  'CM_CD4_T_cells' = '#1B998B',
  'GD_T_cells' = '#E91E63',
  'EM_CD4_T_cells' = '#FFC107',
  'Tregs' = '#D4A017',
  'PLCG2_cells' = '#A8DADC',
  "LECs" = "#FD1486", # pink
  "VECs" = "#99D969" # yellow
)

colors_lookup <- c(
  # B cells (pinks/magentas)
  "Plasma_cells"          = "#A23B72",
  "B_cells"               = "#F032E6",
  "IGFBP2_cells"          = "#F72585",
  "IGFBP2"                = "#F72585",
  "PLCG2_cells"           = "#EC407A",
  
  # Stromal (greens/olives)
  "Mesothelial"           = "#808000",
  "Mesothelial_cells"     = "#808000",
  "Pre-adipocytes"        = "#2D6A4F",
  "Pre_adipocytes"        = "#2D6A4F",
  "FAPs"                  = "#3CB44B",
  "Myofibroblasts"        = "#7B0000",
  
  # Endothelial (yellow-greens)
  "ECs"                   = "#90BE6D",
  "LECs"                  = "#FF5722",
  "VECs"                  = "#1B998B",
  
  # Myeloid / macrophages (reds → oranges)
  "Myeloid_cells"         = "#7B0000",
  "Monocytes"             = "#C62828",
  "Macrophages"           = "#E6194B",
  "LAMs"                  = "#FF5722",
  "PVMs"                  = "#FF9F1C",
  "CD9_ATM"               = "#F58231",
  
  # Dendritic cells (yellows / lime)
  "pDCs"                  = "#D4B800",
  "pDCS"                  = "#D4B800",
  "cDCs"                  = "#E9C46A",
  "cDC1s"                 = "#BFEF45",
  "cDC2Bs"                = "#9A6324",
  
  # Neutrophils
  "Neutrophils"           = "#6D4C41",
  
  # T cells (blues)
  "T_cells"               = "#4363D8",
  "CD4_T_cells"           = "#3A86FF",
  "NA_CD4_T_cells"        = "#7986CB",
  "CM_CD4_T_cells"        = "#0096C7",
  "EM_CD4_T_cells"        = "#00ACC1",
  "Tregs"                 = "#264653",
  "CD8_T_cells"           = "#311B92",
  "Effector_CD8_T_cells"  = "#000075",
  "CD8_TRM"               = "#42D4F4",
  "GD_T_cells"            = "#4B0082",
  
  # NK cells (purples)
  "NK_cells"              = "#7209B7",
  "mNK1_cells"            = "#911EB4",
  "mNK2_cells"            = "#AB47BC",
  "mNK2_cell"             = "#AB47BC",
  "mNK3_cells"            = "#C77DFF"
)


# ---- Visualisation helper ----------------------------------------------------
# contrast_tbl and ligand_target_matrix are accessed from the calling environment.
plot_multinichenet <- function(mnn_out, out_dir, contrast_tbl,
                               lt_matrix, top_n = TOP_N_LR) {

  # Generate one large qualitative palette for (senders + receivers) combined so
  # no colour is shared between the two maps. Senders take the first n_s slots,
  # receivers take the next n_r slots — all from the same evenly-spaced HCL wheel.
  # One base colour per cell type, chosen for maximum perceptual separation.
  # Senders use the full-saturation base colour; receivers use a lighter tint
  # of the same colour (~40% blended toward white), so same cell type = same hue
  # but role is immediately distinguishable by brightness.
  lighten <- function(hex, amount = 0.40) {
    m <- col2rgb(hex) / 255
    m <- m + (1 - m) * amount
    rgb(m[1, ], m[2, ], m[3, ])
  }

  make_sr_colors <- function(tbl) {
    cell_types  <- sort(unique(c(tbl$sender, tbl$receiver)))
    n           <- length(cell_types)

    # hcl.colors with "Dark 3" spaces hues evenly on the HCL wheel at fixed
    # chroma + luminance — gives perceptually uniform, maximally distinct colours.
    base_cols <- colors_lookup #setNames(hcl.colors(n, palette = "Dark 3", fixup = TRUE), cell_types)

    list(
      sender   = base_cols[sort(unique(tbl$sender))],
      receiver = setNames(lighten(base_cols[sort(unique(tbl$receiver))]),
                          sort(unique(tbl$receiver)))
    )
  }

  # --- 1. Cross-contrast circos (top pairs per group so every group has entries) -
  # rank_per_group = TRUE ensures each condition contributes pairs; FALSE can leave
  # one condition with zero entries, which triggers a zero-length unit vector error
  # in circlize when it tries to draw an empty sector.
  prio_all <- get_top_n_lr_pairs(mnn_out$prioritization_tables,
                                  top_n = top_n, rank_per_group = TRUE)
  prio_tbl_all <- mnn_out$prioritization_tables$group_prioritization_tbl %>%
    filter(id %in% prio_all$id) %>%
    distinct(id, sender, receiver, ligand, receptor, group) %>%
    left_join(prio_all)   # no by= → joins on all common columns, avoids .x/.y suffixes
  prio_tbl_all$prioritization_score[is.na(prio_tbl_all$prioritization_score)] <- 0

  n_groups_present <- n_distinct(prio_tbl_all$group)
  if (nrow(prio_tbl_all) > 0 && n_groups_present >= 2) {
    sr_colors <- make_sr_colors(prio_tbl_all)
    tryCatch({
      pdf(file.path(out_dir, "00_circos_all_groups.pdf"), 10, 10)
      multinichenetr::make_circos_group_comparison(prio_tbl_all,
                                   sr_colors$sender, sr_colors$receiver)
      dev.off()
    }, error = function(e) { dev.off(); cat("    Circos (all):", conditionMessage(e), "\n") })
  } else {
    cat("    Circos skipped: only", n_groups_present, "group(s) in top pairs\n")
  }

  # --- Per-contrast plots ------------------------------------------------------
  for (i in seq_len(nrow(contrast_tbl))) {
    group_oi   <- contrast_tbl$group[i]
    contrast_name <- contrast_tbl$contrast[i]
    safe_cname <- gsub("[^A-Za-z0-9]", "_", contrast_name)
    cat("    Plotting contrast:", contrast_name, "\n")

    # Top LR pairs for this group
    prio_group <- get_top_n_lr_pairs(mnn_out$prioritization_tables,
                                      top_n = top_n, groups_oi = group_oi)

    if (nrow(prio_group) == 0) {
      cat("      No prioritized pairs for this contrast, skipping\n"); next
    }

    # 2. Prioritization dot plot
    # get_top_n_lr_pairs() returns a slim table; fetch full columns from the
    # prioritization table so fraction_expressing_ligand_receptor is available.
    prio_group_full <- mnn_out$prioritization_tables$group_prioritization_tbl %>%
      filter(id %in% prio_group$id, group == group_oi)

    p_prio <- prio_group_full %>%
      mutate(lr_pair = paste0(ligand, "→", receptor),
             sr_pair = paste0(sender, "→", receiver)) %>%
      ggplot(aes(x = sr_pair, y = reorder(lr_pair, prioritization_score),
                 color = prioritization_score,
                 size  = fraction_expressing_ligand_receptor)) +
      geom_point() +
      scale_color_gradient(low = "lightgrey", high = "#b2182b",
                           name = "Priority\nscore") +
      scale_size_continuous(range = c(1, 6), name = "Fraction\nexpressing") +
      labs(title    = paste("Top L-R pairs —", contrast_name),
           subtitle = paste("Up in", group_oi),
           x = "Sender → Receiver", y = "Ligand → Receptor") +
      theme_classic(base_size = 9) +
      theme(axis.text.x      = element_text(angle = 45, hjust = 1),
            legend.position  = "right",
            panel.grid.major = element_line(colour = "grey88", linewidth = 0.3))

    ggsave(file.path(out_dir, paste0(safe_cname, "_01_prioritization_dotplot.pdf")),
           p_prio,
           width  = max(10, n_distinct(prio_group_full$sender) * 1.5 + 4),
           height = max(8,  nrow(prio_group_full) * 0.22 + 3))

    # 3. Sample LR expression + activity bubble plot
    tryCatch({
      p_sample <- make_sample_lr_prod_activity_plots(
        mnn_out$prioritization_tables, prio_group)
      ggsave(file.path(out_dir, paste0(safe_cname, "_02_sample_lr_activity.pdf")),
             p_sample, width = 14, height = max(6, nrow(prio_group) * 0.25 + 3))
    }, error = function(e) cat("      Sample-LR plot:", conditionMessage(e), "\n"))

    # 4. Ligand activity + target gene heatmap (per top receiver)
    # mnn_out_targets is a copy of mnn_out whose group_prioritization_tbl is
    # restricted to LR pairs that have at least one non-NA target in any contrast.
    # This is built once here (contrast-agnostic) so the function's internal
    # lookups into prioritization_tables stay consistent.
    ligands_with_any_target <- mnn_out$ligand_activities_targets_DEgenes$ligand_activities %>%
      filter(!is.na(target)) %>%
      distinct(receiver, ligand)

    mnn_out_targets <- mnn_out
    mnn_out_targets$prioritization_tables$group_prioritization_tbl <-
      mnn_out$prioritization_tables$group_prioritization_tbl %>%
      filter(group == group_oi) %>%
      filter(group == top_group) %>%
      semi_join(ligands_with_any_target, by = c("receiver", "ligand"))

    prio_group_targets <- get_top_n_lr_pairs(mnn_out_targets$prioritization_tables,
                                      top_n = 5000, groups_oi = group_oi)

    top_receivers <- prio_group_targets %>%
      dplyr::count(receiver, sort = TRUE) %>%
      # slice_head(n = 3) %>%
      pull(receiver)

    for (recv in top_receivers) {
      prio_recv_full <- prio_group_targets %>%
        filter(receiver == recv)

      if (nrow(prio_recv_full) == 0) next

      safe_recv <- gsub("[^A-Za-z0-9]", "_", recv)
      tryCatch({
        p_act <- make_ligand_activity_target_plot(
          group_oi,
          recv,
          prio_recv_full,
          mnn_out_targets$prioritization_tables,
          mnn_out_targets$ligand_activities_targets_DEgenes,
          contrast_tbl,
          mnn_out_targets$grouping_tbl,
          mnn_out_targets$celltype_info,
          lt_matrix,
          plot_legend = FALSE
        )
        if (!is.null(p_act)) {
          pdf(file.path(out_dir,
                paste0(safe_cname, "_03_ligand_activity_", safe_recv, ".pdf")),
              16, 10)
          # Function returns a list of mixed plot types (Heatmap + ggplot)
          plot_list <- if (is.list(p_act) &&
                           !inherits(p_act, c("Heatmap", "HeatmapList", "gg")))
            p_act else list(p_act)
          for (pl in plot_list) {
            # grid::grid.newpage()
            if (inherits(pl, c("Heatmap", "HeatmapList"))) draw(pl)
            else print(pl)
          }
          dev.off()
        }
      }, error = function(e)
        cat("      Ligand-activity (", recv, "):", conditionMessage(e), "\n"))
    }
  }
}

# ---- Main loop ----------------------------------------------------------------
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

for (tissue_name in names(TISSUES)) {
  cfg <- TISSUES[[tissue_name]]
  cat("\n=== Tissue:", tissue_name, "===\n")

  cat("  Loading h5ad...\n")
  sce <- readH5AD(file.path(BASE_DIR, cfg$h5ad), use_hdf5 = FALSE)

  # Ensure raw integer counts (MultiNicheNet requires counts assay)
  if (!"counts" %in% assayNames(sce))
    stop("counts assay not found in ", cfg$h5ad)

  # Harmonise gene names to match make.names()-processed networks.
  # This converts e.g. "HLA-A" -> "HLA.A" so joins on gene names succeed.
  cat("  Harmonising SCE gene names...\n")
  sce <- alias_to_symbol_SCE(sce, organism = ORGANISM) %>% makenames_SCE()

  # Harmonise participant column → "sample_id"
  if (cfg$part_col %in% colnames(colData(sce))) {
    colData(sce)$sample_id <- colData(sce)[[cfg$part_col]]
  } else {
    stop("Participant column '", cfg$part_col, "' not found in colData")
  }

  # Make colData character columns syntactically valid (required by multinichenetr)
  colData(sce)[[COHORT_COL]]  <- make.names(colData(sce)[[COHORT_COL]])
  colData(sce)[["sample_id"]] <- make.names(colData(sce)[["sample_id"]])

  # Sync CONTRAST_TBL group names in case make.names changed them
  CONTRAST_TBL_safe <- CONTRAST_TBL %>%
    mutate(group = make.names(group))

  # Also harmonise contrasts_oi to match made-names group labels
  # "Unhealthy-Healthy" → groups Unhealthy, Healthy — make.names keeps these unchanged
  CONTRASTS_OI_safe <- CONTRASTS_OI

  cat("  Samples:", n_distinct(colData(sce)$sample_id),
      "| Cohort distribution:\n")
  print(table(colData(sce)[[COHORT_COL]]))

  # Diagnostic: verify gene name overlap after harmonisation
  check_network_sce(sce, ligand_target_matrix, lr_network)

  for (ct_col in CELL_TYPE_COLS) {
    if (!ct_col %in% colnames(colData(sce))) {
      cat("  Column", ct_col, "not found, skipping\n"); next
    }

    # Make cell type labels syntactically valid
    colData(sce)[[ct_col]] <- make.names(colData(sce)[[ct_col]])

    cat("\n  -- Cell type column:", ct_col, "--\n")
    cat("  Cell types:", paste(sort(unique(colData(sce)[[ct_col]])), collapse = ", "), "\n")

    out_dir  <- file.path(OUTPUT_DIR, tissue_name, ct_col)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    rds_path <- file.path(out_dir, "multinichenet_output.rds")

    if (file.exists(rds_path)) {
      cat("  Cached result found, loading:", rds_path, "\n")
      mnn_out <- readRDS(rds_path)
    } else {
      cat("  Running MultiNicheNet analysis...\n")

      senders_oi   <- unique(colData(sce)[[ct_col]])
      receivers_oi <- unique(colData(sce)[[ct_col]])

      mnn_out <- multi_nichenet_analysis(
        sce                  = sce,
        celltype_id          = ct_col,
        sample_id            = "sample_id",
        group_id             = COHORT_COL,
        batches              = NA,
        covariates           = NA,
        senders_oi           = senders_oi,
        receivers_oi         = receivers_oi,
        lr_network           = lr_network,
        ligand_target_matrix = ligand_target_matrix,
        contrasts_oi         = CONTRASTS_OI_safe,
        contrast_tbl         = CONTRAST_TBL_safe,
        min_cells            = MIN_CELLS,
        fraction_cutoff      = FRACTION_CUTOFF,
        min_sample_prop      = MIN_SAMPLE_PROP,
        logFC_threshold      = LOGFC_THRESHOLD,
        p_val_threshold      = P_VAL_THRESHOLD,
        p_val_adj            = TRUE,
        empirical_pval       = FALSE,
        top_n_target         = TOP_N_TARGETS,
        n.cores              = N_CORES,
        scenario             = SCENARIO,
        ligand_activity_down = FALSE,
        verbose              = TRUE
      )

      saveRDS(mnn_out, rds_path)
      cat("  Saved:", rds_path, "\n")
    }

    # Save prioritization table
    prio_tbl <- mnn_out$prioritization_tables$group_prioritization_tbl
    write_csv(prio_tbl, file.path(out_dir, "prioritization_table.csv"))
    cat("  Prioritized L-R pairs:",
        nrow(filter(prio_tbl, prioritization_score > 0.5)), "with score > 0.5\n")

    # Top 10 per contrast
    cat("\n  Top 10 L-R pairs per contrast:\n")
    prio_tbl %>%
      filter(fraction_expressing_ligand_receptor >= FRACTION_CUTOFF) %>%
      group_by(contrast) %>%
      slice_max(prioritization_score, n = 10, with_ties = FALSE) %>%
      dplyr::select(contrast, sender, receiver, ligand, receptor, prioritization_score) %>%
      print(n = 30)

    cat("\n  Generating plots...\n")
    plot_multinichenet(mnn_out, out_dir,
                       contrast_tbl = CONTRAST_TBL_safe,
                       lt_matrix    = ligand_target_matrix)

    cat("  Done:", tissue_name, ct_col, "\n")
  }

  rm(sce); gc()
}

cat("\n=== MultiNicheNet complete ===\n")
cat("Results in:", OUTPUT_DIR, "\n")

# -------------------------------------
# ----- Some specific figures ---------
# -------------------------------------

tissue_name <- "subq"
# tissue_name <- "omentum"
ct_col <- "cell_type_"

out_dir  <- file.path(OUTPUT_DIR, tissue_name, ct_col)
rds_path <- file.path(out_dir, "multinichenet_output.rds")

mnn_out <- readRDS(rds_path)

contrast_tbl <- CONTRAST_TBL#_safe
lt_matrix <- ligand_target_matrix

top_n <- 25



########################################################################################
########################################################################################
# ----- Filter LR pairs by sender-receiver for the mesothelial figure ------
senders <- c("Mesothelial", "IGFBP2_cells")
receivers <- c("Mesothelial", "IGFBP2_cells")


mnn_filtered <- mnn_out
mnn_filtered$prioritization_tables$group_prioritization_tbl <- mnn_filtered$prioritization_tables$group_prioritization_tbl %>% 
  dplyr::filter(sender %in% senders | receiver %in% receivers)

i <- 1
group_oi   <- contrast_tbl$group[i]
contrast_name <- contrast_tbl$contrast[i]
safe_cname <- gsub("[^A-Za-z0-9]", "_", contrast_name)
cat("    Plotting contrast:", contrast_name, "\n")


# Top LR pairs for this group
prio_group <- get_top_n_lr_pairs(mnn_filtered$prioritization_tables,
                                 top_n = top_n, groups_oi = group_oi)


# 3. Sample LR expression + activity bubble plot

sender_receiver_order <- c(
  "Mesothelial --> ECs",
  "Mesothelial --> FAPs",
  "Mesothelial --> Mesothelial",
  "ECs --> Mesothelial",
  "Myeloid_cells --> Mesothelial",
  "NK_cells --> Mesothelial",
  "Pre_adipocytes --> Mesothelial",
  "IGFBP2_cells --> Mesothelial",
  "IGFBP2_cells --> FAPs",
  "IGFBP2_cells --> Myeloid_cells"
)
# sender_receiver_order <- c(
#   "Mesothelial --> FAPs",
#   "Mesothelial --> Myeloid_cells",
#   "Myeloid_cells --> Mesothelial",
#   "IGFBP2_cells --> FAPs",
#   "IGFBP2_cells --> Myeloid_cells"
# )

tryCatch({
  p_sample <- make_sample_lr_prod_activity_plots(
    mnn_out$prioritization_tables,
    prio_group, 
    row_spacing = 0.001,
    sender_receiver_order = sender_receiver_order,
    group_labels = c("Healthy" = "Ob", "Unhealthy" = "Ob+T2D"),
    expr_limits = c(-2.5, 2.5),
    cs_limits = c(0.5, 1.01)
    )
  ggsave(file.path(out_dir, paste0(safe_cname, "_02_sample_lr_activity_meso.pdf")),
         p_sample, width = 10, height = max(7, nrow(prio_group) * 0.11 + 3))
}, error = function(e) cat("      Sample-LR plot:", conditionMessage(e), "\n"))


########################################################################################
########################################################################################
# ----- Filter LR pairs by sender-receiver for the circos figure ------


lighten <- function(hex, amount = 0.40) {
  m <- col2rgb(hex) / 255
  m <- m + (1 - m) * amount
  rgb(m[1, ], m[2, ], m[3, ])
}

make_sr_colors <- function(tbl) {
  cell_types  <- sort(unique(c(tbl$sender, tbl$receiver)))
  n           <- length(cell_types)
  
  # hcl.colors with "Dark 3" spaces hues evenly on the HCL wheel at fixed
  # chroma + luminance — gives perceptually uniform, maximally distinct colours.
  # base_cols <- setNames(hcl.colors(n, palette = "Dark 3", fixup = TRUE), cell_types)
  base_cols <- colors_lookup
  
  list(
    sender   = base_cols[sort(unique(tbl$sender))],
    receiver = setNames(lighten(base_cols[sort(unique(tbl$receiver))]),
                        sort(unique(tbl$receiver)))
  )
}


cts_found <- union(
  unique(mnn_out$prioritization_tables$group_prioritization_tbl$sender),
  unique(mnn_out$prioritization_tables$group_prioritization_tbl$receiver)
  )

lymphoid_re <- regex("(CD8|T_cells|B_cells|Plasma_cells|mNK|Tregs|PLCG2_cells)", ignore_case = TRUE)
myeloid_re <- regex("(Monocytes|Macrophages|DCs|Neutrophils|cDC)", ignore_case = TRUE)
others_re <- regex("(LECs|VECs|Pre_adipocytes|FAPs|Myofibroblasts)", ignore_case = TRUE)

# lymphoid
lymphoid_cts <- cts_found[str_detect(cts_found, lymphoid_re)]
nonlymphoid_cts <- cts_found[!str_detect(cts_found, lymphoid_re)]

# myeloid
myeloid_cts <- cts_found[str_detect(cts_found, myeloid_re)]
nonmyeloid_cts <- cts_found[!str_detect(cts_found, myeloid_re)]

# others (and non mesothelial)
others_cts <- cts_found[str_detect(cts_found, others_re)]
nonothers_cts <- cts_found[!str_detect(cts_found, others_re)]

cts_list <- list(
  "others" = others_cts,
  "myeloid" = myeloid_cts,
  "lymphoid" = lymphoid_cts
)

cts_name <- "others"

cts <- cts_list[[cts_name]]
mnn_filtered <- mnn_out
mnn_filtered$prioritization_tables$group_prioritization_tbl <- mnn_filtered$prioritization_tables$group_prioritization_tbl  %>% 
  dplyr::filter(sender %in% cts | receiver %in% cts) %>%
  dplyr::filter(!(sender %in% c("Mast_cells_") | receiver %in% c("Mast_cells_")))

i <- 1
group_oi   <- contrast_tbl$group[i]
contrast_name <- contrast_tbl$contrast[i]
safe_cname <- gsub("[^A-Za-z0-9]", "_", contrast_name)
cat("    Plotting contrast:", contrast_name, "\n")

top_n <- 50
# Top LR pairs for this group
prio_group <- get_top_n_lr_pairs(mnn_filtered$prioritization_tables,
                                 top_n = top_n, groups_oi = group_oi)


# --- 1. Cross-contrast circos (top pairs per group so every group has entries) -
# rank_per_group = TRUE ensures each condition contributes pairs; FALSE can leave
# one condition with zero entries, which triggers a zero-length unit vector error
# in circlize when it tries to draw an empty sector.
prio_all <- get_top_n_lr_pairs(mnn_filtered$prioritization_tables,
                               top_n = top_n, rank_per_group = TRUE)
prio_tbl_all <- mnn_out$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prio_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>%
  left_join(prio_all)   # no by= → joins on all common columns, avoids .x/.y suffixes
prio_tbl_all$prioritization_score[is.na(prio_tbl_all$prioritization_score)] <- 0
cts_name <- "all"
write.csv(prio_tbl_all, file.path(out_dir, paste0("prio_tbl_all.", tissue_name, ".", cts_name, ".csv")))

n_groups_present <- n_distinct(prio_tbl_all$group)
if (nrow(prio_tbl_all) > 0 && n_groups_present >= 2) {
  sr_colors <- make_sr_colors(prio_tbl_all)
  tryCatch({
    pdf(file.path(out_dir, paste0("00_circos_all_groups_", cts_name, ".pdf")), 10, 10)
    make_circos_group_comparison(prio_tbl_all,
                                 sr_colors$sender, sr_colors$receiver)
    dev.off()
  }, error = function(e) { dev.off(); cat("    Circos (all):", conditionMessage(e), "\n") })
} else {
  cat("    Circos skipped: only", n_groups_present, "group(s) in top pairs\n")
}

# Violin comparison of ligand, receptor and target expression by cell type and
# group for a set of ligands in a MultiNicheNet result object.
#
# Unlike the CellChat version (which filters by pathway_name), this function
# takes ligands directly because the MNN prioritization table has no
# pathway_name column. Ligands → receptors are resolved via the prioritization
# table; target genes are taken from ligand_activities_targets_DEgenes.
#
# Expression values are pseudobulk per sample (pb_df), one dot per sample in
# the jitter layer, which reflects the N used in edgeR DE internally.
#
# Args:
#   mnn_out       : multinichenet result object (readRDS of multinichenet_output.rds)
#   ligands_oi    : character vector of ligands of interest
#   senders_oi    : restrict ligand expression to these sender cell types (NULL = all)
#   receivers_oi  : restrict receptor / target expression to these receivers (NULL = all)
#   contrast_oi   : contrast name for target lookup, e.g. "Unhealthy-Healthy"
#                   (NULL = first contrast in ligand_activities table)
#   top_n_targets : number of top target genes per receiver to include (0 = no targets)
#   groups_oi     : groups to show; default NULL uses all groups in grouping_tbl
#   group_colors  : named color vector for groups; NULL auto-assigns two colors
#   ncol          : columns in facet_wrap
#   point_size    : size of per-sample jitter points
plot_lr_violin_mnn <- function(
    mnn_out,
    ligands_oi,
    senders_oi    = NULL,
    receivers_oi  = NULL,
    contrast_oi   = NULL,
    top_n_targets = 5,
    groups_oi     = NULL,
    group_colors  = NULL,
    ncol          = 3,
    point_size    = 1.5
) {
  prio_tbl <- mnn_out$prioritization_tables$group_prioritization_tbl

  # ---- Resolve ligands → receptors via prioritization table ---------------
  lr_tbl <- prio_tbl %>%
    dplyr::filter(ligand %in% ligands_oi) %>%
    { if (!is.null(senders_oi))   dplyr::filter(., sender   %in% senders_oi)   else . } %>%
    { if (!is.null(receivers_oi)) dplyr::filter(., receiver %in% receivers_oi) else . } %>%
    dplyr::distinct(ligand, receptor, sender, receiver)

  if (nrow(lr_tbl) == 0)
    stop("No L-R pairs found for the supplied ligands / sender / receiver filters.")

  ligands   <- unique(lr_tbl$ligand)
  receptors <- unique(lr_tbl$receptor)

  # ---- Resolve target genes ------------------------------------------------
  act_df <- mnn_out$ligand_activities_targets_DEgenes$ligand_activities

  if (is.null(contrast_oi))
    contrast_oi <- act_df$contrast[!is.na(act_df$contrast)][1]

  if (top_n_targets > 0) {
    targets <- act_df %>%
      dplyr::filter(
        ligand %in% ligands,
        contrast == contrast_oi,
        !is.na(target),
        !is.na(activity_scaled)
      ) %>%
      { if (!is.null(receivers_oi)) dplyr::filter(., receiver %in% receivers_oi) else . } %>%
      dplyr::group_by(receiver) %>%
      dplyr::slice_max(activity_scaled, n = top_n_targets, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::pull(target) %>%
      unique()
  } else {
    targets <- character(0)
  }

  # ---- Build gene-role table (ligands first, then receptors, then targets) -
  gene_roles <- dplyr::bind_rows(
    data.frame(gene = ligands,   role = "Ligand"),
    data.frame(gene = receptors, role = "Receptor"),
    data.frame(gene = targets,   role = "Target")
  ) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::mutate(role = factor(role, levels = c("Ligand", "Receptor", "Target")))

  # ---- Pull pseudobulk expression and join group labels --------------------
  pb <- mnn_out$celltype_info$pb_df %>%
    dplyr::inner_join(mnn_out$grouping_tbl, by = "sample") %>%
    dplyr::filter(gene %in% gene_roles$gene) %>%
    { if (!is.null(groups_oi)) dplyr::filter(., group %in% groups_oi) else . }

  # For ligands show expression in sender cell types; for receptors and targets
  # show expression in receiver cell types (when the user supplied restrictions).
  if (!is.null(senders_oi) || !is.null(receivers_oi)) {
    pb_lig <- pb %>%
      dplyr::filter(gene %in% ligands) %>%
      { if (!is.null(senders_oi)) dplyr::filter(., celltype %in% senders_oi) else . }

    pb_other <- pb %>%
      dplyr::filter(gene %in% c(receptors, targets)) %>%
      { if (!is.null(receivers_oi)) dplyr::filter(., celltype %in% receivers_oi) else . }

    pb <- dplyr::bind_rows(pb_lig, pb_other)
  }

  if (nrow(pb) == 0)
    stop("No expression data remained after filtering. Check cell type names.")

  pb <- pb %>%
    dplyr::left_join(gene_roles, by = "gene") %>%
    dplyr::mutate(
      facet_label = paste0(gene, "\n[", role, "]"),
      facet_label = factor(
        facet_label,
        levels = gene_roles %>%
          dplyr::mutate(fl = paste0(gene, "\n[", role, "]")) %>%
          dplyr::pull(fl)
      )
    )

  # ---- Colors --------------------------------------------------------------
  grps <- if (!is.null(groups_oi)) groups_oi else unique(mnn_out$grouping_tbl$group)

  if (is.null(group_colors)) {
    default_pal <- c("#4dac26", "#d01c8b", "#1f78b4", "#ff7f00")
    group_colors <- setNames(default_pal[seq_along(grps)], grps)
  }

  # ---- Plot ----------------------------------------------------------------
  ggplot(pb, aes(x = celltype, y = pb_sample, fill = group)) +
    geom_violin(
      position  = position_dodge(0.8),
      scale     = "width",
      trim      = TRUE,
      linewidth = 0.3
    ) +
    geom_point(
      aes(color = group),
      position    = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
      size        = point_size,
      alpha       = 0.7,
      show.legend = FALSE
    ) +
    facet_wrap(~ facet_label, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values  = group_colors) +
    scale_color_manual(values = group_colors) +
    labs(
      x    = NULL,
      y    = "Pseudobulk expression",
      fill = NULL,
      subtitle = paste0(
        "Ligands: ", paste(ligands, collapse = ", "),
        if (length(targets) > 0) paste0(" | Targets: contrast ", contrast_oi) else ""
      )
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = 9),
      legend.position  = "top"
    )
}

# Plot pseudobulk expression of an arbitrary gene set in selected cell types,
# split by group. One facet per gene, x-axis = cell type.
#
# Args:
#   mnn_out      : multinichenet result object
#   genes_oi     : character vector of genes to plot
#   celltypes_oi : character vector of cell types to show (NULL = all)
#   groups_oi    : character vector of groups to show (NULL = all)
#   group_colors : named color vector; NULL auto-assigns
#   ncol         : columns in facet_wrap
#   point_size   : size of per-sample jitter points
plot_pb_expression <- function(
    mnn_out,
    genes_oi,
    celltypes_oi = NULL,
    groups_oi    = NULL,
    group_colors = NULL,
    ncol         = 3,
    point_size   = 1.5
) {
  pb <- mnn_out$celltype_info$pb_df %>%
    dplyr::inner_join(mnn_out$grouping_tbl, by = "sample") %>%
    dplyr::filter(gene %in% genes_oi) %>%
    { if (!is.null(celltypes_oi)) dplyr::filter(., celltype %in% celltypes_oi) else . } %>%
    { if (!is.null(groups_oi))    dplyr::filter(., group    %in% groups_oi)    else . } %>%
    dplyr::mutate(gene = factor(gene, levels = genes_oi))

  if (nrow(pb) == 0)
    stop("No data after filtering. Check gene and cell type names.")

  grps <- if (!is.null(groups_oi)) groups_oi else unique(mnn_out$grouping_tbl$group)

  if (is.null(group_colors)) {
    default_pal  <- c("#4dac26", "#d01c8b", "#1f78b4", "#ff7f00")
    group_colors <- setNames(default_pal[seq_along(grps)], grps)
  }

  ggplot(pb, aes(x = celltype, y = pb_sample, fill = group)) +
    geom_violin(
      position  = position_dodge(0.8),
      scale     = "width",
      trim      = TRUE,
      linewidth = 0.3
    ) +
    geom_point(
      aes(color = group),
      position    = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
      size        = point_size,
      alpha       = 0.7,
      show.legend = FALSE
    ) +
    facet_wrap(~ gene, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values  = group_colors) +
    scale_color_manual(values = group_colors) +
    labs(x = NULL, y = "Pseudobulk expression", fill = NULL) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = 9),
      legend.position  = "top"
    )
}

# Plot pseudobulk expression of an arbitrary gene set read directly from a
# pseudobulk AnnData (h5ad). Counts in assay "X" are log1p-CPM normalised
# internally. Gene symbols are resolved via rowData$entrez_name.
#
# Args:
#   h5ad_path     : path to the pseudobulk h5ad file
#   genes_oi      : character vector of gene symbols to plot
#   celltypes_oi  : cell types to include (NULL = all; matched on celltype_col)
#   groups_oi     : groups to include (NULL = all; matched on group_col)
#   group_col     : colData column for group labels (default "cohort")
#   celltype_col  : colData column for cell type labels (default "cell_type2_")
#   sample_col    : colData column for sample ID (default "common_participant")
#   group_colors  : named color vector for groups; NULL auto-assigns
#   ncol          : columns in facet_wrap
#   point_size    : size of per-sample jitter points
plot_pb_expression_h5ad <- function(
    h5ad_path,
    genes_oi,
    celltypes_oi = NULL,
    groups_oi    = NULL,
    group_col    = "cohort",
    celltype_col = "cell_type2_",
    sample_col   = "common_participant",
    group_colors = NULL,
    ncol         = 3,
    point_size   = 1.5
) {
  sce <- zellkonverter::readH5AD(h5ad_path, verbose = FALSE)

  # Map Entrez IDs (rownames) → gene symbols via rowData$entrez_name
  gene_map    <- setNames(rowData(sce)$entrez_name, rownames(sce))
  genes_found <- names(gene_map)[gene_map %in% genes_oi]
  missing     <- setdiff(genes_oi, gene_map[genes_found])
  if (length(missing) > 0)
    message("Genes not found in h5ad: ", paste(missing, collapse = ", "))
  if (length(genes_found) == 0)
    stop("None of the requested genes were found.")

  # log1p-CPM: normalise each observation column by its library size
  counts_full <- assay(sce, "X")
  lib_sizes   <- colSums(counts_full)
  mat         <- counts_full[genes_found, , drop = FALSE]
  rownames(mat) <- gene_map[genes_found]
  mat <- log1p(sweep(mat, 2, lib_sizes / 1e6, "/"))

  # Observation metadata
  cd <- as.data.frame(colData(sce)) |>
    tibble::rownames_to_column("obs_id") |>
    dplyr::select(
      obs_id,
      sample   = dplyr::all_of(sample_col),
      group    = dplyr::all_of(group_col),
      celltype = dplyr::all_of(celltype_col)
    )

  # Pivot to long and attach metadata
  present_genes <- genes_oi[genes_oi %in% gene_map]

  df <- as.data.frame(t(mat)) |>
    tibble::rownames_to_column("obs_id") |>
    tidyr::pivot_longer(cols = -obs_id, names_to = "gene", values_to = "expression") |>
    dplyr::left_join(cd, by = "obs_id") |>
    dplyr::filter(if (is.null(celltypes_oi)) TRUE else celltype %in% celltypes_oi) |>
    dplyr::filter(if (is.null(groups_oi))    TRUE else group    %in% groups_oi)    |>
    dplyr::mutate(gene = factor(gene, levels = present_genes))

  if (nrow(df) == 0)
    stop("No data after filtering. Check cell type and group names.")

  grps <- if (!is.null(groups_oi)) groups_oi else unique(df$group)
  if (is.null(group_colors)) {
    default_pal  <- c("#4dac26", "#d01c8b", "#1f78b4", "#ff7f00")
    group_colors <- setNames(default_pal[seq_along(grps)], grps)
  }

  ggplot(df, aes(x = celltype, y = expression, fill = group)) +
    geom_violin(
      position  = position_dodge(0.8),
      scale     = "width",
      trim      = TRUE,
      linewidth = 0.3,
      alpha     = 0.3,
    ) +
    geom_point(
      aes(color = group),
      position    = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
      size        = point_size,
      alpha       = 1,
      show.legend = FALSE
    ) +
    facet_wrap(~ gene, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values  = group_colors) +
    scale_color_manual(values = group_colors) +
    labs(x = NULL, y = "log1p(CPM)", fill = NULL) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = 9),
      legend.position  = "top"
    )
}

################################################################
### --------- Plotting Signaling PB expression! ----------------
################################################################
contrast_oi <- "Unhealthy-Healthy"
ligands_oi <- c("BMP2", "BMP6", "CSF1", "IL10")
ncols <- 4

p <- plot_lr_violin_mnn(
  mnn_out,
  ligands_oi    = ligands_oi,
  senders_oi    = NULL, #"Mesothelial",
  receivers_oi  = NULL, #"Myeloid_cells", #c("Macrophages", "Monocytes", "Neutrophils", "cDCs", "pDCs"),
  contrast_oi   = contrast_oi,
  top_n_targets = 20,
  ncol          = ncols
)
p


nrows <- length(unique(p$data$facet_label)) %/% ncols + 1

ggsave(file.path(out_dir, paste0(contrast_oi, "_pb_expression_violins_", paste(ligands_oi, collapse="_"), ".pdf")),
       p, width = 12, height = max(6, nrows * 0.5 + 3))


################################################################
### --------- Plotting just some PB expression! ----------------
################################################################

genes_oi <- c("RXRA", "CD300E", "ADGRE1")

p <- plot_pb_expression(
    mnn_out,
    genes_oi,
    celltypes_oi = myeloid_cts,
    groups_oi    = NULL,
    group_colors = NULL,
    ncol         = 3,
    point_size   = 1.5
)
p


################################################################
### --------- Plotting just some PB expression from h5ad! ----------------
################################################################

h5ad_path <- "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq/inputs/anndatas/omentum_reannotated_finely_pb_cell_type2_.h5ad"
genes_oi <- c("RXRA", "CD300E", "ADGRE1")

p <- plot_pb_expression_h5ad(
    h5ad_path,
    genes_oi = genes_oi,
    celltypes_oi = myeloid_cts,
    groups_oi    = NULL,
    group_col    = "cohort",
    celltype_col = "cell_type2_",
    sample_col   = "common_participant",
    group_colors = NULL,
    ncol         = 3,
    point_size   = 1.5
)
p
