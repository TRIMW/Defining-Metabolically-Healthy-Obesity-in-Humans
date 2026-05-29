library(tidyverse)
library(patchwork)
library(ggrepel)
library(ggnewscale)

setwd("/Users/willtrim/Documents/projs/Wills_analysis")

# ── Configuration ─────────────────────────────────────────────────────────────
DE_PATHS <- list(
  Omentum = c("obesity_scRNAseq/results/DE/omentum/MAST_results_per_cell_type2_.csv"),
  Subq    = c("obesity_scRNAseq/results/DE/subq/MAST_results_per_cell_type2_.csv") #"obesity_scRNAseq/results/DE/subq/MAST_results_per_cell_type_merged.csv", 
)
OUT_DIR    <- "obesity_scRNAseq/results/overview_figures"
CELL_TYPES <- c("Myofibroblasts", "VECs", "FAPs", "Pre_adipocytes")
# CELL_TYPES <- c("PVMs", "LAMs", "Monocytes", "DC", "CD4", "CD8", "NK", "B", "GD")
FDR_CUT    <- 0.05
N_LABEL    <- 8     # top genes labelled per panel in scatter

dir_colours    <- c(Up = "#B2182B", Down = "#2166AC")
tissue_colours <- c(Omentum = "#dba065", Subq = "#65db82")
sig_colours    <- c(
  Both    = "#4D4D4D",
  Omentum = "#dba065",
  Subq    = "#65db82",
  Neither = "grey85"
)
sig_colours_dark <- setNames(colorspace::darken(sig_colours, 0.25), names(sig_colours))

cell_types_map <- NULL
#   list(
#   "B_cells" = "B",
#   "CD4_T_cells" = "CD4",
#   "CD8_T_cells" = "CD8",
#   "GD_T_cells"  =  "GD"
# )

# ── Load & filter DE data ──────────────────────────────────────────────────────
de_all <- imap(DE_PATHS, \(paths, tissue) {
  map(paths, \(path)
    read_csv(path, show_col_types = FALSE) |>
      mutate(cell_type = if (length(cell_types_map) > 0)
               dplyr::recode(cell_type, !!!unlist(cell_types_map))
             else cell_type) |>
      filter(cell_type %in% CELL_TYPES) |>
      transmute(
        gene       = primerid,
        entrez_id  = entrez_id,
        ensembl_id = ensembl_id,
        logFC      = coef,
        padj       = FDR,
        cell_type,
        tissue     = tissue
      )
  ) |> list_rbind()
}) |> list_rbind()

# ── Panel A: DEG count barplot ─────────────────────────────────────────────────
counts <- de_all |>
  dplyr::filter(!is.na(padj), padj < FDR_CUT) |>
  dplyr::mutate(direction = if_else(logFC > 0, "Up", "Down")) |>
  dplyr::count(tissue, cell_type, direction) |>
  dplyr::mutate(
    n_plot    = if_else(direction == "Down", -n, n),
    cell_type = factor(cell_type, levels = CELL_TYPES),
    tissue    = factor(tissue, levels = c("Omentum", "Subq"),
                       labels = c("Omentum", "Subcutaneous"))
  )

fig_A <- ggplot(counts, aes(x = cell_type, y = n_plot, fill = direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
  facet_wrap(~ tissue, ncol = 1) +
  scale_fill_manual(values = dir_colours, name = NULL) +
  scale_y_continuous(
    breaks = c(-2000, -500, -100, -20, 0, 20, 100, 500, 2000),
    labels = \(x) abs(x),
    trans  = scales::pseudo_log_trans(base = 10, sigma = 1)
  ) +
  scale_x_discrete(labels = \(x) gsub("_", " ", x)) +
  labs(
    title    = "Significant DEGs  (FDR < 0.05)",
    subtitle = "Unhealthy vs Healthy  |  positive = up in Unhealthy",
    x = NULL, y = "# DEGs"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x        = element_text(angle = 40, hjust = 1),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", colour = "grey20", size=16),
    legend.position    = "top",
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
    axis.text         = element_text(size = 20),
    axis.title        = element_text(size = 20),
    legend.title      = element_text(size = 18),
    legend.text       = element_text(size = 16),
  )

# ── Diagnostic: common genes per cell_type between tissues ────────────────────
de_all |>
  dplyr::filter(!is.na(ensembl_id)) |>
  dplyr::group_by(cell_type) |>
  dplyr::summarise(
    n_omentum = dplyr::n_distinct(ensembl_id[tissue == "Omentum"]),
    n_subq    = dplyr::n_distinct(ensembl_id[tissue == "Subq"]),
    n_common  = length(intersect(
      unique(ensembl_id[tissue == "Omentum"]),
      unique(ensembl_id[tissue == "Subq"])
    )),
    .groups = "drop"
  ) |>
  print()

# ── Panel B: logFC scatter (omentum vs subq) ──────────────────────────────────
de_wide <- de_all |>
  dplyr::filter(!is.na(padj), !is.na(ensembl_id)) |>
  dplyr::group_by(ensembl_id, cell_type, tissue) |>
  dplyr::slice_min(padj, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  pivot_wider(
    id_cols     = c(ensembl_id, cell_type),
    names_from  = tissue,
    values_from = c(logFC, padj, gene)
  ) |>
  dplyr::filter(!is.na(logFC_Omentum), !is.na(logFC_Subq)) |>
  dplyr::mutate(
    sig_cat = case_when(
      padj_Omentum < FDR_CUT & padj_Subq < FDR_CUT  ~ "Both",
      padj_Omentum < FDR_CUT & padj_Subq >= FDR_CUT  ~ "Omentum",
      padj_Omentum >= FDR_CUT & padj_Subq < FDR_CUT  ~ "Subq",
      TRUE                                             ~ "Neither"
    ),
    sig_cat   = factor(sig_cat, levels = c("Both", "Omentum", "Subq", "Neither")),
    cell_type = factor(cell_type, levels = CELL_TYPES)
  )

de_wide$gene <- de_wide$gene_Omentum

ax_lim <- max(abs(c(de_wide$logFC_Omentum, de_wide$logFC_Subq)), na.rm = TRUE) * 1.05

label_dat <- de_wide |>
  filter(sig_cat != "Neither") |>
  mutate(.dist = sqrt(logFC_Omentum^2 + logFC_Subq^2)) |>
  dplyr::group_by(cell_type) |>
  dplyr::slice_max(.dist, n = N_LABEL, with_ties = FALSE) |>
  ungroup()

fig_B <- ggplot(de_wide, aes(x = logFC_Subq, y = logFC_Omentum, colour = sig_cat)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4, linetype = "dashed",
              colour = "grey50") +
  geom_point(data = \(d) filter(d, sig_cat == "Neither"),
             size = 0.5, alpha = 0.25) +
  geom_point(data = \(d) filter(d, sig_cat != "Neither"),
             size = 1.2, alpha = 0.85) +
  scale_colour_manual(values = sig_colours, name = "Significant in", drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  new_scale_colour() +
  geom_text_repel(
    data               = label_dat,
    aes(label          = gene,
        colour         = sig_cat),
    size               = 3.2,
    max.overlaps       = Inf,
    min.segment.length = 0,
    segment.size       = 0.3,
    segment.alpha      = 0.7,
    arrow              = arrow(length = unit(0.008, "npc"), type = "closed"),
    show.legend        = FALSE
  ) +
  scale_colour_manual(values = sig_colours_dark, guide = "none") +
  facet_wrap(~ cell_type, nrow = 1,
             labeller = as_labeller(\(x) gsub("_", " ", x))) +
  coord_cartesian(
    xlim = c(-ax_lim, ax_lim),
    ylim = c(-ax_lim, ax_lim)
  ) +
  labs(
    title    = "logFC concordance: Omentum vs Subcutaneous",
    subtitle = "Genes tested in both depots  |  dashed line = identity",
    x        = "Subcutaneous logFC",
    y        = "Omentum logFC"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = 14),
    axis.text         = element_text(size = 14),
    axis.title        = element_text(size = 18),
    legend.title      = element_text(size = 18),
    legend.text       = element_text(size = 20),
    legend.position   = "right",
    legend.key.size   = unit(1.2, "cm"),
    aspect.ratio      = 1,
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_line(colour = "grey92", linewidth = 0.25),
    panel.spacing     = unit(1.4, "cm")
  )

# ── Save ───────────────────────────────────────────────────────────────────────
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(OUT_DIR, "stromal_nonim_A_deg_counts.pdf"),
       fig_A, width = 4, height = 9, device = cairo_pdf)
message("Saved: stromal_A_deg_counts.pdf")

ggsave(file.path(OUT_DIR, "nonim_B_logfc_scatter.pdf"),
       fig_B, width = 12, height = 4.5, device = cairo_pdf)
message("Saved: stromal_B_logfc_scatter.pdf")

combined <- (fig_A | fig_B) +
  plot_layout(widths = c(1, 3.5)) +
  plot_annotation(
    title = "Stromal/structural cell types – DE overview (MAST, Unhealthy vs Healthy)",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave(file.path(OUT_DIR, "stromal_AB_overview.pdf"),
       combined, width = 20, height = 7, device = cairo_pdf)
message("Saved: stromal_AB_overview.pdf")

# ── Panel C: Broadly DE genes dotplot (logFC + FDR) ──────────────────────────
N_GENES_DOT <- 20

DE_PATHS_broad <- c(
  Omentum = "obesity_scRNAseq/results/DE/omentum/MAST_results_per_cell_type_.csv",
  Subq    = "obesity_scRNAseq/results/DE/subq/MAST_results_per_cell_type_.csv"
)

# Load DE for all cell types (not just the 5 stromal ones)
de_all_cts <- imap_dfr(DE_PATHS_broad, function(path, tissue) {
  read_csv(path, show_col_types = FALSE) |>
    transmute(
      gene      = primerid,
      logFC     = coef,
      padj      = FDR,
      cell_type,
      tissue    = tissue
    )
})

# Top genes ranked by number of cell types they are significantly DE in
# (regardless of direction, across both tissues); median |logFC| breaks ties
intersection_genes <- de_all_cts |>
  filter(!is.na(padj), padj < FDR_CUT,
         !str_detect(gene, "^RP[SL]\\d")) |>
  group_by(gene) |>
  summarise(
    n_cell_types = n_distinct(cell_type),
    med_logFC    = median(abs(logFC)),
    .groups      = "drop"
  ) |>
  arrange(desc(n_cell_types), desc(med_logFC)) |>
  slice_head(n = N_GENES_DOT)

genes_dot <- intersection_genes$gene
message(length(genes_dot), " broadly DE genes: ", paste(genes_dot, collapse = ", "))

dot_dat <- de_all_cts |>
  filter(gene %in% genes_dot, !is.na(padj)) |>
  mutate(
    neg_log_fdr = pmin(-log10(padj), 10),
    tissue      = factor(tissue, levels = c("Omentum", "Subq"),
                         labels = c("Omentum", "Subcutaneous")),
    gene        = factor(gene,
                         levels = intersection_genes |>
                           arrange(desc(n_cell_types), desc(med_logFC)) |>
                           pull(gene)),
    cell_type   = gsub("_", " ", cell_type)
  )

fc_lim <- max(abs(dot_dat$logFC), na.rm = TRUE)

fig_C <- ggplot(dot_dat, aes(x = cell_type, y = gene,
                              size = neg_log_fdr, colour = logFC)) +
  geom_point() +
  facet_wrap(~ tissue, nrow = 1) +
  scale_colour_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = c(-fc_lim, fc_lim),
    name     = "logFC\n(Ob+T2D/Ob)"
  ) +
  scale_size_continuous(
    name  = "-log10(FDR)",
    range = c(0.5, 6)
  ) +
  labs(
    title    = paste0("Top ", N_GENES_DOT,
                      " broadly DE genes  (FDR < ", FDR_CUT, ")"),
    subtitle = "Ranked by number of cell types DE in  |  colour direction can vary per cell type",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold"),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    axis.text         = element_text(size = 15),
    axis.title.y      = element_text(size = 16, margin = margin(r = 12)),
    legend.title      = element_text(size = 18),
  )

ggsave(file.path(OUT_DIR, "stromal_C_intersection_dotplot.pdf"),
       fig_C, width = 10, height = 6, device = cairo_pdf)
message("Saved: stromal_C_intersection_dotplot.pdf")
