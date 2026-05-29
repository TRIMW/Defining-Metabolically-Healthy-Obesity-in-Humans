library(tidyverse)

# ---- Paths ---------------------------------------------------------------
base_dir  <- "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq"
fgsea_dir <- file.path(base_dir, "results/DE/omentum/fgsea/cell_type_/MAST")
go_ref_fp <- "/Users/willtrim/Documents/reference/GO/go_bp_id_to_name.csv"
out_file  <- file.path(base_dir, "figures/mesothelial_fgsea_summary.pdf")

# ---- Load ----------------------------------------------------------------
bp_raw   <- read_csv(file.path(fgsea_dir, "Mesothelial.bp.csv"),      show_col_types = FALSE)
hallmark <- read_csv(file.path(fgsea_dir, "Mesothelial.hallmark.csv"), show_col_types = FALSE)
kegg     <- read_csv(file.path(fgsea_dir, "Mesothelial.kegg.csv"),     show_col_types = FALSE)
reactome <- read_csv(file.path(fgsea_dir, "Mesothelial.reactome.csv"), show_col_types = FALSE)
go_ref   <- read_csv(go_ref_fp, show_col_types = FALSE)

go_ref <- go_ref %>% distinct(msigdb_name, .keep_all = TRUE)


# ---- GO BP: keep level >= 3 ----------------------------------------------
# Level 1-3 are overly broad ("Growth", "Cytokine production", "Cell-cell
# adhesion"). Level 3+ gives biologically informative terms. Change to 3 for
# less aggressive filtering.
bp <- bp_raw |>
  left_join(go_ref |> select(msigdb_name, level), by = c("pathway" = "msigdb_name")) |>
  filter(!is.na(level), level >= 3)

# ---- Helpers -------------------------------------------------------------
clean_label <- function(x, prefix = NULL, to_sentence = TRUE, wrap_width = 38) {
  if (!is.null(prefix)) x <- str_remove(x, fixed(prefix))
  x <- str_replace_all(x, "_", " ")
  if (to_sentence) x <- str_to_sentence(x)
  str_wrap(x, width = wrap_width)
}

# Ribosomal genes heavily overlap viral pathway gene sets → artifactual hits
viral_re <- regex("(sars.cov|covid|coronavirus|influenza infection)", ignore_case = TRUE)

# ---- Custom dotplot for a hand-picked cross-database pathway selection ---
# selection: named list of pathway ID vectors, e.g.
#   list(Hallmark = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
#        KEGG     = c("KEGG_ADIPOGENESIS"),
#        `GO BP`  = c("GOBP_LIPID_METABOLIC_PROCESS"),
#        Reactome = c("REACTOME_FATTY_ACID_METABOLISM"))
# facet_by_db = FALSE → single panel, pathways ordered globally by NES (default)
# facet_by_db = TRUE  → one facet per database (same layout as the auto plot)
plot_selected_pathways <- function(
    selection,
    data_list = list(
      Hallmark = hallmark,
      KEGG     = kegg,
      `GO BP`  = bp,
      Reactome = reactome
    ),
    prefixes = list(
      Hallmark = "HALLMARK_",
      KEGG     = NULL,
      `GO BP`  = "GOBP_",
      Reactome = "REACTOME_"
    ),
    to_sentence = c(Hallmark = TRUE, KEGG = FALSE, `GO BP` = TRUE, Reactome = TRUE),
    facet_by_db = FALSE,
    wrap_width  = 38,
    title       = NULL
) {
  dat <- imap_dfr(selection, function(pathways, db_name) {
    data_list[[db_name]] |>
      filter(pathway %in% pathways) |>
      mutate(
        label = clean_label(
          pathway,
          prefix      = prefixes[[db_name]],
          to_sentence = isTRUE(to_sentence[[db_name]]),
          wrap_width  = wrap_width
        ),
        db = db_name
      )
  }) |>
    mutate(
      direction  = if_else(NES > 0, "Up", "Down"),
      neg_log10p = -log10(padj),
      db         = factor(db, levels = intersect(names(data_list), names(selection)))
    )

  if (facet_by_db && dplyr::n_distinct(dat$db) > 1) {
    dat <- dat |> arrange(db, NES)
  } else {
    dat <- dat |> arrange(NES)
  }
  dat <- dat |> mutate(y_key = factor(formatC(seq_len(n()), width = 4, flag = "0")))

  label_lookup <- setNames(dat$label, as.character(dat$y_key))

  p <- ggplot(dat, aes(x = NES, y = y_key, size = neg_log10p, color = direction)) +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "grey55", linetype = "dashed") +
    geom_point(alpha = 0.85) +
    scale_y_discrete(labels = label_lookup) +
    coord_cartesian(xlim = c(min(dat$NES) * 1.2, max(dat$NES) * 1.2)) +
    scale_color_manual(
      values = c("Up" = "#c0392b", "Down" = "#2980b9"),
      name   = NULL
    ) +
    scale_size_continuous(
      name   = bquote(-log[10](p[adj])),
      range  = c(1.5, 5.5),
      breaks = c(2, 5, 20, 40)
      # breaks = scales::pretty_breaks(n = 4)
    ) +
    labs(x = "Normalized Enrichment Score (NES)", y = NULL, title = title) +
    theme_bw(base_size = 8) +
    theme(
      strip.background   = element_rect(fill = "grey93", color = "grey65", linewidth = 0.4),
      strip.text         = element_text(face = "bold", size = 18),
      axis.text.y        = element_text(size = 14, lineheight = 0.9),
      axis.text.x        = element_text(size = 14),
      axis.title.x       = element_text(size = 14),
      panel.grid.major.x = element_line(linewidth = 0.25, color = "grey85"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(color = "grey65", linewidth = 0.4),
      legend.position    = "bottom",
      legend.box         = "vertical",
      legend.key.size    = unit(0.35, "cm"),
      legend.text        = element_text(size = 14),
      legend.title       = element_text(size = 18),
      plot.margin        = unit(c(0.5, 4, 0, 2), "cm")
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 12), order = 1),
      size  = guide_legend(nrow = 1, order = 2)
    )

  if (facet_by_db && dplyr::n_distinct(dat$db) > 1) {
    p <- p + facet_wrap(~ db, scales = "free_y", ncol = 1)
  }

  p
}

meso_selection <- list(
  Hallmark = c("HALLMARK_INTERFERON ALPHA RESPONSE",
               "HALLMARK_OXIDATIVE PHOSPHORYLATION",
               "HALLMARK_P53_PATHWAY",
               "HALLMARK_MITOTIC_SPINDLE",
               "HALLMARK_KRAS_SIGNALING_UP",
               "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
               "HALLMARK_G2M_CHECKPOINT",
               "HALLMARK_CHOLESTEROL_HOMEOSTASIS"#,
               # "HALLMARK_ANDROGEN_RESPONSE"
  ),
  KEGG     = c("Proteoglycans in cancer", 
               "IgSF CAM signaling", 
               "Tight junction", 
               "Cadherin signaling", 
               "Cytoskeleton in muscle cells"
  ),
  `GO BP`  = c("GOBP_CELL_JUNCTION_ORGANIZATION",
               "GOBP_CELL_CELL_ADHESION",
               "GOBP_WOUND_HEALING"
  ),
  Reactome = c(
    "REACTOME_SELENOAMINO_ACID_METABOLISM",
    "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",
    "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
    "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
    "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
    "REACTOME_NUCLEAR_ENVELOPE_NE_REASSEMBLY_RHO_GTPASES_ACTIVATE_IQGAPS"
  )
)

igfbp2_selection <- list(
  KEGG     = c("Antigen processing and presentation"),
  `GO BP`  = c("GOBP_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE",
               "GOBP_REGULATION_OF_MAP_KINASE_ACTIVITY",
   "REACTOME_Cytokine_signaling_in_immune system",
               "GOBP_REGULATION_OF_BMP_SIGNALING_PATHWAY",
               "GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION"
  ),
  Reactome = c(
    "REACTOME_ELASTIC_FIBRE_FORMATION",
    "REACTOME_DISEASES_OF_GLYCOSYLATION",
    "REACTOME_ECM_PROTEOGLYCANS",
    "REACTOME_EXTRACELLULAR_MATRIX_ORGANISATION",
    "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE SYSTEM",
    "REACTOME_INTERFERON_SIGNALING"
  )
  )
                  
p <- plot_selected_pathways(meso_selection)
out_file_sel <- file.path(base_dir, "figures/Mesothelial_fgsea_summary_together.pdf")
ggsave(out_file_sel, p, width = 7.5, height = 7, device = cairo_pdf)

# ---- Select top 5 significant pathways per direction per database --------
top_n <- 6

select_top <- function(df, prefix = NULL, to_sentence = TRUE) {
  df <- df |>
    filter(padj < 0.05, !str_detect(pathway, viral_re)) |>
    mutate(label = clean_label(pathway, prefix, to_sentence))

  bind_rows(
    df |> filter(NES > 0) |> arrange(padj) |> slice_head(n = top_n),
    df |> filter(NES < 0) |> arrange(padj) |> slice_head(n = top_n)
  )
}

dat <- bind_rows(
  select_top(hallmark, prefix = "HALLMARK_")     |> mutate(db = "Hallmark"),
  select_top(kegg,     to_sentence = FALSE)       |> mutate(db = "KEGG"),
  select_top(bp,       prefix = "GOBP_")          |> mutate(db = "GO BP"),
  select_top(reactome, prefix = "REACTOME_")      |> mutate(db = "Reactome")
) |>
  mutate(
    direction  = if_else(NES > 0, "Up", "Down"),
    neg_log10p = -log10(padj),
    db         = factor(db, levels = c("Hallmark", "KEGG", "GO BP", "Reactome"))
  )

# ---- Per-facet y-axis ordering by NES (low → high = bottom → top) -------
dat <- dat |>
  arrange(db, NES) |>
  dplyr::mutate(y_key = factor(
    paste(as.integer(db), formatC(seq_len(n()), width = 4, flag = "0"), sep = "_")
  ))

label_lookup <- setNames(dat$label, as.character(dat$y_key))

# ---- Plot ----------------------------------------------------------------
p <- ggplot(dat, aes(x = NES, y = y_key, size = neg_log10p, color = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey55", linetype = "dashed") +
  geom_point(alpha = 0.85) +
  scale_y_discrete(labels = label_lookup) +
  scale_color_manual(
    values = c("Up" = "#c0392b", "Down" = "#2980b9"),
    name   = NULL
  ) +
  scale_size_continuous(
    name   = bquote(-log[10](p[adj])),
    range  = c(1.5, 5.5),
    breaks = c(2, 5, 10, 20, 40)
  ) +
  facet_wrap(~ db, scales = "free_y", ncol = 1) +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL) +
  theme_bw(base_size = 8) +
  theme(
    strip.background   = element_rect(fill = "grey93", color = "grey65", linewidth = 0.4),
    strip.text         = element_text(face = "bold", size = 8),
    axis.text.y        = element_text(size = 12, lineheight = 0.9),
    axis.text.x        = element_text(size = 12),
    axis.title.x       = element_text(size = 12),
    panel.grid.major.x = element_line(linewidth = 0.25, color = "grey85"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "grey65", linewidth = 0.4),
    legend.position    = "bottom",
    legend.box         = "vertical",
    legend.key.size    = unit(0.35, "cm"),
    legend.text        = element_text(size = 7),
    legend.title       = element_text(size = 7),
    plot.margin        = unit(c(1,1,1,6), "cm")
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 6), order = 1),
    size  = guide_legend(nrow = 1, order = 2)
  )

ggsave(out_file, p, width = 7.5, height = 18, device = cairo_pdf)
message("Saved: ", out_file)
