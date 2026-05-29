suppressPackageStartupMessages({
  library(pathview)
  library(magick)
  library(tidyverse)
})

# ---- Parameters ---------------------------------------------------------------
# Score: sign(coef) * -log10(FDR)
#   Positive = significantly up in Unhealthy (red)
#   Negative = significantly up in Healthy   (blue)
# GENE_LIMIT: colour scale clipped at ±GENE_LIMIT (-log10 units).
#   -log10(0.05) ≈ 1.3  |  -log10(0.01) ≈ 2  |  -log10(1e-5) = 5

BASE_DIR   <- "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq"
OUTPUT_DIR <- file.path(BASE_DIR, "results/pathview")

TISSUES <- list(
  omentum = list(
    cell_type_  = "results/DE/omentum/MAST_results_per_cell_type_.csv",
    cell_type2_ = "results/DE/omentum/MAST_results_per_cell_type2_.csv"
  ),
  subq = list(
    cell_type_  = "results/DE/subq/MAST_results_per_cell_type_.csv",
    cell_type2_ = "results/DE/subq/MAST_results_per_cell_type2_.csv"
  )
)

# All normalised to "hsa" prefix (= Entrez IDs for human KEGG)
PATHWAYS <- c(
  "hsa04060",  # Cytokine-cytokine receptor interaction
  "hsa04668",  # TNF signaling pathway
  "hsa04064",  # NF-kappa B signaling pathway
  "hsa04350",  # TGF-beta signaling pathway
  "hsa04010",  # MAPK signaling pathway
  "hsa04931",  # Insulin resistance          (from map04931)
  "hsa00190",  # Oxidative phosphorylation   (from map00190)
  "hsa04068",  # FoxO signaling pathway      (from map04068)
  "hsa04024",  # cAMP signaling pathway
  "hsa04657",  # IL-17 signaling pathway
  "hsa04660",  # T cell receptor signaling pathway
  "hsa04662",  # B cell receptor signaling pathway
  "hsa04210",  # Apoptosis
  "hsa04519"   # Cadherin signaling
)

PATHWAY_IDS <- sub("^hsa", "", PATHWAYS)  # 5-digit IDs for pathview
GENE_LIMIT  <- 5                          # ±5 on the -log10(FDR) scale

# ---- Helper: build scored gene vector -----------------------------------------
# score = sign(coef) * -log10(FDR); FDR floored at 1e-300 to avoid Inf.
# If multiple genes share an Entrez ID keep the largest |score|.
make_gene_vec <- function(de_df, ct) {
  de_df %>%
    filter(cell_type == ct,
           !is.na(entrez_id), !is.na(coef), !is.na(FDR)) %>%
    mutate(
      entrez_str = as.character(as.integer(entrez_id)),
      score      = sign(coef) * (-log10(pmax(FDR, 1e-300)))
    ) %>%
    filter(entrez_str != "NA") %>%
    group_by(entrez_str) %>%
    slice_max(abs(score), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    { setNames(.$score, .$entrez_str) }
}

safe_name <- function(s) gsub("[^A-Za-z0-9]", "_", s)

# ---- Helper: combine PNGs for one cell type into a single PDF -----------------
pngs_to_pdf <- function(png_paths, pdf_path, pathway_labels) {
  existing <- file.exists(png_paths)
  if (!any(existing)) {
    cat("    No PNG outputs found to combine\n")
    return(invisible(NULL))
  }

  pages <- lapply(which(existing), function(i) {
    img <- image_read(png_paths[i])
    # Add a small label bar at the top so each page is identifiable
    label <- pathway_labels[i]
    img <- image_border(img, "white", "0x30")
    image_annotate(img, label, size = 22, gravity = "North",
                   location = "+0+5", color = "black", font = "sans-bold")
  })

  combined <- image_join(pages)
  image_write(combined, path = pdf_path, format = "pdf")
  cat("    PDF saved:", basename(pdf_path), "\n")
}

# ---- Main loop ----------------------------------------------------------------
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
orig_wd <- getwd()

for (tissue_name in names(TISSUES)) {
  tissue_cfg <- TISSUES[[tissue_name]]

  for (ct_col in names(tissue_cfg)) {
    de_path <- file.path(BASE_DIR, tissue_cfg[[ct_col]])
    cat("\n=== Tissue:", tissue_name, "| Cell-type column:", ct_col, "===\n")

    if (!file.exists(de_path)) {
      cat("  File not found, skipping:", de_path, "\n")
      next
    }

    de_all <- read_csv(de_path, show_col_types = FALSE)
    cat("  Loaded", nrow(de_all), "rows across",
        n_distinct(de_all$cell_type), "cell types\n")

    out_dir <- file.path(OUTPUT_DIR, tissue_name, ct_col)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    cell_types <- sort(unique(de_all$cell_type))

    for (ct in cell_types) {
      gene_vec <- make_gene_vec(de_all, ct)

      if (length(gene_vec) < 5) {
        cat("  Skipping", ct, "(only", length(gene_vec), "genes with Entrez ID)\n")
        next
      }

      cat("  [", ct, "]", length(gene_vec), "genes |",
          sum(gene_vec > 0), "up-Unhealthy,",
          sum(gene_vec < 0), "up-Healthy\n")
      safe_ct <- safe_name(ct)

      # pathview saves to the current working directory
      setwd(out_dir)

      for (i in seq_along(PATHWAY_IDS)) {
        
        
        pid    <- PATHWAY_IDS[i]
        hsa_id <- PATHWAYS[i]
        
        out_fn <- paste0(hsa_id, ".", safe_ct, ".png")
        
        if(file.exists(out_fn)) {
          next
        }
        
        cat("    ", hsa_id, "... ")

        result <- tryCatch({
          suppressMessages(
            pathview(
              gene.data   = gene_vec,
              pathway.id  = pid,
              species     = "hsa",
              gene.idtype = "entrez",
              out.suffix  = safe_ct,
              kegg.native = TRUE,
              same.layer  = TRUE,
              limit       = list(gene = GENE_LIMIT, cpd = 1),
              low         = list(gene = "steelblue", cpd = "blue"),
              mid         = list(gene = "white",     cpd = "gray"),
              high        = list(gene = "firebrick",  cpd = "red")
            )
          )
          "ok"
        }, error = function(e) {
          paste("ERROR:", conditionMessage(e))
        })

        cat(result, "\n")
      }

      # ---- Combine all pathway PNGs into one PDF per cell type ----------------
      png_paths      <- file.path(out_dir,
                                   paste0(PATHWAYS, ".", safe_ct, ".png"))
      pathway_labels <- paste0(PATHWAYS, " (", safe_ct, ")")
      pdf_path       <- file.path(out_dir, paste0(safe_ct, ".pdf"))

      pngs_to_pdf(png_paths, pdf_path, pathway_labels)
    }
  }
}

setwd(orig_wd)  # safety restore
cat("\n=== Pathview complete ===\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("Per-cell-type PDFs: {cell_type}.pdf  (one page per pathway)\n")
cat("Score: sign(coef) * -log10(FDR)  |  blue = Healthy up, red = Unhealthy up\n")
