suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(NMF)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
})

# ---- Parameters -------------------------------------------------------------
N_PATTERNS         <- 4    # NMF rank (number of communication patterns per condition)
N_TOP_DIFF         <- 3    # top differential interactions shown per sender cell type
MIN_SCORE          <- 0.5  # minimum overall prioritization_score to include
MIN_COMPONENT      <- 0.5  # both scaled_lfc_ligand AND scaled_lfc_receptor must exceed this
MIN_SAMPLE_FRAC    <- 0.75 # fraction of top-group samples where sender+receiver must both be
                           # present (abundance consistency filter)
N_TOP_FEATS        <- 30   # top ligands/receptors shown in H heatmap
NMF_NRUN           <- 10   # NMF restarts for consensus (reduce to 1 for speed)

# ---- Paths ------------------------------------------------------------------
MNN_PATHS <- list(
  omentum = "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq/results/multinichenet/regular/omentum/cell_type_/multinichenet_output.rds",
  subq    = "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq/results/multinichenet/regular/subq/cell_type_/multinichenet_output.rds"
)
OUTPUT_DIR <- "/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq/results/multinichenet/figures"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# # ---- Cell-type colours (same as cellchat_analysis.R) -----------------------
# CT_COLOURS <- c(
#   B_cells          = "#E41A1C",
#   ECs              = "#377EB8",
#   FAPs             = "#4DAF4A",
#   IGFBP2_cells     = "#984EA3",
#   Mesothelial      = "#FF7F00",
#   Myeloid_cells    = "#A65628",
#   Myofibroblasts   = "#F781BF",
#   NK_cells         = "#999999",
#   `Pre-adipocytes` = "#66C2A5",
#   Pre_adipocytes   = "#66C2A5",
#   Plasma_cells     = "#FC8D62",
#   T_cells          = "#8DA0CB"
# )
CT_COLOURS <- c(
  "Plasma_cells" = "#6e5f5b", 
  "B_cells" = "#fc7447", # orange
  "IGFBP2_cells" = "#FD1486", # pink
  "Mesothelial" = "#00546B", # dark green
  "Pre-adipocytes" = "#941DEC", # violet
  "Pre_adipocytes" = "#941DEC", # violet
  "T_cells" = "#FED000", # yellow
  "FAPs" = "#6EB8E7", # light blue
  "NK_cells" = "#408EF7", # dark blue
  "Myeloid_cells" = "#a12e2a", # dark red
  "ECs" = "#7bde18", # light green
  "Myofibroblasts" = "#de1818"
)

# ---- Shared colour scales ---------------------------------------------------
spectral_pal <- function(lo = 0, hi = 1, n = 11) {
  colorRamp2(seq(lo, hi, length.out = n), rev(brewer.pal(n, "Spectral")))
}

# ---- Helper: filter prioritization table ------------------------------------
# Three conservative filters applied in sequence:
#   1. Overall prioritization_score >= min_score
#   2. Both scaled_lfc_ligand AND scaled_lfc_receptor >= min_component
#      (prevents interactions driven by only one partner being DE)
#   3. Sender + receiver both present in >= min_sample_frac of the top-group
#      samples (removes interactions inconsistent across participants)
get_prio_tbl <- function(mnn,
                          min_score       = MIN_SCORE,
                          min_component   = MIN_COMPONENT,
                          min_sample_frac = MIN_SAMPLE_FRAC) {
  group_tbl  <- mnn$prioritization_tables$group_prioritization_tbl
  sample_tbl <- mnn$prioritization_tables$sample_prioritization_tbl

  # 1. Overall score + remove NAs
  tbl <- group_tbl[!is.na(group_tbl$prioritization_score) &
                     group_tbl$prioritization_score >= min_score, ]

  # # 2. Both ligand and receptor LFC components must individually be strong
  # tbl <- tbl[!is.na(tbl$scaled_lfc_ligand) &
  #              !is.na(tbl$scaled_lfc_receptor) &
  #              tbl$scaled_lfc_ligand   >= min_component &
  #              tbl$scaled_lfc_receptor >= min_component, ]

  # 3. Abundance consistency: require sender & receiver present in most samples
  #    of the interaction's top_group
  top_group_map <- unique(tbl[, c("id", "top_group")])

  consistency <- sample_tbl |>
    dplyr::inner_join(top_group_map, by = "id") |>
    dplyr::filter(group == top_group) |>
    dplyr::group_by(id) |>
    dplyr::summarise(
      frac_present = mean(keep_sender_receiver == "Sender & Receiver present",
                          na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::filter(frac_present >= min_sample_frac)

  tbl[tbl$id %in% consistency$id, ]
}

# ---- Helper: build cell-type × feature matrix (rows=ct, cols=ligand/receptor)
build_ct_feat_mat <- function(tbl, ct_col, feat_col) {
  long <- tbl |>
    dplyr::group_by(.data[[ct_col]], .data[[feat_col]]) |>
    dplyr::summarise(score = sum(prioritization_score, na.rm = TRUE), .groups = "drop")
  wide <- tidyr::pivot_wider(long,
                              names_from  = feat_col,
                              values_from = "score",
                              values_fill = 0)
  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- wide[[ct_col]]
  # Drop all-zero rows and columns
  mat <- mat[rowSums(mat) > 0, colSums(mat) > 0, drop = FALSE]
  mat
}

# ---- Helper: cosine similarity (rows of A vs rows of B) → ka × kb matrix --
cos_sim_mat <- function(A, B) {
  An <- A / sqrt(rowSums(A^2) + 1e-12)
  Bn <- B / sqrt(rowSums(B^2) + 1e-12)
  An %*% t(Bn)
}

# ---- Helper: greedy pattern matching by cosine similarity ------------------
# Returns index vector: match_idx[i] = which column of sim best matches row i
greedy_match <- function(sim_mat) {
  used <- rep(FALSE, ncol(sim_mat))
  idx  <- integer(nrow(sim_mat))
  for (i in seq_len(nrow(sim_mat))) {
    avail <- which(!used)
    best  <- avail[which.max(sim_mat[i, avail])]
    idx[i] <- best
    used[best] <- TRUE
  }
  idx
}

# ============================================================================
# 01  Sender / receiver scatter  (analogue of CellChat 07_*)
# ============================================================================
plot_sr_scatter <- function(tbl, title = NULL) {
  out_df <- tbl |>
    dplyr::group_by(group, sender) |>
    dplyr::summarise(outgoing = sum(prioritization_score, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::rename(cell_type = sender)

  in_df <- tbl |>
    dplyr::group_by(group, receiver) |>
    dplyr::summarise(incoming = sum(prioritization_score, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::rename(cell_type = receiver)

  sr <- dplyr::full_join(out_df, in_df, by = c("group", "cell_type")) |>
    tidyr::replace_na(list(outgoing = 0, incoming = 0))

  sr_h <- dplyr::filter(sr, group == "Healthy")
  sr_u <- dplyr::filter(sr, group == "Unhealthy")

  arrows_df <- dplyr::inner_join(
    dplyr::select(sr_h, cell_type, x0 = outgoing, y0 = incoming),
    dplyr::select(sr_u, cell_type, x1 = outgoing, y1 = incoming),
    by = "cell_type"
  )

  all_df     <- dplyr::bind_rows(sr_h, sr_u)
  present    <- unique(all_df$cell_type)
  ct_cols    <- structure(
    ifelse(is.na(CT_COLOURS[present]), "grey70", CT_COLOURS[present]),
    names = present
  )

  ggplot() +
    geom_segment(
      data      = arrows_df,
      aes(x = x0, y = y0, xend = x1, yend = y1, colour = cell_type),
      arrow     = arrow(length = unit(0.35, "cm"), type = "closed"),
      linewidth = 1.4, alpha = 0.25, show.legend = FALSE
    ) +
    geom_point(
      data  = dplyr::filter(all_df, group == "Healthy"),
      aes(x = outgoing, y = incoming, colour = cell_type),
      shape = 21, size = 6, stroke = 1.1, alpha = 0.45, fill = NA
    ) +
    geom_point(
      data  = dplyr::filter(all_df, group == "Unhealthy"),
      aes(x = outgoing, y = incoming, colour = cell_type),
      shape = 16, size = 6, alpha = 0.65
    ) +
    ggrepel::geom_text_repel(
      data  = all_df,
      aes(x = outgoing, y = incoming,
          label  = ifelse(group == "Unhealthy", gsub("_", " ", cell_type), ""),
          colour = cell_type),
      size               = 5.6,
      show.legend        = FALSE,
      max.overlaps       = Inf,
      force              = 12,
      force_pull         = 0.5,
      box.padding        = 0.5,
      point.padding      = 0.4,
      min.segment.length = 1,
      seed               = 12
    ) +
    scale_colour_manual(values = ct_cols) +
    labs(
      title   = title %||% "Sender / receiver signaling strength",
      x       = "Outgoing (sender) total prioritization score",
      y       = "Incoming (receiver) total prioritization score",
      colour  = NULL,
      caption = "Hollow = Healthy  |  Filled = Unhealthy  |  Arrows show direction of change"
    ) +
    theme_bw(base_size = 18) +
    theme(legend.position = "right")
}

# ============================================================================
# 02  Sender × receiver communication heatmap
# ============================================================================
plot_comm_heatmap <- function(tbl, title = NULL) {
  make_sr_mat <- function(grp) {
    m <- tbl |>
      dplyr::filter(group == grp) |>
      dplyr::group_by(sender, receiver) |>
      dplyr::summarise(score = sum(prioritization_score, na.rm = TRUE),
                       .groups = "drop") |>
      tidyr::pivot_wider(names_from = "receiver",
                          values_from = "score", values_fill = 0)
    mat <- as.matrix(m[, -1, drop = FALSE])
    rownames(mat) <- m$sender
    mat
  }

  mat_h <- make_sr_mat("Healthy")
  mat_u <- make_sr_mat("Unhealthy")

  # Align rows and columns to the union of both conditions
  all_snd <- sort(union(rownames(mat_h), rownames(mat_u)))
  all_rcv <- sort(union(colnames(mat_h), colnames(mat_u)))
  pad <- function(m) {
    out <- matrix(0, nrow = length(all_snd), ncol = length(all_rcv),
                  dimnames = list(all_snd, all_rcv))
    out[rownames(m), colnames(m)] <- m
    out
  }
  mat_h <- pad(mat_h)
  mat_u <- pad(mat_u)

  col_max  <- max(c(mat_h, mat_u))
  col_fun  <- spectral_pal(0, col_max)

  # Cluster on Healthy; reuse order for Unhealthy so panels are directly comparable
  row_dend <- as.dendrogram(hclust(dist(mat_h),   method = "ward.D2"))
  col_dend <- as.dendrogram(hclust(dist(t(mat_h)), method = "ward.D2"))

  base_args <- list(
    col             = col_fun,
    row_names_gp    = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 9),
    column_names_rot = 45,
    heatmap_legend_param = list(title = "Total\nscore", direction = "horizontal")
  )

  ht <- do.call(Heatmap, c(list(
    mat_h,
    name         = "Healthy",
    cluster_rows = row_dend,
    cluster_columns = col_dend,
    column_title = "Healthy",
    row_title    = "Sender →"
  ), base_args)) +
  do.call(Heatmap, c(list(
    mat_u,
    name            = "Unhealthy",
    cluster_rows    = row_dend,
    cluster_columns = col_dend,
    column_title    = "Unhealthy",
    show_row_dend   = FALSE,
    show_heatmap_legend = FALSE
  ), base_args))

  draw(ht,
       column_title    = title %||% "Sender → Receiver total prioritization score",
       column_title_gp = gpar(fontsize = 12, fontface = "bold"),
       heatmap_legend_side = "bottom", merge_legend = TRUE)
}

# ============================================================================
# 03/04  NMF pattern comparison  (analogue of CellChat 08b_*)
# ============================================================================
plot_nmf_patterns <- function(tbl, role = c("outgoing", "incoming"),
                               n_patterns = N_PATTERNS,
                               n_top      = N_TOP_FEATS,
                               out_pdf, width = 18, height = 7) {
  role     <- match.arg(role)
  ct_col   <- if (role == "outgoing") "sender"   else "receiver"
  feat_col <- if (role == "outgoing") "ligand"   else "receptor"
  feat_lbl <- if (role == "outgoing") "Ligand"   else "Receptor"

  # Build per-condition cell-type × feature matrices (native; no padding yet)
  mats <- lapply(c("Healthy", "Unhealthy"), function(grp)
    build_ct_feat_mat(dplyr::filter(tbl, group == grp), ct_col, feat_col))
  names(mats) <- c("Healthy", "Unhealthy")

  message("  Running NMF (rank=", n_patterns, ", nrun=", NMF_NRUN,
          ") for Healthy [", role, "]...")
  nmf_h <- NMF::nmf(mats$Healthy,   rank = n_patterns, method = "lee",
                     seed = 42, nrun = NMF_NRUN, .options = "v0")
  message("  Running NMF for Unhealthy [", role, "]...")
  nmf_u <- NMF::nmf(mats$Unhealthy, rank = n_patterns, method = "lee",
                     seed = 42, nrun = NMF_NRUN, .options = "v0")

  W_h <- NMF::basis(nmf_h)   # cell_type_h × k
  H_h <- NMF::coef(nmf_h)    # k × feature_h
  W_u <- NMF::basis(nmf_u)
  H_u <- NMF::coef(nmf_u)

  # Align W to the union of cell types (pad absent ones with 0)
  all_cts <- sort(union(rownames(W_h), rownames(W_u)))
  pad_w <- function(W) {
    out <- matrix(0, nrow = length(all_cts), ncol = ncol(W),
                  dimnames = list(all_cts, colnames(W)))
    out[rownames(W), ] <- W
    out
  }
  W_h <- pad_w(W_h)
  W_u <- pad_w(W_u)

  # Use the intersection of features for cosine similarity matching
  common_feats <- intersect(colnames(H_h), colnames(H_u))
  H_h_common   <- H_h[, common_feats, drop = FALSE]
  H_u_common   <- H_u[, common_feats, drop = FALSE]

  # Column-normalise W (each pattern sums to 1 across cell types)
  nc <- function(m) { cs <- colSums(m); cs[cs == 0] <- 1; sweep(m, 2, cs, "/") }
  # Row-normalise H (each pattern sums to 1 across features)
  nr <- function(m) { rs <- rowSums(m); rs[rs == 0] <- 1; sweep(m, 1, rs, "/") }

  W_h_n <- nc(W_h);  W_u_n <- nc(W_u)
  H_h_n <- nr(H_h);  H_u_n <- nr(H_u)

  # Match Unhealthy patterns to Healthy by cosine similarity on shared features
  H_h_c    <- nr(H_h_common)
  H_u_c    <- nr(H_u_common)
  sim      <- cos_sim_mat(H_h_c, H_u_c)   # k × k
  order_u  <- greedy_match(sim)

  W_u_m <- W_u_n[, order_u, drop = FALSE]
  H_u_m <- H_u_n[order_u, , drop = FALSE]

  pat_names <- paste0("P", seq_len(n_patterns))
  colnames(W_h_n) <- colnames(W_u_m) <- pat_names
  rownames(H_h_n) <- rownames(H_u_m) <- pat_names

  # Top features for H display: top n_top from each Healthy pattern,
  # restricted to features present in both conditions
  top_feats <- unique(unlist(lapply(seq_len(n_patterns), function(i)
    names(sort(H_h_n[i, common_feats], decreasing = TRUE))[seq_len(n_top)])))

  # Colour scales shared across both conditions
  w_max   <- max(c(W_h_n, W_u_m), na.rm = TRUE)
  h_max   <- max(c(H_h_n[, top_feats], H_u_m[, top_feats]), na.rm = TRUE)
  col_w   <- spectral_pal(0, w_max)
  col_h   <- spectral_pal(0, h_max)
  col_sim <- colorRamp2(c(0, 0.5, 1), c("white", "#FDB863", "#B2182B"))

  # Row dendrogram for W: cluster cell types on Healthy W, reuse for Unhealthy
  row_dend_w <- as.dendrogram(hclust(dist(W_h_n), method = "ward.D2"))

  pdf(out_pdf, width = 7, height = height)

  # ── Page 1: cosine similarity overview ────────────────────────────────────
  sim_display <- sim[, order_u, drop = FALSE]
  rownames(sim_display) <- paste0("H_", pat_names)
  colnames(sim_display) <- paste0("U_", pat_names)

  draw(
    Heatmap(
      sim_display,
      name            = "Cosine\nsimilarity",
      col             = col_sim,
      cluster_rows    = FALSE,
      cluster_columns = FALSE,
      cell_fun        = function(j, i, x, y, w, h, fill)
        grid.text(sprintf("%.2f", sim_display[i, j]), x, y, gp = gpar(fontsize = 11)),
      row_names_gp    = gpar(fontsize = 10),
      column_names_gp = gpar(fontsize = 10),
      column_title    = sprintf(
        "%s patterns — Healthy vs Unhealthy similarity (matched)",
        tools::toTitleCase(role)),
      column_title_gp = gpar(fontsize = 12, fontface = "bold"),
      heatmap_legend_param = list(direction = "horizontal")
    ),
    heatmap_legend_side = "bottom"
  )

  # ── Page 2: cell-type loadings (W) ────────────────────────────────────────
  base_w <- list(
    cluster_rows    = row_dend_w,
    cluster_columns = FALSE,
    col             = col_w,
    column_names_rot      = 0,
    column_names_centered = TRUE,
    row_labels      = gsub("_", " ", rownames(W_h_n)),
    row_names_gp    = gpar(fontsize = 20),
    column_names_gp = gpar(fontsize = 16),
    heatmap_legend_param = list(
      title         = "Contribution\n(col-norm)",
      title_gp      = gpar(fontsize = 16, fontface = "bold"),
      labels_gp     = gpar(fontsize = 15),
      legend_height = unit(5, "cm"),
      grid_width    = unit(0.6, "cm")
    )
  )

  ht_w <- do.call(Heatmap, c(list(
    W_h_n,
    name         = "Healthy",
    column_title = "Healthy"
  ), base_w)) +
  do.call(Heatmap, c(list(
    W_u_m,
    name              = "Unhealthy",
    column_title      = "Unhealthy",
    show_row_dend     = FALSE,
    show_heatmap_legend = FALSE
  ), base_w))

  draw(ht_w,
       column_title    = sprintf("%s patterns — cell-type loadings",
                                  tools::toTitleCase(role)),
       column_title_gp = gpar(fontsize = 20, fontface = "bold"),
       row_title_gp = gpar(fontsize = 22),
       heatmap_legend_side = "right", merge_legend = TRUE)

  # ── Page 3: feature loadings (H, top features) ────────────────────────────
  # Cluster features on Healthy H, reuse for Unhealthy
  col_dend_h <- as.dendrogram(
    hclust(dist(t(H_h_n[, top_feats, drop = FALSE])), method = "ward.D2"))

  base_h <- list(
    cluster_rows    = FALSE,
    cluster_columns = col_dend_h,
    col             = col_h,
    column_names_rot = 45,
    row_names_gp    = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 7),
    heatmap_legend_param = list(title = "Contribution\n(row-norm)",
                                 direction = "horizontal")
  )

  ht_h <- do.call(Heatmap, c(list(
    H_h_n[, top_feats, drop = FALSE],
    name         = "Healthy",
    column_title = "Healthy",
    row_title    = "Pattern"
  ), base_h)) +
  do.call(Heatmap, c(list(
    H_u_m[, top_feats, drop = FALSE],
    name              = "Unhealthy",
    column_title      = "Unhealthy",
    show_row_names    = FALSE,
    show_heatmap_legend = FALSE
  ), base_h))

  draw(ht_h,
       column_title    = sprintf(
         "%s patterns — %s loadings (top %d per Healthy pattern)",
         tools::toTitleCase(role), feat_lbl, n_top),
       column_title_gp = gpar(fontsize = 12, fontface = "bold"),
       heatmap_legend_side = "bottom", merge_legend = TRUE)

  dev.off()
  message("  Saved: ", basename(out_pdf))
}

# ============================================================================
# 04  Differential ligand-activity heatmap  (outgoing or incoming)
#     Outgoing: sender × ligand  — which ligands does each sender use more?
#     Incoming: receiver × ligand — which ligands activate each receiver more?
#     Both use max_scaled_activity, same feature space as the NMF patterns.
#     Fill = max_scaled_activity_Unhealthy − max_scaled_activity_Healthy.
# ============================================================================
plot_diff_heatmap <- function(mnn, tissue_name, role = c("outgoing", "incoming"),
                               top_n = N_TOP_DIFF, out_pdf) {
  role      <- match.arg(role)
  ct_col    <- if (role == "outgoing") "sender" else "receiver"
  ct_lbl    <- if (role == "outgoing") "sender" else "receiver"

  group_tbl <- mnn$prioritization_tables$group_prioritization_tbl

  # Aggregate max ligand activity per (cell_type, ligand, group), then diff
  diff_tbl <- group_tbl |>
    dplyr::filter(!is.na(max_scaled_activity)) |>
    dplyr::group_by(.data[[ct_col]], ligand, group) |>
    dplyr::summarise(activity = max(max_scaled_activity, na.rm = TRUE),
                     .groups = "drop") |>
    tidyr::pivot_wider(
      id_cols     = c(ct_col, "ligand"),
      names_from  = group,
      values_from = activity,
      values_fill = 0
    ) |>
    dplyr::mutate(diff = Unhealthy - Healthy)

  # Top N ligands per cell type by |diff|
  top_int <- diff_tbl |>
    dplyr::group_by(.data[[ct_col]]) |>
    dplyr::slice_max(order_by = abs(diff), n = top_n, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data[[ct_col]], dplyr::desc(abs(diff)))

  if (nrow(top_int) == 0) {
    message("  No interactions found — skipping diff heatmap [", role, "]")
    return(invisible(NULL))
  }

  # Unique ligands nominated by any cell type (no repeats)
  top_ligands <- unique(top_int$ligand)

  # Full matrix: ALL cell types × ALL nominated ligands
  cts <- sort(unique(diff_tbl[[ct_col]]))
  mat <- diff_tbl |>
    dplyr::filter(ligand %in% top_ligands) |>
    tidyr::pivot_wider(id_cols     = ct_col,
                       names_from  = "ligand",
                       values_from = "diff",
                       values_fill = 0) |>
    (\(w) {
      m <- as.matrix(w[, -1, drop = FALSE])
      rownames(m) <- w[[ct_col]]
      m[, top_ligands, drop = FALSE]   # keep nomination order initially
    })()

  # Left annotation: cell-type colours
  cts_in_mat <- rownames(mat)
  ct_col_vec <- structure(
    ifelse(is.na(CT_COLOURS[cts_in_mat]), "grey70", CT_COLOURS[cts_in_mat]),
    names = cts_in_mat
  )
  row_ann <- rowAnnotation(
    " " = cts_in_mat,
    col = list(" " = ct_col_vec),
    show_annotation_name = FALSE,
    annotation_legend_param = list(" " = list(
      title     = "Cell type",
      title_gp  = gpar(fontsize = 14, fontface = "bold"),
      labels_gp = gpar(fontsize = 13)
    ))
  )

  diff_max <- max(abs(mat), na.rm = TRUE)
  col_fun  <- colorRamp2(
    c(-diff_max, 0, diff_max),
    c("#4DBBD5", "white", "#E64B35")
  )

  ht <- Heatmap(
    mat,
    name                  = "Activity diff\n(U − H)",
    col                   = col_fun,
    left_annotation       = row_ann,
    cluster_rows          = TRUE,
    clustering_method_rows    = "ward.D2",
    cluster_columns       = TRUE,
    clustering_method_columns = "ward.D2",
    show_row_names        = TRUE,
    show_column_names     = TRUE,
    row_names_gp          = gpar(fontsize = 18),
    column_names_gp       = gpar(fontsize = 18),
    column_names_rot      = 45,
    rect_gp               = gpar(col = "grey85", lwd = 0.5),
    heatmap_legend_param  = list(
      title         = "Ligand activity diff\n(Ob+T2D − Ob)",
      title_gp      = gpar(fontsize = 14, fontface = "bold"),
      labels_gp     = gpar(fontsize = 13),
      legend_height = unit(5, "cm"),
      grid_width    = unit(0.6, "cm")
    )
  )

  width  <- max(8, length(top_ligands) * 0.5 + 4)
  height <- max(4, nrow(mat) * 0.5 + 2)

  pdf(out_pdf, width = width, height = height)
  draw(ht,
       column_title    = paste0(tools::toTitleCase(tissue_name),
                                 " — Top ", top_n, " differential ligands per ",
                                 ct_lbl, " (max scaled activity, ", role, ")"),
       column_title_gp = gpar(fontsize = 12, fontface = "bold"),
       heatmap_legend_side = "right", merge_legend = TRUE)
  dev.off()
  message("  Saved: ", basename(out_pdf))
}

# ============================================================================
# Main: iterate over tissues
# ============================================================================
for (tissue in names(MNN_PATHS)) {
  cat("\n==== Tissue:", tissue, "====\n")
  mnn <- readRDS(MNN_PATHS[[tissue]])
  tbl <- get_prio_tbl(mnn)
  cat("  Interactions kept (score >=", MIN_SCORE, "):", nrow(tbl), "\n")

  tdir <- file.path(OUTPUT_DIR, tissue)
  dir.create(tdir, recursive = TRUE, showWarnings = FALSE)

  # 01: Sender / receiver scatter
  p_sr <- plot_sr_scatter(
    tbl, title = paste0(tools::toTitleCase(tissue),
                         " — Sender/Receiver signaling strength"))
  ggsave(file.path(tdir, "01_sender_receiver_scatter.pdf"),
         p_sr, width = 12, height = 10)
  message("  Saved: 01_sender_receiver_scatter.pdf")

  # 02: Communication topology heatmap
  pdf(file.path(tdir, "02_comm_heatmap.pdf"), width = 14, height = 6)
  plot_comm_heatmap(
    tbl, title = paste0(tools::toTitleCase(tissue),
                         " — Sender × Receiver communication"))
  dev.off()
  message("  Saved: 02_comm_heatmap.pdf")

  # 03: NMF pattern comparison (outgoing + incoming)
  for (role in c("outgoing", "incoming")) {
    out_pdf <- file.path(tdir, sprintf("03_nmf_%s_patterns.pdf", role))
    plot_nmf_patterns(tbl, role = role, out_pdf = out_pdf)
  }

  # 04: Top differential ligands per sender / receiver
  for (role in c("outgoing", "incoming")) {
    plot_diff_heatmap(mnn, tissue_name = tissue, role = role,
                      out_pdf = file.path(tdir,
                                          sprintf("04_diff_%s_heatmap.pdf", role)))
  }
}

cat("\nDone. Figures saved to:", OUTPUT_DIR, "\n")
