# mofa_analysis.R
# ============================================================
# MOFA+ Integration: Flow Cytometry (3 tissues) + Cytokines
# ============================================================
# Input:
#   data/flow_cytometry_perct_nonstim.csv  — all features as % (SubQ/Blood/Omentum)
#   data/multiplex_cytokines.csv           — 11 plasma cytokines
#   data/clinical_data.csv                 — metadata + diagnoses
#
# MOFA+ setup:
#   4 views: SubQ_flow, Blood_flow, Omentum_flow, Cytokines
#   2 groups (multi-group): T2D vs Non-T2D
#   Missing samples handled natively by MOFA+
#
# Output: results/mofa/
# ============================================================

suppressPackageStartupMessages({
  library(MOFA2)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(corrplot)
  library(RColorBrewer)
})

set.seed(42)
OUT_DIR <- "results/mofa_with_stim_INS_variance"
setwd("~/Documents/projs/Wills_analysis")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 1. LOAD DATA ────────────────────────────────────────────

flow_raw  <- read.csv("data/flow_cytometry_perct_nonstim.csv",
                      check.names = FALSE, fileEncoding = "UTF-8-BOM")
names(flow_raw) <- make.unique(names(flow_raw))
cyto_raw  <- read.csv("data/multiplex_cytokines.csv",
                      check.names = FALSE, fileEncoding = "UTF-8-BOM")
names(cyto_raw) <- make.unique(names(cyto_raw))
clinical  <- read.csv("data/clinical_data.csv", check.names = FALSE)
names(clinical) <- make.unique(names(clinical))

# Clean patient IDs
flow_raw$Patient  <- trimws(flow_raw$Patient)
cyto_raw$Patient  <- trimws(cyto_raw$Patient)
clinical$Patient  <- trimws(clinical$Patient)
rownames(clinical) <- clinical$Patient

# ── 2. PREPROCESS FLOW CYTOMETRY ────────────────────────────
# All features are percentages (0–100). Apply arcsine-sqrt transformation,
# which is the standard variance-stabilising transform for compositional data.
# Formula: asin(sqrt(x / 100))  →  maps [0,100] → [0, π/2]

asin_sqrt <- function(x) asin(sqrt(pmax(0, pmin(100, x)) / 100))

subq_cols    <- grep("^SubQ",    names(flow_raw), value = TRUE)
blood_cols   <- grep("^Blood",   names(flow_raw), value = TRUE)
omentum_cols <- grep("^Omentum", names(flow_raw), value = TRUE)

make_view <- function(df, patient_col, feature_cols, transform_fn = identity) {
  # Returns a features × samples matrix (MOFA format), with NaN for missing patients
  m <- df[, c(patient_col, feature_cols), drop = FALSE]
  rownames(m) <- m[[patient_col]]
  m <- m[, feature_cols, drop = FALSE]
  m[m == ""] <- NA
  m <- apply(m, 2, as.numeric)
  rownames(m) <- df[[patient_col]]
  m <- m[rowSums(is.na(m)) != ncol(m),]
  # colnames(m) <- iconv(colnames(m), from = "UTF-8", to = "ASCII", sub = "?")
  # Apply transformation feature-wise (NA-aware)
  m_t <- apply(m, 2, function(col) ifelse(is.na(col), NA, transform_fn(col)))
  return(t(m_t))  # features × samples
}

subq_mat    <- make_view(flow_raw, "Patient", subq_cols,    asin_sqrt)
blood_mat   <- make_view(flow_raw, "Patient", blood_cols,   asin_sqrt)
omentum_mat <- make_view(flow_raw, "Patient", omentum_cols, asin_sqrt)

# ── 3. PREPROCESS CYTOKINES ─────────────────────────────────
# Cytokines are right-skewed plasma concentrations → log1p transform

cyto_cols <- setdiff(names(cyto_raw), "Patient")
cyto_mat  <- make_view(cyto_raw, "Patient", cyto_cols,
                       transform_fn = function(x) log1p(x))

# ── 3b. PRE-MOFA QC: NORMALITY & VARIANCE CHECKS ────────────
# Runs BEFORE feature filtering or MOFA training.
# Outputs → OUT_DIR/00_QC/
#
# Answers two questions:
#   (a) Did the transforms achieve approximate normality?
#       (MOFA assumes Gaussian likelihood)
#   (b) Are there near-zero-variance features worth dropping to reduce noise?
# ─────────────────────────────────────────────────────────────

QC_DIR <- file.path(OUT_DIR, "00_QC")
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

# Helper: raw (untransformed) samples × features matrix
raw_view <- function(df, patient_col, feature_cols) {
  m <- df[, c(patient_col, feature_cols), drop = FALSE]
  rownames(m) <- m[[patient_col]]
  m <- m[, feature_cols, drop = FALSE]
  m[m == ""] <- NA
  m <- apply(m, 2, as.numeric)
  rownames(m) <- df[[patient_col]]
  m <- m[rowSums(is.na(m)) != ncol(m), , drop = FALSE]
  colnames(m) <- iconv(colnames(m), from = "UTF-8", to = "ASCII", sub = "?")
  m  # samples × features, raw values
}

# Collect raw + transformed matrices for each view
# trf matrices are transposed back to samples × features for QC convenience
views_qc <- list(
  SubQ = list(
    raw = raw_view(flow_raw, "Patient", subq_cols),
    trf = t(subq_mat),
    raw_label = "Raw %", trf_label = "asin(sqrt(x/100))"
  ),
  Blood = list(
    raw = raw_view(flow_raw, "Patient", blood_cols),
    trf = t(blood_mat),
    raw_label = "Raw %", trf_label = "asin(sqrt(x/100))"
  ),
  Omentum = list(
    raw = raw_view(flow_raw, "Patient", omentum_cols),
    trf = t(omentum_mat),
    raw_label = "Raw %", trf_label = "asin(sqrt(x/100))"
  ),
  Cytokines = list(
    raw = raw_view(cyto_raw, "Patient", cyto_cols),
    trf = t(cyto_mat),
    raw_label = "Raw concentration", trf_label = "log1p"
  )
)

view_cols <- c(SubQ = "#3182BD", Blood = "#E6550D",
               Omentum = "#31A354", Cytokines = "#756BB1")

# ── 3b-i. POOLED DENSITY: raw vs transformed ─────────────────
# Pooled across all features in the view — gives a gestalt of the
# shift in distribution shape achieved by each transform.

message("QC: plotting pooled raw vs transformed densities...")

pdf(file.path(QC_DIR, "01_pooled_density_raw_vs_transformed.pdf"),
    width = 14, height = 6)

for (vname in names(views_qc)) {
  v <- views_qc[[vname]]
  raw_vals <- as.vector(v$raw); raw_vals <- raw_vals[!is.na(raw_vals)]
  trf_vals <- as.vector(v$trf); trf_vals <- trf_vals[!is.na(trf_vals)]

  par(mfrow = c(1, 2), mar = c(4, 4, 3.5, 1))

  # Raw
  d_raw <- density(raw_vals)
  plot(d_raw, col = view_cols[vname], lwd = 2,
       main = paste0(vname, "  —  ", v$raw_label, "\n(pooled across all features)"),
       xlab = v$raw_label)
  polygon(d_raw, col = adjustcolor(view_cols[vname], 0.2), border = NA)
  rug(sample(raw_vals, min(300, length(raw_vals))),
      col = adjustcolor(view_cols[vname], 0.35), ticksize = 0.03)

  # Transformed
  d_trf <- density(trf_vals)
  plot(d_trf, col = view_cols[vname], lwd = 2,
       main = paste0(vname, "  —  ", v$trf_label, "\n(pooled across all features)"),
       xlab = v$trf_label)
  polygon(d_trf, col = adjustcolor(view_cols[vname], 0.2), border = NA)
  rug(sample(trf_vals, min(300, length(trf_vals))),
      col = adjustcolor(view_cols[vname], 0.35), ticksize = 0.03)
  # Overlay reference normal
  curve(dnorm(x, mean(trf_vals), sd(trf_vals)), add = TRUE,
        col = "grey30", lwd = 1.5, lty = 2)
  legend("topright", c("observed", "normal ref."),
         col = c(view_cols[vname], "grey30"), lwd = c(2, 1.5), lty = c(1, 2),
         bty = "n", cex = 0.8)
}

dev.off()

# ── 3b-ii. PER-FEATURE HISTOGRAMS (sample of 16) ─────────────
# Spot-checks individual feature distributions after transformation.
# Overlaid normal curve highlights remaining skew or bimodality.

message("QC: plotting per-feature histograms (sample of 16)...")

pdf(file.path(QC_DIR, "02_perfeature_histograms_transformed.pdf"),
    width = 16, height = 12)

for (vname in names(views_qc)) {
  v     <- views_qc[[vname]]
  feats <- colnames(v$trf)
  set.seed(42)
  sel   <- if (length(feats) <= 16) feats else sample(feats, 16)

  par(mfrow = c(4, 4), mar = c(3, 3, 2.5, 0.5), oma = c(0, 0, 3, 0))
  for (f in sel) {
    vals <- v$trf[, f]; vals <- vals[!is.na(vals)]
    if (length(vals) < 3) { plot.new(); next }
    h <- hist(vals, breaks = 12, plot = FALSE)
    hist(vals, breaks = 12,
         main = substr(gsub(paste0("^", vname, "_?"), "", f), 1, 32),
         xlab = "", col = adjustcolor(view_cols[vname], 0.5),
         border = "white", cex.main = 0.72, cex.axis = 0.8)
    # Normal curve scaled to histogram
    bin_w <- diff(h$breaks[1:2])
    curve(dnorm(x, mean(vals), sd(vals)) * length(vals) * bin_w,
          add = TRUE, col = "grey25", lwd = 1.5)
  }
  mtext(paste0(vname, ": ", v$trf_label, "  — random sample of ", length(sel), " features"),
        outer = TRUE, cex = 1, font = 2)
}

dev.off()

# ── 3b-iii. QQ PLOTS (sample of 16 features) ─────────────────
# Points on the diagonal = normally distributed.
# Curved tails = skew; S-shape = heavy tails or bimodality.

message("QC: plotting QQ plots (sample of 16 per view)...")

pdf(file.path(QC_DIR, "03_QQ_plots_transformed.pdf"), width = 16, height = 12)

for (vname in names(views_qc)) {
  v     <- views_qc[[vname]]
  feats <- colnames(v$trf)
  set.seed(42)
  sel   <- if (length(feats) <= 16) feats else sample(feats, 16)

  par(mfrow = c(4, 4), mar = c(3, 3, 2.5, 0.5), oma = c(0, 0, 3, 0))
  for (f in sel) {
    vals <- v$trf[, f]; vals <- vals[!is.na(vals)]
    if (length(vals) < 3) { plot.new(); next }
    qqnorm(vals,
           main = substr(gsub(paste0("^", vname, "_?"), "", f), 1, 32),
           pch = 16, cex = 0.7, col = view_cols[vname],
           cex.main = 0.72, cex.axis = 0.8)
    qqline(vals, col = "grey30", lwd = 1.5)
  }
  mtext(paste0(vname, ": QQ plots after ", v$trf_label),
        outer = TRUE, cex = 1, font = 2)
}

dev.off()

# ── 3b-iv. SHAPIRO-WILK TEST (all features) ──────────────────
# SW p-value > 0.05 (after BH FDR correction) = fail to reject normality.
# With n~30-40 SW is sensitive — expect some rejections even for roughly
# normal data. Focus on the overall % and p-value distribution shape
# (right-skewed p-values = most features are approximately normal).

message("QC: running Shapiro-Wilk normality tests...")

sw_results <- lapply(names(views_qc), function(vname) {
  trf_mat <- views_qc[[vname]]$trf   # samples × features
  pvals <- apply(trf_mat, 2, function(col) {
    col <- col[!is.na(col)]
    if (length(col) < 5) return(NA_real_)
    tryCatch(shapiro.test(col)$p.value, error = function(e) NA_real_)
  })
  data.frame(view = vname, feature = names(pvals), sw_pval = pvals,
             stringsAsFactors = FALSE)
})
sw_df <- do.call(rbind, sw_results)

sw_df$sw_padj <- NA_real_
for (vname in unique(sw_df$view)) {
  idx <- sw_df$view == vname & !is.na(sw_df$sw_pval)
  sw_df$sw_padj[idx] <- p.adjust(sw_df$sw_pval[idx], method = "BH")
}
sw_df$normal <- sw_df$sw_padj > 0.05   # TRUE = does not reject normality

sw_summary <- sw_df %>%
  group_by(view) %>%
  summarise(
    n_features     = n(),
    n_normal       = sum(normal, na.rm = TRUE),
    pct_normal     = round(100 * mean(normal, na.rm = TRUE), 1),
    median_sw_pval = round(median(sw_pval, na.rm = TRUE), 4),
    .groups = "drop"
  )

message("  Shapiro-Wilk summary (BH FDR<0.05 = reject normality):")
print(as.data.frame(sw_summary))
write.csv(sw_df,      file.path(QC_DIR, "04_shapiro_wilk_perfeature.csv"), row.names = FALSE)
write.csv(sw_summary, file.path(QC_DIR, "04_shapiro_wilk_summary.csv"),    row.names = FALSE)

pdf(file.path(QC_DIR, "04_shapiro_wilk_plots.pdf"), width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5))

# Panel A: p-value density per view
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 5),
     xlab = "Shapiro-Wilk p-value (after transform)",
     ylab = "Density",
     main = "SW p-value distribution per view\n(right-skewed = approximately normal)")
for (vname in names(views_qc)) {
  pv <- sw_df$sw_pval[sw_df$view == vname & !is.na(sw_df$sw_pval)]
  if (length(pv) > 1) {
    d <- density(pv, from = 0, to = 1, bw = 0.08)
    lines(d, col = view_cols[vname], lwd = 2)
  }
}
abline(v = 0.05, lty = 2, col = "grey40", lwd = 1)
legend("topright", names(view_cols), col = view_cols, lwd = 2, bty = "n", cex = 0.85)

# Panel B: % features passing normality per view
bp <- barplot(sw_summary$pct_normal,
              names.arg = sw_summary$view,
              col       = view_cols[sw_summary$view],
              ylim      = c(0, 115),
              ylab      = "% features not rejecting normality (BH-adjusted)",
              main      = "% approximately normal features\nafter transformation",
              border    = "white")
abline(h = 50, lty = 2, col = "grey50")
text(bp, sw_summary$pct_normal + 4,
     labels = paste0(sw_summary$pct_normal, "%"),
     cex = 0.95, font = 2)

dev.off()

# ── 3b-v. FEATURE VARIANCE PLOTS ─────────────────────────────
# Three panels per view:
#   A) Ranked variance elbow plot — look for a natural elbow
#   B) Variance histogram — bimodal shape → two populations present
#   C) Cumulative variance — how many features cover 80/90% of variance?
#
# Two threshold lines are drawn:
#   - Median (orange dashed): keeps top 50% most variable
#   - 75th percentile (green dotted): keeps top 25% most variable
#
# FILTER_LOW_VARIANCE below lets you act on these thresholds.

message("QC: computing and plotting feature variances...")

var_results <- lapply(names(views_qc), function(vname) {
  trf_mat <- views_qc[[vname]]$trf
  feat_var <- apply(trf_mat, 2, function(col) {
    col <- col[!is.na(col)]
    if (length(col) < 3) return(NA_real_)
    var(col)
  })
  data.frame(view = vname, feature = names(feat_var), variance = feat_var,
             stringsAsFactors = FALSE)
})

var_df <- do.call(rbind, var_results) %>%
  dplyr::filter(!is.na(variance)) %>%
  group_by(view) %>%
  arrange(dplyr::desc(variance), .by_group = TRUE) %>%
  dplyr::mutate(
    rank         = row_number(),
    cum_var_frac = cumsum(variance) / sum(variance)
  ) %>%
  dplyr::ungroup()

var_thresholds <- var_df %>%
  group_by(view) %>%
  summarise(
    n_features     = n(),
    var_min        = round(min(variance), 5),
    var_median     = round(median(variance), 5),
    var_75pct      = round(quantile(variance, 0.75), 5),
    n_above_median = sum(variance > median(variance)),
    n_above_75pct  = sum(variance > quantile(variance, 0.75)),
    k_80pct_cumvar = min(which(cumsum(sort(variance, decreasing=TRUE)) /
                                 sum(variance) >= 0.80)),
    k_90pct_cumvar = min(which(cumsum(sort(variance, decreasing=TRUE)) /
                                 sum(variance) >= 0.90)),
    .groups = "drop"
  )

message("  Variance thresholds per view:")
print(as.data.frame(var_thresholds))
write.csv(var_df,         file.path(QC_DIR, "05_feature_variances.csv"),    row.names = FALSE)
write.csv(var_thresholds, file.path(QC_DIR, "05_variance_thresholds.csv"), row.names = FALSE)

pdf(file.path(QC_DIR, "05_variance_plots.pdf"),
    width = 15, height = 1.5 * length(views_qc))

for (vname in names(views_qc)) {
  vd  <- filter(var_df, view == vname)
  thr <- filter(var_thresholds, view == vname)
  col <- view_cols[vname]

  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 1.5), oma = c(0, 0, 3.5, 0))

  # A: Ranked variance (elbow plot)
  plot(vd$rank, vd$variance, type = "l", lwd = 1.8, col = col,
       xlab = "Feature rank (1 = most variable)",
       ylab = "Variance (after transform)",
       main = "A: Ranked variance\n(look for an elbow)")
  abline(h = thr$var_median, lty = 2, col = "#E6550D", lwd = 1.3)
  abline(h = thr$var_75pct,  lty = 3, col = "#31A354", lwd = 1.3)
  abline(v = thr$n_above_median, lty = 2, col = "#E6550D", lwd = 1)
  abline(v = thr$n_above_75pct,  lty = 3, col = "#31A354", lwd = 1)
  legend("topright", bty = "n", cex = 0.78,
         legend = c(sprintf("Median — keep n=%d", thr$n_above_median),
                    sprintf("75th pct — keep n=%d", thr$n_above_75pct)),
         lty = c(2, 3), col = c("#E6550D", "#31A354"), lwd = 1.5)

  # B: Variance histogram
  hist(vd$variance, breaks = 30,
       col = adjustcolor(col, 0.5), border = "white",
       xlab = "Variance (after transform)",
       main = "B: Variance distribution\n(bimodal = informative vs flat features)")
  abline(v = thr$var_median, lty = 2, col = "#E6550D", lwd = 1.5)
  abline(v = thr$var_75pct,  lty = 3, col = "#31A354", lwd = 1.5)
  rug(vd$variance, col = adjustcolor(col, 0.4), ticksize = 0.03)

  # C: Cumulative variance
  plot(vd$rank, vd$cum_var_frac * 100, type = "l", lwd = 1.8, col = col,
       xlab = "Top-k features (ranked by variance)",
       ylab = "Cumulative % of total variance",
       main = "C: Cumulative variance\n(how many features capture 80 / 90%?)",
       ylim = c(0, 100))
  abline(h = c(80, 90), lty = c(2, 3), col = "grey50", lwd = 1)
  abline(v = thr$k_80pct_cumvar, lty = 2, col = "grey50", lwd = 1)
  abline(v = thr$k_90pct_cumvar, lty = 3, col = "grey50", lwd = 1)
  text(thr$k_80pct_cumvar, 78,
       paste0("80%: k=", thr$k_80pct_cumvar), adj = c(-0.05, 1), cex = 0.8, col = "grey30")
  text(thr$k_90pct_cumvar, 88,
       paste0("90%: k=", thr$k_90pct_cumvar), adj = c(-0.05, 1), cex = 0.8, col = "grey30")

  mtext(paste0(vname, ": variance-based feature assessment"),
        outer = TRUE, cex = 1.1, font = 2)
}

dev.off()

# ── 3b-vi. OPTIONAL VARIANCE FILTERING ───────────────────────
# Set FILTER_LOW_VARIANCE <- TRUE to drop features below a per-view threshold.
# Default threshold = median variance (keeps top 50% most variable features).
#
# Note: with scale_views=TRUE and ARD priors, MOFA+ naturally down-weights
# low-variance features, so filtering is optional. Apply it if:
#   - The elbow plot (panel A) shows a clear break, OR
#   - Many features are near-zero-variance (seen in panel B histogram), OR
#   - You want faster training / sparser, more interpretable weights.

FILTER_LOW_VARIANCE <- TRUE   # ← change to TRUE to activate

if (FILTER_LOW_VARIANCE) {
  message("Applying variance filter (threshold = median variance per view)...")

  filter_by_var <- function(mat, min_var) {
    # mat is features × samples (MOFA orientation)
    feat_var <- apply(mat, 1, function(r) {
      r <- r[!is.na(r)]; if (length(r) < 3) 0 else var(r)
    })
    keep <- feat_var >= min_var
    message(sprintf("  Keeping %d / %d features (var >= %.5f)",
                    sum(keep), length(keep), min_var))
    mat[keep, , drop = FALSE]
  }

  subq_mat    <- filter_by_var(subq_mat,
                   filter(var_thresholds, view == "SubQ")$var_median)
  blood_mat   <- filter_by_var(blood_mat,
                   filter(var_thresholds, view == "Blood")$var_median)
  omentum_mat <- filter_by_var(omentum_mat,
                   filter(var_thresholds, view == "Omentum")$var_median)
  cyto_mat <- filter_by_var(cyto_mat,
                               filter(var_thresholds, view == "Cytokines")$var_median)
}

message("QC complete. Outputs in: ", QC_DIR)
message("  01_pooled_density_raw_vs_transformed.pdf")
message("  02_perfeature_histograms_transformed.pdf")
message("  03_QQ_plots_transformed.pdf")
message("  04_shapiro_wilk_*.pdf / *.csv")
message("  05_variance_plots.pdf / *.csv")

# ── 3b-vii. NON-NORMALITY DIAGNOSIS & INT TRANSFORMATION ─────
#
# Context: ~50% of flow features failing SW is common. There are two
# distinct causes requiring different responses:
#
#   (A) ZERO-INFLATION: cell subset absent in many patients → spike at 0.
#       asin-sqrt of 0% = 0, leaving a point mass SW correctly rejects.
#       These features mislead the MOFA Gaussian likelihood.
#       → Remove features with > ZERO_FRAC_THRESHOLD fraction of exact zeros.
#
#   (B) BIMODALITY: feature genuinely splits patients into two groups.
#       SW rejects this too, but MOFA *will capture it as a factor* even
#       under the Gaussian assumption — this is signal, not noise.
#       → Keep these features; INT will spread them into a normal shape.
#
# Solution applied here:
#   Step 1 — drop zero-inflated features (point A)
#   Step 2 — apply Rank-based Inverse Normal Transformation (INT) to
#             all remaining features. INT maps values to normal quantiles
#             via their ranks. It is guaranteed to produce N(0,1) output
#             regardless of input distribution, and preserves relative order.
#
# INT formula (Blom offset, standard in omics):
#   INT(x) = Φ⁻¹( (rank(x) - 0.5) / n_observed )
#
# USE_INT <- TRUE/FALSE controls whether INT or asin-sqrt is used in MOFA.
# Inspect 06_INT_normality_comparison.pdf before deciding.

ZERO_FRAC_THRESHOLD <- 0.40   # drop features with >40% exact zeros (raw %)
USE_INT             <- TRUE    # ← set FALSE to keep the asin-sqrt transform

# Rank-based inverse normal transformation (applied per feature, NA-aware)
int_transform <- function(x) {
  r <- rank(x, na.last = "keep", ties.method = "average")
  n <- sum(!is.na(x))
  qnorm((r - 0.5) / n)
}

# -- Step 1: identify zero-inflated features (on RAW percentages) -----------
message("QC-INT: diagnosing zero-inflation in raw flow features...")

zero_inf_results <- lapply(
  c("SubQ", "Blood", "Omentum"),
  function(vname) {
    raw_mat <- views_qc[[vname]]$raw   # samples × features, raw %
    zero_frac <- apply(raw_mat, 2, function(col) {
      col <- col[!is.na(col)]
      if (length(col) == 0) return(1)
      mean(col == 0)
    })
    data.frame(view = vname, feature = names(zero_frac),
               zero_frac = zero_frac, stringsAsFactors = FALSE)
  }
)
zero_inf_df <- do.call(rbind, zero_inf_results)
zero_inf_df$zero_inflated <- zero_inf_df$zero_frac > ZERO_FRAC_THRESHOLD

zi_summary <- zero_inf_df %>%
  group_by(view) %>%
  summarise(
    n_features       = n(),
    n_zero_inflated  = sum(zero_inflated),
    pct_zero_inflated = round(100 * mean(zero_inflated), 1),
    .groups = "drop"
  )
message("  Zero-inflation summary (threshold >", ZERO_FRAC_THRESHOLD * 100, "% zeros):")
print(as.data.frame(zi_summary))
write.csv(zero_inf_df, file.path(QC_DIR, "06_zero_inflation.csv"), row.names = FALSE)

# -- Step 2: apply INT and rebuild MOFA matrices ----------------------------
# We apply INT to the (already row-subsetted) raw matrices.
# Zero-inflated features are either dropped (if USE_INT) or kept as-is.
# Variance filter threshold on asin-sqrt scale (applied before INT)
# "bottom quartile" drops the 25% least variable features per view.
# Set to 0 to disable, or use a fixed value from the QC elbow plots.
VAR_QUANTILE_THRESHOLD <- 0.25   # drop bottom 25% by asin-sqrt variance

apply_int_to_view <- function(raw_mat, zi_df, drop_zi = TRUE,
                              var_qthr = VAR_QUANTILE_THRESHOLD) {
  # Step 1: ZI filter
  if (drop_zi) {
    keep_feats <- zi_df$feature[!zi_df$zero_inflated]
    raw_mat    <- raw_mat[, intersect(keep_feats, colnames(raw_mat)), drop = FALSE]
  }
  
  # Step 2: variance filter on asin-sqrt scale
  if (var_qthr > 0) {
    asinsqrt_mat <- apply(raw_mat, 2,
                          function(x) asin(sqrt(pmax(0, pmin(100, x)) / 100)))
    feat_var  <- apply(asinsqrt_mat, 2, var, na.rm = TRUE)
    var_cutoff <- quantile(feat_var, var_qthr, na.rm = TRUE)
    keep      <- feat_var >= var_cutoff
    message(sprintf("    Variance filter (bottom %.0f%%): dropping %d / %d features",
                    var_qthr * 100, sum(!keep), length(keep)))
    raw_mat   <- raw_mat[, keep, drop = FALSE]
  }
  
  # Step 3: INT on the filtered raw values
  m_int <- apply(raw_mat, 2, int_transform)
  rownames(m_int) <- rownames(raw_mat)
  t(m_int)   # features × samples (MOFA orientation)
}

if (USE_INT) {
  message("Applying INT transformation (dropping zero-inflated features)...")

  subq_zi    <- filter(zero_inf_df, view == "SubQ")
  blood_zi   <- filter(zero_inf_df, view == "Blood")
  omentum_zi <- filter(zero_inf_df, view == "Omentum")

  subq_mat    <- apply_int_to_view(views_qc$SubQ$raw,    subq_zi,    drop_zi = TRUE)
  blood_mat   <- apply_int_to_view(views_qc$Blood$raw,   blood_zi,   drop_zi = TRUE)
  omentum_mat <- apply_int_to_view(views_qc$Omentum$raw, omentum_zi, drop_zi = TRUE)

  # Cytokines: no zeros expected (they're concentrations), apply INT without ZI filter
  cyto_raw_mat <- views_qc$Cytokines$raw
  cyto_mat     <- t(apply(cyto_raw_mat, 2, int_transform))
  rownames(cyto_mat) <- colnames(cyto_raw_mat)

  message(sprintf("  SubQ:    %d features after ZI filter + INT",    nrow(subq_mat)))
  message(sprintf("  Blood:   %d features after ZI filter + INT",    nrow(blood_mat)))
  message(sprintf("  Omentum: %d features after ZI filter + INT",    nrow(omentum_mat)))
  message(sprintf("  Cytokines: %d features after INT",              nrow(cyto_mat)))
} else {
  message("Using asin-sqrt / log1p transforms (USE_INT = FALSE)")
}

# -- Step 3: re-run SW test on INT-transformed features & compare -----------
message("QC-INT: comparing normality before and after INT...")

sw_after_int <- lapply(
  list(SubQ = subq_mat, Blood = blood_mat, Omentum = omentum_mat),
  function(mat) {
    apply(mat, 1, function(row) {   # mat is features × samples; row = one feature
      vals <- row[!is.na(row)]
      if (length(vals) < 5) return(NA_real_)
      tryCatch(shapiro.test(vals)$p.value, error = function(e) NA_real_)
    })
  }
)

sw_int_df <- bind_rows(lapply(names(sw_after_int), function(vname) {
  pvals <- sw_after_int[[vname]]
  data.frame(view = vname, feature = names(pvals), sw_pval_int = pvals,
             stringsAsFactors = FALSE)
}))
sw_int_df$sw_padj_int <- NA_real_
for (vname in unique(sw_int_df$view)) {
  idx <- sw_int_df$view == vname & !is.na(sw_int_df$sw_pval_int)
  sw_int_df$sw_padj_int[idx] <- p.adjust(sw_int_df$sw_pval_int[idx], method = "BH")
}
sw_int_df$normal_int <- sw_int_df$sw_padj_int > 0.05

# Join with original SW results for comparison
sw_compare <- sw_df %>%
  select(view, feature, sw_pval, sw_padj, normal) %>%
  left_join(sw_int_df, by = c("view", "feature")) %>%
  left_join(select(zero_inf_df, view, feature, zero_frac, zero_inflated),
            by = c("view", "feature"))

write.csv(sw_compare, file.path(QC_DIR, "06_normality_comparison.csv"), row.names = FALSE)

sw_compare_summary <- sw_compare %>%
  filter(!zero_inflated | !USE_INT) %>%   # only features kept in the model
  group_by(view) %>%
  summarise(
    n_features_original  = n(),
    pct_normal_asinsqrt  = round(100 * mean(normal, na.rm = TRUE), 1),
    pct_normal_INT       = round(100 * mean(normal_int, na.rm = TRUE), 1),
    .groups = "drop"
  )
message("  Normality comparison (% features passing SW BH<0.05):")
print(as.data.frame(sw_compare_summary))

# -- Step 4: visualise the comparison ---------------------------------------
pdf(file.path(QC_DIR, "06_INT_normality_comparison.pdf"), width = 14, height = 10)

for (vname in c("SubQ", "Blood", "Omentum")) {
  vc   <- filter(sw_compare, view == vname)
  col  <- view_cols[vname]

  par(mfrow = c(2, 3), mar = c(4, 4, 3.5, 1.5), oma = c(0, 0, 4, 0))

  # Original SW p-value distribution
  pv_orig <- vc$sw_pval[!is.na(vc$sw_pval)]
  hist(pv_orig, breaks = 20, col = adjustcolor(col, 0.5), border = "white",
       main = "SW p-values: asin-sqrt",
       xlab = "p-value", xlim = c(0, 1))
  abline(v = 0.05, lty = 2, col = "grey40")
  mtext(sprintf("%.0f%% normal", 100 * mean(pv_orig > 0.05)), line = -1.5,
        cex = 0.85, col = "grey30")

  # INT SW p-value distribution
  pv_int <- vc$sw_pval_int[!is.na(vc$sw_pval_int) & (!vc$zero_inflated | !USE_INT)]
  hist(pv_int, breaks = 20, col = adjustcolor("#31A354", 0.5), border = "white",
       main = "SW p-values: INT (after ZI filter)",
       xlab = "p-value", xlim = c(0, 1))
  abline(v = 0.05, lty = 2, col = "grey40")
  mtext(sprintf("%.0f%% normal", 100 * mean(pv_int > 0.05, na.rm = TRUE)),
        line = -1.5, cex = 0.85, col = "grey30")

  # Zero-inflation: fraction of zeros per feature
  zi_v <- filter(zero_inf_df, view == vname)
  hist(zi_v$zero_frac * 100, breaks = 20,
       col = adjustcolor("#E6550D", 0.5), border = "white",
       main = sprintf("Zero fraction per feature\n(%d dropped, threshold=%.0f%%)",
                      sum(zi_v$zero_inflated), ZERO_FRAC_THRESHOLD * 100),
       xlab = "% patients with exact 0%")
  abline(v = ZERO_FRAC_THRESHOLD * 100, lty = 2, col = "#E6550D", lwd = 1.5)

  # QQ plots: 3 examples — pick worst SW p-val features
  worst_feats <- vc %>%
    filter(!zero_inflated | !USE_INT) %>%
    arrange(sw_pval) %>%
    slice_head(n = 3) %>%
    pull(feature)

  for (f in worst_feats) {
    # side-by-side: original vs INT
    raw_vals  <- views_qc[[vname]]$trf[, f]   # asin-sqrt
    int_vals  <- if (USE_INT && f %in% rownames(subq_mat) ||
                                f %in% rownames(blood_mat) ||
                                f %in% rownames(omentum_mat)) {
      mat_name <- switch(vname, SubQ = subq_mat, Blood = blood_mat, Omentum = omentum_mat)
      if (f %in% rownames(mat_name)) mat_name[f, ] else NULL
    } else NULL

    fname_short <- substr(gsub(paste0("^", vname, "_?"), "", f), 1, 28)

    qqnorm(raw_vals[!is.na(raw_vals)],
           main = paste0(fname_short, "\nasin-sqrt"),
           pch = 16, cex = 0.7, col = col,
           cex.main = 0.7, cex.axis = 0.8)
    qqline(raw_vals[!is.na(raw_vals)], col = "grey30", lwd = 1.5)
  }

  mtext(paste0(vname, ": normality before vs after INT + ZI filter"),
        outer = TRUE, cex = 1.1, font = 2)
}

dev.off()

message("QC-INT complete.")
message("  06_zero_inflation.csv          — fraction of zeros per feature")
message("  06_normality_comparison.csv    — SW p-values before and after INT")
message("  06_INT_normality_comparison.pdf")
message(if (USE_INT) "  → INT is ACTIVE: MOFA will use INT-transformed matrices."
        else         "  → asin-sqrt is ACTIVE (USE_INT = FALSE).")

# ── 3c. LOAD & PROCESS STIMULATED FLOW DATA ──────────────────
# Stimulated panel: cytokine production (IFNγ, TNF, IL4, IL10, IL17A/F)
# measured per cell type (MAIT, iNKT, CD4, CD8, γδ T, Lin- T) after
# ex vivo PMA/ionomycin stimulation.
#
# These are SEPARATE views from the nonstim data — they measure functional
# capacity, not abundance. Do NOT merge with nonstim in the same view.
# A MOFA factor loading on both nonstim (cell frequency) AND stim (cytokine
# response) for the same cell type is particularly strong biological evidence.
#
# Sample coverage (much smaller than nonstim):
#   Blood_stim: ~31 patients  |  Omentum_stim: ~28  |  SubQ_stim: ~18
# MOFA handles this via its missing-data mechanism — views with fewer
# samples get higher uncertainty and will be down-weighted by the ARD prior
# if signal is weak.
#
# Compositionality note: quadrant features (e.g. TNF+IFNγ+, TNF+IFNγ-,
# TNF-IFNγ+, TNF-IFNγ-) sum to 100% within each cell type × stimulation
# combination. The ZI filter will drop near-zero quadrants, and INT handles
# the remaining collinearity by working on ranks rather than raw values.

message("Loading stimulated flow cytometry data...")

stim_raw <- read.csv("data/flow_cytometry_perct_stim.csv",
                     check.names = FALSE, fileEncoding = "UTF-8-BOM")
names(stim_raw) <- make.unique(names(stim_raw))
stim_raw$Patient <- trimws(stim_raw$Patient)

subq_stim_cols    <- grep("^SubQ",    names(stim_raw), value = TRUE)
blood_stim_cols   <- grep("^Blood",   names(stim_raw), value = TRUE)
omentum_stim_cols <- grep("^Omentum", names(stim_raw), value = TRUE)

message(sprintf("  Stim features — SubQ: %d | Blood: %d | Omentum: %d",
                length(subq_stim_cols), length(blood_stim_cols), length(omentum_stim_cols)))

# Build raw matrices for QC and ZI assessment
subq_stim_raw    <- raw_view(stim_raw, "Patient", subq_stim_cols)
blood_stim_raw   <- raw_view(stim_raw, "Patient", blood_stim_cols)
omentum_stim_raw <- raw_view(stim_raw, "Patient", omentum_stim_cols)

message(sprintf("  Patients with stim data — SubQ: %d | Blood: %d | Omentum: %d",
                nrow(subq_stim_raw), nrow(blood_stim_raw), nrow(omentum_stim_raw)))

# Zero-inflation in stim data
# Stim data is especially prone to ZI: rare cytokine-producing subsets
# (e.g. IL17F+ MAIT in SubQ) may be undetectable in most patients.
# Apply the same ZERO_FRAC_THRESHOLD used for nonstim.

zero_inf_stim <- lapply(
  list(SubQ_stim    = subq_stim_raw,
       Blood_stim   = blood_stim_raw,
       Omentum_stim = omentum_stim_raw),
  function(raw_mat) {
    apply(raw_mat, 2, function(col) {
      col <- col[!is.na(col)]
      if (length(col) == 0) return(1)
      mean(col == 0)
    })
  }
)

zi_stim_df <- bind_rows(lapply(names(zero_inf_stim), function(vname) {
  data.frame(view = vname, feature = names(zero_inf_stim[[vname]]),
             zero_frac = zero_inf_stim[[vname]], stringsAsFactors = FALSE)
}))
zi_stim_df$zero_inflated <- zi_stim_df$zero_frac > ZERO_FRAC_THRESHOLD

zi_stim_summary <- zi_stim_df %>%
  group_by(view) %>%
  summarise(n_features = n(), n_zi = sum(zero_inflated),
            pct_zi = round(100 * mean(zero_inflated), 1), .groups = "drop")
message("  Stim zero-inflation (threshold >", ZERO_FRAC_THRESHOLD * 100, "% zeros):")
print(as.data.frame(zi_stim_summary))

write.csv(zi_stim_df, file.path(QC_DIR, "07_stim_zero_inflation.csv"), row.names = FALSE)

# Apply INT (with ZI filter) to stim views
apply_int_stim <- function(raw_mat, view_name) {
  zi <- zi_stim_df %>% filter(view == view_name)
  apply_int_to_view(raw_mat, zi, drop_zi = USE_INT)
}

subq_stim_mat    <- apply_int_stim(subq_stim_raw,    "SubQ_stim")
blood_stim_mat   <- apply_int_stim(blood_stim_raw,   "Blood_stim")
omentum_stim_mat <- apply_int_stim(omentum_stim_raw, "Omentum_stim")

message(sprintf("  Stim features after ZI filter + INT:"))
message(sprintf("    SubQ_stim: %d | Blood_stim: %d | Omentum_stim: %d",
                nrow(subq_stim_mat), nrow(blood_stim_mat), nrow(omentum_stim_mat)))

# ── 4. DEFINE GROUPS (T2D STATUS) ───────────────────────────
# MOFA+ multi-group mode: learns shared factors but allows group-specific
# factor activity. Reveals which variation is T2D-specific vs general.

# All samples present across all views (union)


# Union of all samples across all views (nonstim + stim).
# MOFA handles per-view missingness natively — samples absent from a view
# are represented as NA columns in that view's matrix.
all_samples <- Reduce(union, list(
  colnames(subq_mat),    colnames(blood_mat),    colnames(omentum_mat),
  colnames(subq_stim_mat), colnames(blood_stim_mat), colnames(omentum_stim_mat)
))

data_list <- list(
  SubQ         = subq_mat,
  Blood        = blood_mat,
  Omentum      = omentum_mat,
  SubQ_stim    = subq_stim_mat,
  Blood_stim   = blood_stim_mat,
  Omentum_stim = omentum_stim_mat
  # Cytokines  = cyto_mat   # uncomment to add plasma cytokines
)

for(i in 1:length(data_list)) {
  # print(i)
  # print(data_list[[i]])
  new_mat_cols <- setdiff(all_samples, colnames(data_list[[i]]))
  new_mat <- matrix(NA, nrow=nrow(data_list[[i]]), ncol=length(new_mat_cols))
  colnames(new_mat) <- new_mat_cols
  rownames(new_mat) <- rownames(data_list[[i]])
  
  data_list[[i]] <- cbind(data_list[[i]], new_mat)
  data_list[[i]] <- data_list[[i]][,all_samples]
  print(dim(data_list[[i]]))
}

# Build sample → group mapping
t2d_map <- clinical[all_samples[all_samples %in% rownames(clinical)], "T2D"]
names(t2d_map) <- all_samples[all_samples %in% rownames(clinical)]

sample_groups <- sapply(all_samples, function(s) {
  if (!s %in% names(t2d_map) || is.na(t2d_map[s]) || t2d_map[s] == "") {
    "Unknown"
  } else if (as.character(t2d_map[s]) == "1") {
    "T2D"
  } else {
    "Non_T2D"
  }
})

# ── 5. CREATE MOFA OBJECT ───────────────────────────────────

# MOFA2 input: named list of matrices (features × samples)
# Each matrix can have different sample coverage — missingness is fine


# Create multi-group MOFA object
mofa <- create_mofa(data_list) #, groups = sample_groups[all_samples])

message("MOFA object overview:")
print(mofa)

# ── 6. SET OPTIONS ──────────────────────────────────────────

data_opts <- get_default_data_options(mofa)
data_opts$scale_views  <- TRUE   # z-score each view (important: cytokines vs flow have different variance)
# data_opts$scale_groups <- FALSE  # do not scale between T2D/Non-T2D groups

model_opts <- get_default_model_options(mofa)
model_opts$num_factors         <- 10
model_opts$spikeslab_weights   <- TRUE   # sparse weights (ARD prior per factor×view)

train_opts <- get_default_training_options(mofa)
train_opts$seed              <- 42
train_opts$convergence_mode  <- "slow"   # ELBO convergence, more robust
train_opts$maxiter           <- 10000
train_opts$drop_factor_threshold <- 0.01  # drop factors explaining <1% variance

mofa <- prepare_mofa(
  mofa,
  data_options     = data_opts,
  model_options    = model_opts,
  training_options = train_opts
)

# ── 7. TRAIN ────────────────────────────────────────────────
# use_basilisk = TRUE auto-manages the Python mofapy2 environment via basilisk

mofa_file <- file.path(OUT_DIR, "mofa_model_medium_with_stim.hdf5")

mofa <- run_mofa(mofa, outfile = mofa_file, use_basilisk = TRUE)
message("MOFA+ training complete. Model saved to: ", mofa_file)

r2 <- get_variance_explained(mofa)$r2_per_factor[[1]]  # first group
# r2 is factors × views — sum across views per factor
total_r2 <- rowSums(r2)
plot(cumsum(total_r2), type = "b", xlab = "Factor", ylab = "Cumulative R²")



mofa <- load_model(mofa_file)

# ── 8. ATTACH METADATA ──────────────────────────────────────

cyto_raw_ <- cbind(Patient=cyto_raw[,"Patient"], cyto_raw[,rownames(cyto_mat)])

# Prepare a metadata table for all samples
meta <- merge(
  clinical[, c("Patient", "Sex", "Age", "BMI", "Glucose", "HbA1c",
                      "TotChol", "TG", "HDL", "LDL", "CRP", "Insulin",
                      "T2D", "HTN", "Dyslip", "OSA", "Asthma/COPD",
                      "LiverDx", "GORD", "Sleeve", "OAGB", "ASA", "OSMRS")],
  cyto_raw_,
  by="Patient",
  all.x=TRUE,
  all.y=TRUE
)

# Recode categorical variables
meta$Sex_F    <- as.integer(meta$Sex == "F")
meta$T2D      <- as.numeric(meta$T2D)
meta$HTN      <- as.numeric(meta$HTN)
meta$Dyslip   <- as.numeric(meta$Dyslip)
meta$OSA      <- as.numeric(meta$OSA)
meta$`Asthma/COPD` <- as.numeric(meta$`Asthma/COPD`)
meta$LiverDx  <- as.numeric(meta$LiverDx)
meta$GORD     <- as.numeric(meta$GORD)
meta$HbA1c    <- as.numeric(meta$HbA1c)  # stored as mmol/mol in this dataset
meta$BMI      <- as.numeric(meta$BMI)
meta$Age      <- as.numeric(meta$Age)
meta$Glucose  <- as.numeric(meta$Glucose)
meta$CRP      <- as.numeric(meta$CRP)
meta$Insulin  <- as.numeric(meta$Insulin)
meta$ASA  <- as.numeric(meta$ASA)
meta$OSMRS  <- as.numeric(meta$OSMRS)

for (col in colnames(cyto_raw_)[2:length(colnames(cyto_raw_))]) {
  meta[,col] <- as.numeric(meta[,col])
}

# Only keep samples in the MOFA model
mofa_samples <- unlist(samples_names(mofa))
meta_filtered <- meta[meta$Patient %in% mofa_samples, ]

colnames(meta_filtered)[1] <- "sample"
samples_metadata(mofa) <- meta_filtered

# ── 9. VARIANCE EXPLAINED ───────────────────────────────────

message("Plotting variance explained...")

p_r2 <- plot_variance_explained(mofa) +
  scale_fill_gradientn(
    colours = c("white", "#FEE090", "#FC8D59", "#D73027", "#4A0011"),
    limits  = c(0, NA),
    name    = expression(R^2~"(%)")
  ) +
  ggtitle("Variance explained per factor per view") +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(size = 13, face = "bold"),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y     = element_text(size = 15),
    legend.title    = element_text(size = 20),
    legend.text     = element_text(size = 20),
    strip.text       = element_blank(),
    strip.background = element_blank(),
    panel.grid       = element_blank()
  )

p_r2_total <- plot_variance_explained(mofa, plot_total = TRUE)[[2]] +
  scale_fill_manual(values = c(
    "Blood"        = "#E6550D",
    "Blood_stim"   = "#FDAE6B",
    "SubQ"         = "#3182BD",
    "SubQ_stim"    = "#9ECAE1",
    "Omentum"      = "#31A354",
    "Omentum_stim" = "#A1D99B",
    "Cytokines"    = "#756BB1"
  )) +
  ggtitle("Total variance explained per view") +
  theme_classic(base_size = 11) +
  theme(
    plot.title       = element_text(size = 13, face = "bold"),
    axis.text        = element_text(size = 20),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y     = element_text(size = 20),
    legend.title     = element_text(size = 20),
    legend.text      = element_text(size = 20),
    axis.title.y     = element_text(size = 20),
    strip.text       = element_blank(),
    strip.background = element_blank()
  )

ggsave(file.path(OUT_DIR, "01_variance_explained.pdf"),
       p_r2 / p_r2_total, width = 5, height = 10)
ggsave(file.path(OUT_DIR, "01_variance_explained.png"),
       p_r2 / p_r2_total, width = 5, height = 10, dpi = 250)


# ── 10. FACTOR CORRELATIONS WITH CLINICAL METADATA ──────────
# Key question: which factors associate with T2D, HTN, Dyslip, BMI, HbA1c?

message("Correlating factors with clinical covariates...")

clin_covariates <- c(
  "T2D", "HTN", "Dyslip",
  "ASA", "OSMRS",
  "BMI", "Age", "Sex_F",
  "HbA1c", "Glucose", "CRP", "Insulin",
  "TotChol", "TG", "HDL", "LDL",
  colnames(cyto_raw_)[2:length(colnames(cyto_raw_))]
)

# Only keep covariates with sufficient non-NA values
clin_covariates <- clin_covariates[sapply(clin_covariates, function(v) {
  col <- meta_filtered[[v]]
  sum(!is.na(as.numeric(col))) >= 10
})]

pdf(file.path(OUT_DIR, "02_factor_covariate_correlations.pdf"), width = 12, height = 8)
correlate_factors_with_covariates(
  mofa,
  covariates = clin_covariates,
  plot       = "log_pval",        # -log10(p) heatmap
  abs        = FALSE,             # show direction of correlation
  alpha      = 0.2               # mark significant cells
)
dev.off()

# Also save as PNG
png(file.path(OUT_DIR, "02_factor_covariate_correlations_r.png"),
    width = 1800, height = 900, res = 150)
correlate_factors_with_covariates(
  mofa,
  covariates = clin_covariates,
  plot       = "r",
  abs        = FALSE,
  alpha      = 0.1
) + theme_classic(base_size = 30) +
  scale_fill_viridis_c(option = "plasma") +
  scale_color_viridis_c(option = "plasma")
dev.off()

# Signed log-p heatmap: (-log10(p)) * sign(r), blue = negative, red = positive.
# Non-significant cells (p >= sig_alpha) shown at 25% opacity so the gradient
# direction is still visible but significance pops. Covariates in a fixed order.
{
  fac_df    <- get_factors(mofa, as.data.frame = TRUE)
  meta_mofa <- samples_metadata(mofa)
  sig_alpha <- 0.1

  all_factors   <- unique(fac_df$factor)
  factor_levels <- all_factors[order(as.integer(sub("Factor", "", all_factors)))]

  corr_rows <- do.call(rbind, lapply(factor_levels, function(fac) {
    fvals        <- fac_df$value[fac_df$factor == fac]
    names(fvals) <- fac_df$sample[fac_df$factor == fac]

    do.call(rbind, lapply(clin_covariates, function(cov) {
      cvals <- as.numeric(meta_mofa[[cov]])
      samps <- intersect(names(fvals), meta_mofa$sample)
      fv    <- fvals[samps]
      cv    <- cvals[match(samps, meta_mofa$sample)]
      ok    <- !is.na(fv) & !is.na(cv) & length(unique(cv[!is.na(cv)])) > 1
      if (sum(ok) < 5)
        return(data.frame(factor = fac, covariate = cov, r = NA_real_, p = NA_real_))
      ct <- tryCatch(cor.test(fv[ok], cv[ok], method = "pearson"), error = function(e) NULL)
      if (is.null(ct))
        return(data.frame(factor = fac, covariate = cov, r = NA_real_, p = NA_real_))
      data.frame(factor = fac, covariate = cov, r = as.numeric(ct$estimate), p = ct$p.value)
    }))
  }))

  corr_rows$signed_logp <- with(corr_rows, ifelse(is.na(p), 0, (-log10(p)) * sign(r)))
  corr_rows$significant <- !is.na(corr_rows$p) & corr_rows$p < sig_alpha

  corr_rows$covariate <- factor(corr_rows$covariate, levels = clin_covariates)
  corr_rows$factor    <- factor(corr_rows$factor, levels = rev(factor_levels))  # F1 at top

  clim <- max(abs(corr_rows$signed_logp), na.rm = TRUE)

  # Asterisk colour: white when fill is intense (|signed_logp| / clim > 0.45),
  # dark grey when near the white midpoint so it stays readable.
  corr_rows$star_color <- ifelse(
    abs(corr_rows$signed_logp) / clim > 0.45, "white", "#333333"
  )

  p_signed <- ggplot(corr_rows,
                     aes(x = covariate, y = factor, fill = signed_logp)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(data = subset(corr_rows, significant),
              aes(label = "*", color = star_color), size = 9, vjust = 0.75,
              fontface = "bold", lineheight = 0) +
    scale_color_identity() +
    scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "white",
      high     = "#B2182B",
      midpoint = 0,
      limits   = c(-clim, clim),
      name     = expression(-log[10](p) %.% sign(r))
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 90, hjust = 1, size = 28),
      axis.text.y      = element_text(size = 26),
      axis.title       = element_blank(),
      legend.position  = "right",
      legend.title     = element_text(size = 22),
      legend.text      = element_text(size = 20),
      legend.key.height = unit(2, "cm")
    ) +
    ggtitle(sprintf("Factor-covariate correlations  (* = p < %.2f)", sig_alpha))

  ggsave(file.path(OUT_DIR, "02_factor_covariate_correlations_signed.pdf"),
         p_signed, width = 22, height = 7)
  ggsave(file.path(OUT_DIR, "02_factor_covariate_correlations_signed.png"),
         p_signed, width = 22, height = 7, dpi = 250)
}

# ── 11. FACTOR SCATTER PLOTS ────────────────────────────────

message("Plotting factor scatter plots...")

# Plot top 3 factor pairs, coloured by T2D
factor_pairs <- list(c(1, 2), c(1, 3), c(2, 3), c(5, 6), c(5, 7), c(5, 8))

plots_t2d <- lapply(factor_pairs, function(fs) {
  plot_factors(mofa,
               factors  = fs,
               color_by = "T2D",
               dot_size = 3,
               stroke   = 0.3) +
    scale_color_manual(values = c("0" = "#3182BD", "1" = "#E6550D"),
                       labels = c("0" = "Non-T2D", "1" = "T2D"),
                       na.value = "grey70") +
    theme_classic(base_size = 10)
})

p_scatter_t2d <- wrap_plots(plots_t2d, ncol = 3)
ggsave(file.path(OUT_DIR, "03_factors_scatter_T2D.pdf"),
       p_scatter_t2d, width = 10, height = 5)
ggsave(file.path(OUT_DIR, "03_factors_scatter_T2D.png"),
       p_scatter_t2d, width = 10, height = 5, dpi = 150)

library(scales)

# Also colour by BMI and HbA1c
p_bmi <- plot_factors(mofa, factors = c(1, 2), color_by = "BMI", dot_size = 3) +
  scale_fill_viridis_c(option = "plasma") + theme_classic(base_size = 10)
p_hba1c <- plot_factors(mofa, factors = c(5, 7), color_by = "HbA1c", dot_size = 3) + 
  scale_fill_viridis_c(option = "plasma", limits = c(25, 70),  oob = squish) + theme_classic(base_size = 10)
p_hdl <- plot_factors(mofa, factors = c(5, 7), color_by = "HDL", dot_size = 3) +
  scale_fill_viridis_c(option = "plasma") + theme_classic(base_size = 10)
ggsave(file.path(OUT_DIR, "03_factors_scatter_BMI_HbA1c.pdf"),
       p_bmi + p_hba1c + p_hdl, width = 15, height = 5)

# ── 12. TOP FEATURE WEIGHTS PER FACTOR ──────────────────────

message("Plotting feature weights...")

n_factors_plot <- min(8, get_dimensions(mofa)$K)
views_to_plot  <- c("Blood", "Blood_stim", "SubQ", "SubQ_stim", "Omentum",
                    "Omentum_stim")
                    # add "Cytokines" here too if that view is enabled

pdf(file.path(OUT_DIR, "04_feature_weights.pdf"), width = 14, height = 20)
for (f in seq_len(n_factors_plot)) {
  plots <- lapply(views_to_plot, function(v) {
    tryCatch(
      plot_top_weights(mofa, view = v, factor = f, nfeatures = 20,
                       scale = TRUE, abs = FALSE) +
        ggtitle(paste0("Factor ", f, " — ", v)) +
        theme_classic(base_size = 9),
      error = function(e) NULL
    )
  })
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) > 0) {
    print(wrap_plots(plots, ncol = 2))
  }
}
dev.off()

# ── 12b. FOCUSED WEIGHT BARPLOTS: Factors 1 & 6 ─────────────────────────────
# Horizontal barplots, top N positive + negative per view, blue-red by weight.

{
  n_top <- 5

  wdf <- get_weights(mofa, as.data.frame = TRUE)
  wdf$factor_num <- as.integer(sub("Factor", "", wdf$factor))

  make_weight_barplot <- function(fnum) {
    sub_w <- wdf[wdf$factor_num == fnum, ]

    top_df <- do.call(rbind, lapply(unique(sub_w$view), function(v) {
      vd   <- sub_w[sub_w$view == v, ]
      vd   <- vd[order(vd$value), ]
      n    <- nrow(vd)
      keep <- unique(c(seq_len(min(n_top, n)),
                       seq(max(1L, n - n_top + 1L), n)))
      vd[keep, ]
    }))

    # Strip any tissue/view prefix from feature name for display
    top_df$label <- gsub(
      "^(Blood_stim|Blood|SubQ_stim|SubQ|Omentum_stim|Omentum)[_ ]+",
      "", top_df$feature
    )
    # Greek symbol substitutions (specific before general)
    top_df$label <- gsub("Monocgamma_tes", "Monocytes", top_df$label)
    top_df$label <- gsub("gamma_d T", "γδT", top_df$label)
    top_df$label <- gsub("gamma_d", "γδT", top_df$label)
    top_df$label <- gsub("gamma_",    "γ",         top_df$label)

    # Ordered y factor per view; use ||| as separator so view prefixes don't break regex
    top_df <- top_df[order(top_df$view, top_df$value), ]
    top_df$yid <- factor(
      paste(top_df$view, top_df$label, sep = "|||"),
      levels = paste(top_df$view, top_df$label, sep = "|||")
    )

    clim <- max(abs(top_df$value), na.rm = TRUE)

    ggplot(top_df, aes(x = value, y = yid, fill = value)) +
      geom_col(width = 0.75, color = NA) +
      geom_vline(xintercept = 0, linewidth = 0.5, color = "grey40") +
      facet_wrap(~ view, scales = "free_y", nrow = 1) +
      scale_fill_gradient2(
        low      = "#2166AC",
        mid      = "white",
        high     = "#B2182B",
        midpoint = 0,
        limits   = c(-clim, clim),
        name     = "Weight"
      ) +
      scale_x_continuous(limits = c(-clim, clim),
                         expand = expansion(mult = 0.05)) +
      scale_y_discrete(labels = function(x) sub("^.*\\|\\|\\|", "", x)) +
      theme_classic(base_size = 10) +
      theme(
        axis.title.y       = element_blank(),
        axis.text.y        = element_text(size = 26),
        axis.text.x        = element_text(size = 20),
        axis.title.x       = element_text(size = 24),
        strip.text         = element_text(size = 28, face = "bold"),
        strip.background   = element_blank(),
        legend.position    = "right",
        legend.title       = element_text(size = 28, margin = margin(b = 32)),
        legend.text        = element_text(size = 28),
        legend.key.height  = unit(3, "cm"),
        legend.key.width   = unit(0.8, "cm"),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
      ) +
      labs(x     = "Weight",
           title = paste0("Factor ", fnum, "  —  top ", n_top,
                          " feature weights per view"))
  }

  p_f1 <- make_weight_barplot(1)
  p_f6 <- make_weight_barplot(6)

  p_weights_focus <- p_f1 / p_f6

  ggsave(file.path(OUT_DIR, "05_weights_factor1_6.pdf"),
         p_weights_focus, width = 46, height = 16)
  ggsave(file.path(OUT_DIR, "05_weights_factor1_6.png"),
         p_weights_focus, width = 46, height = 16, dpi = 150)
}

library(ggpubr)
plot_data_scatter(mofa,
                  view = "Omentum",         # view of interest
                  factor = 5,             # factor of interest
                  features = 6,           # number of features to plot (they are selected by weight)
                  add_lm = TRUE,          # add linear regression
                  color_by = "T2D"
)
plot_data_scatter(mofa,
                  view = "Omentum",         # view of interest
                  factor = 7,             # factor of interest
                  features =6,           # number of features to plot (they are selected by weight)
                  add_lm = TRUE,          # add linear regression
                  color_by = "T2D"
)


# ── 13. FACTOR BEESWARM BY GROUP ────────────────────────────
# Visualise each factor's distribution across T2D vs Non-T2D

message("Plotting factor distributions by T2D group...")

p_groups <- MOFA2::plot_factors(mofa,
                          factors  = seq_len(n_factors_plot),
                          color_by = "T2D",
                          dot_size = 2.5) + #, stroke   = 0.2) + #, dodge    = TRUE
  scale_color_manual(values = c("0" = "#3182BD", "1" = "#E6550D"),
                     labels = c("0" = "Non-T2D", "1" = "T2D"),
                     na.value = "grey70") +
  theme_classic(base_size = 10)

ggsave(file.path(OUT_DIR, "05_factors_by_T2D_group.pdf"),
       p_groups, width = 14, height = 6)
ggsave(file.path(OUT_DIR, "05_factors_by_T2D_group.png"),
       p_groups, width = 14, height = 6, dpi = 150)

# ── 14. WEIGHT HEATMAP (multi-view overview) ─────────────────

message("Plotting weight heatmaps per view...")

pdf(file.path(OUT_DIR, "06_weight_heatmaps.pdf"), width = 14, height = 8)
for (v in views_to_plot) {
  tryCatch({
    p <- plot_weights_heatmap(
      mofa, view = v,
      factors = seq_len(n_factors_plot),
      show_colnames = FALSE,
      fontsize = 8
    )
    print(p)
  }, error = function(e) message("  Skipping heatmap for view: ", v, " — ", e$message))
}
dev.off()

p <- plot_weights_heatmap(
  mofa, view = v,
  factors = seq_len(n_factors_plot),
  show_colnames = FALSE,
  fontsize = 8
)
p

# ── 15. SAVE FACTOR VALUES & WEIGHTS FOR DOWNSTREAM USE ─────

message("Saving factor values and weights as CSV...")

# Factor (Z) values: samples × factors
Z_list <- get_factors(mofa, as.data.frame = TRUE)
Z_df   <- Z_list %>%
  pivot_wider(names_from = factor, values_from = value) %>%
  rename(Patient = sample)
write.csv(Z_df, file.path(OUT_DIR, "mofa_factor_values.csv"), row.names = FALSE)

# Weights: features × factors, per view
W_list <- get_weights(mofa, as.data.frame = TRUE)
write.csv(W_list, file.path(OUT_DIR, "mofa_weights_long.csv"), row.names = FALSE)

# Variance explained table
r2_df <- get_variance_explained(mofa)$r2_per_factor
if (is.list(r2_df)) {
  r2_df <- bind_rows(lapply(names(r2_df), function(g) {
    as.data.frame(r2_df[[g]]) %>%
      rownames_to_column("factor") %>%
      mutate(group = g)
  }))
}
write.csv(r2_df, file.path(OUT_DIR, "mofa_variance_explained.csv"), row.names = FALSE)

message("\n=== MOFA+ analysis complete ===")
message("Outputs in: ", OUT_DIR)
message("  mofa_model.hdf5                      — trained model")
message("  01_variance_explained.pdf/png        — R² per factor × view")
message("  02_factor_covariate_correlations.*   — correlation with clinical metadata")
message("  03_factors_scatter_*                 — factor scatter plots")
message("  04_feature_weights.pdf               — top loading features per factor")
message("  05_factors_by_T2D_group.*            — factor distributions T2D vs Non-T2D")
message("  06_weight_heatmaps.pdf               — weight overview heatmaps")
message("  mofa_factor_values.csv               — Z scores for downstream analysis")
message("  mofa_weights_long.csv                — all feature weights")

