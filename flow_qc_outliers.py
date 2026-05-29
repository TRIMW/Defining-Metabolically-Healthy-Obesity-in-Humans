#!/usr/bin/env python3
# flow_qc_outliers.py
# ============================================================
# Distribution histograms + outlier flags for all flow features
# Covers both: flow_cytometry_perct_nonstim.csv
#              flow_cytometry_counts_nonstim.csv
#
# Output: results/flow_qc/
#   perct_<view>_histograms.pdf   — one page per feature, percentages
#   counts_<view>_histograms.pdf  — one page per feature, counts
#   outlier_report.csv            — every (patient, feature, value) flagged
#   outlier_summary.csv           — per-feature: n outliers, threshold used
# ============================================================

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf_backend
from pathlib import Path
from scipy import stats

OUT_DIR = Path("results/flow_qc")
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size":   9,
    "axes.titlesize": 8.5,
    "figure.dpi":  120,
})

VIEWS   = ["SubQ", "Blood", "Omentum"]
N_COLS  = 4   # histograms per row on the page
N_ROWS  = 5   # rows per page  → 20 features per page

# Outlier detection: flag values beyond OUTLIER_IQR_K * IQR from Q1/Q3.
# IQR×3 (Tukey "far out") is more conservative than ×1.5 and appropriate
# for biological data where right-skew is normal.
OUTLIER_IQR_K = 3.0

VIEW_COLORS = {
    "perct":  {"SubQ": "#3182BD", "Blood": "#E6550D", "Omentum": "#31A354"},
    "counts": {"SubQ": "#6BAED6", "Blood": "#FD8D3C", "Omentum": "#74C476"},
}

# ── UNITS & PHYSIOLOGICAL REFERENCE RANGES ───────────────────
# Blood:   cells / µL of blood (lymphocyte gate, mononuclear fraction)
# Adipose: cells / g of tissue (stromal vascular fraction)
#
# Reference ranges are shown as green shaded bands on count histograms.
# Sources: standard haematology references + adipose immunology literature.
# These are approximate — use as sanity-check guides, not hard cutoffs.
#
# Structure: { col_substring: (lo, hi, label) }
# Matched by checking if the substring appears anywhere in the column name.

UNITS = {
    "Blood":   "cells / µL blood",
    "SubQ":    "cells / g tissue",
    "Omentum": "cells / g tissue",
}

UNITS_PERCT = {
    "Blood":   "% of parent gate",
    "SubQ":    "% of parent gate",
    "Omentum": "% of parent gate",
}

# Physiological reference bands for count data (lo, hi, description)
# Only shown when a column name contains the key substring.
PHYS_REFS = {
    # Blood — lymphocyte-gate CD45+ (not whole-blood WBC)
    "Blood_CD45+ of Total":    (500,  4000, "CD45⁺ lymph. gate /µL"),
    "Blood_CD3+ of CD45+":     (300,  2500, "CD3⁺ /µL"),
    "Blood_CD4+ % of CD45+":   (150,  1200, "CD4⁺ /µL"),   # these are count cols labelled oddly
    "Blood_CD8+ % of CD45+":   (80,   700,  "CD8⁺ /µL"),
    "Blood_CD3+ Count":        (700,  2100, "CD3⁺ absolute /µL"),
    "Blood_CD4+ Count":        (400,  1600, "CD4⁺ absolute /µL"),
    "Blood_CD8+ Count":        (200,  900,  "CD8⁺ absolute /µL"),
    # Adipose — cells/g SVF (literature: Wentworth 2010, Vohralik 2019)
    "SubQ_CD45+ of Total":     (50000,  1500000, "CD45⁺ /g SVF"),
    "Omentum_CD45+ of Total":  (80000,  2000000, "CD45⁺ /g SVF"),
    "SubQ_CD3+ of CD45+":      (5000,   500000,  "CD3⁺ /g SVF"),
    "Omentum_CD3+ of CD45+":   (8000,   600000,  "CD3⁺ /g SVF"),
    "SubQ_CD3+ Count":         (30000,  600000,  "CD3⁺ /g (count tube)"),
    "Omentum_CD3+ Count":      (50000,  800000,  "CD3⁺ /g (count tube)"),
    "SubQ_CD4+ Count":         (20000,  350000,  "CD4⁺ /g"),
    "Omentum_CD4+ Count":      (30000,  500000,  "CD4⁺ /g"),
    "SubQ_CD8+ Count":         (5000,   150000,  "CD8⁺ /g"),
    "Omentum_CD8+ Count":      (8000,   250000,  "CD8⁺ /g"),
}

# ── 1. LOAD ──────────────────────────────────────────────────

def load_file(path):
    df = pd.read_csv(path, encoding="utf-8-sig")
    df["Patient"] = df["Patient"].str.strip()
    df = df.set_index("Patient")
    # coerce everything to numeric
    df = df.apply(pd.to_numeric, errors="coerce")
    return df

perct_df  = load_file("data/flow_cytometry_perct_nonstim.csv")
counts_df = load_file("data/flow_cytometry_counts_nonstim.csv")

print(f"Percentages: {perct_df.shape[0]} patients, {perct_df.shape[1]} features")
print(f"Counts:      {counts_df.shape[0]} patients, {counts_df.shape[1]} features")

# ── 2. OUTLIER DETECTION ─────────────────────────────────────

def iqr_bounds(series, k=OUTLIER_IQR_K):
    """Return (lo, hi) Tukey fences using k×IQR. lo is clamped to 0."""
    q1, q3 = series.quantile(0.25), series.quantile(0.75)
    iqr = q3 - q1
    return max(0, q1 - k * iqr), q3 + k * iqr

def flag_outliers(df, data_label):
    """Return long-form DataFrame of flagged (patient, feature, value, bounds)."""
    records = []
    for col in df.columns:
        s = df[col].dropna()
        if len(s) < 4:
            continue
        lo, hi = iqr_bounds(s)
        flagged = s[(s < lo) | (s > hi)]
        for patient, val in flagged.items():
            records.append({
                "data":     data_label,
                "feature":  col,
                "patient":  patient,
                "value":    val,
                "fence_lo": lo,
                "fence_hi": hi,
                "iqr_k":    OUTLIER_IQR_K,
                "pct_rank": round(stats.percentileofscore(s, val), 1),
            })
    return pd.DataFrame(records)

outliers_perct  = flag_outliers(perct_df,  "percentages")
outliers_counts = flag_outliers(counts_df, "counts")
outliers_all    = pd.concat([outliers_perct, outliers_counts], ignore_index=True)

print("\nNotes on units:")
print("  Blood counts  → cells / µL of blood (lymphocyte scatter gate)")
print("  Adipose counts→ cells / g of tissue (SVF, bead-normalised)")
print("  *_Count columns (CD3+, CD4+, CD8+, MAIT, iNKT) → from separate")
print("    absolute counting tube (TruCount); lower n than bead-normalised cols")
print(f"\nOutlier flags (IQR×{OUTLIER_IQR_K}):")
print(f"  Percentages: {len(outliers_perct)} values flagged across "
      f"{outliers_perct['feature'].nunique()} features")
print(f"  Counts:      {len(outliers_counts)} values flagged across "
      f"{outliers_counts['feature'].nunique()} features")

outliers_all.to_csv(OUT_DIR / "outlier_report.csv", index=False)

# Per-feature summary
def outlier_summary(df, data_label):
    rows = []
    for col in df.columns:
        s = df[col].dropna()
        if len(s) < 4:
            continue
        lo, hi = iqr_bounds(s)
        n_out = int(((s < lo) | (s > hi)).sum())
        rows.append({
            "data":       data_label,
            "feature":    col,
            "n_obs":      len(s),
            "n_outliers": n_out,
            "pct_outliers": round(100 * n_out / len(s), 1),
            "fence_lo":   round(lo, 3),
            "fence_hi":   round(hi, 3),
            "median":     round(s.median(), 3),
            "max":        round(s.max(), 3),
            "min":        round(s.min(), 3),
        })
    return pd.DataFrame(rows)

summary = pd.concat([
    outlier_summary(perct_df,  "percentages"),
    outlier_summary(counts_df, "counts"),
], ignore_index=True)
summary = summary.sort_values(["data", "n_outliers"], ascending=[True, False])
summary.to_csv(OUT_DIR / "outlier_summary.csv", index=False)

# ── 3. HISTOGRAM PAGES ───────────────────────────────────────
# One page = N_ROWS × N_COLS features.
# Each histogram shows:
#   - distribution of all non-NA values (grey bars)
#   - IQR×k fences as vertical red dashed lines
#   - each outlier patient labelled by ID on the x-axis
#   - n_obs and n_outliers in subtitle

def short_name(col, prefix):
    """Strip the view prefix for a shorter axis title."""
    name = col
    for p in [prefix + "_", prefix]:
        if name.startswith(p):
            name = name[len(p):]
    # further shorten common suffixes
    name = name.replace("% of CD45+", "% CD45")
    name = name.replace("% of CD4+", "% CD4")
    name = name.replace("% of CD8+", "% CD8")
    name = name.replace("% of B cells", "% B")
    name = name.replace("% of MAIT", "% MAIT")
    name = name.replace("% of iNKT", "% iNKT")
    name = name.replace("% of Vd1", "% Vd1")
    name = name.replace("% of ATM", "% ATM")
    return name[:45]   # hard cap

def get_phys_ref(col):
    """Return (lo, hi, label) physiological reference for this column, or None."""
    for key, ref in PHYS_REFS.items():
        if col == key:
            return ref
    return None

def get_unit_label(col, file_key):
    """Return the appropriate unit string for x-axis labelling."""
    view = next((v for v in VIEWS if col.startswith(v)), None)
    if view is None:
        return ""
    if file_key == "perct":
        return UNITS_PERCT.get(view, "%")
    else:
        return UNITS.get(view, "cells / unit")

def plot_feature_hist(ax, series, col, prefix, color, feat_outliers, file_key="perct"):
    vals = series.dropna().values
    if len(vals) == 0:
        ax.set_visible(False)
        return

    lo, hi = iqr_bounds(series)

    # Histogram
    n_bins = min(30, max(8, len(vals) // 2))
    ax.hist(vals, bins=n_bins, color=color, alpha=0.72,
            edgecolor="white", linewidth=0.4)

    # Physiological reference band (counts only, where defined)
    phys = get_phys_ref(col)
    if phys and file_key == "counts":
        plo, phi, plabel = phys
        ax.axvspan(plo, phi, color="#31A354", alpha=0.10, zorder=0)
        ax.axvline(plo, color="#31A354", lw=0.8, ls=":", alpha=0.6)
        ax.axvline(phi, color="#31A354", lw=0.8, ls=":", alpha=0.6)
        # Label the band unobtrusively at top
        ax.text(0.98, 0.97, f"ref: {plabel}",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=5, color="#31A354", style="italic")

    # IQR×k fence lines
    ax.axvline(hi, color="#D73027", lw=1.3, ls="--", alpha=0.85)
    if lo > 0:
        ax.axvline(lo, color="#D73027", lw=1.3, ls="--", alpha=0.85)

    # Annotate outlier patients just above the x-axis
    if len(feat_outliers) > 0:
        ymax = ax.get_ylim()[1]
        for _, row in feat_outliers.iterrows():
            ax.annotate(
                row["patient"],
                xy=(row["value"], 0),
                xytext=(row["value"], ymax * 0.05),
                fontsize=5.5,
                color="#D73027",
                rotation=90,
                ha="center", va="bottom",
                arrowprops=dict(arrowstyle="-", color="#D73027", lw=0.6),
            )

    n_out = len(feat_outliers)
    title = short_name(col, prefix)
    subtitle = f"n={len(vals)}"
    if n_out:
        subtitle += f"  ⚠ {n_out} outlier{'s' if n_out > 1 else ''}"
    ax.set_title(f"{title}\n{subtitle}",
                 fontsize=7.5,
                 color="#D73027" if n_out else "black",
                 fontweight="bold" if n_out else "normal",
                 pad=2)
    ax.tick_params(labelsize=6.5)
    ax.set_xlabel(get_unit_label(col, file_key), fontsize=6, labelpad=1)
    ax.set_ylabel("n", fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)

def make_pdf(df, data_label, file_key, outlier_df):
    """Generate one PDF per view for a given data file."""
    for view in VIEWS:
        feat_cols = [c for c in df.columns if c.startswith(view + "_") or c.startswith(view)]
        # also catch columns without underscore (e.g. 'Blood_...')
        feat_cols = [c for c in df.columns if c.startswith(view)]
        if not feat_cols:
            continue

        color = VIEW_COLORS[file_key][view]
        n_feats = len(feat_cols)
        n_pages = int(np.ceil(n_feats / (N_ROWS * N_COLS)))

        out_path = OUT_DIR / f"{file_key}_{view.lower()}_histograms.pdf"
        with pdf_backend.PdfPages(out_path) as pdf:
            for page in range(n_pages):
                start = page * N_ROWS * N_COLS
                end   = min(start + N_ROWS * N_COLS, n_feats)
                page_feats = feat_cols[start:end]
                n_on_page  = len(page_feats)

                fig, axes = plt.subplots(
                    N_ROWS, N_COLS,
                    figsize=(N_COLS * 3.8, N_ROWS * 2.8),
                    constrained_layout=True
                )
                axes_flat = axes.flatten()

                unit_str = UNITS.get(view, "") if file_key == "counts" else UNITS_PERCT.get(view, "%")

                for i, col in enumerate(page_feats):
                    ax = axes_flat[i]
                    feat_out = outlier_df[outlier_df["feature"] == col]
                    plot_feature_hist(ax, df[col], col, view, color, feat_out, file_key)

                # Hide unused axes on last page
                for j in range(n_on_page, len(axes_flat)):
                    axes_flat[j].set_visible(False)

                fig.suptitle(
                    f"{data_label} | {view}  [{unit_str}]"
                    f"  —  features {start+1}–{end} of {n_feats}"
                    f"  (page {page+1}/{n_pages})"
                    f"  |  red dashed = IQR×{OUTLIER_IQR_K} fence"
                    + ("  |  green band = physiological reference" if file_key == "counts" else ""),
                    fontsize=9.5, fontweight="bold", y=1.01
                )
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

        print(f"  Saved: {out_path.name}  ({n_feats} features, {n_pages} pages)")

print("\nGenerating percentage histograms...")
make_pdf(perct_df,  "Flow percentages", "perct",  outliers_perct)

print("Generating count histograms...")
make_pdf(counts_df, "Flow counts",      "counts", outliers_counts)

# ── 4. OUTLIER HEATMAP (overview across all features) ────────
# One heatmap per file: patients × features, coloured by whether the value
# is flagged. Gives a bird's-eye view of which patients are systematic
# outliers across many features.

def outlier_heatmap(df, data_label, file_key, outlier_df, view):
    feat_cols = [c for c in df.columns if c.startswith(view)]
    if not feat_cols:
        return

    # Binary: 1 = flagged outlier, 0 = not, NaN = missing
    mat = pd.DataFrame(0.0, index=df.index, columns=feat_cols)
    for _, row in outlier_df[outlier_df["feature"].isin(feat_cols)].iterrows():
        if row["patient"] in mat.index and row["feature"] in mat.columns:
            mat.loc[row["patient"], row["feature"]] = 1.0
    # Mark missing as -1 for distinct colour
    for col in feat_cols:
        mat.loc[df[col].isna(), col] = np.nan

    # Shorten feature names
    short_cols = [short_name(c, view) for c in feat_cols]

    n_patients, n_feats = mat.shape
    fig_w = max(12, n_feats * 0.22)
    fig_h = max(5,  n_patients * 0.28)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    cmap = plt.cm.colors if hasattr(plt.cm, "colors") else None

    im = ax.imshow(
        mat.values,
        aspect="auto",
        cmap=plt.cm.RdYlGn_r,
        vmin=-0.05, vmax=1.05,
        interpolation="none"
    )

    # Patient labels
    ax.set_yticks(range(n_patients))
    ax.set_yticklabels(mat.index, fontsize=6.5)

    # Feature labels (rotated)
    ax.set_xticks(range(n_feats))
    ax.set_xticklabels(short_cols, rotation=90, ha="center", fontsize=5.5)

    # Count outliers per patient (row sum)
    n_out_per_patient = (mat == 1).sum(axis=1)
    for i, (pat, n) in enumerate(n_out_per_patient.items()):
        if n > 0:
            ax.text(n_feats + 0.3, i, f"  {n}✕", va="center",
                    fontsize=6, color="#D73027")

    ax.set_title(
        f"{data_label} | {view}: outlier heatmap\n"
        f"Red = flagged (IQR×{OUTLIER_IQR_K}), White = normal, Grey = missing\n"
        f"Numbers on right = total flagged features per patient",
        fontsize=9, fontweight="bold"
    )
    plt.colorbar(im, ax=ax, fraction=0.015, pad=0.01,
                 label="0 = normal  |  1 = outlier")
    plt.tight_layout()

    out_path = OUT_DIR / f"{file_key}_{view.lower()}_outlier_heatmap.png"
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path.name}")

print("\nGenerating outlier heatmaps...")
for view in VIEWS:
    outlier_heatmap(perct_df,  "Flow percentages", "perct",  outliers_perct,  view)
    outlier_heatmap(counts_df, "Flow counts",      "counts", outliers_counts, view)

# ── 5. PRINT TOP OUTLIER PATIENTS ────────────────────────────

print("\n=== TOP 15 MOST FLAGGED PATIENTS (percentages) ===")
top_perct = (outliers_perct.groupby("patient")
             .size().sort_values(ascending=False).head(15))
for pat, n in top_perct.items():
    print(f"  {pat}: {n} features flagged")

print("\n=== TOP 15 MOST FLAGGED PATIENTS (counts) ===")
top_counts = (outliers_counts.groupby("patient")
              .size().sort_values(ascending=False).head(15))
for pat, n in top_counts.items():
    print(f"  {pat}: {n} features flagged")

print("\n=== MOST EXTREME SINGLE VALUES (counts, top 10) ===")
extreme = (outliers_counts
           .assign(excess=lambda d: d["value"] / d["fence_hi"])
           .sort_values("excess", ascending=False)
           .head(10)
           [["patient", "feature", "value", "fence_hi", "excess"]])
extreme["value"]    = extreme["value"].round(1)
extreme["fence_hi"] = extreme["fence_hi"].round(1)
extreme["excess"]   = extreme["excess"].round(1)
print(extreme.to_string(index=False))

print(f"\nAll outputs in: {OUT_DIR}")
print("  perct_<view>_histograms.pdf     — percentage distributions")
print("  counts_<view>_histograms.pdf    — count distributions")
print("  perct/counts_<view>_outlier_heatmap.png — patient×feature overview")
print("  outlier_report.csv             — every flagged (patient, feature, value)")
print("  outlier_summary.csv            — per-feature outlier count + fences")
