import numpy as np
import pandas as pd
import matplotlib

from module_02_immune_landscape import load_flow_data

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.stats import spearmanr, mannwhitneyu, kruskal, chi2_contingency
from scipy.stats import shapiro
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.cluster import KMeans
from pathlib import Path
import warnings

warnings.filterwarnings("ignore")

import sys

sys.path.insert(0, ".")
from COHORTS import (
    COHORT_B3,
    CLINICAL_BLOODS,
    GLYCAEMIA,
    LIPIDS,
    INFLAMMATION,
    DIAGNOSES,
    STAT_GUIDANCE,
    cohort_list,
)

# ── Setup ─────────────────────────────────────────────────────
OUT = Path("results/PCA")
OUT.mkdir(parents=True, exist_ok=True)

PALETTE = {
    "T2D": "#D6604D",
    "No T2D": "#4393C3",
    "Female": "#E07BAD",
    "Male": "#5B8DB8",
    "Unknown": "#AAAAAA",
}

DIAG_COLORS = {
    "T2D": "#D6604D",
    "HTN": "#4393C3",
    "Dyslip": "#74C476",
    "OSA": "#9E9AC8",
    "LiverDx": "#F6A623",
}

DIAG_LABELS = {
    "T2D": "Type 2 Diabetes",
    "HTN": "Hypertension",
    "Dyslip": "Dyslipidaemia",
    "OSA": "Sleep Apnoea",
    "LiverDx": "Liver Disease",
}

VAR_LABELS = {
    "HbA1c": "HbA1c (%)",
    "Glucose": "Fasting glucose (mmol/L)",
    "Insulin": "Fasting insulin (mU/L)",
    "TotChol": "Total cholesterol (mmol/L)",
    "TG": "Triglycerides (mmol/L)",
    "HDL": "HDL cholesterol (mmol/L)",
    "LDL": "LDL cholesterol (mmol/L)",
    "nonHDL": "non-HDL cholesterol (mmol/L)",
    "CRP": "CRP (mg/L)",
    "ASA": "ASA score",
    "OSMRS": "OS-MRS score",
}


def _pc1_corr(pc1, df, min_n=3):
    """Spearman r + p between pc1 and every column of df.
    Returns DataFrame(r, p, annot) sorted by r, dropping columns with < min_n valid pairs.
    """
    rows = {}
    for col in df.columns:
        vals = df[col].values.astype(float)
        mask = ~np.isnan(vals)
        if mask.sum() < min_n:
            continue
        r, p = spearmanr(pc1[mask], vals[mask])
        stars = "***" if p < 0.0001 else "**" if p < 0.001 else "*" if p < 0.05 else ""
        rows[col] = {"r": r, "p": p, "annot": f"{r:.2f}{stars}"}
    return pd.DataFrame(rows).T.sort_values("r")


def plot_pca(sub, clinical, out_name, top_n=5, multiplex=None):
    """
    PCA.
    Colour by HbA1c (continuous severity) and by T2D status.
    Loadings plot shows which variables drive each component.
    This addresses: discrete subgroups vs continuous gradient.
    """

    # Mean imputation for missing values (report % missing per variable)
    missing_pct = sub.isnull().mean() * 100
    print("  PCA variable missingness:")
    for v, pct in missing_pct.items():
        if pct > 0:
            print(f"    {v}: {pct:.1f}% missing — mean-imputed")

    imp = SimpleImputer(strategy="mean")
    X = imp.fit_transform(sub)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    pca_vars = sub.columns

    pca = PCA(n_components=min(len(pca_vars), 4))
    scores = pca.fit_transform(X_scaled)
    loadings = pca.components_

    pc1_series = pd.Series(scores[:, 0], index=sub.index)
    pc2_series = pd.Series(scores[:, 1], index=sub.index)

    n_extra = 2 if multiplex is not None else 1
    fig = plt.figure(figsize=(32, 6 + 3 * n_extra))
    gs = fig.add_gridspec(
        1 + n_extra,
        3,
        height_ratios=[3] + [1.2] * n_extra,
        hspace=0.45,
        wspace=0.35,
    )
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    ax_heatmap = fig.add_subplot(gs[1, :])
    ax_heatmap_mx = fig.add_subplot(gs[2, :]) if multiplex is not None else None

    # Panel A: PC1 vs PC2, coloured by HbA1c
    ax = axes[0]
    hba1c = (
        clinical["HbA1c"].values
        if "HbA1c" in clinical.columns
        else np.zeros(len(clinical))
    )
    sc = ax.scatter(
        scores[:, 0],
        scores[:, 1],
        c=hba1c,
        cmap="RdYlBu_r",
        s=50,
        alpha=0.75,
        edgecolors="white",
        linewidths=0.4,
    )
    plt.colorbar(sc, ax=ax, label="HbA1c (%)")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)", fontsize=10)
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)", fontsize=10)
    ax.set_title(f"{out_name} PCA — coloured by HbA1c", fontsize=10)
    sns.despine(ax=ax)

    cluster_labels = None
    if out_name == "flow":
        km = KMeans(n_clusters=2, random_state=42, n_init=10)
        cluster_labels = km.fit_predict(scores[:, :2])
        cluster_labels_df = pd.DataFrame(
            cluster_labels, index=sub.index, columns=["Cluster"]
        )
        cluster_labels_df.to_csv(OUT / f"{out_name}_pca_clusters.csv")

    # Panel B: PC1 vs PC2, coloured by color_by value (coolwarm)
    ax = axes[1]
    color_by = "HTN"
    if color_by in clinical.columns:
        color_by_vals = clinical[color_by].values.astype(float)
        vmin, vmax = np.nanmin(color_by_vals), np.nanmax(color_by_vals)
        cmap_b = plt.cm.coolwarm
        norm_b = plt.Normalize(vmin=vmin, vmax=vmax)
        colors_pt = [
            cmap_b(norm_b(v)) if not np.isnan(v) else PALETTE["Unknown"]
            for v in color_by_vals
        ]
        ax.scatter(
            scores[:, 0],
            scores[:, 1],
            c=colors_pt,
            s=50,
            alpha=0.75,
            edgecolors="white",
            linewidths=0.4,
        )
        unique_vals = sorted(pd.Series(color_by_vals).dropna().unique())
        legend_handles = [
            mpatches.Patch(color=cmap_b(norm_b(v)), label=f"{color_by} = {v:.0f}")
            for v in unique_vals
        ]
        ax.legend(handles=legend_handles, fontsize=9)
        panel_b_title = f"{out_name} PCA — coloured by {color_by}"
    else:
        panel_b_title = f"{out_name} PCA"
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)", fontsize=10)
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)", fontsize=10)
    ax.set_title(panel_b_title, fontsize=10)
    sns.despine(ax=ax)

    # Panel C: Loadings biplot (PC1 and PC2) — top_n variables only
    ax = axes[2]
    loading_magnitude = np.abs(loadings[0]) + np.abs(loadings[1])
    top_idx = np.argsort(loading_magnitude)[-top_n:]
    for i, var in enumerate(pca_vars):
        if i not in top_idx:
            continue
        ax.arrow(
            0,
            0,
            loadings[0, i],
            loadings[1, i],
            head_width=0.02,
            head_length=0.02,
            fc="#D6604D",
            ec="#D6604D",
            alpha=0.8,
        )
        ax.text(
            loadings[0, i] * 1.12,
            loadings[1, i] * 1.12,
            VAR_LABELS.get(var, var),
            fontsize=8,
            ha="center",
        )
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.axvline(0, color="grey", lw=0.5, ls="--")
    circle = plt.Circle((0, 0), 1, fill=False, color="grey", lw=0.8, ls="--")
    ax.add_patch(circle)
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.set_xlabel("PC1 loadings", fontsize=10)
    ax.set_ylabel("PC2 loadings", fontsize=10)
    ax.set_title("PCA loadings  (variable contributions)", fontsize=10)
    ax.set_aspect("equal")
    sns.despine(ax=ax)

    # Panel D: Spearman correlations between PC1 and clinical variables (horizontal)
    clin_aligned = clinical.loc[sub.index].copy()
    cat_cols = clin_aligned.select_dtypes(
        include=["object", "category"]
    ).columns.tolist()
    num_cols = clin_aligned.select_dtypes(include=[np.number]).columns.tolist()
    clin_encoded = clin_aligned[num_cols].copy()
    for col in cat_cols:
        dummies = pd.get_dummies(clin_aligned[col], prefix=col)
        clin_encoded = pd.concat([clin_encoded, dummies], axis=1)
    clin_corr_pc1 = _pc1_corr(pc1_series.values, clin_encoded)
    clin_cols = clin_corr_pc1.index  # column order fixed by PC1 r sort
    clin_corr_pc2 = _pc1_corr(pc2_series.values, clin_encoded[clin_cols]).reindex(
        clin_cols
    )
    clin_labels = [VAR_LABELS.get(v, v) for v in clin_cols]
    clin_r_mat = np.vstack(
        [
            clin_corr_pc1["r"].values.astype(float),
            clin_corr_pc2["r"].values.astype(float),
        ]
    )
    clin_annot_mat = np.vstack(
        [
            clin_corr_pc1["annot"].values,
            clin_corr_pc2["annot"].values,
        ]
    )
    sns.heatmap(
        clin_r_mat,
        ax=ax_heatmap,
        cmap="RdBu_r",
        vmin=-1,
        vmax=1,
        annot=clin_annot_mat,
        fmt="",
        annot_kws={"size": 7},
        xticklabels=clin_labels,
        yticklabels=["PC1", "PC2"],
        cbar_kws={"label": "Spearman r", "shrink": 0.4},
        linewidths=0.3,
    )
    ax_heatmap.set_title(
        "PC1 & PC2 — clinical correlations (Spearman r; * p<0.05  ** p<0.001  *** p<0.0001)",
        fontsize=10,
    )
    ax_heatmap.tick_params(axis="x", labelsize=8, rotation=90)
    ax_heatmap.tick_params(axis="y", labelsize=9, rotation=0)

    # Panel E: Spearman correlations between PC1 and multiplex cytokines (horizontal)
    if multiplex is not None and ax_heatmap_mx is not None:
        common_idx = sub.index.intersection(multiplex.index)
        mx_aligned = multiplex.loc[common_idx].copy()
        imp_mx = SimpleImputer(strategy="mean")
        mx_enc = pd.DataFrame(
            imp_mx.fit_transform(mx_aligned),
            index=mx_aligned.index,
            columns=mx_aligned.columns,
        )
        pc1_mx = pc1_series.loc[common_idx].values
        pc2_mx = pc2_series.loc[common_idx].values
        mx_corr_pc1 = _pc1_corr(pc1_mx, mx_enc)
        mx_cols = mx_corr_pc1.index  # column order fixed by PC1 r sort
        mx_corr_pc2 = _pc1_corr(pc2_mx, mx_enc[mx_cols]).reindex(mx_cols)
        mx_r_mat = np.vstack(
            [
                mx_corr_pc1["r"].values.astype(float),
                mx_corr_pc2["r"].values.astype(float),
            ]
        )
        mx_annot_mat = np.vstack(
            [
                mx_corr_pc1["annot"].values,
                mx_corr_pc2["annot"].values,
            ]
        )
        sns.heatmap(
            mx_r_mat,
            ax=ax_heatmap_mx,
            cmap="RdBu_r",
            vmin=-1,
            vmax=1,
            annot=mx_annot_mat,
            fmt="",
            annot_kws={"size": 7},
            xticklabels=list(mx_cols),
            yticklabels=["PC1", "PC2"],
            cbar_kws={"label": "Spearman r", "shrink": 0.4},
            linewidths=0.3,
        )
        ax_heatmap_mx.set_title(
            f"PC1 & PC2 — multiplex cytokine correlations (Spearman r; * p<0.05  ** p<0.001  *** p<0.0001, n={len(common_idx)})",
            fontsize=10,
        )
        ax_heatmap_mx.tick_params(axis="x", labelsize=8, rotation=90)
        ax_heatmap_mx.tick_params(axis="y", labelsize=9, rotation=0)

    fig.suptitle(
        f"{out_name} PCA  (n={len(sub)}, " f"missing values mean-imputed)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(OUT / f"{out_name}_pca.pdf", bbox_inches="tight", dpi=150)
    plt.close(fig)

    # Variance explained summary
    var_exp = pd.Series(
        pca.explained_variance_ratio_,
        index=[f"PC{i+1}" for i in range(len(pca.explained_variance_ratio_))],
    )
    var_exp.to_csv(OUT / f"{out_name}_pca_variance_explained.csv")
    loadings_df = pd.DataFrame(
        loadings.T,
        index=pca_vars,
        columns=[f"PC{i+1}" for i in range(loadings.shape[0])],
    )
    loadings_df.to_csv(OUT / f"{out_name}_pca_loadings.csv")
    scores_df = pd.DataFrame(
        scores,
        index=sub.index,
        columns=[f"PC{i+1}" for i in range(scores.shape[1])],
    )
    scores_df.to_csv(OUT / f"{out_name}_pca_scores.csv")
    print(
        f"  Saved: {out_name}_pca.pdf  (PC1={pca.explained_variance_ratio_[0]*100:.1f}%, "
        f"PC2={pca.explained_variance_ratio_[1]*100:.1f}%)"
    )


def plot_pca_bulk(pca_df, pct_var, clinical, out_name, multiplex=None):
    """
    Plot PCA from pre-computed bulkRNAseq coordinates.
    pca_df    : DataFrame with columns PC1, PC2, participant (int)
    pct_var   : Series indexed 1,2 – variance explained (%)
    clinical  : full clinical DataFrame (index P01, P02 …)
    out_name  : e.g. "bulk_omentum"
    """
    pca_df = pca_df.copy()
    pca_df.index = [f"P{int(p):02d}" for p in pca_df["participant"]]

    common = pca_df.index.intersection(clinical.index)
    pca_df = pca_df.loc[common]
    clin = clinical.loc[common]

    pc1 = pca_df["PC1"].values
    pc2 = pca_df["PC2"].values
    pc1_series = pca_df["PC1"]
    pc2_series = pca_df["PC2"]
    pct1 = float(pct_var.iloc[0])
    pct2 = float(pct_var.iloc[1])

    n_extra = 2 if multiplex is not None else 1
    fig = plt.figure(figsize=(21, 6 + 3 * n_extra))
    gs = fig.add_gridspec(
        1 + n_extra,
        3,
        height_ratios=[3] + [1.2] * n_extra,
        hspace=0.45,
        wspace=0.35,
    )
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    ax_heatmap = fig.add_subplot(gs[1, :])
    ax_heatmap_mx = fig.add_subplot(gs[2, :]) if multiplex is not None else None

    # Panel A: coloured by HbA1c
    ax = axes[0]
    hba1c = clin["HbA1c"].values if "HbA1c" in clin.columns else np.zeros(len(clin))
    sc = ax.scatter(
        pc1,
        pc2,
        c=hba1c,
        cmap="RdYlBu_r",
        s=50,
        alpha=0.75,
        edgecolors="white",
        linewidths=0.4,
    )
    plt.colorbar(sc, ax=ax, label="HbA1c (%)")
    ax.set_xlabel(f"PC1 ({pct1:.1f}%)", fontsize=10)
    ax.set_ylabel(f"PC2 ({pct2:.1f}%)", fontsize=10)
    ax.set_title(f"{out_name} PCA — coloured by HbA1c", fontsize=10)
    sns.despine(ax=ax)

    # Panel B: coloured by OS-MRS (continuous coolwarm)
    ax = axes[1]
    colorby = "ASA"
    if colorby in clin.columns:
        osmrs_vals = clin[colorby].values.astype(float)
        vmin_o, vmax_o = np.nanmin(osmrs_vals), np.nanmax(osmrs_vals)
        cmap_b = plt.cm.coolwarm
        norm_b = plt.Normalize(vmin=vmin_o, vmax=vmax_o)
        colors_pt = [
            cmap_b(norm_b(v)) if not np.isnan(v) else PALETTE["Unknown"]
            for v in osmrs_vals
        ]
        ax.scatter(
            pc1, pc2, c=colors_pt, s=50, alpha=0.75, edgecolors="white", linewidths=0.4
        )
        unique_vals = sorted(pd.Series(osmrs_vals).dropna().unique())
        legend_handles = [
            mpatches.Patch(color=cmap_b(norm_b(v)), label=f"{colorby} = {v:.0f}")
            for v in unique_vals
        ]
        ax.legend(handles=legend_handles, fontsize=9)
        panel_b_title = f"{out_name} PCA — coloured by {colorby}"
    else:
        panel_b_title = f"{out_name} PCA"
    ax.set_xlabel(f"PC1 ({pct1:.1f}%)", fontsize=10)
    ax.set_ylabel(f"PC2 ({pct2:.1f}%)", fontsize=10)
    ax.set_title(panel_b_title, fontsize=10)
    sns.despine(ax=ax)

    # Panel C: coloured by T2D (discrete)
    ax = axes[2]
    if "T2D" in clin.columns:
        t2d_vals = clin["T2D"].values
        colors_c = [
            (
                PALETTE["T2D"]
                if t == 1
                else PALETTE["No T2D"] if t == 0 else PALETTE["Unknown"]
            )
            for t in t2d_vals
        ]
        ax.scatter(
            pc1, pc2, c=colors_c, s=50, alpha=0.75, edgecolors="white", linewidths=0.4
        )
        handles_c = [
            mpatches.Patch(color=PALETTE["T2D"], label="T2D"),
            mpatches.Patch(color=PALETTE["No T2D"], label="No T2D"),
        ]
        ax.legend(handles=handles_c, fontsize=9)
        panel_c_title = f"{out_name} PCA — coloured by T2D"
    else:
        panel_c_title = f"{out_name} PCA"
    ax.set_xlabel(f"PC1 ({pct1:.1f}%)", fontsize=10)
    ax.set_ylabel(f"PC2 ({pct2:.1f}%)", fontsize=10)
    ax.set_title(panel_c_title, fontsize=10)
    sns.despine(ax=ax)

    # Heatmap: PC1 & PC2 vs clinical variables
    clin_cat = clin.select_dtypes(include=["object", "category"]).columns.tolist()
    clin_num = clin.select_dtypes(include=[np.number]).columns.tolist()
    clin_encoded = clin[clin_num].copy()
    for col in clin_cat:
        dummies = pd.get_dummies(clin[col], prefix=col)
        clin_encoded = pd.concat([clin_encoded, dummies], axis=1)
    clin_corr_pc1 = _pc1_corr(pc1_series.values, clin_encoded)
    clin_cols = clin_corr_pc1.index
    clin_corr_pc2 = _pc1_corr(pc2_series.values, clin_encoded[clin_cols]).reindex(
        clin_cols
    )
    clin_labels = [VAR_LABELS.get(v, v) for v in clin_cols]
    clin_r_mat = np.vstack(
        [
            clin_corr_pc1["r"].values.astype(float),
            clin_corr_pc2["r"].values.astype(float),
        ]
    )
    clin_annot_mat = np.vstack(
        [
            clin_corr_pc1["annot"].values,
            clin_corr_pc2["annot"].values,
        ]
    )
    sns.heatmap(
        clin_r_mat,
        ax=ax_heatmap,
        cmap="RdBu_r",
        vmin=-1,
        vmax=1,
        annot=clin_annot_mat,
        fmt="",
        annot_kws={"size": 7},
        xticklabels=clin_labels,
        yticklabels=["PC1", "PC2"],
        cbar_kws={"label": "Spearman r", "shrink": 0.4},
        linewidths=0.3,
    )
    ax_heatmap.set_title(
        "PC1 & PC2 — clinical correlations (Spearman r; * p<0.05  ** p<0.001  *** p<0.0001)",
        fontsize=10,
    )
    ax_heatmap.tick_params(axis="x", labelsize=8, rotation=90)
    ax_heatmap.tick_params(axis="y", labelsize=9, rotation=0)

    # Heatmap: PC1 & PC2 vs multiplex cytokines
    if multiplex is not None and ax_heatmap_mx is not None:
        common_mx = pca_df.index.intersection(multiplex.index)
        mx_aligned = multiplex.loc[common_mx].copy()
        imp_mx = SimpleImputer(strategy="mean")
        mx_enc = pd.DataFrame(
            imp_mx.fit_transform(mx_aligned),
            index=mx_aligned.index,
            columns=mx_aligned.columns,
        )
        pc1_mx = pc1_series.loc[common_mx].values
        pc2_mx = pc2_series.loc[common_mx].values
        mx_corr_pc1 = _pc1_corr(pc1_mx, mx_enc)
        mx_cols = mx_corr_pc1.index
        mx_corr_pc2 = _pc1_corr(pc2_mx, mx_enc[mx_cols]).reindex(mx_cols)
        mx_r_mat = np.vstack(
            [
                mx_corr_pc1["r"].values.astype(float),
                mx_corr_pc2["r"].values.astype(float),
            ]
        )
        mx_annot_mat = np.vstack(
            [
                mx_corr_pc1["annot"].values,
                mx_corr_pc2["annot"].values,
            ]
        )
        sns.heatmap(
            mx_r_mat,
            ax=ax_heatmap_mx,
            cmap="RdBu_r",
            vmin=-1,
            vmax=1,
            annot=mx_annot_mat,
            fmt="",
            annot_kws={"size": 7},
            xticklabels=list(mx_cols),
            yticklabels=["PC1", "PC2"],
            cbar_kws={"label": "Spearman r", "shrink": 0.4},
            linewidths=0.3,
        )
        ax_heatmap_mx.set_title(
            f"PC1 & PC2 — multiplex cytokine correlations "
            f"(Spearman r; * p<0.05  ** p<0.001  *** p<0.0001, n={len(common_mx)})",
            fontsize=10,
        )
        ax_heatmap_mx.tick_params(axis="x", labelsize=8, rotation=90)
        ax_heatmap_mx.tick_params(axis="y", labelsize=9, rotation=0)

    fig.suptitle(f"{out_name} PCA  (n={len(pca_df)})", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / f"{out_name}_pca.pdf", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_name}_pca.pdf  (PC1={pct1:.1f}%, PC2={pct2:.1f}%)")


def run(
    multiplex_path="data/multiplex_cytokines.csv",
    flow_path="data/flow_cytometry_perct_nonstim.csv",
    stim_path="data/flow_cytometry_perct_stim.csv",
    clinical_path="data/clinical_data.csv",
):

    flow_cytometry = pd.read_csv(flow_path, index_col=0)
    multiplex_cytokines = pd.read_csv(multiplex_path, index_col=0)
    clinical = pd.read_csv(clinical_path, index_col=0)

    flow_cohort = COHORT_B3
    multiplex_cohort = multiplex_cytokines.index

    # plot_pca(
    #     flow_cytometry.loc[flow_cohort, :],
    #     clinical.loc[flow_cohort],
    #     "flow",
    #     multiplex=multiplex_cytokines,
    # )

    # plot_pca(
    #     multiplex_cytokines.loc[multiplex_cohort, :],
    #     clinical.loc[multiplex_cohort],
    #     "multiplex",
    # )

    # # ── bulkRNAseq PCA (pre-computed coordinates) ─────────────────
    # bulk_dir = Path("data/bulk/reanalysis/deseq2_samples_ok")
    # bulk_tissues = {
    #     "bulk_omentum": ("O.pca_data.csv", "O.pca_percentVar.csv"),
    #     "bulk_subq": ("S.pca_data.csv", "S.pca_percentVar.csv"),
    #     # "bulk_liver": ("L.pca_data.csv", "L.pca_percentVar.csv"),
    # }
    # for out_tag, (data_file, var_file) in bulk_tissues.items():
    #     bulk_df = pd.read_csv(bulk_dir / data_file, index_col=0)
    #     pct_var = pd.read_csv(bulk_dir / var_file, index_col=0)["x"]
    #     plot_pca_bulk(
    #         bulk_df,
    #         pct_var,
    #         clinical,
    #         out_tag,
    #         multiplex=None,
    #     )

    flow_tissues = load_flow_data(flow_path)
    stim_tissues = load_flow_data(stim_path)

    # ── Nonstim PCA (per tissue) ──────────────────────────────
    plot_pca(
        flow_tissues["Blood"].loc[flow_cohort, :],
        clinical.loc[flow_cohort],
        "flow_blood",
        multiplex=multiplex_cytokines,
    )
    plot_pca(
        flow_tissues["Omentum"].loc[flow_cohort, :],
        clinical.loc[flow_cohort],
        "flow_omentum",
        multiplex=multiplex_cytokines,
    )
    plot_pca(
        flow_tissues["SubQ"].loc[flow_cohort, :],
        clinical.loc[flow_cohort],
        "flow_subq",
        multiplex=multiplex_cytokines,
    )

    # ── Stimulated PCA (per tissue) ───────────────────────────
    # Stim data has fewer patients than nonstim; intersect with COHORT_B3
    # and with the clinical index to avoid KeyErrors.
    for tissue, tag in [
        ("Blood", "flow_blood_stim"),
        ("Omentum", "flow_omentum_stim"),
        ("SubQ", "flow_subq_stim"),
    ]:
        stim_df = stim_tissues[tissue]
        stim_cohort = stim_df.index.intersection(flow_cohort).intersection(
            clinical.index
        )
        if len(stim_cohort) < 5:
            print(f"  Skipping {tag}: only {len(stim_cohort)} patients with stim data")
            continue
        plot_pca(
            stim_df.loc[stim_cohort, :],
            clinical.loc[stim_cohort],
            tag,
            multiplex=multiplex_cytokines,
        )


if __name__ == "__main__":
    run()
