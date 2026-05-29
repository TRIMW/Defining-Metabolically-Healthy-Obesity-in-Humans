"""
Step 1: Preprocessing, Harmony integration, UMAP comparison, and pseudotime seed.

Outputs (in outputs/pseudotime/):
  - umap_integration_comparison.pdf   – 4x2 grid: non-integrated vs Harmony
  - umap_pseudotime.pdf               – DPT pseudotime on Harmony UMAP
  - pseudotime_distributions.pdf      – histograms by cell type / cohort / depot
  - fap_preadip_processed.h5ad        – combined object with Harmony PCA + DPT

Usage:
  pip install scanpy harmonypy matplotlib pandas
  python pseudotime_01_preprocessing.py
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

try:
    import harmonypy as hm
except ImportError:
    raise ImportError("Install harmonypy:  pip install harmonypy")

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_IN  = "obesity_scRNAseq/inputs/anndatas"
BASE_OUT = "obesity_scRNAseq/outputs/pseudotime"

# ── 1. Load, label, concatenate ───────────────────────────────────────────────
print("Loading data...")
om = sc.read_h5ad(f"{BASE_IN}/omentum_reannotated_finely_truncated.h5ad")
sq = sc.read_h5ad(f"{BASE_IN}/sq_reannotated_finely_truncated.h5ad")

om.obs["depot"] = "Omentum"
sq.obs["depot"] = "Subcutaneous"

# Store per-dataset UMAPs before concatenation (to show original embeddings)
om.obsm["X_umap_orig"] = om.obsm["X_umap"].copy()
sq.obsm["X_umap_orig"] = sq.obsm["X_umap"].copy()

adata_full = sc.concat(
    [om, sq],
    keys=["Om", "SQ"],
    label="dataset",
    merge="same",
    uns_merge="same",
)
print(f"Combined dataset: {adata_full.shape}")

# ── 2. Subset to FAPs + Pre_adipocytes ────────────────────────────────────────
adata = adata_full[
    adata_full.obs["cell_type_"].isin(["FAPs", "Pre_adipocytes"])
].copy()
print(f"After subset: {adata.shape}")
print(adata.obs.groupby(["depot", "cell_type_"]).size().to_string())

# ── 3. Reprocess from raw counts ──────────────────────────────────────────────
# Keep raw counts for tradeSeq export
adata.layers["counts"] = adata.layers["counts"]  # already present

adata.X = adata.layers["counts"].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["logcounts"] = adata.X.copy()

# HVG selection, batch-aware
sc.pp.highly_variable_genes(
    adata, n_top_genes=3000, batch_key="common_participant", flavor="seurat_v3",
    layer="counts",
)
print(f"HVGs selected: {adata.var['highly_variable'].sum()}")

# Scale and PCA on HVGs
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=50, svd_solver="arpack")
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"].copy()
adata.uns["pca"] = adata_hvg.uns["pca"].copy()

# ── 4. Non-integrated neighbors + UMAP ───────────────────────────────────────
print("Computing non-integrated UMAP...")
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30, key_added="neighbors_noint")
sc.tl.umap(adata, neighbors_key="neighbors_noint")
adata.obsm["X_umap_noint"] = adata.obsm["X_umap"].copy()

# ── 5. Harmony integration on participant ────────────────────────────────────
print("Running Harmony integration on common_participant...")
pca_mat = adata.obsm["X_pca"].copy()
meta    = adata.obs[["common_participant"]].copy()

ho = hm.run_harmony(
    pca_mat, meta,
    vars_use=["common_participant"],
    max_iter_harmony=30,
    random_state=42,
)
adata.obsm["X_pca_harmony"] = ho.Z_corr.T
print("Harmony done.")

sc.pp.neighbors(
    adata, n_neighbors=30, n_pcs=30,
    use_rep="X_pca_harmony", key_added="neighbors_harmony",
)
sc.tl.umap(adata, neighbors_key="neighbors_harmony")
adata.obsm["X_umap_harmony"] = adata.obsm["X_umap"].copy()

# ── 6. UMAP comparison: non-integrated vs Harmony ────────────────────────────
print("Plotting UMAP comparison...")
PALETTE = {
    "cell_type_": None,
    "cohort":     {"Healthy": "#2196F3", "Unhealthy": "#F44336"},
    "depot":      {"Omentum": "#FF9800", "Subcutaneous": "#9C27B0"},
    "sex":        {"F": "#E91E63", "M": "#00BCD4"},
}

fig, axes = plt.subplots(2, 4, figsize=(22, 11))
fig.suptitle("UMAP: Non-integrated  vs  Harmony-integrated (FAPs + Pre_adipocytes)",
             fontsize=13, y=1.01)

panels = [
    ("cell_type_", "umap_noint",   "Non-integrated — Cell type"),
    ("cohort",     "umap_noint",   "Non-integrated — Cohort"),
    ("depot",      "umap_noint",   "Non-integrated — Depot"),
    ("sex",        "umap_noint",   "Non-integrated — Sex"),
    ("cell_type_", "umap_harmony", "Harmony — Cell type"),
    ("cohort",     "umap_harmony", "Harmony — Cohort"),
    ("depot",      "umap_harmony", "Harmony — Depot"),
    ("sex",        "umap_harmony", "Harmony — Sex"),
]

for ax, (color_by, basis, title) in zip(axes.flat, panels):
    palette = PALETTE[color_by]
    sc.pl.embedding(
        adata, basis=basis, color=color_by, ax=ax, show=False,
        title=title, frameon=False,
        palette=palette,
        legend_loc="on data" if color_by == "cell_type_" else "right margin",
        size=8,
    )

plt.tight_layout()
plt.savefig(f"{BASE_OUT}/umap_integration_comparison.pdf", bbox_inches="tight", dpi=150)
plt.savefig(f"{BASE_OUT}/umap_integration_comparison.png", bbox_inches="tight", dpi=150)
print(f"  Saved: umap_integration_comparison.pdf/png")
plt.close()

# ── 7. Diffusion pseudotime (DPT) ─────────────────────────────────────────────
# Using Harmony-integrated neighbors for the diffusion map
print("Computing diffusion map + DPT...")
sc.tl.diffmap(adata, neighbors_key="neighbors_harmony", n_comps=15)

# Root: FAP cell closest to centroid of all FAPs in diffusion space
fap_mask  = adata.obs["cell_type_"] == "FAPs"
fap_dc    = adata.obsm["X_diffmap"][fap_mask]
centroid  = fap_dc.mean(axis=0)
dists     = np.linalg.norm(fap_dc - centroid, axis=1)
root_idx  = int(np.where(fap_mask)[0][np.argmin(dists)])

adata.uns["iroot"] = root_idx
sc.tl.dpt(adata, n_dcs=10)
print(f"  Root cell index: {root_idx}  "
      f"(cell type: {adata.obs['cell_type_'].iloc[root_idx]})")

# Per-group pseudotime stats
for grp in ["cell_type_", "cohort", "depot"]:
    stats = adata.obs.groupby(grp)["dpt_pseudotime"].describe()[["mean","std","50%"]]
    print(f"\nPseudotime by {grp}:\n{stats.to_string()}")

# ── 8. Pseudotime UMAP plots ──────────────────────────────────────────────────
fig, axes = plt.subplots(1, 4, figsize=(22, 5))
fig.suptitle("DPT Pseudotime — Harmony UMAP", fontsize=12)

for ax, (color, title) in zip(axes, [
    ("dpt_pseudotime", "Pseudotime"),
    ("cell_type_",     "Cell type"),
    ("cohort",         "Cohort"),
    ("depot",          "Depot"),
]):
    sc.pl.embedding(
        adata, basis="umap_harmony", color=color, ax=ax, show=False,
        title=title, frameon=False,
        color_map="viridis" if color == "dpt_pseudotime" else None,
        palette=PALETTE.get(color),
        legend_loc="on data" if color == "cell_type_" else "right margin",
        size=8,
    )

plt.tight_layout()
plt.savefig(f"{BASE_OUT}/umap_pseudotime.pdf", bbox_inches="tight", dpi=150)
plt.savefig(f"{BASE_OUT}/umap_pseudotime.png", bbox_inches="tight", dpi=150)
print(f"\n  Saved: umap_pseudotime.pdf/png")
plt.close()

# ── 9. Pseudotime distribution plots ─────────────────────────────────────────
COLORS = {
    "FAPs":            "#FF9800",
    "Pre_adipocytes":  "#4CAF50",
    "Healthy":         "#2196F3",
    "Unhealthy":       "#F44336",
    "Omentum":         "#FF9800",
    "Subcutaneous":    "#9C27B0",
}

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
fig.suptitle("DPT Pseudotime distributions", fontsize=12)

for ax, grp in zip(axes, ["cell_type_", "cohort", "depot"]):
    for cat in sorted(adata.obs[grp].unique()):
        vals = adata.obs.loc[adata.obs[grp] == cat, "dpt_pseudotime"]
        ax.hist(vals, bins=50, density=True, alpha=0.6,
                label=cat, color=COLORS.get(cat))
    ax.set_xlabel("Pseudotime")
    ax.set_ylabel("Density")
    ax.set_title(f"By {grp}")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

plt.tight_layout()
plt.savefig(f"{BASE_OUT}/pseudotime_distributions.pdf", bbox_inches="tight", dpi=150)
plt.savefig(f"{BASE_OUT}/pseudotime_distributions.png", bbox_inches="tight", dpi=150)
print(f"  Saved: pseudotime_distributions.pdf/png")
plt.close()

# ── 10. Save processed object for R ──────────────────────────────────────────
# Restore logcounts to .X so zellkonverter sees it as normalized expression
adata.X = adata.layers["logcounts"].copy()

# Drop scaled matrix (stored in HVG subset, not in adata)
out_path = f"{BASE_OUT}/fap_preadip_processed.h5ad"
adata.write_h5ad(out_path)
print(f"\nSaved processed object → {out_path}")
print("  obsm keys:", list(adata.obsm.keys()))
print("  layers   :", list(adata.layers.keys()))
print("\nDone. Next: run pseudotime_02_tradeseq.R")
