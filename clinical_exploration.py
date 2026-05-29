#!/usr/bin/env python3
# clinical_exploration.py
# ============================================================
# Clinical Data Exploration
#   1. Comorbidity UpSet plot (T2D, HTN, Dyslip, OSA, Asthma, GORD)
#   2. Continuous variable distributions by T2D status
#   3. Logistic regression: how well do basic clinical variables
#      (Age, Sex, BMI, blood measurements) predict T2D and other
#      diagnoses — serves as baseline to compare MOFA factors against
#   4. Association heatmap (all binary diagnoses × continuous vars)
# ============================================================
# Requirements:
#   pip install pandas numpy matplotlib seaborn statsmodels scikit-learn upsetplot scipy

import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from pathlib import Path
from scipy import stats
from scipy.stats import mannwhitneyu, chi2_contingency, spearmanr
from statsmodels.formula.api import logit
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import LeaveOneOut, cross_val_score
import upsetplot

OUT_DIR = Path("results/clinical_exploration")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── STYLE ────────────────────────────────────────────────────
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "figure.dpi": 150,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
    }
)
PALETTE_T2D = {"Non-T2D": "#3182BD", "T2D": "#E6550D"}

# ── 1. LOAD DATA ─────────────────────────────────────────────

clin = pd.read_csv("data/clinical_data.csv")
clin["Patient"] = clin["Patient"].str.strip()
clin = clin.set_index("Patient")

# Recode Sex
clin["Sex_F"] = (clin["Sex"] == "F").astype(float)
clin.loc[clin["Sex"].isna(), "Sex_F"] = np.nan

# Numeric coercion
numeric_cols = [
    "Age",
    "BMI",
    "Glucose",
    "HbA1c",
    "TotChol",
    "TG",
    "HDL",
    "LDL",
    "CRP",
    "Insulin",
]
for col in numeric_cols:
    clin[col] = pd.to_numeric(clin[col], errors="coerce")

# Binary diagnoses (0/1, treat empty as NaN)
binary_cols = ["T2D", "HTN", "Dyslip", "OSA", "Asthma/COPD", "LiverDx", "GORD"]
for col in binary_cols:
    clin[col] = pd.to_numeric(clin[col], errors="coerce")

# T2D group label
clin["T2D_label"] = clin["T2D"].map({1.0: "T2D", 0.0: "Non-T2D"})

print(f"Loaded {len(clin)} patients")
print(f"T2D: {clin['T2D'].value_counts().to_dict()}")
print(f"Sex: {clin['Sex'].value_counts().to_dict()}")

# ── 2. UPSET PLOT: COMORBIDITY CO-OCCURRENCE ─────────────────
# Shows how diagnoses cluster together in the cohort.
# Each bar = patients with exactly that combination of diagnoses.

print("\nGenerating UpSet plot...")

COMORBIDITIES = ["T2D", "HTN", "Dyslip", "OSA", "Asthma/COPD", "GORD"]

# Restrict to patients with complete comorbidity data
upset_data = clin[COMORBIDITIES].dropna()
print(f"  Patients with complete comorbidity data: {len(upset_data)}")

# upsetplot expects a boolean MultiIndex
upset_bool = upset_data.astype(bool)
upset_series = upsetplot.from_indicators(indicators=COMORBIDITIES, data=upset_bool)

fig_upset = plt.figure(figsize=(14, 7))
upset_obj = upsetplot.UpSet(
    upset_series,
    subset_size="count",
    show_counts=True,
    min_subset_size=1,
    sort_by="cardinality",  # largest subsets first
    sort_categories_by="-cardinality",
    facecolor="#4292C6",
    shading_color=0.08,
    totals_plot_elements=3,
)
upset_obj.style_subsets(
    present="T2D", absent="OSA", facecolor="#E6550D", label="T2D without OSA"
)
upset_obj.plot(fig=fig_upset)
fig_upset.suptitle(
    "Comorbidity co-occurrence in the bariatric surgery cohort",
    y=1.01,
    fontsize=13,
    fontweight="bold",
)
fig_upset.savefig(OUT_DIR / "01_comorbidity_upset.png", dpi=200, bbox_inches="tight")
fig_upset.savefig(OUT_DIR / "01_comorbidity_upset.pdf", bbox_inches="tight")
plt.close()

# Count table
combo_counts = (
    upset_bool.apply(lambda row: tuple(row[COMORBIDITIES]), axis=1)
    .value_counts()
    .reset_index()
)
combo_counts.columns = ["combination", "count"]
combo_counts.to_csv(OUT_DIR / "01_comorbidity_combinations.csv", index=False)

# ── 3. DISTRIBUTIONS OF CONTINUOUS VARIABLES BY T2D ─────────
# Violin + strip plots for each continuous variable split by T2D status.
# Annotated with Mann-Whitney U p-values.

print("Plotting continuous variable distributions by T2D...")

CONT_VARS = [
    ("Age", "Age (years)"),
    ("BMI", "BMI (kg/m²)"),
    ("Glucose", "Fasting Glucose (mmol/L)"),
    ("HbA1c", "HbA1c (mmol/mol)"),
    ("TotChol", "Total Cholesterol (mmol/L)"),
    ("TG", "Triglycerides (mmol/L)"),
    ("HDL", "HDL Cholesterol (mmol/L)"),
    ("LDL", "LDL Cholesterol (mmol/L)"),
    ("CRP", "CRP (mg/L)"),
    ("Insulin", "Fasting Insulin (pmol/L)"),
]

ncols = 5
nrows = int(np.ceil(len(CONT_VARS) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.5, nrows * 4))
axes = axes.flatten()

for ax, (col, label) in zip(axes, CONT_VARS):
    sub = clin[["T2D_label", col]].dropna()
    sns.violinplot(
        data=sub,
        x="T2D_label",
        y=col,
        ax=ax,
        palette=PALETTE_T2D,
        order=["Non-T2D", "T2D"],
        inner=None,
        alpha=0.6,
        linewidth=1,
    )
    sns.stripplot(
        data=sub,
        x="T2D_label",
        y=col,
        ax=ax,
        palette=PALETTE_T2D,
        order=["Non-T2D", "T2D"],
        size=4,
        jitter=True,
        alpha=0.8,
        linewidth=0.3,
        edgecolor="white",
    )
    # Mann-Whitney test
    g1 = sub.loc[sub["T2D_label"] == "Non-T2D", col]
    g2 = sub.loc[sub["T2D_label"] == "T2D", col]
    if len(g1) >= 3 and len(g2) >= 3:
        stat, pval = mannwhitneyu(g1, g2, alternative="two-sided")
        ptext = f"p={pval:.3f}" if pval >= 0.001 else f"p<0.001"
        star = (
            "***"
            if pval < 0.001
            else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
        )
        ymax = sub[col].max()
        ax.text(
            0.5,
            0.97,
            f"{star} ({ptext})",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=8.5,
            color="black",
        )
    ax.set_title(label, fontsize=9.5, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(labelsize=8)

for ax in axes[len(CONT_VARS) :]:
    ax.set_visible(False)

fig.suptitle(
    "Clinical variables: T2D vs Non-T2D", fontsize=13, fontweight="bold", y=1.01
)
plt.tight_layout()
fig.savefig(OUT_DIR / "02_distributions_by_T2D.png")
fig.savefig(OUT_DIR / "02_distributions_by_T2D.pdf")
plt.close()

# ── 4. LOGISTIC REGRESSION: PREDICTING EACH DIAGNOSIS ────────
# For each binary outcome (T2D, HTN, Dyslip, OSA, Asthma/COPD),
# fit a logistic regression using basic clinical predictors.
# Reports: OR, 95% CI, p-value, AUC (LOO cross-validated).
# This forms the clinical baseline — MOFA factors should explain
# variance BEYOND what these simple variables already capture.

print("Fitting logistic regression models...")

PREDICTORS = ["Age", "Sex_F", "BMI", "HbA1c", "Glucose", "TotChol", "TG", "HDL", "CRP"]

OUTCOMES = {
    "T2D": "Type 2 Diabetes",
    "HTN": "Hypertension",
    "Dyslip": "Dyslipidaemia",
    "OSA": "Obstructive Sleep Apnoea",
    "Asthma/COPD": "Asthma / COPD",
    "GORD": "GORD",
}

results_all = []
auc_summary = {}
coef_matrix = pd.DataFrame(index=PREDICTORS)

for outcome, outcome_label in OUTCOMES.items():
    cols_needed = [outcome] + PREDICTORS
    sub = clin[cols_needed].dropna()

    if sub[outcome].nunique() < 2:
        print(f"  Skipping {outcome}: only one class present")
        continue
    if len(sub) < 15:
        print(f"  Skipping {outcome}: n={len(sub)} too small")
        continue

    y = sub[outcome].astype(int)
    X_raw = sub[PREDICTORS]

    # ─ Univariate logistic regression for OR table ────────────
    for pred in PREDICTORS:
        sub_uni = sub[[outcome, pred]].dropna()
        if len(sub_uni) < 10 or sub_uni[pred].nunique() < 2:
            continue
        y_uni = sub_uni[outcome].astype(int)
        x_uni = sm.add_constant(sub_uni[pred])
        try:
            m = sm.Logit(y_uni, x_uni).fit(disp=False, maxiter=100)
            coef = m.params[pred]
            ci_lo = m.conf_int()[0][pred]
            ci_hi = m.conf_int()[1][pred]
            pval = m.pvalues[pred]
            results_all.append(
                {
                    "Outcome": outcome_label,
                    "Predictor": pred,
                    "OR": np.exp(coef),
                    "OR_lo": np.exp(ci_lo),
                    "OR_hi": np.exp(ci_hi),
                    "p_value": pval,
                    "n": len(sub_uni),
                }
            )
        except Exception:
            pass

    # ─ Multivariate model + LOO-AUC ──────────────────────────
    sub_multi = sub[[outcome] + PREDICTORS].dropna()
    if len(sub_multi) < 15:
        continue
    y_multi = sub_multi[outcome].astype(int).values
    X_multi = sub_multi[PREDICTORS].values

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_multi)

    clf = LogisticRegression(penalty="l2", C=0.5, max_iter=500, random_state=42)

    # LOO AUC (appropriate for small cohorts)
    loo = LeaveOneOut()
    y_prob = np.zeros(len(y_multi))
    for train_idx, test_idx in loo.split(X_scaled):
        clf.fit(X_scaled[train_idx], y_multi[train_idx])
        y_prob[test_idx] = clf.predict_proba(X_scaled[test_idx])[:, 1]

    try:
        loo_auc = roc_auc_score(y_multi, y_prob)
    except Exception:
        loo_auc = np.nan

    auc_summary[outcome_label] = loo_auc
    print(
        f"  {outcome_label}: LOO-AUC = {loo_auc:.3f}  (n={len(sub_multi)}, "
        f"cases={y_multi.sum()})"
    )
    coef_matrix.loc[PREDICTORS, outcome] = clf.coef_[0]

# ── 5. FOREST PLOT: ORs PER OUTCOME ─────────────────────────

print("Plotting forest plots...")

results_df = pd.DataFrame(results_all)

# FDR correction per outcome
results_df["p_fdr"] = np.nan
for outcome in results_df["Outcome"].unique():
    idx = results_df["Outcome"] == outcome
    _, p_corr, _, _ = multipletests(results_df.loc[idx, "p_value"], method="fdr_bh")
    results_df.loc[idx, "p_fdr"] = p_corr

results_df["sig"] = results_df["p_fdr"].apply(
    lambda p: "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
)
results_df.to_csv(OUT_DIR / "03_logistic_regression_ORs.csv", index=False)

# Forest plot — one panel per outcome
outcomes_plot = [o for o in OUTCOMES.values() if o in results_df["Outcome"].values]
n_outcomes = len(outcomes_plot)
fig, axes = plt.subplots(1, n_outcomes, figsize=(n_outcomes * 4, 6), sharey=False)
if n_outcomes == 1:
    axes = [axes]

PRED_LABELS = {
    "Age": "Age",
    "Sex_F": "Sex (F)",
    "BMI": "BMI",
    "HbA1c": "HbA1c",
    "Glucose": "Glucose",
    "TotChol": "Total Chol",
    "TG": "TG",
    "HDL": "HDL",
    "CRP": "CRP",
}

for ax, outcome in zip(axes, outcomes_plot):
    sub = results_df[results_df["Outcome"] == outcome].copy()
    sub["pred_label"] = sub["Predictor"].map(PRED_LABELS)
    sub = sub.sort_values("OR")
    y_pos = range(len(sub))

    colors = ["#E6550D" if p < 0.05 else "#AAAAAA" for p in sub["p_fdr"]]
    # errorbar doesn't accept a list for markerfacecolor — draw CI lines and
    # markers separately so each point can have its own colour
    ax.errorbar(
        sub["OR"],
        y_pos,
        xerr=[sub["OR"] - sub["OR_lo"], sub["OR_hi"] - sub["OR"]],
        fmt="none",
        capsize=3,
        ecolor="dimgrey",
        elinewidth=1.2,
    )
    for xi, yi, ci in zip(sub["OR"], y_pos, colors):
        ax.plot(
            xi,
            yi,
            "o",
            color=ci,
            markeredgecolor="dimgrey",
            markersize=7,
            markeredgewidth=0.6,
        )
    ax.axvline(1, color="black", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(sub["pred_label"].tolist(), fontsize=9)
    ax.set_xlabel("Odds Ratio (95% CI)", fontsize=9)
    ax.set_title(outcome, fontsize=10, fontweight="bold", wrap=True)
    ax.tick_params(axis="x", labelsize=8)

    # Annotate significance
    for i, (_, row) in enumerate(sub.iterrows()):
        if row["sig"]:
            ax.text(
                row["OR_hi"] * 1.02,
                i,
                row["sig"],
                va="center",
                ha="left",
                fontsize=8,
                color="#E6550D",
            )

fig.suptitle(
    "Univariate logistic regression: odds ratios for each diagnosis\n"
    "(FDR-corrected; coloured = FDR<0.05)",
    fontsize=12,
    fontweight="bold",
)
plt.tight_layout()
fig.savefig(OUT_DIR / "03_forest_plots_ORs.png")
fig.savefig(OUT_DIR / "03_forest_plots_ORs.pdf")
plt.close()

# ── 6. AUC BAR CHART ─────────────────────────────────────────
# Shows how well basic clinical features predict each diagnosis.
# MOFA factors should improve on these AUCs.

if auc_summary:
    fig, ax = plt.subplots(figsize=(8, 4))
    labels = list(auc_summary.keys())
    aucs = list(auc_summary.values())
    colors_bar = ["#E6550D" if l == "Type 2 Diabetes" else "#3182BD" for l in labels]
    bars = ax.barh(labels, aucs, color=colors_bar, edgecolor="white", height=0.55)
    ax.axvline(
        0.5,
        color="grey",
        linestyle="--",
        linewidth=0.9,
        alpha=0.7,
        label="Random (AUC=0.5)",
    )
    ax.axvline(
        0.7,
        color="#31A354",
        linestyle=":",
        linewidth=0.9,
        alpha=0.7,
        label="Good (AUC=0.7)",
    )
    for bar, auc in zip(bars, aucs):
        ax.text(
            auc + 0.01,
            bar.get_y() + bar.get_height() / 2,
            f"{auc:.2f}",
            va="center",
            ha="left",
            fontsize=9,
        )
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("LOO cross-validated AUC", fontsize=11)
    ax.set_title(
        "Predictive performance of basic clinical variables\n"
        "(Age, Sex, BMI, HbA1c, Glucose, Cholesterols, CRP)\n"
        "→ Baseline to compare MOFA factor models against",
        fontsize=11,
        fontweight="bold",
    )
    ax.legend(fontsize=8)
    plt.tight_layout()
    fig.savefig(OUT_DIR / "04_AUC_clinical_baseline.png")
    fig.savefig(OUT_DIR / "04_AUC_clinical_baseline.pdf")
    plt.close()

# ── 7. ASSOCIATION HEATMAP: DIAGNOSES × CONTINUOUS VARS ──────
# Spearman correlation of each continuous variable with each binary
# diagnosis. Colour = Spearman rho; annotated with FDR-corrected p-values.

print("Computing association heatmap...")

assoc_rho = pd.DataFrame(index=OUTCOMES.keys(), columns=PREDICTORS, dtype=float)
assoc_pval = pd.DataFrame(index=OUTCOMES.keys(), columns=PREDICTORS, dtype=float)

for outcome in OUTCOMES:
    for pred in PREDICTORS:
        sub = clin[[outcome, pred]].dropna()
        sub[outcome] = pd.to_numeric(sub[outcome], errors="coerce")
        sub = sub.dropna()
        if len(sub) < 8 or sub[outcome].nunique() < 2:
            continue
        rho, pval = spearmanr(sub[pred], sub[outcome])
        assoc_rho.loc[outcome, pred] = rho
        assoc_pval.loc[outcome, pred] = pval

# FDR correction across all tests
flat_pvals = assoc_pval.values.flatten()
valid_mask = ~np.isnan(flat_pvals.astype(float))
flat_pvals_valid = flat_pvals[valid_mask].astype(float)
_, p_fdr, _, _ = multipletests(flat_pvals_valid, method="fdr_bh")
flat_fdr = np.full_like(flat_pvals, np.nan, dtype=float)
flat_fdr[valid_mask] = p_fdr
assoc_fdr = pd.DataFrame(
    flat_fdr.reshape(assoc_pval.shape),
    index=assoc_pval.index,
    columns=assoc_pval.columns,
)

# Annotation: * for FDR<0.05, ** for FDR<0.01
annot = assoc_fdr.copy().astype(str)
for i in range(assoc_fdr.shape[0]):
    for j in range(assoc_fdr.shape[1]):
        p = assoc_fdr.iloc[i, j]
        if pd.isna(p):
            annot.iloc[i, j] = ""
        elif p < 0.01:
            annot.iloc[i, j] = "**"
        elif p < 0.05:
            annot.iloc[i, j] = "*"
        else:
            annot.iloc[i, j] = ""

row_labels = [OUTCOMES[o] for o in assoc_rho.index]
col_labels = [PRED_LABELS.get(p, p) for p in assoc_rho.columns]

fig, ax = plt.subplots(figsize=(12, 5))
sns.heatmap(
    assoc_rho.astype(float),
    ax=ax,
    annot=annot,
    fmt="s",
    cmap="RdBu_r",
    center=0,
    vmin=-0.8,
    vmax=0.8,
    linewidths=0.4,
    linecolor="white",
    cbar_kws={"label": "Spearman ρ", "shrink": 0.7},
    xticklabels=col_labels,
    yticklabels=row_labels,
    annot_kws={"size": 11, "weight": "bold"},
)
ax.set_title(
    "Association of clinical variables with each diagnosis\n"
    "Spearman ρ  (* FDR<0.05  ** FDR<0.01)",
    fontsize=11,
    fontweight="bold",
)
ax.set_xticklabels(ax.get_xticklabels(), rotation=35, ha="right", fontsize=9)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=9)
plt.tight_layout()
fig.savefig(OUT_DIR / "05_association_heatmap.png")
fig.savefig(OUT_DIR / "05_association_heatmap.pdf")
plt.close()

# ── 8. PAIRPLOT OF KEY METABOLIC VARIABLES ────────────────────
# Scatter matrix of BMI, HbA1c, Glucose, TG, HDL, CRP coloured by T2D.
# Gives a quick gestalt of metabolic clustering.

print("Plotting pairplot...")

pairplot_vars = ["BMI", "HbA1c", "Glucose", "TG", "HDL", "CRP"]
pp_data = clin[pairplot_vars + ["T2D_label"]].dropna(subset=["T2D_label"])

g = sns.pairplot(
    pp_data.dropna(),
    hue="T2D_label",
    palette=PALETTE_T2D,
    hue_order=["Non-T2D", "T2D"],
    diag_kind="kde",
    plot_kws={"alpha": 0.6, "s": 40, "edgecolor": "white", "linewidth": 0.3},
    diag_kws={"fill": True, "alpha": 0.4},
)
g.figure.suptitle(
    "Pairplot of key metabolic variables by T2D status",
    y=1.01,
    fontsize=12,
    fontweight="bold",
)
g.figure.savefig(OUT_DIR / "06_pairplot_metabolic.png")
g.figure.savefig(OUT_DIR / "06_pairplot_metabolic.pdf")
plt.close()

# ── 9. SEX × T2D × COMORBIDITY SUMMARY TABLE ─────────────────

print("Generating summary table...")

summary_rows = []
for outcome in OUTCOMES:
    if outcome == "T2D":
        continue  # skip T2D as an outcome since it's the main stratifier
    sub = clin[[outcome, "T2D", "Sex"]].dropna(subset=[outcome])
    sub[outcome] = pd.to_numeric(sub[outcome], errors="coerce")
    sub["T2D"] = pd.to_numeric(sub["T2D"], errors="coerce")
    sub = sub.dropna()

    n_total = len(sub)
    n_pos = sub[outcome].sum()
    pct_pos = 100 * n_pos / n_total if n_total > 0 else np.nan

    # Stratified by T2D
    for t2d_val, label in [(0.0, "Non-T2D"), (1.0, "T2D")]:
        grp = sub[sub["T2D"] == t2d_val]
        n_grp = len(grp)
        n_grp_pos = grp[outcome].sum() if n_grp > 0 else 0
        pct_grp = 100 * n_grp_pos / n_grp if n_grp > 0 else np.nan
        summary_rows.append(
            {
                "Comorbidity": OUTCOMES[outcome],
                "T2D_status": label,
                "N": n_grp,
                "N_positive": int(n_grp_pos),
                "Pct_positive": round(pct_grp, 1),
            }
        )

    # Chi-square T2D vs diagnosis
    if sub["T2D"].nunique() == 2 and sub[outcome].nunique() == 2:
        ct = pd.crosstab(sub["T2D"], sub[outcome])
        if ct.shape == (2, 2):
            chi2, pchi, _, _ = chi2_contingency(ct, correction=True)
            summary_rows[-1]["Chi2_pval_T2D"] = round(pchi, 4)
            summary_rows[-2]["Chi2_pval_T2D"] = round(pchi, 4)

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(OUT_DIR / "07_comorbidity_T2D_summary.csv", index=False)

# ── 10. PRINT SUMMARY ────────────────────────────────────────

print("\n" + "=" * 60)
print("CLINICAL EXPLORATION COMPLETE")
print("=" * 60)
print(f"Outputs in: {OUT_DIR}")
print()
print("Files generated:")
for f in sorted(OUT_DIR.glob("*")):
    print(f"  {f.name}")
print()
print("AUC summary (LOO cross-validated, multivariate model):")
print("  — baseline for comparison with MOFA factor models —")
for outcome, auc in auc_summary.items():
    print(f"  {outcome:<35} AUC = {auc:.3f}")
print()
print("Key notes:")
print("  • HbA1c values are in mmol/mol (not %)")
print("    (e.g. 48 mmol/mol ≈ 6.5% = T2D threshold)")
print("  • OSA prevalence is very high in this cohort (~80%)")
print("    — low discriminating power; interpret OSA analyses cautiously")
print("  • LiverDx has only ~4 cases — excluded from regression models")
