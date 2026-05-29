# COHORTS.py
# ============================================================
# SINGLE SOURCE OF TRUTH — patient cohort definitions
# Generated from Patient_samples.xlsx (verified March 2026)
# Import this file in every analysis script.
#
# scRNAseq column is noted here but IGNORED in all analyses.
# Those 10 patients are fully included under every other
# dataset they contributed to (flow, clinical bloods, etc.)
# The scRNAseq data itself is not used anywhere in this pipeline.
# ============================================================

# ── All 59 working patients (P47 excluded — sex only, no data) ──
ALL_PATIENTS = [
    "P01","P02","P03","P04","P05","P06","P07","P08","P09","P10",
    "P11","P12","P13","P14","P15","P16","P17","P18","P19","P20",
    "P21","P22","P23","P24","P25","P26","P27","P28","P29","P30",
    "P31","P32","P33","P34","P35","P36","P37","P38","P39","P40",
    "P41","P42","P43","P44","P45","P46","P48","P49","P50","P51",
    "P52","P53","P54","P55","P56","P57","P58","P59","P60",
]

# ── Cohort A: Clinical characterisation (n=52) ───────────────
# ≥3 clinical blood measures available
# Data: HbA1c(50), Glucose(45), Lipids(51-52), CRP(42), Insulin(26)
#       T2D(20), HTN(19), Dyslip(18), OSA(37), Liver(4)
#       ASA/OS-MRS(46), sex known for 50/52
COHORT_A = [
    "P02","P04","P05","P06","P07","P08","P09","P10","P11","P12",
    "P13","P14","P15","P16","P17","P18","P19","P20","P21","P22",
    "P23","P24","P25","P26","P27","P28","P30","P31","P32","P34",
    "P35","P36","P37","P38","P39","P40","P41","P42","P43","P44",
    "P45","P46","P48","P49","P50","P51","P52","P53","P54","P55",
    "P56","P57",
]  # n=52

# ── Cohort B: Core immune phenotyping (n=41) ─────────────────
# Blood Flow + Omentum Flow (both required)
# Data: all clinical bloods ~80-93%, diagnoses ~29-68%
#       SubQ flow(34), stim flows(18-31), multiplex(26)
#       H&E(16), RNAseq(10), sex known for all 41
COHORT_B = [
    "P01","P02","P03","P04","P05","P06","P07","P08","P09","P10",
    "P11","P12","P13","P14","P15","P16","P17","P18","P19","P20",
    "P21","P22","P23","P24","P25","P26","P27","P28","P29","P30",
    "P31","P32","P33","P34","P35","P37","P38","P39","P40","P41",
    "P42",
]  # n=41

# ── Cohort B3: Triple-tissue flow (n=34) ─────────────────────
# Blood Flow + OM Flow + SubQ Flow (all three required)
# Primary cohort for unsupervised stratification
# Data: HbA1c(33/34), clinical bloods ~85-97%, diagnoses~35-74%
#       stim flows(17-27), multiplex(22), H&E(15), RNAseq(10)
#       sex known for all 34
COHORT_B3 = [
    "P03","P04","P05","P06","P07","P08","P09","P12","P13","P14",
    "P15","P17","P18","P19","P20","P21","P22","P23","P24","P25",
    "P26","P27","P28","P30","P31","P32","P33","P34","P35","P38",
    "P39","P40","P41","P42",
]  # n=34

# ── Cohort B_stim: All 3 stimulated flows (n=17) ─────────────
# Blood Stim + OM Stim + SubQ Stim (all three required)
# For within-patient cross-tissue functional comparisons
COHORT_B_STIM = [
    "P04","P05","P06","P08","P12","P13","P17","P20","P21","P22",
    "P23","P24","P26","P28","P30","P31","P33",
]  # n=17

# ── Cohort C: Bulk RNAseq deep phenotyping (n=10) ────────────
# Omental + SubQ + liver bulk RNAseq
# ALL 10 also have: all 3 flow panels, HbA1c, blood stim(8/10)
#                   OM stim(7/10), H&E(8/10), most clinical bloods
# NOTE: n=10, hypothesis-generating for multivariate; use
#       descriptive stats + effect sizes, report p-values with caution
COHORT_C = [
    "P19","P20","P22","P25","P27","P30","P33","P34","P35","P40",
]  # n=10

# ── H&E Histology cohort (n=17) ──────────────────────────────
COHORT_HE = [
    "P19","P20","P21","P23","P24","P27","P28","P29","P30","P31",
    "P33","P34","P35","P36","P38","P39","P40",
]  # n=17

# ── TLR stimulation cohort (n=10) ────────────────────────────
# Omental adipose TLR stimulations
# 8/10 also have flow (B), 3/10 also have RNAseq (C)
COHORT_TLR = [
    "P12","P20","P27","P29","P35","P36","P37","P38","P42","P60",
]  # n=10

# ── LPS/PolyIC stimulation cohort (n=10) ─────────────────────
# Omental adipose LPS & PolyIC stimulations
# Only 2/10 overlap with flow (P39, P40)
# 8/10 have multiplex cytokines; 8/10 have clinical bloods
COHORT_LPS = [
    "P39","P40","P44","P45","P46","P52","P53","P57","P58","P59",
]  # n=10

# ── Pre-computed sub-cohort intersections ─────────────────────
# RNAseq integrations
C_x_B        = [p for p in COHORT_C if p in COHORT_B]        # n=10
C_x_B3       = [p for p in COHORT_C if p in COHORT_B3]       # n=10
C_x_HE       = [p for p in COHORT_C if p in COHORT_HE]       # n=8
C_x_OM_STIM  = ["P19","P20","P22","P25","P27","P30","P33"]    # n=7
C_x_SQ_STIM  = ["P20","P22","P30","P33"]                      # n=4
C_x_B_STIM   = ["P19","P20","P22","P25","P27","P30","P33","P34"] # n=8
C_x_MULTI    = ["P20","P22","P25","P40"]                       # n=4
C_x_TLR      = ["P20","P27","P35"]                             # n=3

# Flow integrations
B3_x_HE      = [p for p in COHORT_B3 if p in COHORT_HE]      # n=15
B3_x_MULTI   = [p for p in COHORT_B3 if p in
                ["P02","P04","P05","P06","P08","P09","P11","P12",
                 "P13","P15","P16","P17","P18","P20","P21","P22",
                 "P23","P24","P25","P26","P37","P38","P39","P40",
                 "P41","P42","P52","P53","P57","P59","P60",
                 "P45","P46","P48","P49"]]                     # n=22

# Stim pairings
STIM_B_OM    = ["P04","P05","P06","P07","P08","P09","P10","P11",
                "P12","P13","P14","P16","P17","P18","P19","P20",
                "P21","P22","P23","P24","P25","P26","P27","P28",
                "P29","P30","P31","P33"]                       # n=28 Blood+OM stim
STIM_B_SQ    = ["P04","P05","P06","P08","P12","P13","P17","P20",
                "P21","P22","P23","P24","P26","P28","P30","P31",
                "P32","P33"]                                   # n=18 Blood+SubQ stim
STIM_OM_SQ   = ["P04","P05","P06","P08","P12","P13","P17","P20",
                "P21","P22","P23","P24","P26","P28","P30","P31",
                "P33"]                                         # n=17 OM+SubQ stim

# ── Column name mapping (Excel → short code) ─────────────────
COL_SHORT = {
    "Blood Multiplex":                     "BMultiplex",
    "Omental Adipose TLR Stimulations":    "OM_TLR",
    "Omental Adipose LPS & PolyIC Stimulations": "OM_LPS",
    "Omental and SubQ adipose, and liver RNAseq": "RNAseq",
    "Adipose H&E Imaging":                 "HE",
    "Omentum adipose Flow":                "OM_Flow",
    "Blood Flow":                          "B_Flow",
    "SubQ Adipose Flow":                   "SQ_Flow",
    "SubQ Adipose Stimulated flow":        "SQ_Stim",
    "Blood Stimulated flow":               "B_Stim",
    "Omentum adipose Stimulated flow":     "OM_Stim",
    "Blood Glucose":                       "Glucose",
    "Blood Total Cholesterol":             "TotChol",
    "Blood Triglyceride":                  "TG",
    "Blood HDL Cholesterol":               "HDL",
    "Blood non-HDL Cholesterol":           "nonHDL",
    "Blood LDL Cholesterol":               "LDL",
    "Blood CRP":                           "CRP",
    "Blood Insulin":                       "Insulin",
    "Blood HbA1c":                         "HbA1c",
    "Type 2 Diabetes Diagnosis":           "T2D",
    "American Society of Anesthesiologists (ASA) Physical Status Classification System": "ASA",
    "Obesity Surgery Mortality Risk Score (OS-MRS)": "OSMRS",
    "Hypertension  Diagnosis":             "HTN",
    "Dyslipidemia Diagnosis":              "Dyslip",
    "Obstructive sleep apneoa Diagnosis":  "OSA",
    "Liver Disease Diagnosis":             "LiverDx",
    "Female": "Female",
    "Male":   "Male",
}

# ── Clinical blood variable groups ───────────────────────────
CLINICAL_BLOODS = ["HbA1c","Glucose","TotChol","TG","HDL","LDL","CRP","Insulin","nonHDL"]
GLYCAEMIA       = ["HbA1c","Glucose","Insulin"]
LIPIDS          = ["TotChol","TG","HDL","LDL","nonHDL"]
INFLAMMATION    = ["CRP"]
DIAGNOSES       = ["T2D","HTN","Dyslip","OSA"]  # LiverDx n=4, descriptive only

# ── Statistical guidance ──────────────────────────────────────
# Liver Disease: n=4 total, 0 in Cohort C — DESCRIPTIVE ONLY
# Cohort C n=10: hypothesis-generating; no standalone p-values
# Sex: fully recorded in B, B3, C; 2 unknown in A (P33, P43)
# Primary clinical stratifier: T2D (n=17-20 across cohorts)
# Secondary stratifiers: HTN, Dyslip, OSA (sufficient n in B/B3)
STAT_GUIDANCE = {
    "A":  "n=52. Parametric tests valid for continuous variables. "
          "FDR correction for multiple comparisons. Sex as covariate.",
    "B":  "n=41. UMAP/clustering valid. Non-parametric group comparisons. "
          "Spearman correlations. Sex fully recorded — use as covariate.",
    "B3": "n=34. Primary stratification cohort. Sufficient for NMF, PSN, UMAP. "
          "T2D comparisons: 17 vs 17. Sex fully recorded.",
    "B_stim": "n=17. Paired within-patient tests only (Wilcoxon signed-rank). "
              "No multivariate. Effect sizes reported.",
    "C":  "n=10. HYPOTHESIS-GENERATING. Spearman r with confidence intervals. "
          "No p-values reported in isolation. Effect sizes + visualisation primary.",
    "TLR": "n=10. Descriptive + paired comparisons within patients only.",
    "LPS": "n=10. Descriptive. Only 2 overlap with flow — no flow integration.",
}

def cohort_list(cohort):
    """Return sorted list from set or list input."""
    return sorted(set(cohort), key=lambda x: int(x[1:]))
