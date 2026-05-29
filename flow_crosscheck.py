"""
flow_crosscheck.py
==================
For each "X % of Y" column in the percentages file:
  - Parse child (X) and parent (Y) from the column name
  - Find exactly-matching columns in the counts file by name
    (stripping tissue prefix and trailing " #" / " Count" / " ATM" etc.)
  - If both found: compute Spearman r between stated % and child/parent*100
  - Report what was found and what was not

Columns with "% of Tot" (Total cells — not in counts file) are skipped.

Output: results/flow_crosscheck/reconstruction_results.csv
"""

import warnings

warnings.filterwarnings("ignore")

import re
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
from pathlib import Path

# ── Config ─────────────────────────────────────────────────────
PCT_PATH = "data/flow_cytometry_perct_nonstim.csv"
CNT_PATH = "data/flow_cytometry_counts_nonstim.csv"
OUT = Path("results/flow_crosscheck")
OUT.mkdir(parents=True, exist_ok=True)

MIN_VALID = 8
R_WARN = 0.90
TISSUES = ["SubQ", "Blood", "Omentum"]

# ── Load ───────────────────────────────────────────────────────
pct = pd.read_csv(PCT_PATH, index_col=0, encoding="utf-8-sig")
cnt = pd.read_csv(CNT_PATH, index_col=0, encoding="utf-8-sig")
pct.index = pct.index.str.strip()
cnt.index = cnt.index.str.strip()

common = pct.index.intersection(cnt.index)
print(f"Patients — pct: {len(pct)}  cnt: {len(cnt)}  common: {len(common)}")
pct = pct.loc[common].apply(pd.to_numeric, errors="coerce")
cnt = cnt.loc[common].apply(pd.to_numeric, errors="coerce")


# ── Normalise a counts column name to a bare population token ──
# Strips: tissue prefix, trailing " #", "COUNT"/"Count", " ATM", "% of ..." , extra spaces
def normalise_cnt(col, tissue):
    s = col
    if s.startswith(tissue + "_"):
        s = s[len(tissue) + 1 :]
    s = re.sub(r"\s*#\s*$", "", s)
    s = re.sub(r"\s*Count\s*$", "", s, flags=re.IGNORECASE)
    s = re.sub(r"\s*%\s+of\s+.*$", "", s, flags=re.IGNORECASE)
    s = re.sub(r"\s*%\s*.*$", "", s, flags=re.IGNORECASE)
    s1 = s.strip().lower()
    s = re.sub(r"\s+ATM\s*$", "", s)
    s2 = s.strip().lower()
    return s1, s2


# Build lookup: tissue → {normalised_name → original_col}
cnt_lookup = {}
for tissue in TISSUES:
    cnt_lookup[tissue] = {}
    for col in cnt.columns:
        if col.startswith(tissue + "_"):
            keys = normalise_cnt(col, tissue)
            for key in keys:
                # keep first occurrence if duplicates (make.unique may have added suffix)
                if key not in cnt_lookup[tissue]:
                    cnt_lookup[tissue][key] = col

CD45_TOKENS = {"cd45", "cd45+", "cd45.1"}
CD45_COL = "CD45+ of Total"  # bare name without tissue prefix


def find_child_col(child_token, parent_token, tissue):
    """
    Find the count column for the child population.
    Strategy:
      1. CD45 parent → special column (handled at parent level, not here)
      2. Exact match of child_token after normalisation
      3. If not found, try "child_token parent_token" (e.g. "CD9+" + "ATM" → "CD9+ ATM")
         because count columns often append the gate name to the marker.
    Returns (col_name_or_None, method_string).
    """
    key = child_token.strip().lower()
    if key in cnt_lookup[tissue]:
        return cnt_lookup[tissue][key], "exact"
    # try child + parent appended
    key2 = (child_token + " " + parent_token).strip().lower()
    if key2 in cnt_lookup[tissue]:
        return cnt_lookup[tissue][key2], "child+parent"
    return None, ""


def complement_token(token):
    """
    Return the +/- complement of a marker token.
    e.g. "CD9+" → "CD9-",  "CD206-" → "CD206+".
    Swaps only the trailing + or -.
    """
    t = token.strip()
    if t.endswith("+"):
        return t[:-1] + "-"
    if t.endswith("-"):
        return t[:-1] + "+"
    return None  # no +/- to swap


def find_parent_col_or_complement(child_token, child_col, parent_token, tissue):
    """
    Find or synthesise the parent population count.
    Strategy:
      1. CD45 tokens → special column "CD45+ of Total".
      2. Exact match of parent_token after normalisation.
      3. Complement synthesis: take the child count column name, substitute
         the child marker with its +/- complement to get the sibling column,
         then sum child + sibling to get the parent total.
         e.g. child "CD9+ ATM #" → sibling "CD9- ATM #" → total ATM.
    Returns (array_or_None, col_description_string).
    """
    pt_lower = parent_token.strip().lower()

    # Rule 1: CD45 special case
    if pt_lower in CD45_TOKENS:
        full = f"{tissue}_{CD45_COL}"
        if full in cnt.columns:
            return cnt[full].values.astype(float), full
        return None, f"{full} (NOT FOUND)"

    # Rule 2: exact match
    if pt_lower in cnt_lookup[tissue]:
        col = cnt_lookup[tissue][pt_lower]
        return cnt[col].values.astype(float), col

    # Rule 3a: complement synthesis using the child count column
    # Works for simple markers: "CD9+ ATM" → complement "CD9- ATM"
    if child_col is not None:
        comp = complement_token(child_token)
        if comp is not None:
            child_norm = normalise_cnt(child_col, tissue)[0]
            comp_norm = child_norm.replace(child_token.strip().lower(), comp.lower(), 1)
            sibling_col = cnt_lookup[tissue].get(comp_norm)
            if sibling_col:
                c1 = cnt[child_col].values.astype(float)
                c2 = cnt[sibling_col].values.astype(float)
                arr = c1 + c2
                arr[np.isnan(c1) & np.isnan(c2)] = np.nan
                return arr, f"SUM({child_col} + {sibling_col})"

    # Rule 3b: fallback for multi-marker children (e.g. "CD63+CD9+ % of ATM")
    # Find any X+/X- complementary pair in the counts file where both columns
    # end with " {parent_token}", and sum that pair as the parent total.
    suffix = " " + pt_lower  # e.g. " atm"
    parent_cands = {
        key: col for key, col in cnt_lookup[tissue].items() if key.endswith(suffix)
    }
    for key, col in parent_cands.items():
        marker = key[: -len(suffix)]  # strip trailing " parent"
        comp = complement_token(marker)
        if comp is None:
            continue
        comp_key = comp.lower() + suffix
        if comp_key in parent_cands:
            c1 = cnt[col].values.astype(float)
            c2 = cnt[parent_cands[comp_key]].values.astype(float)
            arr = c1 + c2
            arr[np.isnan(c1) & np.isnan(c2)] = np.nan
            return arr, f"SUM({col} + {parent_cands[comp_key]}) [fallback]"

    return None, ""


# ── Parse and check ────────────────────────────────────────────
parse_re = re.compile(r"^(.+?)\s+%\s+of\s+(.+)$", re.IGNORECASE)

rows = []
for pct_col in pct.columns:
    tissue = next((t for t in TISSUES if pct_col.startswith(t + "_")), None)
    if tissue is None:
        continue

    short = pct_col[len(tissue) + 1 :]
    m = parse_re.match(short)
    if not m:
        continue

    child_token = m.group(1).strip()
    parent_token = m.group(2).strip()

    # Skip "% of Tot" — total cell counts not in counts file
    if parent_token.lower() in ("tot", "total"):
        continue

    child_col, child_method = find_child_col(child_token, parent_token, tissue)
    parent_arr, parent_desc = find_parent_col_or_complement(
        child_token, child_col, parent_token, tissue
    )

    r = n = None
    status = ""

    if child_col and parent_arr is not None:
        y_stated = pct[pct_col].values.astype(float)
        num = cnt[child_col].values.astype(float)
        den = parent_arr
        with np.errstate(divide="ignore", invalid="ignore"):
            recon = np.where((den > 0) & ~np.isnan(den), num / den * 100, np.nan)
        mask = ~(np.isnan(y_stated) | np.isnan(recon))
        n = int(mask.sum())
        if n >= MIN_VALID:
            r, _ = spearmanr(y_stated[mask], recon[mask])
            r = round(float(r), 4)
            status = "OK" if r >= R_WARN else "MISMATCH"
        else:
            status = "TOO_FEW_PAIRS"
    elif not child_col and parent_arr is None:
        status = "BOTH_NOT_FOUND"
    elif not child_col:
        status = "CHILD_NOT_FOUND"
    else:
        status = "PARENT_NOT_FOUND"

    rows.append(
        {
            "tissue": tissue,
            "pct_col": pct_col,
            "child_token": child_token,
            "parent_token": parent_token,
            "child_cnt_col": child_col or "",
            "child_method": child_method,  # "exact" or "child+parent"
            "parent_desc": parent_desc,  # col name, SUM(...), or ""
            "n_valid": n if n is not None else "",
            "spearman_r": r if r is not None else "",
            "status": status,
        }
    )

df = pd.DataFrame(rows)
df.to_csv(OUT / "reconstruction_results.csv", index=False)

# ── Print summary ───────────────────────────────────────────────
print(f"\nTotal '% of' columns checked : {len(df)}")
for status, grp in df.groupby("status"):
    print(f"  {status:<22} {len(grp)}")

print("\nMismatches (r < {:.2f}):".format(R_WARN))
mm = df[df["status"] == "MISMATCH"]
if mm.empty:
    print("  None")
else:
    print(
        mm[["pct_col", "child_cnt_col", "parent_desc", "spearman_r"]].to_string(
            index=False
        )
    )

print("\nNot found:")
nf = df[df["status"].str.contains("NOT_FOUND")]
if nf.empty:
    print("  None")
else:
    print(
        nf[["pct_col", "child_token", "parent_desc", "status"]].to_string(index=False)
    )

print("\nColumns resolved via child+parent name or SUM (for info):")
non_exact = df[
    df["child_method"].isin(["child+parent"])
    | df["parent_desc"].str.startswith("SUM", na=False)
]
if non_exact.empty:
    print("  None")
else:
    print(
        non_exact[
            ["pct_col", "child_method", "parent_desc", "spearman_r", "status"]
        ].to_string(index=False)
    )

print(f"\nSaved: {OUT}/reconstruction_results.csv")
