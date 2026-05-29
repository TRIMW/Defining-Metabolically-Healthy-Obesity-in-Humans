import os

import pandas as pd
import scanpy as sc
import pertpy as pt

cell_type_col = "cell_type2_"

results_per_ct = {}
res_list = []
skipped = []


adata_pb = sc.read_h5ad(
    f"/Users/willtrim/Documents/projs/Wills_analysis/obesity_scRNAseq/inputs/anndatas/sq_reannotated_finely_pb_{cell_type_col}.h5ad"
)


for ct in sorted(adata_pb.obs[cell_type_col].unique()):
    pdata_ct = adata_pb[adata_pb.obs[cell_type_col] == ct].copy()

    # sc.pp.filter_genes(pdata_ct, min_cells=1)
    sc.pp.filter_genes(pdata_ct, min_counts=10)

    # Need both cohort levels and enough samples for sex+cohort design
    cohorts_present = pdata_ct.obs["cohort"].nunique()
    sexes_present = pdata_ct.obs["sex"].nunique()
    n = pdata_ct.n_obs

    if cohorts_present < 2:
        skipped.append((ct, f"n={n}, only one cohort"))
        continue
    if n < 6:
        skipped.append((ct, f"n={n}, too few samples for sex+cohort design"))
        continue

    design = "~sex+cohort" if sexes_present > 1 else "~cohort"

    # try:
    pds2_ct = pt.tl.PyDESeq2(adata=pdata_ct, design=design)
    pds2_ct.fit()
    res_ct = pds2_ct.test_contrasts(
        pds2_ct.contrast(
            column="cohort", baseline="Healthy", group_to_compare="Unhealthy"
        )
    )
    res_ct["cell_type"] = ct
    results_per_ct[ct] = res_ct
    print(
        f"  {ct:35s}  n={n:3d}  design={design}  top gene: {res_ct.iloc[0]['variable']}  padj={res_ct.iloc[0]['adj_p_value']:.2e}"
    )
    res_list.append(res_ct)
    # except Exception as e:
    # skipped.append((ct, str(e)))

# print("\nSkipped:")
for ct, reason in skipped:
    print(f"  {ct}: {reason}")

res_all = pd.concat(res_list, axis=0)
res_all.to_csv(
    f"obesity_scRNAseq/results/DE/subq/deseq2_results_per_{cell_type_col}.csv"
)
