"""
06_de_marker_check.py

不靠marker清單,反過來用差異表現分析(Wilcoxon rank-sum test,one-vs-rest)
驗證05判定的cell_type合不合理,邏輯跟
scripts/atlas_annotation/03_de_marker_check.py完全一樣。

Input:  results/atlas_aging/wt_aging_annotated.h5ad (05的輸出)
Output: results/atlas_aging/cluster_top_de_genes.xlsx (老師要看的表)
"""

import os
import scanpy as sc

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from load_utils import get_project_paths

N_TOP_GENES = 20


def main():
    paths = get_project_paths(__file__)
    out_dir = os.path.join(paths["script_dir"], "..", "..", "results", "atlas_aging")

    input_path = os.path.join(out_dir, "wt_aging_annotated.h5ad")
    print(f"Loading: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"Input: {adata.shape[0]} cells x {adata.shape[1]} genes, "
          f"{adata.obs['leiden'].nunique()} Leiden clusters")

    cluster_to_celltype = adata.obs.drop_duplicates("leiden").set_index("leiden")["cell_type"].to_dict()

    print("\nRunning rank_genes_groups (Wilcoxon test, per Leiden cluster)...")
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", use_raw=True)

    de_df = sc.get.rank_genes_groups_df(adata, group=None)
    de_df["gene_symbol"] = de_df["names"].str.split("|").str[-1]
    de_df["assigned_cell_type"] = de_df["group"].map(cluster_to_celltype)

    top_de = de_df.groupby("group", sort=False).head(N_TOP_GENES).copy()
    top_de["rank"] = top_de.groupby("group").cumcount() + 1

    top_de = top_de[[
        "group", "assigned_cell_type", "rank", "gene_symbol",
        "scores", "logfoldchanges", "pvals_adj",
    ]].rename(columns={"group": "leiden"})

    print(f"\nTop {N_TOP_GENES} differentially expressed genes per cluster "
          f"(vs. 05's marker-based cell_type):")
    print(top_de.to_string(index=False))

    excel_path = os.path.join(out_dir, "cluster_top_de_genes.xlsx")
    top_de.to_excel(excel_path, index=False)
    print(f"\nSaved: {excel_path}")


if __name__ == "__main__":
    main()