"""
05_marker_annotation.py

用marker基因幫04算出來的Leiden分群結果貼上cell_type標籤,邏輯、marker清單
跟scripts/atlas_annotation/02_marker_annotation.py完全一樣(細胞類型的
marker基因不會因為年齡而改變,3個月跟16個月的Adipocytes都還是同一組
marker基因)。

每個cluster的判定方式:對cluster裡所有細胞的每種marker score取平均,
分數最高的那個細胞類型當作這個cluster的標籤。

Input:  results/atlas_aging/wt_aging_clustered.h5ad (04的輸出)
Output: results/atlas_aging/wt_aging_annotated.h5ad
        results/atlas_aging/cluster_marker_scores.xlsx (老師要看的表)
        results/atlas_aging/figures/umap_celltype.png
"""

import os
import scanpy as sc
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from load_utils import get_project_paths

POSTER_DPI = 300
plt.rcParams.update({
    "font.size": 14, "axes.titlesize": 18, "axes.labelsize": 14, "legend.fontsize": 10,
})

# 跟scripts/atlas_annotation/02_marker_annotation.py同一份清單(Loft et al.,
# Nature Metabolism 2025,表3),原封不動沿用——marker基因是細胞類型的生物學
# 定義,不會因為老鼠年齡而改變。
MARKER_GENES = {
    "Adipocytes": ["Adipoq", "Plin4"],
    "ASPCs": ["Pdgfra", "Dcn"],
    "Pericytes": ["Steap4", "Enpep"],
    "SMCs": ["Myocd", "Myh11", "Acta2"],
    "Endothelial cells": ["Pecam1", "Cdh5"],
    "Mesothelial cells": ["Msln", "Krt19"],
    "B cells": ["Ms4a1", "Cd79a", "Cd79b"],
    "T cells": ["Cd3d", "Cd3e"],
    "NK cells": ["Klrd1", "Xcl1"],
    "ILCs": ["Il7r", "Gata3", "Arg1"],
    "DCs": ["Flt3"],
    "Neutrophils": ["S100a8", "Csf3r", "S100a9"],
    "Monocytes": ["Plac8", "Lyz1", "Lyz2"],
    "Macrophages": ["Mrc1", "Adgre1"],
    "Mast cells": ["Kit", "Cpa3", "Mcpt4"],
    "Schwann cells": ["Mpz"],
}


def main():
    paths = get_project_paths(__file__)
    out_dir = os.path.join(paths["script_dir"], "..", "..", "results", "atlas_aging")
    figures_dir = os.path.join(out_dir, "figures")
    os.makedirs(figures_dir, exist_ok=True)

    input_path = os.path.join(out_dir, "wt_aging_clustered.h5ad")
    print(f"Loading: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"Input: {adata.shape[0]} cells x {adata.shape[1]} genes, "
          f"{adata.obs['leiden'].nunique()} Leiden clusters")

    adata_scoring = adata.raw.to_adata()

    gene_symbols = adata_scoring.var_names.str.split("|").str[-1]
    symbol_to_varname = dict(zip(gene_symbols, adata_scoring.var_names))

    print("\nScoring marker genes...")
    for label, genes in MARKER_GENES.items():
        genes_present = [symbol_to_varname[g] for g in genes if g in symbol_to_varname]
        if not genes_present:
            print(f"  Skipping '{label}': none of {genes} found in var_names.")
            continue
        if len(genes_present) < len(genes):
            missing = set(genes) - set(symbol_to_varname)
            print(f"  '{label}': only {genes_present} found (missing {missing}).")
        sc.tl.score_genes(adata_scoring, genes_present, score_name=f"score_{label}")

    score_cols = [c for c in adata_scoring.obs.columns if c.startswith("score_")]

    cluster_scores = adata_scoring.obs.groupby(adata.obs["leiden"].values)[score_cols].mean()
    cluster_scores.index.name = "leiden"
    cluster_scores.columns = [c.replace("score_", "") for c in cluster_scores.columns]
    cluster_scores["assigned_cell_type"] = cluster_scores.idxmax(axis=1)

    cluster_sizes = adata.obs["leiden"].value_counts()
    cluster_scores.insert(0, "n_cells", cluster_sizes.reindex(cluster_scores.index))

    print("\nCluster-level marker scores and assigned cell types:")
    print(cluster_scores)

    excel_path = os.path.join(out_dir, "cluster_marker_scores.xlsx")
    cluster_scores.to_excel(excel_path)
    print(f"\nSaved: {excel_path}")

    cluster_to_celltype = cluster_scores["assigned_cell_type"].to_dict()
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_to_celltype).astype(str)
    print(f"\nAssigned cell_type to all cells based on cluster-level marker scores:")
    print(adata.obs["cell_type"].value_counts())

    fig = sc.pl.umap(
        adata, color="cell_type", show=False, return_fig=True,
        size=8, frameon=False, title="3mo+16mo WT atlas: cell type (marker-based)",
        legend_fontsize=9,
    )
    fig_path = os.path.join(figures_dir, "umap_celltype.png")
    fig.savefig(fig_path, dpi=POSTER_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")

    output_path = os.path.join(out_dir, "wt_aging_annotated.h5ad")
    adata.write_h5ad(output_path)
    print(f"\nSaved annotated atlas to: {output_path}")


if __name__ == "__main__":
    main()