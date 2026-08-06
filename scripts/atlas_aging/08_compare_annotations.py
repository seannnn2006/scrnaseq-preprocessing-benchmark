"""
08_compare_annotations.py

把07的最終標籤(cell_type_final)跟PanSci原始論文自帶的獨立標註
(main_cell_type)拿來對照,算cell-type similarity matrix,邏輯跟
scripts/atlas_annotation/04b_compare_annotations.py完全一樣(老師確認過
只做到cluster層級的最終標註、跟PanSci比對即可,不需要再往下做
sub-clustering)。

Input:  results/atlas_aging/wt_aging_final_annotated.h5ad (07的輸出)
Output: results/atlas_aging/annotation_comparison.xlsx (老師要看的交叉表)
        results/atlas_aging/figures/celltype_similarity_matrix.png
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from load_utils import get_project_paths

POSTER_DPI = 300


def main():
    paths = get_project_paths(__file__)
    out_dir = os.path.join(paths["script_dir"], "..", "..", "results", "atlas_aging")
    figures_dir = os.path.join(out_dir, "figures")

    input_path = os.path.join(out_dir, "wt_aging_final_annotated.h5ad")
    print(f"Loading: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"Input: {adata.shape[0]} cells, "
          f"{adata.obs['cell_type_final'].nunique()} cell_type_final (ours) x "
          f"{adata.obs['main_cell_type'].nunique()} main_cell_type (PanSci)")

    counts = pd.crosstab(adata.obs["cell_type_final"], adata.obs["main_cell_type"])
    proportions = counts.div(counts.sum(axis=1), axis=0)

    best_match = proportions.idxmax(axis=1)
    best_match_frac = proportions.max(axis=1)
    summary = pd.DataFrame({
        "cell_type_final": counts.index,
        "n_cells": counts.sum(axis=1).values,
        "best_matching_pansci_label": best_match.values,
        "best_match_fraction": best_match_frac.values,
    }).sort_values("best_match_fraction")

    print("\n各cell_type_final跟PanSci最佳匹配的比例(由低到高排序,越低越值得注意):")
    print(summary.to_string(index=False))

    excel_path = os.path.join(out_dir, "annotation_comparison.xlsx")
    with pd.ExcelWriter(excel_path) as writer:
        summary.to_excel(writer, sheet_name="best_match_summary", index=False)
        counts.to_excel(writer, sheet_name="raw_counts")
        proportions.to_excel(writer, sheet_name="row_proportions")
    print(f"\nSaved: {excel_path}")

    g = sns.clustermap(
        proportions, cmap="viridis", vmin=0, vmax=1,
        linewidths=0.5, linecolor="white",
        figsize=(max(10, 0.5 * proportions.shape[1]), max(8, 0.5 * proportions.shape[0])),
        dendrogram_ratio=0.12,
        cbar_kws={"label": "Proportion of cell_type_final's cells"},
    )
    g.ax_heatmap.set_xlabel("PanSci original annotation (main_cell_type)")
    g.ax_heatmap.set_ylabel("Our final annotation (cell_type_final)")
    g.fig.suptitle("Cell-type similarity matrix: ours (3mo+16mo WT) vs. PanSci", y=1.02)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")

    fig_path = os.path.join(figures_dir, "celltype_similarity_matrix.png")
    g.savefig(fig_path, dpi=POSTER_DPI, bbox_inches="tight")
    plt.close(g.fig)
    print(f"Saved: {fig_path}")


if __name__ == "__main__":
    main()