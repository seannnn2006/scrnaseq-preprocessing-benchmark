"""
08c_refined_comparison.py

08b把ASPCs拆成20個subcluster後,人工核對每個subcluster裡PanSci標註
(main_cell_type)的組成,找到幾個乾淨獨立、能對應到08原本抓不到的稀有
類型的subcluster:

    subcluster 18+19 (共2822個細胞) -> 95.8%/99.0%純的Rdh16 positive cells
    subcluster 11    (1068個細胞)   -> 92.6%純的Neural cells
    subcluster 10    (1063個細胞)   -> 91.3%純的Mesothelial cells
    subcluster 14    (108個細胞)    -> 76.9%的Epididymal cells(細胞數少,標tentative)

subcluster 17(Lymphatic/Vascular endothelial混在一起)跟15(Neural/Fcgbp
positive混在一起)不夠乾淨,沒有獨立標籤,維持ASPCs不變。其餘12個
subcluster都是99%左右純的ASPCs,也維持不變。

這支腳本把上面4組標籤合併回主資料,產生新的cell_type_refined欄位
(9種cluster層級類型,ASPCs細分出4種之後變成12種),再重新畫一次
similarity matrix,看之前全紫色的欄位是不是解決了。

Input:  results/atlas_aging/wt_aging_final_annotated.h5ad (07的輸出)
        results/atlas_aging/subclusters/subcluster_ASPCs.h5ad (08b的輸出)
Output: results/atlas_aging/wt_aging_refined_annotated.h5ad
        results/atlas_aging/annotation_comparison_refined.xlsx
        results/atlas_aging/figures/celltype_similarity_matrix_refined.png
        results/atlas_aging/figures/umap_celltype_refined.png
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
plt.rcParams.update({
    "font.size": 14, "axes.titlesize": 18, "axes.labelsize": 14, "legend.fontsize": 9,
})

# ASPCs的subcluster編號(字串,08b存的subcluster欄位是category字串)
# -> 新標籤。沒列在這裡的subcluster維持ASPCs不變。
ASPC_SUBCLUSTER_LABELS = {
    "18": "Rdh16+ epithelial-like cells",
    "19": "Rdh16+ epithelial-like cells",
    "11": "Neural cells",
    "10": "Mesothelial cells",
    "14": "Epididymal cells (tentative)",
}


def main():
    paths = get_project_paths(__file__)
    out_dir = os.path.join(paths["script_dir"], "..", "..", "results", "atlas_aging")
    figures_dir = os.path.join(out_dir, "figures")
    subclusters_dir = os.path.join(out_dir, "subclusters")

    main_path = os.path.join(out_dir, "wt_aging_final_annotated.h5ad")
    print(f"Loading: {main_path}")
    adata = sc.read_h5ad(main_path)

    sub_path = os.path.join(subclusters_dir, "subcluster_ASPCs.h5ad")
    print(f"Loading: {sub_path}")
    sub = sc.read_h5ad(sub_path, backed="r")

    # 每個ASPCs細胞(用barcode對應)-> subcluster編號 -> 新標籤(沒對到的維持None)
    subcluster_series = sub.obs["subcluster"].astype(str)
    new_label_for_aspc_cell = subcluster_series.map(ASPC_SUBCLUSTER_LABELS)  # 對不到的變NaN

    adata.obs["cell_type_refined"] = adata.obs["cell_type_final"].astype(str)
    # 用barcode(obs_names)對齊,只更新ASPCs裡有新標籤的那些細胞
    relabel_map = new_label_for_aspc_cell.dropna()
    matched = adata.obs_names.isin(relabel_map.index)
    adata.obs.loc[matched, "cell_type_refined"] = relabel_map.reindex(adata.obs_names[matched]).values

    print("\ncell_type_refined分布(ASPCs細分後):")
    print(adata.obs["cell_type_refined"].value_counts())

    output_path = os.path.join(out_dir, "wt_aging_refined_annotated.h5ad")
    adata.write_h5ad(output_path)
    print(f"\nSaved: {output_path}")

    fig = sc.pl.umap(
        adata, color="cell_type_refined", show=False, return_fig=True,
        size=8, frameon=False, title="3mo+16mo WT atlas: cell type (ASPCs refined)",
        legend_fontsize=9,
    )
    fig_path = os.path.join(figures_dir, "umap_celltype_refined.png")
    fig.savefig(fig_path, dpi=POSTER_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")

    # -------------------------------------------------------------------
    # 重新畫similarity matrix(邏輯跟08一樣),看紫色缺口是否補上
    # -------------------------------------------------------------------
    counts = pd.crosstab(adata.obs["cell_type_refined"], adata.obs["main_cell_type"])
    proportions = counts.div(counts.sum(axis=1), axis=0)

    best_match = proportions.idxmax(axis=1)
    best_match_frac = proportions.max(axis=1)
    summary = pd.DataFrame({
        "cell_type_refined": counts.index,
        "n_cells": counts.sum(axis=1).values,
        "best_matching_pansci_label": best_match.values,
        "best_match_fraction": best_match_frac.values,
    }).sort_values("best_match_fraction")

    print("\n各cell_type_refined跟PanSci最佳匹配的比例:")
    print(summary.to_string(index=False))

    excel_path = os.path.join(out_dir, "annotation_comparison_refined.xlsx")
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
        cbar_kws={"label": "Proportion of cell_type_refined's cells"},
    )
    g.ax_heatmap.set_xlabel("PanSci original annotation (main_cell_type)")
    g.ax_heatmap.set_ylabel("Our refined annotation (cell_type_refined)")
    g.fig.suptitle("Cell-type similarity matrix (ASPCs refined): ours vs. PanSci", y=1.02)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")

    fig_path = os.path.join(figures_dir, "celltype_similarity_matrix_refined.png")
    g.savefig(fig_path, dpi=POSTER_DPI, bbox_inches="tight")
    plt.close(g.fig)
    print(f"Saved: {fig_path}")


if __name__ == "__main__":
    main()
