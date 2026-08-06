"""
08i_final_umap_v2_adipocytes.py

老師要求再對Adipocytes做一次subclustering(邏輯跟08b對ASPCs/Macrophages
做的一樣)。08b跑完後(見subcluster_Adipocytes.h5ad),人工核對每個subcluster
裡PanSci原始標註(main_cell_type)的組成,發現:

    subcluster 1,4,5,6,7,8 (共94,764個細胞) -> 97.3~99.1%純的Brown adipocytes
    subcluster 15 (279個細胞)               -> 70.6%的Vascular endothelial cells,
                                                併入既有的Endothelial cells
                                                (跟08e併Macrophages subcluster 8
                                                進DCs用同一套邏輯:比例過半、
                                                對應到既有類別就直接併入)
    subcluster 0 (7,456個細胞)  -> 52.8% Brown adipocytes但混了24.4%骨骼肌,
                                    不夠乾淨,維持Adipocytes不變
    subcluster 2 (13,551個細胞) -> 49.2% Adipocytes vs 49.1% Brown adipocytes,
                                    幾乎五五對分,無法判斷,維持Adipocytes不變
    其餘subcluster(3,9,10,11,12,13,14)本來就是97~99.9%純的Adipocytes,不變

這支腳本接在08e之後,在08e已經套用過的ASPCs+Macrophages修正基礎上,
再疊加這次Adipocytes的修正,直接覆蓋08e的輸出(wt_aging_v2_annotated.h5ad、
umap_celltype_v2_final.png),並且順便重算annotation_comparison_v2.xlsx
(邏輯跟08_compare_annotations.py一樣,只是改用cell_type_v2)——這樣
08f/08g兩支下游腳本不用改,直接重跑就會讀到含Brown adipocytes的新結果。

Input:  results/atlas_aging/wt_aging_final_annotated.h5ad (07的輸出)
        results/atlas_aging/subclusters/subcluster_ASPCs.h5ad
        results/atlas_aging/subclusters/subcluster_Macrophages.h5ad
        results/atlas_aging/subclusters/subcluster_Adipocytes.h5ad
Output: results/atlas_aging/wt_aging_v2_annotated.h5ad (覆蓋08e的版本)
        results/atlas_aging/figures/umap_celltype_v2_final.png (覆蓋)
        results/atlas_aging/annotation_comparison_v2.xlsx (覆蓋)
"""

import os
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from load_utils import get_project_paths

POSTER_DPI = 300
plt.rcParams.update({
    "font.size": 14, "axes.titlesize": 18, "axes.labelsize": 14, "legend.fontsize": 9,
})

ASPC_SUBCLUSTER_LABELS = {
    "18": "Rdh16+ epithelial-like cells",
    "19": "Rdh16+ epithelial-like cells",
    "11": "Neural cells",
    "10": "Mesothelial cells",
    "14": "Epididymal cells (tentative)",
    "17": "Lymphatic endothelial cells (tentative)",
}

MACROPHAGE_SUBCLUSTER_LABELS = {
    "8": "DCs",
}

# 只列出跟原本"Adipocytes"標籤不同的subcluster;沒列到的(3,9,10,11,12,13,14)
# 本來就是97%以上純Adipocytes,維持不變。0跟2太混雜(見上面docstring),
# 也刻意不列進來,維持Adipocytes不變。
ADIPOCYTE_SUBCLUSTER_LABELS = {
    "1": "Brown adipocytes",
    "4": "Brown adipocytes",
    "5": "Brown adipocytes",
    "6": "Brown adipocytes",
    "7": "Brown adipocytes",
    "8": "Brown adipocytes",
    "15": "Endothelial cells",
}


def apply_subcluster_labels(adata, cell_type_final_value, subcluster_h5ad_path, label_map):
    """讀一個subcluster h5ad,依label_map把符合的細胞的cell_type_v2改成新標籤"""
    sub = sc.read_h5ad(subcluster_h5ad_path, backed="r")
    subcluster_series = sub.obs["subcluster"].astype(str)
    new_label = subcluster_series.map(label_map).dropna()
    matched = adata.obs_names.isin(new_label.index)
    adata.obs.loc[matched, "cell_type_v2"] = new_label.reindex(adata.obs_names[matched]).values
    print(f"  {cell_type_final_value}: 更新了{matched.sum()}個細胞的標籤")


def main():
    paths = get_project_paths(__file__)
    out_dir = os.path.join(paths["script_dir"], "..", "..", "results", "atlas_aging")
    figures_dir = os.path.join(out_dir, "figures")
    subclusters_dir = os.path.join(out_dir, "subclusters")

    main_path = os.path.join(out_dir, "wt_aging_final_annotated.h5ad")
    print(f"Loading: {main_path}")
    adata = sc.read_h5ad(main_path)

    adata.obs["cell_type_v2"] = adata.obs["cell_type_final"].astype(str)

    print("\n套用sub-clustering人工核對後的修正:")
    apply_subcluster_labels(
        adata, "ASPCs",
        os.path.join(subclusters_dir, "subcluster_ASPCs.h5ad"),
        ASPC_SUBCLUSTER_LABELS,
    )
    apply_subcluster_labels(
        adata, "Macrophages",
        os.path.join(subclusters_dir, "subcluster_Macrophages.h5ad"),
        MACROPHAGE_SUBCLUSTER_LABELS,
    )
    apply_subcluster_labels(
        adata, "Adipocytes",
        os.path.join(subclusters_dir, "subcluster_Adipocytes.h5ad"),
        ADIPOCYTE_SUBCLUSTER_LABELS,
    )

    print("\n最終cell_type_v2分布:")
    print(adata.obs["cell_type_v2"].value_counts())

    output_path = os.path.join(out_dir, "wt_aging_v2_annotated.h5ad")
    adata.write_h5ad(output_path)
    print(f"\nSaved: {output_path}")

    fig = sc.pl.umap(
        adata, color="cell_type_v2", show=False, return_fig=True,
        size=8, frameon=False,
        title="3mo+16mo WT atlas: cell type (v2, ASPCs+Macrophages+Adipocytes refined)",
        legend_fontsize=9,
    )
    fig_path = os.path.join(figures_dir, "umap_celltype_v2_final.png")
    fig.savefig(fig_path, dpi=POSTER_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")

    # -------------------------------------------------------------------
    # 重算annotation_comparison_v2.xlsx,邏輯跟08_compare_annotations.py
    # 一樣,只是改用cell_type_v2,這樣08f/08g不用改就能讀到新結果
    # -------------------------------------------------------------------
    counts = pd.crosstab(adata.obs["cell_type_v2"], adata.obs["main_cell_type"])
    proportions = counts.div(counts.sum(axis=1), axis=0)

    best_match = proportions.idxmax(axis=1)
    best_match_frac = proportions.max(axis=1)
    summary = pd.DataFrame({
        "cell_type_v2": counts.index,
        "n_cells": counts.sum(axis=1).values,
        "best_matching_pansci_label": best_match.values,
        "best_match_fraction": best_match_frac.values,
    }).sort_values("best_match_fraction")

    print("\n各cell_type_v2跟PanSci最佳匹配的比例(由低到高排序):")
    print(summary.to_string(index=False))

    excel_path = os.path.join(out_dir, "annotation_comparison_v2.xlsx")
    with pd.ExcelWriter(excel_path) as writer:
        summary.to_excel(writer, sheet_name="best_match_summary", index=False)
        counts.to_excel(writer, sheet_name="raw_counts")
        proportions.to_excel(writer, sheet_name="row_proportions")
    print(f"Saved: {excel_path}")


if __name__ == "__main__":
    main()