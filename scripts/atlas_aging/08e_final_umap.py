"""
08e_final_umap.py

整合兩輪sub-clustering的人工核對結果,產生最終版本的cell type標註跟UMAP:

1. ASPCs(08b/08c已經做過)-> 拆出Rdh16+ epithelial-like cells、Neural cells、
   Mesothelial cells、Epididymal cells (tentative)

2. Macrophages(這輪新做的)-> subcluster 8(1,985個細胞)過半數(62.7%)是
   PanSci標註的Myeloid cells_Dendritic cells,改標成DCs,併入原本的DCs類別

   Endothelial cells、Pericytes、T cells、DCs這輪也做了sub-clustering,
   但每個subcluster都跟原本的類型高度一致(95%以上都是同一個PanSci
   標籤),沒有發現新的、值得獨立標示的稀有類型,維持原樣不變。

   B cells的subcluster 0有17.4%是Lymphoid cells_Plasma cells,但主體
   仍是81.6%的B cells(沒有過半),不夠格獨立拆開,維持B cells不變
   (值得跟老師口頭提一下這個訊號,但不改標籤)。

Input:  results/atlas_aging/wt_aging_final_annotated.h5ad (07的輸出)
        results/atlas_aging/subclusters/subcluster_ASPCs.h5ad
        results/atlas_aging/subclusters/subcluster_Macrophages.h5ad
Output: results/atlas_aging/wt_aging_v2_annotated.h5ad
        results/atlas_aging/figures/umap_celltype_v2_final.png
"""

import os
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

ADIPOCYTE_SUBCLUSTER_LABELS = {
    "2": "Brown adipocytes",
    "6": "Brown adipocytes",
    "7": "Brown adipocytes",
    "8": "Brown adipocytes",
    "9": "Brown adipocytes",
    "10": "Brown adipocytes",
    "12": "Brown adipocytes",
    "11": "Endothelial cells (tentative)",
    # subcluster 0(46.7%Brown/29.7%骨骼肌)、1(56.7%Adipocytes/41.5%Brown)
    # 都沒有過半數的單一PanSci類型,不強行標籤,維持Adipocytes不變
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
        title="3mo+16mo WT atlas: cell type (v2, ASPCs+Macrophages refined)",
        legend_fontsize=9,
    )
    fig_path = os.path.join(figures_dir, "umap_celltype_v2_final.png")
    fig.savefig(fig_path, dpi=POSTER_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")


if __name__ == "__main__":
    main()
