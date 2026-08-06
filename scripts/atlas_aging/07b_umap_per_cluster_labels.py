"""
07b_umap_per_cluster_labels.py

老師的要求:07最終的UMAP(umap_celltype_final.png)把同一種cell type的所有
cluster都合併成同一個顏色,如果一種類型底下其實有好幾個cluster(例如
Adipocytes有8個cluster),合併著色會看不出這些cluster之間的差異。

這支腳本改成每個Leiden cluster各自一個顏色/標籤,但標籤不是單純的cluster
編號(0、1、2...),而是「cell type + 序號」,例如Adipocytes底下的8個cluster
分別叫Adipocytes 1、Adipocytes 2...Adipocytes 8,序號依leiden編號由小到大排。
只有1個cluster的類型(例如Pericytes)也一樣編號成"Pericytes 1",全部類型
用同一套規則,不特別區分。

這張新圖不會覆蓋掉07原本存的umap_celltype_final.png,存成新檔名;另外多存
一張左右並排對照圖,方便直接比較「合併後」vs.「拆開後」的差異。

Input:  results/atlas_aging/wt_aging_final_annotated.h5ad (07的輸出)
Output: results/atlas_aging/figures/umap_celltype_final_per_cluster.png (新圖,23個cluster各自標籤)
        results/atlas_aging/figures/umap_celltype_final_side_by_side.png (兩張並排對照)
"""

import os
from collections import defaultdict
import scanpy as sc
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from load_utils import get_project_paths

POSTER_DPI = 300
plt.rcParams.update({
    "font.size": 14, "axes.titlesize": 18, "axes.labelsize": 14, "legend.fontsize": 9,
})


def main():
    paths = get_project_paths(__file__)
    out_dir = os.path.join(paths["script_dir"], "..", "..", "results", "atlas_aging")
    figures_dir = os.path.join(out_dir, "figures")

    input_path = os.path.join(out_dir, "wt_aging_final_annotated.h5ad")
    print(f"Loading: {input_path}")
    adata = sc.read_h5ad(input_path)

    # -------------------------------------------------------------------
    # 每個leiden cluster -> "cell_type_final 序號"。序號依leiden編號由小到大
    # 排,同一個cell_type_final底下才重新從1開始編號。
    # -------------------------------------------------------------------
    cluster_to_celltype = adata.obs.drop_duplicates("leiden").set_index("leiden")["cell_type_final"].to_dict()

    celltype_to_clusters = defaultdict(list)
    for cl, ct in cluster_to_celltype.items():
        celltype_to_clusters[ct].append(cl)

    cluster_to_label = {}
    for ct, clusters in celltype_to_clusters.items():
        for i, cl in enumerate(sorted(clusters, key=lambda x: int(x)), start=1):
            cluster_to_label[cl] = f"{ct} {i}"

    adata.obs["cell_type_final_per_cluster"] = (
        adata.obs["leiden"].map(cluster_to_label).astype(str)
    )

    # 圖例順序:依cell_type_final字母排序,同一類型底下再依序號排,
    # 不然預設會照字串排序把"Adipocytes 10"排到"Adipocytes 2"前面
    # (這裡最多到8,暫時不會出現這個問題,但排序邏輯還是照嚴謹的方式寫)
    ordered_labels = []
    for ct in sorted(celltype_to_clusters.keys()):
        n = len(celltype_to_clusters[ct])
        ordered_labels.extend([f"{ct} {i}" for i in range(1, n + 1)])
    adata.obs["cell_type_final_per_cluster"] = (
        adata.obs["cell_type_final_per_cluster"].astype("category").cat.set_categories(ordered_labels)
    )

    print("\n每個leiden cluster對應的新標籤:")
    for cl in sorted(cluster_to_label, key=lambda x: int(x)):
        print(f"  cluster {cl}: {cluster_to_label[cl]}")

    # -------------------------------------------------------------------
    # 新圖1:23個cluster各自一個標籤/顏色(單獨存檔,不覆蓋umap_celltype_final.png)
    # -------------------------------------------------------------------
    fig = sc.pl.umap(
        adata, color="cell_type_final_per_cluster", show=False, return_fig=True,
        size=8, frameon=False,
        title="3mo+16mo WT atlas: cell type (final, per-cluster labels)",
        legend_fontsize=9,
    )
    fig_path = os.path.join(figures_dir, "umap_celltype_final_per_cluster.png")
    fig.savefig(fig_path, dpi=POSTER_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")

    # -------------------------------------------------------------------
    # 新圖2:跟07原本合併著色的版本並排比較
    # -------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(24, 10))
    sc.pl.umap(
        adata, color="cell_type_final", show=False, ax=axes[0],
        size=8, frameon=False, legend_fontsize=8,
        title="Merged by cell type (9 types)",
    )
    sc.pl.umap(
        adata, color="cell_type_final_per_cluster", show=False, ax=axes[1],
        size=8, frameon=False, legend_fontsize=8,
        title="Split by cluster (23 clusters, labelled)",
    )
    fig.suptitle("3mo+16mo WT atlas: merged vs. per-cluster cell type labels", y=1.02, fontsize=18)
    fig.tight_layout()
    fig_path = os.path.join(figures_dir, "umap_celltype_final_side_by_side.png")
    fig.savefig(fig_path, dpi=POSTER_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")


if __name__ == "__main__":
    main()
