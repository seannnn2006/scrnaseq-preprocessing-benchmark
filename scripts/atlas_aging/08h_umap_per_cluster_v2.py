"""
08h_umap_per_cluster_v2.py

跟07b同樣的需求(老師要求:同一種cell type如果底下有好幾個群,每群要各自
標一個號碼,不要合併成一個顏色,不然看不出群跟群之間的差異),但這次用
最新的v2結果(08e/08g,ASPCs跟Macrophages都做過sub-clustering後的版本)。

分組依據:
- 如果這個細胞原本(cell_type_final)屬於ASPCs或Macrophages,用它在
  08b算出來的subcluster編號分組(因為這兩種類型的細分是靠subcluster
  才分出來的,不是原本04算的leiden)
- 其他類型(Adipocytes、Endothelial cells、T cells、B cells、Pericytes、
  DCs、Skeletal muscle cells)維持用04算的leiden cluster編號分組

同一個cell_type_v2底下,依編號由小到大重新從1開始標號,例如Adipocytes
的8個leiden cluster變成Adipocytes 1~8;ASPCs細分後剩下的12個subcluster
變成ASPCs 1~12。

Input:  results/atlas_aging/wt_aging_v2_annotated.h5ad (08e的輸出)
        results/atlas_aging/subclusters/subcluster_ASPCs.h5ad
        results/atlas_aging/subclusters/subcluster_Macrophages.h5ad
Output: results/atlas_aging/figures/umap_celltype_v2_per_cluster.png
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
    "font.size": 14, "axes.titlesize": 18, "axes.labelsize": 14, "legend.fontsize": 8,
})

SUBCLUSTERED_TYPES = {
    "ASPCs": "subcluster_ASPCs.h5ad",
    "Macrophages": "subcluster_Macrophages.h5ad",
}

# 這幾個cell_type_v2的值,是這輪對ASPCs做sub-clustering才新發現、獨立標示
# 出來的稀有類型(08e/08g的成果),原本04/07的cluster層級判斷完全沒有這幾
# 個標籤。在圖上額外標註出來,跟老師強調「這些是子分群才找到的,不是原本
# 分群結果就有的」。
NEWLY_DISCOVERED_LABELS = [
    "Rdh16+ epithelial-like cells",
    "Neural cells",
    "Mesothelial cells",
    "Epididymal cells (tentative)",
    "Lymphatic endothelial cells (tentative)",
]


def main():
    paths = get_project_paths(__file__)
    out_dir = os.path.join(paths["script_dir"], "..", "..", "results", "atlas_aging")
    figures_dir = os.path.join(out_dir, "figures")
    subclusters_dir = os.path.join(out_dir, "subclusters")

    input_path = os.path.join(out_dir, "wt_aging_v2_annotated.h5ad")
    print(f"Loading: {input_path}")
    adata = sc.read_h5ad(input_path)

    # group_id: 每個細胞用來分組的編號,格式是"leiden_<n>"或
    # "<cell_type>_sub<n>",確保ASPCs/Macrophages的細分群不會跟其他類型的
    # leiden編號搞混
    group_id = "leiden_" + adata.obs["leiden"].astype(str)
    group_id = group_id.astype(object)

    for ct, fname in SUBCLUSTERED_TYPES.items():
        sub_path = os.path.join(subclusters_dir, fname)
        print(f"Loading: {sub_path}")
        sub = sc.read_h5ad(sub_path, backed="r")
        sub_ids = ct + "_sub" + sub.obs["subcluster"].astype(str)
        matched = adata.obs_names.isin(sub.obs_names)
        group_id.loc[matched] = sub_ids.reindex(adata.obs_names[matched]).values

    adata.obs["group_id"] = group_id

    # 每個cell_type_v2底下有哪些group_id,依leiden/subcluster編號由小到大
    # 排序後,重新從1開始編號
    celltype_to_groups = defaultdict(set)
    for gid, ct in zip(adata.obs["group_id"], adata.obs["cell_type_v2"]):
        celltype_to_groups[ct].add(gid)

    def sort_key(gid):
        # "leiden_9" -> 9; "ASPCs_sub12" -> 12,純粹取數字部分排序
        return int("".join(c for c in gid if c.isdigit()))

    group_to_label = {}
    for ct, groups in celltype_to_groups.items():
        for i, gid in enumerate(sorted(groups, key=sort_key), start=1):
            group_to_label[gid] = f"{ct} {i}"

    adata.obs["cell_type_v2_per_cluster"] = adata.obs["group_id"].map(group_to_label)

    ordered_labels = []
    for ct in sorted(celltype_to_groups.keys()):
        n = len(celltype_to_groups[ct])
        ordered_labels.extend([f"{ct} {i}" for i in range(1, n + 1)])
    adata.obs["cell_type_v2_per_cluster"] = (
        adata.obs["cell_type_v2_per_cluster"].astype("category").cat.set_categories(ordered_labels)
    )

    print(f"\n總共{len(ordered_labels)}個分群標籤:")
    for lbl in ordered_labels:
        print(f"  {lbl}")

    fig = sc.pl.umap(
        adata, color="cell_type_v2_per_cluster", show=False, return_fig=True,
        size=8, frameon=False,
        title="3mo+16mo WT atlas: cell type (v2, per-cluster labels)",
        legend_fontsize=8,
    )

    # -------------------------------------------------------------------
    # 額外標註:這5種是這輪對ASPCs做sub-clustering才新發現的稀有類型,
    # 04/07的cluster層級分析完全沒有這幾個標籤。在圖上用文字+箭頭指出它們
    # 的位置,並加上"(子分群新發現)"的說明,跟老師強調這是sub-clustering
    # 額外做出來的成果,不是原本的分群結果就有的。
    # -------------------------------------------------------------------
    ax = fig.axes[0]
    umap_coords = adata.obsm["X_umap"]
    x_min, x_max = umap_coords[:, 0].min(), umap_coords[:, 0].max()
    y_min, y_max = umap_coords[:, 1].min(), umap_coords[:, 1].max()

    # 把5個標註框固定排在圖的最上方,依序水平排開,每個各自的x位置錯開,
    # 全部框的y都在圖的正上方(同一水平線),彼此之間不會疊在一起;
    # 再各自用箭頭指回它實際的cluster中心點
    present_labels = [
        label for label in NEWLY_DISCOVERED_LABELS
        if (adata.obs["cell_type_v2"] == label).sum() > 0
    ]
    import textwrap
    n_labels = len(present_labels)
    text_y = y_max + 0.55 * (y_max - y_min)
    slot_w = (x_max - x_min) / n_labels
    for i, label in enumerate(present_labels):
        mask = adata.obs["cell_type_v2"] == label
        cx, cy = umap_coords[mask.values, 0].mean(), umap_coords[mask.values, 1].mean()
        tx = x_min + (i + 0.5) * slot_w
        wrapped_label = textwrap.fill(label, width=14)
        ax.annotate(
            f"{wrapped_label}\n(new: sub-cluster)",
            xy=(cx, cy), xytext=(tx, text_y),
            fontsize=6.5, color="black", fontweight="bold", ha="center", va="bottom",
            bbox=dict(boxstyle="round,pad=0.25", fc="yellow", ec="black", alpha=0.9),
            arrowprops=dict(arrowstyle="->", color="black", lw=1),
        )
    ax.set_ylim(top=text_y + 0.55 * (y_max - y_min))
    fig.set_size_inches(fig.get_size_inches()[0], fig.get_size_inches()[1] * 1.35)

    fig_path = os.path.join(figures_dir, "umap_celltype_v2_per_cluster.png")
    fig.savefig(fig_path, dpi=POSTER_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {fig_path}")


if __name__ == "__main__":
    main()