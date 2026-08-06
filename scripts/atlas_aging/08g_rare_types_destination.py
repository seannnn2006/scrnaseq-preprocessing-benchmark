"""
08g_rare_types_destination.py

08f的相似度矩陣是「依我們的類型」正規化(我們每一列各自加總100%),PanSci
裡4種真正稀有的類型(Fcgbp positive cells、Lymphoid cells_Plasma cells、
Erythroblasts、Circulating hepatoblasts)因為細胞數太少,稀釋在我們的大
cluster裡,那個方向看永遠會是接近0的紫色,不是分析沒做完整。

這支腳本反過來,只挑這4種稀有類型,改成「依PanSci稀有類型」的方向正規化
(每一列改成這4種稀有類型,各自加總100%),這樣才能看出它們自己的細胞
實際上被我們歸到哪裡去了,顏色不會是紫色一片。是08f主圖的補充說明,
不是取代主圖。

Input:  results/atlas_aging/wt_aging_v2_annotated.h5ad
Output: results/atlas_aging/figures/rare_pansci_types_destination.png
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
    "font.size": 13, "axes.titlesize": 16, "axes.labelsize": 13,
})

RARE_TYPES = [
    "Fcgbp positive cells",
    "Lymphoid cells_Plasma cells",
    "Erythroblasts",
    "Circulating hepatoblasts",
]


def main():
    paths = get_project_paths(__file__)
    out_dir = os.path.join(paths["script_dir"], "..", "..", "results", "atlas_aging")
    figures_dir = os.path.join(out_dir, "figures")

    input_path = os.path.join(out_dir, "wt_aging_v2_annotated.h5ad")
    print(f"Loading: {input_path}")
    adata = sc.read_h5ad(input_path, backed="r")

    obs = adata.obs[adata.obs["main_cell_type"].isin(RARE_TYPES)][["main_cell_type", "cell_type_v2"]]
    counts = pd.crosstab(obs["main_cell_type"], obs["cell_type_v2"])
    counts = counts.reindex(RARE_TYPES)

    # 只留有出現過的欄位,不然14種我們的類型裡大部分是0,圖會很空
    counts = counts.loc[:, (counts.sum(axis=0) > 0)]
    proportions = counts.div(counts.sum(axis=1), axis=0)

    n_cells = obs["main_cell_type"].value_counts().reindex(RARE_TYPES)
    row_labels = [f"{rt} (n={n_cells[rt]})" for rt in RARE_TYPES]
    proportions.index = row_labels
    counts.index = row_labels

    print("\n細胞數:")
    print(counts.to_string())
    print("\n比例:")
    print(proportions.round(3).to_string())

    # 欄數不固定(依實際有出現的cell_type_v2而定),寬度依欄數留夠空間給4位數
    # 的annot數字(例如1006)跟旋轉45度的x軸標籤,不然會像之前一樣互相疊字。
    # 高度也加大,標題才不會跟colorbar的刻度("1.0")擠在一起。
    n_cols = proportions.shape[1]
    fig, ax = plt.subplots(figsize=(max(11, 1.3 * n_cols), 6.5))
    sns.heatmap(
        proportions, cmap="viridis", vmin=0, vmax=1, annot=counts, fmt="d",
        annot_kws={"size": 11},
        linewidths=0.5, linecolor="white", ax=ax,
        cbar_kws={"label": "Proportion of this PanSci rare type's cells", "shrink": 0.8},
    )
    ax.set_xlabel("Our final annotation (cell_type_v2)")
    ax.set_ylabel("PanSci rare category\n(too few cells to form its own cluster)")
    ax.set_title("Where do these 4 rare PanSci categories end up in our classification?", pad=20)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    plt.setp(ax.get_yticklabels(), rotation=0)
    fig.tight_layout()

    fig_path = os.path.join(figures_dir, "rare_pansci_types_destination.png")
    fig.savefig(fig_path, dpi=POSTER_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {fig_path}")


if __name__ == "__main__":
    main()
