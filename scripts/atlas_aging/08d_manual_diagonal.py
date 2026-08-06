"""
08d_manual_diagonal.py

08c的clustermap是用階層式聚類自動排序,雖然同一群顏色會排得比較靠近,
但不保證每一列的最佳匹配欄位剛好落在同一個對角線位置上(列跟欄是各自
獨立聚類的)。這支腳本改成手動指定順序:每一列(我們的類型)後面接的
第一欄,就是它在08c算出來的best_matching_pansci_label,強制排成真正
的左上到右下對角線,對不到任何一列的PanSci欄位放在最後面。

Input:  results/atlas_aging/annotation_comparison_refined.xlsx (08c的輸出)
Output: results/atlas_aging/figures/celltype_similarity_matrix_diagonal.png
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from load_utils import get_project_paths

POSTER_DPI = 300
plt.rcParams.update({
    "font.size": 14, "axes.titlesize": 18, "axes.labelsize": 14,
})

# 列的顯示順序(自訂,把生物學上相近的類型放在一起,免疫細胞一群、
# 血管/間質一群、脂肪細胞放最後,方便閱讀)
ROW_ORDER = [
    "T cells", "B cells", "Macrophages", "DCs",
    "Pericytes", "Endothelial cells", "Neural cells", "Mesothelial cells",
    "ASPCs", "Epididymal cells (tentative)",
    "Skeletal muscle cells", "Rdh16+ epithelial-like cells", "Adipocytes",
]


def main():
    paths = get_project_paths(__file__)
    out_dir = os.path.join(paths["script_dir"], "..", "..", "results", "atlas_aging")
    figures_dir = os.path.join(out_dir, "figures")

    excel_path = os.path.join(out_dir, "annotation_comparison_refined.xlsx")
    print(f"Loading: {excel_path}")
    proportions = pd.read_excel(excel_path, sheet_name="row_proportions", index_col=0)
    summary = pd.read_excel(excel_path, sheet_name="best_match_summary")
    best_match = dict(zip(summary["cell_type_refined"], summary["best_matching_pansci_label"]))

    proportions = proportions.reindex(ROW_ORDER)

    # 每一列的最佳匹配欄,依ROW_ORDER的順序排在最前面;沒被任何列選中的
    # PanSci欄位(例如Erythroblasts這種對不到的稀有類型)排在最後
    matched_cols = [best_match[r] for r in ROW_ORDER]
    remaining_cols = [c for c in proportions.columns if c not in matched_cols]
    col_order = matched_cols + remaining_cols
    proportions = proportions[col_order]

    fig, ax = plt.subplots(figsize=(max(10, 0.55 * proportions.shape[1]), max(8, 0.5 * proportions.shape[0])))
    sns.heatmap(
        proportions, cmap="viridis", vmin=0, vmax=1,
        linewidths=0.5, linecolor="white", ax=ax,
        cbar_kws={"label": "Proportion of cell_type_refined's cells"},
    )
    ax.set_xlabel("PanSci original annotation (main_cell_type)")
    ax.set_ylabel("Our refined annotation (cell_type_refined)")
    ax.set_title("Cell-type similarity matrix (manually ordered to diagonal)")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.tight_layout()

    fig_path = os.path.join(figures_dir, "celltype_similarity_matrix_diagonal.png")
    fig.savefig(fig_path, dpi=POSTER_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")


if __name__ == "__main__":
    main()
