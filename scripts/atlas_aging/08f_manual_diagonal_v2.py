"""
08f_manual_diagonal_v2.py

跟08d邏輯一樣(手動排成對角線),但改用08e產生的最終版cell_type_v2
(ASPCs+Macrophages都修正過)重新畫一次。

Input:  results/atlas_aging/annotation_comparison_v2.xlsx
Output: results/atlas_aging/figures/celltype_similarity_matrix_v2_diagonal.png
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

ROW_ORDER = [
    "T cells", "B cells", "Macrophages", "DCs",
    "Pericytes", "Endothelial cells",
    "Lymphatic endothelial cells (tentative)",
    "Neural cells", "Mesothelial cells",
    "ASPCs", "Epididymal cells (tentative)",
    "Skeletal muscle cells", "Rdh16+ epithelial-like cells",
    "Adipocytes", "Brown adipocytes",
]
# "Endothelial cells (tentative)"(279個細胞,71%)故意不放進這張圖——它的
# 最佳匹配欄位(Vascular endothelial cells)跟上面"Endothelial cells"共用,
# 兩列搶同一欄會讓對角線在那個位置歪掉。這個類別數量少、訊號也不夠乾淨,
# 用文字/Excel另外說明就好,不用勉強塞進這張要求乾淨對角線的圖。

# PanSci欄位裡完全對不到任何一列(全部細胞比例都接近0、整欄呈現紫色)的
# 那幾個,老師覺得留著看起來不完整,直接從圖上移除。這4個已經在對話記錄
# 裡查過去向(Circulating hepatoblasts/Erythroblasts數量個位數,
# Fcgbp positive cells/Lymphoid cells_Plasma cells在任何subcluster裡都
# 沒有過半數訊號,go查詳見rare_pansci_types_destination.png這張補充圖)。
COLUMNS_TO_DROP = [
    "Circulating hepatoblasts",
    "Erythroblasts",
    "Fcgbp positive cells",
    "Lymphoid cells_Plasma cells",
]


def main():
    paths = get_project_paths(__file__)
    out_dir = os.path.join(paths["script_dir"], "..", "..", "results", "atlas_aging")
    figures_dir = os.path.join(out_dir, "figures")

    excel_path = os.path.join(out_dir, "annotation_comparison_v2.xlsx")
    print(f"Loading: {excel_path}")
    proportions = pd.read_excel(excel_path, sheet_name="row_proportions", index_col=0)
    summary = pd.read_excel(excel_path, sheet_name="best_match_summary")
    best_match = dict(zip(summary["cell_type_v2"], summary["best_matching_pansci_label"]))

    proportions = proportions.reindex(ROW_ORDER)

    # 有些列(例如細胞數很少的"Endothelial cells (tentative)")最佳匹配欄位
    # 跟另一列共用(都是Vascular endothelial cells),用dict.fromkeys去重,
    # 保留第一次出現的順序,不然同一欄會被選兩次、圖裡出現重複欄位
    matched_cols = list(dict.fromkeys(best_match[r] for r in ROW_ORDER))
    remaining_cols = [
        c for c in proportions.columns
        if c not in matched_cols and c not in COLUMNS_TO_DROP
    ]
    col_order = matched_cols + remaining_cols
    proportions = proportions[col_order]

    fig, ax = plt.subplots(figsize=(max(10, 0.55 * proportions.shape[1]), max(8, 0.5 * proportions.shape[0])))
    sns.heatmap(
        proportions, cmap="viridis", vmin=0, vmax=1,
        linewidths=0.5, linecolor="white", ax=ax,
        cbar_kws={"label": "Proportion of cells"},
    )
    ax.set_xlabel("PanSci original annotation")
    ax.set_ylabel("Our final cell-type annotation")
    ax.set_title("Validating our cell-type annotation against PanSci's original labels")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.tight_layout()

    fig_path = os.path.join(figures_dir, "celltype_similarity_matrix_v2_diagonal.png")
    fig.savefig(fig_path, dpi=POSTER_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")


if __name__ == "__main__":
    main()
