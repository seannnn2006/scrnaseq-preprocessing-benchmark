"""
03b_filter_doublets.py

套用doublet_score門檻、濾掉doublet,秒級操作,邏輯跟
scripts/atlas/05b_filter_doublets.py完全一樣。門檻先沿用同一個0.2
(等03實際跑完、看過doublet_score_distribution.png的分布後,可能需要
依這份資料實際的百分位數重新調整,不用重跑03/scrublet)。

Input:  results/atlas_aging/qc_filtered_predoublet.h5ad
Output: results/atlas_aging/qc_filtered_atlas.h5ad
"""

import os
import scanpy as sc

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from load_utils import get_project_paths

DOUBLET_SCORE_THRESHOLD = 0.2


def main():
    paths = get_project_paths(__file__)
    atlas_dir = os.path.join(paths["script_dir"], "..", "..", "results", "atlas_aging")

    input_path = os.path.join(atlas_dir, "qc_filtered_predoublet.h5ad")
    print(f"Loading: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"Before doublet filtering: {adata.shape[0]} cells x {adata.shape[1]} genes")

    if DOUBLET_SCORE_THRESHOLD is None:
        is_doublet = adata.obs["predicted_doublet"]
        print("使用scrublet自動判斷的predicted_doublet欄位。")
    else:
        is_doublet = adata.obs["doublet_score"] > DOUBLET_SCORE_THRESHOLD
        print(f"使用手動門檻 doublet_score > {DOUBLET_SCORE_THRESHOLD}。")

    n_doublets = int(is_doublet.sum())
    print(f"標記 {n_doublets} 個doublet ({n_doublets / adata.shape[0]:.1%})。")

    adata = adata[~is_doublet].copy()
    print(f"After doublet removal: {adata.shape[0]} cells x {adata.shape[1]} genes")

    output_path = os.path.join(atlas_dir, "qc_filtered_atlas.h5ad")
    adata.write_h5ad(output_path)
    print(f"Saved to: {output_path}")


if __name__ == "__main__":
    main()