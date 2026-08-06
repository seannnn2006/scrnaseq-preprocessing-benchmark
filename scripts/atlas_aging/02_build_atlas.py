"""
02_build_atlas.py

讀01篩選、匯出的3個組織(WT-only,3個月+16個月),合併成單一個AnnData物件。
邏輯跟scripts/atlas/04_build_atlas.py完全一樣,只是輸入換成01的輸出資料夾。

Input:  data/mouse-pansci-filtered-3-16month-wt/<tissue>/
Output: results/atlas_aging/combined_atlas.h5ad
"""

import os
import gzip
import anndata as ad
import pandas as pd
from scipy.io import mmread

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from load_utils import get_project_paths

TISSUES = ["BAT", "iWAT", "gWAT"]


def load_filtered_tissue(filtered_dir: str, tissue: str) -> ad.AnnData:
    tissue_dir = os.path.join(filtered_dir, tissue)
    mtx_path = os.path.join(tissue_dir, "matrix.mtx.gz")
    barcodes_path = os.path.join(tissue_dir, "barcodes.tsv")
    features_path = os.path.join(tissue_dir, "features.tsv")
    meta_path = os.path.join(tissue_dir, "meta.tsv")

    for p in (mtx_path, barcodes_path, features_path, meta_path):
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"[{tissue}] Missing {p}. Run 01_filter_export_mtx.py first."
            )

    print(f"[{tissue}] Loading filtered mtx...")
    with gzip.open(mtx_path, "rt") as f:
        X = mmread(f).tocsr().T.tocsr()  # mtx is genes x cells -> cells x genes

    with open(barcodes_path) as f:
        barcodes = [line.strip() for line in f]
    with open(features_path) as f:
        features = [line.strip() for line in f]

    adata = ad.AnnData(X=X)
    adata.obs_names = barcodes
    adata.var_names = features

    meta = pd.read_csv(meta_path, sep="\t", index_col=0)
    adata.obs = meta.loc[adata.obs_names]
    adata.obs["tissue"] = tissue

    print(f"[{tissue}] {adata.shape[0]} cells x {adata.shape[1]} genes")
    return adata


def main():
    paths = get_project_paths(__file__)
    filtered_dir = os.path.join(paths["script_dir"], "..", "..", "data", "mouse-pansci-filtered-3-16month-wt")
    results_dir = os.path.join(paths["script_dir"], "..", "..", "results", "atlas_aging")
    os.makedirs(results_dir, exist_ok=True)

    adata_list = [load_filtered_tissue(filtered_dir, t) for t in TISSUES]

    print("\nMerging tissues...")
    adata_combined = ad.concat(adata_list, label="batch", keys=TISSUES, join="inner")
    print(f"Combined atlas: {adata_combined.shape[0]} cells x {adata_combined.shape[1]} genes")

    if "genotype" in adata_combined.obs.columns:
        print("\nGenotype distribution (should be 100% WT):")
        print(adata_combined.obs["genotype"].value_counts())

    if "age_group" in adata_combined.obs.columns:
        print("\nAge group distribution by tissue:")
        print(adata_combined.obs.groupby("tissue")["age_group"].value_counts())

    output_path = os.path.join(results_dir, "combined_atlas.h5ad")
    adata_combined.write_h5ad(output_path)
    print(f"\nSaved combined atlas to: {output_path}")


if __name__ == "__main__":
    main()