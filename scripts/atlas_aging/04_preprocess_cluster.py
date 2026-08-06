"""
04_preprocess_cluster.py

前處理+分群,邏輯跟scripts/atlas_annotation/01_filter_wt_cluster.py一樣
(normalize->log1p->HVG->scale->PCA->Harmony->neighbors->Leiden->UMAP),
差別只在這裡不用再篩WT——01已經在原始資料階段就先篩過genotype=WT了。

Harmony批次校正的batch_key沿用"tissue"(跟atlas/atlas_annotation兩條既有
管線一致)。age_group刻意不當作Harmony的batch——Harmony是用來拉近「技術性」
差異(不同組織處理批次造成的差異),而年齡(3個月 vs 16個月)在這裡是我們
真正想觀察的生物學變數,不該被當成雜訊拉近,否則會把老化造成的細胞組成
差異也一起校正掉,整個分析就沒意義了。

Input:  results/atlas_aging/qc_filtered_atlas.h5ad
Output: results/atlas_aging/wt_aging_clustered.h5ad
        results/atlas_aging/figures/umap_leiden.png
        results/atlas_aging/figures/umap_tissue.png
        results/atlas_aging/figures/umap_age_group.png
"""

import os
import scanpy as sc
import harmonypy
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from load_utils import get_project_paths

N_TOP_GENES = 2000
N_PCS = 30
LEIDEN_RESOLUTION = 1.0
POSTER_DPI = 300
plt.rcParams.update({
    "font.size": 14, "axes.titlesize": 18, "axes.labelsize": 14, "legend.fontsize": 10,
})


def main():
    paths = get_project_paths(__file__)
    out_dir = os.path.join(paths["script_dir"], "..", "..", "results", "atlas_aging")
    figures_dir = os.path.join(out_dir, "figures")
    os.makedirs(figures_dir, exist_ok=True)

    input_path = os.path.join(out_dir, "qc_filtered_atlas.h5ad")
    print(f"Loading: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"Input: {adata.shape[0]} cells x {adata.shape[1]} genes")

    print("\nAge group distribution:")
    print(adata.obs["age_group"].value_counts())
    print("\nTissue x age_group breakdown:")
    print(adata.obs.groupby("tissue")["age_group"].value_counts())

    adata.layers["counts"] = adata.X.copy()

    print("\nNormalizing...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    print("Finding highly variable genes...")
    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES, batch_key="tissue")

    adata.raw = adata
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    print("Scaling...")
    sc.pp.scale(adata_hvg, max_value=10)

    print("Running PCA...")
    sc.tl.pca(adata_hvg, n_comps=N_PCS, svd_solver="arpack")

    print("Running Harmony integration (batch=tissue)...")
    harmony_out = harmonypy.run_harmony(adata_hvg.obsm["X_pca"], adata_hvg.obs, "tissue")
    adata_hvg.obsm["X_pca_harmony"] = harmony_out.Z_corr

    print("Computing neighbors and UMAP...")
    sc.pp.neighbors(adata_hvg, n_pcs=N_PCS, use_rep="X_pca_harmony")
    sc.tl.umap(adata_hvg, min_dist=0.1)

    print(f"Clustering (Leiden, resolution={LEIDEN_RESOLUTION})...")
    sc.tl.leiden(
        adata_hvg,
        resolution=LEIDEN_RESOLUTION,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    print(f"Found {adata_hvg.obs['leiden'].nunique()} clusters.")

    adata.obs["leiden"] = adata_hvg.obs["leiden"].values
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    adata.obsm["X_pca_harmony"] = adata_hvg.obsm["X_pca_harmony"]

    plt_kwargs = dict(show=False, return_fig=True, size=8, frameon=False, legend_fontsize=9)

    fig = sc.pl.umap(adata, color="leiden", title="3mo+16mo WT atlas: Leiden clusters", **plt_kwargs)
    fig.savefig(os.path.join(figures_dir, "umap_leiden.png"), dpi=POSTER_DPI, bbox_inches="tight")
    plt.close(fig)

    fig = sc.pl.umap(adata, color="tissue", title="3mo+16mo WT atlas: tissue", **plt_kwargs)
    fig.savefig(os.path.join(figures_dir, "umap_tissue.png"), dpi=POSTER_DPI, bbox_inches="tight")
    plt.close(fig)

    fig = sc.pl.umap(adata, color="age_group", title="3mo+16mo WT atlas: age group", **plt_kwargs)
    fig.savefig(os.path.join(figures_dir, "umap_age_group.png"), dpi=POSTER_DPI, bbox_inches="tight")
    plt.close(fig)

    output_path = os.path.join(out_dir, "wt_aging_clustered.h5ad")
    adata.write_h5ad(output_path)
    print(f"\nSaved clustered atlas to: {output_path}")


if __name__ == "__main__":
    main()