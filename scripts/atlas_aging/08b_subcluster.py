"""
08b_subcluster.py

老師覺得similarity matrix(08的輸出)裡,PanSci幾種稀有類型(Rdh16+細胞、
Neural cells、Mesothelial cells等)在我們這邊完全對不到、整欄呈現紫色,
這樣不完整。往回查(見對話記錄裡的分析)發現:這幾種稀有類型有80~99%的
細胞,都被我們歸進籠統的"ASPCs"這個大類別裡,是ASPCs這個類別把它們都
吸收掉了,不是真的不存在。

這支腳本把ASPCs單獨挑出來,重新從原始counts做一次HVG/PCA/分群
(邏輯跟scripts/atlas_annotation/05_subcluster.py完全一樣),看能不能在
ASPCs底下找到對應到這些稀有類型的獨立subcluster。

用法(預設處理ASPCs,也可以指定其他cell_type_final的值):
    python 08b_subcluster.py
    python 08b_subcluster.py ASPCs
    python 08b_subcluster.py Adipocytes 0.5

Input:  results/atlas_aging/wt_aging_final_annotated.h5ad (07的輸出)
Output: results/atlas_aging/subclusters/subcluster_<cell_type>.h5ad
        results/atlas_aging/figures/umap_subcluster_<cell_type>.png
"""

import os
import sys
import scanpy as sc
import harmonypy
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from load_utils import get_project_paths

N_TOP_GENES = 2000
N_PCS = 30
DEFAULT_LEIDEN_RESOLUTION = 1.0
DEFAULT_CELL_TYPE = "ASPCs"
POSTER_DPI = 300
plt.rcParams.update({
    "font.size": 14, "axes.titlesize": 18, "axes.labelsize": 14, "legend.fontsize": 10,
})


def main():
    args = sys.argv[1:]
    if not args:
        leiden_resolution = DEFAULT_LEIDEN_RESOLUTION
        target_cell_type = DEFAULT_CELL_TYPE
    else:
        try:
            leiden_resolution = float(args[-1])
            cell_type_args = args[:-1]
        except ValueError:
            leiden_resolution = DEFAULT_LEIDEN_RESOLUTION
            cell_type_args = args
        target_cell_type = " ".join(cell_type_args) if cell_type_args else DEFAULT_CELL_TYPE

    paths = get_project_paths(__file__)
    out_dir = os.path.join(paths["script_dir"], "..", "..", "results", "atlas_aging")
    subclusters_dir = os.path.join(out_dir, "subclusters")
    figures_dir = os.path.join(out_dir, "figures")
    os.makedirs(subclusters_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    input_path = os.path.join(out_dir, "wt_aging_final_annotated.h5ad")
    print(f"Loading: {input_path}")
    adata = sc.read_h5ad(input_path)

    available_types = adata.obs["cell_type_final"].unique().tolist()
    if target_cell_type not in available_types:
        raise SystemExit(
            f"'{target_cell_type}' 不在cell_type_final裡。可用的選項:\n  "
            + "\n  ".join(sorted(available_types))
        )

    adata = adata[adata.obs["cell_type_final"] == target_cell_type].copy()
    print(f"'{target_cell_type}': {adata.shape[0]} cells x {adata.shape[1]} genes")

    print("\n這個子集裡PanSci原始標註(main_cell_type)的分布(方便對照subcluster結果):")
    print(adata.obs["main_cell_type"].value_counts())

    # 關鍵:adata.X此時是04存檔前就已經normalize+log1p過的數值,要先從
    # layers["counts"]接回原始counts,才能正確地重新走一次完整流程
    adata.X = adata.layers["counts"].copy()

    print("\nNormalizing...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    print("Finding highly variable genes (recomputed within this subset)...")
    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES, batch_key="tissue")

    adata.raw = adata
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    print("Scaling...")
    sc.pp.scale(adata_hvg, max_value=10)

    print("Running PCA (recomputed within this subset)...")
    sc.tl.pca(adata_hvg, n_comps=N_PCS, svd_solver="arpack")

    print("Running Harmony integration (batch=tissue)...")
    harmony_out = harmonypy.run_harmony(adata_hvg.obsm["X_pca"], adata_hvg.obs, "tissue")
    adata_hvg.obsm["X_pca_harmony"] = harmony_out.Z_corr

    print("Computing neighbors and UMAP...")
    sc.pp.neighbors(adata_hvg, n_pcs=N_PCS, use_rep="X_pca_harmony")
    sc.tl.umap(adata_hvg, min_dist=0.1)

    print(f"Clustering (Leiden, resolution={leiden_resolution})...")
    sc.tl.leiden(
        adata_hvg,
        resolution=leiden_resolution,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    n_subclusters = adata_hvg.obs["leiden"].nunique()
    print(f"Found {n_subclusters} subclusters within '{target_cell_type}'.")

    adata.obs["subcluster"] = adata_hvg.obs["leiden"].values
    adata.obsm["X_umap_subcluster"] = adata_hvg.obsm["X_umap"]
    adata.obsm["X_pca_subcluster"] = adata_hvg.obsm["X_pca"]
    adata.obsm["X_pca_harmony_subcluster"] = adata_hvg.obsm["X_pca_harmony"]

    print("\n每個subcluster裡,PanSci原始標註的分布(交叉比對用):")
    print(adata.obs.groupby("subcluster")["main_cell_type"].value_counts().to_string())

    fig = sc.pl.embedding(
        adata, basis="X_umap_subcluster", color="subcluster", show=False, return_fig=True,
        size=8, frameon=False, title=f"{target_cell_type}: subclusters", legend_fontsize=9,
    )
    safe_name = target_cell_type.replace(" ", "_")
    fig_path = os.path.join(figures_dir, f"umap_subcluster_{safe_name}.png")
    fig.savefig(fig_path, dpi=POSTER_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")

    output_path = os.path.join(subclusters_dir, f"subcluster_{safe_name}.h5ad")
    adata.write_h5ad(output_path)
    print(f"\nSaved subclustered atlas to: {output_path}")


if __name__ == "__main__":
    main()
