"""
08j_top50_marker_umap.py

老師的要求:針對08i最終版的cell_type_v2(15種類型),對每一種類型做差異
表現分析(跟03/06一樣用sc.tl.rank_genes_groups),取每種類型前50個
top DE基因,取聯集後只用這批基因重新做一次HVG-free的PCA/Harmony/UMAP
(不再用原本2000個HVG),驗證單靠這些marker基因就能把15種類型重新分開,
是對目前分群/標註結果的獨立驗證。

跟03/06的差異:03/06是用rank_genes_groups的結果去反過來檢查cluster的
標籤對不對(標籤驗證);這支是反過來,用這些marker基因當作重新降維/畫UMAP
的特徵,看embedding本身穩不穩(特徵選擇驗證)。

Input:  results/atlas_aging/wt_aging_v2_annotated.h5ad (08i的輸出)
Output: results/atlas_aging/top50_marker_genes_per_celltype.xlsx (老師要看的表)
        results/atlas_aging/figures/umap_top50marker_celltype_v2.png
"""

import os
import pandas as pd
import scanpy as sc
import harmonypy
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from load_utils import get_project_paths

N_TOP_GENES = 50
N_PCS = 30
POSTER_DPI = 300
plt.rcParams.update({
    "font.size": 14, "axes.titlesize": 18, "axes.labelsize": 14, "legend.fontsize": 9,
})


def main():
    paths = get_project_paths(__file__)
    out_dir = os.path.join(paths["script_dir"], "..", "..", "results", "atlas_aging")
    figures_dir = os.path.join(out_dir, "figures")

    input_path = os.path.join(out_dir, "wt_aging_v2_annotated.h5ad")
    print(f"Loading: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"Input: {adata.shape[0]} cells x {adata.shape[1]} genes, "
          f"{adata.obs['cell_type_v2'].nunique()} cell_type_v2 types")

    print("\nRunning rank_genes_groups (Wilcoxon test, per cell_type_v2)...")
    sc.tl.rank_genes_groups(adata, groupby="cell_type_v2", method="wilcoxon", use_raw=True)
    de_df = sc.get.rank_genes_groups_df(adata, group=None)
    de_df["gene_symbol"] = de_df["names"].str.split("|").str[-1]

    top_de = de_df.groupby("group", sort=False).head(N_TOP_GENES).copy()
    top_de["rank"] = top_de.groupby("group").cumcount() + 1
    top_de = top_de[["group", "rank", "gene_symbol", "names", "scores", "logfoldchanges", "pvals_adj"]]
    top_de = top_de.rename(columns={"group": "cell_type_v2", "names": "var_name"})

    excel_path = os.path.join(out_dir, "top50_marker_genes_per_celltype.xlsx")
    top_de.to_excel(excel_path, index=False)
    print(f"Saved: {excel_path}")

    # 15種類型各50個top基因取聯集(會有重疊,總數比15*50少)
    marker_genes = sorted(top_de["var_name"].unique().tolist())
    print(f"\n{adata.obs['cell_type_v2'].nunique()}種cell_type_v2,"
          f"各取top {N_TOP_GENES}後聯集共{len(marker_genes)}個獨立基因。")

    # -------------------------------------------------------------------
    # 只用這批marker基因重新走一次PCA/Harmony/UMAP,不再用原本的2000個HVG。
    # adata.X此時是normalize+log1p過的完整基因集(跟adata.raw一樣,04存檔前
    # 沒有把X本身subset成HVG,只有暫時的adata_hvg拿去做PCA),可以直接subset。
    # -------------------------------------------------------------------
    adata_markers = adata[:, marker_genes].copy()

    print("\nScaling...")
    sc.pp.scale(adata_markers, max_value=10)

    print("Running PCA (marker-gene-only feature set)...")
    sc.tl.pca(adata_markers, n_comps=min(N_PCS, len(marker_genes) - 1), svd_solver="arpack")

    print("Running Harmony integration (batch=tissue)...")
    harmony_out = harmonypy.run_harmony(adata_markers.obsm["X_pca"], adata_markers.obs, "tissue")
    adata_markers.obsm["X_pca_harmony"] = harmony_out.Z_corr

    print("Computing neighbors and UMAP...")
    sc.pp.neighbors(adata_markers, n_pcs=N_PCS, use_rep="X_pca_harmony")
    sc.tl.umap(adata_markers, min_dist=0.1)

    fig = sc.pl.umap(
        adata_markers, color="cell_type_v2", show=False, return_fig=True,
        size=8, frameon=False, legend_fontsize=9,
        title=f"3-Month + 16-Month WT Adipose Atlas: Final Cell-Type Annotation\n"
              f"(UMAP re-embedded using top {N_TOP_GENES} marker genes/type, "
              f"{len(marker_genes)} genes total)",
    )
    fig_path = os.path.join(figures_dir, "umap_top50marker_celltype_v2.png")
    fig.savefig(fig_path, dpi=POSTER_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")


if __name__ == "__main__":
    main()