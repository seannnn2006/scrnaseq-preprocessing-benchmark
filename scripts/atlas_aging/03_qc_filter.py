"""
03_qc_filter.py

QC過濾+doublet分數計算,邏輯、門檻都跟scripts/atlas/05_qc_filter.py一致
(同一份snRNA-seq資料,QC標準沒有理由變)。這裡資料量已經先在01篩成
WT-only+3個月/16個月,細胞數比全genotype的版本少很多,scrublet(doublet
偵測,整條管線最慢的一步)要處理的細胞數也跟著變少。

跟05一樣,這裡只算doublet_score,不在這裡濾細胞——真正套用門檻的步驟
在03b_filter_doublets.py,那支是秒級操作,可以放心反覆調整門檻重跑。

Input:  results/atlas_aging/combined_atlas.h5ad
Output: results/atlas_aging/qc_filtered_predoublet.h5ad
        results/atlas_aging/qc_plots/doublet_score_distribution.png
"""

import os
import numpy as np
import scanpy as sc
import scanpy.preprocessing._utils as sc_pp_utils
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from load_utils import get_project_paths


def _sample_without_replacement_int64(n_population, n_samples, method="auto", random_state=None):
    """跟scripts/atlas/05_qc_filter.py同一個修正:scrublet在大樣本下用舊版
    numpy RandomState.randint在Windows上會誤判成32-bit整數運算而溢位,
    改用新版numpy Generator.choice避開這個問題。"""
    rng = np.random.default_rng(random_state)
    return rng.choice(n_population, size=n_samples, replace=False).astype(np.int64)


sc_pp_utils.sample_without_replacement = _sample_without_replacement_int64

MIN_GENES_PER_CELL = 100
MIN_CELL_PER_GENE = 3
MAX_PCT_MITO = 15.0


def main():
    paths = get_project_paths(__file__)
    atlas_dir = os.path.join(paths["script_dir"], "..", "..", "results", "atlas_aging")
    plots_dir = os.path.join(atlas_dir, "qc_plots")
    os.makedirs(plots_dir, exist_ok=True)

    input_path = os.path.join(atlas_dir, "combined_atlas.h5ad")
    print(f"Loading: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"Before QC: {adata.shape[0]} cells x {adata.shape[1]} genes")

    gene_symbols = adata.var_names.str.split("|").str[-1]
    adata.var["mt"] = gene_symbols.str.lower().str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    sc.pl.violin(adata, "n_genes_by_counts", ax=axes[0], stripplot=False, show=False)
    sc.pl.violin(adata, "total_counts", ax=axes[1], stripplot=False, show=False)
    sc.pl.violin(adata, "pct_counts_mt", ax=axes[2], stripplot=False, show=False)
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, "qc_violin_before_filtering.png"), dpi=150)
    plt.close(fig)

    n_before = adata.shape[0]
    sc.pp.filter_cells(adata, min_genes=MIN_GENES_PER_CELL)
    sc.pp.filter_genes(adata, min_cells=MIN_CELL_PER_GENE)
    adata = adata[adata.obs["pct_counts_mt"] < MAX_PCT_MITO].copy()

    n_after = adata.shape[0]
    print(f"Filtered out {n_before - n_after} cells ({(n_before - n_after) / n_before:.1%}).")
    print(f"After QC: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # batch_key="ID":跟05一樣,doublet是同一次定序run才會發生的技術性問題,
    # 要照個別老鼠樣本分開偵測
    print("Detecting doublets (Scrublet, per-sample)...")
    sc.pp.scrublet(adata, batch_key="ID", random_state=0)

    n_auto = int(adata.obs["predicted_doublet"].sum())
    print(f"scrublet自動門檻會標記 {n_auto} 個doublet "
          f"({n_auto / adata.shape[0]:.1%})——實際門檻留到03b再決定。")

    for pct in (90, 95, 97, 99, 99.5):
        val = adata.obs["doublet_score"].quantile(pct / 100)
        print(f"  第{pct}百分位數的doublet_score = {val:.4f}")

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(adata.obs["doublet_score"], bins=100)
    ax.set_xlabel("doublet_score")
    ax.set_ylabel("cell count")
    fig.tight_layout()
    fig.savefig(os.path.join(plots_dir, "doublet_score_distribution.png"), dpi=150)
    plt.close(fig)

    if "scrublet" in adata.uns and "batches" in adata.uns["scrublet"]:
        adata.uns["scrublet"]["batches"] = {
            str(k): v for k, v in adata.uns["scrublet"]["batches"].items()
        }

    output_path = os.path.join(atlas_dir, "qc_filtered_predoublet.h5ad")
    adata.write_h5ad(output_path)
    print(f"Saved (doublet-filtering not yet applied) to: {output_path}")


if __name__ == "__main__":
    main()