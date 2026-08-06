"""
01_filter_export_mtx.py

比較3個月vs16個月(WT-only)的細胞組成差異。

這支對應scripts/02_filter_and_export_mtx.py的邏輯(從data/mouse-pansci/的
原始超大exprMatrix.tsv.gz逐行流式讀取,只留符合條件的欄位),但篩選條件
換成:
    age_group in {03_months, 16_months}  且  genotype == WT

刻意在這一步(讀最原始資料的階段)就先篩掉Rag/Prkdc基因剔除鼠,而不是
像scripts/atlas_annotation/01那樣等全部genotype都跑完QC+doublet偵測後才篩
——因為doublet偵測(scrublet)是整條管線最慢的一步(80萬細胞規模要跑數
小時),如果Rag/Prkdc的細胞根本不會被用到,提早篩掉可以讓scrublet少處理
一大批用不到的細胞,省下大量時間。

Output per tissue (data/mouse-pansci-filtered-3-16month-wt/<tissue>/):
    - matrix.mtx.gz, barcodes.tsv, features.tsv, meta.tsv
"""

import os
import gzip
import time
from operator import itemgetter
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.io import mmwrite

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from load_utils import get_project_paths

TISSUES = ["BAT", "iWAT", "gWAT"]

# mouse-pansci這份資料每個組織的基因數固定是42,332(跟已經跑過的其他管線
# 輸出的h5ad shape一致)。這裡只拿來換算進度百分比給人看,不是精確地重新
# 掃過檔案算出來的——如果哪個組織的基因數剛好不一樣,百分比會有一點誤差,
# 但不影響實際的過濾/匯出邏輯,純粹是顯示用。
TOTAL_GENES_ESTIMATE = 42332
AGE_COLUMN = "age_group"
AGE_FILTER = ["03_months", "16_months"]
GENOTYPE_COLUMN = "genotype"
GENOTYPE_FILTER = ["WT"]


def filter_and_export_tissue(data_dir: str, out_dir: str, tissue: str):
    src_path = os.path.join(data_dir, tissue)
    meta_path = os.path.join(src_path, "meta.tsv")
    expr_path = os.path.join(src_path, "exprMatrix.tsv.gz")

    if not os.path.exists(meta_path) or not os.path.exists(expr_path):
        print(f"[{tissue}] SKIPPED: missing meta.tsv or exprMatrix.tsv.gz in {src_path}")
        return

    print(f"\n=== {tissue} ===")
    print(f"[{tissue}] Reading metadata...")
    meta = pd.read_csv(meta_path, sep="\t", index_col=0)

    for col in (AGE_COLUMN, GENOTYPE_COLUMN):
        if col not in meta.columns:
            raise ValueError(
                f"[{tissue}] '{col}' column not found. Available columns: "
                f"{list(meta.columns)}."
            )

    age_mask = meta[AGE_COLUMN].astype(str).isin({str(v) for v in AGE_FILTER})
    genotype_mask = meta[GENOTYPE_COLUMN].astype(str).isin({str(v) for v in GENOTYPE_FILTER})
    meta_filtered = meta[age_mask & genotype_mask]
    print(f"[{tissue}] {len(meta_filtered)} / {len(meta)} cells match "
          f"{AGE_COLUMN}={AGE_FILTER} AND {GENOTYPE_COLUMN}={GENOTYPE_FILTER}")

    if len(meta_filtered) == 0:
        print(f"[{tissue}] WARNING: no cells matched. Actual {AGE_COLUMN} values: "
              f"{sorted(meta[AGE_COLUMN].dropna().unique().tolist(), key=str)}; "
              f"actual {GENOTYPE_COLUMN} values: "
              f"{sorted(meta[GENOTYPE_COLUMN].dropna().unique().tolist(), key=str)}")
        return

    print(f"[{tissue}] Breakdown by {AGE_COLUMN}:")
    print(meta_filtered[AGE_COLUMN].value_counts())

    selected_barcodes = set(meta_filtered.index.astype(str))

    # -------------------------------------------------------------------
    # 掃過header,決定要保留哪些欄(細胞)
    # -------------------------------------------------------------------
    print(f"[{tissue}] Scanning header...")
    with gzip.open(expr_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")

    matrix_cells = header[1:]
    keep_col_positions = [i for i, c in enumerate(matrix_cells) if c in selected_barcodes]
    kept_cell_ids = [matrix_cells[i] for i in keep_col_positions]
    print(f"[{tissue}] {len(kept_cell_ids)} / {len(matrix_cells)} matrix columns kept.")

    if not kept_cell_ids:
        print(f"[{tissue}] WARNING: none of the filtered barcodes matched matrix "
              f"columns. Check that meta.tsv row names match the matrix header.")
        return

    data, row_idx, col_idx = [], [], []
    genes = []

    # 逐行(逐基因)讀取,避免整份1.4M欄寬的檔案一次載入記憶體
    # (理由跟scripts/02_filter_and_export_mtx.py完全一樣)
    PROGRESS_EVERY = 500
    getter = itemgetter(*keep_col_positions)

    print(f"[{tissue}] Streaming expression matrix line-by-line "
          f"({len(kept_cell_ids)} of {len(matrix_cells)} columns kept)...")

    start_time = time.time()
    with gzip.open(expr_path, "rt") as f:
        f.readline()

        for gene_idx, line in enumerate(f):
            fields = line.rstrip("\n").split("\t")
            genes.append(fields[0])

            row_vals = np.array(getter(fields[1:]), dtype=np.float32)
            nz_cols = np.flatnonzero(row_vals)
            if nz_cols.size:
                data.extend(row_vals[nz_cols])
                row_idx.extend([gene_idx] * nz_cols.size)
                col_idx.extend(nz_cols)

            if (gene_idx + 1) % PROGRESS_EVERY == 0:
                # 百分比是用TOTAL_GENES_ESTIMATE概略換算的,僅供顯示;
                # 用已經花的時間/目前進度,推算大概還要多久跑完(ETA)
                pct = (gene_idx + 1) / TOTAL_GENES_ESTIMATE
                elapsed = time.time() - start_time
                eta_min = (elapsed / (gene_idx + 1)) * (TOTAL_GENES_ESTIMATE - (gene_idx + 1)) / 60
                print(f"[{tissue}]   Progress: {pct:.0%} ({gene_idx + 1}/{TOTAL_GENES_ESTIMATE} genes), "
                      f"elapsed {elapsed / 60:.1f} min, ETA ~{eta_min:.1f} min")

    n_genes = len(genes)
    n_cells = len(kept_cell_ids)
    X = csr_matrix((data, (row_idx, col_idx)), shape=(n_genes, n_cells), dtype=np.float32)
    print(f"[{tissue}] Final sparse matrix: {n_genes} genes x {n_cells} cells "
          f"({X.nnz} nonzero entries)")

    tissue_out = os.path.join(out_dir, tissue)
    os.makedirs(tissue_out, exist_ok=True)

    mtx_path = os.path.join(tissue_out, "matrix.mtx")
    print(f"[{tissue}] Writing {mtx_path}.gz ...")
    mmwrite(mtx_path, X)
    with open(mtx_path, "rb") as f_in, gzip.open(mtx_path + ".gz", "wb") as f_out:
        f_out.writelines(f_in)
    os.remove(mtx_path)

    with open(os.path.join(tissue_out, "barcodes.tsv"), "w") as f:
        f.write("\n".join(kept_cell_ids) + "\n")

    with open(os.path.join(tissue_out, "features.tsv"), "w") as f:
        f.write("\n".join(genes) + "\n")

    meta_filtered.loc[kept_cell_ids].to_csv(os.path.join(tissue_out, "meta.tsv"), sep="\t")

    print(f"[{tissue}] Done. Output written to: {tissue_out}")


def main():
    paths = get_project_paths(__file__)
    # paths["data_dir"]是以script_dir/../data算的,只適合放在scripts/正下方
    # 的腳本(例如scripts/02_filter_and_export_mtx.py)。這支腳本在
    # scripts/atlas_aging/,比scripts/多一層,要多跳一層".."才會到專案根目錄
    # 的data/,否則會指到不存在的scripts/data/mouse-pansci。
    data_dir = os.path.join(paths["script_dir"], "..", "..", "data", "mouse-pansci")
    out_dir = os.path.join(paths["script_dir"], "..", "..", "data", "mouse-pansci-filtered-3-16month-wt")
    os.makedirs(out_dir, exist_ok=True)

    for tissue in TISSUES:
        filter_and_export_tissue(data_dir, out_dir, tissue)

    print(f"\nAll tissues processed. Filtered data in: {os.path.abspath(out_dir)}")


if __name__ == "__main__":
    main()