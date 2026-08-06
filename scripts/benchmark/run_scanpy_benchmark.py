import os

#限制執行緒數量的環境變數，一定要在 import numpy/scanpy 之前設定才會生效
#這樣三個工具(Scanpy/Seurat/BPCells)才是在「單執行緒」的公平基準下比較，
#不會有些工具偷偷用多核心運算而佔便宜
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMBA_NUM_THREADS"] = "1"

import sys #存取命令列參數與系統功能
import gzip #讀取壓縮的mtx檔案
import time #計算執行時間
import csv #將結果寫入csv檔

import numpy as np
import scanpy as sc #單細胞資料分析
import anndata as ad #單細胞的資料格式
from scipy.io import mmread #從 scipy 匯入 mmread，用於讀取Matrix Market格式檔案

N_PCS = 50  #主成分數量，對齊參考論文的設定(50 個成分)
N_HVG = 2000  #高變異基因數量。改成每個工具各自用自己的原生方法選，
              #不再讀取一份三個工具共用的清單(見檔案下方run_pipeline的說明)


def get_peak_memory_mb() -> float:
    """取得該程序的峰值記憶體，單位為 MB
    概念上對應論文用「GNU time -v 的 Max RSS」做的量測方式(程序歷史最高記憶體用量)"""
    #resource 模組只存在於 macOS/Linux，Windows 上沒有這個模組，
    #所以 Windows 改用 psutil 讀取 peak_wset（Windows 對等的峰值記憶體欄位）
    if sys.platform == "win32":
        import psutil
        return psutil.Process(os.getpid()).memory_info().peak_wset / (1024 ** 2)
    import resource
    peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss #獲取系統報告最大駐留集大小
    if sys.platform == "darwin": #若作業系統是 macOS
        return peak_kb / (1024 ** 2) # macOS 的 ru_maxrss 單位是 byte ，除以1024的二次方轉MB
    return peak_kb / 1024 #Linux的ru_maxrss單位KB，除以1024轉MB


#負責把03腳本輸出的一組子樣本檔案，讀回成Scanpy慣用的AnnData物件
#將三個分開存的檔案重組成Scanpy能直接用的單一物件且處理好行列方向
def load_mtx_dataset(subset_dir: str) -> ad.AnnData:
    mtx_path = os.path.join(subset_dir, "matrix.mtx.gz")
    barcodes_path = os.path.join(subset_dir, "barcodes.tsv")
    features_path = os.path.join(subset_dir, "features.tsv")

    #mmread(f):讀取Matrix Market格式，還原成矩陣(原:基因*細胞)
    #.tocsr:先轉成CSR(壓縮稀疏列)格式
    #第二次:轉置後矩陣內部會變成CSC，再轉一次確保運算效率
    with gzip.open(mtx_path, "rt") as f:
        X = mmread(f).tocsr().T.tocsr()  # mtx is genes x cells -> cells x genes

    #兩者各自讀成一個list，strip()換掉換行符
    with open(barcodes_path) as f:
        barcodes = [line.strip() for line in f]
    with open(features_path) as f:
        #Seurat 會強制把基因名稱裡的 "_" 跟 "|" 換成 "-"（無法關閉這個行為），
        #這裡統一套用同一個規則，確保跟 variable_genes.txt(由 Seurat 產生)對得上名稱
        features = [line.strip().replace("_", "-").replace("|", "-") for line in f]

    #用剛讀到的矩陣建立一個AnnData物件，內部結構是X存數值矩陣，obs存細胞資訊，var基因資訊(欄)
    adata = ad.AnnData(X=X)
    adata.obs_names = barcodes #把細胞ID指派給列的名稱
    adata.var_names = features #把基因ID指派給欄的名稱
    return adata


#執行到 PCA 為止(不含分群)。高變異基因改成Scanpy自己內建的方法選，
#不再讀取外部共用清單——三個工具各自用自己原生的方式選HVG，比較貼近
#真實使用情境(使用者實際用某個工具時，本來就是用它自己內建的HVG演算法，
#不會手動注入別的工具算出來的清單)。代價是三個工具選出來的基因集合
#不會100%一樣，但矩陣維度(都是約2000個基因)一致，對效能量測的影響很小
def run_pipeline(adata: ad.AnnData) -> ad.AnnData:
    # 步驟1: 只保留在所有細胞裡至少被偵測到 1 次的基因
    sc.pp.filter_genes(adata, min_counts=1)

    # 步驟2: 正規化，把每個細胞的總表達量縮放成 10000
    sc.pp.normalize_total(adata, target_sum=1e4)

    # 步驟3: log(x+1) 轉換
    sc.pp.log1p(adata)

    # 步驟4: Scanpy 自己選高變異基因(取代原本讀取共用清單的步驟)
    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG)
    adata = adata[:, adata.var["highly_variable"]].copy()

    # 步驟5: z-score 標準化(平均0、標準差1)。這裡不做極端值裁切(max_value=None)，
    # 對齊論文單純的 z-score 步驟，不加上限
    sc.pp.scale(adata, max_value=None)

    # 步驟6: PCA降維，取 50 個主成分
    n_comps = min(N_PCS, min(adata.shape) - 1)
    sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack")

    return adata


def main():
    #sys.argv是執行程式時，終端機傳進來的命令列參數清單
    if len(sys.argv) != 4:
        print("Usage: python run_scanpy_benchmark.py <subset_dir> <label> <results_csv>")
        sys.exit(1)

    subset_dir, label, results_csv = sys.argv[1], sys.argv[2], sys.argv[3]

    print(f"[Scanpy] Loading {subset_dir} ...")
    adata = load_mtx_dataset(subset_dir)
    n_cells_full, n_genes_full = adata.shape
    print(f"[Scanpy] Loaded {n_cells_full} cells x {n_genes_full} genes")

    t0 = time.perf_counter() #記下開始計時的時間點
    print("[Scanpy] Running pipeline (filter -> normalize -> log1p -> HVG "
          "-> scale -> PCA)...")
    adata = run_pipeline(adata)
    elapsed = time.perf_counter() - t0

    peak_mb = get_peak_memory_mb()
    n_cells, n_genes = adata.shape

    print(f"[Scanpy] Done. elapsed={elapsed:.2f}s peak_mem={peak_mb:.1f}MB "
          f"(final: {n_cells} cells x {n_genes} genes used in PCA)")

    row = {
        "tool": "Scanpy",
        "label": label,
        "n_cells": n_cells_full,
        "n_genes": n_genes_full,
        "elapsed_sec": round(elapsed, 3),
        "peak_memory_MB": round(peak_mb, 1),
        "status": "OK",
    }

    file_exists = os.path.exists(results_csv)
    os.makedirs(os.path.dirname(results_csv), exist_ok=True)
    with open(results_csv, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)

    print(f"[Scanpy] Result appended to: {results_csv}")


if __name__ == "__main__":
    main()