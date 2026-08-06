# run_seurat_benchmark.R
# 用Seurat(標準的記憶體內稀疏矩陣格式)載入一個mtx子樣本
#
# 對齊 Parks & Greenleaf (2025, bioRxiv 2025.03.27.645853) 論文
# 「RNA Normalization & PCA benchmark」方法。計時的流程只做到 PCA(50個主成分)為止，
# 分群(clustering)跟UMAP刻意不包含在內(論文認為這些是非瓶頸的下游步驟)
#
# 關於 ScaleData 的備註:論文提到除了 Seurat 以外，PCA 過程中都隱含做了置中(mean-centering)，
# 所以 Seurat 這裡保留原本明確呼叫 ScaleData() 的步驟，是刻意對齊論文的 Seurat 設定，不是不一致
#
# 高變異基因改成用 Seurat 自己內建的 FindVariableFeatures() 選，不再讀取外部
# 共用清單——三個工具各自用自己原生的方式選HVG，比較貼近真實使用情境

#相當於python的import
#把訊息隱藏
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(peakRAM)
})

# 單執行緒，對齊論文主要的工具間比較(多執行緒是另外獨立測的，不是這次主要比較的一部分)
# Seurat 核心的 Normalize/Scale/PCA 步驟預設就是單執行緒，除非有另外設定 future::plan()，
# 所以這裡不用額外改設定，這行只是註記說明

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript run_seurat_benchmark.R <subset_dir> <label> <results_csv>")
}
subset_dir  <- args[1]
label       <- args[2]
results_csv <- args[3]

N_PCS <- 50  # 主成分數量，對齊論文設定(50個成分)

load_mtx_dataset <- function(subset_dir) {
  mtx_path      <- file.path(subset_dir, "matrix.mtx.gz")
  barcodes_path <- file.path(subset_dir, "barcodes.tsv")
  features_path <- file.path(subset_dir, "features.tsv")

  mat <- Matrix::readMM(gzfile(mtx_path))
  barcodes <- readLines(barcodes_path)
  features <- readLines(features_path)
  rownames(mat) <- make.unique(features)
  colnames(mat) <- barcodes
  mat
}

cat(sprintf("[Seurat] Loading %s ...\n", subset_dir))
mat <- load_mtx_dataset(subset_dir)
cat(sprintf("[Seurat] Loaded %d cells x %d genes\n", ncol(mat), nrow(mat)))

bench <- peakRAM::peakRAM({
  # 步驟1: 只保留在所有細胞裡至少被偵測到1次的基因
  mat_f <- mat[Matrix::rowSums(mat) >= 1, ]

  obj <- CreateSeuratObject(counts = mat_f, min.cells = 0, min.features = 0)

  # 步驟2+4(合併): Seurat 的 NormalizeData() 會在同一次呼叫裡，
  # 先做總表達量縮放成10000、再做log1p(這是預設方法 LogNormalize)——
  # 這是 Seurat 的標準寫法，也是最貼近論文分開的步驟2跟4的實際作法
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                        scale.factor = 10000, verbose = FALSE)

  # 步驟3: Seurat 自己選高變異基因(取代原本讀取共用清單的步驟)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  genes_present <- VariableFeatures(obj)

  # 步驟5: z-score標準化(平均0、標準差1)
  obj <- ScaleData(obj, features = genes_present, verbose = FALSE)

  # 步驟6: PCA降維，取50個主成分
  n_pcs <- min(N_PCS, ncol(obj) - 1, length(genes_present) - 1)
  obj <- RunPCA(obj, features = genes_present, npcs = n_pcs, verbose = FALSE)
})

n_cells <- ncol(mat)
n_genes <- nrow(mat)
elapsed <- bench$Elapsed_Time_sec
peak_mb <- bench$Peak_RAM_Used_MiB

cat(sprintf("[Seurat] Done. elapsed=%.2fs peak_mem=%.1fMB\n", elapsed, peak_mb))

row <- data.frame(
  tool = "Seurat",
  label = label,
  n_cells = n_cells,
  n_genes = n_genes,
  elapsed_sec = round(elapsed, 3),
  peak_memory_MB = round(peak_mb, 1),
  status = "OK"
)

dir.create(dirname(results_csv), showWarnings = FALSE, recursive = TRUE)
write.table(
  row,
  file = results_csv,
  sep = ",",
  row.names = FALSE,
  col.names = !file.exists(results_csv),
  append = file.exists(results_csv)
)

cat(sprintf("[Seurat] Result appended to: %s\n", results_csv))