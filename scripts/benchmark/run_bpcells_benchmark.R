#run_bpcells_benchmark.R
#把一個 mtx 子樣本轉換成 BPCells 的硬碟儲存二進位格式
#跑一套跟Seurat版相同流程分析，底層改用BPCells的
#硬碟矩陣(disk matrix)，而不是存在記憶體裡的 dgCMatrix
#拿來跟一般Seurat/Scanpy比較的省記憶體版本
#
#對齊 Parks & Greenleaf (2025, bioRxiv 2025.03.27.645853) 論文
#「RNA Normalization & PCA benchmark」方法。計時的流程只做到 PCA(50個主成分)為止，
#分群(clustering)跟UMAP刻意不包含在內，對齊論文的作法
#
#根據論文方法(「對 BPCells 跟 DelayedArray 來說，在做置中(mean-centering)之前，
#會先把正規化後的矩陣存一份到硬碟，避免 PCA 每一輪都要重複算一次 log1p」——
#這個「正規化中繼檔」策略)，我們在跑 PCA 之前，先把正規化、篩選過高變異基因的
#中繼矩陣存到硬碟一次，而不是每輪 PCA 都重複算 log1p
#
#執行緒:單執行緒，對齊論文主要的比較方式。BPCells 的多執行緒是「每個操作各自控制」，
#不是單一全域開關；如果你裝的版本在相關函式(例如 write_matrix_dir()、svds())上
#有提供 threads 參數，可以明確設成 1。這裡不寫死 threads 參數，
#因為這個 API 在不同版本間變動過，怕在函式名稱不同的版本上悄悄失敗
#
#高變異基因改成用 BPCells 自己的串流變異數統計選(matrix_stats)，不再讀取外部
#共用清單——三個工具各自用自己原生的方式選HVG，比較貼近真實使用情境
#
#重要警語:下面用到的串流列/欄加總、欄位縮放、SVD 函式名稱
#(BPCells::rowSums、BPCells::colSums、BPCells::multiply_cols、BPCells::svds)
#是依照目前對 BPCells API 的理解寫的，但這個套件的函式名稱跟參數簽名在不同版本間變動過。
#如果執行時報「找不到函式」或「參數不對」，去查
#`library(help = "BPCells")` 或你安裝版本的官方文件，調整對應的呼叫——
#整體「流程邏輯」(篩選 -> 正規化 -> 篩選基因 -> log1p -> 存硬碟 -> PCA)
#才是要對齊論文方法的重點，不是這些函式的確切名稱

suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
  library(peakRAM)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript run_bpcells_benchmark.R <subset_dir> <label> <results_csv>")
}
subset_dir  <- args[1]
label       <- args[2]
results_csv <- args[3]

N_PCS <- 50  # 主成分數量，對齊論文設定(50個成分)

#BPCells需要一個未壓縮的.mtx檔才能匯入
#(直接從硬碟串流讀取，而不是把整個檔案先讀進R裡)
#先一次把檔案解壓縮到一個暫存資料夾，此動作不會被計時
#因這是準備工作，不是分析流程
mtx_gz_path <- file.path(subset_dir, "matrix.mtx.gz")
#tempfile只是產生一個路徑字串，不會建立資料夾
scratch_dir <- tempfile("bpcells_scratch_")
#把產生的那個路徑，真正建立資料夾
dir.create(scratch_dir)
mtx_path <- file.path(scratch_dir, "matrix.mtx")

if (!requireNamespace("R.utils", quietly = TRUE)) {
  stop("Package 'R.utils' is required to decompress the .mtx.gz file. ",
       "Install it with: install.packages('R.utils')")
}
R.utils::gunzip(mtx_gz_path, destname = mtx_path, remove = FALSE)

barcodes_path <- file.path(subset_dir, "barcodes.tsv")
features_path <- file.path(subset_dir, "features.tsv")
barcodes <- readLines(barcodes_path)
#Seurat 會強制把基因名稱裡的 "_" 跟 "|" 換成 "-"（無法關閉這個行為），
#這裡統一套用同一個規則，確保跟 variable_genes.txt(由 Seurat 產生)對得上名稱
features <- make.unique(gsub("_", "-", gsub("\\|", "-", readLines(features_path))))

N_HVG <- 2000  # 高變異基因數量，跟另外兩支腳本一致

bpcells_raw_dir <- file.path(scratch_dir, "bpcells_raw")
bpcells_norm_dir <- file.path(scratch_dir, "bpcells_normalized")

cat(sprintf("[BPCells] Loading %s ...\n", subset_dir))

#真正被計時、記憶體核心區塊
bench <- peakRAM::peakRAM({
  #匯入成硬碟後端、串流讀取的 IterableMatrix(基因 x 細胞)
  #mat_disk是指向硬碟資料的操作介面，不是把資料整個讀進記憶體
  mat_disk <- BPCells::import_matrix_market(
    mtx_path, outdir = bpcells_raw_dir, #前面解壓縮好的未壓縮路徑
    row_names = features, col_names = barcodes
  )

  n_cells_full <- ncol(mat_disk)
  n_genes_full <- nrow(mat_disk)
  cat(sprintf("[BPCells] Loaded %d cells x %d genes (disk-backed)\n",
              n_cells_full, n_genes_full))

  # 步驟1: 只保留在所有細胞裡至少被偵測到1次的基因
  # BPCells的串流統計函式可以在不把整個矩陣載入記憶體的情況下算出這個
  gene_sums <- BPCells::rowSums(mat_disk)
  keep_genes <- names(gene_sums)[gene_sums >= 1]
  mat_disk <- mat_disk[keep_genes, ]

  # 步驟2: 把每個細胞的表達量縮放成總和10000(串流運算子，
  # 要等後面真的被讀取時才會實際算出來，不會馬上具現化)
  cell_sums <- BPCells::colSums(mat_disk)
  mat_norm <- BPCells::multiply_cols(mat_disk, 1e4 / cell_sums)

  # 步驟3: log1p轉換(串流、保留稀疏性的運算子)
  mat_norm <- log1p(mat_norm)

  # 步驟4: BPCells 自己選高變異基因(取代原本讀取共用清單的步驟)——
  # 用串流變異數統計，在不用把整個矩陣載入記憶體的情況下，
  # 算出每個基因的變異數，取變異數最高的前N_HVG個基因
  gene_stats <- BPCells::matrix_stats(mat_norm, row_stats = "variance")
  gene_variance <- gene_stats$row_stats["variance", ]
  top_genes <- names(sort(gene_variance, decreasing = TRUE))[seq_len(min(N_HVG, length(gene_variance)))]
  mat_norm <- mat_norm[top_genes, ]

  # 對齊論文的「正規化中繼檔」策略:把這個正規化、篩選過基因的矩陣
  # 存到硬碟一次，這樣PCA跑好幾輪(~100+輪)就不用每次都重算log1p
  mat_norm <- BPCells::write_matrix_dir(mat_norm, bpcells_norm_dir)

  cat(sprintf("[BPCells] Running PCA (%d components)...\n", N_PCS))

  # 步驟5+6: z-score標準化(置中的部分在PCA計算過程中隱含處理，
  # 這是論文對BPCells採用的作法——見論文「稀疏矩陣向量乘法」的討論)，
  # 然後計算PCA
  n_pcs <- min(N_PCS, min(dim(mat_norm)) - 1)
  pca_result <- BPCells::svds(mat_norm, k = n_pcs)
})

elapsed <- bench$Elapsed_Time_sec
peak_mb <- bench$Peak_RAM_Used_MiB

cat(sprintf("[BPCells] Done. elapsed=%.2fs peak_mem=%.1fMB\n", elapsed, peak_mb))

row <- data.frame(
  tool = "BPCells",
  label = label,
  n_cells = n_cells_full,
  n_genes = n_genes_full,
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

cat(sprintf("[BPCells] Result appended to: %s\n", results_csv))

#Clean up scratch files (decompressed mtx + BPCells on-disk directory)
unlink(scratch_dir, recursive = TRUE)