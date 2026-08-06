#R package requirements for scripts/benchmark/*.R

options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages(c(
  "Seurat", # standard single-cell workflow
  "SeuratObject",
  "Matrix",
  "peakRAM", # peak memory + elapsed time measurement
  "R.utils" # used by run_bpcells_benchmark.R to decompress .mtx.gz
))


# BPCells install from GitHub
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("bnprks/BPCells/r")