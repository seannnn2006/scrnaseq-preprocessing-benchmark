# scRNA-seq Preprocessing Benchmark: Scanpy vs. Seurat vs. BPCells

Benchmarks Scanpy (Python), Seurat (R), and BPCells (R) for scRNA-seq preprocessing (filter → normalize → log1p → HVG → scale → PCA) on mouse adipose tissue single-cell RNA-seq data (PanSci dataset, BAT/iWAT/gWAT tissues), across dataset sizes from 1,000 to ~2,047,431 cells, measuring runtime and peak memory on a single CPU thread.

## Structure

```
scripts/
  benchmark/              benchmark runner scripts (Scanpy/Seurat/BPCells) + summarizing/averaging
  requirements.txt         Python dependencies
  requirement_r.R          R dependencies

results/
  benchmark_comparison_linear.png / _loglog.png   runtime & memory vs. cell count
  benchmark_results*.csv, benchmark_summary.xlsx   raw + averaged benchmark timings
```

## Data

Raw data and the per-scale benchmark input matrices (`results/benchmark_data/n_*/`) are **not included** — together they run to well over 10 GB, beyond what's practical for git. The PanSci mouse adipose dataset is the original source; the benchmark inputs are reproducible from it by running the scripts in `scripts/benchmark/` in order.

A related aging-atlas annotation analysis on the same dataset lives in a separate repo: `adipose-aging-atlas`.