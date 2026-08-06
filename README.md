# Mouse Adipose scRNA-seq: Preprocessing Benchmark + Aging Atlas Annotation

This repo covers two related internship analyses on mouse adipose tissue single-cell RNA-seq data (PanSci dataset, BAT/iWAT/gWAT tissues):

1. **`scripts/benchmark/`** — Benchmarks Scanpy (Python), Seurat (R), and BPCells (R) for scRNA-seq preprocessing (filter → normalize → log1p → HVG → scale → PCA) across dataset sizes from 1,000 to ~2,047,431 cells, measuring runtime and peak memory on a single CPU thread.
2. **`scripts/atlas_aging/`** — Builds a combined 3-month + 16-month WT mouse adipose atlas (518,332 cells), annotates cell types via marker-gene scoring + differential-expression validation + targeted sub-clustering, and validates the final annotation against the PanSci dataset's own `main_cell_type` labels.

## Structure

```
scripts/
  benchmark/            benchmark runner scripts (Scanpy/Seurat/BPCells) + summarizing/averaging
  atlas_aging/           01-08j: filter → QC/doublets → cluster → annotate → sub-cluster → validate
  load_utils.py          shared path-resolution helper used by atlas_aging scripts
  requirements.txt        Python dependencies
  requirement_r.R         R dependencies

results/
  benchmark_comparison_linear.png / _loglog.png   runtime & memory vs. cell count
  benchmark_results*.csv, benchmark_summary.xlsx   raw + averaged benchmark timings
  atlas_aging/
    figures/            UMAPs, similarity-matrix comparisons vs. PanSci labels
    qc_plots/           QC/doublet-score distributions
    *.xlsx              cluster marker scores, DE genes, annotation-comparison tables

reports/
  atlas_aging_meeting_0819/   lab-meeting slides (atlas_aging_report.pptx) + speaker notes
```

## Data

Raw data, intermediate `.h5ad` AnnData objects, and the per-scale benchmark input matrices are **not included** in this repo — the full set is on the order of 150+ GB, well beyond what's practical for git. This includes:

- `data/` — raw and filtered PanSci mouse adipose export (WT/age-filtered subsets)
- `results/atlas_aging/*.h5ad`, `results/atlas_aging/subclusters/` — clustered/annotated AnnData objects at each pipeline stage
- `results/benchmark_data/n_*/` — per-scale subsampled matrices used as benchmark input
- `results/atlas/`, `results/atlas_annotation/` — outputs from a separate (earlier, 3-month-only) analysis track, out of scope for this repo

The PanSci mouse adipose dataset is the original source; every file above is reproducible from it by running the numbered scripts in `scripts/atlas_aging/` and `scripts/benchmark/` in order.