"""
load_utils.py

Shared helper functions used by multiple scripts in the pipeline
(combine_atlas.py, benchmark_memory.py, qc_filter.py, etc.).

NOTE ON DATA SOURCE (per Zhang et al., Science 2025 / PanSci):
- This mouse-pansci dataset is snRNA-seq (single-nucleus), not scRNA-seq.
  Mitochondrial reads should be near-zero in nuclear preparations, so QC
  mitochondrial thresholds should be stricter than typical whole-cell
  scRNA-seq defaults (see 03_qc_filter.py).
- Wild-type mice were profiled at five discrete ages: 3, 6, 12, 16, and
  23 months. "age" in meta.tsv is expected to take one of these values
  (confirm exact format/dtype with 00_inspect_metadata.py before filtering,
  since it may be stored as int, float, or a string like "3" or "3mo").
- Sex is sex-balanced (male/female) in the original study.
"""

import os
import gzip
import scanpy as sc
import pandas as pd
import anndata as ad


def get_project_paths(this_file: str):
    """
    Compute standard project directories relative to the calling script's
    own location, so paths work regardless of the current working directory.

    Parameters
    ----------
    this_file : str
        Pass __file__ from the calling script.

    Returns
    -------
    dict with keys: script_dir, data_dir, results_dir
    """
    script_dir = os.path.dirname(os.path.abspath(this_file))
    return {
        "script_dir": script_dir,
        "data_dir": os.path.join(script_dir, "..", "data", "mouse-pansci"),
        "results_dir": os.path.join(script_dir, "..", "results"),
    }


def _select_matching_barcodes(meta: pd.DataFrame, filter_criteria: dict, tissue: str):
    """
    Given a metadata DataFrame and a dict of {column_name: [allowed_values]},
    return the subset of meta matching ALL criteria, and print a warning
    (with the available unique values) if a column produces zero matches.

    Values are compared both as-is and as strings, since numeric columns
    like 'age' may be stored as int, float, or string depending on how
    the metadata file was exported.
    """
    mask = pd.Series(True, index=meta.index)

    for col, allowed_values in filter_criteria.items():
        if col not in meta.columns:
            print(f"[{tissue}] WARNING: filter column '{col}' not found in "
                  f"metadata. Available columns: {list(meta.columns)}")
            continue

        allowed_str = {str(v) for v in allowed_values}
        col_mask = meta[col].isin(allowed_values) | meta[col].astype(str).isin(allowed_str)
        n_match = col_mask.sum()
        print(f"[{tissue}] Filter '{col}' in {allowed_values}: "
              f"{n_match} / {len(meta)} cells match.")

        if n_match == 0:
            print(f"[{tissue}] WARNING: no cells matched '{col}' in {allowed_values}. "
                  f"Actual unique values found: "
                  f"{sorted(meta[col].dropna().unique().tolist(), key=str)}")

        mask &= col_mask

    return meta[mask]


def load_tissue_data(
    data_dir: str,
    tissue: str,
    use_sparse: bool = False,
    filter_criteria: dict = None,
):
    """
    Load expression matrix and metadata for a single tissue, optionally
    pre-filtering cells by metadata BEFORE reading the (large) expression
    matrix into memory. This is the key memory-saving step when datasets
    are too large to load in full.

    Parameters
    ----------
    data_dir : str
        Path to the mouse-pansci data folder (contains one sub-folder per tissue).
    tissue : str
        Tissue folder name (e.g. "BAT", "iWAT", "gWAT").
    use_sparse : bool
        If True and matrix.mtx.gz is available, load using the sparse
        10x-style format (sc.read_mtx) instead of the dense CSV reader.
    filter_criteria : dict, optional
        e.g. {"age": [3]} to keep only 3-month-old cells (both sexes),
        or {"age": [3], "sex": ["M", "F"]} to be explicit about keeping
        both sexes. Cells are matched against meta.tsv BEFORE the
        expression matrix is read; only matching columns of
        exprMatrix.tsv.gz are loaded, which drastically reduces memory
        use on large datasets. Run 00_inspect_metadata.py first to
        confirm exact column names/values if unsure.
        If None, all cells are loaded (original behavior).

    Returns
    -------
    AnnData object with expression data, metadata, and a 'tissue' column
    added to .obs.
    """
    path = os.path.join(data_dir, tissue)
    meta_path = os.path.join(path, "meta.tsv")

    if not os.path.exists(meta_path):
        raise FileNotFoundError(f"[{tissue}] Metadata file not found: {meta_path}")

    print(f"[{tissue}] Loading metadata...")
    meta = pd.read_csv(meta_path, sep="\t", index_col=0)

    # -------------------------------------------------------------------
    # Apply metadata-based pre-filter (age / sex / etc.) BEFORE touching
    # the expression matrix. This is what keeps memory usage down.
    # -------------------------------------------------------------------
    if filter_criteria:
        meta = _select_matching_barcodes(meta, filter_criteria, tissue)
        print(f"[{tissue}] {len(meta)} cells remain after metadata filtering.")
        if len(meta) == 0:
            raise ValueError(
                f"[{tissue}] No cells left after filtering with {filter_criteria}. "
                f"Check column names/values with 00_inspect_metadata.py."
            )

    selected_barcodes = set(meta.index.astype(str))

    if use_sparse:
        # NOTE: column-level pre-filtering for sparse .mtx format is more
        # involved (requires filtering the barcodes file and re-indexing
        # the sparse matrix). For now, sparse mode loads all cells and
        # filters afterward — still far more memory-efficient than dense,
        # but if the sparse matrix itself is too large, filter the
        # barcodes.tsv file directly before loading.
        mtx_path = os.path.join(path, "matrix.mtx.gz")
        barcodes_path = os.path.join(path, "barcodes.tsv")
        features_path = os.path.join(path, "features.tsv")

        for p in (mtx_path, barcodes_path, features_path):
            if not os.path.exists(p):
                raise FileNotFoundError(
                    f"[{tissue}] Sparse format requested but missing: {p}"
                )

        print(f"[{tissue}] Loading SPARSE matrix (matrix.mtx.gz)...")
        adata = sc.read_mtx(mtx_path).T  # genes x cells -> cells x genes
        barcodes = pd.read_csv(barcodes_path, header=None)[0].astype(str).values
        features = pd.read_csv(features_path, header=None, sep="\t")[0].values
        adata.obs_names = barcodes
        adata.var_names = features

        if filter_criteria:
            adata = adata[adata.obs_names.isin(selected_barcodes)].copy()

        adata.obs = meta.loc[adata.obs_names]

    else:
        expr_matrix_path = os.path.join(path, "exprMatrix.tsv.gz")
        if not os.path.exists(expr_matrix_path):
            raise FileNotFoundError(
                f"[{tissue}] Expression matrix not found: {expr_matrix_path}"
            )

        if filter_criteria:
            # Read only the header row to find which columns (cell barcodes)
            # are present, then load just the matching columns. This avoids
            # ever materializing the full dense matrix in memory.
            print(f"[{tissue}] Reading header to match barcodes...")
            with gzip.open(expr_matrix_path, "rt") as f:
                header = f.readline().rstrip("\n").split("\t")

            gene_col = header[0]
            matching_cols = [c for c in header[1:] if c in selected_barcodes]
            print(f"[{tissue}] {len(matching_cols)} / {len(header) - 1} "
                  f"columns in exprMatrix.tsv.gz match the filtered barcodes.")

            if not matching_cols:
                raise ValueError(
                    f"[{tissue}] None of the filtered barcodes were found as "
                    f"columns in exprMatrix.tsv.gz. Check that meta.tsv row "
                    f"names match the matrix header exactly."
                )

            print(f"[{tissue}] Loading filtered DENSE matrix "
                  f"({len(matching_cols)} cells only)...")
            usecols = [gene_col] + matching_cols
            df = pd.read_csv(
                expr_matrix_path, sep="\t", usecols=usecols, index_col=0
            )
            # df is genes x cells -> transpose to cells x genes for AnnData
            adata = ad.AnnData(df.T)
            meta = meta.loc[adata.obs_names]
        else:
            print(f"[{tissue}] Loading DENSE matrix (exprMatrix.tsv.gz)...")
            adata = sc.read_csv(expr_matrix_path, delimiter="\t").T

        adata.obs = meta

    adata.obs["tissue"] = tissue

    if "sex" in adata.obs.columns:
        print(f"[{tissue}] Final 'sex' value counts:")
        print(adata.obs["sex"].value_counts(dropna=False))
    if "age" in adata.obs.columns:
        print(f"[{tissue}] Final 'age' value counts:")
        print(adata.obs["age"].value_counts(dropna=False))

    return adata
