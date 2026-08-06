"""
average_results.py

讀取 run_all.sh 產生的多輪結果(benchmark_results_round1.csv、round2.csv、
round3.csv……)，把同一個 (tool, label) 組合、狀態是 OK 的那幾輪取平均，
輸出成一份 benchmark_results_averaged.csv，讓 summarize_results.py 讀取。

只平均「成功(status=OK)」的輪數，失敗的輪不列入平均。如果某個
(tool, label) 全部輪都失敗，最終這筆的 status 就是 FAILED，elapsed_sec、
peak_memory_MB 留空。輸出多一欄 n_success_rounds/n_total_rounds，記錄
「這個數字是幾輪裡面成功幾輪平均出來的」，方便之後在報告裡註明可信度。

Usage:
    python average_results.py
"""

import glob
import os

import pandas as pd


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.join(script_dir, "..", "..", "results")

    round_csvs = sorted(glob.glob(os.path.join(results_dir, "benchmark_results_round*.csv")))
    if not round_csvs:
        raise FileNotFoundError(
            f"No benchmark_results_round*.csv found in {results_dir}. Run run_all.sh first."
        )

    n_total_rounds = len(round_csvs)
    print(f"Found {n_total_rounds} round file(s): {[os.path.basename(p) for p in round_csvs]}")

    dfs = []
    for path in round_csvs:
        df = pd.read_csv(path)
        dfs.append(df)
    all_rounds = pd.concat(dfs, ignore_index=True)

    rows = []
    for (tool, label), group in all_rounds.groupby(["tool", "label"]):
        n_cells = group["n_cells"].iloc[0]
        n_genes = group["n_genes"].iloc[0]

        ok_rows = group[group["status"] == "OK"]
        n_success = len(ok_rows)

        if n_success == 0:
            rows.append({
                "tool": tool,
                "label": label,
                "n_cells": n_cells,
                "n_genes": n_genes,
                "elapsed_sec": None,
                "peak_memory_MB": None,
                "status": "FAILED",
                "n_success_rounds": 0,
                "n_total_rounds": n_total_rounds,
            })
        else:
            rows.append({
                "tool": tool,
                "label": label,
                "n_cells": n_cells,
                "n_genes": n_genes,
                "elapsed_sec": round(ok_rows["elapsed_sec"].mean(), 3),
                "peak_memory_MB": round(ok_rows["peak_memory_MB"].mean(), 1),
                "status": "OK",
                "n_success_rounds": n_success,
                "n_total_rounds": n_total_rounds,
            })

    averaged = pd.DataFrame(rows).sort_values(["tool", "n_cells"])

    out_path = os.path.join(results_dir, "benchmark_results_averaged.csv")
    averaged.to_csv(out_path, index=False)

    n_failed = (averaged["status"] == "FAILED").sum()
    n_partial = ((averaged["status"] == "OK") & (averaged["n_success_rounds"] < n_total_rounds)).sum()
    print(f"\nSaved averaged results to: {out_path}")
    print(f"  {len(averaged)} (tool, label) combinations total")
    print(f"  {n_failed} fully FAILED (0/{n_total_rounds} rounds succeeded)")
    print(f"  {n_partial} partially successful (some but not all rounds succeeded)")


if __name__ == "__main__":
    main()
