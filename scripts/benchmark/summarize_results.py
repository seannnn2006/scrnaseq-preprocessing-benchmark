"""
summarize_results.py

Reads results/benchmark_results_averaged.csv (produced by
average_results.py, which combines the multiple rounds written by
run_all.sh) and generates comparison plots: peak memory and elapsed
time vs. cell counts, one line per tool (Scanpy, Seurat, BPCells).

"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
 
 
def _first_missing_size(df, tool):
    """Smallest n_cells where `tool` has no result even though it succeeded
    at every smaller size that was attempted -- i.e. the point where it
    likely crashed (e.g. Seurat running out of memory) rather than a size
    that was simply never part of the benchmark plan for that tool."""
    all_sizes = sorted(df["n_cells"].unique())
    tool_sizes = set(df.loc[df["tool"] == tool, "n_cells"])
    if not tool_sizes:
        return None
    max_tool_size = max(tool_sizes)
    missing_beyond = [s for s in all_sizes if s not in tool_sizes and s > max_tool_size]
    return missing_beyond[0] if missing_beyond else None


def make_plot(df, results_dir, log_x: bool, log_y: bool, suffix: str, mem_unit: str = "GB",
              mem_zoom_ylim=None, label_min_n=None,
              crossover_xytext=(0.36, 0.78), foldiff_xytext=(0.74, 0.24)):
    mem_col = "peak_memory_GB" if mem_unit == "GB" else "peak_memory_MB"

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6))

    tool_colors = {}
    for tool, sub in df.groupby("tool"):
        line, = axes[0].plot(sub["n_cells"], sub[mem_col], marker="o", label=tool)
        axes[1].plot(sub["n_cells"], sub["elapsed_sec"], marker="o", label=tool, color=line.get_color())
        tool_colors[tool] = line.get_color()

    # Sized for a poster/slide, not just an on-screen read -- default
    # matplotlib sizes are noticeably small once printed or projected.
    axes[0].set_xlabel("Number of cells", fontsize=13)
    axes[0].set_ylabel(f"Peak memory ({mem_unit})", fontsize=13)
    axes[0].set_title("Peak memory vs. cell count", fontsize=15)
    axes[0].legend(fontsize=11)
    axes[0].tick_params(labelsize=11)

    axes[1].set_xlabel("Number of cells", fontsize=13)
    axes[1].set_ylabel("Elapsed time (sec)", fontsize=13)
    axes[1].set_title("Runtime vs. cell count", fontsize=15)
    axes[1].legend(fontsize=11)
    axes[1].tick_params(labelsize=11)

    # x linear (default): small sizes compress toward the left, large sizes
    # spread out across the middle/right -- more intuitive at a glance than
    # log-x, at the cost of the small-size points being hard to distinguish
    # from each other. y stays log-able independently so the ~4 orders of
    # magnitude spread in memory/time (BPCells vs. Seurat) stays visible even
    # when x is linear.
    if log_x:
        axes[0].set_xscale("log")
        axes[1].set_xscale("log")
    if log_y:
        axes[0].set_yscale("log")
        axes[1].set_yscale("log")
        # Plain numbers ("100", "10,000") instead of matplotlib's default
        # power-of-ten notation ("10^2", "10^4") on a log axis.
        plain_fmt = mticker.FuncFormatter(lambda v, _: f"{v:,.0f}")
        axes[0].yaxis.set_major_formatter(plain_fmt)
        axes[1].yaxis.set_major_formatter(plain_fmt)

    # Mark tools that stopped appearing at some cell count while other
    # tools kept going further -- a dashed vertical line + label at the
    # size where that tool is believed to have crashed (e.g. OOM).
    for tool, color in tool_colors.items():
        fail_at = _first_missing_size(df, tool)
        if fail_at is None:
            continue
        for ax in axes:
            ymax = ax.get_ylim()[1]
            ax.axvline(fail_at, color=color, linestyle="--", alpha=0.5)
            ax.text(fail_at, ymax, f" {tool} failed here ", color=color,
                    rotation=90, va="top", ha="right", fontsize=9.5)

    # Short reference tick (not a full-height line -- that got visually
    # noisy crossing through the data) at each tested size, labeled with
    # the exact cell count -- on a log x-axis most of the 10 sizes don't
    # land on a major gridline (1e3/1e4/1e5/1e6), so the exact value at
    # each point is otherwise unreadable. On a linear x-axis the small
    # sizes are deliberately squeezed together near the left edge, so their
    # labels just overlap into an unreadable smudge -- label_min_n skips
    # labeling (but keeps the tick for) sizes below the given threshold.
    all_sizes = sorted(df["n_cells"].unique())
    label_sizes = [n for n in all_sizes if label_min_n is None or n >= label_min_n]
    for ax, ycol in zip(axes, [mem_col, "elapsed_sec"]):
        ymin, ymax = ax.get_ylim()
        tick_top = ymin * (ymax / ymin) ** 0.045 if log_y else ymin + (ymax - ymin) * 0.045
        for n in all_sizes:
            at_n = df.loc[df["n_cells"] == n, ycol]
            if at_n.empty:
                continue
            ax.plot([n, n], [ymin, tick_top], color="gray", linestyle=":",
                    linewidth=0.9, alpha=0.6, zorder=0)
            if n in label_sizes:
                ax.text(n, ymin, f"{n:,}  ", rotation=90, fontsize=7.5,
                        color="gray", va="bottom", ha="right", alpha=0.9)
        ax.set_ylim(ymin, ymax)

    # Two call-out annotations, both anchored to n=1,000,000 -- the largest
    # size where BPCells, Scanpy, and Seurat all still have data, and (per
    # the actual measured results) the point where BPCells overtakes Scanpy
    # on speed. xytext uses axes-fraction coordinates so the same relative
    # placement works whether the surrounding axes are log or linear.
    callout_n = 1_000_000
    bpc_time = df.loc[(df["tool"] == "BPCells") & (df["n_cells"] == callout_n), "elapsed_sec"]
    scanpy_time = df.loc[(df["tool"] == "Scanpy") & (df["n_cells"] == callout_n), "elapsed_sec"]
    if not bpc_time.empty and not scanpy_time.empty:
        axes[1].annotate(
            "BPCells overtakes Scanpy\nhere (~1M cells)",
            xy=(callout_n, bpc_time.iloc[0]), xycoords="data",
            xytext=crossover_xytext, textcoords="axes fraction",
            fontsize=9.5, color="#0b6e63", ha="center",
            arrowprops=dict(arrowstyle="->", color="#0b6e63", lw=1.3),
        )

    bpc_mem = df.loc[(df["tool"] == "BPCells") & (df["n_cells"] == callout_n), mem_col]
    seurat_mem_mb = df.loc[(df["tool"] == "Seurat") & (df["n_cells"] == callout_n), "peak_memory_MB"]
    bpc_mem_mb = df.loc[(df["tool"] == "BPCells") & (df["n_cells"] == callout_n), "peak_memory_MB"]
    if not bpc_mem.empty and not seurat_mem_mb.empty:
        fold = seurat_mem_mb.iloc[0] / bpc_mem_mb.iloc[0]
        axes[0].annotate(
            f"~{fold:.0f}x less memory\nthan Seurat here",
            xy=(callout_n, bpc_mem.iloc[0]), xycoords="data",
            xytext=foldiff_xytext, textcoords="axes fraction",
            fontsize=9.5, color="#0b6e63", ha="center",
            arrowprops=dict(arrowstyle="->", color="#0b6e63", lw=1.3),
        )

    fig.tight_layout()

    # On a linear y-axis, BPCells' curve (tens to low-hundreds of MB) is
    # dwarfed by Seurat/Scanpy (tens of thousands of MB) and looks like a
    # flat line pinned to zero -- easy to misread as "BPCells uses no
    # memory at all". A zoomed-in inset shows it's actually climbing. It's
    # placed entirely OUTSIDE axes[0], in a blank margin opened up above
    # both subplots (via subplots_adjust), so it can't overlap any line,
    # marker, or label inside the chart -- indicate_inset_zoom() draws the
    # connector lines bridging the gap back to the zoomed-in region.
    if mem_zoom_ylim is not None:
        fig.subplots_adjust(top=0.62)
        pos = axes[0].get_position()
        axins = fig.add_axes([pos.x0 + pos.width * 0.18, 0.82, pos.width * 0.64, 0.14])
        for tool, sub in df.groupby("tool"):
            axins.plot(sub["n_cells"], sub[mem_col], marker="o", markersize=3,
                       linewidth=1.2, color=tool_colors[tool])
        axins.set_xlim(df["n_cells"].min(), df["n_cells"].max())
        axins.set_ylim(*mem_zoom_ylim)
        axins.set_title(f"Zoomed: 0-{mem_zoom_ylim[1]:g} {mem_unit}", fontsize=8.5)
        axins.tick_params(labelsize=7)
        axes[0].indicate_inset_zoom(axins, edgecolor="gray")

    plot_path = os.path.join(results_dir, f"benchmark_comparison_{suffix}.png")
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    print(f"Comparison plot ({suffix}) saved to: {plot_path}")
 
 
def main():
    # This script lives in scripts/benchmark/, so the project's results/
    # folder is two levels up from here (scripts/benchmark/ -> scripts/ -> root).
    script_dir = os.path.dirname(os.path.abspath(__file__))
    results_dir = os.path.join(script_dir, "..", "..", "results")
    csv_path = os.path.join(results_dir, "benchmark_results_averaged.csv")

    if not os.path.exists(csv_path):
        raise FileNotFoundError(
            f"{csv_path} not found. Run benchmark/run_all.sh then "
            f"benchmark/average_results.py first."
        )

    df = pd.read_csv(csv_path)
    print("Loaded averaged benchmark results:")
    print(df)

    failed = df[df["status"] == "FAILED"]
    if len(failed):
        print("\nThe following (tool, label) combinations failed in every round "
              "and are excluded from the plots/tables below:")
        print(failed[["tool", "label", "n_cells", "n_total_rounds"]])

    # Only plot/tabulate combinations that succeeded in at least one round --
    # rows with status == "FAILED" have no elapsed_sec/peak_memory_MB to plot.
    df = df[df["status"] == "OK"].copy()

    # benchmark_results_round*.csv keeps MB (finer-grained, matches what the
    # three benchmark scripts record). GB is only for readability in the
    # plots and summary table, since values now range into the tens of GB.
    df["peak_memory_GB"] = df["peak_memory_MB"] / 1024

    df = df.sort_values(["tool", "n_cells"])

    make_plot(df, results_dir, log_x=True, log_y=True, suffix="loglog", mem_unit="MB")
    make_plot(df, results_dir, log_x=False, log_y=False, suffix="linear", mem_unit="MB",
              mem_zoom_ylim=(0, 1000), label_min_n=300000,
              crossover_xytext=(0.30, 0.60), foldiff_xytext=(0.87, 0.52))

    # Also save a pivoted summary table for quick reading.
    summary_mem = df.pivot(index="n_cells", columns="tool", values="peak_memory_GB").round(2)
    summary_time = df.pivot(index="n_cells", columns="tool", values="elapsed_sec")
    summary_path = os.path.join(results_dir, "benchmark_summary.xlsx")
    with pd.ExcelWriter(summary_path) as writer:
        summary_mem.to_excel(writer, sheet_name="peak_memory_GB")
        summary_time.to_excel(writer, sheet_name="elapsed_sec")
    print(f"Summary tables saved to: {summary_path}")
 
 
if __name__ == "__main__":
    main()
 