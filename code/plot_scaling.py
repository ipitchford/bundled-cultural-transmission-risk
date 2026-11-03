"""Plot V_N(b)/V_N(1) against b from results/scaling.json (reads only).

Usage: python3 -I -B code/plot_scaling.py [--results DIR] [--out DIR]
Writes figures/variance_scaling.png: one panel per beta, one curve per N.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse
import json
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="culture_plot_"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from certify import DEFAULT_FIGURES, DEFAULT_RESULTS, output_dir, check_writable

INK, MUTED, GRID = "#263333", "#586261", "#dfe5e3"
# sequential, colour-blind-safe (viridis samples, light to dark with N)
COLS = ["#bddf26", "#6ece58", "#25ab82", "#1f988b", "#2c728e", "#433e85", "#440154"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", help="directory holding scaling.json (default: results/)")
    parser.add_argument("--out", help="figure directory (default: figures/)")
    args = parser.parse_args()
    res = Path(args.results).resolve() if args.results else DEFAULT_RESULTS
    data = json.loads((res / "scaling.json").read_text())
    out = output_dir(args.out, DEFAULT_FIGURES)
    target = check_writable(out / "variance_scaling.png", DEFAULT_FIGURES)
    runs = [r for r in data["results"] if not r.get("failed")]
    betas = sorted({r["beta"] for r in runs})
    Ns = sorted({r["N"] for r in runs})
    col = {N: COLS[min(i, len(COLS) - 1)] for i, N in enumerate(Ns)}
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.edgecolor": "#85928e",
                         "savefig.facecolor": "white", "figure.facecolor": "white"})
    fig, axes = plt.subplots(1, len(betas), figsize=(3.1 * len(betas), 3.6), sharey=True, squeeze=False)
    axes = axes[0]
    for ax, beta in zip(axes, betas):
        for r in sorted([r for r in runs if r["beta"] == beta], key=lambda r: r["N"]):
            bs = [x["b"] for x in r["by_b"]]
            ys = [x["ratio_to_b1"] for x in r["by_b"]]
            ax.plot(bs, ys, color=col[r["N"]], lw=1.9, marker="o", ms=3, label=f"N = {r['N']}",
                    markeredgecolor="#33333355", markeredgewidth=0.4)
        ax.set_title(rf"$\beta/u = {beta:g}$", fontsize=12, color=INK)
        ax.set_xlim(0, 1.02)
        ax.set_ylim(0, 1.03)
        ax.grid(color=GRID, lw=0.7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlabel(r"Copying probability $b$")
    axes[0].set_ylabel(r"$V_N(b)\,/\,V_N(1)$")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center right", frameon=False, fontsize=10.5, title="Population",
               bbox_to_anchor=(1.0, 0.5))
    fig.tight_layout(rect=(0, 0, 0.9, 1))
    fig.savefig(target, dpi=300, bbox_inches="tight", pad_inches=0.1,
                metadata={"Title": "Variance relative to whole-repertoire copying",
                          "Description": "Exact pair-system solutions, u=1, mu=0.1; not empirical data."})
    plt.close(fig)
    print(target)


if __name__ == "__main__":
    main()
