"""Plot previously verified cultural-transmission calculations.

Run from any working directory:
    python cultural_transmission_research/code/plot_transmission.py
    python cultural_transmission_research/code/plot_transmission.py --paper

Inputs are read only. The N=2 variance curve uses the analytical expression
implemented in transmission.py; the larger-population points are the saved
finite-state linear-system solutions. Error bars are Monte Carlo standard
errors of estimated stationary summaries, not uncertainty for single counts.
"""
from __future__ import annotations

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import json
import os
import argparse
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="culture_plot_"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np

from transmission import n2_population_variance


ROOT = Path(__file__).resolve().parents[1]
PLUM = "#78395b"
OCHRE = "#bd7631"
TEAL = "#377f78"
INK = "#263333"
MUTED = "#586261"
GRID = "#dfe5e3"


def main(paper: bool = False, out: Path | None = None, results: Path | None = None) -> None:
    results = ROOT / "results" if results is None else results
    large = json.loads((results / "large_N_moments.json").read_text())
    verification = json.loads((results / "verification.json").read_text())
    n50 = sorted(
        [d for d in large if d["N"] == 50 and d["n"] == 50],
        key=lambda d: d["b"],
    )
    forward = sorted(verification["forward"], key=lambda d: d["b"])

    # Guard against silently plotting a changed parameter block or incomplete run.
    if [d["b"] for d in n50] != [0, .1, .5, .9, .99, 1]:
        raise ValueError("Expected the six prespecified N=50 numerical points.")
    if [d["b"] for d in forward] != [.1, .25, .5, 1]:
        raise ValueError("Expected four prespecified forward-simulation runs.")
    for d in n50:
        if (d["beta"], d["u"], d["mu"]) != (1.3, 1.0, .1):
            raise ValueError("Unexpected parameters in the N=50 calculation.")
    for d in forward:
        if (d["N"], d["beta"], d["u"], d["mu"], d["cycles"]) != (
            2, 2, 1, .5, 100000
        ):
            raise ValueError("Unexpected forward-simulation parameters.")
    if np.ptp([d["mean"] for d in n50]) > 1e-9:
        raise ValueError("The expected repertoire should be constant in b.")

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11.3 if paper else 11,
        "axes.titlesize": 11.5,
        "axes.labelsize": 11.2 if paper else 10.5,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelsize": 10.8 if paper else 9.5,
        "ytick.labelsize": 10.8 if paper else 9.5,
        "axes.edgecolor": "#85928e",
        "axes.linewidth": .8,
        "mathtext.fontset": "dejavusans",
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
    })

    if paper:
        fig, axes = plt.subplots(3, 1, figsize=(6.8, 8.25))
        fig.subplots_adjust(left=.17, right=.97, bottom=.17, top=.865, hspace=.60)
        fig.suptitle(
            "Shared learning encounters change variability",
            x=.17, y=.98, ha="left", fontsize=13.0, fontweight="semibold",
        )
        fig.text(.17, .94, "Stationary neutral model · Whole-population counts",
                 fontsize=10.8, color=MUTED)
        fig.text(.17, .913, r"Effective transmission rate $\beta$ fixed in each comparison",
                 fontsize=10.8, color=MUTED)
    else:
        fig, axes = plt.subplots(1, 3, figsize=(12.1, 4.5))
        fig.subplots_adjust(left=.066, right=.988, bottom=.205, top=.775, wspace=.34)
        fig.suptitle(
            "Shared learning encounters change variability at fixed mean",
            x=.066, y=.965, ha="left", fontsize=15, fontweight="semibold",
        )
        fig.text(
            .066, .898,
            "Stationary neutral model  ·  Whole-population trait counts  ·  "
            r"Effective transmission rate $\beta$ held fixed within each comparison",
            fontsize=10.0, color=MUTED,
        )

    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color=GRID, linewidth=.7, zorder=0)
        ax.set_axisbelow(True)
        ax.set_xlim(-.045, 1.055)
        ax.set_xticks([0, .25, .5, .75, 1])
        ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])
        ax.tick_params(length=3.5, width=.7)
        ax.set_xlabel(r"Per-trait copying probability, $b$", labelpad=8)
    if paper:
        for ax in axes[:2]:
            ax.set_xlabel("")

    # A: closed-form population variance and independent typed forward simulation.
    ax = axes[0]
    ax.set_title("A   Two-person variance" + ("  |  " if paper else "\n") +
                 r"$N=2,\ \beta=2,\ u=1,\ \mu=0.5$", loc="left", pad=11 if paper else 13)
    dense_b = np.linspace(0, 1, 501)
    exact_v = np.array([n2_population_variance(2, 1, .5, b)[1] for b in dense_b])
    ax.plot(dense_b, exact_v, color=PLUM, lw=2.25, zorder=2)
    ax.scatter([0], [exact_v[0]], s=43, facecolors="white", edgecolors=TEAL,
               linewidths=1.5, zorder=5)
    sim_b = np.array([d["b"] for d in forward])
    sim_v = np.array([d["variance"] for d in forward])
    sim_v_error = 1.96 * np.array([d["variance_se"] for d in forward])
    ax.errorbar(sim_b, sim_v, yerr=sim_v_error, fmt="o", ms=5.0,
                color=OCHRE, ecolor=OCHRE, elinewidth=1.25, capsize=3,
                markeredgecolor="white", markeredgewidth=.65, zorder=4)
    ax.set_ylim(2.56, 3.39)
    ax.yaxis.set_major_locator(MultipleLocator(.2))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.set_ylabel((r"Variance, $\mathrm{Var}(C_b)$" if paper else
                   r"Trait-count variance, $\mathrm{Var}(C_b)$"), labelpad=7)
    ax.text(.04, .96, r"$E[C_b]=1.5$ for every $b$", transform=ax.transAxes,
            va="top", fontsize=10.8 if paper else 9.8)
    ax.legend(handles=[
        Line2D([], [], color=PLUM, lw=2.25, label="Exact formula"),
        Line2D([], [], color=OCHRE, marker="o", lw=0, markersize=5,
               label="Forward simulation"),
    ], loc="upper left", bbox_to_anchor=(.01, .82 if paper else .88), frameon=False,
       fontsize=10.8 if paper else 9.2, handlelength=2.05, labelspacing=.55,
       ncols=2 if paper else 1, columnspacing=1.0)

    # B: saved numerical solutions of the exact finite-state moment equations.
    ax = axes[1]
    ax.set_title("B   Population of 50" + ("  |  " if paper else "\n") +
                 r"$N=50,\ \beta=1.3,\ u=1,\ \mu=0.1$", loc="left", pad=11 if paper else 13)
    large_b = np.array([d["b"] for d in n50])
    large_v = np.array([d["variance"] for d in n50])
    ax.plot(large_b, large_v, color=PLUM, lw=1.6, ls=(0, (3, 2)), zorder=2)
    ax.scatter(large_b[1:], large_v[1:], s=33, facecolors=PLUM,
               edgecolors="white", linewidths=.65, zorder=4)
    ax.scatter([0], [large_v[0]], s=43, facecolors="white", edgecolors=TEAL,
               linewidths=1.5, zorder=5)
    ax.set_ylim(0, 1025)
    ax.yaxis.set_major_locator(MultipleLocator(200))
    ax.set_ylabel((r"Variance, $\mathrm{Var}(C_b)$" if paper else
                   r"Trait-count variance, $\mathrm{Var}(C_b)$"), labelpad=7)
    ax.text(.04, .96, r"$E[C_b]=39.1408$ for every $b$", transform=ax.transAxes,
            va="top", fontsize=10.8 if paper else 9.8)
    ratio = large_v[-1] / large_v[0]
    ax.text(.04, .75 if paper else .82, f"{ratio:.2f}× variance from" +
            (" " if paper else "\n") + r"$b\to0$ to $b=1$", transform=ax.transAxes,
            va="top", fontsize=11.0 if paper else 10.1, color=PLUM, linespacing=1.4)
    ax.text(.04, .49 if paper else .59, "Dots: moment-equation solutions\n"
            "Dashes: guide between points", transform=ax.transAxes,
            va="top", fontsize=10.8 if paper else 8.8, color=MUTED, linespacing=1.35 if paper else 1.45)

    # C: independent simulation of stationary absence; one exact endpoint check.
    # For b=1 the two-line ancestry Laplace equations at mu=.5 give 4/11.
    ax = axes[2]
    ax.set_title("C   Probability of no traits" + ("  |  " if paper else "\n") +
                 r"$N=2,\ \beta=2,\ u=1,\ \mu=0.5$", loc="left", pad=11 if paper else 13)
    zero = np.array([d["p_zero"] for d in forward])
    zero_error = 1.96 * np.array([d["p_zero_se"] for d in forward])
    ax.errorbar(sim_b, zero, yerr=zero_error, fmt="o", ms=5.0,
                color=OCHRE, ecolor=OCHRE, elinewidth=1.25, capsize=3,
                markeredgecolor="white", markeredgewidth=.65, zorder=4)
    exact_zero = 4 / 11
    ax.hlines(exact_zero, .88, 1.04, color=TEAL, lw=1.4, zorder=3)
    ax.scatter([1], [exact_zero], marker="D", s=63, facecolors="none",
               edgecolors=TEAL, linewidths=1.35, zorder=5)
    ax.annotate(r"Exact $b=1$: $4/11$", xy=(.94, exact_zero),
                xytext=(.49 if paper else .36, .3682), fontsize=10.8 if paper else 9.4, color=TEAL,
                arrowprops=dict(arrowstyle="-", color=TEAL, lw=.8),
                ha="left", va="center")
    ax.set_ylim(.325, .3735)
    ax.yaxis.set_major_locator(MultipleLocator(.01))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.set_ylabel((r"Absence, $P(C_b=0)$" if paper else
                   r"Stationary probability, $P(C_b=0)$"), labelpad=7)

    if paper:
        fig.text(.17, .018,
                 r"Error bars: estimate $\pm\,1.96$ Monte Carlo SE."
                 "\n100,000 regenerative cycles per setting."
                 "\n" + r"Open teal circles: $b=0$ limiting model.",
                 fontsize=10.8, color=MUTED, linespacing=1.45)
    else:
        fig.text(
            .066, .055,
            r"Error bars: estimate $\pm\,1.96$ Monte Carlo SE from 100,000 "
            "regenerative cycles per setting.  "
            r"Open teal circles: $b=0$ limiting model.",
            fontsize=9.4, color=MUTED,
        )

    default = ROOT / "figures"
    out = default if out is None else out
    out.mkdir(parents=True, exist_ok=True)
    path = out / ("transmission_variability_paper.png" if paper else "transmission_variability.png")
    if path.exists() and path.resolve().parent != default.resolve():
        raise SystemExit(f"Refusing to overwrite existing file outside figures/: {path}")
    fig.savefig(path, dpi=360, bbox_inches="tight", pad_inches=.14,
                metadata={"Title": "Stationary cultural counts under partial transmission",
                          "Description": "Exact model moments and independent forward Monte Carlo; not empirical data."})
    plt.close(fig)
    print(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper", action="store_true", help="Use a stacked 6.8-by-8.25-inch layout.")
    parser.add_argument("--out", help="figure directory (default: figures/)")
    parser.add_argument("--results", help="directory holding the input JSON (default: results/)")
    args = parser.parse_args()
    main(paper=args.paper,
         out=None if args.out is None else Path(args.out).expanduser().resolve(),
         results=None if args.results is None else Path(args.results).expanduser().resolve())
