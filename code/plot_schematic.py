"""Schematic of bundled versus partial copying (no data; illustration only).

Usage: python3 -I -B code/plot_schematic.py [--out DIR]

Left: bundled copying (b=1), few encounters, each passing every trait the
donor holds. Right: partial copying (b=1/4), four times as many encounters,
each passing each trait with probability 1/4. Two traits (circle / square
markers, two colour-blind-safe colours) start at the same individual. On the
left their transmission paths use the same encounters; on the right they
mostly use different ones. A cross marks a reset (death and empty
replacement). Writes figures/schematic.png.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="culture_plot_"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from certify import DEFAULT_FIGURES, output_dir, check_writable

INK, MUTED, GRID = "#263333", "#586261", "#c9d1ce"
TRAIT_A, TRAIT_B = "#0072B2", "#D55E00"   # Okabe-Ito blue and vermillion
Y = {"1": 2.0, "2": 1.0, "3": 0.0}


def lane(ax, title):
    for name, y in Y.items():
        ax.hlines(y, 0, 10, color=GRID, lw=6, zorder=0, capstyle="round")
        ax.text(-0.35, y, f"Individual {name}", ha="right", va="center", fontsize=11, color=INK)
    ax.set_xlim(-3.2, 10.4)
    ax.set_ylim(-0.75, 3.05)
    ax.axis("off")
    ax.text(-3.2, 2.62, title, ha="left", va="bottom", fontsize=12.5, fontweight="semibold", color=INK)
    ax.annotate("", xy=(10.2, -0.55), xytext=(0, -0.55),
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=1))
    ax.text(5, -0.72, "time", ha="center", va="top", fontsize=10.5, color=MUTED)


def encounter(ax, t, src, dst, a=False, bb=False):
    y0, y1 = Y[src], Y[dst]
    ax.annotate("", xy=(t, y1), xytext=(t, y0),
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.1, shrinkA=4, shrinkB=4,
                                mutation_scale=11), zorder=2)
    off = 0.11
    if a:
        ax.plot([t - off, t - off], [y0, y1], color=TRAIT_A, lw=3.0, zorder=3, solid_capstyle="butt")
    if bb:
        ax.plot([t + off, t + off], [y0, y1], color=TRAIT_B, lw=3.0, zorder=3, solid_capstyle="butt")


def held(ax, t0, t1, who, trait):
    y = Y[who] + (0.13 if trait == "a" else -0.13)
    col = TRAIT_A if trait == "a" else TRAIT_B
    ax.hlines(y, t0, t1, color=col, lw=3.2, zorder=3)


def reset(ax, t, who):
    ax.plot(t, Y[who], marker="x", ms=13, mew=2.6, color=INK, zorder=5)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", help="figure directory (default: figures/)")
    args = parser.parse_args()
    out = output_dir(args.out, DEFAULT_FIGURES)
    target = check_writable(out / "schematic.png", DEFAULT_FIGURES)
    plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.facecolor": "white",
                         "figure.facecolor": "white"})
    fig, (L, R) = plt.subplots(1, 2, figsize=(12.4, 4.3))
    fig.subplots_adjust(left=0.01, right=0.995, top=0.9, bottom=0.08, wspace=0.08)

    # ---- left: bundled copying, b = 1 ------------------------------------
    lane(L, "Bundled copying (b = 1)\nfew encounters, each copies every trait")
    reset(L, 6.3, "1")                          # person 1 dies: both traits lost there
    for tr in "ab":
        held(L, 0.3, 6.3, "1", tr)
    encounter(L, 2.4, "1", "2", True, True)     # one shared encounter carries both
    for tr in "ab":
        held(L, 2.4, 10, "2", tr)
    encounter(L, 4.6, "3", "1")                 # donor 3 holds nothing yet
    encounter(L, 7.6, "2", "3", True, True)
    for tr in "ab":
        held(L, 7.6, 10, "3", tr)
    L.text(5.0, 2.36, "both traits share the same encounters", ha="center", fontsize=10, color=MUTED)

    # ---- right: partial copying, b = 1/4 ---------------------------------
    lane(R, "Partial copying (b = 1/4)\n4\u00d7 as many encounters, each copies a trait with probability 1/4")
    reset(R, 6.3, "1")
    for tr in "ab":
        held(R, 0.3, 6.3, "1", tr)
    # potential encounters (grey) and which trait (if any) passes
    events = [(0.9, "1", "2", False, False), (1.6, "1", "3", False, False),
              (2.4, "1", "2", True, False), (3.1, "2", "3", False, False),
              (3.8, "1", "3", False, True), (4.6, "3", "1", False, False),
              (5.2, "2", "1", False, False), (5.8, "1", "2", False, False),
              (7.0, "2", "3", True, False), (7.7, "3", "2", False, True),
              (8.4, "1", "2", False, False), (9.1, "3", "1", False, False)]
    for t, s, d, a, bb in events:
        encounter(R, t, s, d, a, bb)
    held(R, 2.4, 10, "2", "a")
    held(R, 7.0, 10, "3", "a")
    held(R, 3.8, 10, "3", "b")
    held(R, 7.7, 10, "2", "b")
    R.text(5.0, 2.36, "the two traits mostly travel through different encounters",
           ha="center", fontsize=10, color=MUTED)

    handles = [Line2D([], [], color=TRAIT_A, lw=3.2, label="trait A held / copied"),
               Line2D([], [], color=TRAIT_B, lw=3.2, label="trait B held / copied"),
               Line2D([], [], color=MUTED, lw=1.1, marker=">", ms=6, label="encounter (donor → learner)"),
               Line2D([], [], color=INK, lw=0, marker="x", ms=10, mew=2.4, label="reset (death; repertoire lost)")]
    fig.legend(handles=handles, loc="lower center", ncols=4, frameon=False, fontsize=10.5,
               bbox_to_anchor=(0.5, -0.04))
    fig.savefig(target, dpi=300, bbox_inches="tight", pad_inches=0.12,
                metadata={"Title": "Bundled versus partial copying", "Description": "Schematic; not data."})
    plt.close(fig)
    print(target)


if __name__ == "__main__":
    main()
