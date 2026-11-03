#!/usr/bin/env python3
"""Stretch exponent c by shell group (manuscript Fig. 5).

Reads results/abc_fits.json, written by stretched_exponential_reporter.py, so the
plotted exponents are exactly the ones tabulated -- the fits are never redone
here.

Usage: python scripts/plot_stretch_exponent.py
Output: results/stretch_exponent_by_group.png
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# All text sizes are multiplied by TEXT_SCALE so the labels stay legible at
# journal column width. Raise it for larger text.
TEXT_SCALE = 1.5

# Plot style: sans-serif text, outward ticks, open (unboxed) axes, 300 dpi.
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9 * TEXT_SCALE, "axes.labelsize": 9 * TEXT_SCALE,
    "axes.titlesize": 9 * TEXT_SCALE, "legend.fontsize": 8 * TEXT_SCALE,
    "xtick.labelsize": 7 * TEXT_SCALE, "ytick.labelsize": 7 * TEXT_SCALE,
    "axes.linewidth": 0.6,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": 300, "figure.dpi": 300,
})

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import result_path

GROUP_ORDER = ["capped", "wiffle_str", "wiffle_exp"]
GROUP_LABEL = {"capped": "Capped\n(intact shell)",
               "wiffle_str": "Wiffle,\nstretched",
               "wiffle_exp": "Wiffle,\nsingle\nexponential"}
COLORS = {"capped": "C0", "wiffle_str": "C1", "wiffle_exp": "C2"}
MARKERS = {"capped": "o", "wiffle_str": "s", "wiffle_exp": "^"}
META_GREY = "#666666"


def main():
    with open(result_path("abc_fits.json")) as fh:
        rows = json.load(fh)
    by_group = {g: [r["c"] for r in rows if r["group"] == g] for g in GROUP_ORDER}

    fig, ax = plt.subplots(figsize=(4.4, 3.4))
    rng = np.random.default_rng(1)
    for i, g in enumerate(GROUP_ORDER):
        vals = by_group[g]
        if not vals:
            continue
        x = i + rng.uniform(-0.07, 0.07, size=len(vals))
        ax.scatter(x, vals, s=42, color=COLORS[g], marker=MARKERS[g],
                   zorder=3, edgecolor="white", linewidth=0.6)
        ax.hlines(np.median(vals), i - 0.22, i + 0.22,
                  color=COLORS[g], lw=2.0, zorder=2)

    ax.axhline(1.0, color=META_GREY, ls=":", lw=1.0, zorder=1)

    ax.set_xticks(range(len(GROUP_ORDER)))
    ax.set_xticklabels([GROUP_LABEL[g] for g in GROUP_ORDER])
    ax.set_xlim(-0.5, 2.42)
    ax.set_ylabel("Stretch exponent $c$")
    ax.margins(y=0.12)

    fig.tight_layout()
    out = result_path("stretch_exponent_by_group.png")
    fig.savefig(out, dpi=300)
    print("wrote", out)


if __name__ == "__main__":
    main()
