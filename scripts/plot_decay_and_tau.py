#!/usr/bin/env python3
"""Two-panel figure: log-absorbance decay curves (left) and characteristic
decay time by group (right), both colored consistently by shell group
(manuscript Fig. 4).

Reads the group classification from results/fit_groups_plotdata.json (written
by fit_groups.py) so the grouping here always matches what Table 2 uses -- it
is never redefined independently.

NOTE ON .pkl vs .json: results/cache/fit_groups_results.pkl embeds lmfit
MinimizerResult objects, which pull in scipy-internal modules at unpickle time,
so reading it back requires the same scipy/lmfit versions that wrote it. This
script only needs the group assignment and per-replicate model-free numbers
(c, tau, ...), so it reads the .json instead: plain dicts/lists/floats, no
third-party classes, portable across environments.

Usage: python scripts/plot_decay_and_tau.py
Output: results/decay_and_tau_by_group.png and .pdf
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
from common import load_alldata, result_path

GROUP_LABEL = {"capped": "Capped", "wiffle_str": "Wiffle, stretched",
              "wiffle_exp": "Wiffle, single exp."}
GROUP_ORDER = ["capped", "wiffle_str", "wiffle_exp"]
COLORS = {"capped": "C0", "wiffle_str": "C1", "wiffle_exp": "C2"}
MARKERS = {"capped": "o", "wiffle_str": "s", "wiffle_exp": "^"}


def main():
    alldata = load_alldata()
    with open(result_path("fit_groups_plotdata.json")) as fh:
        D = json.load(fh)
    groups = {g: [tuple(pair) for pair in lst] for g, lst in D["groups"].items()}
    per_rep = D["per_replicate"]
    tau_by_group = {g: [per_rep["%s|%s" % (date, shell)]["tau"] for date, shell in groups[g]]
                    for g in groups}

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(8.6, 3.8))

    for g in GROUP_ORDER:
        for date, shell in groups[g]:
            t, y = alldata[date][shell]
            axL.plot(t, np.log(y / y[0]), color=COLORS[g], lw=1.2, alpha=0.8)
    axL.set_xlim(left=0)
    axL.set_xlabel("Time (s)")
    axL.set_ylabel(r"$\ln$(normalized absorbance)")
    axL.margins(x=0.02, y=0.06)

    rng = np.random.default_rng(3)
    for i, g in enumerate(GROUP_ORDER):
        vals = tau_by_group[g]
        x = i + rng.uniform(-0.08, 0.08, size=len(vals))
        axR.scatter(x, vals, s=42, color=COLORS[g], marker=MARKERS[g],
                    edgecolor="white", linewidth=0.6, zorder=3)
        axR.hlines(np.median(vals), i - 0.22, i + 0.22, color=COLORS[g],
                   lw=2.0, zorder=2)
    axR.set_yscale("log")
    axR.set_xticks(range(len(GROUP_ORDER)))
    axR.set_xticklabels([GROUP_LABEL[g] for g in GROUP_ORDER],
                        rotation=12, ha="right")
    axR.set_xlim(-0.5, len(GROUP_ORDER) - 0.5)
    axR.set_ylabel(r"Characteristic decay time $\tau$ (s)")
    axR.margins(y=0.15)

    handles = [plt.Line2D([0], [0], color=COLORS[g], marker=MARKERS[g],
                          linestyle="-", markeredgecolor="white",
                          markeredgewidth=0.6, label=GROUP_LABEL[g])
              for g in GROUP_ORDER]
    axL.legend(handles=handles, loc="lower left", frameon=False, fontsize=7 * TEXT_SCALE)

    fig.tight_layout()
    out = result_path("decay_and_tau_by_group.png")
    fig.savefig(out, dpi=300)
    fig.savefig(result_path("decay_and_tau_by_group.pdf"), dpi=300)
    print("wrote", out)


if __name__ == "__main__":
    main()
