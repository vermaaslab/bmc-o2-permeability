#!/usr/bin/env python3
"""Summary statistics quoted in the manuscript's Results text, computed from the
per-replicate stretched-exponential fits (results/abc_fits.json, written by
stretched_exponential_reporter.py).

Reports
  * the number of replicates per shell type and per group,
  * the geometric-mean characteristic decay time tau for capped and wiffle shells
    (tau spans a multiplicative range, so the geometric mean is the natural
    average) and their ratio,
  * Welch's t-test on log10(tau) and the Mann-Whitney U test on tau, capped vs
    wiffle (two-sided),
  * whether every capped replicate returns c < 1, and how many replicates
    return c = 1 (within the 0.999 classification threshold),
  * the range of the per-fit standard errors on c.

Usage:
    python scripts/tau_statistics.py
Output:
    results/tau_statistics.json
"""
import json
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import result_path


def geomean(x):
    return float(10 ** np.mean(np.log10(x)))


def main():
    with open(result_path("abc_fits.json")) as fh:
        rows = json.load(fh)
    capped = [r for r in rows if r["shell"] == "Capped"]
    wiffle = [r for r in rows if r["shell"] == "Wiffle"]
    tau_c = [r["tau"] for r in capped]
    tau_w = [r["tau"] for r in wiffle]

    welch = stats.ttest_ind(np.log10(tau_c), np.log10(tau_w), equal_var=False)
    mw = stats.mannwhitneyu(tau_c, tau_w, alternative="two-sided")
    c_err = [r["c_err"] for r in rows]
    out = dict(
        n_replicates=len(rows), n_capped=len(capped), n_wiffle=len(wiffle),
        n_by_group={g: sum(r["group"] == g for r in rows)
                    for g in ("capped", "wiffle_str", "wiffle_exp")},
        geomean_tau_capped_s=geomean(tau_c),
        geomean_tau_wiffle_s=geomean(tau_w),
        geomean_tau_ratio_capped_over_wiffle=geomean(tau_c) / geomean(tau_w),
        welch_t_log10_tau=dict(t=float(welch.statistic), p=float(welch.pvalue)),
        mann_whitney_tau=dict(U=float(mw.statistic), p=float(mw.pvalue)),
        all_capped_c_below_1=bool(all(r["c"] < 1 for r in capped)),
        capped_c_max=float(max(r["c"] for r in capped)),
        capped_c_mean=float(np.mean([r["c"] for r in capped])),
        capped_c_sd=float(np.std([r["c"] for r in capped], ddof=1)),
        n_replicates_with_c_equal_1=int(sum(r["group"] == "wiffle_exp" for r in rows)),
        c_standard_error_min=float(min(c_err)), c_standard_error_max=float(max(c_err)),
    )
    with open(result_path("tau_statistics.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    for k, v in out.items():
        print("%-40s %s" % (k, v))
    print("wrote results/tau_statistics.json")


if __name__ == "__main__":
    main()
