#!/usr/bin/env python3
"""Leave-one-date-out jackknife of the capped-shell subdiffusive exponent alpha.

For each of the seven capped replicates in turn, drop that replicate, refit the
joint model (P and alpha shared; k and s per replicate) to the remaining six
with the same multi-start procedure fit_groups.py uses, and record alpha and P.
The spread of alpha across exclusions shows how much any single replicate
drives the exponent.

Scope: each exclusion is a multi-start point fit (fit_groups.fit_group); it
does NOT repeat the profile-likelihood polish that fit_groups.py applies to the
full-data fit. The full-data refit reported first is the like-for-like
reference for these numbers.

Runtime: ~20 minutes on 32 cores (each fit is a parallel multi-start).

Usage:
    python scripts/jackknife_capped.py
Output:
    results/jackknife_capped.json
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fit_groups as fg
from common import load_alldata, result_path


def main():
    alldata = load_alldata()
    groups, _ = fg.classify_replicates(alldata)
    capped_keys = groups["capped"]
    print("Capped replicates:", capped_keys, flush=True)

    t0 = time.time()
    full = fg.fit_group([alldata[d][s] for d, s in capped_keys], alpha_fixed=None)
    out = dict(
        full_dataset=dict(alpha=float(full.params["alpha"].value),
                          P_nm_per_s=float(np.exp(full.params["logP"].value)),
                          rss=float(full.chisqr)),
        exclusions=[])
    print("Full-dataset (all %d) multi-start refit: alpha=%.4f  P=%.3g nm/s  (%.0f s)"
          % (len(capped_keys), out["full_dataset"]["alpha"],
             out["full_dataset"]["P_nm_per_s"], time.time() - t0), flush=True)

    for key in capped_keys:
        subset = [alldata[d][s] for d, s in capped_keys if (d, s) != key]
        t0 = time.time()
        r = fg.fit_group(subset, alpha_fixed=None)
        rec = dict(excluded=key[0], alpha=float(r.params["alpha"].value),
                   P_nm_per_s=float(np.exp(r.params["logP"].value)),
                   rss=float(r.chisqr))
        out["exclusions"].append(rec)
        print("excl=%-11s alpha=%.4f  P=%.3g nm/s  (%.0f s)"
              % (key[0], rec["alpha"], rec["P_nm_per_s"], time.time() - t0), flush=True)

    alphas = [e["alpha"] for e in out["exclusions"]]
    out["alpha_min"], out["alpha_max"] = float(min(alphas)), float(max(alphas))
    with open(result_path("jackknife_capped.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print("\nJackknife alpha range: %.3f to %.3f" % (out["alpha_min"], out["alpha_max"]))
    print("wrote results/jackknife_capped.json")


if __name__ == "__main__":
    main()
