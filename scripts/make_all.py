#!/usr/bin/env python3
"""Regenerate every table, figure and number in the manuscript from the data
in data/.

Pipeline (each step reads the previous steps' outputs from results/):

  0. parse_raw_spectra.py              raw spectra -> data/*.csv    (optional,
                                       needs the raw files; --raw-dir)
  1. stretched_exponential_reporter.py per-trace fits  -> Table 1, abc_fits.json
  2. plot_stretch_exponent.py          Fig. 5          -> stretch_exponent_by_group.png
  3. fit_groups.py                     ODE fits        -> Table 2, fit_groups_*.json
  4. plot_decay_and_tau.py             Fig. 4          -> decay_and_tau_by_group.png/.pdf

Robustness checks (optional, slow; --checks), which need step 3 first:

  5. jackknife_capped.py               leave-one-date-out alpha -> jackknife_capped.json
  6. constrain_ks.py                   bounded/shared k and s   -> constrain_ks.json

Step 3 is the expensive one (multi-start + profiling; ~8 minutes on 32 cores).
Use --skip-fits to reuse the committed ODE-fit outputs and redo only the cheap
steps.

Usage:
    python scripts/make_all.py                    # steps 1-4
    python scripts/make_all.py --skip-fits        # steps 1, 2, 4
    python scripts/make_all.py --checks           # steps 1-6
    python scripts/make_all.py --raw-dir /path/to/inputdata   # also step 0
"""
import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common import DATA_DIR, RESULTS_DIR, TRACES_CSV  # noqa: E402

# (script, description, files it produces, relative to results/)
STEPS = [
    ("stretched_exponential_reporter.py", "per-trace stretched-exponential fits (Table 1)",
     ["abc_fits.json", "abctable.tex", "ctautable.tex"]),
    ("plot_stretch_exponent.py", "Fig. 5, stretch exponent by group",
     ["stretch_exponent_by_group.png"]),
    ("fit_groups.py", "grouped ODE permeability fits (Table 2)",
     ["fit_groups_plotdata.json", "fit_groups_summary.json", "table2_odeparams.tex"]),
    ("plot_decay_and_tau.py", "Fig. 4, decay curves and decay times",
     ["decay_and_tau_by_group.png", "decay_and_tau_by_group.pdf"]),
]
CHECKS = [
    ("jackknife_capped.py", "leave-one-date-out jackknife of capped alpha", ["jackknife_capped.json"]),
    ("constrain_ks.py", "bounded / shared k and s sensitivity", ["constrain_ks.json"]),
]


def run(script, description, extra=()):
    print("\n=== %s: %s ===" % (script, description), flush=True)
    result = subprocess.run([sys.executable, os.path.join(HERE, script), *extra], cwd=HERE)
    if result.returncode != 0:
        sys.exit("FAILED: %s (exit %d)" % (script, result.returncode))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--raw-dir", help="also re-parse the raw spectra in this directory")
    ap.add_argument("--skip-fits", action="store_true",
                    help="reuse existing ODE fit outputs (skip fit_groups.py)")
    ap.add_argument("--checks", action="store_true",
                    help="also run the slow robustness checks (jackknife, k/s bounds)")
    args = ap.parse_args()

    if args.raw_dir:
        run("parse_raw_spectra.py", "parse raw spectra", ["--raw-dir", args.raw_dir])
    if not os.path.exists(TRACES_CSV):
        sys.exit("missing %s -- run with --raw-dir, or restore data/" % TRACES_CSV)

    for script, description, _ in STEPS:
        if args.skip_fits and script == "fit_groups.py":
            print("--- skipping %s ---" % script)
            continue
        run(script, description)
    if args.checks:
        for script, description, _ in CHECKS:
            run(script, description)

    print("\nOutputs in %s:" % RESULTS_DIR)
    for _, _, produces in STEPS + (CHECKS if args.checks else []):
        for p in produces:
            ok = os.path.exists(os.path.join(RESULTS_DIR, p))
            print("  %s %s" % ("ok     " if ok else "MISSING", p))


if __name__ == "__main__":
    main()
