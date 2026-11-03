"""Shared paths and data loading for the analysis scripts.

Every script in this directory imports its input/output locations from here, so
the repository can be cloned anywhere and run from any working directory.

    data/      inputs (tidy CSV of the 780 nm traces + a provenance manifest)
    results/   everything the scripts write (tables, figures, JSON summaries)
    results/cache/   large intermediate files (pickles); git-ignored
"""
import csv
import os

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(REPO, "data")
RESULTS_DIR = os.path.join(REPO, "results")
CACHE_DIR = os.path.join(RESULTS_DIR, "cache")

TRACES_CSV = os.path.join(DATA_DIR, "ethanol_780nm_traces.csv")
MANIFEST_CSV = os.path.join(DATA_DIR, "raw_file_manifest.csv")

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(CACHE_DIR, exist_ok=True)


def result_path(name):
    """Path of an output file inside results/."""
    return os.path.join(RESULTS_DIR, name)


def load_alldata(path=TRACES_CSV):
    """Read the 780 nm traces.

    Returns {replicate: {shell: array}} where `replicate` is the date key used
    throughout the manuscript (e.g. '11.12.25', or '11.12.25a' for a second
    sample of the same shell type on the same day), `shell` is 'Capped' or
    'Wiffle', and the array has shape (2, n): row 0 is time in seconds, row 1
    is the raw absorbance at 780 nm.
    """
    cols = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            key = (row["replicate"], row["shell"])
            t, a = cols.setdefault(key, ([], []))
            t.append(float(row["time_s"]))
            a.append(float(row["absorbance_780nm"]))
    alldata = {}
    for (rep, shell), (t, a) in cols.items():
        alldata.setdefault(rep, {})[shell] = np.vstack([t, a])
    return alldata
