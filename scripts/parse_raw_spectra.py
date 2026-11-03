#!/usr/bin/env python3
"""Extract the 780 nm absorbance traces from the raw spectrometer files.

Input : a directory laid out as  <raw-dir>/<date>/<NAME>_formatted_output.txt
        (one text file per laser experiment; first row = time stamps, first
        column = wavelength, entries = absorbance).
Output: data/ethanol_780nm_traces.csv   -- the 11 traces used in every analysis
        data/raw_file_manifest.csv      -- what happened to EVERY raw file

The raw files are ~16 MB each (~240 MB for the full set) and are not stored in
this repository; the CSV written here is what every downstream script reads.
You only need to run this script if you have the raw files and want to confirm
that the CSV reproduces from them.

Selection rules (identical to those described in the Methods):
  1. The row whose wavelength is exactly 780 nm is extracted together with the
     time-stamp row.
  2. Repeated time stamps (the spectrometer sometimes writes the same stamp
     twice) are removed by keeping only strictly increasing times.
  3. A trace containing ANY negative absorbance is rejected as unphysical, not
     fit.
  4. File names decide the shell type: 'CAPPED' -> Capped, 'WIFFLE' -> Wiffle.
     Anything else (e.g. the ICG-only controls) is recorded but not analyzed.
  5. The first sample of a shell type on a date is keyed by the date
     (e.g. '11.12.25'); a second sample of the SAME shell type on the SAME date
     is keyed '<date>a'. Files are visited in sorted order, so this assignment
     is deterministic.

Usage:
    python scripts/parse_raw_spectra.py --raw-dir /path/to/inputdata
"""
import argparse
import csv
import glob
import os

import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import DATA_DIR, MANIFEST_CSV, TRACES_CSV


def parse(raw_dir):
    files = sorted(glob.glob(os.path.join(raw_dir, "*", "*txt")))
    if not files:
        raise SystemExit("no <date>/*.txt files found under %s" % raw_dir)
    alldata, manifest = {}, []
    for f in files:
        date = os.path.basename(os.path.dirname(f))
        name = os.path.basename(f)
        data = np.loadtxt(f)
        row = np.nonzero(data[:, 0] == 780)[0][0]
        d780 = np.vstack([data[0][1:], data[row][1:]])
        n_raw = d780.shape[1]
        keep = np.hstack([0, 1 + np.nonzero(np.diff(d780[0]) > 0)[0]])
        d780 = d780[:, keep]
        entry = dict(date=date, raw_file=name, shell="", replicate="",
                     n_points_raw=n_raw, n_points_used=d780.shape[1])
        if np.any(d780 < 0):
            entry.update(status="rejected: negative absorbance")
        elif "CAPPED" in name:
            shell = "Capped"
        elif "WIFFLE" in name:
            shell = "Wiffle"
        else:
            shell = None
            entry.update(status="not analyzed: control (not a capped/wiffle shell sample)")
        if "status" not in entry:
            rep = date
            if rep in alldata and shell in alldata[rep]:
                rep = date + "a"
            alldata.setdefault(rep, {})[shell] = d780
            entry.update(status="analyzed", shell=shell, replicate=rep)
        manifest.append(entry)
    return alldata, manifest


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--raw-dir", required=True,
                    help="directory containing <date>/*_formatted_output.txt")
    args = ap.parse_args()

    alldata, manifest = parse(args.raw_dir)
    os.makedirs(DATA_DIR, exist_ok=True)
    with open(TRACES_CSV, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["replicate", "shell", "time_s", "absorbance_780nm"])
        for rep in sorted(alldata):
            for shell in sorted(alldata[rep]):
                t, a = alldata[rep][shell]
                for ti, ai in zip(t, a):
                    w.writerow([rep, shell, repr(float(ti)), repr(float(ai))])
    with open(MANIFEST_CSV, "w", newline="") as fh:
        cols = ["date", "raw_file", "status", "replicate", "shell",
                "n_points_raw", "n_points_used"]
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(manifest)
    n = sum(len(v) for v in alldata.values())
    print("analyzed traces: %d  (capped %d, wiffle %d)" % (
        n, sum("Capped" in v for v in alldata.values()),
        sum("Wiffle" in v for v in alldata.values())))
    for m in manifest:
        print("  %-9s %-34s %s %s" % (m["date"], m["raw_file"], m["status"],
                                      ("-> " + m["replicate"] + " " + m["shell"]) if m["replicate"] else ""))
    print("wrote", TRACES_CSV)
    print("wrote", MANIFEST_CSV)


if __name__ == "__main__":
    main()
