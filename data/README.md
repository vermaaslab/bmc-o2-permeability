# Data

## `ethanol_780nm_traces.csv`

The 11 absorbance traces used in every analysis: the raw absorbance at 780 nm
of ICG (added from an ethanol stock to 5% of the final volume) in samples of
SOPP3-conjugated shells under continuous 473 nm laser irradiation. One row per
time point (7,988 rows).

| column | meaning |
|---|---|
| `replicate` | date key used throughout the manuscript, `MM.DD.YY`; a trailing `a` marks a second sample of the same shell type on the same date |
| `shell` | `Capped` (BMC-P pentamers added) or `Wiffle` (open shells) |
| `time_s` | spectrometer time stamp (s), starting at 0 |
| `absorbance_780nm` | raw absorbance at 780 nm (no baseline subtraction or normalization) |

Values are written with round-trip precision, so loading this file
(`scripts/common.py: load_alldata`) reproduces the originally parsed arrays
exactly.

| replicate | shell | points | t range (s) | median dt, first 100 pts (s) | median dt, last 100 pts (s) | A780 at t0 | A780 at end |
|---|---|---:|---|---:|---:|---:|---:|
| 01.22.26 | Capped | 727 | 0.0-1201.8 | 0.58 | 2.08 | 0.2782 | 0.1599 |
| 01.22.26a | Capped | 727 | 0.0-1201.9 | 0.58 | 2.08 | 0.2298 | 0.1255 |
| 10.02.25 | Capped | 726 | 0.0-1201.5 | 0.58 | 2.08 | 0.2745 | 0.1734 |
| 11.12.25 | Capped | 726 | 0.0-1201.5 | 0.58 | 2.08 | 0.2741 | 0.1823 |
| 11.12.25a | Capped | 726 | 0.0-1200.9 | 0.58 | 2.08 | 0.1420 | 0.0550 |
| 12.10.25 | Capped | 726 | 0.0-1201.6 | 0.58 | 2.08 | 0.2539 | 0.2062 |
| 12.10.25a | Capped | 726 | 0.0-1201.6 | 0.58 | 2.08 | 0.0762 | 0.0389 |
| 01.22.26 | Wiffle | 726 | 0.0-1201.3 | 0.58 | 2.08 | 0.1241 | 0.0337 |
| 10.02.25 | Wiffle | 726 | 0.0-1201.4 | 0.58 | 2.08 | 0.0674 | 0.0238 |
| 11.12.25 | Wiffle | 726 | 0.0-1200.6 | 0.58 | 2.08 | 0.0390 | 0.0055 |
| 12.10.25 | Wiffle | 726 | 0.0-1201.5 | 0.58 | 2.08 | 0.1645 | 0.0627 |

## `raw_file_manifest.csv`

What happened to **every** raw spectrometer file, including those that were
not analyzed. Written by `scripts/parse_raw_spectra.py`.

| date | raw file | status | replicate key | shell |
|---|---|---|---|---|
| 01.22.26 | `CAPPED_1_formatted_output.txt` | analyzed | 01.22.26 | Capped |
| 01.22.26 | `CAPPED_2_formatted_output.txt` | analyzed | 01.22.26a | Capped |
| 01.22.26 | `WIFFLE_1_formatted_output.txt` | rejected: negative absorbance | - | - |
| 01.22.26 | `WIFFLE_2_formatted_output.txt` | analyzed | 01.22.26 | Wiffle |
| 10.02.25 | `CAPPED_STIR_TEMP_formatted_output.txt` | analyzed | 10.02.25 | Capped |
| 10.02.25 | `WIFFLE_STIR_TEMP_formatted_output.txt` | analyzed | 10.02.25 | Wiffle |
| 11.12.25 | `CAPPED_10X_formatted_output.txt` | analyzed | 11.12.25 | Capped |
| 11.12.25 | `CAPPED_2X_formatted_output.txt` | analyzed | 11.12.25a | Capped |
| 11.12.25 | `ICG_Only_formatted_output.txt` | not analyzed: control (not a capped/wiffle shell sample) | - | - |
| 11.12.25 | `WIFFLE_formatted_output.txt` | analyzed | 11.12.25 | Wiffle |
| 12.10.25 | `CAPPED_1_formatted_output.txt` | analyzed | 12.10.25 | Capped |
| 12.10.25 | `CAPPED_2_formatted_output.txt` | analyzed | 12.10.25a | Capped |
| 12.10.25 | `ICG_Only_2_formatted_output.txt` | not analyzed: control (not a capped/wiffle shell sample) | - | - |
| 12.10.25 | `ICG_Only_formatted_output.txt` | not analyzed: control (not a capped/wiffle shell sample) | - | - |
| 12.10.25 | `WIFFLE_2_DYE_Spiked_formatted_output.txt` | analyzed | 12.10.25 | Wiffle |

* **Rejected:** a trace containing any negative absorbance is rejected as
  unphysical rather than fit (Methods, "Computational Modeling"). One trace
  was rejected this way.
* **Not analyzed:** files named `ICG_Only*` do not contain `CAPPED` or `WIFFLE`
  in their names and are therefore not treated as shell samples.
* **Keys:** the first sample of a shell type on a date takes the date as its
  key; a second sample of the same shell type on the same date takes
  `<date>a`. Files are visited in sorted order, so the assignment is
  deterministic.

## Raw spectra (not in this repository)

The 15 raw files (`<date>/<name>_formatted_output.txt`, about 16 MB each, about
240 MB in total) hold full spectra; only the 780 nm row is used. They are not
stored here because of their size. To confirm that the CSV reproduces from
them:

```
python scripts/parse_raw_spectra.py --raw-dir /path/to/inputdata
git diff --stat data/        # no change expected
```

## License

The files in this directory are licensed under CC BY 4.0; see
[`../LICENSE-CC-BY-4.0.txt`](../LICENSE-CC-BY-4.0.txt).
