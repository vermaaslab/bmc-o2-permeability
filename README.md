# bmc-o2-permeability

[![DOI](https://zenodo.org/badge/1444898684.svg)](https://doi.org/10.5281/zenodo.23933032)

Analysis code and data for the manuscript *Construction and characterization of
synthetic bacterial microcompartment shells containing oxygen-sensing
reporters*. The analysis estimates how permeable synthetic bacterial
microcompartment (BMC) shells are to a photobleachable dye, from the kinetics of
dye bleaching by SOPP3-generated singlet oxygen.

## The experiment and the analysis in brief

SOPP3, an FMN-binding photosensitizing protein, is tethered inside the lumen of
*wiffle ball* shells (open vertices) and, for *capped* samples, the vertices are
closed with BMC-P pentamers. Under 473 nm laser light SOPP3 turns oxygen into
singlet oxygen, which bleaches indocyanine green (ICG) loaded into the sample.
The absorbance of ICG at 780 nm is recorded for 20 minutes. In the model, oxygen is
in large excess everywhere but the singlet oxygen is made, and consumed by dye,
inside the shell, so dye must cross the shell to be bleached and the decay
rate reports on shell permeability.

The analysis has three layers:

1. **Model-free curve shape.** Each trace is fit on its own with a stretched
   exponential, $`A(t)=a\,e^{-bt^{c}}`$. $`c<1`$ means the decay is slower at late
   times than a single exponential predicts. The characteristic decay time is
   $`\tau=b^{-1/c}`$. (Table 1, Figure 5.)
2. **Mechanistic model.** A two-compartment (inside/outside the shell)
   diffusion-reaction model with a time-dependent permeability
   $`P_\mathrm{eff}(t)=P\,t^{\alpha-1}`$ is fit to each shell group, giving
   $`P`$ and $`\alpha`$ and honest uncertainty bands. (Table 2, Figure 4.)
3. **Robustness checks (not discussed in the current manuscript text).** A
   leave-one-date-out jackknife of $`\alpha`$, and a test of how much the
   per-replicate nuisance parameters matter.

The equations are in [`docs/model.md`](docs/model.md).

## Quick start

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt        # or: conda env create -f environment.yml

make quick     # seconds: Table 1, Figures 4-5, reusing the committed ODE-fit outputs
make           # ~10 min on 32 cores: also refits the ODE model (Table 2)
make checks    # ~30 min: everything, plus the jackknife and nuisance-parameter checks
```

Everything the scripts write goes to `results/`, and the committed copies there
are the outputs behind the manuscript, so a reviewer can read them without
running anything. All scripts can be run from any working directory.

## Where each manuscript item comes from

| Manuscript item | Script (in `scripts/`) | Output (in `results/`) |
|---|---|---|
| Table 1: per-replicate $`c`$ and $`\tau`$ | `stretched_exponential_reporter.py` | `ctautable.tex` (also `abctable.tex` with $`a,b,c`$; `abc_fits.json`) |
| Quoted statistics: geometric-mean $`\tau`$, Welch and Mann-Whitney tests, "$`c<1`$ for every capped replicate" | `tau_statistics.py` | `tau_statistics.json` |
| Figure 4: decay curves and decay times | `plot_decay_and_tau.py` | `decay_and_tau_by_group.png`, `.pdf` |
| Figure 5: stretch exponent by group | `plot_stretch_exponent.py` | `stretch_exponent_by_group.png` |
| Table 2: ODE model parameters, near-optimal bands, classical-vs-subdiffusive RSS | `fit_groups.py` | `table2_odeparams.tex`, `fit_groups_summary.json` |
| Not in the current manuscript: leave-one-date-out jackknife of $`\alpha`$ (capped) | `jackknife_capped.py` | `jackknife_capped.json` |
| Not in the current manuscript: effect of bounding or sharing the per-replicate $`k`$ and $`s`$ | `constrain_ks.py` | `constrain_ks.json` |
| SI Eqs. S1-S7 | (see [`docs/model.md`](docs/model.md)) | |

`scripts/make_all.py` runs them in dependency order. The `.tex` tables are
generated and meant to be `\input{}` into the manuscript, not edited by hand.

## Results at a glance

All values below are read from `results/*.json`; nothing is typed by hand.

| Group | α | P at t = 1 s (nm/s) | P_eff at t = 1200 s (nm/s) | Near-optimal band for P (nm/s) | RSS improvement over α = 1 |
|---|---|---:|---:|---|---:|
| Capped (n=7) | 0.438 ± 0.008 | 146.7 | 2.73 | 21.5–464 | 23.1× |
| Wiffle, stretched (n=2) | 0.705 ± 0.010 | 39.8 | 4.92 | 10–464 | 22.6× |
| Wiffle, single exponential (n=2) | 1 (fixed) | 4.0 | 4.00 | 8.6–10⁷ | – |

* Characteristic decay time (geometric mean): capped 2891 s, wiffle 980 s (ratio 2.95); Welch t-test on log τ, p = 0.0052; Mann-Whitney, p = 0.0061.
* Stretch exponent c: below 1 for all 7 capped replicates (mean 0.77 ± 0.04 s.d.); exactly 1 for 2 of 11 traces, both wiffle.
* Leave-one-date-out jackknife of the capped exponent (not in the current manuscript): α = 0.325–0.525 across the 7 exclusions, versus 0.438 from the full-data multi-start refit; P ranges 80–483 nm/s (`results/jackknife_capped.json`).
* For the single-exponential wiffle group P is bounded only from below (the band runs to the upper search limit), and its "P_eff" is just P because α = 1 by construction.

## Repository layout

```
data/
  ethanol_780nm_traces.csv      the 11 analyzed traces (raw absorbance at 780 nm)
  raw_file_manifest.csv         what happened to every raw spectrometer file
  README.md                     column definitions, per-trace summary, exclusions
scripts/
  common.py                     paths and the CSV loader shared by all scripts
  parse_raw_spectra.py          raw spectra -> data/*.csv (needs the raw files)
  stretched_exponential_reporter.py   model-free per-trace fits (Table 1)
  tau_statistics.py             statistics quoted in the Results text
  plot_stretch_exponent.py      Figure 5
  fit_groups.py                 ODE model: multi-start fit + profiling (Table 2)
  plot_decay_and_tau.py         Figure 4
  jackknife_capped.py           leave-one-date-out check on alpha
  constrain_ks.py               bounded / shared k and s check
  make_all.py                   runs the above in order
results/                        committed outputs (tables, figures, JSON summaries)
docs/model.md                   the model, its parameters, and how it is fit
requirements.txt  environment.yml  Makefile
```

## Reading the results

* **What the data determine.** The decay is not a single exponential for 9 of
  the 11 traces, and the stretched exponent $`c<1`$ is returned for every capped
  replicate. In the ODE model the exponent $`\alpha<1`$ is robustly supported for
  the capped and stretched-wiffle groups (the classical $`\alpha=1`$ fit is far
  worse; see the table above).
* **What they do not.** $`P`$ and $`\alpha`$ trade off, so the point estimate of
  $`P`$ is much better constrained than its covariance standard error suggests.
  Quote the near-optimal band (Table 2), not the standard error. The capped and
  stretched-wiffle bands overlap almost entirely, and the single-exponential
  wiffle group yields only a lower bound on $`P`$. Details in
  [`docs/model.md`](docs/model.md), section 6.

## Reproducibility

* No fit uses random numbers: the multi-start grids are fixed. (The small
  horizontal jitter of points in Figures 4 and 5 uses a fixed seed.) Re-running
  `fit_groups.py` here reproduced the committed tables and JSON values exactly.
  The jackknife agreed with an earlier run to within 5×10⁻⁵ in α, so expect
  agreement to roughly that level rather than bit-for-bit for the jackknife.
* Results were produced with Python 3.11, numpy 2.4.6, scipy 1.17.1,
  lmfit 1.3.4 and matplotlib 3.11.0 (`requirements.txt`). Other versions should
  reproduce the findings but may change the trailing digits of the ODE fits.
* The raw spectra (about 240 MB) are not stored here; `data/README.md` explains
  how to confirm that the CSV reproduces from them.
* `results/cache/` holds a pickle of the full lmfit objects. It is git-ignored,
  needs matching scipy/lmfit versions to read, and is regenerated by
  `fit_groups.py`; nothing in the committed results depends on it except the
  optional `constrain_ks.py`.

## Status

The manuscript is in preparation and this repository will be updated with it.

## License

| Path | License | Text |
|---|---|---|
| `data/`, `results/` (traces, tables, figures, JSON summaries), `docs/`, and this README | Creative Commons Attribution 4.0 (CC BY 4.0) | [`LICENSE-CC-BY-4.0.txt`](LICENSE-CC-BY-4.0.txt) |
| everything else (`scripts/`, `Makefile`, `requirements.txt`, `environment.yml`, `.gitignore`) | GNU General Public License, version 3 | [`LICENSE`](LICENSE) |

If you reuse the data, figures or documentation, please credit the authors of the associated
manuscript (Dwyer, Leslie, Yang, Ferlez, Beck, Vermaas and Kerfeld),
*Construction and characterization of synthetic bacterial microcompartment
shells containing oxygen-sensing reporters*, and cite the paper once it is
published.
