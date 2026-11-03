#!/usr/bin/env python3
"""Fit the subdiffusive two-compartment permeability model to each shell group
INDEPENDENTLY, assess how well each parameter is actually identified, and emit
the LaTeX for Table 2 of the manuscript.

DESIGN
------
Each shell group (capped; wiffle, stretched; wiffle, single exponential) goes
through the SAME procedure, and Table 2 is generated from the saved fit
objects rather than transcribed. The columns of Table 2 are therefore directly
comparable to each other, and to the model-free stretched exponents c in
Table 1 (see the "implied c" diagnostic below).

The model is Eqs. S5-S7 of the Supplementary Information; see docs/model.md
for the equations and for where each one lives in this file.

IMPLEMENTATION NOTES
--------------------
1. USE method='least_squares', NOT 'leastsq'. On this problem lmfit's default
   leastsq returns its own starting values essentially unchanged (nfev ~ 35 for
   28 parameters, ier=2, "Fit succeeded"). It is not fitting; it declares
   convergence on the parameter-tolerance test before exploring. Switching to
   least_squares (trust-region reflective) improved RSS by ~49x on an otherwise
   identical problem. A "successful" leastsq fit here is not evidence of
   anything.

2. FIT log(P), NOT P. The permeabilities span many orders of magnitude between
   groups (~1 nm/s for capped, ~1e5+ for the single-exponential wiffles). A
   large linear P has a near-zero residual gradient, which stalls the solver
   globally -- not just for that one parameter. Same reasoning as the log_s
   reparameterization of s (log_s_i) used below.

3. MULTI-START IS MANDATORY, AND STILL NOT SUFFICIENT. The RSS surface is rough:
   profiling P (re-optimizing alpha at each fixed P) routinely finds points
   BELOW the best multi-start minimum. This script therefore treats the profile
   as part of the fit -- it polishes from the best profile point and iterates --
   and reports profiled RANGES, not just point estimates. If you only report the
   point estimate from this model you will overstate your precision.

Outputs:
    results/cache/fit_groups_results.pkl  -- dict of lmfit MinimizerResult per
                                             group + profiles (git-ignored; needs
                                             matching scipy/lmfit to unpickle)
    results/fit_groups_summary.json       -- plain-text summary of every number
                                             in Table 2
    results/fit_groups_plotdata.json      -- group assignment + per-replicate
                                             model-free numbers (read by the
                                             plotting scripts)
    results/table2_odeparams.tex          -- \\input{}-able LaTeX for Table 2

Usage:
    python scripts/fit_groups.py            # fit, profile, write everything
    python scripts/fit_groups.py --quick    # fewer starts/profile points, smoke test
"""

import argparse
import json
import multiprocessing as mp
import os
import pickle
import sys
import time

import numpy as np
from lmfit import Parameters, minimize
from scipy.integrate import odeint
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import CACHE_DIR, load_alldata, result_path

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
SHELL_RADIUS_NM = 19.0
# s = V_outside/V_inside is fit as log_s with a floor derived from the shell
# stock concentration the collaborator reported (150 nM before desalting and
# dilution into the 1 mL reaction). Dilution only ever LOWERS concentration,
# which only ever RAISES s, so 150 nM sets a hard floor on s and nothing else.
SHELL_STOCK_CONC_NM = 150.0
_V_SHELL_L = (4.0 / 3.0) * np.pi * SHELL_RADIUS_NM ** 3 * 1e-24  # nm^3 -> L
_N_A = 6.022e23
S_FLOOR = 1.0 / (SHELL_STOCK_CONC_NM * 1e-9 * _N_A * _V_SHELL_L)

# EPS_SUBDIFF must be 1e-3, NOT a "more accurate" tiny value like 1e-12.
# P_eff(t) = P*(t+eps)^(alpha-1) diverges as t->0 for alpha<1, so eps sets the
# size of an unphysical flux spike in the very first ODE step. At alpha=0.3,
# eps=1e-12 gives (1e-12)^-0.7 ~ 1e8.4 -- an enormous spurious flux that
# caused ODEintWarning stiffness failures and NaN residuals during
# development. eps=1e-3 keeps that prefactor O(10) while staying far below the
# first measurement time (t[0] >= 1 s), so it cannot affect the fitted region.
# The profile seeds below deliberately probe alpha as low as 0.3, which is
# exactly the regime that breaks with a smaller eps.
EPS_SUBDIFF = 1.0e-3
ABSERR, RELERR = 1.0e-12, 1.0e-10


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------
def subdiffusive_ode(t, wvec, pvec):
    """Two-compartment dye model with a time-dependent permeability.

    P_eff(t) = P * t^(alpha-1); alpha=1 recovers constant-P Fickian transport.
    The 3/19nm prefactor is the shell's surface-area-to-volume ratio (sphere,
    r=19nm): S/V_i = 3/r. It converts a permeability (nm/s) into a rate (1/s).
    """
    mi, mo = wvec
    k, p, log_s, alpha = pvec
    s = np.exp(log_s)
    p_eff = p * (t + EPS_SUBDIFF) ** (alpha - 1.0)
    return [-k * mi + 3.0 * p_eff * (mo - mi) / SHELL_RADIUS_NM,
            -3.0 * p_eff * (mo - mi) / (s * SHELL_RADIUS_NM)]


def model_curve(t, A, alpha, k, p, log_s, minit):
    """Absorbance predicted by the model: volume-weighted mean of the two
    compartment concentrations, A * ([M_i] + s[M_o]) / (s+1)."""
    y = odeint(subdiffusive_ode, [minit, minit], t,
               args=([k, p, log_s, alpha],), atol=ABSERR, rtol=RELERR, tfirst=True)
    s = np.exp(log_s)
    return A * (y[:, 0] + s * y[:, 1]) / (s + 1.0)


def stretched_exp(t, a, b, c):
    return a * np.exp(-b * t ** c)


# ---------------------------------------------------------------------------
# Group classification -- derived from the data, not hardcoded
# ---------------------------------------------------------------------------
def classify_replicates(alldata):
    """Assign every replicate to a group by fitting a stretched exponential to
    the raw trace. Wiffle replicates whose unconstrained fit returns c=1 (to
    within tol) are 'wiffle_exp'; the rest are 'wiffle_str'. This classification
    is an OUTPUT of the model-free fits, never an input to the ODE fits."""
    groups = {"capped": [], "wiffle_str": [], "wiffle_exp": []}
    per_replicate = {}
    for date in sorted(alldata):
        for shell in sorted(alldata[date]):
            t, y = alldata[date][shell]
            popt, pcov = curve_fit(stretched_exp, t, y, p0=[y[0], 1.0, 0.9],
                                   bounds=(0, [np.inf, 100.0, 1.0]), maxfev=40000)
            c, c_err = popt[2], float(np.sqrt(np.diag(pcov))[2])
            tau = popt[1] ** (-1.0 / c)
            if shell == "Capped":
                g = "capped"
            else:
                g = "wiffle_exp" if c > 0.999 else "wiffle_str"
            groups[g].append((date, shell))
            per_replicate[(date, shell)] = dict(group=g, c=c, c_err=c_err, tau=tau,
                                                a=popt[0], b=popt[1])
    return groups, per_replicate


# ---------------------------------------------------------------------------
# Group fitting
# ---------------------------------------------------------------------------
def build_params(traces, logP0, alpha0, alpha_fixed):
    """One shared logP (and alpha, unless fixed) per group; k and log_s per
    replicate as nuisance parameters."""
    p = Parameters()
    p.add("A", value=float(np.mean([y[0] for _, y in traces])), min=1e-9)
    p.add("logP", value=np.log(logP0), min=np.log(1e-9), max=np.log(1e12))
    if alpha_fixed is None:
        p.add("alpha", value=alpha0, min=0.05, max=1.0)
    else:
        p.add("alpha", value=alpha_fixed, vary=False)
    for i, (t, y) in enumerate(traces):
        p.add("k_%d" % i, value=1e-3, min=0)
        p.add("log_s_%d" % i, value=np.log(max(1000.0, S_FLOOR * 2)),
              min=np.log(S_FLOOR), max=np.log(1e7))
    return p


def group_objective(params, traces):
    A = params["A"].value
    P = np.exp(params["logP"].value)
    alpha = params["alpha"].value
    out = []
    for i, (t, y) in enumerate(traces):
        fit = model_curve(t, A, alpha, params["k_%d" % i].value, P,
                          params["log_s_%d" % i].value, y[0] / A)
        out.append(fit - y)
    return np.concatenate(out)


def _run_one_start(job):
    """Worker for the parallel multi-start. Returns (rss, params_dict) rather than
    the MinimizerResult because lmfit results do not pickle cleanly across
    processes; the caller re-fits briefly from the winning values."""
    traces, P0, a0, alpha_fixed = job
    p = build_params(traces, P0, a0 if a0 else 0.8, alpha_fixed)
    try:
        r = minimize(group_objective, p, args=(traces,),
                     method="least_squares", max_nfev=60000)
        return (float(r.chisqr), {k: v.value for k, v in r.params.items()})
    except Exception:
        return (float("inf"), None)


def fit_group(traces, alpha_fixed=None, P_starts=(1, 10, 100, 1e3, 1e5),
              alpha_starts=(0.4, 0.6, 0.8, 0.95), label="", nproc=None):
    """Multi-start fit, parallelized over starts.

    Always least_squares -- see note 1 in the module docstring. The starts are
    independent, so they run across cores; on a 32-core box this turns a
    multi-hour serial sweep into minutes.
    """
    a_seeds = [None] if alpha_fixed is not None else list(alpha_starts)
    jobs = [(traces, P0, a0, alpha_fixed) for P0 in P_starts for a0 in a_seeds]
    nproc = nproc or min(len(jobs), max(1, (os.cpu_count() or 2) - 1))

    if nproc > 1 and len(jobs) > 1:
        with mp.Pool(nproc) as pool:
            out = pool.map(_run_one_start, jobs)
    else:
        out = [_run_one_start(j) for j in jobs]

    rss, vals = min(out, key=lambda x: x[0])
    if vals is None or not np.isfinite(rss):
        raise RuntimeError("all starts failed for group %s" % label)

    # Re-fit from the winning start's converged values to recover a full
    # MinimizerResult (with covariance) in this process.
    p = build_params(traces, 1.0, 0.8, alpha_fixed)
    for k, v in vals.items():
        if k in p:
            p[k].value = v
    return minimize(group_objective, p, args=(traces,),
                    method="least_squares", max_nfev=60000)


def _run_one_profile(job):
    """Worker for the parallel profile scan."""
    traces, base_vals, Pf, a0, alpha_fixed = job
    p = build_params(traces, 1.0, 0.8, alpha_fixed)
    for k, v in base_vals.items():
        if k in p:
            p[k].value = v
    p["logP"].value = np.log(Pf)
    p["logP"].vary = False
    if a0 is not None:
        p["alpha"].value = a0
    try:
        r = minimize(group_objective, p, args=(traces,),
                     method="least_squares", max_nfev=40000)
        return (float(Pf), float(r.chisqr), float(r.params["alpha"].value))
    except Exception:
        return (float(Pf), float("inf"), float("nan"))


def profile_P(base, traces, P_vals, alpha_fixed=None,
              alpha_seeds=(0.3, 0.55, 0.78, 0.95), nproc=None):
    """Fix P, re-optimize everything else (alpha included), in parallel.

    This is the honest uncertainty on P: the covariance stderr at the minimum
    badly understates it because P and alpha trade off against each other.
    """
    seeds = [None] if alpha_fixed is not None else list(alpha_seeds)
    base_vals = {k: v.value for k, v in base.params.items()}
    jobs = [(traces, base_vals, Pf, a0, alpha_fixed) for Pf in P_vals for a0 in seeds]
    nproc = nproc or min(len(jobs), max(1, (os.cpu_count() or 2) - 1))

    if nproc > 1 and len(jobs) > 1:
        with mp.Pool(nproc) as pool:
            out = pool.map(_run_one_profile, jobs)
    else:
        out = [_run_one_profile(j) for j in jobs]

    by_P = {}
    for Pf, rss, al in out:
        if Pf not in by_P or rss < by_P[Pf][0]:
            by_P[Pf] = (rss, al)
    return [dict(P=float(Pf), rss=float(by_P[Pf][0]), alpha=float(by_P[Pf][1]))
            for Pf in sorted(by_P)]


def fit_with_profile_polish(traces, alpha_fixed=None, rounds=3, quick=False, label=""):
    """Fit, profile, and if the profile finds a better point than the fit did,
    restart from there and repeat. Necessary because the surface is rough enough
    that multi-start alone lands above the profile floor (note 3)."""
    P_starts = (1, 100, 1e4) if quick else (1, 10, 100, 1e3, 1e5)
    a_starts = (0.6, 0.9) if quick else (0.4, 0.6, 0.8, 0.95)
    npts = 8 if quick else 16
    a_seeds = (0.6, 0.9) if quick else (0.3, 0.55, 0.78, 0.95)

    best = fit_group(traces, alpha_fixed, P_starts, a_starts, label)
    prof = None
    for _ in range(rounds):
        P_hi = 1e7 if alpha_fixed is not None else 1e5
        prof = profile_P(best, traces, np.logspace(0, np.log10(P_hi), npts),
                         alpha_fixed, a_seeds)
        finite = [d for d in prof if np.isfinite(d["rss"])]
        if not finite:
            break
        floor = min(finite, key=lambda d: d["rss"])
        if floor["rss"] >= best.chisqr * (1 - 1e-6):
            break  # profile found nothing better; converged
        p = best.params.copy()
        p["logP"].value = np.log(floor["P"])
        p["logP"].vary = True
        if alpha_fixed is None and np.isfinite(floor["alpha"]):
            p["alpha"].value = floor["alpha"]
        r = minimize(group_objective, p, args=(traces,),
                     method="least_squares", max_nfev=60000)
        if r.chisqr < best.chisqr:
            best = r
        else:
            break
    return best, prof


def implied_c(result, traces):
    """Forward-simulate the fitted model and fit a stretched exponential to the
    RESULT, giving the c the model predicts. This is the check that matters for
    the paper: alpha and c are different parameterizations of the same curvature,
    so a fit that reproduces the data must imply a c close to the measured one.

    A fit with low RSS but an implied c far from the data's c means alpha is not
    identified -- the per-replicate k and s are absorbing curvature that alpha
    would otherwise have to explain. The ODE model has three shape knobs per
    replicate (alpha, k, s) against the stretched exponential's two (b, c), so
    this degeneracy is structural, not a fitting artifact.
    """
    A = result.params["A"].value
    P = np.exp(result.params["logP"].value)
    alpha = result.params["alpha"].value
    out = []
    for i, (t, y) in enumerate(traces):
        curve = model_curve(t, A, alpha, result.params["k_%d" % i].value, P,
                            result.params["log_s_%d" % i].value, y[0] / A)
        try:
            popt, _ = curve_fit(stretched_exp, t, curve, p0=[curve[0], 1.0, 0.9],
                                bounds=(0, [np.inf, 100.0, 1.0]), maxfev=40000)
            out.append(float(popt[2]))
        except Exception:
            out.append(float("nan"))
    return out


def near_optimal_band(prof, rss_min, tol=1.25):
    """Points whose RSS is within tol of the minimum. With a rough surface and a
    strong P-alpha trade-off this band, not the stderr, is what the data
    constrain."""
    band = [d for d in prof if np.isfinite(d["rss"]) and d["rss"] <= rss_min * tol]
    if not band:
        return None
    return dict(P_lo=min(d["P"] for d in band), P_hi=max(d["P"] for d in band),
                alpha_lo=min(d["alpha"] for d in band),
                alpha_hi=max(d["alpha"] for d in band), n=len(band))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="fewer starts and profile points; smoke test only")
    args = ap.parse_args()

    alldata = load_alldata()

    groups, per_replicate = classify_replicates(alldata)
    print("Replicate classification (from per-replicate stretched-exponential fits):")
    for g in ("capped", "wiffle_str", "wiffle_exp"):
        for date, shell in groups[g]:
            info = per_replicate[(date, shell)]
            print("  %-11s %-7s -> %-11s c=%.4f tau=%7.0f s"
                  % (date, shell, g, info["c"], info["tau"]))

    traces = {g: [alldata[d][s] for d, s in groups[g]] for g in groups}

    results, profiles, bands, implied = {}, {}, {}, {}
    for g in ("capped", "wiffle_str", "wiffle_exp"):
        if not traces[g]:
            print("\n[%s] no replicates; skipping" % g)
            continue
        # The single-exponential wiffles are classical BY CONSTRUCTION -- their
        # defining property is c=1, so fitting alpha there would be circular.
        alpha_fixed = 1.0 if g == "wiffle_exp" else None
        t0 = time.time()
        print("\n[%s] fitting %d replicate(s)%s..."
              % (g, len(traces[g]), " (alpha fixed at 1)" if alpha_fixed else ""))
        best, prof = fit_with_profile_polish(traces[g], alpha_fixed,
                                             quick=args.quick, label=g)
        results[g], profiles[g] = best, prof
        bands[g] = near_optimal_band(prof, best.chisqr) if prof else None
        P = np.exp(best.params["logP"].value)
        print("  RSS=%.7f  P=%.4g nm/s  alpha=%.4f  (%.1f s, %d varying params)"
              % (best.chisqr, P, best.params["alpha"].value, time.time() - t0,
                 best.nvarys))
        if bands[g]:
            b = bands[g]
            print("  within 25%% of min RSS: P %.3g-%.3g nm/s, alpha %.3f-%.3f (%d pts)"
                  % (b["P_lo"], b["P_hi"], b["alpha_lo"], b["alpha_hi"], b["n"]))

        # Does this fit actually reproduce the observed curve shape?
        ic = implied_c(best, traces[g])
        dc = [per_replicate[key]["c"] for key in groups[g]]
        implied[g] = ic
        print("  implied c from fitted model: %s" % ", ".join("%.3f" % v for v in ic))
        print("  measured c on same traces:   %s" % ", ".join("%.3f" % v for v in dc))
        gap = abs(float(np.nanmean(ic)) - float(np.mean(dc)))
        print("  mean |implied c - measured c| = %.3f%s" % (
            gap, "   <-- MISMATCH: alpha is not identified" if gap > 0.05 else "   (consistent)"))

        # classical comparison for the two subdiffusive groups
        if alpha_fixed is None:
            cl, _ = fit_with_profile_polish(traces[g], 1.0, quick=args.quick,
                                            label=g + "-classical")
            results[g + "_classical"] = cl
            print("  classical (alpha=1) RSS=%.7f  -> improvement %.1fx"
                  % (cl.chisqr, cl.chisqr / best.chisqr))

    write_outputs(results, profiles, bands, per_replicate, groups, implied)


def write_outputs(results, profiles, bands, per_replicate, groups, implied):
    with open(os.path.join(CACHE_DIR, "fit_groups_results.pkl"), "wb") as fh:
        pickle.dump(dict(results=results, profiles=profiles, bands=bands,
                         per_replicate=per_replicate, groups=groups,
                         implied_c=implied), fh)

    # fit_groups_results.pkl embeds lmfit MinimizerResult objects, which pull in
    # scipy's asteval symbol table (scipy.special._special_ufuncs and friends).
    # Unpickling it therefore requires the SAME scipy internals to exist in
    # whatever environment reads it back -- it broke across a scipy version
    # bump. Anything that only needs the group assignment
    # and the per-replicate model-free numbers (e.g. plot_decay_and_tau.py)
    # should read THIS file instead: plain dicts/lists/floats, no third-party
    # classes, portable across any scipy/lmfit version and even across
    # languages.
    plotdata = dict(
        groups={g: [[d, s] for d, s in lst] for g, lst in groups.items()},
        per_replicate={"%s|%s" % (d, s): dict(v) for (d, s), v in per_replicate.items()},
    )
    with open(result_path("fit_groups_plotdata.json"), "w") as fh:
        json.dump(plotdata, fh, indent=2)

    summary = {}
    for g, r in results.items():
        summary[g] = dict(rss=float(r.chisqr), nvarys=int(r.nvarys),
                          P_nm_per_s=float(np.exp(r.params["logP"].value)),
                          alpha=float(r.params["alpha"].value),
                          alpha_stderr=(float(r.params["alpha"].stderr)
                                        if r.params["alpha"].stderr else None))
        if g in bands and bands[g]:
            summary[g]["near_optimal_band"] = bands[g]
    # model-free exponents, for the alpha-vs-c comparison the table now makes
    for g in ("capped", "wiffle_str", "wiffle_exp"):
        cs = [per_replicate[key]["c"] for key in groups[g]]
        if cs:
            summary.setdefault(g, {})["data_c_mean"] = float(np.mean(cs))
            summary[g]["data_c_sd"] = float(np.std(cs, ddof=1)) if len(cs) > 1 else None
            summary[g]["data_c_values"] = [float(c) for c in cs]
            if g in implied:
                summary[g]["implied_c_values"] = implied[g]
                summary[g]["implied_c_mean"] = float(np.nanmean(implied[g]))
                summary[g]["implied_vs_measured_c_gap"] = abs(
                    float(np.nanmean(implied[g])) - float(np.mean(cs)))
    with open(result_path("fit_groups_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)

    with open(result_path("table2_odeparams.tex"), "w") as fh:
        fh.write(render_table(results, bands, summary))
    print("\nWrote results/cache/fit_groups_results.pkl, results/fit_groups_summary.json, "
          "results/fit_groups_plotdata.json, results/table2_odeparams.tex")


def _fmt(v, prec=2):
    if v is None or not np.isfinite(v):
        return "---"
    if abs(v) >= 1e4 or (0 < abs(v) < 1e-2):
        exp = int(np.floor(np.log10(abs(v))))
        return r"$%.1f{\times}10^{%d}$" % (v / 10 ** exp, exp)
    return "$%.*f$" % (prec, v)


def render_table(results, bands, summary):
    cols = [("capped", "Capped"), ("wiffle_str", "Wiffle, stretched"),
            ("wiffle_exp", "Wiffle, single exp.")]
    present = [(g, h) for g, h in cols if g in results]

    L = [r"% Generated by fit_groups.py -- do not edit by hand.",
         r"\begin{table}", r"\centering",
         r"\caption{\label{tab:odeparams}Parameters of the subdiffusive "
         r"permeability model (Eqs.~\ref{eq:concentrations}--\ref{eq:peff}), fit "
         r"to each shell group independently by the same procedure: $P$ and "
         r"$\alpha$ shared within the group, $k$ and $s$ fit per replicate, "
         r"multi-start over $P$ and $\alpha$, followed by profiling $P$ and "
         r"restarting from the profile floor. Note that $\alpha$ is \emph{not} "
         r"the same quantity as the stretched-exponential exponent $c$ of "
         r"Table~\ref{tab:perreplicate} and the two should not be compared "
         r"directly: $\alpha$ is the exponent of the permeation \emph{rate law}, "
         r"$P_{\mathrm{eff}}\propto t^{\alpha-1}$, whereas $c$ describes the "
         r"resulting \emph{absorbance decay}, and the map between them depends "
         r"on the per-replicate $k$ and $s$. The meaningful comparison is the "
         r"third row: we forward-simulate each fitted model and fit a stretched "
         r"exponential to the result, giving the $c$ the model predicts, which "
         r"should match the measured $c$. Because $P$ and $\alpha$ trade off, "
         r"the near-optimal band is a better guide to what the data constrain "
         r"than the covariance standard error.}",
         r"\begin{tabular}{l" + "c" * len(present) + "}", r"\toprule",
         " & " + " & ".join(h for _, h in present) + r" \\",
         " & " + " & ".join("($n=%d$)" % len(summary[g]["data_c_values"])
                            for g, _ in present) + r" \\",
         r"\midrule"]

    row = [r"$\alpha$ (fitted)"]
    for g, _ in present:
        s = summary[g]
        if g == "wiffle_exp":
            row.append(r"$1$ (by construction)")
        elif s.get("alpha_stderr"):
            row.append(r"$%.3f\pm%.3f$" % (s["alpha"], s["alpha_stderr"]))
        else:
            row.append("$%.3f$" % s["alpha"])
    L.append(" & ".join(row) + r" \\")

    row = [r"$c$ implied by fitted model"]
    for g, _ in present:
        s = summary[g]
        row.append("$%.3f$" % s["implied_c_mean"] if s.get("implied_c_mean") is not None
                   else "---")
    L.append(" & ".join(row) + r" \\")

    row = [r"$c$ measured (Table~\ref{tab:perreplicate})"]
    for g, _ in present:
        s = summary[g]
        if s.get("data_c_sd"):
            row.append(r"$%.3f\pm%.3f$" % (s["data_c_mean"], s["data_c_sd"]))
        else:
            row.append("$%.3f$" % s["data_c_mean"])
    L.append(" & ".join(row) + r" \\")

    L.append(" & ".join([r"$P$ at $t=1$\,s (nm/s)"] +
                        [_fmt(summary[g]["P_nm_per_s"]) for g, _ in present]) + r" \\")

    row = [r"$P_{\mathrm{eff}}$ at $t=1200$\,s (nm/s)"]
    for g, _ in present:
        s = summary[g]
        row.append(_fmt(s["P_nm_per_s"] * 1200.0 ** (s["alpha"] - 1.0)))
    L.append(" & ".join(row) + r" \\")

    L += [r"\midrule",
          (r"\multicolumn{%d}{l}{\emph{Near-optimal band (RSS within 25\%% of the "
           r"minimum, $\alpha$ re-optimized at each $P$):}} \\") % (len(present) + 1)]

    row = [r"$P$ (nm/s)"]
    for g, _ in present:
        b = bands.get(g)
        row.append("%s--%s" % (_fmt(b["P_lo"], 1), _fmt(b["P_hi"], 1)) if b else "---")
    L.append(" & ".join(row) + r" \\")

    row = [r"$\alpha$"]
    for g, _ in present:
        b = bands.get(g)
        row.append(("$%.2f$--$%.2f$" % (b["alpha_lo"], b["alpha_hi"]))
                   if b and g != "wiffle_exp" else "---")
    L.append(" & ".join(row) + r" \\")

    L += [r"\midrule",
          " & ".join([r"RSS, subdiffusive"] +
                     ["$%.5f$" % results[g].chisqr for g, _ in present]) + r" \\"]

    row = [r"RSS, classical ($\alpha=1$)"]
    for g, _ in present:
        key = g + "_classical"
        row.append("$%.5f$" % results[key].chisqr if key in results
                   else "$%.5f$" % results[g].chisqr)
    L.append(" & ".join(row) + r" \\")

    row = [r"Improvement from $\alpha<1$"]
    for g, _ in present:
        key = g + "_classical"
        row.append(r"$%.1f\times$" % (results[key].chisqr / results[g].chisqr)
                   if key in results else "none (by construction)")
    L.append(" & ".join(row) + r" \\")

    L += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    sys.exit(main())
