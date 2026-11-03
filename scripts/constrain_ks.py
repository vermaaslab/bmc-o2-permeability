"""Ask whether the per-replicate nuisance parameters k and s really need the
freedom the unconstrained fit gives them.

Three variants per group, same model and same multi-start procedure otherwise:
  free      -- k_i and s_i fully independent per replicate (what Table 2 uses)
  bounded   -- k_i = k_base * m_i and s_i = s_base * n_i, with m_i, n_i in
               [0.5, 1.5]: every replicate within +/-50% of a shared group value
               (so at most a 3x spread between any two replicates)
  shared    -- one k and one s for the whole group (the strictest bracket)

Bounding does NOT reduce the number of free parameters, so this is not an
F-test; the question is simply how much RSS degrades and whether the bounded
parameters pile up against their limits (which would mean the constraint is
fighting the data rather than regularizing it).

Run AFTER fit_groups.py: the 'free' reference values are read from its saved
fit objects (results/cache/fit_groups_results.pkl).

Usage:
    python scripts/constrain_ks.py
Output:
    results/constrain_ks.json
"""
import json
import multiprocessing as mp
import os
import pickle
import sys

import numpy as np
from lmfit import Parameters, minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fit_groups as fg
from common import CACHE_DIR, load_alldata, result_path

SPREAD = 0.5  # +/-50%

alldata = load_alldata()
with open(os.path.join(CACHE_DIR, 'fit_groups_results.pkl'), 'rb') as fh:
    D = pickle.load(fh)
groups, free_results = D['groups'], D['results']
traces = {g: [alldata[d][s] for d, s in groups[g]] for g in groups}


def build(mode, tr, P0, a0, alpha_fixed):
    p = Parameters()
    p.add('A', value=float(np.mean([y[0] for _, y in tr])), min=1e-9)
    p.add('logP', value=np.log(P0), min=np.log(1e-9), max=np.log(1e12))
    if alpha_fixed is None:
        p.add('alpha', value=a0, min=0.05, max=1.0)
    else:
        p.add('alpha', value=alpha_fixed, vary=False)

    lo, hi = 1.0 - SPREAD, 1.0 + SPREAD
    if mode == 'free':
        for i in range(len(tr)):
            p.add('k_%d' % i, value=1e-3, min=0)
            p.add('log_s_%d' % i, value=np.log(max(1000.0, fg.S_FLOOR * 2)),
                  min=np.log(fg.S_FLOOR), max=np.log(1e7))
    elif mode == 'bounded':
        p.add('k_base', value=1e-3, min=0)
        # keep every derived log_s_i above the physical floor even at the
        # loosest allowed offset
        p.add('log_s_base', value=np.log(max(1000.0, fg.S_FLOOR * 2)),
              min=np.log(fg.S_FLOOR) - np.log(lo), max=np.log(1e7) - np.log(hi))
        for i in range(len(tr)):
            p.add('k_mult_%d' % i, value=1.0, min=lo, max=hi)
            p.add('k_%d' % i, expr='k_base*k_mult_%d' % i)
            p.add('s_mult_%d' % i, value=1.0, min=lo, max=hi)
            p.add('log_s_%d' % i, expr='log_s_base+log(s_mult_%d)' % i)
    elif mode == 'shared':
        p.add('k_base', value=1e-3, min=0)
        p.add('log_s_base', value=np.log(max(1000.0, fg.S_FLOOR * 2)),
              min=np.log(fg.S_FLOOR), max=np.log(1e7))
        for i in range(len(tr)):
            p.add('k_%d' % i, expr='k_base')
            p.add('log_s_%d' % i, expr='log_s_base')
    return p


def one(job):
    mode, g, P0, a0 = job
    tr = traces[g]
    alpha_fixed = 1.0 if g == 'wiffle_exp' else None
    try:
        p = build(mode, tr, P0, a0, alpha_fixed)
        r = minimize(fg.group_objective, p, args=(tr,),
                     method='least_squares', max_nfev=60000)
        return (mode, g, float(r.chisqr), {k: v.value for k, v in r.params.items()})
    except Exception as exc:
        return (mode, g, float('inf'), None)




def main():
    jobs = [(m, g, P0, a0)
            for m in ('bounded', 'shared')
            for g in ('capped', 'wiffle_str', 'wiffle_exp')
            for P0 in (4.0, 40.0, 150.0, 1000.0)
            for a0 in (0.45, 0.7, 0.9)]

    with mp.Pool(min(len(jobs), 30)) as pool:
        out = pool.map(one, jobs)

    best = {}
    for mode, g, rss, vals in out:
        if vals is None:
            continue
        if (mode, g) not in best or rss < best[(mode, g)][0]:
            best[(mode, g)] = (rss, vals)

    report = {}
    for g in ('capped', 'wiffle_str', 'wiffle_exp'):
        tr = traces[g]
        n = len(tr)
        free_rss = float(free_results[g].chisqr)
        entry = dict(n=n, free_rss=free_rss,
                     free_k=[free_results[g].params['k_%d' % i].value for i in range(n)],
                     free_s=[float(np.exp(free_results[g].params['log_s_%d' % i].value))
                             for i in range(n)])
        for mode in ('bounded', 'shared'):
            if (mode, g) not in best:
                continue
            rss, vals = best[(mode, g)]
            d = dict(rss=rss, ratio_vs_free=rss / free_rss,
                     P=float(np.exp(vals['logP'])), alpha=float(vals['alpha']),
                     k=[float(vals['k_%d' % i]) for i in range(n)],
                     s=[float(np.exp(vals['log_s_%d' % i])) for i in range(n)])
            if mode == 'bounded':
                km = [float(vals['k_mult_%d' % i]) for i in range(n)]
                sm = [float(vals['s_mult_%d' % i]) for i in range(n)]
                tol = 1e-3
                d['k_mult'] = km
                d['s_mult'] = sm
                d['k_at_bound'] = sum(1 for v in km
                                      if abs(v - (1 - SPREAD)) < tol or abs(v - (1 + SPREAD)) < tol)
                d['s_at_bound'] = sum(1 for v in sm
                                      if abs(v - (1 - SPREAD)) < tol or abs(v - (1 + SPREAD)) < tol)
            entry[mode] = d
        report[g] = entry

    with open(result_path('constrain_ks.json'), 'w') as fh:
        json.dump(report, fh, indent=2)
    print("wrote results/constrain_ks.json")


if __name__ == '__main__':
    main()
