# Model and fitting procedure

This document states the mathematical model fit in `scripts/fit_groups.py` and
maps each equation to the code. Equation numbers (S1-S7) are those of the
Supplementary Information (SI).

## 1. Observable

The dye (indocyanine green, ICG) is monitored by its absorbance at 780 nm,
which is proportional to the monomer concentration $[M]$ (main text, Eq. 1).
Traces are fit as recorded: the raw absorbance, not a background-subtracted or
normalized version.

## 2. Simplest model: exponential decay (S1-S2)

If singlet oxygen is produced at a constant rate and consumes dye
bimolecularly,

```math
\frac{d[M]}{dt} = -k\,[M]\,[{}^{1}\mathrm{O}_2] \tag{S1}
```

and, with $[{}^{1}\mathrm{O}_2]$ constant (sample open to air, constant laser
power),

```math
\frac{d[M]}{dt} = -C\,[M] \quad\Rightarrow\quad [M](t) = [M]_0\,e^{-Ct} \tag{S2}
```

The prediction is a single exponential. Only 2 of the 11 traces are consistent
with it (`results/tau_statistics.json`; `results/abc_fits.json`).

## 3. Two-compartment permeation model (S3-S6)

Dye outside the shell crosses the shell with permeability $P$ (length/time) and
is consumed only inside, where SOPP3 is tethered. Writing Fick's first law,
$J = P\,A_\mathrm{shell}\,\Delta C$, for the number of dye molecules $M_i$,
$M_o$ inside and outside (the SI writes the shell area as $A$; it is written
$A_\mathrm{shell}$ here to keep it distinct from the absorbance scale factor
$A$ of S6):

```math
\begin{aligned}
\frac{dM_i}{dt} &= -k\,M_i - P A_\mathrm{shell}\left([M_i]-[M_o]\right)\\
\frac{dM_o}{dt} &= \phantom{-k\,M_i}\;P A_\mathrm{shell}\left([M_i]-[M_o]\right)
\end{aligned}\tag{S3}
```

Dividing by the compartment volumes, $V_i$ inside and $V_o = s\,V_i$ outside
($s \gg 1$), and treating the shell as a sphere of radius $r = 19$ nm,

```math
\frac{A_\mathrm{shell}}{V_i} = \frac{4\pi r^2}{\tfrac{4}{3}\pi r^3} = \frac{3}{r}
\approx 0.158\ \mathrm{nm^{-1}} \tag{S4}
```

gives the concentration equations that are actually integrated:

```math
\begin{aligned}
\frac{d[M_i]}{dt} &= -k\,[M_i] + \frac{3P}{19\ \mathrm{nm}}\left([M_o]-[M_i]\right)\\
\frac{d[M_o]}{dt} &= -\frac{3P}{s\cdot 19\ \mathrm{nm}}\left([M_o]-[M_i]\right)
\end{aligned}\tag{S5}
```

The spectrometer sees the whole cuvette, so the measured absorbance is the
volume-weighted mean of the two concentrations,

```math
A_{780\,\mathrm{nm}} = A\,\frac{[M_i] + s\,[M_o]}{s+1} \tag{S6}
```

with both compartments initially at the same concentration,
$[M_i](0) = [M_o](0) = A_{780}(t_0)/A$. For constant $P$ the solutions of S5
are again exponential, so a constant-$P$ model cannot reproduce the stretched
decays (SI, "Model Building").

## 4. Time-dependent (subdiffusive) permeability (S7)

To allow single-file, subdiffusive transport through the shell pores the
permeability is made time dependent,

```math
P_\mathrm{eff}(t) = P\,t^{\alpha-1} \tag{S7}
```

and inserted for $P$ in S5. $\alpha = 1$ recovers constant $P$; $\alpha < 1$ means
the effective permeability falls with time. $P$ is therefore the permeability
at $t = 1$ s, and the permeability at the end of a 1200 s run is
$P\,1200^{\alpha-1}$ (these two values are the "P at t = 1 s" and
"P_eff at t = 1200 s" rows of Table 2).

## 5. How the model is fit

| Parameter | Meaning | Scope | Bounds | In code |
|---|---|---|---|---|
| $A$ | concentration-to-absorbance scale (S6) | shared within a group | $> 10^{-9}$ | `A` |
| $P$ | permeability at $t=1$ s (nm/s) | shared within a group | $10^{-9}$ to $10^{12}$, fit as $\ln P$ | `logP` |
| $\alpha$ | permeability exponent (S7) | shared within a group | 0.05 to 1 (fixed at 1 for the single-exponential wiffle group) | `alpha` |
| $k$ | dye-consumption rate inside the shell (1/s) | per replicate | $\ge 0$ | `k_i` |
| $s$ | $V_o/V_i$ | per replicate | $s_\mathrm{floor}$ to $10^7$, fit as $\ln s$ | `log_s_i` |

**Groups.** Replicates are assigned to groups by the data, not by hand:
each trace is first fit with the stretched exponential
$a\,e^{-b t^{c}}$ (`classify_replicates`; Table 1). Capped replicates form the
`capped` group; wiffle replicates whose own unconstrained fit returns
$c > 0.999$ form `wiffle_exp`, the rest `wiffle_str`. The ODE model is then fit
to each group independently by the same procedure.

**Volume ratio.** $s$ is fit as $\ln s$ because $s/(s+1)$ is insensitive to $s$
once $s$ is large, so a linear parameterization drifts over orders of
magnitude. The lower bound is physical: shells were assembled at 150 nM before
desalting and dilution into the reaction, and dilution can only lower the
concentration (raise $s$), so
$s \ge s_\mathrm{floor} = 1/(150\ \mathrm{nM}\cdot N_A\cdot V_\mathrm{shell}) = 385.3$
for $r = 19$ nm (`S_FLOOR`).

**Regularization of the $t \to 0$ singularity.** For $\alpha < 1$,
$t^{\alpha-1}$ diverges at $t=0$. The code evaluates
$P\,(t+\epsilon)^{\alpha-1}$ with $\epsilon = 10^{-3}$ s (`EPS_SUBDIFF`), far
below the first sample time ($\ge 1$ s), so it cannot affect the fitted region
but avoids a spurious flux spike in the first integration step.

**Integration and optimization.** `scipy.integrate.odeint`
(`atol=1e-12`, `rtol=1e-10`) integrates S5 at each recorded time point; the
residuals of S6 against the raw trace are minimized with lmfit's
trust-region-reflective `method="least_squares"`. The default
Levenberg-Marquardt solver terminated at its starting values on this problem
(it satisfied its parameter-tolerance test without exploring), so it is not
used.

**Multi-start and profiling.** The surface in $(P,\alpha)$ is rough and $P$
and $\alpha$ trade off, so:

1. `fit_group` starts from a grid of $P_0 \in \{1, 10, 10^2, 10^3, 10^5\}$ nm/s
   and $\alpha_0 \in \{0.4, 0.6, 0.8, 0.95\}$ and keeps the lowest RSS.
2. `profile_P` fixes $P$ on a logarithmic grid (16 points, 1 to $10^5$ nm/s;
   to $10^7$ for the single-exponential group), re-optimizes everything else
   (including $\alpha$) from four $\alpha$ seeds, and records the RSS.
3. `fit_with_profile_polish` restarts from the best profile point if it beats
   the multi-start fit and repeats (at most 3 rounds).
4. `near_optimal_band` reports the range of $P$ and $\alpha$ over profile
   points whose RSS is within 25% of the minimum. **This band, not the
   covariance standard error, is what the data constrain**: $P$ and $\alpha$ are
   strongly anticorrelated, so the standard error understates the uncertainty
   in $P$ (see Table 2 and `results/fit_groups_summary.json`).

**Classical comparison.** Each subdiffusive group is refit with $\alpha=1$
(`<group>_classical` in the saved results); the ratio of the two RSS values is
the "Improvement from $\alpha<1$" row of Table 2.

**Consistency check against the model-free fits (`implied_c`).** $\alpha$ is
the exponent of the permeability *rate law* and $c$ is the exponent of the
resulting *absorbance decay*; the map between them depends on the per-replicate
$k$ and $s$, so the two are not interchangeable. For each fitted model the
code therefore forward-simulates S5-S6 at the fitted parameters, fits a
stretched exponential to that simulated curve, and reports the $c$ the model
implies. If the ODE fit reproduces the data, this implied $c$ matches the
measured $c$ of Table 1 (third row of Table 2).

## 6. Limits of what the fits determine

* The data determine robustly that $\alpha < 1$ for the capped and
  stretched-wiffle groups: the classical ($\alpha=1$) fit is worse by about
  23x in both groups, and the model-free exponents $c$ agree.
* The absolute value of $P$ is **not** tightly determined, because of the
  $P$-$\alpha$ trade-off; report and use the near-optimal band.
* For the single-exponential wiffle group, $\alpha=1$ by construction and
  only a lower bound on $P$ is obtained: the RSS stays flat out to the largest
  value tested.
* Because the capped and stretched-wiffle $P$ bands overlap almost entirely,
  this model and dataset cannot resolve a permeability difference between
  those two groups.
* The per-replicate $k$ and $s$ absorb some of the curvature that $\alpha$
  would otherwise explain; `scripts/constrain_ks.py` quantifies how much the
  fits degrade when $k$ and $s$ are bounded or shared.

## License

This document is licensed under CC BY 4.0; see
[`../LICENSE-CC-BY-4.0.txt`](../LICENSE-CC-BY-4.0.txt).
