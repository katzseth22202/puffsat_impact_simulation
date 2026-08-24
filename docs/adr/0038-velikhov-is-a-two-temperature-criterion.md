# The Velikhov criterion is a two-temperature criterion, and only the cold leg meets it

`electrothermal_screen` tested two links — Hall drive `beta > BETA_CRIT` and ionisation gain
`S > 1` — against a `BETA_CRIT = 2` that Q-F(b) recorded as "engineering practice taken on
authority", noting that "the real criterion is a linearised dispersion relation (Velikhov 1962;
Kerrebrock 1964) whose sources are not available here". The sources were obtained. **Both halves
of the screen were wrong**, in opposite directions, and the corrected verdict is narrower and
firmer than either: the hot legs are stable by three to four orders of magnitude, and the cold
leg's second half is unstable with a microsecond e-folding time.

The criterion is Petit & Geffray, *Non-Equilibrium Plasma Instabilities*, Acta Physica Polonica A
**115** (2009) 1170, restating Velikhov [1] via Petit & Valensi (1969):

```
s        = 2 k T_e^2 / [E_i (T_e - T_g)] * 1/(1 + 1.5 k T_e / E_i)
beta_cr  = 1.935 f + 0.065 + s          f = -(d mu/mu)/(d n_e/n_e)
g        = sigma E*^2 / [n_e (E_i + 1.5 k T_e)(1 + beta^2)] * (beta - beta_cr)
```

with the stated limits `beta_cr ~ (s^2 + 2s)^{1/2}` (weakly ionised) and `beta_cr ~ 2 + s` (fully
ionised). Implemented as `critical_hall_parameter`, `mobility_sensitivity`, `ionisation_energy`
and `electrothermal_loop`.

## Error 1 — the gain differentiated only the seed

`ionisation_sensitivity` claimed to be `S = d ln n_e / d ln T_e`, the gain of the feedback loop,
but `electron_density` evaluates water's contribution at the **gas** temperature and only the
potassium Saha constant at `t_e`. Above ~7000 K water supplies most of the electrons, so the
derivative saw almost nothing move and reported `S ~ 0`.

| `T` [K] | 2000 | 3000 | 5000 | 8000 | 11 000 | 15 000 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `S`, seed only (was) | 13.3 | 9.05 | 3.51 | **0.11** | **0.01** | **0.00** |
| `S`, all sources (now) | 13.3 | 9.05 | 3.52 | **2.18** | **6.48** | **5.68** |

The old column produced Q-F(b)'s "the gain collapses above ~5000 K … the hot end is safe on
**both** counts", which was **the headline stabilising result and it is not real**. Water at
13.6 eV against a hotter gas is as sensitive as the seed at 4.34 eV against a cool one:
`chi/(2kT) ~ 10` either way. The gain is not monotone — it dips where potassium is spent and
water has not started, then recovers. Two tests asserted the artifact and now assert the shape.

`electron_density` itself is unchanged: holding water at the gas temperature is the right
contract for `sigma`, whose non-equilibrium use here is the 2000–5000 K band where the seed
supplies the electrons. The gain gets a private `_n_e_all_sources_at` instead.

## Error 2 — `BETA_CRIT = 2` is the `s -> 0` limit, and `s` is never small here

**`(T_e - T_g)` sits in the denominator of `s`, so `beta_cr` diverges as the plasma approaches
thermal equilibrium.** Two is what `beta_cr ~ 2 + s` gives for a *strongly* two-temperature
plasma; the paper states the instability as "occurring in `T_e > T_g` regimes", and its Fig. 2
plots `beta_cr` running to infinity as `T_e -> T_g`. It also scopes itself to "low magnetic
Reynolds non-thermal plasmas". This plume is neither: Q-M established equilibrium, and `Rm` runs
~40–900.

So the whole question became: *how far out of equilibrium is the plume actually driven?* That is
Q-F's open item — "deciding needs an electron energy balance (Joule heating against elastic loss),
which needs an E-field and is a real piece of work" — and Q-F(b) noted it became writable once
`B` and `u` were known. Both are, so `electron_energy_balance` writes it:

```
Q_joule = j^2/sigma,  j = |curl B|/mu0 ~ B/(mu0 L_eff)
Q_loss  = 1.5 n_e k_B (delta nu) (T_e - T_g)
L_eff   = min(field-gradient scale, sqrt(t_transit/(mu0 sigma)))
```

Algebraic, not an ODE, on Q-F(b)'s own finding that electron relaxation is ~1e-8 s against a
~1e-3 s transit. Solved by damped fixed point: hotter electrons raise `sigma`, which both thickens
the current skin and cuts `j^2/sigma`, so the feedback is negative and it contracts.

**The driving current is the field-gradient current, not `sigma u B`.** At `Rm` of 40–900 the
field is largely frozen to the plasma, so the plasma-frame field is small. Taking the low-`Rm`
generator form `j = sigma u B` overstates the heating by ~3 orders of magnitude — it gives
`ΔT ~ 5e5 K` — and is simply the wrong regime.

## The verdict

`x_K = 1%`, local `B = 20/(A/A*)` T, equilibrium branch. **This is the elastic-only pass**; the
addendum below adds the inelastic channels, which move the numbers and not the verdict.

| leg | `T_g` at exit | `T_e - T_g` | `beta_cr` | `beta` | e-folding | verdict |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| 75 km/s | 16 224 K | **0.4 K** | 6 535 | 0.38 | — | stable, 4 decades |
| 65 | 14 151 | **1.2 K** | 1 851 | 0.58 | — | stable, 3 decades |
| 56.53 | 11 681 | **7.9 K** | 204 | 1.20 | — | stable, 2 decades |
| 45.58 | 4 597 | **813.5 K** | **2.07** | **4.63** | **3.2 µs** | **UNSTABLE\*** |

**Reproduce: `make analysis-electrothermal`** -> `data/results/electrothermal_scan.csv`.
`electrothermal.py` walks `electrothermal_loop` along `expansion.history()` under **both** readings
of the current-layer thickness, so Addendum 3's disagreement is an output rather than a paragraph.
The table above is the elastic-plus-inelastic (shipped) balance at the exit station of each leg;
the elastic-only column that the addendum contrasts against is in that addendum, not here.

**The resolution trade runs opposite to `expansion.history`'s.** The deliverable is a *threshold
crossing*, so the station spacing — not the integration step — is what has to be converged: at
`expansion.history`'s own default stride the crossing lands anywhere in `A/A*` 1.97–2.28 and the
dwell reads 1.8–2.1 ms. At `steps=160, stride=2` the crossing is pinned at `A/A* = 1.968`,
`T = 7569 K`, `t = 1.012 ms`, unchanged to four digits against `steps=320`.

**\*Qualified by Addendum 3: the cold-leg verdict is not robust.** It rests on the current-layer
thickness, and the alternative defensible choice makes every unstable station stable.

The cold leg crosses at `A/A* ~ 1.97` (7569 K) and is unstable to the exit — **1.70 ms of a
2.71 ms transit**, at e-folding times of 3–54 µs. Hundreds of e-foldings; brevity is not a
defence. The hot legs are not marginal: their `beta_cr` exceeds `beta` by 170× and 4000×.

**Q-F closes as a by-product.** The elevation is computed rather than assumed: +0 K on the hot
legs, +1055 K at the cold exit (+814 K once the addendum's inelastic channels are included).
Q-F's own table puts the conductivity cliff at 1841–2378 K for a +1000 K elevation and
2341–2884 K for +500 K — below every exit temperature on either figure — so **the cliff is not
reached on any leg**.

## Consequences

1. **Q-F(b)'s "the loop is CLOSED at 2000–5000 K … a live concern in this design"** — superseded.
   Right that the cool end is where risk lives, wrong that the bag is exposed: at bag density the
   plume is in equilibrium and `beta_cr` is enormous.
2. **Q-K's three levers are re-scoped.** Exit field, bag geometry and seed fraction all move
   `beta`, and `beta` is not what binds. What binds is `(T_e - T_g)`, i.e. **how cold the plume
   gets before it leaves the field**. Q-K's `BETA_CRIT` caveat — "the design sits on a boundary
   whose location is unknown" — is resolved: the boundary is known and it is not near the hot legs.
3. **Q-M's "no leg is exposed" holds for the hot legs and fails for the cold one**, against a
   criterion Q-M did not have.
4. **The fix is to keep the cold leg hot**, which is the `k` lever: `k ~ 4–5` on the cold tail
   puts the exit above ~11 000 K, where `beta_cr` is ~200 against `beta ~ 1.2`. Isp cost 2–7%
   (bare-plate ballistic model, and that is the tamper study's — it needs checking in `aim`).
   The mirror-stagnation enthalpy is worth ~1000 K and is **not** sufficient on its own, because
   the plume must traverse the unstable band regardless of where it ends.

## Considered and rejected: a loop-gain scaling argument

Before the source was obtained, the same physics was expressed as
`loop gain = [beta^2/(1+beta^2)] * S * (T_e - T_g)/T_e`, on the reasoning that a heating
perturbation can only move `T_e` by the fraction Joule heating supplies. It peaks at **0.176**
and pronounces every station stable by ≥5.7×.

**It is wrong, and it is wrong in the dangerous direction.** The direction was right — near
equilibrium is stabilising — but the O(1) coefficient was invented, and the published criterion
puts the cold leg's `beta` a factor 2.5 *above* `beta_cr` where the scaling put it 5.7× below.
Recorded because it is an easy derivation to re-invent, and because it is a clean case of a
plausible scaling argument being unable to decide a margin near unity.

## What would falsify this

- **`Q_en` and `ln Lambda`.** Both feed `nu`, hence `beta` and `sigma`. `sigma`'s docstring already
  says a factor-of-two claim is undecidable by this model; `beta` at the cold exit is 4.64 against
  `beta_cr = 1.84`, a factor 2.5, so a factor-2 error in `nu` is enough to flip that station. The
  *leg-level* verdict is safer than the station-level one — the crossing at 7569 K has
  `beta/beta_cr` = 1.04 and moves with any of this.
- **`L_eff` — this is the dominant lever, not a weak one. See Addendum 3.**
- ~~**The inelastic channels.**~~ **Closed below: they help and they are not enough.**
- **`f`.** Interpolated by Coulomb collision fraction between the paper's two stated limits rather
  than differentiated. It contributes at most 1.935 of a `beta_cr` near 2, so it matters.

## Addendum — the inelastic channels, and why they cannot be a rescue

The elastic-only relaxation rate was the cheapest open item above and the only one pointing the
helpful way: molecular and alkali excitation channels are far faster than elastic transfer, so
including them lowers `T_e - T_g` and raises `beta_cr`. Implemented as `inelastic_loss` and
`INELASTIC_CHANNELS`, with `electron_energy_balance` now closing on
`Q_joule = Q_elastic + Q_inelastic`.

**They are real, they are helpful, and they do not change the verdict.**

| | elevation `T_e - T_g` | `beta_cr` | `beta` | `beta/beta_cr` |
| --- | ---: | ---: | ---: | ---: |
| cold-leg exit, elastic only | 1055 K | 1.84 | 4.64 | 2.52 |
| cold-leg exit, **+ inelastic** | **814 K** | **2.07** | 4.63 | **2.23** |
| crossing (7569 K), elastic only | 466 K | 3.79 | 3.93 | 1.037 |
| crossing (7569 K), **+ inelastic** | **456 K** | **3.85** | 3.93 | **1.020** |

The unstable stretch is unchanged: `A/A* ~ 1.97` to the exit, 1.70 ms of a 2.71 ms transit. The
crossing station's e-folding lengthens 54 -> 99 µs, still ~17 e-foldings in its own dwell alone,
with 3-4 µs downstream.

**Channel ranking, and it is set by composition rather than by cross-sections:**

| channel | at the crossing (7569 K) | at the exit (4596 K) |
| --- | ---: | ---: |
| K 4s-4p resonance (1.61 eV) | **2.26% of Joule** | **18.46%** |
| H2O stretch `nu1/nu3` (0.453 eV) | 0.00% | 4.82% |
| H2O bend `nu2` (0.198 eV) | 0.00% | 2.01% |
| H2O rotation (~0.005 eV) | 0.00% | 0.29% |

Water's channels have **almost no target**: the plume is 99.7-99.9998% dissociated across the
band (Q-M), so `n_H2O` is 1.9e-6 of the heavies at the crossing and 3.3e-3 at the exit. Rotation
is doubly suppressed — the net transfer scales as `(dE/k T_e)^2` near equilibrium, and 0.005 eV
against `k T_e ~ 0.5 eV` costs four orders. **The alkali resonance is the only inelastic channel
that matters**, which is the standard result for seeded plasmas.

**The structural reason no inelastic channel could have rescued this.** Detailed balance requires
the net transfer to vanish at `T_e == T_gas`, so the bracket
`[1 - exp((dE/k)(1/T_e - 1/T_gas))]` expands to `(dE/k) dT / T^2` for small elevations — **linear
in the elevation, exactly as the elastic channel is.** The inelastic-to-elastic ratio is therefore
a property of the state, not a quantity that grows as the plasma approaches equilibrium. Adding
channels rescales the loss by a fixed factor; it cannot change the *scaling* that sets how far
`T_e` sits above `T_gas`. Anything that removes the drive has to act on `Q_joule` or on how cold
the plume gets, not on the loss channels.

**Consequence for the design: unchanged.** The `k` lever stands as the fix — keep the cold leg's
exit above ~11 000 K, where `beta_cr ~ 200` against `beta ~ 1.2`.

**What is now the binding uncertainty.** The crossing station sits at `beta/beta_cr = 1.02`, and
`Q_en` and `ln Lambda` carry factor-2 uncertainty that `sigma`'s own docstring admits. So *where*
the leg crosses into instability is not decidable by this model; that it crosses, and that the
exit is unstable by a factor >2, is. The alkali cross-section (`3e-19 m^2`, a step model for an
energy-dependent curve) is exposed on `INELASTIC_CHANNELS` for the same reason `Q_EN` is exposed
on `sigma` — but at 18% of the budget it would take an implausible factor of ~10 to matter.

## Addendum 2 — does a *cooling* plasma decouple the way a *driven* one does?

The criterion is validated on generators, where a cold gas is heated by an imposed current. This
plume is the reverse: a hot equilibrium plasma expanding and cooling. Two mechanisms could
separate `T_e` from `T_gas` here that have no generator analogue, and neither was in the balance.

**1. The expansion outrunning electron-heavy equilibration.** It does not, by six orders:

| `A/A*` | `T_gas` | `tau_eps` (e⁻↔heavy) | `tau_exp` (cooling) | ratio |
| ---: | ---: | ---: | ---: | ---: |
| 1.00 | 12 948 K | 0.0004 µs | 1 170 µs | **3.3e6** |
| 1.97 | 7 569 | 0.0037 | 2 317 | **6.2e5** |
| 4.00 | 4 596 | 0.0115 | 5 380 | **4.7e5** |

Ions, protons and neutrals never separate at all — comparable masses equilibrate in ~one
collision. Only the electron population can run, and the expansion is far too slow to make it.

**2. Three-body recombination depositing the store into electrons.** `X+ + e + e -> X + e` puts
the released energy in the *third body*, which is an electron — a genuine recombining-plasma
effect. Upper bound, with **all** of it into electrons and none into heavies:

| `A/A*` | `T_gas` | `dT` from recombination | `dT` from Joule | ratio |
| ---: | ---: | ---: | ---: | ---: |
| 1.00 | 12 948 K | **0.1 K** | 3.2 K | 0.040 |
| 1.97 | 7 569 | **0.3 K** | 456 K | 0.001 |
| 4.00 | 4 596 | **0.1 K** | 814 K | 0.000 |

Fractions of a kelvin. The release is real (4.7e9 W/m^3 at the throat) but spread over
milliseconds against a nanosecond redistribution time. **Both mechanisms raise `T_e`, so both
would have *lowered* `beta_cr` — they run against the design, and they are simply too small.**

So the field-gradient current remains the only thing separating `T_e` from `T_gas` anywhere in
this plume, which is what `electron_energy_balance` already models.

**The extrapolation is not load-bearing where it is largest.** The near-equilibrium end
(`dT/T_e ~ 1e-5`) is far outside the regime the criterion was validated in, so the hot-leg
verdicts were re-checked at `s = 0` — the maximally non-equilibrium limit, i.e. squarely inside
that regime:

| leg | `dT/T_e` | `beta` | floor at `s = 0` (`1.935 f + 0.065`) | needs `s`? |
| --- | ---: | ---: | ---: | --- |
| 75 km/s, all stations | 0.00000–0.00003 | 0.24–0.38 | 1.88 | **no** |
| 56.53, all stations | 0.00002–0.00067 | 0.39–1.20 | 1.58–1.70 | **no** |
| 45.58, throat | 0.00025 | 0.92 | 1.30 | **no** |
| 45.58, unstable band | 0.006–0.150 | 2.06–4.63 | 0.63–1.02 | yes |

**Every stable verdict survives with `s` deleted entirely.** The `beta_cr` values of 200–74 000 are
the formula's asymptote and are not defended as measurements; they are not needed. Conversely the
stations that *do* rely on `s` are the cold leg's unstable band at `dT/T_e` = 0.06–0.15, within a
factor of a few of the generator regime. **Where the transfer is most questionable it does not
matter, and where it matters we are in regime.**

Tightest of the hot-leg numbers: the 56.53 exit clears the `s = 0` floor by only 1.32×
(`beta` 1.20 vs 1.584), against factor-2 uncertainty in `beta` from `Q_en`/`ln Lambda`. Robust in
practice — `s` would have to fall from 203 to ~1, i.e. the elevation wrong by ~200× — but it is
the one hot-leg margin that is not enormous.

## Addendum 3 — the cold-leg verdict rests on one unsolved quantity, and it is not robust

Raised by Seth 2026-08-24, questioning where the `T_e > T_g` gap comes from at all in a plasma
that is neither driven nor strongly two-temperature. The question is well founded and this
section **downgrades the cold-leg finding from a result to a model-dependent possibility.**

Addenda 1 and 2 closed every candidate source of the temperature gap except one: the
field-gradient current. So `T_e - T_g` is set entirely by `Q_joule = j^2/sigma` with
`j = B/(mu0 L_eff)`, and `L_eff` is **estimated by a scaling, not solved**:

| station | `L_eff` choice | `L` [m] | `T_e - T_g` | `beta_cr` | `beta` | verdict |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| crossing (7569 K) | skin, `sqrt(t/(mu0 sigma))` — shipped | 1.72 | 456 K | 3.86 | 3.93 | **UNSTABLE** |
| crossing | full flux-tube radius | 4.21 | **77 K** | **18.13** | 3.93 | **stable** |
| exit (4596 K) | skin — shipped | 2.12 | 699 K | 2.23 | 4.63 | **UNSTABLE** |
| exit | full flux-tube radius | 6.00 | **79 K** | **10.58** | 4.74 | **stable** |

**A factor 3 in thickness flips every unstable station**, because it is 9x in `Q_joule`, hence 9x
in `T_e - T_g`, hence 9x in `s`. The "What would falsify this" entry above originally called
`L_eff` a weak lever on the verdict. That was wrong and is corrected.

**A physical argument points the stable way.** Plasma beta at the exit is `2 mu0 p / B^2` ~ **0.016**
— magnetic pressure exceeds gas pressure ~60x. A plume that cannot appreciably distort the field
does not concentrate its current into a thin resistive skin; the diamagnetic current is small and
distributed over the pressure-gradient scale, i.e. the flux-tube radius. That is the *stable*
column above. The skin estimate is the conservative one and is retained as shipped, but it should
not be reported as the answer.

`Q_en` is a secondary lever: doubling it flips the crossing to stable, while the exit survives ×3.

**What survives unqualified:**

1. **The hot legs are stable under every variant tested** — they clear even the `s = 0` floor
   (Addendum 2), so no choice of `L_eff` or `Q_en` reaches them.
2. **`BETA_CRIT = 2` is wrong for this plume** and the seed-only gain was a real bug. Both
   corrections stand independently of any of this.
3. **The cold leg is the only place the question can bite**, which is a genuine narrowing over
   Q-F(b)'s "live concern in this design".

**What it would take to settle it.** The current distribution in a low-plasma-beta magnetic nozzle
is an MHD problem this repository does not solve and `expansion.py` does not carry. Until it is
solved, the honest statement is: *the cold leg is the only exposed candidate, and whether it is
actually exposed is undecided.* The `k` lever remains the fix **if** it is exposed, and costs
2-7% Isp; sizing it should wait on this.
