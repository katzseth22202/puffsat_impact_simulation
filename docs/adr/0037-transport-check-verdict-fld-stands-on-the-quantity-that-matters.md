# The transport check finds FLD's escape channel wrong and the deliverable unaffected

ADR-0012 gated a coupled M1 solver behind a cheaper question: run a post-hoc Sₙ /
short-characteristics diagnostic first, and escalate **only if** the escape flux differs by >10%.

The diagnostic was built, it fails that gate decisively — 24 of 32 states above 10%, worst 116% —
and **the escalation is nonetheless declined**, because neither cause is a closure error and M1
would fix neither. What earns FLD its reprieve is not accuracy on the escape channel, which it does
not have, but a *measured* insensitivity of `e_eff` to that channel.

**Reported radiation-model error bar: `Δf ≤ 0.0033`.** Below the opacity bracket's 0.0055
(ADR-0036), so it is not the dominant uncertainty.

## What the comparison actually compares

FLD reports the escape to space through its Marshak surface as `F = (c/2)·E_surface` — the
*half-range-isotropic* answer, presuming the emergent intensity is the same in every outward
direction. Transport computes the first moment of the actual angular distribution,
`F = 2π ∫₀¹ I⁺(μ) μ dμ`. The two agree when `I⁺` really is isotropic and part company when it is
not — limb-brightened at small `τ`, forward-peaked across a transparent gap.

The audit is **one-way**: it observes states the FLD run produced without feeding back. A test
asserts the audited bounce is bit-identical to the unaudited one, so what it yields is a bias
estimate on a loss channel, never a corrected `e_eff`.

## Cause 1 — a free-surface defect, real and verified

FLD's Marshak boundary converts the space-facing cell's *average* `E` into an emergent flux, which
is only meaningful while that cell resolves the Planck mean free path. Verified in a controlled
isothermal slab against the exact `σT⁴`:

| `Δτ_P` (surface cell) | 1.0 | 0.5 | 0.25 | 0.125 | 0.031 |
|---|---|---|---|---|---|
| FLD flux / `σT⁴` | 1.43 | 1.27 | 1.17 | 1.11 | 1.07 |

Converging to 1.0, as it must. Production runs carry an escape-weighted `Δτ_P` from 0.2 to **48**.
The signature is exact: every dense row reports `FLD/blackbody = 1.99–2.01`, i.e. `(c/2)aT⁴` against
the true `(c/4)aT⁴`, because emission pins the cell at equilibrium instead of letting it drain to
`aT⁴/2`. **A mesh/BC defect, not a closure error** — and one no M1 addresses, since M1 applies a
surface condition to the same cell.

## Cause 2 — the gray two-mean model, which nothing short of multigroup fixes

FLD here emits and absorbs with `χ_P` and diffuses with `χ_R` (ADR-0006). A gray transport solve has
one extinction coefficient, so **no single-mean Sₙ run is "FLD minus the closure"** — an exact
closure-only comparison does not exist against a two-mean model. Running both brackets the ambiguity,
and the two tallies disagree with *each other* by up to **164%**, more than either disagrees with
FLD. Per ADR-0036 that spread is spectral structure, not model error.

## Why the deliverable survives: measured, not argued

Across a **16× mesh refinement** (300 → 4800 cells):

- the escape channel climbs **+8.5% to +13.8% and is still rising** at 4800 cells;
- `e_eff` settles to better than **0.7%**.

The two are nearly decoupled, and the reason is physical: the escape comes from the outer
re-expansion tail — gas that has already stopped pushing on the wall — so misplacing its energy
barely touches the momentum integral that defines `e_eff`.

| state | `e_eff` @300 | @4800 | Richardson limit | Δ from 300 | Δf |
|---|---|---|---|---|---|
| 69 km/s ρ=0.02 | 0.449977 | 0.444638 | 0.444579 | +0.0054 | 0.0026 |
| 45 km/s ρ=0.04 | 0.628710 | 0.626881 | 0.626760 | +0.0020 | 0.0010 |
| 28 km/s ρ=0.01 | 0.525947 | 0.519861 | 0.519205 | **+0.0067** | **0.0033** |
| 69 km/s ρ=0.30 | 0.703994 | 0.702917 | 0.701433 | +0.0026 | 0.0013 |

**A correction on the record.** The first estimate of the impact was `Δf = 0.058`, a first-order
energy-book calculation assuming the misplaced escape energy would otherwise have become rebound
momentum. It does not. That estimate was high by ~25×, and the refinement study is what replaced it
with a measurement. Production stays at 300 cells: 0.0033 is below the opacity bracket, 600 cells
would cost 4×, and 300 keeps continuity with the published 16–28 km/s results.

## Consequences

- **No coupled M1.** It changes the closure, which is not what is wrong. The `hydro1d` refactor seam
  that was contingent on it stays closed.
- **No free-surface fix either.** It would repair a genuine defect that is not moving `f`, could not
  address cause 2, and — since the natural implementation *is* the in-cell formal solution the audit
  computes — would make the Sₙ escape check circular thereafter.
- **`Δf ≤ 0.0033` is reported as the radiation-model error bar.**

## Recorded, unresolved

An isolated **+46% `peak_wall_pressure` artifact** at 28 km/s, ρ = 0.01, 2400 cells (1.382e7 Pa
against ~9.4e6 at every other resolution). The other three states converge smoothly. It feeds the
facesheet survivability gate (ADR-0010/0027), not the restitution, so this verdict is unaffected —
but it should be explained before that gate is quoted at the dilute end.
