# `f` is factored as `eta_capture · (1 + e_eff)/2`, with `eta_capture` a lossless 2D/1D ratio

`f` is computed as a product of two independently-run solvers — a 1D rad-hydro kernel for
`e_eff` and a radiation-free 2D Euler kernel for `eta_capture` — rather than one monolithic 2D
radiation-hydrodynamics simulation swept across the frontier. This is justified by `tau >> 1`
making radiation local-diffusive (effectively 1D) at the wall, so a full 2D rad-hydro sweep
(hundreds of coupled runs) is unnecessary.

**The subtlety that drives the definition:** the gas-dynamic re-expansion appears in *both*
factors. Even a perfectly lossless gas has `e_eff < 1`, because the rebound is a rarefaction fan
with a spread of velocities — by Cauchy-Schwarz, `p_rebound = ∫ρu dz < (∫ρ dz · ∫ρu² dz)^½ = p_in`,
with equality only for a uniform outflow that never occurs. If `eta_capture` were defined as the
raw 2D axial throughput `J_wall^{2D}/(2·p_in)`, it would carry that same re-expansion penalty, and
multiplying by `(1 + e_eff)/2` would double-count it, understating `f`.

**Decision:** define `eta_capture` as the **lossless 2D/1D wall-impulse ratio**

```
eta_capture ≡ J_wall^{2D, no-loss} / J_wall^{1D, no-loss}
```

both runs adiabatic, same EOS, same cloud, same `v`. The common re-expansion divides out, leaving
only the geometric effect — sideways loss on a flat plate, re-collimation recovery on a concave
one. This makes the 1D case `eta_capture = 1` by construction, which is the correct ceiling: 1D
*is* perfect collimation.

**Contract between the tracks:** the 1D `e_eff` run inherits the 2D operating-point footprint's
column density `Σ = m_engaged / (π·r_foot²)`, because column density sets `tau` and hence radiative
trapping. The tracks are separable in their *outputs* but coupled in this one *input*.

**The combined-physics cross-code at the optimum (rung F) validates this factorization** — that
geometry and thermophysics do not cross-couple appreciably — not merely the individual kernels.

## Amendment (2026-06): the `eta_capture` extraction, operationalized (flat plate)

The 2D kernel (`crates/euler2d`, ADR-0023) landed and the flat-plate `eta_capture` is now extracted.
How the lossless 2D/1D ratio above is realized in code:

- **`J_wall^{2D}` = `Σ_{plate cells} p(r,t)·(n̂·ẑ)·dA`, trapezoid-integrated in time** to the same
  `10⁻³`-of-peak-force decay cutoff as ADR-0001, with the "don't stop in a trough before a focused
  secondary peak" guard carried forward (flat plate does not focus, but the concave follow-on will).
  Gas past the plate edge (`r > r_plate`) is never counted — the plate cells are `r ≤ r_plate`, the
  rest of the `z=0` boundary is transmissive (§7).
- **The `1D, no-loss` denominator is a *confined-2D* run of the same kernel** — same cloud mass / `v`
  / EOS, but reflecting outer `r`-wall and the cloud filling the radius (the perfectly-collimated
  plane-wave limit). Forming the ratio from two runs of the *same* solver makes the scheme error (not
  just the EOS error) common-mode, so it cancels — `eta_capture = J_wall^{free}/J_wall^{confined}`.
  This avoids any per-area normalization and realizes "1D ⇒ `eta_capture = 1`" by construction.
- **Cross-check (independent kernel):** the confined-2D `1+e_eff` reproduces the 1D `hydro1d` kernel's
  `run_bounce` to **~5 %** (Eulerian-Godunov vs Lagrangian-AV) — the strongest single check that the
  2D bounce physics is right, since the two solvers share no code.
- **Result (flat, M ≈ 5, γ = 1.4):** `eta_capture` rises **0.81 → 0.92** as the footprint widens over
  `r_foot/L = 0.5 → 2.0` (wider disks re-collimate toward the 1D ceiling 1; narrow slugs splat with
  more radial relief). Flat is the conservative hemispherical-rebound floor (ADR-0021); shallow
  concave can only raise it. The parametric sweep and the concave plate are the follow-on rung.

## Amendment (2026-06): the concave `eta_capture` and the `f`-reconciliation + `Σ` contract, operationalized

The shallow-concave plate landed (ghost-cell IBM, ADR-0023 amended) and the `f = eta_capture·(1+e_eff)/2`
product is now formed in code (`analysis.py --axis geometry`):

- **Both flat and concave run through the same IBM boundary**, so the same-kernel cancellation above
  extends to the curved case — the curvature gain carries no grid-alignment artifact (a D4c consistency
  gate ties the IBM flat wall to the verified grid-aligned flat `eta_capture`, rel < 0.10). Shallow
  concave can push `eta_capture > 1` (over-collimation past the flat plane-wave limit, ADR-0021); the
  "1D ⇒ `eta_capture = 1`" identity above fixes the *flat* plane-wave denominator and is not a ceiling.
- **The `Σ = m_engaged/(π r_foot²)` contract reduces, for a uniform cylinder, to `Σ/ρ = L = 2·(L/D)·r_foot`** —
  the footprint `r_foot` cancels (mass scales with the same `π r_foot²`), so `r_foot/R` is purely the
  `eta_capture` lever and `L/D` alone sets the 1D column density that fixed `e_eff`. In the sweep's
  normalized units (`ρ = 1`, `r_foot = 1`) this is `sigma_over_rho = 2·(L/D)`.
- **First `f` bracket (the deliverable of this rung):** pair each geometry case with the two 1D `e_eff`
  anchors — the transitional worst case `e_eff = 0.57` (≈ 11 km/s, the conservative floor, ADR-0012)
  and `e_eff = 0.63` (16 km/s). Flat floor `f ≈ 0.696` at the dip, shallow concave **`f ≈ 0.83` (> the
  0.8 useful gate)**; `f ≈ 0.86` at 16 km/s. The full `Σ`-resolved `e_eff(ρ(r_foot))` lookup that feeds
  the dual-curve `f(v)` deliverable (ADR-0013) is the noted refinement, a later rung.

  *Superseded numbers (2026-07):* the bracket above was computed on the original 56×40 grid at
  Mach anchors 5/10, which the 2026-07 audit showed was not grid-converged for the deep-dish/tight-
  footprint corner. On the converged 112×80 grid with physical Mach anchors 10/20 the concave
  maxima are `f ≈ 0.792` at the dip and `f ≈ 0.822` at 16 km/s (before survivability). The
  factorization decision itself is unchanged.

## Considered Options

- **Raw-throughput `eta_capture`** (2D axial momentum ÷ incident). Rejected: double-counts the
  gas-dynamic re-expansion already present in `e_eff`, understating `f`.
- **Monolithic 2D rad-hydro swept across the frontier.** Rejected: hundreds of coupled runs,
  impractical and unnecessary — `tau >> 1` localizes the radiation to a 1D wall problem.
- **Per-area-normalized 2D/1D ratio.** Rejected in favor of the confined-2D denominator: a separate
  1D-equivalent normalization would reintroduce scheme-error mismatch the same-kernel ratio cancels.
