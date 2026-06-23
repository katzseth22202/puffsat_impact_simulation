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

## Considered Options

- **Raw-throughput `eta_capture`** (2D axial momentum ÷ incident). Rejected: double-counts the
  gas-dynamic re-expansion already present in `e_eff`, understating `f`.
- **Monolithic 2D rad-hydro swept across the frontier.** Rejected: hundreds of coupled runs,
  impractical and unnecessary — `tau >> 1` localizes the radiation to a 1D wall problem.
- **Per-area-normalized 2D/1D ratio.** Rejected in favor of the confined-2D denominator: a separate
  1D-equivalent normalization would reintroduce scheme-error mismatch the same-kernel ratio cancels.
