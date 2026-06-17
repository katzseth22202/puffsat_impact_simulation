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

## Considered Options

- **Raw-throughput `eta_capture`** (2D axial momentum ÷ incident). Rejected: double-counts the
  gas-dynamic re-expansion already present in `e_eff`, understating `f`.
- **Monolithic 2D rad-hydro swept across the frontier.** Rejected: hundreds of coupled runs,
  impractical and unnecessary — `tau >> 1` localizes the radiation to a 1D wall problem.
