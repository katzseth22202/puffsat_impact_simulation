# Restitution is measured as a momentum ratio from wall impulse

We define the effective restitution as a **momentum** ratio read directly off the
wall-pressure history:

```
e_eff = J_wall / p_in − 1
```

where `J_wall = ∫ P_wall(t)·A dt` is the time-integrated axial force the gas exerts on the
plate during one pulse, and `p_in = m_pulse · v` is the incident axial momentum in the plate
frame. Equivalently `J_wall = p_in · (1 + e_eff)`: the wall receives the incoming momentum
plus the rebound, so `e_eff` ranges 0 (dead stick) to 1 (elastic). `J_wall` is integrated
until the wall force decays to `10⁻³` of its peak, and `f`'s sensitivity to that cutoff is
reported alongside the result so the long re-expansion tail cannot silently move the number.

This is the only definition under which the smoke-test limits hold — `f → 1` for a lossless
elastic wall (in the precise sense of the Amendment below) and `f → 0.5` for a dead inelastic
stick — and under which `e_eff` plugs directly into the paper's
`f = eta_capture · (1 + e_eff) / 2`. It also sidesteps the need to define a
"rebound velocity," which is ill-posed: the gas leaves the wall as a rarefaction fan with every
fluid element at a different speed.

## Considered Options

- **Energy-based restitution** (rebound KE / incident KE). Rejected: it breaks the `f → 0.5`
  inelastic limit and double-counts thermal energy that stays in the (optically thick) gas
  rather than reappearing as useful axial momentum.
- **Rebound-velocity definitions** (some characteristic outflow speed ÷ `v`). Rejected:
  ill-posed for a rarefaction fan — there is no single rebound speed to divide by.

## Amendment (2026-06-17): the `f → 1` elastic limit, precisely

The phrase "`f → 1` for a lossless elastic wall" above is an idealized limit, not the value a
finite gas slug actually reaches. Design §9 is the governing statement: a real re-expanding
slug has `e_eff < 1` **even with zero losses**, because the rebound is a rarefaction fan with a
velocity spread (by Cauchy–Schwarz the coherent rebound momentum is below the incident). That
finite-amplitude **lossless gas-dynamic `e_eff`** — some value `< 1`, set by γ — is the *true
bounce ceiling*, not 1.

So the rung-A elastic smoke test is **two checks of different character** (neither asserts a
finite-amplitude `e_eff == 1`):

1. **Conservation bookkeeping** (runs every build, value-independent): with a perfectly
   reflecting wall, the impulse measured from the wall-pressure history, `J_wall = ∫P_wall dt`,
   must equal the slug's momentum change read off the mesh, `Δp_slug`, to round-off. This
   catches gross hydro and impulse-accounting bugs regardless of the `e_eff` value.
2. **Low-Mach limit**: `e_eff → 1` as `M → 0`. In the acoustic limit a pulse reflects off a
   rigid wall with full momentum reversal and negligible rarefaction spread, so this is the
   genuine physical limit that makes `f → 1` more than an accounting identity.

The finite-amplitude lossless `e_eff` (the ceiling) is then a **measured, refinement-checked
output**, reported — not an asserted constant. The `f → 0.5` dead-stick limit is unaffected: a
sticking wall returns no momentum, so `J_wall = p_in`, `e_eff = 0`, `f = 0.5` exactly.

**Configuration (rung A).** Both limits run on one harness: a **finite cold gas slug
re-expanding into vacuum** — uniform `ρ₀`, velocity `v` toward a wall at `x=0`, a free-surface
(`p=0`) outer boundary — not a sustained semi-infinite column (a column's wall force never
decays, so the `10⁻³` cutoff is meaningless for it). The slug is parameterized by incident Mach
`M = v/c_s₀`: `M → 0` is the elastic low-Mach limit, finite `M` maps the lossless ceiling. A
**tail guard** stops the run at the `10⁻³` cutoff rather than fully evacuating the domain (the
hot, low-density re-expansion cells would otherwise dominate the CFL timestep while adding
almost nothing to `J_wall`); the residual in-flight momentum at stop-time is reported as a
closure check (`J_wall + p_residual ≈ p_in·(1+e_eff)`). Kernel scheme is ADR-0022.
