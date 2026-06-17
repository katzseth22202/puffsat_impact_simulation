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
elastic wall and `f → 0.5` for a dead inelastic stick — and under which `e_eff` plugs directly
into the paper's `f = eta_capture · (1 + e_eff) / 2`. It also sidesteps the need to define a
"rebound velocity," which is ill-posed: the gas leaves the wall as a rarefaction fan with every
fluid element at a different speed.

## Considered Options

- **Energy-based restitution** (rebound KE / incident KE). Rejected: it breaks the `f → 0.5`
  inelastic limit and double-counts thermal energy that stays in the (optically thick) gas
  rather than reappearing as useful axial momentum.
- **Rebound-velocity definitions** (some characteristic outflow speed ÷ `v`). Rejected:
  ill-posed for a rarefaction fan — there is no single rebound speed to divide by.
