# PuffSat Impact Simulation

A physics simulation to compute the per-collision coefficient of restitution of PuffSat gas against a pusher plate — the fudge factor `f` in the paper *Aim Is All You Need: A Speculative White Paper on PuffSat Pulsed Propulsion*.

## What This Computes

The core output is `f(v)`, the fudge factor as a function of impact speed across the mission envelope (3.2–16 km/s):

```
f = eta_capture * (1 + e_eff) / 2
```

- `eta_capture` — fraction of axial momentum that lands and rebounds usefully (geometry)
- `e_eff` — effective restitution after radiative, conductive, and condensation losses (thermophysics)

This backs paper §3.2, which currently assumes a constant `f = 0.8`. The sim replaces that with a defensible `f(v)` curve and a full loss budget.

## Architecture

Two physics tracks feed a combined result:

| Track | Tool | Output |
|---|---|---|
| Restitution + wall load | 1D Lagrangian rad-hydro (Rust) | `e_eff`, peak flux, peak pressure |
| Geometric capture | 2D axisymmetric Euler, radiation-free (Rust) | `eta_capture` vs plate shape |
| Confirmation | FLASH (deferred) | independent cross-check |

Python handles EOS/opacity table generation (Cantera/CoolProp/CEA) and all analysis/plotting. The boundary between Rust and Python is a file format (HDF5 or Parquet/CSV).

## Build Order

- **A.** 1D ideal-gas smoke test — Sod shock tube, elastic/inelastic momentum limits
- **B.** 1D high-v package — equilibrium EOS, real opacity tables, flux-limited diffusion at 16 km/s
- **C.** 1D low-v package — cool-gas, optically thin, condensation loss; worst case is water at 3.2 km/s
- **D.** 2D Euler geometry — `eta_capture` for flat and shallow-concave plates
- **E.** Ablating wall and levers — transpiring wall, dark-oil opacity seed study
- **F.** Validation — FLASH cross-check and Project Orion impulse/ablation reproduction

## Reference

Full physics rationale, gas model, plate construction, design variables, and validation plan: [`puffsat_impact_sim_design.md`](puffsat_impact_sim_design.md)
