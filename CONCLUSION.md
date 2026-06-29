# Conclusion: `f ≈ 0.8` is realistic

**Headline.** Within this study's models, a per-collision fudge factor of **`f ≈ 0.8` is
physically realistic across the full 3.2–16 km/s impact envelope** — achievable at the worst
case and reachable at the top of the envelope under ordinary design choices. This is a
**single-code result**; an independent hydrocode cross-check is the one open validation gate
(see *Status* below).

> Per [ADR-0009](docs/adr/0009-useful-f-gate.md), no `f` is quoted *externally* (i.e. in the
> paper) until that cross-check lands. The numbers here are internal modeling results recorded
> for the project.

## What was asked

Paper §3.2 assumes a constant fudge factor `f = 0.8`
(`f = eta_capture · (1 + e_eff) / 2`, the axial momentum a gas pulse delivers to the pusher
plate as a fraction of full-capture-perfect-bounce; see [`CONTEXT.md`](CONTEXT.md)). The goal of
this repo was narrow and specific: **show that `0.8` is realistic** — not to derive a single
optimal design, and not to claim `0.8` is guaranteed with margin.

## The result

All three momentum-loss mechanisms are measured independently:

| Mechanism | Symbol | Track |
|---|---|---|
| Geometric capture | `eta_capture` | 2D axisymmetric Euler |
| Effective restitution (radiative/conductive/condensation loss) | `e_eff` | 1D Lagrangian rad-hydro |
| Survivability (peak facesheet pressure) | — | stagnation `≈ 2·ρv²` ([ADR-0010](docs/adr/0010-facesheet-damage-regime.md)/[0011](docs/adr/0011-sic-ti-shock-reflection.md)) |

**The worst case is interior, not an endpoint.** The conservative (equilibrium-EOS) `e_eff`
does not bottom out at either velocity limit — it dips in the *transition* near **~11 km/s** to
`e_eff ≈ 0.57`, below both the 3.2 km/s (`≈ 0.74`) and 16 km/s (`≈ 0.63`) endpoints. The worst
case for the whole study is therefore this dip, and that is the case the `0.8` claim is tested
against (Rung T; [ADR-0012](docs/adr/0012-transitional-anchor.md)).

**Best survivable `f`** (conservative floor: rigid wall, worst-case `e_eff`, concave focusing
penalty applied):

- **At the ~11 km/s dip:** survivable `f ≈ 0.804` on a moderate-footprint shallow-concave disk,
  peaking under the 400 MPa baseline pressure limit. **Clears 0.8.**
- **At 16 km/s:** survivable `f ≈ 0.784` at the *reference* plate (R = 5.0 m, m = 25) under the
  strict 400 MPa baseline — just under 0.8 — rising to **`≈ 0.835`** once the plate is given
  modest radius/mass headroom *at that same conservative limit*, and again under the relaxed
  900 MPa SiC+Ti limit. Two independent levers, either of which clears 0.8.

**`0.8` is a forgiving design family, not a knife-edge.** The margin sweep
(`data/results/frontier_margin.csv`) shows the 16 km/s baseline `f` climbing from `0.784`
(R = 5.0, m = 25) through `0.835` (R ≈ 6.0, m = 20) as plate radius rises relative to pulse mass.
A wider or heavier plate buys margin directly, so reaching `0.8` does not depend on hitting one
exact geometry.

**Loss budget at 16 km/s** (`data/results/frontier.csv`): of the momentum *not* returned,
≈ 76–80 % is lost to the radiative wall flux and the remaining ≈ 20–24 % escapes to space;
conductive loss is negligible on the pulse timescale (Rung C). The ablating wall recovers a real
but small slice of this (`e_eff` 0.63 → up to ~0.68 at 16 km/s; injection-dominated, not
shielding-dominated), but the ~11 km/s dip is an EOS energy sink and is **not** radiatively
fillable ([ADR-0014](docs/adr/0014-ablating-wall.md)).

## Status & open items

**One validation gate (blocks the externally-quoted number):**

- **Independent hydrocode cross-check (FLASH).** Everything above is from this study's own
  solvers. ADR-0009 reserves any external `f` until a second, independent code reproduces the
  result. This is the single remaining gate.

**Two non-blocking refinements (sharpen a known floor; cannot move the verdict):**

- **Real opacity table.** Never pulled (data source firewall-blocked). The dip is EOS-sink-limited,
  not opacity-limited, and `e_eff` was shown insensitive to opacity (≤ 1.6 % over a 100× sweep in
  κ, Rung B). A real table refines the radiative channel; it does not move `f` across 0.8.
- **High-velocity plasma gas conductivity (`k_gas`).** Deferred at the top of the envelope; a
  refinement to the conductive channel, which is already negligible on the pulse timescale.

## Pointers

- Full physics rationale and method: [`puffsat_impact_sim_design.md`](puffsat_impact_sim_design.md)
- The `0.8` gate and its survivable-frontier amendment: [ADR-0009](docs/adr/0009-useful-f-gate.md)
- Dual-curve deliverable (conservative floor vs. best estimate): [ADR-0013](docs/adr/0013-dual-curve-deliverable.md)
- Numbers and figures: `data/results/frontier*.csv` and the `*.png` plots
