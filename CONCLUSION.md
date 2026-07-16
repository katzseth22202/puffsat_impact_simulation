# Conclusion: `f ≈ 0.8` is realistic

**Headline.** Within this study's models, a per-collision fudge factor of **`f ≈ 0.8` is
physically realistic across the full 3.2–16 km/s impact envelope** — the best survivable `f`
lands at `≈ 0.77–0.82` over the envelope, clearing the gate at the top (16 km/s, with the
ablating wall) and sitting just under it (within the numerics band) at the ~11 km/s worst-case
dip. This is a
**single-code result**; an independent hydrocode cross-check is the one open validation gate
(see *Status* below).

> Per [ADR-0009](docs/adr/0009-useful-f-gate.md) (2026-06 paper-draft amendment), a **preliminary**
> `f` *is* quotable externally (i.e. in the white paper) ahead of that cross-check, provided it is
> framed as preliminary: a single-code, conservative-floor result with the FLASH cross-check named
> as the outstanding independent validation. The cross-check upgrades the number from *preliminary*
> to *validated*; it no longer blocks quotation.

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
| Survivability (peak facesheet pressure) | — | reflected-shock stagnation `≈ 1.2·ρv²` ([ADR-0010](docs/adr/0010-facesheet-damage-regime.md)/[0011](docs/adr/0011-sic-ti-shock-reflection.md); 2026-07 correction — the earlier `≈ 2·ρv²` was an artificial-viscosity artifact) |

**The worst case is interior, not an endpoint.** The conservative (equilibrium-EOS) `e_eff`
does not bottom out at either velocity limit — it dips in the *transition* near **~11 km/s** to
`e_eff ≈ 0.57`, below both the 3.2 km/s (`≈ 0.74`) and 16 km/s (`≈ 0.63`) endpoints. The worst
case for the whole study is therefore this dip, and that is the case the `0.8` claim is tested
against (Rung T; [ADR-0012](docs/adr/0012-transitional-anchor.md)).

**Best survivable `f`** (conservative floor: rigid wall, worst-case `e_eff`, concave focusing
penalty applied; 2026-07 numbers — grid-converged 2D `eta_capture` at 112×80 with physical Mach
anchors, and the physical stagnation coefficient `c_stag ≈ 1.2` in place of the earlier
artificial-viscosity artifact `≈ 2.0`):

- **At the ~11 km/s dip:** survivable `f ≈ 0.768` on a moderate-footprint shallow-concave disk
  (`≈ 0.773` with the small ablating-wall recovery) — **just under 0.8**, within the study's
  ±0.03 numerics band of the gate; the plate-radius/pulse-mass margin map reaches `≈ 0.780`.
- **At 16 km/s:** survivable `f ≈ 0.798` at the *reference* plate (R = 5.0 m, m = 25 kg) under
  the strict 400 MPa baseline — on the line — and the ablating wall lifts it to **`0.802–0.821`**
  across its full Q\*/τ bracket, so the top of the envelope clears the gate. With the corrected
  (lower) peak pressure, survivability is barely binding at 16 km/s: relaxing to 900 MPa adds
  nothing at the reference plate, and the margin map plateaus at `≈ 0.810`.

**`0.8` is a realistic family, not a knife-edge either way.** The two 2026-07 audit corrections
moved the two anchors in opposite directions (dip `0.804 → 0.777`, 16 km/s `0.784 → 0.805`), and
a third correction (2026-07-10: two 2D-kernel defects found while chasing an M = 40 blow-up —
[ADR-0023](docs/adr/0023-2d-axisymmetric-euler-kernel-numerics.md) correction — had left the
concave `eta_capture` corners ~1 % high) lowered both by ~0.01 (dip `0.777 → 0.768`, 16 km/s
`0.805 → 0.798`). All are third-decimal shifts inside the ±0.03 numerics band. The honest
statement is unchanged: the best survivable `f` lands at `≈ 0.77–0.82` across the envelope,
centered on the paper's `0.8`, with the worst case (the dip) just under the line and the top of
the envelope reaching above it.

**Freeze-timing bracket (frozen recombination) — the largest quantified physics uncertainty**
([ADR-0026](docs/adr/0026-frozen-recombination-bracket.md)). The equilibrium EOS *returns* the
banked dissociation/ionization energy during the rebound — the one assumption stacked
optimistically, and one a FLASH cross-check with the same equilibrium EOS would silently share.
Bounding the freeze timing both ways: freezing the composition at turnaround (maximal,
"freeze-after-the-plate" bound) drops the dip `e_eff` 0.570 → 0.398 (best survivable `f`
≈ 0.68 dip / ≈ 0.74 at 16 km/s — **below the gate**); chemistry-free pure H₂O
("freeze-before-the-plate") raises it to 0.661 (`f` ≈ 0.81–0.82). The equilibrium curve stays
the headline because three-body recombination at the probed turnaround densities
(`n ~ 10²⁰–10²¹ cm⁻³`) is ~10²–10³× faster than the ~µs rebound, so the gas tracks equilibrium
through the dense phase where the momentum is exchanged — the sudden-freeze end is a bound,
not an expectation. But any external quotation of `f` should carry this bracket.

**Pulse-shape sensitivity — slight shape changes cost only slight impulse** (design §13,
[ADR-0028](docs/adr/0028-shape-sensitivity-fixed-design-protocol.md);
`data/results/shape_sensitivity.csv`/`.png`). Raw `f(shape)` at the fixed baseline design
(`d/D = 0.1`, `L/D = 0.3`, `r_foot/R = 0.5`), perturbing only the arriving pulse over an
**assumed** shape box (±20% footprint and aspect, edge taper to 30% of `r_foot`, radial
divergence to `α = 0.1`): the normalized sensitivity `S = (Δf/f)/(Δx/x)` maxes at **0.26**
(footprint, concave; ≈ 0.1 aspect, ≤ 0.02 taper/divergence per full box) — far under the
"≲ a few" bound, so a 10% delivery-shape error costs ≲ 2.6% of `f`, linearly correctable by the
pushed vehicle's guidance. No cliff survives grid refinement, the taper Σ-profile bound is small
(Δf ≤ 0.003), and a three-point frozen-chemistry spot-check at the dip shows the same gentle
slope (Δe_eff across the box: equilibrium +0.010, sudden-freeze +0.018) — the smoothness is not
an equilibrium-chemistry artifact. Off-axis modes (offset, tilt, drift) are bounded analytically
as linear with a `0.5·R` rim-clip margin (§13). The separate **survivability margin check**:
comfortable at the dip (+28% under 400 MPa over the whole box), but at 16 km/s the nominal
already sits near the line, so footprint *concentration* beyond ~−5% exceeds the 400 MPa
baseline (peak `∝ 1/r_foot³`); the whole box clears the relaxed 900 MPa limit. The shape box is
an **assumption standing in for undetermined delivery dispersion** (the deferred cloud-schedule
study owns the real numbers) — every quote of these `S` values must carry that caveat.

**Loss budget at 16 km/s** (`data/results/frontier.csv`): of the momentum *not* returned,
≈ 76–80 % is lost to the radiative wall flux and the remaining ≈ 20–24 % escapes to space;
conductive loss is negligible on the pulse timescale (Rung C). The ablating wall recovers a real
but small slice of this (`e_eff` 0.63 → up to ~0.68 at 16 km/s; injection-dominated, not
shielding-dominated), but the ~11 km/s dip is an EOS energy sink and is **not** radiatively
fillable ([ADR-0014](docs/adr/0014-ablating-wall.md)).

## Status & open items

**One validation gate (upgrades a preliminary number to validated; no longer blocks quotation):**

- **Independent hydrocode cross-check (FLASH).** Everything above is from this study's own
  solvers. Per ADR-0009 (2026-06 amendment), a preliminary `f` is quotable in the white paper
  ahead of this, framed as single-code; the cross-check is what promotes it from *preliminary* to
  *validated*. This is the single remaining gate. **Note:** a cross-check run with an equilibrium
  EOS does *not* discharge the freeze-timing caveat below — it shares the assumption
  ([ADR-0026](docs/adr/0026-frozen-recombination-bracket.md)).

**Three non-blocking refinements (sharpen a known floor; cannot move the verdict):**

- **Finite-rate recombination chemistry.** The freeze-timing bracket above is two-sided and its
  physical end is argued from the recombination timescale (~10²–10³× faster than the rebound at
  turnaround densities); a partial-equilibrium/freeze-out-density rung would interpolate the
  bracket rather than move the headline (ADR-0026). Named here because the pessimistic bound —
  unlike every other refinement — does cross the gate. The **69 km/s Jupiter special scenario**
  now carries its own, wider bracket (`data/results/frontier_frozen_jupiter.csv`): the multi-charge
  ~140–170 kK turnaround needed the extended-grid frozen pipeline with the full 8-stage O ionization
  ladder, giving equilibrium headline `f = 0.78` → **sudden-freeze `f = 0.69`** (pure-H₂O 0.75).
  There the recombination-tracks-equilibrium defence is *weaker* — the turnaround is dilute
  (`n ~ 10¹⁹ cm⁻³`, ~10–100× below the dip) with an 8-stage ladder to unwind — so this caveat is
  more load-bearing for that scenario ([ADR-0026](docs/adr/0026-frozen-recombination-bracket.md)
  amendment).

- **Real opacity table.** Pulled 2026-07-11 for the **69 km/s Jupiter scenario only** (TOPS/OPLIB
  gray means, ADR-0007 amendment) — there it mattered (~2000× the interim Kramers at stagnation,
  `e_eff` 0.42 → 0.65, best survivable `f` 0.69 → 0.78). For **this envelope** the interim table
  stands: the dip is EOS-sink-limited, not opacity-limited, and `e_eff` was shown insensitive to
  opacity (≤ 1.6 % over a 100× sweep in κ, Rung B) because the slab here is already `τ ≫ 1`. The
  molecular/low-v half of the seam (HITEMP/ExoMol) remains unpulled.
- **High-velocity plasma gas conductivity (`k_gas`).** Deferred at the top of the envelope; a
  refinement to the conductive channel, which is already negligible on the pulse timescale.
- **Coupled-radiation stability window (found 2026-07-16, shape study).** The 1D coupled
  radiation operator has a **resolution-onset radiative collapse**: past a critical grid
  refinement that coarsens with `ρv²`, the thin wall cell's radiative drain zeroes its energy,
  the slab loses pressure support, and the run dies mid-infall with unphysical `e_eff` (the
  radiative sibling of the documented conductive over-drain). At 16 km/s the onset is 1200 cells
  at the production ρ = 0.64 and ~200–300 cells by ρ ≈ 1.2–1.7. All quoted numbers sit on the
  verified **stable plateau** (flat to < 0.002 in `e_eff` across a 4× cell range, and the
  EOS-only curve is smooth everywhere, so the collapse is numerics, not physics); the shape
  study's few affected samples use a two-resolution validity protocol (300/150 cells, flagged in
  `shape_sensitivity.csv`). The named fix is an implicit matter–radiation exchange at the wall —
  a kernel robustness refinement, not a result change.

## Pointers

- Full physics rationale and method: [`puffsat_impact_sim_design.md`](puffsat_impact_sim_design.md)
- The `0.8` gate and its survivable-frontier amendment: [ADR-0009](docs/adr/0009-useful-f-gate.md)
- Dual-curve deliverable (conservative floor vs. best estimate): [ADR-0013](docs/adr/0013-dual-curve-deliverable.md)
- Numbers and figures: `data/results/frontier*.csv` and the `*.png` plots
