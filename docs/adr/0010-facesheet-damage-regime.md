# Facesheet atomic/displacement damage is out of scope: per-atom energy is sub-threshold

A reviewer will ask whether hypervelocity gas atoms crystal-damage the SiC facesheet. They do not,
anywhere in the envelope, for a quantitative reason — so displacement/penetration/sputtering damage
is explicitly out of scope, and the survivability frontier is thermal/structural-limited only.

**Per-atom energy is eV-scale, at or below the SiC displacement threshold.** The damaging species is
the heaviest abundant gas atom — oxygen (16 amu). Its kinetic energy `½mv²`:

| v (km/s) | O atom KE |
|---|---|
| 3.2 | 0.85 eV |
| 11 | 10 eV |
| 16 | 21 eV |

(Hydrogen is ~16× lower; 1.3 eV even at 16 km/s.) The SiC displacement threshold is `E_d ≈ 20–35 eV`
(C ~20, Si ~35). Elastic transfer O→C/Si is near-maximal (`γ = 4m₁m₂/(m₁+m₂)² ≈ 0.93–0.98`), so even
at 16 km/s the transferred energy (~20 eV) only marginally reaches the carbon threshold and never the
silicon one; below ~13 km/s it is sub-threshold outright. Gas atoms cannot crystal-damage the SiC.

**Contrast with Project Orion.** Orion's tungsten jets (184 amu) at ~100 km/s carry ~9.5 keV per atom
— ~10³× the PuffSat per-atom energy, far above any displacement/sputter threshold. Orion was genuinely
in a sputtering/displacement regime; PuffSat is not. Both lower atomic weight (H, O vs W) and lower
velocity contribute, the `v²` term dominating, so `m·v²` lands ~3 orders of magnitude below Orion on
the benign side of the threshold.

**A few-micron ablator stops the atoms with enormous margin.** A ~20 eV atom has a stopping range
< 1 nm; even a 21 keV atom (O at the out-of-scope near-Sun ~500 km/s) ranges only ~tens of nm. An atom
would need ~MeV (oxygen at ~thousands of km/s) to penetrate microns of ablator. Penetration is never
the limiting factor in or beyond the envelope; the binding constraint on the few-micron ablator is
**thermal burn-through per pulse** (ablation depth < layer thickness) — the ablation-per-pulse
requirement reported at the frontier — not atomic penetration.

**The real SiC threats are chemical and thermal, both already handled (§5).** Oxidation by dissociated
atomic O (→ the SiC choice and SiO₂ passivation) and transient thermal shock (→ the pressure-limit
table). The ablator's role is vapor/thermal shielding plus a sacrificial oxidation surface, not
atom-stopping.

## Consequence

Crystal displacement, sputtering, and penetration damage are not modeled. The survivability frontier
is limited by peak facesheet pressure (structural) and ablation burn-through (thermal/logistical) only.

## Amendment (2026-06, Rung S): the peak-pressure frontier, operationalized

The pressure-limited frontier this ADR scopes is now computed (design §7). Peak facesheet pressure is
the stagnation pressure **`≈ 2.0·ρv²`**, verified physical from the 1D kernel's `peak_wall_force` (the
coefficient is 2.0 at both 11 and 16 km/s, i.e. `peak ≈ 2·ρv²`, the cold-cloud ram pressure
recompressed at the wall). The `Σ = m/(π r_foot²) = ρL` contract (ADR-0003) maps each cloud shape to a
density, so the frontier is `peak(L/D, r_foot/R, v)` against the **`P_limit = 400 MPa` baseline** (the
§5 conservative floor), swept to 700/900 MPa at 16 km/s. The damage regime is unchanged — thermal /
structural, atomic damage out of scope — and the binding limit is the compressive facesheet pressure
(the reflected-tensile spall stays sub-dominant, ADR-0011). Ablation-per-pulse is reported as a
back-propagated MEMS-replenishment *requirement*, not gated here. **Result:** the `f`-maximizing
short-disk / tight-footprint corner fails by a wide margin (~2.3 GPa at 16 km/s); the best *survivable*
`f` is ≈ 0.80 (dip) / ≈ 0.78 (16 km/s baseline) / ≈ 0.84 (16 km/s relaxed 900 MPa).

## Amendment (2026-06): the closed-form `f`-margin map over plate radius `R` and pulse mass `m`

Peak facesheet pressure is **intensive** — the local stagnation stress `≈ 2·ρv²`, set entirely by the
gas at the wall — so it is blind to every plate *geometry* change (facesheet thickness, plate width as
empty acreage, total force). Only two things move it: the **gas density** `ρ` (cloud shape, already
pushed to its `eta_capture` limit at the frontier above) and the **material allowable** `P_limit`
(low-leverage, design §5/line 105). This amendment records the one remaining handle — the two **scale
knobs**, plate radius `R` and pulse mass `m` — as a closed-form margin map, since a reviewer naturally
asks "why not a wider plate / a smaller pulse?"

**The scaling is analytic, so this needs no kernel reruns.** Via the Σ contract
`ρ = m / (2π·(L/D)·(r_foot/R)³·R³)` (ADR-0003), the peak scales as **`ρ ∝ m/R³`**, while
`eta_capture(r_foot/R)` is a pure *ratio* and so is scale-invariant. A wider plate or a smaller pulse
therefore only **relaxes the pressure ceiling** by a `headroom = (R/R₀)³·(m₀/m)` factor — admitting a
denser, higher-`eta` shape that failed at the baseline, and buying `f` back. The margin map is a pure
rescaling of the Rung-S frontier (`analysis.py --axis margin`, `margin_map`; `make analysis-margin`):

| 16 km/s, 400 MPa baseline | headroom | best survivable `f` |
|---|---|---|
| `R = 5 m, m = 25 kg` (pinned baseline) | 1.0× | **0.784** (just under the 0.8 gate) |
| e.g. `R = 5 m, m = 15 kg` | 1.7× | **0.814** (clears the gate) |
| e.g. `R = 6.5 m, m = 25 kg`  or  `R = 5.5 m, m = 15 kg` | 2.2× | **0.835** (comfortable margin) |
| `R = 7 m, m = 15 kg` (grid corner) | 4.6× | 0.835 (plateaus) |

(The dip behaves the same: `0.804 → 0.829`.) Two honest bounds: the gain is **real but modest and
stepped** (limited to the discrete sampled cloud shapes, `~+0.05` in `f`), and it **plateaus at
≈0.835** — the absolute `f`-max corner (≈0.86) needs `~5.75×` headroom, beyond even the grid's `4.6×`,
because that corner is the densest case the frontier already forecloses.

**These are not levers the study pulls — they are pinned by *external* budgets, and this is only the
`f`-side of a system trade.** `R` is set by the vehicle dry-mass budget (a wider pusher plate is the
mass the architecture exists to minimize; `R = 5 m` fixed, design §2/§7 — and *shrinking* it is
catastrophic, `1/R³` the wrong way). `m` is set by the per-pulse thrust × pulse-rate (a smaller pulse
lowers the peak but lowers delivered impulse `(1+e_eff)·m·v` in lockstep, so it needs proportionally
*more* pulses — a throughput knob, not an efficiency one; `f` itself is scale-invariant in `m`). The
cost-side (plate-mass(`R`), pulse-count(`m`)) lives outside this per-collision study, so the map
yields trade *curves*, not an optimum. The `m`-leg additionally leans on the deferred Σ-resolved
`e_eff(ρ)` lookup (ADR-0013) to be quantitatively sharp. **Decision:** recorded as a de-risking margin
map (the passing `f ≈ 0.8` does not *need* it); not pursued as a sweep rung, and not a change to the
deliverable `f(v)`. Plate *bending* remains a separate, cloud-shape-invariant thrust-structure problem
(impulse-driven, ms-scale) handled by the impedance-Ti + truss backing under the rigid-wall condition
(ADR-0011), and buys back no `f`.
