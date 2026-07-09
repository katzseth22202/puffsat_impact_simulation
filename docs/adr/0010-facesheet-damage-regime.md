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
the stagnation pressure ~~**`≈ 2.0·ρv²`**, verified physical from the 1D kernel's `peak_wall_force`~~
**`≈ 1.2·ρv²`** (see the 2026-07 correction below — the 2.0 coefficient was an artificial-viscosity
artifact, not physical). The `Σ = m/(π r_foot²) = ρL` contract (ADR-0003) maps each cloud shape to a
density, so the frontier is `peak(L/D, r_foot/R, v)` against the **`P_limit = 400 MPa` baseline** (the
§5 conservative floor), swept to 700/900 MPa at 16 km/s. The damage regime is unchanged — thermal /
structural, atomic damage out of scope — and the binding limit is the compressive facesheet pressure
(the reflected-tensile spall stays sub-dominant, ADR-0011). Ablation-per-pulse is reported as a
back-propagated MEMS-replenishment *requirement*, not gated here. **Result** (as corrected 2026-07):
the `f`-maximizing short-disk / tight-footprint corner still fails by a wide margin (~1.6 GPa at
16 km/s); the best *survivable* `f` is ≈ 0.78 (dip) / ≈ 0.81 (16 km/s baseline) / ≈ 0.81 (16 km/s
relaxed 900 MPa — with the corrected lower peak, pressure is barely binding at 16 km/s and the relaxed
limit buys almost nothing).

## Correction (2026-07): the stagnation coefficient is ≈ 1.2, not 2.0 — the 2.0 was an AV artifact

The Rung-S amendment above originally claimed `peak ≈ 2.0·ρv²`, "verified physical from the 1D
kernel's `peak_wall_force`," on the grounds that the coefficient was 2.0 at both 11 and 16 km/s
(velocity-independence read as physicality). That inference was wrong. `peak_wall_force` is the wall
cell's **total** pressure `p + q`, and its first-impact spike is dominated by the von Neumann–Richtmyer
artificial-viscosity term `q ≈ c_q·ρ·Δu² ≈ 2·ρv²` (production `c_q = 2.0`) — which is *also* ∝ ρv² and
so velocity-independent in coefficient. Varying `c_q` moves the "measured" coefficient in lockstep
(c_q = 2.0 → 2.02, c_q = 1.0 → ~1.2), the signature of a numerical artifact. The **physical** wall
pressure — the EOS `p(0, t)` alone — converges under grid refinement to the reflected-shock stagnation
value `(γ_eff+1)/2·ρv²` ≈ **1.20·ρv² at 11 km/s and 1.24·ρv² at 16 km/s** (water EOS, γ_eff ≈ 1.1–1.2).

**Fix:** the kernel now reports `peak_wall_pressure` (EOS `p` only, AV excluded) alongside
`peak_wall_force`, and the survivability/margin analyses back `c_stag` out of that (a stale pre-fix
JSONL fails loudly). The impulse/`e_eff` bookkeeping is untouched (the AV term belongs in the
*integrated* wall force; only the *peak* was misattributed). Direction: the old number was
**conservative** — survivability was ~1.7× too pessimistic. Folded together with the same audit's 2D
grid-convergence fix (the geometry sweep's 56×40 grid was not converged for the deep-dish/tight-
footprint corner; now 112×80 with physical Mach anchors 10/20), the best-survivable `f` moved
`0.804 → 0.777` at the dip and `0.784 → 0.805` at 16 km/s — opposite-sign third-decimal shifts inside
the study's ±0.03 numerics band. See CONCLUSION.md for the corrected headline numbers.

## Amendment (2026-06): the closed-form `f`-margin map over plate radius `R` and pulse mass `m`

Peak facesheet pressure is **intensive** — the local stagnation stress `≈ 1.2·ρv²` (2026-07
correction above), set entirely by the gas at the wall — so it is blind to every plate *geometry*
change (facesheet thickness, plate width as
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

(Numbers as corrected 2026-07 — physical `c_stag ≈ 1.2` and the converged 112×80 geometry grid:)

| 16 km/s, 400 MPa baseline | headroom | best survivable `f` |
|---|---|---|
| `R = 5 m, m = 25 kg` (pinned baseline) | 1.0× | **0.805** (already clears the gate) |
| mid-grid (e.g. `R = 6 m, m = 25 kg`) | ~1.7× | **0.806** (flat — pressure barely binds) |
| `R = 7 m, m = 15 kg` (grid corner) | 4.6× | **0.822** (plateaus) |

(The dip: `0.777 → 0.792` across the same grid.) Two honest bounds: with the corrected (lower) peak
pressure the baseline already survives its best shapes, so the headroom gain is **small and stepped**
(`~+0.02` in `f`, limited to the discrete sampled cloud shapes), and it **plateaus at ≈0.822** — the
absolute `f`-max corner is the densest case the frontier forecloses at any sampled headroom.

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
