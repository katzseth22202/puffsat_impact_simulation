# PuffSat Pusher-Plate Restitution Simulation — Design Document

**Status:** Draft for review
**Date:** 2026-06-16
**Author context:** Seth Katz (paper author). Design distilled from a grilling session on how to simulate the energy-loss "fudge factor" `f` in *Aim Is All You Need: A Speculative White Paper on PuffSat Pulsed Propulsion*.
**Intended home:** a separate implementation repository (not this LaTeX paper repo). Companion to `puffsat_control_sim_design.md`; the two share only the vehicle, not the physics.
**Backs:** paper §3.2 ("Radiative Differences From Project Orion"), which currently states "detailed computer analysis is needed for full confirmation." This sim is that analysis.

---

## 1. Purpose

Compute one quantity: the **per-collision coefficient of restitution** of the PuffSat gas against the pusher plate, expressed as the paper's fudge factor

```
f = eta_capture * (1 + e_eff) / 2
```

as a function of relative impact speed across the mission envelope. `eta_capture` is the fraction of axial momentum that lands on the plate and rebounds usefully (geometry). `e_eff` is the effective restitution that survives radiative, conductive, and condensation losses (thermophysics).

The paper's mass-ratio law `m_r/m_p = 2f / ln((v_p - v_ri)/(v_p - v_rf))` (Appendix, eq. PuffSat_ratio) currently assumes a single constant `f = 0.8`. This sim replaces that assumption with a defensible `f(v)`.

**This is a per-collision study, not a performance study.** The result is `f(v)` and its loss budget. How those pulses are then smoothed into vehicle acceleration (the buffer invariant) is explicitly out of scope (§9).

---

## 2. Reference case

| Quantity | Value | Note |
|---|---|---|
| Pulse mass | 25 kg | gas delivered per PuffSat |
| Plate radius | <= 5 m | sets footprint, eta_capture, column density |
| Plate mass | 3-4 t | the **plate**, not the vehicle; only used to confirm rigid-wall (§5) |
| Impact speed | 3.2 to 16 km/s | swept; see §3 |
| Baseline gas | water (H2O) | covers icy/off-world PuffSats; see §4 |
| Chemical energy | **0** | conservative; see §4 |

Velocity envelope reasoning: during LEO insertion the target accelerates, so the relative impact speed falls from ~11 km/s toward ~3.2 km/s (the target reaches ~7.8 km/s orbital). Decelerating shots can close at up to ~16 km/s. Near-Sun work (hundreds of km/s) uses solid projectiles and a different chamber; it is a separate study (§9).

---

## 3. The physics that shaped the design

A scaling pass fixed the regime before any tool was chosen. The numbers below are order-of-magnitude and drove the decisions; they are not the deliverable.

**Bulk motion is the energy source, not the chemistry.** At 11 km/s the specific kinetic energy is 60.5 MJ/kg; total 1.5 GJ for a 25 kg pulse, about 13x the chemical energy of the same mass of TNT. So composition matters for the gas *species* (EOS, opacity, condensation), not for the energy budget. We do not model detonation; we model a cold gas arriving at `v`.

**The collision crosses three regimes across the velocity envelope:**

| v (km/s) | KE (MJ/kg) | KE/chem | T (naive cap) | Regime |
|---|---|---|---|---|
| 3.2 | 5.1 | 1.3 | < 1700 K | cool neutral gas, **optically thin** |
| 6.0 | 18 | 4.5 | ~6000 K | dissociating / partially ionized |
| 11 | 60 | 15 | ~20 kK | strong plasma, optically thick |
| 16 | 128 | 32 | ~43 kK | strong plasma, optically thick |

**At high speed, opacity is "won" with margin.** At plausible cloud densities (0.16-0.64 kg/m^3) the stagnated plasma is deeply optically thick (tau ~ 1e2 to 1e5 even on a pessimistic free-free estimate). Project Orion's opacity argument holds. Radiation does not free-stream into the plate; it diffuses slowly. So `e_eff` is set mainly by gas dynamics and geometry, not by radiative loss, and "ensuring a good bounce" is the dominant lever precisely because opacity is favorable. Radiation re-enters only as the survivability constraint.

**At low speed the opacity argument fails, and a different one takes over.** At 3.2 km/s the stagnated gas is ~1700 K, neutral, and optically thin. But it is so cold it radiates almost nothing (flux scales ~v^8, so ~5e-5 of the 11 km/s level). The low-v bounce is governed by classical gas-dynamic re-expansion minus condensation/recombination loss. With chemical energy zeroed, this is the **worst case in the whole study: water at 3.2 km/s.** If it clears a useful `f`, every other case is easier.

**The survivability-defining case is the opposite end: 16 km/s.** Ram pressure scales as v^2 (2.1x the 11 km/s value) and radiative flux as ~v^8 (~20x). A disk-like shot at 16 km/s is ~393 MPa stagnation pressure. So the worst case for `f` and the worst case for survivability sit at opposite ends of the velocity range.

---

## 4. Gas model (settled)

- **Baseline gas: water.** This covers the icy and LOX off-world PuffSats directly, and with chemical energy zeroed it also covers terrestrial explosive PuffSats, so one `f(v)` curve serves all three. Water is also the conservative choice for restitution: it re-condenses readily (latent heat 2.26 MJ/kg), and condensed mass that deposits at the wall sticks instead of bouncing, which is the dominant low-v loss now that chemistry is gone.
- **Confirmatory spot-checks:** oxygen (LOX) and an ANFO-product mix (N2/H2O/CO2/CO), run at the optimum rather than swept.
- **Chemical/explosion energy set to zero in the bounce.** The chemistry still gasifies the cloud (it sets the cloud-at-impact state) but contributes nothing to the rebound. Consequence: `f <= 1` strictly. The pessimism is a ~40% energy haircut at 3.2 km/s (where KE ~= chemical energy) and negligible above ~8 km/s (6% at 11, 3% at 16). This is the price of one curve valid for explosive and icy PuffSats alike.
- **EOS:** equilibrium with dissociation and ionization (Saha / CEA-style). An ideal-gamma EOS is wrong here because it skips the ionization enthalpy that buffers the temperature (the paper's §3.2 "higher specific heat capacity" point) and so over-predicts T and the ~T^4 flux.
- **Opacity:** real tabulated opacities for the actual species, **not Kramers.** Kramers free-free assumes full ionization and is orders of magnitude wrong at the 10-20 kK partial-ionization conditions. The survivability flux that defines the frontier rides entirely on getting this right.

---

## 5. Plate model (settled)

**Construction: a three-layer stack.**

1. **Sacrificial ablator** (silicone or ionic liquid), renewed each pulse by space-qualified inkjet-style MEMS (paper §4.x): vapor-shields and self-heals the surface.
2. **SiC-based hot face** (CVD or monolithic SiC; C/SiC composite as an alternative), permanent: oxidation resistance plus high-temperature strength. This is the layer that meets the gas.
3. **Ductile metallic backing (titanium)**: carries the impulse into the shock absorber and gives fracture toughness so the brittle ceramic cannot spall catastrophically. It must stay below ~700 K, so it sits behind the ceramic and ablator and never sees the hot face directly.

**Why SiC, not carbon-carbon.** Every PuffSat gas we care about is oxidizing: LOX directly, and water and CO2 both dissociate to atomic oxygen at impact temperatures. Bare carbon-carbon burns (C + O -> CO/CO2) and is consumed. SiC passivates, growing a protective SiO2 layer, and survives. So C-C is out as the hot face for our chemistry. This sharpens what the paper already gestures at (§4.x: oxygen attacking PICA phenolic, adding B4C/HfC; and the LOX "non-combustible lining" note). The real engineering risk in the stack is the SiC-Ti bond: their thermal expansion differs by ~2x (SiC ~4 ppm/K, Ti ~9), so the join needs a graded interlayer or compliant braze. That is a Phase-2 structural and thermal-shock task, not a gas-dynamics one.

**Rigid wall is a good approximation.** A 3-4 t plate recoils only `2*m_gas/M_plate` ~= 1.4% of the impact speed per pulse (velocity-independent, ~230 m/s at 16 km/s). That 1.4% is the entire error in the rigid-wall assumption. We start rigid (Phase 1) and add the ablating, transpiring wall (Phase 2) for the real frontier. If the 1% recoil correction is ever wanted, the 1D solver takes the plate as a heavy piston for free.

**Plate shape: model flat and shallow concave (depth/diameter ~0.1-0.15).** Flat is the conservative `eta_capture` floor (hemispherical rebound, large sideways loss). Shallow concave re-collimates the rebound toward the axis without a deep dish's penalties. A deep dish was considered for the highest `f` ceiling and dropped to an unrun upper-bound note: its focal hot spot is geometric (not removable by reshaping the incoming cloud) and it couples radiation into the otherwise radiation-free 2D geometry track.

**Transient pressure-limit table** (engineering placeholders; impulsive load with thermal shock):

| Construction | Transient limit | Note |
|---|---|---|
| ablator / light substrate | 50-150 MPa | floor, non-viable |
| carbon-carbon face | 150-300 MPa | **oxidizes -> out** |
| C/SiC composite face | 200-400 MPa | tough; fibers vulnerable if matrix cracks |
| **SiC (CVD/monolithic) + Ti** | **400-700 MPa** | baseline; best oxidation; Ti carries load |
| UHTC (HfC/ZrC/HfB2-SiC) + Ti | 600-900 MPa | highest T and oxidation; heavy, costly |

At 16 km/s a stretched pulse is ~100 MPa and a disk-like one ~393 MPa, so the design is feasible across C/SiC and up. We are above the feasibility threshold, not fighting to clear it.

**Pressure-limit sensitivity (prediction to confirm, not a result).** Optimizing the pressure limit is predicted to be **low leverage** on `e_eff`, closer to "+5%, not worth it" than "doubles it," for three reasons: (1) `e_eff`'s radiative trapping is set by *column* density (mass / footprint area) and opacity, not by volumetric density or pressure, so a stronger plate only lets you run denser/sharper at a fixed footprint, which barely moves the column density; (2) past the tau >> 1 knee the residual loss is condensation and conduction, which pressure does not fix; (3) the local pressure is the only structural limit in scope, and even it sits above the feasibility threshold. **Decision rule:** baseline SiC+Ti for oxidation and feasibility margin, sweep `e_eff(P_limit)` and `f(P_limit)` at the 16 km/s anchor to find the knee, and push to UHTC only if the sweep shows a steep gain past SiC+Ti (doubted). The genuine `e_eff` levers are the dark-oil opacity boost (§7) and the footprint/curvature, not facesheet strength.

---

## 6. The fudge factor and its loss budget

```
f = eta_capture * (1 + e_eff) / 2
```

- `eta_capture` (geometry): axial-momentum fraction that lands and rebounds usefully. Set by plate radius, curvature, and cloud footprint/divergence. Computed by the 2D Euler track (§7).
- `e_eff` (thermophysics): `1 - (radiative + conductive + condensation losses)`. Computed by the 1D rad-hydro track (§7).

**Explicit output: the loss decomposition.** `(1 - f)` is reported split into four channels, per velocity anchor:

1. radiative-to-plate (the **optical** fraction, called out as a headline number),
2. conductive-to-plate,
3. condensation / recombination (mass that sticks),
4. sideways escape / non-capture (geometric).

Expectation: the optical fraction is largest at 16 km/s (but trapped by tau >> 1) and small at 3.2 km/s (cold gas barely radiates; the low-v loss is condensation- and conduction-dominated). Reporting it directly is what lets us decide whether the opacity-boost lever (§7) is worth pulling.

`f = 0.8` from the paper, for reference, implies e.g. `eta_capture ~ 0.95` and `e_eff ~ 0.68`. The sim's job is to put real numbers on both, across `v`.

---

## 7. Velocity sweep, design variables, and the opacity-seed study

**Two physics packages, anchored at ~3 velocities** (3.2 / ~8 / 16 km/s):

- **Low-v package:** cool water vapor, optically thin, condensation/recombination loss, no chemistry. The condensation model is the key piece, since it is the dominant low-v sink.
- **High-v package:** water plasma, equilibrium EOS, real opacity tables, flux-limited radiation diffusion.
- A transitional anchor between them.

The opacity argument (Orion / §3.2) justifies only the high-v end. The low-v end rests on a separate argument: a cold gas radiates little, so little energy is lost radiatively regardless of opacity, and the bounce is gas-dynamic minus condensation.

**Plate fixed, cloud scheduled.** The plate is one built object (curvature fixed, optimized once). The cloud shape can vary per shot, so it is a schedule `shape(v)`: stretch hard at 16 km/s to survive, run sharper/denser at 3.2 km/s to beat condensation.

**Design variables (sweep primitive = the cloud-at-impact state, not deployment hardware):**

- `rho_impact` (impact density). There is no independent "impact pressure" knob: a coasting cloud is internally cold, so the plate feels the ram pressure `rho*v^2` (~2.4x at stagnation). The optimal density is sandwiched between an opacity floor (need tau >~ a few; already met at ~0.16 kg/m^3) and a survivability ceiling, and within that band you want the lowest density that keeps tau >> 1 (the longest survivable pulse). Past the tau >> 1 knee, more density barely helps `e_eff`.
- `L/D` (the disk <-> cylinder axis).
- footprint coverage `r_foot / R`.
- radial divergence.
- plate curvature (flat, shallow).

Standoff, initial size, and detonation energy are backed out afterward; they are not modeled.

**Survivability frontier (local limits only):** for each shape, find the densest/sharpest cloud that keeps peak facesheet pressure under the SiC+Ti limit and ablation under budget; record `f` there. The total-force/buffer-invariant leg is dropped (out of scope). Disk-like shots fail on pressure, over-concentrated cylinders fail on local burn-through and waste plate area, so the optimum is intermediate: an elongated, mildly diverging cloud fanned to roughly plate diameter on a shallow concave plate.

**Dark-oil opacity seed (study point at the cool anchor).** Seeding ~0.5% mass of a dark ablative oil into the cloud (distinct from the plate-surface dark ablative already in the paper) pyrolyzes to soot. A rough calc shows 0.5% soot raises the cool 3.2 km/s gas from tau ~0.05 to ~11-45 (a ~375x jump), flipping it from optically thin to thick. The study point asks: does seeding recover optical loss and improve `f` (and survivability) at the cool end? It is **gated** on the loss decomposition showing the optical fraction is a meaningful part of the loss there. Likely the bigger payoff is at the transitional/high-v end and on survivability, not on `f` at 3.2 km/s. A 0.5% contamination is judged acceptable against the biocompatibility goal.

---

## 8. Tooling and implementation

### Architecture: factored, because tau >> 1 makes the radiation nearly 1D

A monolithic 2D rad-hydro swept across the frontier is impractical (hundreds of runs) and unnecessary: tau >> 1 makes radiation local-diffusive rather than a global transport problem, and the radiative-vs-hydro timescale competition is a 1D problem at the wall.

| Track | Tool | Produces |
|---|---|---|
| restitution + wall load | **1D Lagrangian rad-hydro** (equilibrium EOS, real opacity tables, flux-limited diffusion, rigid wall then ablating wall) | `e_eff`, peak flux, peak pressure vs (rho, pulse shape) |
| geometric capture | **2D axisymmetric Euler**, radiation-free | `eta_capture` vs plate curvature and cloud shape |
| confirmation | **FLASH** at the chosen optimum (deferred) | independent cross-check |

The 2D track runs radiation-free because the gas is optically thick (energy stays in the gas) and the capture question is adiabatic geometry. Exception: a deep dish would couple radiation into the focus, which is one reason it was dropped.

### Implementation language: Rust core, Python at the table and analysis boundary

The Rust lean is correct, but the decisive reasons are kernel speed and a trustworthy from-scratch solver, **not** the GIL. The sweep is embarrassingly parallel (independent runs), and even in Python that parallelizes cleanly via multiprocessing, which dodges the GIL entirely; the GIL only bites in-process threading, which a parameter sweep does not need.

- **Rust for the hot path:** the 1D rad-hydro kernel, the 2D Euler kernel, and the sweep driver (rayon for data-parallel runs across all cores). This is where the hundreds of runs live. Rust gives machine speed, fearless parallelism, and a solver whose correctness you can defend.
- **Python for the cold path:** generate the EOS and opacity tables once as offline preprocessing (Cantera / CoolProp / CEA for equilibrium chemistry; opacity from tables or atomic-physics interfaces), and do all plotting, frontier extraction, and loss-decomposition analysis on the solver's output. This ecosystem is the real reason to keep Python in the loop, and it never touches the hot loop.
- **Boundary = a file format** (HDF5 or Parquet/CSV): Python writes tables, Rust reads them and writes results, Python plots results. No FFI or PyO3 binding needed unless later desired.

Two worries answered directly:

- *Do we need Python optimization libraries?* No. The "optimization" is a grid sweep plus knee-finding (the problem is constraint-dominated and `f` saturates near 1, so there is no smooth interior maximum to chase). Even if a real optimizer is wanted later, it is a cheap outer loop wrapping the expensive kernel, so its language is almost irrelevant. Rust has `argmin` if needed.
- *Are there Python hydrocodes we should use?* None fit this problem. Custom EOS, opacity, ablating-wall, and condensation models mean you write the kernel yourself regardless, so Rust is the better choice for it. (PyClaw/Clawpack and astro teaching codes exist but would all need to be rebuilt for this physics.)

**Optional de-risking path:** prototype the 1D kernel in NumPy first to pass the Sod / Marshak / elastic-limit smoke tests (§10), then port the validated kernel to Rust for the production sweep. This avoids debugging the physics and the borrow checker at the same time. Going Rust-first with `ndarray` is fine if preferred.

---

## 9. Validation and verification

**Deferred (frontier first), with a live smoke test.** The decision was to map the frontier first and do full V&V after, but the two analytic momentum limits are kept as a continuous smoke test during development: `f -> 1` for a lossless elastic wall and `f -> 0.5` for an inelastic stick. They cost nothing and catch gross hydro bugs. **No `f` value is quoted externally until the validation below lands.**

- **Verification (solver solves its equations):** Sod, Noh, Sedov for the Euler track; a Marshak wave for the flux-limited diffusion; the `f -> 1` / `f -> 0.5` momentum limits.
- **Validation (physics is right), keystone = Orion:** reproduce Project Orion's published impulse-per-pulse and ablation-per-pulse from their stated inputs (Balcomb 1970, already cited in the paper). The entire `f` argument descends from Orion's opacity finding, so reproducing their numbers and then extrapolating the same code to PuffSat conditions (lower-Z, lower-T, ~1/5 the velocity) is the §3.2 argument made quantitative. No other single check buys as much credibility.
- **Cross-code:** FLASH at the chosen optimum.
- **Provenance:** opacity tables from a citable source (TOPS / OPLIB), not hand-rolled.

---

## 10. Build order (rungs)

- **A. 1D ideal-gas smoke test.** Lagrangian hydro, effective-gamma EOS. Pass Sod and the `f -> 1` / `f -> 0.5` limits. (NumPy prototype acceptable; this is where the physics gets debugged.)
- **B. 1D high-v package.** Add equilibrium EOS + real opacity + flux-limited diffusion + rigid wall. Produce `e_eff(rho)` and wall flux/pressure at 16 km/s.
- **C. 1D low-v package.** Cool-gas EOS, optically thin, condensation/recombination. Produce the worst case: water at 3.2 km/s.
- **D. 2D Euler geometry.** `eta_capture` for flat and shallow-concave plates vs cloud shape.
- **E. Ablating wall + levers.** Phase-2 transpiring wall; the dark-oil opacity-seed study; the pressure-limit sensitivity sweep.
- **F. Validation.** FLASH cross-check; Orion impulse/ablation reproduction. Only now is `f(v)` quotable.

---

## 11. Out of scope

- **Total-performance / vehicle-acceleration analysis** (the buffer invariant `m*s ~= M a T^2/4` and its inputs: vehicle mass, target acceleration, pulse cadence, absorber stroke). `f(v)` is a per-collision property reusable regardless of how pulses are smoothed into motion.
- **Spatially-varying `f` across the plate** (off-center hits load the RCS). The first pass takes `eta_capture` as a footprint-averaged quantity.
- **Near-Sun regime** (hundreds of km/s), which the paper handles with solid projectiles and a separate chamber.
- **Composition sweep beyond water + the O2/ANFO spot-checks.**
- **Multi-pulse plate thermal accumulation between shots** (depends on cadence, which is a performance input).
