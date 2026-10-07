# Walled nozzle: the blast inside the near-term chamber, and what fixes it

Status: analysis note, 2026-10-07. Follows [`walled_nozzle_chamber_fatigue.md`](walled_nozzle_chamber_fatigue.md),
which found that a reflected blast spike, not the mean pressure, decides whether the chamber wall
survives. This note solves that spike. The chamber is the parent's near-term one: the 2.5 kg rod
at 75 km/s (7.03 GJ) in 20 m^3, methane at 7000 K (496 bar) and hydrogen at 5500 K (818 bar).

Reproduce:
- `make walled-nozzle-chamber-blast`: 1-D spherical runs, `euler2d::sphere`.
- `make walled-nozzle-vessel-blast`: 2-D axisymmetric chamber shapes, `euler2d::vessel`.
- `PYTHONPATH=python uv run python -m puffsat.walled_nozzle.volume_trade`.

## 1. Bottom line

1. **The paper's plug sends a 5-6 GPa hammer into the nozzle end on every pulse.** The 4.4 kg
   foam plug of `sec:material_chamber_plug` stops the rod but leaves 36% of its energy, 2.5 GJ,
   as the merged 6.9 kg moving on at 27 km/s. Almost nothing stands in its path: the gas on the
   axis is ~0.2 kg. In the 2-D chamber it strikes the convergent nose at **5.6-5.9 GPa** and the
   throat at **~2.8 GPa**. Steel spalls at ~2.9 GPa. The paper's own estimate, "about 0.2 GPa
   lasting microseconds", counts only the spherical part of the blast.
2. **A longer, lighter plug removes it.** Take ~10 kg of plug spread over 1.75-2.5 m of the
   axis: the paper's 4.4 kg plus ~5.6 kg, 8% of the methane charge, as foam or frozen methane.
   The rod's energy then becomes heat along the axis and no fast remnant survives. The nose
   drops to 0.4-0.8 GPa, and every wall region sits at 0.4-0.8 GPa at 1 cm cells, about half the
   point-blast level and roughly uniform. Mass matters less than length: 20 and 40 kg plugs do
   no better, and short heavy plugs load the port end.
3. **Shape the chamber like a rocket chamber, curved at both ends.** A 2:1 elliptical (domed)
   port head and an elliptical convergence give the lowest loads. A flat port face doubles the
   spike in its corner (2.8 against 1.4 GPa at 1 cm). A straight cone focuses the blast
   mid-length (1.6-1.8 GPa at 2 cm, against 0.7-0.8 on a cylinder).
4. **No carbon overwrap survives at 20 m^3.** The spike puts 0.36-0.6 of itself through the
   wall as tension. The resin holds 64 MPa across the fibres. Even the long-plug loads are about
   5-10x over that.
   - More steel ahead of the carbon does not shield it.
   - A spike longer than a layer's transit passes straight through, so the 5 mm aluminium layer
     does nothing at these pulse lengths.
   - Aluminium-alloy and magnesium backings spall earlier than steel.
   - The working wall at 20 m^3 is steel: multilayer (shrink-fitted) Cr-Mo for methane, about
     35-47 t, with ~2-3x margin on spall under the long plug.
5. **A bigger chamber trades Isp for a gentler spike.** The spike falls as 1/V at fixed pulse
   energy. Holding V/A* (so blowdown stays ~55 ms) and taking the nozzle's freeze debit, for
   methane:
   - 60-80 m^3 costs ~4-5% Isp and brings the long-plug spike to ~0.3-0.4 GPa. That is the
     range for a through-thickness-reinforced or dry Vectran wrap.
   - ~120-160 m^3 costs ~6-7% Isp and reaches ~0.16-0.21 GPa, the plain overwrap's threshold.
   - Hydrogen pays about twice the Isp for the same volume, because recombination freezes more
     at lower pressure.
6. **The liquid-share spray does not soften the wall shocks.** That is the paper's stated intent
   for it (`sec:carbon_overwrap`).
   - Swept-up mist is methane, and so is the gas around it. Both pay the same dissociation toll,
     with Γ_eff 1.13-1.27.
   - Drops of 0.5-1 mm shatter and match the gas in ~3 µs, while the front moves ~3 cm.
   - Inert extra mass at the wall raised the spike 10-20%.
   - The liquid is better spent as plug mass on the axis.
7. **Hydrogen at this pulse size is still unresolved.** Its centred-blast spike, ~2.95 GPa, is
   at steel's spall strength. A heavy axial plug either does not fit as frozen hydrogen or costs
   too much exhaust speed as polyethylene.

## 2. The solvers

- **`euler2d::sphere`** is 1-D spherical MUSCL-Hancock/HLLC, conservative in mass and energy,
  with a well-balanced pressure source. Its tests: gas at rest stays at rest; the closed sphere
  conserves mass and energy; the blast radius matches Sedov (ξ0 = 1.0328 at γ 1.4) within 1%.
- **`euler2d::vessel`** is a closed chamber wall `r_w(z)` (domed head, cylinder, cone-to-elliptical
  convergence, open throat) as a ghost-cell immersed boundary in the 2-D kernel. Its tests: gas
  at rest stays at rest in a curved chamber; a closed-cylinder blast holds its mass within 2%.
- **Gas.** The gas is γ-law at 1.2. Equilibrium methane has Γ_eff = 1 + p/(ρe) of 1.13-1.27 over
  the swept states (`eos_methane`), and 1.4 is the bracket. Each history is read against its
  own quasi-static pressure (QSP) and scaled to the solved real-gas 496 / 818 bar.
- **Deposition.** The rod's 7.03 GJ (rod + plug 6.9 kg) goes in as heat and, where stated,
  axial momentum, with energy and momentum exact on every grid.

## 3. Results

**Centred blast, 1-D sphere, grid series** (methane, γ 1.2, model units, QSP 70.7 MPa):

| cell | 2 mm | 1 mm | 0.5 mm | 0.25 mm | extrapolated |
|---|---|---|---|---|---|
| wall spike | 1.72 GPa | 2.02 | 2.23 | 2.36 | ~2.57 (~36x QSP) |

That is ~1.8 GPa for methane and ~2.95 GPa for hydrogen at the real QSP. The breathing overshoot
D is 2.0-2.3 here, against the paper's 1.7, because the reflections refocus at the centre.

**Chamber shapes, 2-D, 2 cm cells** (methane, γ 1.2, peak in GPa):

| contour (20 m^3) | port-face corner | body | convergence |
|---|---|---|---|
| cylinder + 45° cone | 1.97 | 0.71-0.80 | 0.37-0.76 |
| cylinder + elliptical nose | 1.99 | 0.68-0.81 | 0.33-0.71 |
| curved cone | 1.11 | 0.27-0.47 | 0.83-1.00 toward the throat |
| straight cone | 1.12 | n/a | **1.55-1.82 mid-length (focusing)** |

**Domed head + elliptical nose, grid series** (methane, point blast; peak with median):

| region | 2 cm | 1 cm | 0.5 cm |
|---|---|---|---|
| port face | 0.81 | 1.37 | 2.00 |
| head | 1.01 | 1.41 | 1.93 |
| cylinder | 0.82 (0.72) | 1.18 (1.02) | 1.68 (1.33) |
| nose | 0.60 | 0.69 | 1.03 |

The 2-D peaks rise ~40% per halving of the cell. Read 1 cm values as roughly half to two-thirds
of converged. **Ratios between cases are robust; absolute levels are not.**

**Deposition geometry, domed chamber, 1 cm** (methane, peak GPa; median in brackets):

| deposition | port | head | cylinder | nose | throat | worst |
|---|---|---|---|---|---|---|
| point blast | 1.37 | 1.41 | 1.18 (1.02) | 0.69 | 0.87 | 1.41 |
| **paper's plug + 27 km/s merged plug** | 0.93 | 1.22 | 0.74 (0.24) | **5.56** | **2.81** | **5.56** |
| + 15% eroded along the entry track | 1.42 | 0.67 | 0.54 | **5.86** | **2.92** | 5.86 |
| long column, all heat (ideal) | 0.85 | 0.55 | 0.83 | 0.98 | 0.29 | 0.98 |
| 2 / 5 kg port-side disk | 1.63 / 1.40 | 1.36 / 1.41 | 0.74 | 5.56 | 2.81 | 5.56 |
| 3 / 10 kg shell around the plug | 1.43 / 0.91 | 1.70 / 2.03 | 0.71-0.74 | 5.54 / 5.25 | 2.8 | 5.5 |

Sacrificial polyethylene placed beside or around the plug does nothing for the nose. Only mass
on the merged plug's path matters.

**The long plug, domed chamber, 1 cm** (methane). The plug's mass includes the paper's 4.4 kg,
and the extra comes out of the charge. The rod accretes it slice by slice: each slice releases
`p²/2 (1/m_before - 1/m_after)` as heat where it sits, and the merged body leaves the far end at
`p/(m_rod + m_plug)`.

| plug | port | head | cylinder | nose | throat | worst |
|---|---|---|---|---|---|---|
| 10 kg, 1.75 m | 0.43 | 0.39 | 0.56 (0.51) | 0.77 | 0.44 | 0.77 |
| **10 kg, 2.5 m, back-weighted** | 0.49 | 0.52 | 0.56 (0.46) | **0.50** | 0.61 | **0.61** |
| 20 kg, 1.75 m | 0.65 | 0.48 | 0.62 (0.55) | 0.72 | 0.41 | 0.72 |
| 40 kg, 1.75 m | 0.69 | 0.56 | 0.66 (0.57) | 0.65 | 0.34 | 0.69 |
| 40 kg, 2.5 m | 1.04 | 1.11 | 0.65 (0.56) | 0.75 | 0.50 | 1.11 |

18 cases in `vessel_blast_long_plug.jsonl`. One, 20 kg over 1 m back-weighted, put 1.82 GPa at
the throat, probably where this model places the whole merged body. Treat it as unconfirmed.

**Grid series for the best long plug** (10 kg, 2.5 m, back-weighted; `--long-plug-fine`):

| | 2 cm | 1 cm | 0.5 cm |
|---|---|---|---|
| worst on the wall | 0.38 GPa | 0.61 | 0.92 |
| cylinder (median) | 0.38 (0.34) | 0.56 (0.46) | 0.76 (0.61) |
| point blast, worst, same grid | (1.01) | 1.41 | ~1.9-2.0 |

The 2-D peaks are **not yet converging**: they rise ~50% per halving, and the increments have
not started to shrink. The 1-D sphere's did converge, at ~0.65x per halving. The long plug holds
at ~0.45x the point blast on every grid, so **its advantage is robust; its absolute level is
not**. Scaled from the converged sphere to the cylinder's closer radius, the converged long-plug
spike is estimated at **~1-1.5 GPa, possibly more**. Every absolute 2-D number below carries that
uncertainty.

**Volume trade** (`volume_trade.py`). Throat grown with volume at V/A* = 134 m. The Isp change
is `sqrt(eta_net ratio)`, where eta_net is the equilibrium blowdown efficiency less the stranded
H-atom store at the Bray freeze, at A/A* = 300.

| methane 7000 K | 20 m^3 | 40 | 60 | 80 | 120 | 160 | 240 |
|---|---|---|---|---|---|---|---|
| chamber pressure (bar) | 496 | 246 | 164 | 122 | 81 | 61 | 41 |
| eta_net | 0.538 | 0.511 | 0.496 | 0.485 | 0.470 | 0.460 | 0.446 |
| Isp vs 20 m^3 | 1 | -2.6% | -4.0% | -5.0% | -6.5% | -7.5% | -9.0% |
| long-plug spike (GPa, est. converged) | ~1.25 | 0.62 | 0.42 | 0.31 | 0.21 | 0.16 | 0.10 |
| exit diameter at A/A* 300 (m) | 7.5 | 10.7 | 13.1 | 15.1 | 18.5 | 21.3 | 26.1 |
| wall area factor | 1 | 1.6 | 2.1 | 2.5 | 3.3 | 4.0 | 5.2 |

| hydrogen 5500 K | 20 m^3 | 40 | 60 | 80 | 120 | 160 | 240 |
|---|---|---|---|---|---|---|---|
| eta_net | 0.858 | 0.783 | 0.723 | 0.681 | 0.631 | 0.599 | 0.557 |
| Isp vs 20 m^3 | 1 | -4.5% | -8.2% | -10.9% | -14.2% | -16.4% | -19.4% |
| centred spike (GPa) | 2.95 | 1.48 | 0.98 | 0.74 | 0.49 | 0.37 | 0.25 |

The chamber temperature does not move with volume: it is set by energy per kilogram. The
pressure does, and the loss is the recombination race in the nozzle. Blowdown stays 54-60 ms
(methane) and 39-44 ms (hydrogen), well inside the 250 ms pulse period.

**Dry wrap over a 10 mm steel liner** (`dry_wrap.py`; radial layer chain, wound layers in contact
only, each layer a hoop ring; the long plug's cylinder history scaled to each volume, grid factor
1-2):
- **Sizing is by strain compatibility.** The liner is bonded in strain to the wrap and may swing
  only +/-0.75 Y (0.42%). The wrap must stiffen the wall to that limit, so low-modulus fibre is
  heavy.
  - Vectran: 0.33 / 0.22 / 0.15 t/m^2 at 40 / 80 / 160 m^3.
  - Kevlar 49: 0.22 / 0.16 / 0.12 t/m^2.
  - Dry carbon tow: 0.17 / 0.13 / 0.10 t/m^2.
  - At 80 m^3 (~90 m^2 of wall) that is roughly 11-20 t.
- **The wrap itself survives.** Peak hoop strain is 0.8-0.9% (grid factor 1) or 1.4-1.55%
  (factor 2), against 2.1-3.3% to break. Lift-off gaps are 0.01-5 mm, layers fling at 30-200 m/s,
  and the steel's radial tension is negligible.
- **Hoop overshoot is settled by the shell model below, not by this ring.** This table's grid
  factor 2 scales the whole history, which doubles the spike's impulse. The hoop response depends
  on that impulse, and it is nearly grid-converged (next item). So the "factor 2" liner stresses
  (~2.1 GPa) overstate. The factor-1 column (770-815 MPa, D ~2.0 against the sizing's 1.7) is
  the realistic one.

**Hoop overshoot with the wall connected** (`shell_response.py`). The axisymmetric shell is a
beam on an elastic foundation along the meridian, clamped at the port ring and throat insert,
loaded by the 0.5 cm long-plug history at every station:

| wall | free ring D | connected shell D | sized for |
|---|---|---|---|
| monolithic steel, 20 m^3 | 1.74 | **1.66** | 1.7 |
| multilayer steel (8 shells), 20 m^3 | 1.74 | **1.63** | 1.7 |
| Kevlar dry wrap + 10 mm steel, 80 m^3 | 2.01 | **1.93** | 1.7 |

The hoop response depends on the spike's **impulse**, because the wall's 1-2 ms period is far
longer than the spike. The impulse converges where the peak does not:

| grid | 2 cm | 1 cm | 0.5 cm |
|---|---|---|---|
| cylinder peak | 0.38 GPa | 0.56 | 0.76 |
| cylinder 400 µs overpressure impulse (max) | 11.5 kPa s | 12.5 | 13.5 |

- **The paper's D = 1.7 holds for steel walls** under the long plug in the rocket-shaped chamber.
- **Dry wraps need D ~2.0-2.2.** They are stiffer for their mass, so they ring faster. That makes
  the wrap ~20-30% thicker: ~13-25 t at 80 m^3.
- **The earlier "D 3-6 at every volume" was an error** from doubling the impulse. Corrected here.
- **Most of the fix is the shape and the plug.** Coupling to the neighbours adds only ~5%, because
  under the long plug the spike is already spread along the wall.
- **Through-thickness effects (spall, delamination) follow the peak,** which is not converged. The
  overwrap conclusions stand.

**Pitch against volume** (methane; `near_term.wall_heat` at each volume, pitch scaled from the
20 m^3 anchors of 1.4-5.6 kg per pulse):

| | 20 m^3 | 40 | 80 | 160 | 240 |
|---|---|---|---|---|---|
| wall heat per pulse (MJ) | 184-644 | 222-894 | 242-1258 | 228-1785 | 204-2194 |
| pitch per pulse (kg) | 1.4-5.6 | 1.7-7.8 | 1.8-10.9 | 1.7-15.5 | 1.6-19.1 |

The low edge (H- radiation floor plus frozen-c_p convection) barely moves. The high edge (opaque
blackbody) grows with wall area.
- **Mass:** at 80 m^3 the extra pitch is +0.3-4.2 t per 780-pulse departure, against ~15-25 t of
  wall saved. Both are spent each departure, since the chamber is expended.
- **Isp:** the pitch leaves in the exhaust. 2-5.4 kg costs 6-15 s of effective Isp (near-term
  note), so ~11 kg at the high edge costs ~30 s (~4%). With the -5% from lower pressure, 80 m^3
  costs **~5-6% Isp at the low flux edge and ~8-9% at the high**.

**Wall materials** (wave model, the solved spike, `wall_waves.py`):
- **Thicker steel, less carbon** (methane, 25 -> 100 mm steel): carbon tension stays 6-7x the
  resin. The steel's own radial tension rises from 0.9 to 1.5 GPa. A no-overwrap wall of ~125 mm
  (35 t at the paper's D = 1.7, ~47 t at the solved D) keeps ~2x on spall.
- **Interlayers** (5-50 mm aluminium, titanium, elastomer, syntactic foam between steel and
  carbon): none lowers carbon tension for spikes of 10 µs or longer. Soft layers cut it 3x only
  for ~1 µs spikes, and raise the steel's tension.
- **Estimates, not sourced or simulated:**
  - Through-thickness reinforcement (z-pins, stitching, 3-D weaves) bridges delaminations; it
    does not stop the resin cracking each pulse. Carrying 0.3-1.3 GPa would take roughly 15-60%
    of the fibre running through the thickness.
  - A dry (unimpregnated) Vectran or aramid wrap cannot delaminate. Instead the reflected wave
    flings its outer layers outward at ~2p/Z (~400 m/s per GPa), which strains the hoop
    ~1-2% at 0.2-0.3 GPa and ~6% at 1 GPa, against ~3-3.6% to break.
  - Aluminium (5083 class) and magnesium alloys spall at roughly 1-1.5 and 0.6-1 GPa. Magnesium
    also forms a hydride in hot hydrogen and loses strength above ~150-250 °C.

## 3b. Caveat added 2026-10-07: the long plug's result rests on where the remnant's energy goes

The long-plug model leaves each collision's heat where the collision happened, and sends the
merged remnant on as a cold, 10 cm-wide body. Two bounding runs show that choice decides the
answer:

| plug (back-weighted) | ends at | cold remnant, worst | hot remnant*, worst |
|---|---|---|---|
| 4.4 kg (paper's mass), 2.5 m | 3.5 m (near the throat) | 0.48-0.67 GPa | ~7.9 GPa |
| 10 kg, 2.5 m | 3.5 m | 0.61 GPa | ~5.4 GPa |
| 10 kg, 1.5 m | 2.3 m (in the cylinder) | ~3.5 GPa | ~7.5 GPa |
| 40 kg, 1.5 m | 2.3 m | ~1.8 GPa | ~2.1 GPa |

\*The heat of the back half of the collisions travels with the remnant and is released at the
plug's end (`--long-plug-hot`, `--plug-clear-of-nose`; values scaled to the real chamber).

- **The good long-plug numbers came largely from the remnant leaving through the throat.** It
  starts close to the throat and, in the cold model, stays narrow. End the plug in the cylinder
  and the remnant crosses the chamber and strikes the nose, as the paper's short plug does. Only a
  very heavy plug (40 kg, 58% of the methane charge) slows it to ~2 GPa.
- **Momentum cannot be dissipated, only routed.** The rod's 187.5 kN s ends on the far wall or
  leaves through the throat. A plug converts kinetic energy to heat (p²/2M left moving), but not
  momentum.
- **Neither bound is the physics.** **Owed: a direct rod-into-plug simulation** in the 2-D solver.
  The 2.5 kg rod as dense material at 75 km/s penetrating a dense plug column, both
  hydrodynamic at these pressures, with no deposition model. It decides whether a long plug
  aimed at the throat lets the remnant leave cleanly, or whether a heavy plug is needed.
- **Until then, the §1 and §5 claim that "a ~10 kg long plug brings the wall to ~0.6 GPa" is
  unconfirmed.** The paper's short plug is still shown to hammer the nose at 5-6 GPa.

## 3c. The rod and plug as material (added 2026-10-07): harder than any deposition model

`--rod-plug` drops the deposition models. The 2.5 kg rod (polyethylene, 950 kg/m^3, 3.5 cm
radius, 0.7 m) moves in at 75 km/s and penetrates a plug column of material at rest. At these
pressures both are hydrodynamic. Peaks are in GPa, scaled to the real chamber; "out" is the share
of the rod's 187.5 kN s that leaves through the throat; the impulse is the worst wall station's,
over 400 µs:

| plug | grid | port | cylinder | nose | throat | out | max impulse |
|---|---|---|---|---|---|---|---|
| paper's 4.4 kg, 1.4-2.38 m | 1 cm | 1.41 | 0.87 | 2.53 | 3.69 | 10% | 53 kPa s |
| | 0.5 cm | 1.51 | 0.86 | **4.89** | 3.80 | 6% | 66 |
| 10 kg, 0.8-2.3 m, back-weighted | 1 cm | 0.36 | 1.32 | 1.26 | 2.02 | 41% | 43 |
| | 0.5 cm | 1.42 | 2.45 | 3.04 | **5.49** | 12% | **134** |
| 10 kg, 1.4-3.15 m | 0.5 cm | 1.61 | 1.33 | 5.40 | 2.02 | 3% | 83 |
| 10 kg, 1.0-3.5 m | 0.5 cm | 1.24 | 1.25 | 3.63 | 3.91 | 1% | 77 |

- **At 1 cm the rod is 3-4 cells wide.** Numerical diffusion smears it into a slow blob that
  spreads early and partly leaves through the throat. At 0.5 cm it stays coherent, penetrates as
  a jet, hits harder and more locally, and little escapes. **The 1 cm "spreading helps" result was
  largely an artifact.**
- **Every 0.5 cm case reaches 3-5.5 GPa somewhere,** at or past steel's ~2.9 GPa spall strength.
  The trend with resolution is unfavourable, so the 0.25 cm run decides the level.
- **The shell model's D ~1.6-1.7 was for the deposition models' loads** (impulse 12-14 kPa s). The
  direct rod-plug impulses are 5-10x larger and not grid-converged. **Hoop overshoot for the real
  rod-and-plug load is open again.**

The port column excludes the port opening (r < 7.5 cm). Its axis cell carries the rod's own ram
pressure, not a load on steel; an earlier version of this table included it.

**Grid convergence.** The 10 kg, 0.8-2.3 m case at 0.25 cm: port 1.20, cylinder 3.06, nose
3.90, throat 6.19 GPa; 4% out; 92 kPa s. **Peaks rise 15-25% per halving and are not
converged; impulse falls.** Read every 0.5 cm peak below as a floor.

### Nose, volume and dense layers (`--nose-study`, `--nose-study-40`, `--membrane`; 0.5 cm)

| case | port | cylinder | nose | throat | out | impulse |
|---|---|---|---|---|---|---|
| 20 m^3 long 12° cone (cone starts at z 1.29) | 0.98 | 0.39 | 4.82 | 2.07 | 21% | 161 |
| 20 m^3 narrow, r_c 1.2 m | 0.88 | 3.53 | 1.98 | 5.77 | 5% | 133 |
| 20 m^3 domed, rod enters through the port, plug 0.4-1.9 m | 0.85 | 2.13 | 3.61 | 3.49 | 2% | 103 |
| 20 m^3 long cone, rod enters | 0.98 | 0.94 | 4.18 | 1.31 | 29% | 132 |
| **40 m^3 long 12° cone** (throat 0.298 m^2), plug 0.8-2.3 m | 0.41 | **2.19** | **0.81** | **0.66** | 36% | 86 |
| 40 m^3 long cone, rod enters, plug 0.4-1.9 m | 0.57 | 2.57 | 1.41 | 0.34 | 34% | 60 |
| 20 m^3, 50% / 80% of charge in the last metre | 1.32 / 0.93 | 2.16 / 2.88 | 4.11 / 4.14 | 3.25 / 6.66 | 11-12% | 122 / 109 |
| 40 m^3, 50% / 80% of charge in 1 m before the cone | 0.54 / 0.40 | 3.22 / 4.62 | 0.61 / 0.84 | 0.43 / 0.26 | 62% / 25% | 48 / 52 |

- **The wall nearest where the jet stops takes 2-5 GPa, whatever the shape.** At 20 m^3 every
  nose shape leaves the plug near the convergence.
- **The 40 m^3 long cone clears the nose and throat** (0.3-1.4 GPa). Its 4.5 m cylinder keeps
  the plug 2.25 m from the cone. The remaining hot spot is a ~1 m band of cylinder beside the
  plug (z 1.8-2.8 m, 2.2-2.6 GPa at ~80 µs).
- **Dense end layers do not help the 20 m^3 nose.** Moving gas away from the plug makes the
  40 m^3 band worse.

### The plug band is a line blast: cushions fail, stand-off works (`--cushion-40`, `--shape-40`)

Scored as the worst back-face tension in the steel under the solved histories (`wall_waves`,
1-D elastic, free back), as a fraction of the 2.9 GPa spall strength:

| 40 m^3 long cone, 10 kg plug | volume | 105 mm steel | 10 mm liner + dry wrap | impulse |
|---|---|---|---|---|
| uniform fill (reference) | 40 | 0.68-0.72 | 0.26 | 86 |
| 50% / 80% of charge in a ring plug-to-wall, z 0.5-3 m | 40 | 0.76 / 0.75 | | 73 / 55 |
| 50% / 80% of charge in a layer at the wall (r > 0.9 m) | 40 | 1.12 / 1.58 | | 65 / 55 |
| wider, r_c 1.7 m, 20° cone | 40 | 0.70 | 0.20 | 60 |
| plug 2.5 m long (0.8-3.3 m) | 40 | 0.91 | 0.37 | 59 |
| bulge to r 1.8 m at z 1.4-2.8 (0.5 m ramps) | 47.6 | 0.85 | 0.25 | 35 |
| bulge to r 2.2 m at z 1.4-2.8 | 56.9 | 0.99 | 0.35 | 38 |
| bulge r 2.2 m at z 1.4-3.6, plug 2.5 m | 64.2 | 0.66 | 0.21 | 34 |
| bulge r 1.8 m at z 1.4-3.6, plug 2.5 m | 50.8 | 0.80 | 0.27 | 37 |

- **Gas or mist placed around the plug does not cushion it. A dense layer at the wall makes it
  worse.** The band is a strong line blast. Its shock pressure at radius R is
  `p ~ (E/L)/R^2`, independent of the ambient density. ~5.6 GJ stopped over ~1.5 m gives
  ~2 GPa at 1.4 m, as solved. A dense layer at the wall is driven into it like a slab.
  Small droplets follow the gas within millimetres, so mist behaves as this dense gas.
- **Only stand-off (R) and the stopping length (L) act.**
  - The bulge plateau follows `1/R^2`: ~0.5-1.1 GPa at r 2.2 m.
  - The bulge's steep closing shoulder (0.8 m over 0.5 m, ~68°) faces the blast. The merged
    body still moves downstream, so the shoulder takes 2.3-3 GPa.
  - A longer plug moves the stop downstream, not along. The rod stops where the back-weighted
    plug is densest.
- **Every bulge halves the impulse** (86 -> 34-38 kPa s), which is what the hoop overshoot
  follows.
- **A thin liner suffers less than a thick wall.** In 10 mm of steel the round trip is
  3.4 µs, about the spike's 2-4 µs, so the reflected tension overlaps the arriving compression
  and cancels. In 105 mm the round trip is 36 µs, so the whole spike reflects. This is the
  plate's rise-time rule.
  - Caveats: the spikes are sampled every 2 µs and sharpen with resolution. The liner's
    0.6-1.0 GPa per pulse is at or above Cr-Mo's static yield, so fatigue, not spall, decides.
- **Hoop overshoot under the rod-and-plug load** (`shell_response` on these histories):
  - 20 m^3 baseline: D 2.75 monolithic, 3.62 multilayer;
  - 40 m^3 long cone: D 4.1-4.7;
  - all against the 1.7 used for sizing.
  - It is concentrated in the plug band. The bulges' halved impulse should bring it down;
    rerun it on the final shape.

### A tapered bulge clears most of the band (`--shape-40b`; 0.5 cm)

The bulge closes over 2-3 m instead of 0.5 m (steepest wall ~32° or ~23°). It sits on the
40 m^3 layout, so the volume grows. The throat stayed at 0.298 m^2; it should be rescaled for
2 Hz (§3e).

| case | volume | cylinder peak | 105 mm steel | 10 mm liner + dry wrap | impulse | hoop D, shell (mono / x8) |
|---|---|---|---|---|---|---|
| **r 2.2 m, z 1.4-2.8, 2 m taper, plug 1.5 m** | 63.3 | 1.21 GPa | **0.45** | **0.14** | 26 | **2.09 / 2.35** |
| r 2.2 m, z 1.4-2.8, 3 m taper, plug 1.5 m | 67.6 | 1.45 | 0.51 | 0.21 | 37 | |
| **r 2.2 m, z 1.4-3.6, 2 m taper, plug 2.5 m** | 70.5 | 1.05 | **0.35** | **0.12** | 27 | 2.75 / 2.82 |
| r 1.8 m, z 1.4-2.8, 2 m taper, plug 1.5 m | 50.5 | 1.86 | 0.63 | 0.24 | 54 | |

The tension columns are back-face tension as a fraction of the 2.9 GPa spall strength.

- **Tapering the shoulder removes the shoulder hit.**
  - The r 2.2 m bulges bring the hottest wall to ~1-1.2 GPa, half the 40 m^3 cylinder's 2.2.
  - Steel tension falls to 0.35-0.45x spall, and 0.12-0.14x for the thin liner.
  - The impulse falls to 26-27 kPa s, from 86.
- **A 3 m taper is no better than 2 m.** Stand-off, not wall angle, is now the lever.
- **r 1.8 m is not enough.**
- **Hoop D falls to 2.1-2.8,** from 4.1-4.7, against 1.7 sized. This uses a uniform wall
  sized at the bulge radius, so the zone-by-zone wall study sizes it properly.
- **Chosen shape for the wall study: r 2.2 m bulge, z 1.4-2.8, 2 m taper, 1.5 m plug,
  ~63 m^3.** The shape:
  - a 1.4 m domed port head;
  - widening to 4.4 m across around the plug;
  - tapering over 2 m back to 1.4 m;
  - a 12° cone to the throat.
  The volume costs ~3% impulse per rod against 40 m^3 (volume trade).

## 3f. Which wall: zone by zone (`wall_zones.py`; added 2026-10-07)

`python -m puffsat.walled_nozzle.wall_zones [--histories PATH --case LABEL]`.
- **Hoop:** each zone (head, bulge ramp, plateau, taper, cone) is sized for D = 1.7 at its own
  radius. The allowable is the autofrettaged swing of 0.5625%, `dry_wrap.SWING`; earlier steel
  masses here used half of it. The coupled shell (`shell_response.respond_profile`, now with
  zoned stiffness and mass) is then solved under the real history, and any zone past the swing
  is thickened.
- **Spall:** the zone's hardest hits run through `wall_waves.simulate`, now with contact
  (unbonded) interfaces that pass compression, open in tension and track the gap. It runs over
  each spike only (-30/+100 µs), since the 1-D model has no hoop restoring force.
- **Cracks:** a 0.5 mm penny crack parallel to the face, on the air and CC2938 hydrogen curves,
  failing at K = 45.
- **Fling:** `dry_wrap.simulate`'s contact chain gives the hoop strain of layers thrown off by
  the spike. It is a free ring, so overstated.
- The requirement is **~4,200 pulses**: one departure of a 500-600 t stack.

| bulge (volume) | all steel | steel band, wrap elsewhere | multilayer (any) | **all dry wrap** |
|---|---|---|---|---|
| r 2.2 m (63 m^3) | 43.9 t: H2 cracks, 30 pulses | 35.3 t: H2, 131 | fling 1.3-1.4% vs 0.56% | 27.3 t: fling 3.95% vs 2.4% |
| r 2.6 m (78.5 m^3) | 45.2 t: H2, 301 | 35.4 t: H2, 232 | fling | 27.1 t: fling 2.48% vs 2.4% |
| **r 3.0 m (106 m^3)** | 47.1 t: H2, 662 (air 7,780) | 39.3 t: H2, 437 (air 5,600) | fling | **27.3 t: passes** |

Maraging 300 behind the alumina barrier, scored on the air curve with spall 4.1 GPa, `K_Ic`
66.5 and an assumed 1.2 GPa hoop swing:

| bulge | maraging band, wrap elsewhere | all maraging |
|---|---|---|
| r 2.2 m | 33.4 t, crack life 4,220 (1.0x) | fails, 1,802 |
| r 2.6 m | **34.1 t, 5,960 (1.4x)** | 42.5 t, passes |
| r 3.0 m | 37.2 t, 8,270 (2.0x) | 45.0 t, passes |

**Coupling the chain along the wall does not reduce the dry wrap's fling.**
`dry_wrap.simulate_coupled` couples the liner's bending along the meridian; a narrow struck band
flings ~5% harder. The spike hands its momentum to the wrap in microseconds, against ~0.4 ms for
bending, and lifted dry layers are free rings. So r 2.6 m's miss is real.

- **One spike, two fates.** The spike's few kPa s, delivered in microseconds, either reflects
  as tension in a solid wall (cracks) or throws an unbonded wall's outer layers outward
  (fling). The wall type moves the problem; only a smaller spike removes it.
- **The dry wrap wins: a 10 mm liner under Kevlar, 72-199 mm by zone.**
  - Its liner sees almost no spall (≤0.03x), because the spike is about its 3.4 µs round trip.
  - Its wrap tolerates fling up to 2.4%.
  - At r 3.0 m the worst fling is 1.69%. r 2.6 m misses by 3% on a free-ring estimate that
    shell coupling lowers by 25-40%, so it probably passes.
- **Steel fails mainly in hydrogen.** On the air curve (the alumina barrier) the steel band
  passes at r 3.0 m with a thin margin.
- **The taper is the hot zone in every design.** The blast still runs downstream into the
  closing wall.
- **Cost:** 106 m^3 is ~-7% impulse per rod against 40 m^3 (78.5 m^3: ~-5%), ~6 m across.
  27 t per chamber against the paper's 19 t, which is a bonded overwrap that delaminates.
- **Two designs survive:**
  - **Baseline:** all dry wrap at r 3.0 m (106 m^3, 27.3 t). It does not depend on hydrogen;
    the margin is fling, 1.4x.
  - **Fallback and smaller:** a maraging band with wrap elsewhere at r 2.6 m (78.5 m^3,
    34.1 t, -5% instead of -7% impulse per rod). Its margin is crack life, 1.4x, and it rests on
    an alumina barrier untested under GPa spikes.
- **Open:**
  - 0.5 cm peaks (+15-25% per halving);
  - the wrap's through-thickness stiffness (3 GPa assumed);
  - Kevlar's breaking strain under cycling;
  - maraging's hoop swing (assumed).
  - Run a 0.25 cm check of the r 3.0 m shape.

## 3g. One 5 kg chamber or two 2.5 kg chambers (`wall_zones --energy 2`; added 2026-10-07)

Inviscid flow has no length scale. Doubling the pulse and growing every length and time by
`2^(1/3)` leaves every pressure unchanged, so the solved r 3.0 m history serves the 5 kg rod
exactly, and the r 2.6 m one likewise. Only the fixed 10 mm liner and the sizing differ.

| wall, r 3.0 m bulge | one 2.5 kg chamber (106 m^3) | two of them | **one 5 kg chamber (212 m^3)** |
|---|---|---|---|
| all dry wrap | 27.3 t, fling 1.69% | 54.6 t | **52.9 t, passes, fling 1.97% (1.2x)** |
| maraging band (barrier), wrap elsewhere | 37.2 t | 74.4 t | 72.9 t, passes |

- **The same wall mass, ~3% lighter.** The steady load is `p V`, proportional to the energy.
- **21% less wall area** (`2^(2/3)/2`), so 21% less pitch, barrier and inspection. There is
  also one port, plug feed and membrane set instead of two.
- **The same chemistry and impulse per rod:** the pressure is the same.
- **The fling margin shrinks from 1.4x to 1.2x.** The 10 mm liner does not scale, so the wrap
  carries more of the hoop.
- **Size:** 212 m^3, ~7.6 m across at the bulge. That fits a 9 m Starship fairing but not
  smaller ones; ring segments either way (§3d).
- **The trade is redundancy.** Two chambers survive one failure at half thrust; one 5 kg
  chamber does not.
- At r 2.6 m the 5 kg dry wrap still fails on fling, as the 2.5 kg one does.

## 3d. Volume, manufacture and in-orbit assembly (added 2026-10-07; scaling, not a cost model)

**The steady-load steel does not grow with volume.** The charge's energy fixes `p V`
(`p ≈ (γ-1) E / V`), and a pressure vessel's minimum wall mass is `∝ ρ p V / σ`. So 20, 40
and 80 m^3 need about the same ~50 t of steel for the quasi-static load. The wall just gets
thinner, `t ∝ p r`:

| | 20 m^3 | 80 m^3, shape scaled up | 80 m^3, same radius, 4x longer |
|---|---|---|---|
| radius | 1.4 m | 2.2 m | 1.4 m |
| monolithic Cr-Mo wall | ~210 mm | ~85 mm | ~52 mm |
| wetted area | ~32 m^2 | ~80 m^2 (x2.5) | ~110-120 m^2 (x3.5) |
| steady-load steel | ~50 t | ~50 t | ~50 t |

- **What grows is area.** Pitch respray, the alumina barrier, inspection and wall heat all
  scale with it. Scaling the whole shape costs less area than lengthening at fixed radius.
- **The throat grows with volume** at fixed `V/A*`: x4 the area at 80 m^3, at a quarter of the
  pressure.
- **Estimated vessel cost at 80 m^3: ~1-1.5x if heavy steel and forging dominate; ~2-3.5x if
  area costs dominate. Never 4x,** since no major cost scales with volume itself.
- **Thin walls should be cheaper per tonne to make.** 50-85 mm rings are ordinary
  pressure-vessel and rocket-case practice. A 210 mm autofrettaged Cr-Mo forging is specialist
  work. *Unsourced here; owed to the parent (§5).*
- **The hammer does not shrink with volume.** The spike beside the plug depends on the plug's
  stand-off from the wall (2-2.6 GPa at 0.5 cm, r_c 1.4 m, at both 20 and 40 m^3). A larger
  chamber therefore becomes a **thin vessel with one armoured, possibly replaceable ring** at the
  plug station, as solid-rocket cases and gun chambers carry thick rings locally.

**In-orbit assembly: ring segments with circumferential joints only.**
- **Split the chamber into full hoops,** never into staves. A hoop carries the hoop stress
  `p r / t` without a joint. A circumferential joint carries only the axial load `p π r^2`,
  i.e. half the hoop stress. A longitudinal seam would carry the full hoop load and should be
  avoided.
- **Every ring fits a launch fairing.** The 20 m^3 chamber is 2.8 m across and the 80 m^3
  scaled one 4.4 m, so no ring needs splitting along its length.
- **Joint load.** The static axial load is `p π r^2`. At 20 m^3 that is ~305 MN; at 80 m^3
  scaled, ~190 MN, because the load falls as `V^(-1/3)` when the shape scales. A dynamic factor
  ~2 applies on top.
- **Joint type: a preloaded joint, not a separation clamp.**
  - Candidates are a bolted flange, a hub clamp (high-pressure process-piping practice), or an
    interrupted-thread "breech lock". The breech lock suits snap-together assembly: rotate a
    fraction of a turn and it is locked.
  - V-band ("Marman") clamps are separation joints. They hold by wedge friction and could slip
    and fret under 3,000 GPa-hammer pulses.
  - Preload above the separating load keeps the bolts' cyclic stress to a fraction of the
    pulse, which is the standard defence against joint fatigue.
  - Preload can be set and checked robotically: hydraulic tensioners, ultrasonic bolt-load
    measurement.
- **Where the joints go.** Put them in the quiet stations, never in the hammer band. In the
  40 m^3 cone the wall between the plug band and the cone sees 0.2-0.4 GPa. The armoured
  plug-band ring can itself be one segment, bolted between thin segments and swapped as a
  maintenance item.
- **Seals.** Each joint needs a hydrogen-tight metal seal (C-ring or similar) under the pitch,
  shielded from the hot gas. The seal, not the joint's strength, is the likelier weak point,
  and it is a test to list.

## 3e. Pulse rate, emptying, and two chambers (added 2026-10-07; estimates)

**The chamber must empty before the next rod.** A choked throat passes gas at the sound speed,
~3.7 km/s at 7,000 K. So the blowdown e-fold is `τ ≈ V / (0.6 A* a) ≈ 60 ms` at `V/A* = 134 m`.
The gas cools only slowly as it empties: recombination keeps the effective γ near 1.1-1.15. The
residue is still ~4,200 K at 2.5% of the charge (no wall loss in the model, so an upper bound).
`eta_blowdown` and the quoted Isp already follow this cooling.

Fraction of the charge left when the next rod arrives (methane 7,000 K, `near_term` blowdown):

| chamber | throat | at 0.25 s (4 Hz) | at 0.5 s (2 Hz) | nozzle exit, area ratio 100 |
|---|---|---|---|---|
| 20 m^3 | 0.149 m^2 | 2.7% | 0.2% | 4.4 m |
| 40 m^3 | 0.30 m^2 (`V/A*` 134 m) | 2.4% | 0.1% | 6.2 m |
| 40 m^3 | 0.25 m^2 (`V/A*` 160 m) | 4.0% | 0.3% | 5.6 m |
| 40 m^3 | 0.149 m^2 | 13% | 2.4% | 4.4 m |

- **The binding limit at 4 Hz is the refill, not the blowdown.** In what remains of the interval
  the chamber must be recharged (~60 kg), the plug placed, and the throat membrane resealed.
  The plug (frozen methane, below 90.7 K) would meet ~4,000 K residue. **Baseline: 2 Hz per
  chamber**, throat ~0.25 m^2 at 40 m^3, leaving ~200 ms for the refill.
- **A smaller throat costs little chemistry** (`eta_net` -0.5-1%). It erodes the throat only
  +10-20% per pulse, because the heat flux scales as `p^0.8` and the pressure falls. **But it
  multiplies the wall heat** by the longer dwell (x2-4 at a fixed 0.149 m^2), and `eta_net`
  does not debit that heat.
- **A hotter chamber does not pay.** At 10,000 K the charge falls to ~45 kg: +11-13% Isp,
  -13-17% impulse per rod, and 1.6-3x the wall heat (6-15 kg of pitch per pulse at 40 m^3).
  Counted per kilogram consumed, pitch included, the two temperatures are about even. Keep
  7,000 K.

**The pulse rate sets the departure burn's Oberth loss.** A finite burn centred on the 600 km
periapsis of the turnaround ellipse (613,000 km apoapsis), 5.43 km/s, ~644 kN s per pulse,
constant thrust along the velocity, on the paper's 500-600 t departing stacks:

| total rate | thrust | burn | extra Δv against an impulsive burn |
|---|---|---|---|
| 1 Hz | 0.64 MN | 50-60 min | +800-965 m/s (+15-18%) |
| 2 Hz | 1.29 MN | 25-30 min | +335-433 m/s (+6-8%) |
| 4 Hz | 2.58 MN | 13-15 min | +108-149 m/s (+2-3%) |

The periapsis passage takes about `r_p / v_p ≈ 11 min`. A burn much longer than that spends its
impulse where the ship is slower. **Two chambers at 2 Hz each give the 4 Hz loss, with each
chamber keeping the 0.5 s it needs to empty and refill.**

**Whether a second chamber pays depends on the wall type.**
- **The paper's 19 t methane wall is a bonded carbon overwrap** (18 cm), and §3 finds it
  delaminates under the spike. Its 130 t wall (10 kg hydrogen pulse, 80 m^3) is the same
  autofrettaged overwrap; steel alone there was 385 t. Wall mass scales about linearly with the
  pulse energy (`p V`), so there is no hidden size penalty.
- **On the survivable walls:**
  - all steel: multilayer Cr-Mo, ~35-50 t per 2.5 kg chamber at any volume;
  - dry Kevlar over a 10 mm liner: ~11-20 t at ~80 m^3, sized at D 2.0-2.2.
- **Second chamber, all steel:** +35-50 t on a 500-600 t stack (7-9%), against ~4-5% of Δv
  saved. About a wash or worse.
- **Second chamber, dry wrap:** +13-22 t (2-4%), with extension. It pays.
- **One chamber at 5 kg per rod** is the same total wall as two 2.5 kg chambers (`p V` scales
  with energy). Cube-root scaling leaves the hammer, spall and hoop overshoot unchanged. It has
  2^(2/3)/2 = 79% of the wall area of two chambers, so less pitch and wall heat, and one set of
  port, plug and membrane hardware. It gives up redundancy.
- **The stack mass per chamber is inconsistent in the parent.**
  - The chamber requirement is 630-780 pulses per departure over 159-194 s at 4 Hz
    (`sec:methane_7000_near_term`; used by `chamber_fatigue.py`).
  - The 5.43 km/s burn at ~780 s consumes about half the stack, ~65 kg per pulse with charge,
    rod and pitch. So 780 pulses is a ~100 t craft. The 500-600 t stacks of the growth ledger
    need ~4,000-4,300 pulses (~17 min at 4 Hz).
  - On ~100 t the Oberth loss above is negligible at any rate.
  - On 500-600 t it applies, and the fatigue budget is 5-6x the 780 pulses
    `chamber_fatigue.py` assumed.
  - The parent should state which.
- **Unchecked for the dry wrap: spall of its 10 mm liner.** Its back face rides on an unbonded
  wrap, close to a free surface, and the plug-band spikes (2-4 µs) are about its 3.4 µs round
  trip.

## 4. Limits

- **γ-law gas.** Real-gas equilibrium enters only through Γ_eff and the QSP scaling. A
  real-EOS run is owed.
- **Grid.** 2-D spikes are not grid-converged; the sphere's are. 2-D levels are probably 1.5-2x
  higher than the 1 cm tables show.
- **Deposition.** The long-plug model places the hot accreted mass where it was struck. A
  penetration model (rod erosion and crater flow) would refine where along the axis the heat
  goes.
- **Wall response.** The per-station ring in `vessel_blast.py` (thin wall, point blast) overstates
  the overshoot. `shell_response.py` settles it for the long plug: D 1.63-1.93. It neglects
  meridional curvature of the head and nose.
- **Hydrogen.** Not run with a long plug.
- **Volume trade.** Spikes are scaled as 1/V from the 20 m^3 runs, not re-solved. Nozzle
  extension mass grows with the exit diameter squared.

## 5. What the paper should change (owed to `Balloon-Pulse-Propulsion`)

- **`sec:material_chamber_plug`, the stand-off paragraph** ("the 0.2 GPa at the wall lasts
  microseconds ... the shell's response to it has not been checked").
  - Replace it: the plug as sized leaves 36% of the rod's energy (2.5 GJ) in the merged 6.9 kg
    at 27 km/s. Solved in 2-D, that strikes the convergent nose at 5.6-5.9 GPa and the throat at
    ~2.8 GPa on every pulse, about twice steel's ~2.9 GPa spall strength.
  - The fix: a longer plug of ~10 kg over 1.75-2.5 m of the axis, denser at the back. Make it
    **frozen methane inside a polyethylene container** rather than polyethylene foam:
    - CH4 carries twice polyethylene's hydrogen per carbon, so in the methane chamber the core
      matches the charge and costs essentially no exhaust speed. A polyethylene plug costs ~0.9%.
    - The container is the skin the paper already sizes to carry the plug's sideways pull as the
      port window tracks the rod (1.3-5 mm for the 2.5 kg rod). Solid methane alone is weak and
      would break up.
    - At solid methane's ~500 kg/m^3, 10 kg over 2.5 m is a column ~5 cm in radius, close to the
      paper's 5.5 cm plug width.
    - It must stay below methane's 90.7 K triple point until impact.
    - In the hydrogen chamber a methane core still raises the mean molecular mass (+8.8% for the
      paper's plug width, against +11% for polyethylene). It is the lesser cost, not a free one. It turns the rod's
    energy into heat along the axis, so the nose falls below ~1 GPa and the wall loads become
    roughly uniform at about half the centred-blast level.
  - The penetration-depth rule sizes the minimum plug. The plug should be longer than that
    minimum, so the merged body keeps colliding after the rod is consumed.
- **`sec:carbon_overwrap` and `tab:layered_wall_mass`.**
  - At 20 m^3 the reflected blast spike (centred ~36x the quasi-static pressure; ~1-1.5 GPa
    under the long plug) exceeds the overwrap's through-thickness strength (64 MPa resin) by
    about 5-20x everywhere. The aluminium layer and more steel do not change that.
  - The near-term wall should be steel: multilayer shrink-fitted Cr-Mo, about 35-47 t for
    methane. Otherwise the chamber must grow.
  - A **dry (unbonded) wrap over a thin steel liner** cannot delaminate. Its flung layers stay at
    0.8-1.5% hoop strain against 2.1-3.3% to break. At ~80 m^3 it weighs ~11-20 t (Kevlar 49 or
    dry carbon lighter than Vectran, because sizing is by stiffness). The cost is ~5-9% Isp,
    including the extra pitch.
  - A bonded overwrap needs ~120-160 m^3.
  - Hoop sizing: D = 1.7 holds for steel walls (shell model 1.63-1.66); dry wraps need ~2.0-2.2.
  - State that the "softens the local shocks" role given to the liquid share does not hold. The
    drops equilibrate in microseconds and pay the same toll as the gas. The liquid serves better
    as axial plug mass.
- **Chamber geometry.**
  - Draw the chamber as a rocket chamber: a 2:1 elliptical port head and an elliptical
    convergence, not a sphere or a straight cone. A flat port face doubles the corner spike; a
    straight cone focuses mid-length.
  - `fig:rod_port_plug` should show the longer plug.
- **`sec:steel_chamber_service`, `eq:ramp_dlf`.** The blast's reverberation gives D = 2.0-2.3 in
  a centred sphere, not 1.7. Size with the solved D, or with the long plug's lower, spread load.
- **Hydrogen barrier on the chamber steel** (`docs/hydrogen_and_wall_steels.md`). Dissociated
  methane puts atomic hydrogen on the methane chamber's wall during every pulse, and the pitch is
  not a designed barrier. Propose a thin dense alumina film on the Cr-Mo, aluminized and
  oxidized before the final temper, under the pitch, and under the GRCop-84 liner in the hydrogen
  chamber. Its permeation reduction under pulsed plasma is a test to list. Keep the parent's
  sub-950 MPa Cr-Mo qualified on CC2938. Maraging is "extreme" in hydrogen and is not proposed
  for the chambers.
- **Chamber size, manufacture and assembly (§3d). Please find supporting references.**
  - The steady-load steel is ~constant with volume (`p V` fixed by the charge energy). A larger
    chamber trades a thick forging for thin rings plus more surface area. Estimated
    ~1-1.5x vessel cost at 80 m^3 if steel and forging dominate, ~2-3.5x if area costs
    dominate. The Isp debit from area is separate.
  - **References wanted:**
    1. that thick (~200 mm) high-strength forgings, autofrettaged or heat-treated, cost more per
       tonne than ~50-85 mm rings, and how much more (forging-size limits, supplier base,
       through-hardening, inspection);
    2. segmented pressure-vessel or rocket-case manufacture: SRM case segments and their
       clevis joints, multilayer vessels;
    3. preloaded circumferential joints under cyclic pressure: hub clamps, bolted flanges,
       interrupted-thread breech closures, and metal seals for hydrogen.
  - State that the chamber can be assembled in orbit from full-hoop rings joined
    circumferentially. That means preloaded joints placed away from the plug's hammer band,
    with the hammer-band ring as a replaceable armoured segment. No longitudinal seams.
- **Pulse rate and chamber count (§3e).**
  - State the chamber's pulse rate. It must empty and refill between rods.
  - 4 Hz per chamber is refill-limited: plug placement in ~4,000 K residue, recharge, and the
    throat membrane reseal. Propose **2 Hz per chamber**, with the throat scaled to the volume
    (`V/A*` ~134-160 m).
  - **Use two chambers on the departing stack.** The departure burn's finite-burn loss is
    +15-18% Δv at 1 Hz, +6-8% at 2 Hz and +2-3% at 4 Hz on a 500-600 t stack. Two chambers at
    2 Hz recover the 4 Hz figure, and each still has 0.5 s to empty.
  - **One 5 kg chamber matches two 2.5 kg chambers on wall mass** (52.9 against 54.6 t,
    dry wrap, §3g), with 21% less area and one set of mechanisms. It gives up redundancy and
    needs a ~7.6 m diameter. Ask: which does the mission need?
  - Recheck "on the 500-600 t stacks that depart, one chamber is best" with the finite-burn
    loss charged. The ledger does not appear to include it.
    - The answer turns on the wall: a second chamber costs ~35-50 t all steel (about a wash),
      or ~13-22 t as a dry wrap (it pays).
    - One chamber at 5 kg per rod is the alternative: the same total wall, ~21% less wall area,
      no redundancy.
    - The paper's 19 t and 130 t walls are bonded overwraps, which this study finds delaminate.
      Scale from the survivable walls.
  - List the refill sequence at the chosen rate as an open engineering item: plug dwell in the
    residue, mist timing, membrane reseal.
  - **Reconcile the departing mass per chamber.** 630-780 pulses (159-194 s at 4 Hz) is a
    ~100 t craft, while the growth ledger departs 500-600 t stacks (~4,000-4,300 pulses). The
    finite-burn loss, the chamber count and the fatigue budget all depend on which.
- **Chamber wall (§3f).** Make the near-term methane chamber a **dry Kevlar wrap over a
  10 mm Cr-Mo liner**, ~27 t, on a bulged rocket shape:
  - an r 2.6-3.0 m bulge around the plug, closed by a 2-3 m taper;
  - 78-106 m^3, throat scaled for 2 Hz;
  - solid or multilayer steel at the plug band fails (hydrogen cracks, or flung shells).
  - Replace the 19 t bonded-overwrap figure.
  - Name the fallback: a maraging 300 band behind a thin alumina barrier, wrap elsewhere,
    at r 2.6 m (~34 t). List the barrier's spike test as its gating item.
- **Hydrogen at the 2.5 kg rod.** No wall under ~50 t survives the centred spike (~2.95 GPa, at
  steel's spall strength). Say so. The levers are a larger chamber (~8-11% Isp at 60-80 m^3) or a
  smaller rod per pulse.

## 6. Sources

- [KT07] J. R. Kamm and F. X. Timmes, "On efficient generation of numerically robust Sedov
  solutions", LA-UR-07-2849 (2007): α = 0.851 (γ 1.4, spherical).
- [PE87] M. Pilch and C. A. Erdman, Int. J. Multiphase Flow 13, 741 (1987): drop breakup time
  of 5-6 d/u √(ρ_l/ρ_g). Cited by the parent as `pilch1987_breakup`.
- Parent `templateArxiv.tex`: `sec:material_chamber_plug`, `sec:carbon_overwrap`,
  `sec:batch_filled_steel`, `tab:drop_survival`, `sec:steel_chamber_service`.
- Spall and fatigue data: as in
  [`walled_nozzle_chamber_fatigue.md`](walled_nozzle_chamber_fatigue.md) and
  [`spray_plate_fatigue_and_spall.md`](spray_plate_fatigue_and_spall.md).
