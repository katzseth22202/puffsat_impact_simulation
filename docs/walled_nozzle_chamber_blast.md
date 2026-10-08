# Walled nozzle: the blast inside the near-term chamber, and what fixes it

Status: analysis note, 2026-10-07; bottom line and paper asks updated 2026-10-08. Follows [`walled_nozzle_chamber_fatigue.md`](walled_nozzle_chamber_fatigue.md),
which found that a reflected blast spike, not the mean pressure, decides whether the chamber wall
survives. This note solves that spike. The chamber is the parent's near-term one: the 2.5 kg rod
at 75 km/s (7.03 GJ) in 20 m^3, methane at 7000 K (496 bar) and hydrogen at 5500 K (818 bar).

Reproduce:
- `make walled-nozzle-chamber-blast`: 1-D spherical runs, `euler2d::sphere`.
- `make walled-nozzle-vessel-blast`: 2-D axisymmetric chamber shapes, `euler2d::vessel`.
- `PYTHONPATH=python uv run python -m puffsat.walled_nozzle.volume_trade`.

## 1. Bottom line (current design, 2026-10-08)

**The design that survives one departure (~4,200 pulses), methane at 7000 K, 2.5 kg rod:**

| part | choice |
|---|---|
| plug | 10 kg of frozen methane in a polyethylene container, a 1.5 m column at 0.8-2.3 m from the port, 3x denser at the back |
| shape | domed port head, r 1.4 m; a bulge to **r 3.0 m** around the plug (z 1.4-2.8 m), opened over 0.5 m and closed by a **3 m taper**; a long 12° cone to the throat. **~106 m^3, ~6 m across** |
| wall | a **10 mm Cr-Mo shell** (gas-tight boundary, pitch substrate) under a **dry (unbonded) Kevlar 49 wrap**, 72-199 mm by zone. **~25-27 t.** Pitch on the steel; a thin alumina film as defence in depth |
| rate | **2 Hz per chamber**, throat scaled to the volume (`V/A*` ~134-160 m) |
| count | **one 5 kg chamber** (decided 2026-10-08): 212 m^3, ~7.6 m across, wall 48.2 t, ~53-63 t with the nozzle extension (§3h) |
| stack | 500-600 t departing; ~2,000-2,650 pulses, 17-21 min; finite-burn loss +3-5% Δv |
| cost | ~-7% impulse per rod against 40 m^3 (~-11% against 20 m^3); effective Isp ~flat |

**How the study got there:**
1. **The paper's 4.4 kg plug sends a 5-6 GPa hammer into the nose** (§3), against steel's
   ~2.9 GPa spall strength.
2. **With the rod and plug solved as material, every 20 m^3 shape takes 3-6 GPa** where the
   jet stops (§3c). The peaks rise 15-25% per grid halving.
   - The earlier deposition models' "long plug fixes it" (§3, §3b) was an artifact of where
     they put the heat.
3. **The plug band is a strong line blast,** `p ~ (E/L)/R^2`. It does not depend on the gas
   density, so no arrangement of charge, mist or membrane-held gas cushions it (§3c). A dense
   layer at the wall makes it worse.
   - Only stand-off acts: a wider bulge at the plug, closed gently. A steep wall facing
     downstream takes the remnant's push.
   - The long 12° cone keeps the nose and throat under ~1 GPa.
4. **The wall decides between tension and fling** (§3d). The spike's few kPa s, delivered in
   microseconds, either reflects as tension in a solid wall or throws an unbonded wall's layers
   outward.
   - Solid Cr-Mo cracks in hydrogen within 30-660 pulses.
   - Multilayer steel shells fling past their 0.56% swing.
   - Bonded carbon overwrap delaminates everywhere (64 MPa resin).
   - **The dry wrap passes at r 3.0 m:** its thin shell cancels the spike (≤0.03x spall), and
     its flung layers reach 1.68% against 2.4% (1.43x). This is confirmed at 0.25 cm.
   - **Fallback:** a maraging 300 band behind an alumina barrier, at r 2.6 m (~34 t). It is
     gated on a barrier that is untested under GPa spikes.
5. **Volume costs chemistry, not steel** (§3f). The steady-load steel is fixed by `p V`, so it
   does not grow with volume; the wall thins. What grows is area: pitch, barrier, inspection.
   Early escape of unthermalized gas costs only ~0.4% (§3g).
6. **Rate and count** (§3g, §3h): 4 Hz is refill-limited, so use 2 Hz per chamber.
   - One 5 kg chamber matches two 2.5 kg chambers on wall mass, with 21% less area and one set
     of hardware.
   - On the 500-600 t stack, at net thrust, its finite-burn loss is +3-5% Δv.
   - 10,000 K does not pay; keep 7000 K.
7. **The liquid-share spray does not soften the wall shocks** (§3). The drops equilibrate in
   microseconds and pay the gas's toll. Spend the liquid as plug mass.
8. **The hydrogen chamber at this pulse size is unresolved.** It has not been run with the rod
   and plug as material or with the bulge.

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

*§3 and §3b are the first, deposition-model results, kept as the record. Where §3c onward
differ (the plug's effect, the shape, the wall), the later sections and §1 stand.*

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
2 Hz (§3g).

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

## 3d. Which wall: zone by zone (`wall_zones.py`; added 2026-10-07)

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
- **The r 3.0 m dry wrap holds at 0.25 cm** (`--shape-40c-fine`; `wall_zones ... --suffix fine`).
  - The band peak falls 0.89 -> 0.67 GPa and the impulse 10 -> 9 kPa s.
  - The cone near the throat picks up a later 0.68 GPa hit at ~380 µs (spall 0.04x, fling
    1.01%).
  - The worst fling is unchanged, 1.69% -> 1.68% (1.43x against 2.4%), and the wrap sizes to
    24.9 t.
  - The maraging-band fallback passes at 35.2 t.
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

## 3e. One 5 kg chamber or two 2.5 kg chambers (`wall_zones --energy 2`; added 2026-10-07)

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
  smaller ones; ring segments either way (§3f).
- **The trade is redundancy.** Two chambers survive one failure at half thrust; one 5 kg
  chamber does not.
- At r 2.6 m the 5 kg dry wrap still fails on fling, as the 2.5 kg one does.

## 3f. Volume, manufacture and in-orbit assembly (added 2026-10-07; scaling, not a cost model)

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

## 3g. Pulse rate, emptying, and two chambers (added 2026-10-07; estimates)

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

**Early escape costs ~0.4%** (`--efficiency-r3`: the r 3.0 m chamber, throat 0.66 m^2 for 2 Hz,
5 ms at 0.5 cm, recording each outflow parcel's stagnation enthalpy `h0`).
- In the first 5 ms, 5.5 kg (7% of the charge) leaves carrying 0.79 GJ (11% of the rod's
  energy).
- Its uniformity `(sum dm sqrt(2 h0))^2 / (2 M E)` is 0.956.
- Adding the rest of the charge's ideal isentropic blowdown gives 1,023 kN s, against 1,027 kN s
  had everything thermalized first: a ratio of 0.996 (gamma-law, ideal nozzle).
- So assuming full thermalization in the volume trade is sound. A bigger chamber gains nothing
  from thermalizing faster; its cost remains the lower-pressure chemistry.

**The pulse rate sets the departure burn's Oberth loss.** *(Corrected in §3h: the table below
uses gross exhaust momentum, ~644 kN s per pulse. The head-on rod debit leaves ~450-490 kN s, so
the thrust is ~30% lower and the losses higher; one 5 kg chamber at 2 Hz gives +3-5%.)* A finite burn centred on the 600 km
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
  - dry Kevlar over a 10 mm liner: ~11-20 t at ~80 m^3 sized at D 2.0-2.2 (early estimate);
    solved zone by zone, **~25-27 t at 106 m^3** (§3d).
- **Second chamber, all steel:** +35-50 t on a 500-600 t stack (7-9%), against ~4-5% of Δv
  saved. About a wash or worse.
- **Second chamber, dry wrap:** +~27-29 t with its extension, ~5% of the stack, against the
  ~4-5% of Δv saved. Roughly a wash; the parent's ledger decides. (An early 80 m^3 estimate of
  13-22 t said it paid; the solved r 3.0 m wall is heavier.) One 5 kg chamber gets the 4 Hz loss
  without a second set of hardware (§3e).
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

## 3h. The chosen configuration and its masses (decided 2026-10-08)

**Decision:** a 500-600 t departing stack and **one 5 kg chamber** at 2 Hz. It is the r 3.0 m
design scaled by `2^(1/3)`: 212 m^3, ~7.6 m across at the bulge, throat 1.33 m^2
(`V/A* = 160 m`).

| item | mass | basis |
|---|---|---|
| wall: dry Kevlar 49 over a 10 mm Cr-Mo shell | **48.2 t** (52.9 t on the 0.5 cm history) | `wall_zones --energy 2 --suffix fine`; worst fling 1.91% against 2.4% (1.26x) |
| fallback wall: maraging band (barrier) + wrap | 68.8 t | same |
| nozzle extension, area ratio 100 / 300 | 4.9 t (13.0 m exit) / 14.7 t (22.5 m exit) | the ledger's 2.19 t at 8.7 m, scaled by exit area |
| throat insert, port window, plug feed, membranes | not estimated | |
| **carried per pulse** | **~115-117 kg** methane, plug included, **+3.5-24 kg pitch** | `near_term` at 106 m^3 x2 (~0.55 kg/m^3, ~95 bar at 7000 K); pitch scaled from the 20 m^3 wall heat |
| rod per pulse | 5 kg, delivered | |
| net impulse per pulse | 906 / 971 kN s (AR 100 / 300) | `Isp_eff` 788 / 845 s on the carried charge and plug, head-on debit taken |

The departure burn, 5.43 km/s at the 600 km periapsis, constant thrust at 2 Hz, finite burn
charged:

| stack | thrust | burn | extra Δv (finite burn) | consumed | pulses |
|---|---|---|---|---|---|
| 500 t | 1.81-1.94 MN | 16-18 min | +167-196 m/s (3.1-3.6%) | 251-291 t | ~2,000-2,200 |
| 600 t | 1.81-1.94 MN | 19-21 min | +225-262 m/s (4.1-4.8%) | 304-352 t | ~2,400-2,650 |

- **The chamber hardware is ~53-63 t:** wall plus extension, ~10% of the stack. The paper's
  ledger has a 19 t bonded-overwrap wall and a 2.19 t extension per 2.5 kg chamber.
- **~2,000-2,650 pulses per departure,** about half the 4,200 the wall was checked for. The
  fatigue and crack margins grow.
- **The finite-burn loss is 3-5%, not the 2-3% quoted earlier in §3g.** That figure used gross
  exhaust momentum (~644 kN s per 2.5 kg pulse). The rod arrives head-on, so the net is
  ~450-490 kN s per 2.5 kg pulse.
- **Pitch is the widest uncertainty:** 7-60 t per departure between the low and high
  wall-heat edges. Its high edge comes from the gas's radiation, which the solved wall heat does
  not yet resolve.
- **Area ratio 300 is worth ~5% impulse but needs a ~22 m exit**, against ~13 m at 100.

## 4. Limits

- **Gamma-law gas.** Real-gas equilibrium enters only through `Gamma_eff` and the QSP scaling. A
  real-EOS run is owed.
- **Grid.** 0.5 cm peaks rose 15-25% per halving. The chosen shape is checked at 0.25 cm, where
  the fling held; others are not.
- **Wall waves are 1-D and elastic.** The spike's footprint would spread in 2-D, so the tensions
  are upper bounds. The fling chain is a free ring, and liner coupling was shown not to reduce
  it.
- **Material inputs assumed:**
  - the dry wrap's through-thickness stiffness (3 GPa);
  - Kevlar's 2.4% break, with no knockdown for cycling or fretting;
  - maraging's hoop swing (1.2 GPa).
  These are the bench tests of §5.
- **Not modelled:** the throat insert; ring joints and seals; meridional curvature of head,
  bulge and nose in the shell; hydrogen uptake through pitch (the CC2938 curve at chamber pressure
  stands in).
- **Hydrogen chamber.** Not run with the rod and plug as material.
- **The departing mass per chamber is inconsistent in the parent** (§3g). The fatigue budget
  here uses 4,200 pulses (500-600 t).

## 5. What the paper should change (owed to `Balloon-Pulse-Propulsion`)

- **`sec:material_chamber_plug`: the plug and the stand-off paragraph** ("the 0.2 GPa at the
  wall lasts microseconds ... the shell's response to it has not been checked").
  - Replace it. The 4.4 kg plug leaves 36% of the rod's energy (2.5 GJ) in a merged body at
    27 km/s, which strikes the nose at 5-6 GPa.
  - Even a 10 kg plug only moves the problem. The rod stops where the plug is densest, and the
    wall beside that point takes a line blast of ~2 GPa at 1.4 m radius. The shape and wall
    below are what fix it.
  - Make the plug **10 kg of frozen methane inside a polyethylene container**: a 1.5 m column
    at 0.8-2.3 m from the port, denser at the back, ~6.5 cm in radius at 500 kg/m^3.
    - CH4 carries twice polyethylene's hydrogen per carbon, so in the methane chamber the core
      matches the charge and costs essentially no exhaust speed. A polyethylene plug costs ~0.9%.
    - The container is the skin the paper already sizes to carry the plug's sideways pull as the
      port window tracks the rod (1.3-5 mm). Solid methane alone is weak and would break up.
    - It must stay below methane's 90.7 K triple point until impact. At 2 Hz it meets ~4,000 K
      residue (2.5% of the charge at ~3 bar), so it goes in at the last moment.
    - In the hydrogen chamber a methane core still raises the mean molecular mass (+8.8% against
      +11% for polyethylene): the lesser cost, not a free one.
  - If a heavier plug is ever needed, water is the candidate, not more methane.
- **Chamber geometry (`fig:rod_port_plug`, the chamber sketch).** Draw:
  - a domed port head at r 1.4 m;
  - a **bulge to r 3.0 m around the plug, closed by a ~3 m taper**;
  - a long 12° cone to the throat;
  - ~106 m^3, ~6 m across.

  Not a sphere, a flat port face (it doubles the corner spike), a short convergent nose (it
  takes the remnant), or any steep wall facing downstream.
  - The bulge works because the band's shock pressure falls as `1/R^2`. The taper keeps the
    downstream-moving remnant from striking a shoulder.
- **`sec:carbon_overwrap`, `tab:layered_wall_mass`, `sec:steel_chamber_service`: the wall.**
  - Replace the bonded overwrap with a **dry (unbonded) Kevlar 49 wrap over a 10 mm Cr-Mo
    shell**: ~25-27 t, 72-199 mm of wrap by zone, pitch on the steel.
  - It passes hoop, spall, crack growth and fling for one 4,200-pulse departure. The worst
    fling is 1.68% against 2.4%, and the shell sees ≤0.03x spall, so hydrogen barely matters
    to it.
  - The bonded overwrap delaminates under the reflected spike everywhere: 64 MPa resin against
    0.1-1 GPa spikes, even in the quiet zones.
  - Solid Cr-Mo, the paper's steel, cracks in hydrogen at the plug band within 30-660 pulses.
  - Multilayer shrink-fitted steel flings its outer shells past the swing.
  - The paper's 19 t (methane) and 130 t (10 kg hydrogen) walls are bonded overwraps; replace
    them. Wall mass scales about linearly with pulse energy.
  - **Size the hoop with the solved shell,** zone by zone, not `eq:ramp_dlf`'s D = 1.7. The
    solved overshoot varies by zone and wall type (§3d).
  - **Fallback:** a maraging 300 band at the plug behind the alumina barrier, dry wrap
    elsewhere, at r 2.6 m (78.5 m^3, ~34 t, ~2% more impulse per rod). Its crack-life margin
    (1.4x) holds only if the barrier works under GPa spikes.
- **Hydrogen and the steel (`docs/hydrogen_and_wall_steels.md`).**
  - Atomic hydrogen from dissociated methane, and H2 from the pitch's own pyrolysis, reach the
    steel every pulse. The steel sits at ~350-420 K, where embrittlement is strongest. The gas's
    temperature does not protect it.
  - Put a sub-micron alumina film on the Cr-Mo shell under the pitch, as defence in depth. It is
    **required** for the maraging fallback.
  - List the barrier's uptake under pulsed hot dissociated gas as a test.
  - Keep sub-950 MPa Cr-Mo, qualified on CC2938, for the shell.
- **Chamber size, manufacture and assembly (§3f). Please find supporting references.**
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
    circumferentially, with preloaded joints placed in the quiet zones (the cylinder between head
    and bulge, and the cone), the bulge as its own segment, and no longitudinal seams.
- **Pulse rate and chamber count (§3g).**
  - State the chamber's pulse rate. It must empty and refill between rods.
  - 4 Hz per chamber is refill-limited: plug placement in ~4,000 K residue, recharge, and the
    throat membrane reseal. Propose **2 Hz per chamber**, with the throat scaled to the volume
    (`V/A*` ~134-160 m).
  - **Decided: one 5 kg chamber on a 500-600 t stack (§3h).** The masses are there: ~53-63 t of
    hardware, ~117 kg carried per pulse plus pitch, ~2,000-2,650 pulses, +3-5% finite-burn
    loss at net thrust.
  - **The carried mass is energy-sized at 7000 K and depends on the chamber's density.**
    - 117 kg is the `near_term` run with its 4.4 kg polyethylene plug per rod: 108 kg methane
      plus 8.8 kg plug. With the §1 frozen-methane plug, so that only the rod is polyethylene, it
      is 115 kg. The pulse energy fixes the total either way.
    - The paper's 20 m^3 chamber needs **~146 kg** for the same 5 kg pulse (144 kg with a methane
      plug), at ~496 bar. At ~0.55 kg/m^3 more of the gas is dissociated and ionized, so each
      kilogram holds 115-120 MJ at 7000 K, against 93 MJ at 3.8 kg/m^3. Both figures are 7000 K.
  - (Earlier analysis.) **Use two chambers on the departing stack.** The departure burn's finite-burn loss is
    +15-18% Δv at 1 Hz, +6-8% at 2 Hz and +2-3% at 4 Hz on a 500-600 t stack. Two chambers at
    2 Hz recover the 4 Hz figure, and each still has 0.5 s to empty.
  - **One 5 kg chamber matches two 2.5 kg chambers on wall mass** (52.9 against 54.6 t,
    dry wrap, §3e), with 21% less area and one set of mechanisms. It gives up redundancy and
    needs a ~7.6 m diameter. Ask: which does the mission need?
  - Recheck "on the 500-600 t stacks that depart, one chamber is best" with the finite-burn
    loss charged. The ledger does not appear to include it.
    - The answer turns on the wall: a second chamber costs ~35-50 t all steel (about a wash),
      or ~27-29 t as the solved dry wrap with its extension (~5% of the stack, roughly a wash
      against the 4-5% of Δv saved).
    - One chamber at 5 kg per rod is the alternative: the same total wall, ~21% less wall area,
      no redundancy.
    - The paper's 19 t and 130 t walls are bonded overwraps, which this study finds delaminate.
      Scale from the survivable walls.
  - List the refill sequence at the chosen rate as an open engineering item: plug dwell in the
    residue, mist timing, membrane reseal.
  - **Reconcile the departing mass per chamber** (resolved here as 500-600 t; the paper's
    630-780 pulses should become ~2,000-2,650 five-kilogram pulses). 630-780 pulses (159-194 s at 4 Hz) is a
    ~100 t craft, while the growth ledger departs 500-600 t stacks (~4,000-4,300 pulses). The
    finite-burn loss, the chamber count and the fatigue budget all depend on which.
- **Confirmations the paper can cite.**
  - Unthermalized gas escaping early costs only ~0.4% of the pulse's ideal momentum (§3g), so
    assuming full thermalization is sound.
  - A 10,000 K methane chamber gives +11-13% Isp but -13-17% impulse per rod and 1.6-3x the
    wall heat. Counting pitch it is no better per kilogram; keep 7000 K.
  - The liquid share does not soften the wall shocks, so state that the "softens the local
    shocks" role does not hold. The drops equilibrate in microseconds and pay the gas's toll;
    the liquid serves better as plug mass.
- **The tests that gate the dry-wrap chamber.** No full-scale article is needed.
  - Fling strain `v / sqrt(E/rho)`, with `v` the impulse per area over the mass per area, does
    not depend on size.
  - The blast obeys cube-root scaling: at 1/20 scale, 1/8,000 of the energy (~0.9 MJ) gives the
    same pressures over 1/20 the time.
  - Separate-effects first:
    1. Kevlar strand and ring tension fatigue at 1-2% strain for 5,000+ cycles.
    2. A fretting rig for dry layers slapping together.
    3. Pressure cycling a subscale autofrettaged liner + dry-wrap vessel (COPV practice).
    4. **Ring fling:** a 10-30 cm liner + wrap ring impulse-loaded from inside
       (electromagnetic ring expansion, exploding wire or sheet explosive).
    5. Gas-gun flyer-plate spall of the thin liner at µs pulse lengths.
    6. **Hydrogen uptake** through pitch, with and without the alumina film, in a
       reflected-shock tube or pulsed arcjet (hot dissociated gas at hundreds of bar for ms);
       thermal desorption and crack growth afterwards.
  - Last, to validate the blast solution rather than the materials: a subscale chamber with
    a line-blast source (detonating cord or an exploding wire) standing in for the stopped
    rod, with wall pressure gauges. The 75 km/s rod itself is beyond lab guns (~7-10 km/s), but
    the wall sees only the energy per length and where it is released.
- **Hydrogen at the 2.5 kg rod.** Unresolved.
  - The centred-blast spike (~2.95 GPa) is at steel's spall strength.
  - The hydrogen chamber has not been run with the rod and plug as material, the bulge or the
    dry wrap.
  - Say so. The levers are the same bulge and wrap, a larger chamber, or a smaller rod per
    pulse.

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
