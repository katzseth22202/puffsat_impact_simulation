# Argon spray plate: running summary for the companion repos

Started 2026-10-05 in a grill with Seth. Decisions are in
[ADR-0055](adr/0055-argon-arm-of-the-spray-plate.md), and terms are in `CONTEXT.md` under
"spray plate". This file is written so it can be copied into `Balloon-Pulse-Propulsion` (parent)
and `puffsats_for_datacenters` (child). Update it as the argon arm runs.

**Status of every number here: provisional.** They come from a 1-D planar ideal-gas Lagrangian
estimate with no ionization, radiation or rim spill (`todos/argon_plate/mix1d.py`, not
committed), or from back-of-envelope algebra. The exception is the step 1 table, which is a
repository solve (`make water-plate-argon-radiation`).

## The questions

1. Does the PuffSat gas heat up and mix with the argon spray, or does it bounce back?
2. With 100 t of steel plate at `k = 10`, how much PuffSat and argon can one pulse carry?
3. How close to the overtake ceiling `m w (1 + √(1 + k))` can ideal mixing get?
4. The same questions for water, for comparison.

## Answers so far

**Heating and mixing.** The PuffSat gas cannot stream through the argon. The argon column is about
10^6 collision lengths thick, and fine liquid-argon droplets break up and boil in nanoseconds to
microseconds. Bounce-back is a question of layout. A thin, dense pooled layer acts as a wall
(**cushion bounce**). A diffuse **spray cloud** does not.

**Do the fluid models hold at these speeds? Yes, by about six decades (R3, 2026-10-06;
`make water-plate-collisionality`).** At the 150 t design point (k = 8.52, 12 MN s, 4 m cloud on
a 5 m footprint) the argon is 2.5-3.8e25 atoms/m^3. Cold neutral water meets cold neutral argon at
46-68 km/s, so O-Ar collides at 123-274 eV and H-Ar at 11-24 eV in the centre-of-mass frame.
Momentum-transfer cross-sections from the ZBL screened-Coulomb potential (Ziegler, Biersack and
Littmark 1985), integrated exactly, are 1.3-2.1e-20 m^2 for O-Ar and 2.1-3.4e-20 for H-Ar. The
hard-sphere value at that speed (`pi r0^2`, the generous case) is 1.0-1.6e-20.

| | 45.58 km/s | 65.13 km/s | 68 km/s |
|---|---|---|---|
| O stopping length in the spray | 1.8 µm | 4.0 µm | 4.4 µm |
| Kn against 4 m, generous neutral | 4.5e-7 | 1.0e-6 | 1.1e-6 |
| Momentum-relaxation lengths in the argon column | 2.2e6 | 1.0e6 | 0.9e6 |
| Ion-neutral (Ar+-Ar reference, 1e-18 m^2) | 0.03 µm | 0.04 µm | 0.04 µm |
| Coulomb drag, O+ in Ar+ at Te 3 eV | 3.7 µm | 19 µm | 23 µm |

- **Coulomb drag is not shorter.** Its cross-section falls as `v^-4`, so at these speeds a fully
  ionized spray stops an ion over a *longer* path than the neutrals do. That is still five
  decades under the cloud depth. The earlier line here, that Coulomb collisions take over with
  larger cross sections, was wrong at 50-70 km/s.
- **Leading edge.** Kn reaches 0.1 against 4 m at a PuffSat density of 1.5-2.6e20 m^-3, which is
  (1.6-3.9)e-5 of the pulse's mean. The mass arriving below it is 2.5e-6 to 6.7e-6 for a 1-D
  Gaussian, 1.6e-5 to 3.9e-5 for an exponential tail and 0.6e-4 to 1.5e-4 for a 3-D Gaussian. That
  is its *own* collisionality, though. A thin leading-edge atom still meets the spray at full
  density and stops within microns. It is deposited, not streamed through.
- **Droplets are the one escape.** If the argon is still liquid when the gas arrives, an atom
  travels `4 r rho_l / (3 rho_spray)` between droplets. That reaches Kn = 0.1 at a droplet radius of
  0.36-0.54 mm. Sub-0.1 mm spray is safe. Millimetre drops that have not boiled are not.
- **PLX benchmark.** Argon plasma jets show collisional stagnation at stopping length / layer
  width 0.3-0.6 (Merritt et al., Phys. Plasmas 21, 055703 (2014), oblique, ~40 km/s). They
  interpenetrate at 40 and turn to stagnation by about 5 (Moser and Hsu, Phys. Plasmas 22, 055707
  (2015), head-on, v_rel ≈ 90 km/s). This study sits at ~1e-6, about six decades on the
  collisional side of a measured boundary. PLX checks the criterion, not our cross-sections:
  its collisions are Coulomb at ~1e20 m^-3, and ours are neutral at ~3e25 m^-3.
- Spreading the cloud deeper does not thin the column, which `k` fixes. Only spray thrown outside
  the footprint, or spray still in large drops, can leak.

**Share of the overtake ceiling at `k = 10`** (ceiling 4.32 `m w`; a bare elastic bounce is 2, or
0.46). The stratified and pearl rows are extrapolated from 10/20/40 cells per layer. The
premixed row converges directly.

| Spray layout | Argon | Water (γ = 1.25 stand-in, before chemistry) |
|---|---|---|
| Thin dense pooled layer | ~0.57 | ~0.55 |
| One unmixed layer, as deep as the pulse | ~0.79 | ~0.73 |
| 16 interleaved sub-puffs (axial) | ~0.85 | ~0.80 |
| Fully premixed (ideal mixing) | **0.87** (3.75 `m w`, η_jet ≈ 0.83) | 0.81 (3.50 `m w`) |

The parent's requirement is η_jet = 0.775, or 3.57 `m w` (0.83). Ideal argon mixing clears it in
1-D with about 5% to spare, before any rim spill or radiation loss. Even perfect mixing stops near
0.87, because gas expanding off a wall leaves at a spread of speeds. Interleaving axial sub-puffs
was the first model of the "string of pearls". The real design is a **necklace charge**: a ring of
spheres shaped explosively to the ideal mixing shape. The premixed row is the 1-D bound that the
necklace tries to reach.

**Water's chemistry** comes on top of that. If the bonds never re-form, the child repo's 0.844
toll applies, which takes water to about η_jet 0.64 (3.1 `m w`). This repository's finite-rate run
returned 23-31% of the bond energy by 250 µs. That is preliminary. The study will count only the
bonds that re-form while the gas still presses on the plate.

**Step 1 solved (2026-10-05): the premixed pulse with real EOS and radiation.**
`make water-plate-argon-radiation` runs the merged pulse (the necklace's target state) into the
plate. It uses the repository's 1-D Lagrangian kernel with flux-limited radiation diffusion. Argon
uses a new atomic H/O/Ar Saha table (NIST ionization energies) with TOPS gray opacities for the
`k = 10` mixture. Water uses the equilibrium water table with TOPS water opacities. The footprint
is 78.5 m^2, with 47 kg of PuffSat at 45.58 km/s and 33 kg at 65.13. The spray-cloud depth is
swept from 0.5 to 8 m. All 20 cases converged.

| | Argon 45.58 | Argon 65.13 | Water 45.58 | Water 65.13 |
|---|---|---|---|---|
| Starting temperature | 31-40 kK | 47-61 kK | 13.6-14.9 kK | 22-30 kK |
| Share of overtake ceiling | 0.69-0.74 | 0.69-0.74 | 0.67-0.68 | 0.72-0.75 |
| Radiation onto the face, share of pulse KE | **1.7-11%** | 1.5-8.7% | 0.17-0.74% | 0.42-2.7% |
| Peak face pressure (0.5 m -> 8 m) | 763 -> 42 MPa | 1094 -> 60 MPa | 744 -> 43 MPa | 1148 -> 62 MPa |

What it changes:

- **The heat gate fails, and argon fails it by 10x more than water.** Steel plus a 5 µm film
  take about 0.07%. Argon at 1 m depth puts 2.6% on the face, about 1.3 GJ per pulse. That is
  roughly 130 kg of film per pulse at 10 MJ/kg, about three times the PuffSat's mass. Water at
  0.24% needs about 12 kg. Argon runs hotter because it has fewer particles per kilogram, and
  radiation grows as T^4.
- **With real chemistry and radiation, argon's impulse advantage is small or gone.** Both sit near
  0.7 of the ceiling, which is η_jet ≈ 0.58-0.68, below the parent's 0.775 requirement before any
  rim spill. The 1-D ideal-gas bound of 0.87 was optimistic. The comparison also favours water:
  argon here pays water's bond energy permanently, while the equilibrium water table returns it.
- **Depth trade.** Deeper clouds radiate more and press less. 0.5 m exceeds the steel ladder at
  65 km/s, so the working depth is about 1-2 m.
- **Caveats.** 1-D planar, premixed, gray flux-limited diffusion. The face is a cold black
  absorber with no vapor shielding. The repository's ablating-wall solver (Rung E) models the
  film's vapor absorbing radiation before it reaches the steel, and is the next run.

**Step 1b solved (2026-10-05): the same pulse on an ablating, vapor-shielded face.**
`make water-plate-argon-ablating` runs the Rung E ablating wall: the film boils at
`q_in / Q*`, and its vapor forms a curtain of optical depth `κ_vapor × ablated mass`. Swept over Q*
2/5/10 MJ/kg (ADR-0014) and κ_vapor 0 / 200 (the Rung E calibration) / 1,500 (the child repo's
20%-carbon film at Bond & Bergstrom's 7.5 m^2/g, 550 nm) / 10^4 (an EUV photoabsorption estimate,
not sourced), at depths of 1 and 2 m. All 96 cases converged.

Film per pulse at 1 m depth, Q* = 5 MJ/kg:

| | unshielded | κ 200 | κ 1,500 | κ 10^4 |
|---|---|---|---|---|
| Argon 45.58 km/s | 124 kg | 27 kg | 13 kg | 6 kg |
| Argon 65.13 km/s | 118 kg | 33 kg | 18 kg | 9 kg |
| Water 45.58 km/s | 14 kg | 4.7 kg | 2.2 kg | 0.9 kg |
| Water 65.13 km/s | 27 kg | 10 kg | 5.3 kg | 2.5 kg |

- **Shielding makes argon affordable.** At κ 1,500-10^4, argon costs 6-18 kg of film per pulse,
  1-4% of its 470 kg of spray, against 1-5 kg for water. The film must be about 75-230 µm thick for
  argon and 10-70 µm for water. The child repo's 5 µm is too thin, and Orion's 150 µm was about
  right.
- **Argon's advantage returns at the cold end.** On the ablating wall argon reaches 0.74-0.75 of
  the ceiling at 45.58 km/s against water's 0.68, about 10% ahead. At 65.13 km/s they tie at
  about 0.745. Water is still favoured, because its table returns all of its bond energy.
- Boiled film joins the exhaust, which raises the impulse slightly (0.77-0.80 unshielded).
- **Still premixed.** Poor mixing leaves unmixed PuffSat gas stagnating near 1 GJ/kg, far hotter
  than the merged ~85 MJ/kg. Those pockets radiate as T^4, partly screened from the face by the
  argon. That is step 2, and it needs the 1-D kernel's radiation to handle two materials.

**PuffSat droplets (estimate).** Taking a 1 m cloud (ρ ≈ 6.6 kg/m^3, ~6 kg/m^2 of argon):
- An unbroken drop of radius `r` stops after ~2.7 ρ_water r of column, so up to ~2 mm stop inside
  the cloud.
- Aerodynamic breakup comes first, after ~10 r √(ρ_water ρ_gas) ≈ 810 r, which shatters drops up
  to ~7 mm.
- A stopped fragment needs ~2.5 MJ/kg to boil, against ~85 MJ/kg around it.

The requirement is that **the largest surviving fragment has a radius of at most 1 mm**, which
stops it with 2x margin even without breakup. The parent's foam grid (~25 mm^2 of water per strand)
caps fragments near 2.5 mm radius. That passes with ~2.8x margin if breakup is credited and misses
by ~1.1x if not. Carried to the parent in `docs/argon_plate_owed_to_companions.md` (P1).
Droplets that stop deep in the argon also help mixing, because they deposit there.

**The necklace is the parent's ring plus hub bag** (`templateArxiv.tex:533-547`). The parent notes
that filling a ring's center takes radial expansion of order the ring radius, and the hub fills
it. So the delivered shape is a disk, prescribed with its areal mass matched to the tapered plate.

**Step 2a solved (2026-10-05): mixing quality, two materials with radiation.**
`make water-plate-argon-mixing` puts the PuffSat's water gas and the cold spray in separate
Lagrangian layers, ordered from the plate outward: spray, PuffSat, and so on. It runs 1, 4 or 16
layer pairs on the new per-cell-material kernel (pure argon table with TOPS argon opacities). Each
layout runs on a bare face and on ablating faces at κ_vapor 5,000 and 100. Cells never exchange
mass, so finer interleaving is the only mixing represented. All 72 cases converged.
`--spray-mixing-convergence` repeats the extreme cases at 192-1,536 cells.

| Layered, k = 10, 1-2 m | Argon | Water |
|---|---|---|
| Share of overtake ceiling (premixed 0.69-0.75) | **0.64-0.67** | 0.59-0.65 |
| Change from 1 to 16 layer pairs | +0.01 at most | up to +0.04 |
| Face radiation, bare (premixed 2.6-4.2%) | 1.0-2.3% | 0.05-0.4% |
| Film per pulse, κ 5,000 | 5.5-12.5 kg | 0.4-2.9 kg |
| Film per pulse, κ 100 (soot lost) | 26-48 kg | 2-13 kg |
| Converged peak face pressure, stratified, 1 m, 45.58 km/s | **~4.7 GPa** | **~7 GPa** |

- **Poor mixing costs 0.05-0.10 of the ceiling,** converged to ±0.005. Interleaving hardly helps
  in 1-D with real chemistry and radiation, unlike the ideal-gas estimate.
- **Peak pressure becomes the binding limit, ahead of heat.** An unmerged PuffSat layer drives a
  shock through the spray and into the steel, at 5-10x the 400-900 MPa ladder. That is about 3x
  the ~1.5 GPa the same pulse would put on a bare plate (1.2 ρ w^2), because the heavy spray layer
  compresses and reflects. Interleaving only appears to help on a coarse grid: at 16 pairs the
  peak climbs 0.9 -> 3.7 GPa from 192 to 1,536 cells. The premixed solve (0.17-0.53 GPa) hides
  this, because there the merge has already happened away from the plate.
- **Face radiation is overestimated at production resolution.** It falls ~12% per grid doubling
  (argon stratified 1.5% -> 0.74% from 192 to 1,536 cells). The heat and film numbers in steps 1,
  1b and 2a are therefore conservative by roughly 1.5-2x.
- **The vapor opacity matters 4-5x for film.** Argon needs 5.5-12.5 kg per pulse with sooty vapor,
  but 26-48 kg if the curtain loses its soot. Water stays at 0.4-13 kg either way.

**Step 2a' solved (2026-10-05): keeping an unmerged pulse off the plate.**
`make water-plate-argon-levers` takes step 2a's worst case (stratified, bare face, 768 cells) and
sweeps the spray-cloud depth and the arriving pulse length independently at 1, 4 and 16 m. The
column masses are fixed by `k`, so a deeper cloud is more dilute and adds distance between contact
and plate. The 1 m / 1 m row reproduces step 2a's 4.627 GPa exactly.

Peak face pressure, argon / water, at 45.58 km/s:

| | 1 m cloud | 4 m cloud | 16 m cloud |
|---|---|---|---|
| 1 m pulse | 4.6 / 6.8 GPa | 1.3 / 2.1 | **0.32 / 0.54** |
| 4 m pulse | 2.7 / 5.4 | 1.2 / 1.8 | **0.33 / 0.55** |
| 16 m pulse | 0.46 / 1.0 | 0.72 / 1.4 | 0.31 / 0.46 |

- **A deep spray cloud is the lever that works.** The PuffSat's shock decays crossing it. At
  16 m, the converged peaks (384-1,536 cells, rising 2-4% per doubling) extrapolate to ~0.34 /
  0.52 GPa for argon and ~0.57 / 0.77 GPa for water, at 45.58 / 65.13 km/s. Argon sits inside the
  700 MPa rung. Water is marginal at the fast end.
- **A long pulse alone is a trap.** On a 1 m cloud, a 16 m pulse lowers the peak but drops the
  impulse to ~0.47 of the ceiling, and radiates up to 18% of the pulse's energy to space.
- **The cost of depth.** The impulse falls from ~0.65 to ~0.58-0.60 of the ceiling. Argon's face
  radiation rises to 3.8-5.2%, and unlike the 1 m case it does not fall with resolution. Water's
  stays at ~0.1-0.4%.

**Design point: 16 m cloud, 4 m pulse,** on shielded faces
(`--spray-levers-shielded`):

| | Argon | Water |
|---|---|---|
| Share of overtake ceiling | 0.616-0.626 | 0.575-0.607 |
| η_jet equivalent | ~0.50 | ~0.45 |
| Film per pulse, κ 5,000 | 10-15 kg (2-3% of spray) | 1.1-3.3 kg |
| Film per pulse, κ 100 | 51-61 kg (11-13% of spray) | 5.8-14 kg |
| Peak face pressure | 0.35-0.64 GPa | 0.56-0.77 GPa |

Argon leads water by 2-7% in impulse and has more pressure margin, at 3-10x the film. Both are
well below the parent's η_jet = 0.775 in every survivable 1-D configuration, before rim spill. The
premixed solve's ~0.65 is the ceiling if the merge could be made to happen upstream of the plate.
In 2-D a 16 m column over a 10 m footprint will spill sideways, and lateral relief will also lower
the peak, so step 2b must sweep cloud depth rather than inherit 16 m.

**Step 2b solved (2026-10-05): rim spill and lateral relief in 2-D.**
`make water-plate-argon-2d` uses the ADR-0003/0008 geometry method: the resting spray cloud and the
arriving 4 m slug as one effective-γ gas (γ 1.4, Mach 20), 5 m footprint radius on the 10 m plate.
Each depth runs free and confined. `eta_capture` = free/confined impulse, and the relief factor =
free/confined peak facesheet pressure. Both multiply the 1-D real-physics column at the same depth.
At 2x resolution `eta_capture` moves 1.5-2.5%; relief at 16 m moves 0.36 -> 0.40.

| Spray cloud depth | 4 m | 8 m | 12 m | 16 m |
|---|---|---|---|---|
| `eta_capture` (impulse kept on the 10 m plate) | **0.79** | **0.61** | **0.44** | **0.30** |
| Pressure relief factor | 1.00 | 1.00 | 0.68 | 0.36-0.40 |
| Argon delivered share (η_jet) | 0.50 (0.35) | 0.38 (0.19) | 0.27 (0.05) | 0.18 (< 0) |
| Argon peak, 45.58 / 65.13 km/s | 1.2 / 1.7 GPa | 0.64 / 0.94 | 0.29 / 0.44 | 0.12 / 0.18 |
| Water peak, 45.58 / 65.13 km/s | 1.8 / 2.7 GPa | 1.0 / 1.4 | 0.48 / 0.63 | 0.20 / 0.26 |

- **The deep cloud that fixed pressure in 1-D throws the impulse away in 2-D.** A 20 m column over
  a 10 m-wide footprint spills sideways past the plate. At 16 m the plate gets less than the
  PuffSat's own momentum.
- **Shallow clouds keep the impulse but get no pressure relief.** The peak lands on axis before
  any relief from the cloud's edge arrives (relief factor exactly 1 at 4 and 8 m).
- **No depth gives both in this geometry.** The best survivable case is argon at 8 m and
  45.58 km/s: 0.64 GPa at η_jet ~0.2. Water fails pressure at every depth that keeps useful
  impulse.
- **Caveats.** Effective-γ geometry with no radiation or real EOS in 2-D; flat plate; fixed 5 m
  footprint radius and 10 m plate; unmerged (stratified) pulse. A wider footprint, a dished plate,
  or merging the pulse before it reaches the plate could each change this, and are untested.

**Merge before the plate, 1-D attempt (2026-10-05): the method failed; spike widths recovered.**
`--spray-standoff` placed the spray cloud 0.5-4 m off the plate, with a near-empty gap
(2e-4 kg/m^3) between them. **Every gap > 0 case failed numerically:** peaks of 10^11-10^12 Pa,
near-zero impulse, and most never converged. The arriving shock crushes the near-empty gap cells
against the plate far beyond the table, and the time step collapses. A Lagrangian grid cannot
carry a near-vacuum gap, so the standoff needs the 2-D Eulerian kernel, whose ambient gas
represents the gap natively. Those rows are discarded.

The gap = 0 rows are valid (they reproduce step 2a'), and give the **spike width** (Q29), the time
the face spends above each allowable:

| Unmerged pulse, no gap | Peak (45.58 / 65.13 km/s) | Time > 0.7 GPa | Time > 1.5 GPa |
|---|---|---|---|
| Argon, 1 m cloud | 2.7 / 4.6 GPa | 28 / 24 µs | 16 / 14 µs |
| Argon, 2 m cloud | 2.2 / 3.2 GPa | 26 / 27 µs | 4 / 9 µs |
| Argon, 4 m cloud | 1.2 / 1.7 GPa | 13 / 19 µs | 0 / 1 µs |
| Water, 4 m cloud | 1.8 / 2.7 GPa | 14 / 17 µs | 0 / 4 µs |

The overloads last 1-30 µs, a few steel stress-wave transit times (~7 µs). On maraging's
2.5 GPa allowable, every 4 m cloud passes except water at 65 km/s (2.7 GPa for a few µs), and the
1-2 m clouds exceed it for ~5-16 µs a pulse.

**Merge before the plate, in 2-D (2026-10-05).** `--spray-2d-standoff` lifts the 4 m cloud off
the 10 m plate by a gap carried as ambient gas (the Eulerian kernel holds it natively). Each gap
runs free and confined, compared to the confined no-gap run, which is the 2-D twin of the valid
1-D column. Those ratios scale step 2a's real-physics 4 m numbers. Gap 0 reproduces step 2b
exactly.

| Gap | 0 m | 1 m | 2 m |
|---|---|---|---|
| Impulse vs 1-D column | 0.790 | 0.768 | 0.727 (2x grid 0.721) |
| Peak vs 1-D column | 1.00 | 0.76 | 0.55 (2x grid **0.37**) |
| Argon share (η_jet); peak 45.58 / 65.13 | 0.50 (0.35); 1.2 / 1.7 GPa | 0.48 (0.33); 0.90 / 1.3 | 0.46 (0.29); 0.65 / 0.96 |
| Water peak 45.58 / 65.13 | 1.8 / 2.7 GPa | 1.4 / 2.1 | 1.0 / 1.5 |

- **A 1-2 m gap cuts the peak 25-45%** (more on the finer grid) for 3-8% less impulse. A 1 m gap
  brings water at 65 km/s under maraging's 2.5 GPa.
- **It does not recover the premixed ~3 mw.** Confined, the gap raises the impulse 7-12% (a
  partial-merge gain), but spill past the rim takes more. The delivered impulse stays near
  η_jet 0.3, ~2.1 mw against the bare plate's 1.63.
- **The 4 m and 8 m gap rows are not trusted.** Their confined runs return 4.4-4.7x the no-gap
  impulse, above the overtake ceiling, so something in them is wrong and their free rows are
  suspect too. Under investigation.
- **The peak relief is not grid-converged.** It moves 0.55 -> 0.37 at 2 m. The impulse is.

**Containing the spill (2026-10-05): dish, skirt and flared wall in 2-D.** `--spray-2d-cup` puts
the 4 m cloud 1 m off the 10 m plate, shaped as a dish, a cup with a straight skirt, or a cup with
a flared wall. The kernel's new `PlateShape::Cup` puts a wall on the dish rim, and the axial force
sums every vertical gas-solid contact. Each shape runs free and is compared to the confined flat
no-gap run, and the ratios scale the 1-D real-physics 4 m column.

| Shape | Impulse vs 1-D | Argon η_jet (45.58) | Water peak (65.13) | Added steel |
|---|---|---|---|---|
| open plate | 0.77 | 0.33 | 2.1 GPa | none |
| dish 0.10 / 0.15 | 0.92 / 0.97 | 0.45 / 0.49 | 2.5 / 2.3 | none |
| straight skirt 2 / 4 / 8 m | 0.96 / 1.02 / 1.04 | 0.48 / 0.53 / 0.55 | 2.1 | 9.9 / 19.7 / 39.5 t per cm |
| flare 0.3 / 0.6, 4 m | 1.01 / 0.98 | 0.53 / 0.51 | 2.1 | 21.8 / 25.8 t per cm |
| **dish 0.10 + skirt 4 m** | **1.07** (2x grid 1.06) | **0.58** | 2.5 (2x grid **1.8**) | 19.7 t per cm |
| dish 0.10 + flare 0.3, 4 m | 1.07 | 0.57 | 2.5 | 21.8 t per cm |

- **Containment doubles η_jet, 0.33 -> ~0.58.** Spill was the largest loss, and the walls recover
  it. With the gap, a contained column now beats the 1-D column (x1.07), because the merge
  happens before the plate.
- **The dish is the best value:** η_jet 0.49 at d/D 0.15 for no added steel.
- **A straight skirt beats a flare.** Keeping the gas in matters, and turning it does not.
- **The dish focuses the peak ~18% higher,** but its relief tightens on the finer grid
  (0.90 -> 0.67). Water at 65 km/s sits near 1.8 GPa, under maraging's 2.5 GPa.
- **A 4 m skirt costs ~20 t per cm of wall thickness,** a real share of the 100 t. The wall sees
  far lower pressure than the face, so it can be thinner. It needs structural sizing, and film and
  cooling of its own.
- The paper's η_jet = 0.775 is still not met. The remaining levers are the injection ratio k,
  premixing, and the radiation loss.

**Film choice and the injection ratio (2026-10-05).**

*Pitch against oil* (`--spray-film`, then the thermal cycles):
- **Film use.** Shielded at κ_vapor 5e3, pitch burns 1.8-2.8 kg per pulse (argon) and 0.2-0.7 kg
  (water). Carbon-loaded oil (1.6 MJ/kg; `docs/spray_plate_film_properties.md`) burns 12-18 kg
  and 1-4 kg.
- **Steel heat.** Standing thickness, not consumption, insulates the steel. Pitch kept ~150 µm
  thick holds the steel at ~507 K over 1,500 pulses for both arms, with no cooling water. Oil
  runs it at ~300-490 K, but costs 4-7x the film. Pitch is kept (ADR-0055 decision 11).
- **The thermal solver has no vapor shielding,** so its film-removed figures do not apply. The
  shielded run's do.

*Injection ratio k* (`--spray-k`). The cup (dish 0.10 + 4 m skirt, 1 m gap) at k = 4, 6, 8.5, 10
and 14, both mixing bounds (stratified layers and premixed slab, on each k's mixture table):

| Argon, stratified -> premixed | k = 4 | 6 | 8.5 | 10 | 14 |
|---|---|---|---|---|---|
| η_jet | 0.60 -> 0.69 | 0.59 -> 0.68 | 0.58 -> 0.68 | 0.58 -> 0.68 | 0.54-0.57 -> 0.68-0.69 |
| β = J/(m w) | 2.33 -> 2.53 | 2.57 -> 2.80 | 2.80 -> 3.10 | 2.91 -> 3.27 | 3.1-3.2 -> 3.6-3.7 |
| J per kg consumed (x w) | 0.46 -> 0.50 | 0.36 -> 0.40 | 0.29 -> 0.32 | 0.26 -> 0.30 | 0.21 -> 0.24 |

- **η_jet is nearly flat in k.** Argon holds ~0.58 unmixed and ~0.68 mixed; water's unmixed value
  slides 0.63 -> 0.49. Even premixed, nothing reaches the paper's 0.775. **Mixing is worth ~0.1 in
  η_jet,** which is why the plug branch, which forces the merge, matters.
- **k is a mission trade,** between impulse per scarce PuffSat and impulse per carried kilogram.
  The parent's launch-mass budget should set it (ADR-0055 decision 12).
- Every stratified peak stays at or under maraging's 2.5 GPa (water k = 10, 65 km/s: 2.46).
  Premixed peaks are 0.06-0.16 GPa.

**An inward lip on the skirt (2026-10-06): it lowers the impulse.** `--spray-2d-lip` leans the
working cup's 4 m skirt inward (flare 0 to -1.0, mouth 20 -> 12 m), with its own stepping loop for
wall pressure. The straight case reproduces step 2d (x1.071).

| Skirt lean (mouth) | 0 (20 m) | -0.2 (18.4 m) | -0.4 (16.8 m) | -0.6 (15.2 m) | -1.0 (12 m) |
|---|---|---|---|---|---|
| Impulse vs 1-D | 1.071 | 0.973 | 0.898 | 0.812 | 0.770 |
| Argon η_jet | 0.577 | 0.497 | 0.435 | 0.364 | 0.330 |

- **A lip is the converging half of a nozzle without the diverging half.** Gas pressing up on the
  lip's underside acts against the thrust, and the jet leaves near sonic with pressure still in it,
  then expands outside where nothing catches it. For γ = 1.4 a converging-only nozzle gives ~70% of
  the ideal thrust coefficient, below the open cup's ~87%, because free expansion off the floor
  already turns nearly all the heat into motion.
- **The straight-skirt cup stays the best open design.** A lip pays only with a diverging bell
  after it, which is the heavy chamber end. Even then a pulsed chamber's decaying pressure spreads
  its exit speed, so the ~13% speed-spread loss is close to fundamental for pulsed operation.
- **The first lip run stopped at once:** the ambient gas under an overhang makes the starting net
  force negative, and the kernel's tail guard then fires. The lip loop now arms its guard only
  after the pulse arrives. `run_bounce` in `euler2d` has the same latent bug for any overhanging
  cup, which needs a test-first fix. Step 2d's outward flares were not affected (the bug gives
  near-zero impulse, which they did not show).
- The "time in cup" measure did not discriminate (~19 crossings either way), because the spray's
  own mass dominates it.

**Steel temperature over a push (2026-10-05): ablative film on maraging, propellant-cooled.**
`make water-plate-plate-thermal` drives the coated-wall solver (`walled_nozzle/wall_layers.py`,
extended with a substrate, between-pulse cooling and a cold respray) for 1,500 pulses at 4 Hz.
Each pulse's radiant energy comes from the solved 4 m cloud, 4 m pulse rows (argon 11.1-12.9
MJ/m^2, water 0.5-2.5 MJ/m^2). The film is pitch as a stand-in for the carbon-loaded oil:
150 µm for argon, 30 µm for water. The substrate is maraging 300 with a 753 K aging limit. The
pulse lasts 1 ms (bracketed 0.3-3 ms).

| Cycle | Argon steel peak | Water steel peak |
|---|---|---|
| no cooling | 1,295 K (fails) | 1,346-1,348 K (fails) |
| cold film respray only | **507 K** | 991-993 K (fails) |
| + water spray, h 300 (floating) | 488 K, 0.3 kg/pulse | 915-916 K (fails) |
| + water spray, h 3,000 (wetting) | 418 K, 1.1 kg | **678-679 K, 1.5 kg** |
| + water spray, h 30,000 | 375 K, 1.6 kg | 552-553 K, 1.9-2.0 kg |

- **Cooling is mandatory.** The uncooled steel ratchets to ~1,290 K.
- **Argon: the cold film respray alone is enough.** No water is needed. The pulse-length bracket
  gives 408-664 K, still under the limit.
- **Water needs ~1.5-2 kg of wetting spray per pulse,** despite 5-23x less radiation on the face.
  The steel's heat is set by conduction through the film from its ablating surface, so it tracks
  the film's thickness, not the radiant load.
- **Caveat: the film is nearly consumed every pulse,** ~95% for argon (14 of ~15 kg) and all of
  it for water at 65 km/s (3.0 of 3.1 kg). The solver keeps the film at full thickness as it
  ablates (ADR-0053's fixed grid), so it over-credits the insulation. The film should be specified
  at no less than ~2x its consumption, and the run repeated at that thickness.
- The water spray is a coolant consumable that boils off before the pulse, not propellant. The
  sequence is spray water 10-100 ms after a pulse, then the film at 200 ms, with nothing liquid
  left at 250 ms.

**Pulse size for 100 t of steel.** The gas spring sets this, not steel strength. With a preloaded,
near-constant-force spring, `J ≈ 8 m_plate ν s`. This reproduces the child repo's 3.6 g peak on a
2.5 m stroke. At 4 Hz and 2.5 m, J is 8 MN·s per pulse, enough to carry a vehicle of about 1,090 t
at 3 g.

| Closing speed | PuffSat per pulse | Argon per pulse |
|---|---|---|
| 45.58 km/s | 47-49 kg | 470-490 kg |
| 65.13 km/s | ~33 kg | ~330 kg |

The limit scales with stroke × pulse rate.

## Plate design (ADR-0055)

- **A tapered plate, 20 m in loaded diameter**, about 4 cm mean thickness, thicker at the center.
  A uniform plate fails at any thickness. The struck footprint gets four times the mean kick, and
  steel stores about 100x too little elastic energy to absorb the difference.
- **Aim to 10 cm** is 2% of the footprint radius, against the taper's ~15-20% tolerance. The paper
  should give precise aim as a second reason the taper works.
- **The gas-spring bag's fiber mass is the same at any plate area** (force × stroke, about
  80 MJ). That is ~70 kg of Vectran ideal, or ~300 kg with margin. Bag pressure, about 1 atm at
  20 m, can be raised freely, and a stiffer gas foundation damps the plate's low bending modes.
- **Face pressure is not the binding limit**: about 51 MPa mean at 20 m against a 400-900 MPa
  ladder, assuming a 0.5 ms pulse.
- **Heat is the binding limit.** Each pulse carries about 49 GJ, about 200 GW at 4 Hz. 100 t of
  steel absorbs about 25 GJ over a 500 K rise. A 5 µm ablative film absorbs about 16 MJ per pulse,
  assuming a heat of ablation near 10 MJ/kg (not yet sourced). Over the ~1,500 pulses of the Earth
  push, the steel and film together can take only about 0.07% of each pulse's energy. Above that,
  film mass grows with the heat. At 1% it would equal the PuffSat mass. Diameter sets only the
  film minimum needed for the infrared readout: about 1.6 kg at 20 m, 3.5 kg at 30 m.

## Changes owed to the paper and the child repo

The full, copy-ready list is [`argon_plate_owed_to_companions.md`](argon_plate_owed_to_companions.md).
The short form follows.


- `sec:plate_construction` should say the plate is **tapered** to the impulse profile, and why.
  Precise aim is a second reason.
- `sec:plate_liquid_spray` should state the **spray-cloud** requirement and the cushion-bounce
  failure.
- `sec:plate_aiming`: impact offset is capped at about 0.5 m by the taper.
- `sec:plate_ablative_film`: heat sets the film's thickness, not the readout. That is about
  75-230 µm for argon and 10-70 µm for water, with vapor shielding. The 5 µm readout film is too
  thin, and Orion's 150 µm was about right.
- Parent PuffSat design: state the ≤1 mm largest-fragment requirement on the foam-grid spacing.
- Once the argon arm has solved values, replace the provisional η_jet ≈ 0.83 bound.

- **Argon's case in `sec:plate_liquid_spray` needs a heat caveat.** The premixed solve puts about
  10x more of the pulse onto the face for argon than for water. The 8-11% faster doubling the parent
  credits argon with does not yet pay for that film.

## Open

- **The plug alternative (agreed 2026-10-05, after the k sweep and film comparison).** Each
  necklace sphere arrives compact and buries itself in a **plug**: dense, narrow carried mass on
  tethers beyond the blast radius, steered in the last ~100 ms to millimetre accuracy. This is the
  parent's head-on rod-and-plug idea without the sliding door. The plug forces the merge that the
  spray leaves to chance. Estimates to check first:
  - A ~10 cm water sphere has ~140 kg/m^2 of areal mass, so its 47 kg plug wants ~400+ kg/m^2 of
    column: ~0.4 m across and ~0.4 m long as solid polyethylene (~14 m long as 30 kg/m^3 foam).
  - The impact makes a near-isotropic fireball (~10+ km/s expansion against ~4 km/s drift), so
    capture depends on how much the cup encloses. A bell may pay here where it did not for spray.
  - Layout: one plug per sphere (overlapping cones) against one ring plug. Material: polyethylene
    first (with its ~84 MJ/kg bond toll and soot radiation bracketed, an estimate), then argon ice.
  - Reuse the tamped-nozzle study's snowplow and column-density tools, and the 2-D standoff/cup.
  - **Step 1 done (closed form, 2026-10-05).** The tamped-nozzle ballistic fireball, mirrored for
    the overtake: the blob drifts toward the plate at V = w/(1+k) and expands isotropically at
    u = w√k/(1+k), and elements moving toward the plate reflect specularly. Catching all of them,
    `β = (1+k)[u(1-μ0²)/2 + V(1-μ0)]/w` with `μ0 = -1/√k`, which is **β 2.74, η_jet 0.52 at
    k = 10**. That is below the spray cup's 0.58 unmixed. A plug 4 m off the floor gives 0.45, and
    8 m gives 0.27. Sizing: ten 4.7 kg spheres (10 cm radius, ~140 kg/m^2) want plugs 15-19 cm in
    radius: 0.44-0.73 m long as solid polyethylene, 0.26-0.43 m as argon ice, 14-23 m as foam.
  - **Why it trails the spray:** merging at a point makes an isotropic fireball, and ~34% of it
    leaves away from the plate. The spray's broad slab expands mostly along the axis. The plug
    trades the mixing risk for a geometry loss. Two things could recover it: pressure-mediated
    loading beating the ballistic model (the tamped study's hypothesis), and enclosing the plug
    deep in the cup, so its walls turn the fireball toward the mouth (a short chamber).
  - **Step 2, 2-D (2026-10-05; converged to ~2% on a 2x grid).** One 10 cm sphere into a
    k = 10 plug on its 3.2 m plate share (`make water-plate-plug-2d`). The first pass was invalid.
    Its ambient gas (10^-3 of the sphere's density) weighed ~5 times the plug on this wide domain.
    It tamped the fireball and held the plate above the quiet cutoff, so case 0 ran to its step cap
    at β 52, 12 times the ceiling. The plug runs now use an ambient of 10^-6, which reads the same
    at 10^-7 (β 1.361 against 1.357). The spray cup moves 0.16% under the same change, so the
    committed spray numbers stand (`crates/euler2d/tests/plug_ambient.rs`).

    | case | β | of ceiling |
    |---|---|---|
    | plug 1 m off, dish only | 1.361 | 0.315 |
    | plug 1 m off, skirt 4 m | 2.950 | 0.683 |
    | plug 2 m off, skirt 4 m | 3.373 | 0.781 |
    | plug 1 m off, skirt 8 m | 2.958 | 0.685 |
    | plug 2 m off, skirt 8 m | 3.498 | 0.810 |
    | plug 4 m off, skirt 8 m | 3.515 | 0.814 |
    | spray cup (reference) | 3.065 | 0.710 |

    In a cup, a plug 2-4 m off the floor beats the spray cup and the ballistic model (0.635). The
    open dish loses most of the fireball. On a grid twice as fine (the sphere 10 cells in radius)
    case 1 gives β 3.016 (+2.2%), case 4 3.548 (+1.4%) and case 5 3.485 (-0.9%), so the plug in an
    8 m skirt holds at **0.81-0.82 of ceiling**, η_jet ~0.75-0.77 on the effective-γ method, against
    the spray cup's 0.71 (η_jet 0.62) on the same method. The spray's real-physics η_jet (0.58-0.68)
    has not been carried over to the plug; the ratio suggests ~0.70, an estimate.

- **A trailing pearl in formation does not help (2026-10-05, Q54).** A second PuffSat of the same
  architecture trailing the lead by 5-40 m (tail to front), holding 5-30% of its mass, on the
  working cup (`make water-plate-argon-2d-pearl`). Every case trails a single PuffSat that carries
  the same total mass into the same cloud. The loss grows with share and gap: η_jet −0.03 at
  5% and 5 m, −0.18 at 30% and 40 m. The pearl is 0.5-3% of the rebounding mass, too light to
  equalize anything. Its own mass arrives without spray to mix into, and it meets gas moving away
  from the plate head-on, which turns its energy into heat. At long gaps it reaches an emptied cup
  and makes a bare bounce. The spray cloud is already the tamper, in place before the pulse.

- **Bare plates: tampers do not help, a deep bowl does (scratch, 2026-10-05).** Containment
  against the 1-D column, effective-γ 2-D, production grid:

  | plate | no tamper | trailing pearl (5 m) | thin layer at stand-off |
  |---|---|---|---|
  | flat, open | 0.769 | 0.722 (10%), 0.693 (30%) | 0.779 (1 m), 0.726 (4 m) |
  | dish d/D 0.10 | 0.922 | 0.866 (10%), 0.856 (30%); 0.899 (30%, 2 m) | 0.875 (2 m), 0.837 (4 m) |

  A bare plate loses its impulse to sideways spill, and a tamper holds the hot layer down longer,
  so more escapes. Deeper bare dishes, no wall: d/D 0.15 0.967, 0.20 1.015, 0.25 1.060, 0.30
  1.086, against the dish-and-skirt cup's 1.071. A 6 m bowl matches the cup with ~95 m^2 more
  face instead of 251 m^2 of skirt. Next: the wall-load measurement (Q56) prices both.

- **The skirt's real mass (2026-10-05, Q56; `make water-plate-argon-2d-skirt`).** The 2-D cup
  records the gas pressure on the skirt's inner face, row by row, scaled so the floor carries the
  real 8 MN s pulse. The skirt is sized as a ring (a maraging liner wrapped in fibre) whose
  breathing mode is a single oscillator, stepped exactly through each row's pressure history.

  | | dish 0.10 + 4 m skirt | dish 0.10 + 8 m skirt | flat + 4 m skirt |
  |---|---|---|---|
  | wall impulse at the base / mean / top [kPa s] | 19.8 / 12.4 / 5.0 | 19.9 / 9.6 / 4.0 | 24.3 / 15.9 / 6.9 |
  | outward impulse on the whole wall | 3.1 MN s | 4.8 MN s | 4.0 MN s |
  | aramid wrap, 2 mm liner, ε 0.5 / 0.7 / 1.0% | 60 / 43 / 31 t | 87 / 64 / 46 t | 77 / 56 / 40 t |
  | carbon wrap, 2 mm liner, ε 0.5 / 0.7 / 1.0% | 44 / 32 / 23 t | 64 / 47 / 35 t | 57 / 41 / 30 t |

  The floor averages 25.5 kPa s, and the skirt's base sees ~20. The merged layer pushes sideways
  nearly as hard as it pushes down, and the wall must stop that push in hoop strain alone. The
  wall pulse (~1 ms) is short against the ring's ~7 ms period, so the ring takes it as an impulse.
  The wrap is ~12-26 cm thick at the base. **A 4 m skirt costs ~23-45 t**, not a few tonnes. It
  buys +16% impulse per pulse over the bare d/D 0.10 dish. Caveats: effective-γ 2-D scaled to the
  real pulse; every row is a free ring, which ignores the plate rim's support of the lowest rows;
  the allowable is set by the liner's ~1% yield strain (an autofrettaged liner could run further).
  The deep bowl's steep band will see a similar sideways push, so it needs the same measurement.

- **A 70 t floor in a fixed skirt (2026-10-05, Q57-Q58).** The skirt mounts to the vehicle and
  the floor slides inside it as a piston, so the skirt never takes the floor's kick (a rigid joint
  would see a stress wave of ρcΔv, ~1-3 GPa). Thermal (`make water-plate-plate-thermal`, argon,
  150 µm pitch, cold respray): the steel peaks at 507 K on 30 mm, 509 K on 28 mm and 518 K on
  15 mm, all under maraging's 753 K aging limit. The film does the work, so a thinner floor costs
  little. Pulse: the floor alone is the sprung mass, and `J = 8 m ν s` gives 5.6 MN s for 70 t at
  4 Hz and 2.5 m. Keeping 8 MN s needs a 3.6 m stroke (floor at ±57 m/s) or 5.7 Hz. The annular gap
  leaks hot gas (estimate, ~85 MJ/kg merged gas choked through the gap at the skirt-base pressure):
  ~1.6 kg and ~140 MJ per pulse at 10 mm, ~0.16 kg and ~14 MJ at 1 mm. A 20 m floor grows ~22 mm in
  diameter per 100 K, so a tight gap needs a sprung, segmented seal or a labyrinth vented outward
  away from the gas bags.

- **The 150 t plate's shape (2026-10-05, Q59-Q63; `make water-plate-argon-2d-shape`).** 12 MN s
  per pulse at 4 Hz, k = 8.52, the skirt fixed to the vehicle and the floor sliding inside it.
  Each shape's wall histories (skirt rows and the dish's riser rings) are sized for hoop strain
  0.7% in carbon; the floor takes the rest of 150 t and sets the stroke.

  | shape | η_jet unmixed-premixed | skirt | band wrap | floor (thickness) | stroke |
  |---|---|---|---|---|---|
  | dish 0.30 + 2 m skirt | 0.607-0.707 | 16 t | 36 t | 98 t (30 mm) | 2.8 m |
  | dish 0.25 + 2 m skirt | 0.598-0.697 | 18 t | 23 t | 109 t (36 mm) | 2.8 m |
  | dish 0.30, no skirt | 0.593-0.692 | 0 | 33 t | 117 t (36 mm) | 2.5 m |
  | dish 0.10 + 4 m skirt (old) | 0.583-0.681 | 48 t | 1 t | 101 t (39 mm) | 3.7 m |
  | dish 0.10, no skirt | 0.460-0.544 | 0 | 0 | 150 t (57 mm) | 2.5 m |

  A deep bowl with a short skirt beats the shallow dish with a tall one, by ~0.02-0.025 in η_jet,
  and every shape fits 150 t. The bowl's steep band needs its own hoop wrap (23-40 t at d/D
  0.25-0.30): it is pushed outward much as the skirt is. Caveats: effective-γ 2-D on the
  production grid (the cup dropped ~1% on a finer grid), the cold-end time scale, every ring free.

- **The face at 12 MN s (2026-10-05, Q64; `make water-plate-plate-thermal`).** Argon, respray,
  30 mm floor, ~92 kg PuffSat at 45.58 km/s and ~65 kg at 65.13 km/s (face 21.9 and 25.3 MJ/m^2):

  | standing pitch | steel peak | film removed per pulse (unshielded model) |
  |---|---|---|
  | 150 µm | 507 K | 28-33 kg |
  | 225 µm | 400 K | 28-33 kg |
  | 300 µm | 342 K | 28-33 kg |

  The steel peak matches the 8 MN s runs exactly: the pitch surface sits at its cap, so the
  conduction through the standing film sets the steel temperature, and a bigger pulse only burns
  more film. The unshielded removal exceeds the ~15 kg a 150 µm film holds over the footprint (the
  solver's grid does not thin as it removes, so it does not see burn-through). Shielded by its
  vapor, consumption scales to ~4-6 kg per pulse. 225-300 µm standing is the recommendation.

- **Gravity during the push (estimate, 2026-10-05).** 12 MN s at 4 Hz is 48 MN. On 1,500 t at
  400 km (local g 8.7 m/s^2) that is T/W 3.7 at the start, rising as argon is spent. A planar
  two-body run from rest at 400 km to 10.786 km/s takes ~250 s and spends ~700 t of argon at
  η_jet 0.6. Pushed horizontally the craft falls ~200 km, to ~195 km, while the stream's own
  escape hyperbola rises ~70 km over the ~1,350 km downrange. Holding 400 km costs only ~0.13 km/s
  (1.2%), but needs the thrust tilted up to ~16° at the start, and a plate cannot steer (the
  ledger's ADR-0009: the push axis is the stream's velocity). The likely answer is a lob that meets
  the stream still rising at ~1-1.2 km/s, so the booster pays the gravity loss. The ledger treats
  the push as impulsive; owed to `aim_is_all_you_need`.

- **A wide sheet at stand-off (scratch, 2026-10-05).** The k = 10 reaction mass as a thin layer
  (25 cm, 2 cells) spanning the footprint, 1-8 m off the cup's floor. A tethered sheet cannot leave
  the incoming path in the ~0.2 ms the pulse takes to pass, so it acts as a stand-off reaction
  layer, not a rebound tamper. It trails the deep cloud: 0.61-0.67 of ceiling against 0.71. It
  lowers the peak face pressure (0.78-0.95 against 1.52 in the cloud's units at 6-8 m off).

- The vapor opacity. κ_vapor moves argon's film between 6 and 27 kg per pulse, and only the 200
  and 1,500 values have a source. The EUV absorption of the film's vapor is the number to pin.
- **No survivable, useful configuration yet in 2-D.** Candidate levers: a footprint much wider
  than the cloud is deep (a larger plate), a dished plate (ADR-0021 raises capture), side
  confinement, and merging the pulse upstream so a shallow merged slab arrives.
- Step 2b: rim spill for the delivered disk in 2-D, with the areal mass matched to the taper.
- Step 3: **how good the bounce is**, end to end. Combine 2a's mixing with 2b's capture into one
  delivered impulse per PuffSat, quoted as the parent's η_jet. Compare it with this repository's
  bare-plate fudge factor `f(v)` at the same speeds.

- The necklace geometry: ring radius, sphere count, and the target delivered shape.
- The fraction of pulse energy radiated onto the face, argon versus water. This decides film mass
  and diameter.
- The 2-D losses the 1-D bound omits: rim spill and the radial profile.
- Water bonds that re-form while still pressure-coupled.
