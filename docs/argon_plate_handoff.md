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

**Do the fluid models hold at these speeds?** Atoms meet at about 10-130 eV in the
center-of-mass frame, where momentum-transfer cross sections shrink to roughly 10^-20 m^2 (an
estimate, not sourced). The argon column is fixed by `k` and the pulse's areal mass, about 9e25
atoms/m^2 over the footprint. That still gives about 10^6 collisions, and about 10^4 even with a
cross section 100x smaller. Spreading the cloud deeper does not thin the column. Only spray thrown
outside the footprint does. The fluid treatment fails only in a free-path-thick edge layer,
microns thick. Once the gas ionizes, Coulomb collisions take over, with larger cross sections.

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
