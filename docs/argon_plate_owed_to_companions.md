# Spray plate: changes owed to the companion repos

From `katzseth22202/puffsat_impact_simulation`, started 2026-10-05. This lists every change the
argon arm of the spray plate (ADR-0055) asks of the parent paper (`Balloon-Pulse-Propulsion`) and
the child paper (`puffsats_for_datacenters`). It is written to be copied into the parent as
`docs/spray_plate_answers_from_impact_sim.md`, following the parent's `*_answers_from_impact_sim.md`
convention. The running numbers and their evidence are in this repository's
`docs/argon_plate_handoff.md`.

Section labels and line numbers below are from parent `f15d43c` and child `templateArxiv.tex` as
cloned on 2026-10-05.

**Status of the numbers.** Steps 1, 1b, 2a, 2a' and 2b are repository solves:
`make water-plate-argon-radiation`, `make water-plate-argon-ablating`,
`make water-plate-argon-mixing`, `make water-plate-argon-levers` and `make water-plate-argon-2d`. They are 1-D with gray flux-limited radiation, premixed (1, 1b) or
layered (2a). Everything else here is an estimate and is labelled as one. Rim spill (step 2b) and
the end-to-end bounce (step 3) are still to come, and may move several items.

## The setting

On the overtake, the vehicle sprays `k = 10` kilograms of argon or water per kilogram of PuffSat
into a cloud ahead of the plate. The working case is:

- 100 t of steel plate, tapered, 20 m in loaded diameter, held to 10 cm aim.
- Pulses at 4 Hz on a 2.5 m stroke, giving 8 MN·s per pulse.
- 47 kg of PuffSat at 45.58 km/s, or 33 kg at 65.13 km/s, onto a 78.5 m^2 footprint.

## Parent (`Balloon-Pulse-Propulsion`)

### P1. Largest PuffSat water fragment: `sec:icy_puffsat`, the foam-scaffold paragraphs

Any PuffSat water fragment must stop inside the spray column, about 6 kg/m^2. A fragment that
reaches the plate whole strikes the film as a projectile.

- **Stopping without breakup.** A drop of radius `r` stops after a column of about `2.7 ρ_w r`.
  That allows `r` up to about 2.2 mm.
- **Breakup first.** At 45 km/s a drop shatters after a column of about `10 r √(ρ_w ρ_gas)`, which
  is about `810 r` in a 1 m cloud (6.6 kg/m^3). That allows `r` up to about 7 mm. The fragments
  then stop within a small fraction of the column.
- **Boiling.** A stopped fragment needs about 2.5 MJ/kg to boil, against about 85 MJ/kg around it.

The scaffold gives each strand about 25 mm^2 of water, so the grid caps fragments near 5 mm across
(2.5 mm radius). That passes with about 2.8x margin if breakup is credited. It misses by about 1.1x
if not.

**Change.** State the requirement that the largest surviving fragment has a radius of at most
1 mm, which stops fragments without relying on breakup. Alternatively, keep the ~5 mm grid and say
the design relies on aerodynamic breakup. These are drag and breakup estimates, not a simulation.

### P2. A tapered plate: `sec:segmented_steel_plate`

A uniform plate cannot be used at any thickness. The footprint covers half the radius, so the
struck region takes four times the mean kick: about 320 m/s against 80 m/s per 8 MN·s pulse on
100 t. That is about 29 kJ/kg of relative motion, against the 100-300 J/kg high-strength steel
stores elastically. Steel stiffness cannot spread the load within a ~0.5 ms pulse either. A
10 m, 16 cm plate rings at about 14 Hz.

**Change.** The petals should be tapered so areal mass follows the delivered impulse profile, so
every element takes the same kick. The working plate is 20 m loaded diameter and about 4 cm mean
thickness, thicker at the center. The taper must match the impulse to about 15-20%. These are
estimates.

### P3. Aim as a second reason for the taper: `sec:segmented_steel_plate` or `sec:pusher_plate_down_kick`

A taper is only as good as where the footprint lands. 10 cm aim is 2% of the 5 m footprint radius,
well inside the taper's tolerance.

**Change.** Add that precise aim is what lets the taper work.

### P4. Steering by impact offset is capped: `sec:pusher_plate_down_kick`, the tilt-and-offset passage (lines ~983-988)

An offset footprint breaks the uniform kick the taper relies on.

**Change.** Cap the offset at about 0.5 m (about five aim-sized steps), and let tilt plus thrusters
do the rest.

### P5. The ablative film's thickness is set by heat: `sec:lightweight_pusher_plates` (line ~876)

Solved, step 1b. With the film's vapor shielding the steel, film consumed per pulse is:

| Fluid | Film per pulse | Film thickness |
|---|---|---|
| Argon | 6-18 kg | ~75-230 µm |
| Water | 1-5 kg | ~10-70 µm |

These are at vapor opacities of 1,500-10^4 m^2/kg and Q* = 5 MJ/kg. Without shielding, argon would
need 120+ kg per pulse. The sourced vapor opacity is 5e3 m^2/kg (range 2e3-3e4), from this
repository's `docs/spray_plate_vapor_opacity.md`. That holds while the vapor keeps its soot. If the
plasma side of the curtain loses it, the opacity falls to ~100, and argon needs 26-48 kg per pulse
(step 2a). Heat numbers at production resolution are high by ~1.5-2x, so these are conservative.

**Change.** Keep the film about **150 µm thick** (Orion's ~6 mil, 152 µm, is right), and respray
what each pulse burns: ~2-3 kg (argon) or ~0.2-0.7 kg (water) of pitch. At that thickness the steel
stays near 507 K over a push with no cooling water. A 5 µm readout-sized film is far too thin.
Carbon-loaded oil keeps the steel cooler but burns 4-7x more, so the paper's pitch stands.

### P6. Argon against water: `sec:water_injected_overtake`, and the argon loading near line 1384

Solved, steps 1 and 1b, premixed.

- On a bare cold face, argon puts 1.5-11% of the pulse's kinetic energy on the plate, against
  water's 0.17-2.7%. Argon runs hotter (31-61 kK against 14-30 kK) because it has fewer particles
  per kilogram, and radiation grows as T^4.
- With vapor shielding, argon delivers 0.74-0.75 of the overtake ceiling `1 + √(1+k)` at
  45.58 km/s against water's 0.68, and they tie near 0.745 at 65.13 km/s. That is η_jet about
  0.65 at the cold end, below the 0.775 requirement before rim spill.
- These favour water. The argon table pays water's bond energy permanently, while the water table
  returns all of it.

**Change.** The section opens by saying its numbers are for water. Add that argon's case now has a
film cost (P5) and a premixed η_jet near 0.65. Say also that the 8-11% faster doubling credited to
argon is on the parent's own chemistry ledger, not yet on these solves.

### P7. The spray must be a diffuse cloud: `sec:water_injected_overtake`

Provisional 1-D estimate. A thin, dense pooled layer acts as a wall. The PuffSat gas reflects off
it (cushion bounce) at about 0.57 of the ceiling, near a bare bounce at 0.46. A cloud at least as
deep as the arriving pulse holds about 0.79 even unmixed. The premixed solve puts the depth trade
near 1-2 m: 0.5 m exceeds the steel's pressure ladder at 65 km/s, and deeper clouds radiate more.

**Change.** State that the spray must be diffuse, 1-2 m deep over the footprint, never pooled.

### P8. The necklace is the ring plus hub bag: formation section (lines ~533-547)

There is nothing to change. This repository models the delivered PuffSat as a disk whose areal mass
matches the taper, which is what the existing ring plus hub bag produces.

### P9. The plate survives only a merged pulse: `sec:water_injected_overtake` and `sec:segmented_steel_plate`

Solved, step 2a. When a PuffSat layer reaches the plate before merging with the spray, it drives
a shock through the spray into the steel. The converged peak is ~4.7 GPa (argon) and ~7 GPa
(water) at 45.58 km/s with a 1 m cloud. That is 5-10x the 400-900 MPa ladder, and ~3x what the
same pulse would put on a bare plate. Interleaving the layers does not help once resolved. The
premixed case stays at 0.17-0.53 GPa.

Solved, step 2a'. A deep spray cloud is the fix: at 16 m deep with a 4 m pulse, the converged
peak falls to ~0.34-0.52 GPa (argon) and ~0.57-0.77 GPa (water). A longer pulse on its own lowers
the peak but collapses the impulse to ~0.47 of the ceiling.

**Withdrawn by step 2b.** In 2-D a cloud that deep spills sideways off the 10 m plate and keeps
only 30% of the impulse. Shallow clouds keep the impulse but get no pressure relief.

**Change.** State that survival of an unmerged pulse is unresolved. No cloud depth on the 10 m
plate gives both a survivable peak and a useful impulse. A wider footprint, a dished plate, side
confinement or upstream merging are the open options.

### P10. η_jet = 0.775 is not reached in any survivable 1-D configuration: `sec:water_injected_overtake`, `tab:mass_interest_growth`

Provisional, pending the plug study. The k sweep shows η_jet nearly flat in k on the cup: argon
~0.58 unmixed to ~0.68 premixed from k = 4 to 14, so no choice of k reaches 0.775. On an open 10 m plate the best
survivable case reaches η_jet ~0.33. Shaping the plate as a dish (d/D 0.10-0.15) with a 4 m
straight skirt contains the spill and reaches **η_jet ~0.58** (argon, 45.58 km/s; ~0.56 water),
within maraging's 2.5 GPa. The parent's growth tables are run at η_jet = 0.775 as a requirement.

A paper consequence follows. The overtake plate becomes a shallow cup, between the paper's
"plate" and "nozzle". `sec:water_injected_overtake` and `fig:plate_stack` should show the dish
and skirt.

**Change.** None yet. Step 3 will give the end-to-end number. Flagged now because the gap is large
enough to move the doubling times.

## Child (`puffsats_for_datacenters`)

The child carries summaries of the parent, so each item follows its parent change.

- **C1, `sec:plate_ablative_film`.** The 5 µm film is set by the infrared readout, and heat now sets
  a far thicker film (P5). Keep the readout argument for remapping between pulses, but drop the
  claim that the thin film suffices, and the "0.1% of the pulse" film mass.
- **C2, `sec:plate_construction`.** Add the taper (P2) and aim (P3).
- **C3, `sec:plate_aiming`.** Cap the offset (P4).
- **C4, `sec:plate_liquid_spray`.** Add the diffuse-cloud requirement (P7). Also note that argon
  "runs hotter at the plate" is now quantified: about 10x the face radiation of water (P6).
- **C6, `sec:plate_liquid_spray`.** Add that the plate survives only a merged pulse (P9).
- **C5, `sec:plate_pulse_rate`.** "Of order tens of kilograms per square meter" as the survival
  floor is not the binding limit. Heat is, through the film (P5).
