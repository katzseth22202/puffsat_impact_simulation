# Spray plate: handoff to the companion repos (Draft 1)

**Draft 1, 2026-10-05. Pending the plug study, which may move P11.** From
`katzseth22202/puffsat_impact_simulation`. This lists every change the spray-plate study (ADR-0055)
asks of the parent paper (`Balloon-Pulse-Propulsion`) and the child paper
(`puffsats_for_datacenters`). Copy it into the parent as `docs/spray_plate_answers_from_impact_sim.md`,
following the parent's `*_answers_from_impact_sim.md` convention. The evidence for each number is in
this repository's `docs/argon_plate_handoff.md`. Section labels are from parent `f15d43c` and the
child's `templateArxiv.tex` as cloned on 2026-10-05.

**Status of the numbers.** Each item says whether it is a repository solve (and names the `make`
target) or an estimate. The solves are 1-D real-physics columns (Saha EOS, TOPS opacities, gray
flux-limited radiation), scaled by 2-D effective-γ geometry for spill and containment.

## The design that holds

On the overtake, the vehicle sprays `k` kilograms of argon or water per kilogram of PuffSat into a
cloud ahead of the plate. The working case is 100 t of plate, pulses at 4 Hz on a 2.5 m stroke
(8 MN·s per pulse), and 47 kg of PuffSat at 45.58 km/s or 33 kg at 65.13 km/s.

- **The plate is a shallow cup.** A dish (d/D 0.10-0.15) with a 4 m straight skirt at the 10 m rim,
  20 m across. The spray cloud is ~4 m deep and sits ~1 m off the floor (P10).
- **Maraging steel, tapered, aimed to 10 cm.** The face allowable is 2.5 GPa (P2, P3, P12).
- **A thick ablative film, ~150 µm, preferably pitch,** resprayed each cycle and mapped thermally
  (P5).
- **η_jet ≈ 0.58 (unmixed) to 0.68 (premixed)** for argon, against the paper's 0.775 (P11).
- **The injection ratio k is a mission trade,** not a plate optimum (P13).
- **Plate mass sets pulse size** through the gas spring (P14).

The fluid models hold at these speeds (estimate). Atoms meet at ~10-130 eV, where momentum-transfer
cross sections shrink to ~10^-20 m^2, yet the spray column, fixed by `k` and the pulse's areal mass,
still gives each atom ~10^6 collisions. Only a micron-thick edge layer is non-continuum.

## Parent (`Balloon-Pulse-Propulsion`)

### P1. Largest PuffSat water fragment: `sec:icy_puffsat`, the foam-scaffold paragraphs

Estimate. Any PuffSat water fragment must stop inside the spray column, ~6 kg/m^2.
- An unbroken drop of radius `r` stops after a column of ~`2.7 ρ_w r`, so up to ~2.2 mm.
- At 45 km/s a drop shatters after a column of ~`10 r √(ρ_w ρ_gas)`, so up to ~7 mm, and the
  fragments then stop quickly. A stopped fragment needs ~2.5 MJ/kg to boil, against ~85 MJ/kg
  around it.

The scaffold gives each strand ~25 mm^2 of water, capping fragments near 2.5 mm radius. That passes
with ~2.8x margin if breakup is credited, and misses by ~1.1x if not.

**Change.** Require the largest surviving fragment to have a radius of at most 1 mm, or keep the
~5 mm grid and state that the design relies on aerodynamic breakup.

### P2. A tapered plate: `sec:segmented_steel_plate`

Estimate. A uniform plate fails at any thickness. The footprint covers half the radius, so the
struck region takes ~4x the mean kick (~320 m/s against 80 m/s per pulse on 100 t). That is
~29 kJ/kg of relative motion, against the 100-300 J/kg steel stores elastically. Bending is
impulse-governed (the plate's bending period of tens of ms dwarfs the ~1 ms pulse), so the taper
is sized from the impulse profile.

**Change.** Taper the petals so areal mass follows the delivered impulse, matched to ~15-20%. The
plate is ~4 cm mean thickness, thicker at the center.

### P3. Aim as a second reason for the taper

Estimate. 10 cm aim is 2% of the 5 m footprint radius, well inside the taper's tolerance.

**Change.** Add that precise aim is what lets the taper work.

### P4. Steering by impact offset is capped: `sec:pusher_plate_down_kick` (lines ~983-988)

Estimate. An offset footprint breaks the uniform kick.

**Change.** Cap the offset at ~0.5 m, and use tilt plus thrusters for the rest.

### P5. A thick ablative film, mapped thermally: `sec:lightweight_pusher_plates` (line ~876)

Solved (`make water-plate-argon-film`, `make water-plate-plate-thermal`). Two numbers are set
separately. **Standing thickness** insulates the steel. **Consumption** is what each pulse costs,
and the respray replaces only that.

- **Thickness.** At ~150 µm standing, the maraging stays near 507 K over a 1,500-pulse push for
  both argon and water, with no cooling water. That is under its 753 K aging limit. A film at twice
  its consumption (36-55 µm argon, 3-13 µm water) is too thin: the steel reaches 800-2,500 K.
  Orion's ~6 mil (152 µm, GA-5009 Vol. III) is the right scale.
- **Consumption.** Shielded by its own vapor (κ_vapor ~5e3 m^2/kg, sourced in
  `docs/spray_plate_vapor_opacity.md`), pitch burns ~2-3 kg per pulse (argon) and ~0.2-0.7 kg
  (water). If the curtain loses its soot (κ ~100), argon's film rises ~4x.
- **Oil against pitch.** A carbon-loaded oil (PAO, 20% carbon) cracks at ~1,300 K (estimate)
  rather than holding 3,900 K, so it conducts less heat into the steel (~300-490 K). But it takes
  only ~1.6 MJ/kg to remove, against pitch's 59.3, so it burns 4-7x more film (12-18 kg per pulse for
  argon). Pitch at 150 µm already meets the aging limit, so pitch is preferred. Oil is the fallback
  if the steel runs hotter than modelled.

| Film | Burned per pulse (argon / water) | Steel peak over a push |
|---|---|---|
| Pitch, 150 µm | ~2-3 / ~0.2-0.7 kg | ~507 K, no coolant |
| Oil, 2x consumption | ~12-18 / ~1-4 kg | ~300-490 K, no coolant |

- **No metal interlayer** between film and steel (sourced, `docs/spray_plate_layer_properties.md`).
  A copper alloy would cut the face's temperature spike to ~30%, but it yields on every 2.5 GPa
  pulse (NASA's copper rocket liners thinned each cycle). Inconel 718 stays elastic but runs the
  face 1.2-1.6x hotter.
- **The cycle's order.** Within each 250 ms: an optional water spray to cool (~1.1 kg per pulse
  drops the steel from ~507 K to ~417 K, if it wets), then the cold film respray, then nothing
  liquid on the face when the next pulse arrives. Liquid left pooled acts as a wall (P7). Spray on
  hot steel can float on its own vapor (Leidenfrost) and cool poorly.
- **Wear mapping.** The child paper's near-infrared reflectance readout works only through films a
  few microns thick. A 150 µm film needs **thermal imaging** instead: map the face's temperature
  after each pulse, and respray where the film has thinned. This is how the parent already senses
  the chamber's pitch.

**Change.** Specify a ~150 µm ablative film, pitch preferred. Replace the thin-film optimization
with thermal-imaging wear mapping, and add the oil-versus-pitch trade.

### P6. Argon against water: `sec:water_injected_overtake`, and the argon loading near line 1384

Solved (`make water-plate-argon-radiation`, `-ablating`, `-k`).
- Argon puts ~10x more radiation on the face than water. It runs hotter (31-61 kK against
  14-30 kK) because it has fewer particles per kilogram, and radiation grows as T^4. Shielded, its
  extra film cost is small (~2-3 kg per pulse).
- On the cup at k = 10, argon delivers η_jet ~0.58 unmixed against water's ~0.50-0.56. Premixed,
  both reach ~0.63-0.72. Argon's edge is real but modest, and largest at the slow end.

**Change.** The section says its numbers are for water. Add argon's radiation and film cost, and
that its impulse edge is ~5-15%. The 8-11% faster doubling the parent credits argon with comes from
its own chemistry ledger, not these solves.

### P7. The spray must be a diffuse cloud: `sec:water_injected_overtake`

Solved (`make water-plate-argon-levers`, `-2d`). A thin, dense pooled layer acts as a wall, and
the PuffSat gas reflects off it near a bare bounce ("cushion bounce").

**Change.** State that the spray must form a diffuse cloud, ~4 m deep over the footprint, standing
~1 m off the plate floor. It must never pool on the face.

### P8. The necklace is the ring plus hub bag: formation section (lines ~533-547)

Nothing to change for the spray plate. The study models the delivered PuffSat as a disk matched to
the taper. A **plug** alternative is under study: compact necklace spheres burying themselves in
dense plugs on tethers, steered in the last ~100 ms, as in the parent's head-on rod-and-plug
idea without the sliding door.

### P9. Unmerged pulses and face pressure: `sec:water_injected_overtake`

Solved (`make water-plate-argon-mixing`, `-standoff`, `-2d-standoff`). A PuffSat layer that reaches
the plate before merging drives a shock through the spray into the steel, ~3x harder than on a bare
plate. That's because a dense, low-γ gas reflects shocks strongly. It reaches 4.7-7 GPa with a 1 m
cloud. With a 4 m cloud standing ~1 m off the floor, the peak is ~0.9-2.5 GPa, lasting 1-30 µs,
within maraging's 2.5 GPa.

**Change.** State that the cloud's depth and standoff, not cushioning, keep the face pressure
survivable.

### P10. The plate becomes a shallow cup: `sec:water_injected_overtake`, `fig:plate_stack`

Solved (`make water-plate-argon-2d-cup`). The merged gas is hot and expands sideways as readily as
backward, so an open plate loses ~23-30% of the impulse past its rim. Containment fixes it:

| Shape (10 m plate, 4 m cloud, 1 m gap) | Argon η_jet | Added steel |
|---|---|---|
| open plate | ~0.33 | none |
| dish d/D 0.15 | ~0.49 | none (reshaped) |
| straight skirt 4 m | ~0.53 | ~20 t per cm of wall |
| dish 0.10 + skirt 4 m | **~0.58** | ~20 t per cm of wall |

A flared, bell-like wall did no better than a straight skirt. The skirt's wall pressure is far
below the face's, so it can be thin, but it needs structural sizing and its own film and cooling.

**Change.** Show the plate as a shallow cup: dish plus straight skirt. Note that it sits between
the paper's "plate" and "nozzle". The skirt's mass (~20 t per cm of wall) must come from the 100 t
plate budget, or be added to it.

### P11. η_jet = 0.775 is not reached: `sec:water_injected_overtake`, `tab:mass_interest_growth`

Solved, pending the plug study. On the cup, argon reaches η_jet ~0.58 if PuffSat and spray do not
mix, and ~0.68 if they mix fully. Neither is established, and the difference is worth ~0.1. The
parent's growth tables use 0.775 as a requirement.

**Change.** None yet. If the plug study does not close the gap, rerun the growth tables at
η_jet ~0.6.

### P12. Plate material: maraging steel at a 2.5 GPa face allowable

Sourced (`docs/spray_plate_steel_face_limits.md`). **A short pulse rescues the plate's bending,
not its face.** The plate bends over tens of ms, so bending sees only the ~1 ms pulse's impulse.
A stress wave crosses the steel in ~7 µs, though, so the face sees the full pressure. Orion's lesson
that short pulses are survivable applies to the bulk structure. The face is loaded in uniaxial strain, so its
limit is the Hugoniot elastic limit (measured HY-100 1.89 GPa; maraging 350 4.8 ± 2.0 GPa).
In-plane prestress could raise it to ~3.5 GPa, but a 20 m, 4 cm plate buckles under 1-4 MPa of
in-plane load, so a rim band cannot apply it. The bulk must stay below maraging's ~480 °C aging
temperature. Aramid suits only tension (gas bags, tendons): its compressive strength is
0.14-0.23 GPa. Ceramic tiles do not help, because the steel behind them sees the full pressure.

**Change.** Name maraging 300/350 as the plate steel, with a 2.5 GPa face allowable and a ~480 °C
bulk limit.

### P13. The injection ratio k is a mission trade: `sec:water_injected_overtake`, the growth ledger

Solved (`make water-plate-argon-k`). On the cup, η_jet is nearly flat from k = 4 to 14. Raising k
buys impulse per PuffSat (β = J/mw from 2.3 to 3.7) and costs impulse per carried kilogram (0.46 w
to 0.21 w).

**Change.** Say that k is set by the launch-mass budget's weighting of PuffSat against carried mass.
The plate physics does not pick it.

### P14. Plate mass sets pulse size through the gas spring: `sec:pusher_plate_mass_stroke`

Estimate. With a preloaded gas spring (nearly constant force), the impulse one pulse can deliver is
`J ≈ 8 m_plate ν s`. That reproduces the child paper's 3.6 g peak on a 2.5 m stroke. 100 t at 4 Hz
and 2.5 m gives 8 MN·s per pulse, enough for a ~1,090 t vehicle at 3 g. At k = 10 that is ~47 kg of
PuffSat plus ~470 kg of argon per pulse at 45.58 km/s. The bag's fibre mass is set by force times
stroke (~80 MJ, ~70 kg of Vectran ideal, ~300 kg with margin), not by plate area, so bag pressure
is a free choice.

**Change.** Give the pulse-size relation, and note that plate area is free of the bag's mass.

## Child (`puffsats_for_datacenters`)

The child carries summaries of the parent, so each item follows its parent change.

- **C1, `sec:plate_ablative_film`.** Replace the 5 µm readout-sized film with a ~150 µm film, and
  the near-infrared readout with thermal-imaging wear mapping (P5). Drop the "0.1% of the pulse"
  film mass.
- **C2, `sec:plate_construction`.** Add the taper and aim (P2, P3), maraging steel (P12) and the
  cup with its skirt (P10).
- **C3, `sec:plate_aiming`.** Cap the offset (P4).
- **C4, `sec:plate_liquid_spray`.** Add the diffuse-cloud requirement with its depth and standoff
  (P7, P9). Quantify "argon runs hotter" as ~10x water's face radiation (P6). Give η_jet ~0.58-0.68
  (P11).
- **C6, `sec:plate_pulse_rate`.** Add `J ≈ 8 m_plate ν s` (P14).
- **C5, `sec:plate_pulse_rate`.** Survival is set by heat through the film and the 2.5 GPa face
  allowable, not by "tens of kilograms per square meter".
