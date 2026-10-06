# Spray plate: handoff to the companion repos (Draft 2)

**Draft 2, 2026-10-05.** From
`katzseth22202/puffsat_impact_simulation`. This lists every change the spray-plate study (ADR-0055)
asks of the parent paper (`Balloon-Pulse-Propulsion`) and the child paper
(`puffsats_for_datacenters`). Copy it into the parent as `docs/spray_plate_answers_from_impact_sim.md`,
following the parent's `*_answers_from_impact_sim.md` convention. The evidence for each number is in
this repository's `docs/argon_plate_handoff.md`. Section labels are from parent `f15d43c` and the
child's `templateArxiv.tex` as cloned on 2026-10-05.

**Status of the numbers.** Each item says whether it is a repository solve (and names the `make`
target) or an estimate. The solves are 1-D real-physics columns (Saha EOS, TOPS opacities, gray
flux-limited radiation), scaled by 2-D effective-γ geometry for spill and containment.

**Order of work.** The numbers flow in one direction, so the changes land in this order.

1. **`aim_is_all_you_need` first.** Recompute the growth ledger and the cost model for both plate
   designs at the paper's k = 8.52: the spray cup at η_jet = 0.6 (with 0.57 as the downside) and
   the plug at η_jet ~0.70 (P11). Keep 0.775 only as a reference row. Two conventions to watch in
   `src/growth_ledger.py` (checked 2026-10-05):
   - **Its `plate_efficiency` is the energy efficiency, η = η_jet².** η_jet 0.6 is η = 0.36, the
     0.57 downside is 0.32, and the plug's 0.70 is 0.49. The current grid `PLATE_EFFICIENCIES =
     (0.50, 0.70, 1.00)` is η_jet 0.71, 0.84 and 1.0, so the spray cup sits below all of it and
     needs its own row. (The paper's 0.775 is η = 0.60.)
   - **These η_jet values already charge the PuffSat's water bonds as lost.** This repository's
     argon table carries the PuffSat's water as atoms and charges its 50.9 MJ/kg before the gas
     reaches the plate. Run them with `NO_BONDS` (the "all argon" option), or the ledger's "argon
     on ice" option charges the same bonds a second time.
   - **Its plate is 150 t** (`PLATE_MASS`), not the 100 t this study first assumed. The design is
     re-sized to 150 t at 12 MN s per pulse (`docs/argon_plate_handoff.md`).
   - **Its growth push is impulsive.** At 12 MN s and 4 Hz (48 MN on 1,500 t, T/W 3.7 at 400 km)
     the push lasts ~250 s and spends ~700 t of argon (estimate). Unsupported, the craft falls
     ~200 km while the stream's hyperbola climbs ~70 km, and a plate cannot tilt its thrust to hold
     altitude (ADR-0009). Charge a lob that meets the stream still rising at ~1-1.2 km/s, or the
     ~0.13 km/s an altitude hold would cost, and check the craft stays on the stream's path.
2. **Then the parent paper.** Apply P1-P15, quoting the new ledger and cost figures where P11 and
   `tab:mass_interest_growth` need them.
3. **Then the child paper.** Apply C1-C8 once the parent's figures are in, since the child
   summarizes the parent and recomputes its own ledger and cost budget at η_jet = 0.6 (C8).

## The design that holds

On the overtake, the vehicle sprays `k` kilograms of argon per kilogram of PuffSat into a cloud
ahead of the plate. The working case is the ledger's **150 t plate** and its 1,500 t launch unit,
**12 MN·s per pulse at 4 Hz**, and k = 8.52: ~92 kg of PuffSat and ~790 kg of argon per pulse at
45.58 km/s, ~74 kg at 56.53 km/s (P14). The first studies used 100 t and 8 MN·s; their thermal and
face-pressure results carry over per unit area except where marked.

- **The plate is a deep bowl with a short skirt.** A dish of d/D 0.25-0.30 (5-6 m deep, 20 m
  across) with a 2 m straight skirt. The spray cloud is ~4 m deep and sits ~1 m off the floor
  (P10).
- **The floor slides inside a skirt fixed to the vehicle,** like a piston, so the skirt never takes
  the floor's kick. The skirt and the bowl's steep band carry hoop tension in a carbon wrap (P10).
- **Maraging steel, tapered, aimed to 10 cm.** The face allowable is 2.5 GPa (P2, P3, P12).
- **A thick ablative film, ~150 µm, preferably pitch,** resprayed each cycle and mapped thermally
  (P5). Its check at the 12 MN·s pulse is still owed.
- **Two designs (P11).** The spray cup flies first at η_jet = 0.6 (0.57 unmixed to 0.67 fully
  mixed). The plug, once perfected for the nozzle, reaches η_jet ~0.70 but needs last-minute
  aiming. Neither reaches the paper's 0.775.
- **The injection ratio k is a mission trade,** not a plate optimum (P13).
- **Plate mass sets pulse size** through the gas spring (P14).
- **Tamping the rebound does not work.** A trailing pearl or tethered sheet lowers η_jet. The spray
  cloud is already the tamper (P15).

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
struck region takes ~4x the mean kick (~320 m/s against 80 m/s per pulse on the first 100 t
study plate; the 150 t design's ~100-120 t floor takes ~100-120 m/s). That is
~29 kJ/kg of relative motion, against the 100-300 J/kg steel stores elastically. Bending is
impulse-governed (the plate's bending period of tens of ms dwarfs the ~1 ms pulse), so the taper
is sized from the impulse profile.

**Change.** Taper the petals so areal mass follows the delivered impulse, matched to ~15-20%. The
150 t design's floor is ~30-36 mm mean thickness (P10), thicker at the center.

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
  both argon and water, with no cooling water, at the first study's 8 MN·s pulse. A floor as thin
  as 15 mm changes this by ~10 K (518 K). At the 12 MN·s pulse each pulse carries ~2x the energy
  per area, and that rerun is still owed. That is under its 753 K aging limit. A film at twice
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
the taper. The **plug** design that follows it (P11) uses compact spheres instead.

### P9. Unmerged pulses and face pressure: `sec:water_injected_overtake`

Solved (`make water-plate-argon-mixing`, `-standoff`, `-2d-standoff`). A PuffSat layer that reaches
the plate before merging drives a shock through the spray into the steel, ~3x harder than on a bare
plate. That's because a dense, low-γ gas reflects shocks strongly. It reaches 4.7-7 GPa with a 1 m
cloud. With a 4 m cloud standing ~1 m off the floor, the peak is ~0.9-2.5 GPa, lasting 1-30 µs,
within maraging's 2.5 GPa.

**Change.** State that the cloud's depth and standoff, not cushioning, keep the face pressure
survivable.

### P10. The plate becomes a deep bowl with a short skirt: `sec:water_injected_overtake`, `fig:plate_stack`

Solved (`make water-plate-argon-2d-cup`, `-2d-skirt`, `-2d-shape`). The merged gas is hot and
expands sideways as readily as backward, so an open plate loses ~23-30% of the impulse past its
rim. Containment fixes it, and it is not light. The hot layer pushes the wall outward nearly as
hard as it pushes the floor down, ~20 kPa s per pulse at the skirt's base against the floor's 25
(first study's pulse). The wall can only stop that push in hoop tension.

Sized inside the ledger's 150 t at 12 MN·s, k = 8.52, with carbon hoop wrap at 0.7% strain:

| Shape | Argon η_jet (unmixed-premixed) | Skirt | Bowl-band wrap | Floor | Stroke |
|---|---|---|---|---|---|
| **dish 0.30 + 2 m skirt** | **0.607-0.707** | 16 t | 36 t | 98 t (30 mm) | 2.8 m |
| dish 0.25 + 2 m skirt | 0.598-0.697 | 18 t | 23 t | 109 t (36 mm) | 2.8 m |
| dish 0.30, no skirt | 0.593-0.692 | 0 | 33 t | 117 t (36 mm) | 2.5 m |
| dish 0.10 + 4 m skirt | 0.583-0.681 | 48 t | 1 t | 101 t (39 mm) | 3.7 m |
| dish 0.10, no skirt | 0.460-0.544 | 0 | 0 | 150 t (57 mm) | 2.5 m |
| open flat plate | ~0.33-0.40 | 0 | 0 | 150 t | 2.5 m |

- **A deep bowl with a short skirt wins.** It beats the shallow dish with a tall skirt by ~0.02 in
  η_jet and cuts the skirt from 48 t to 16-18 t. The bowl's steep band is pushed outward too, so it
  needs 23-36 t of its own hoop wrap. Every shape fits 150 t. d/D 0.25-0.30 are within the grid's
  ~1% error of each other.
- **The skirt is fixed to the vehicle and the floor slides inside it.** A skirt bolted to the floor
  would have to take the floor's ~100 m/s kick within the ~1 ms pulse, a stress wave of ρcΔv,
  ~1-3 GPa at the joint. Fixed to the vehicle, it never moves. The cost is a hot sliding seal at the
  floor's rim (estimate: ~1.6 kg and ~140 MJ of gas per pulse through a 10 mm gap, ~14 MJ through
  1 mm). A 20 m floor grows ~22 mm in diameter per 100 K, so the gap needs a sprung, segmented seal
  or a labyrinth vented outward, away from the gas bags. Only the floor is sprung (P14).
- A flared, bell-like wall did no better than a straight skirt, an inward lip did worse, and a
  taller skirt adds ~1%.
- **A bare plate cannot use a tamper instead (solved).** On a flat plate or a d/D 0.10 dish, a
  trailing pearl or a thin stand-off layer lowers the impulse or leaves it unchanged. A bare plate
  loses its impulse to sideways spill, and a tamper pushes the gas back down, so it holds the hot
  layer under pressure longer and more of it escapes past the rim (see P15).

Caveats: effective-γ 2-D on the production grid, scaled to the real pulse; each wall ring sized as
if free; the cold-end time scale.

**Change.** Show the plate as a deep bowl (d/D 0.25-0.30) with a 2 m skirt fixed to the vehicle,
the floor sliding inside it. Note that it sits between the paper's "plate" and "nozzle", and give
the hoop wraps and the sliding seal.

### P11. Two plate designs: the spray cup now, the plug later: `sec:water_injected_overtake`, `tab:mass_interest_growth`

Solved for the spray cup (`make water-plate-argon-k`, `make water-plate-argon-2d-cup`), solved in
2-D only for the plug (`make water-plate-plug-2d`). Neither design reaches η_jet = 0.775. The plug
comes close and the spray cup does not.

**The spray cup is the initial design.** The vehicle sprays argon into a cloud ahead of the deep
bowl and its skirt (P10). Argon reaches η_jet ~0.57-0.61 if PuffSat and spray do not mix and
~0.67-0.71 if they mix fully, across the shallow cup and the deep bowl. The necklace charges mix them partly, and how far is not established. **η_jet = 0.6
is the baseline**, a little above the unmixed floor. 0.57 is the downside.

**The plug is the design to move to once it is perfected for the nozzle.** Each PuffSat sphere
buries itself in a plug, a dense, narrow body of carried mass hung on tethers ahead of the plate.
The plug forces the merge that the spray leaves to chance. In the cup it reaches 0.81-0.82 of the
overtake ceiling against the spray cup's 0.71 on the same 2-D method (converged to ~2% on a grid
twice as fine). That is η_jet ~0.75-0.77 on that method. Carried over the way the spray cup's
full-physics result scales, it is **η_jet ~0.70** (an estimate).

What makes the plug easier:

- **It needs no exploding liquid PuffSat.** The spray cup needs a necklace of liquid charges,
  shaped explosively to the ideal mixing shape. The plug needs only two solid pieces meeting,
  the PuffSat and the plug. There is no droplet size to hold and no shaped charge.
- **It needs no spray system or cloud.** The reaction mass is a solid body, not a mist held in
  place at the right depth.

What makes it harder:

- **Last-minute aiming.** The plug must sit in the sphere's path to millimetres. That needs the
  tethered-plug trick the paper uses for the head-on nozzle: plugs hung beyond the blast radius,
  steered in the last ~100 ms to the measured position of each incoming sphere. It is easier here
  than head-on, because no sliding door has to open. It still needs a fresh tethered plug per
  sphere per pulse, and a miss costs most of that pulse.
- **It needs more enclosure.** The merge makes a nearly round fireball. With no skirt the plug
  delivers only 0.32 of the ceiling. Its best cases were on a shallow dish with an 8 m skirt; it
  has not yet been run in the deep bowl.
- **The plug is argon ice in a cage of stronger material.** Argon has no bonds to break, which
  is why the paper already freezes its head-on plug from argon. Polyethylene was considered and set
  aside. Breaking all its bonds costs ~84 MJ/kg, and at 45.58 km/s and k = 10 the collision
  dissipates only ~86 MJ per kilogram of merged gas. The hydrogen share (~31 MJ/kg) could come
  back as H2 forms while the gas still presses on the plate. The carbon share (~53 MJ/kg) comes
  back only below ~4000 K, after the gas has left (estimates). Argon ice is soft and creeps, so
  the tethers grip a thin cage (aramid net or aluminium sleeve, a few percent of the plug's mass)
  rather than the ice. At 10 g a 0.4 m plug sees ~0.06 MPa of inertial stress. 10 g moves it
  ~25 cm in the last 100 ms (accelerate then stop), at ~4.6 kN of tether pull for a 47 kg plug.
  Argon ice sublimates fast in vacuum unless cold (~1 mm per minute at 40 K, ~1 mm per 10 s at
  50 K), so plugs are stored cold and shaded aboard and deployed seconds before each pulse.
- **One thing is not yet in its number:** the real-physics radiation onto the face, and argon's
  ionization, which the 2-D effective-γ runs do not carry.

**Change.** In `sec:water_injected_overtake`, present the spray cup as the initial design at
η_jet = 0.6 and the plug as the later design at η_jet ~0.70, with the upsides and downsides above.
**Compute the growth ledger and the cost model for both** in the companion calculation repository.
Use η_jet = 0.6 for the spray cup and 0.70 for the plug, at the paper's k = 8.52, and replace the
0.775 requirement in `tab:mass_interest_growth` with these two rows. Keep 0.775 only as a
reference.

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
`J ≈ 8 m ν s` for sprung mass `m`, rate `ν` and stroke `s`. That reproduces the child paper's 3.6 g
peak on a 2.5 m stroke. With the skirt fixed to the vehicle (P10), only the floor is sprung. The
150 t design's ~98-117 t floor gives **12 MN·s per pulse at 4 Hz on a 2.5-2.8 m stroke**. At
k = 8.52 and η_jet = 0.6 (β = J/mw = 2.85) that is **~92 kg of PuffSat plus ~790 kg of argon at
45.58 km/s**, and ~74 kg at 56.53 km/s. The average thrust is 48 MN, 3.3 g on the 1,500 t launch
unit. The bag's fibre mass is set by force times stroke (~80 MJ, ~70 kg of Vectran ideal, ~300 kg
with margin at 8 MN·s), not by plate area, so bag pressure is a free choice.

**Change.** Give the pulse-size relation, the 12 MN·s working pulse with its PuffSat mass, and note
that plate area is free of the bag's mass.

### P15. Tamping the rebound was explored and does not work: `sec:water_injected_overtake`

Solved for the formation case (`make water-plate-argon-2d-pearl`), estimated for the rest. Most of
the gap between η_jet and 1 is speed spread, not heat. Gas leaving a wall leaves at a spread of
speeds, and a linear spread caps a single pulse at √3/2 ≈ 0.87 of the ceiling. A tamper could in
principle close that gap. It is a mass that meets the rebounding gas, slows its fast front and lets
the whole charge leave at one speed. We tried three ways to place one.

- **A slice shaped off the PuffSat itself.** The explosive would throw a small part of the gas
  backward so it arrives late. The delay at the plate is about `D·Δv/w²`. A 1 ms delay from release
  1.8 km out needs ~1.2 km/s of backward kick. The same kick spreads the slice sideways by tens of
  metres, wider than the plate. One event lacks both the time and the aim.
- **A second, lighter PuffSat flying in formation behind the lead.** Same architecture, 5-30% of
  the lead's mass, trailing by 5-40 m. Every case did worse than one PuffSat carrying the same total
  mass into the same cloud. η_jet fell by 0.03 at 5% and 5 m, and by 0.18 at 30% and 40 m.
- **A thin sheet on tethers, such as 100 µm Kapton.** The pulse passes in ~0.2 ms, and no tether
  can move a sheet out of its path that fast. The incoming gas hits the sheet first, vaporizes it
  and carries it to the plate. The sheet becomes part of the incoming charge before any rebound
  exists. Used instead as the reaction mass, a thin sheet 1-8 m off the floor gave 0.61-0.67 of the
  ceiling, against 0.71 for the deep spray cloud.

Why tamping fails here:

- **The tamper is too light.** To slow gas, a tamper must weigh about as much as the gas. The
  rebounding gas is the PuffSat plus k times its mass in spray. A 5-30% pearl is 0.5-3% of it.
- **Collision turns its mass into heat.** The pearl meets gas that is already leaving the plate.
  They close head-on at ~60 km/s, the most dissipative collision there is. Most of the pearl's
  energy becomes heat rather than push.
- **Its mass misses the spray's leverage.** Mass that arrives with the lead shares its energy with
  spray at rest, and that sharing is what raises the impulse. The pearl's mass gets no spray.
- **Late arrivals find an empty cup.** At long gaps the pearl makes a bare bounce, worth about 2
  instead of about 3. So the loss grows with the gap.

The spray cloud is already the tamper. It is a large mass, at rest, in place before the pulse
arrives. That is the only arrangement that worked.

**Change.** In `sec:water_injected_overtake`, add a short paragraph after the η_jet discussion.
Say that tamping the rebound to narrow its speed spread was tested and rejected, and give the
reasons above. The collision point is the one to keep. Name the spray cloud as the tamper.

## Child (`puffsats_for_datacenters`)

The child carries summaries of the parent, so each item follows its parent change.

- **C1, `sec:plate_ablative_film`.** Replace the 5 µm readout-sized film with a ~150 µm film, and
  the near-infrared readout with thermal-imaging wear mapping (P5). Drop the "0.1% of the pulse"
  film mass.
- **C2, `sec:plate_construction`.** Add the taper and aim (P2, P3), maraging steel (P12), and the
  deep bowl with its 2 m skirt fixed to the vehicle and the floor sliding inside it (P10).
- **C3, `sec:plate_aiming`.** Cap the offset (P4).
- **C4, `sec:plate_liquid_spray`.** Add the diffuse-cloud requirement with its depth and standoff
  (P7, P9). Quantify "argon runs hotter" as ~10x water's face radiation (P6). Give η_jet = 0.6 as
  the baseline, 0.57-0.71 as the range (P11).
- **C6, `sec:plate_pulse_rate`.** Add `J ≈ 8 m ν s` and the 12 MN·s pulse on the 150 t plate,
  ~92 kg of PuffSat at the cold end (P14).
- **C7, `sec:plate_liquid_spray`.** One sentence. Trailing tampers were tested and lower the
  impulse, because the late mass meets the rebound head-on and turns to heat (P15).
- **C8, `sec:plate_liquid_spray` and `sec:exponential_mass_growth`.** Say that the parent finds a
  plug design more efficient (η_jet ~0.70, P11) and cite it. Explain why this paper does not use
  it: each plug must be steered into a sphere's path to millimetres in the last ~100 ms, and the
  design is not yet proven for the nozzle. Then **recompute the growth ledger and the cost budget at
  η_jet = 0.6**, the spray cup's baseline, in place of 0.775.
- **C5, `sec:plate_pulse_rate`.** Survival is set by heat through the film and the 2.5 GPa face
  allowable, not by "tens of kilograms per square meter".
