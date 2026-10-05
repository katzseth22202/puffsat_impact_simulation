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

- `sec:plate_construction` should say the plate is **tapered** to the impulse profile, and why.
  Precise aim is a second reason.
- `sec:plate_liquid_spray` should state the **spray-cloud** requirement and the cushion-bounce
  failure.
- `sec:plate_aiming`: impact offset is capped at about 0.5 m by the taper.
- `sec:plate_ablative_film`: the film's thickness is set by heat once more than about 0.07% of
  pulse energy reaches the face, not only by the readout.
- Once the argon arm has solved values, replace the provisional η_jet ≈ 0.83 bound.

- **Argon's case in `sec:plate_liquid_spray` needs a heat caveat.** The premixed solve puts about
  10x more of the pulse onto the face for argon than for water. The 8-11% faster doubling the parent
  credits argon with does not yet pay for that film.

## Open

- Vapor shielding: rerun step 1 on the ablating wall. If shielding does not cut face radiation by
  one to two decades, the film mass decides argon against water.

- The necklace geometry: ring radius, sphere count, and the target delivered shape.
- The fraction of pulse energy radiated onto the face, argon versus water. This decides film mass
  and diameter.
- The 2-D losses the 1-D bound omits: rim spill and the radial profile.
- Water bonds that re-form while still pressure-coupled.
