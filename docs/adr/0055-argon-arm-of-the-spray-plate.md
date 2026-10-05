# ADR-0055: An argon arm for the spray plate, on a tapered 20 m plate

Status: proposed, 2026-10-05 (grill with Seth). Extends the water-injected overtake plate study
([`water_injected_overtake_plate_study.md`](../water_injected_overtake_plate_study.md)).

## Context

`puffsats_for_datacenters` (`templateArxiv.tex:234-246`, §`plate_liquid_spray`) sprays water or
argon into the Jupiter growth wave ahead of the overtake plate. It quotes the **overtake ceiling**
`J_max = m w (1 + √(1 + k))` and the parent's requirement `η_jet = 0.775`, neither of them simulated.
This repository has a water arm only. Nothing here or in either companion repo prices argon at the
plate, caps plate mass, or says how the plate survives repeated pulses.

## Decision

1. **The argon arm lives in this repository, beside the water arm.** It shares the injection ratio
   `k`, the four overtake speeds (45.58 / 56.53 / 61.83 / 65.13 km/s) and the overtake geometry. It
   reuses the shared hydrocode kernels where they fit, and adds its own where they do not.
2. **Working case: `k = 10` and 100 t of steel plate.**
3. **The plate is a tapered plate, 20 m in loaded diameter.** Its areal mass follows the delivered
   impulse, so every element takes the same kick, about 80 m/s per 8 MN·s pulse. A uniform plate
   cannot be used at any thickness. A footprint covering half the radius gives the struck region
   four times the mean kick. The resulting ~29 kJ/kg of relative motion is about 100x what steel
   stores elastically, and steel stiffness cannot spread the load within a ~0.5 ms pulse. Diameter
   stays a working value, to be re-optimized against heat absorbed and ablative mass.
4. **Aim is held to 10 cm.** That is 2% of the 5 m footprint radius, well inside the ~15-20%
   impulse-to-mass mismatch the taper tolerates. Precise aim is a second reason the taper works,
   and the paper should say so. Steering by impact offset stays within about 0.5 m.
5. **The spray must form a spray cloud, never a pooled layer.** A provisional 1-D estimate puts a
   thin dense layer at about 0.57 of the ceiling (cushion bounce). A diffuse cloud as deep as the
   pulse holds about 0.79 even unmixed.
6. **PuffSats arrive as necklace charges.** The explosive shaping is assumed to succeed, so the
   delivered shape is an input to optimize for mixing, not a formation process to simulate.
7. **Heat is the binding survivability criterion, ahead of face pressure.** Face pressure at 20 m
   is about 51 MPa, against a 400-900 MPa ladder. Heat is not comfortable. 100 t of steel absorbs
   about 25 GJ over a 500 K rise, and a 5 µm film absorbs about 16 MJ per pulse. Together they
   cover only about 0.07% of each 49 GJ pulse across the ~1,500 pulses of the Earth push. The
   fraction of pulse energy radiated onto the face therefore sets the film mass, and so the diameter.
8. **Water is charged only for bonds that re-form while the gas still presses on the plate.** The
   31% return by 250 µs from Step 11 is preliminary. The argon arm's flow solve should replace it.

9. **The plate is maraging steel (300/350 class), with no in-plane prestress.** The face allowable
   is **2.5 GPa** peak normal stress. It is set by the Hugoniot elastic limit, because the ~0.1-1
   ms pulse outlasts the ~7 µs stress-wave transit (impulse governs bending, but pressure governs
   the face). The source is `docs/spray_plate_steel_face_limits.md` (measured HELs: HY-100 1.89
   GPa, Thomas et al. 2018; maraging 350 4.8 ± 2.0 GPa, Gust & Royce 1970). Prestress would allow
   ~3.5 GPa, but a 20 m, 4 cm plate buckles under 1-4 MPa of in-plane load, so a rim band cannot
   apply it. The bulk must stay well below maraging's ~480 °C aging temperature, which tightens
   decision 7's 500 K bulk rise. Aramid serves only in tension (gas bags, tendons). Ceramic tiles
   do not raise the steel's limit, because the steel behind them sees the full pressure.

10. **The plate is a shallow cup: a dish (d/D 0.10-0.15) with a 4 m straight skirt at the
    10 m rim.** The spray cloud sits ~1 m off the floor. Spill past the rim was the largest single
    loss. Containing it raises η_jet from ~0.33 to ~0.58 (argon, 45.58 km/s; step 2d), and the
    gap lowers the face peak. The skirt is straight: a flared wall gave no more impulse for more
    steel. The skirt (~20 t per cm of wall thickness) is sized for its own wall pressure, which is
    far below the face's, and gets film and cooling like the face.

11. **The film is pitch, kept about 150 µm thick and resprayed cold each cycle.** Standing
    thickness insulates the steel; consumption is what each pulse costs, and the two are set
    separately. Shielded by its vapor, pitch burns ~2-3 kg per pulse (argon) and ~0.2-0.7 kg
    (water). At 150 µm the steel peaks at ~507 K over 1,500 pulses for both arms with no cooling
    water. Carbon-loaded oil keeps the steel cooler (~300-490 K) but burns 4-7x more film (12-18 kg
    per pulse for argon), so it is not adopted.
12. **The injection ratio k is a mission trade, not a plate optimum.** On the cup, η_jet is nearly
    flat in k: argon ~0.58 unmixed to ~0.68 premixed from k = 4 to 14. Raising k buys impulse per
    PuffSat (β 2.3 -> 3.7) at the cost of impulse per carried kilogram (0.46 -> 0.21 w). k stays with
    the parent's launch-mass budget (8.5-10), whose cost weights decide it.

## Considered options

- **A small, thick uniform plate.** Rejected, because its mismatch between impulse and mass is
  independent of thickness (decision 3).
- **A Kevlar (aramid) face under an ablative.** Rejected: 0.14-0.23 GPa compressive strength and
  a ~200 °C limit.
- **In-plane prestress by a rim band.** Rejected: the thin plate buckles long before the useful
  prestress (0.5 Y) is reached.
- **A flared bell wall instead of a straight skirt.** Rejected: 0.98-1.01 against the straight
  skirt's 1.02 of the 1-D impulse, for 10-30% more steel.
- **An open flat plate.** Rejected: η_jet ~0.33, because hot merged gas spills past the rim.
- **A carbon-loaded oil film.** Rejected for 4-7x the film consumption of pitch, when pitch at
  150 µm already holds the steel under its aging limit without coolant.
- **Steering by impact offset alone, as the child repo describes.** Limited, because large offsets
  break the uniform kick. Tilt and thrusters do the rest.

## Consequences

- Plate mass sets the impulse per pulse only through the gas spring: `J ≈ 8 m_plate ν s`, which is
  8 MN·s at 4 Hz and a 2.5 m stroke. That is about 47-49 kg of PuffSat plus 470-490 kg of argon at
  45.58 km/s. Fiber mass in the gas spring tracks force times stroke, not area, so bag pressure is
  a free choice.
- The first solve (premixed, 1-D, cold black face) puts 1.5-11% of pulse kinetic energy onto the
  face for argon and 0.17-2.7% for water, against the ~0.07% that steel plus film absorb. Decision
  7 is therefore confirmed as binding. The impulse lands near 0.7 of the ceiling for both, so argon
  is not yet shown to beat water.
- On an ablating, vapor-shielded face (step 1b), heat becomes film mass. That is 6-18 kg per pulse
  for argon and 1-5 kg for water at the sourced-to-estimated vapor opacities, so the film is
  tens to hundreds of microns, not 5 µm. Argon then leads water by ~10% at 45.58 km/s and ties it at
  65.13.
- Q9 resolved: the necklace is the parent's ring plus hub bag, delivered as a disk whose areal mass
  is matched to the taper (decision 6 stands). Q8 resolved: the necklace matches the taper.
- Step 2a (layered, two materials with radiation): poor mixing costs 0.05-0.10 of the ceiling, and
  an unmerged PuffSat layer drives a converged ~4.7-7 GPa shock into the steel. **Peak pressure,
  not heat, binds when the pulse reaches the plate unmerged.** Decision 7 is amended accordingly.
  Face radiation at production resolution is high by ~1.5-2x, so the film figures are
  conservative.
- Step 2a': a **deep spray cloud (16 m) with a 4 m pulse** brings the unmerged peak to ~0.34-0.52
  GPa for argon and ~0.57-0.77 GPa for water (converged), at 0.62 / 0.58-0.60 of the ceiling,
  η_jet ~0.50 / ~0.45. A long pulse alone collapses the impulse. Every survivable 1-D
  configuration sits well below the parent's η_jet = 0.775.
- Step 2b (2-D geometry): on the 10 m plate with a 5 m footprint, `eta_capture` falls from 0.79
  (4 m cloud) to 0.30 (16 m). Lateral relief arrives only for clouds deeper than ~8 m. **No cloud
  depth gives both a survivable peak and a useful impulse in this geometry.** The best is argon at
  8 m and 45.58 km/s: 0.64 GPa at η_jet ~0.2. The design point of step 2a' is withdrawn.
- Steps 2c-2d: a 1-2 m standoff cuts the peak 25-45%, and **containing the spill doubles η_jet**.
  A dish (d/D 0.10-0.15) plus a 4 m straight skirt reaches ~0.58 (argon, 45.58 km/s), against
  0.33 open. A flare does no better than a straight skirt. The skirt costs ~20 t per cm of wall.
- Apart from those solves, every number above is a provisional estimate from 1-D ideal-gas or back-of-envelope work. None
  comes from a repository solve, and each must be replaced as the arm runs. The summary for the
  companion repos is [`argon_plate_handoff.md`](../argon_plate_handoff.md).
