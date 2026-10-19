# The push plate's backward plume: answers (S36)

Answers 2026-10-19 to the parent's S36 (its `docs/plume_cadence_asks_for_impact_sim.md`). **Written
to be copied into the parent.** Evidence: Result 18 of
[`xray_collision_study.md`](xray_collision_study.md) and
[ADR-0069](adr/0069-backward-plume-and-floor-wide-sheet.md); `make xray-collision-plume`
(`plume_*.csv`). **Status: proposed.** It assumes ADR-0068's floor-wide sheet (the sheet covers
the 9.5 m floor; the skirt rides the frame).

## Short answers

1. **The exhaust is a broad lobe.** Most of it leaves 10-40 deg off the stream axis at 0.2-0.4 of
   the closing speed. 4-5% of the mass leaves within 5 deg.
2. **Your u, fitted: 8.0 km/s at the 9.86 R_sun headline on the floor-wide sheet**, and
   9.2-11.2 km/s at 115-137 km/s. 100 um Kapton holds: the face peaks at 570 K, and the thinnest
   film that holds is 61 um. **The flown ~8 m sheet fits only u ~ 5.4 km/s, and Kapton fails
   there.** A sheet narrower than the floor focuses its own rebound onto the axis.
3. **The gas is cold and neutral by the time the sheet meets it.** The ions recombine and H2
   re-forms before the plate. Even if every carbon atom stayed atomic and recombined on the face,
   it would add ~0.3% to the dose. The dose is kinetic, ~50 eV/amu.
4. **Kapton, as wide as the floor, holds at every printed closure through 137 km/s**, and at
   208 km/s on maraging only. The wider sheet gives the same impulse within 0.2%. This reverses
   our earlier "no Kapton window at 137 km/s" (Result 15), which assumed your Gaussian at u = 8.
5. **The plume does not blow the sheet off course.** It pushes straight back and evenly across the
   sheet. The film slows 3-5 m/s and arrives 40-145 us late, predictably.
6. **Canting or a cone does not clear the axis.** A convex cone throws the main lobe off axis but
   focuses a jet onto it, and loses 20-40% of the impulse. A 10-20 deg tilt will not save thin
   sheets at 4 Hz.

## The exhaust (item 1)

Mass share by angle from the axis and terminal speed (w = closing speed), on the flat skirted plate
under the floor-wide sheet, k = 7.07:

| angle | < 0.2 w | 0.2-0.3 | 0.3-0.4 | 0.4-0.5 | > 0.5 w |
|---|---|---|---|---|---|
| 0-5 deg | 0.001 | 0.005 | 0.013 | 0.007 | 0.014 |
| 5-10 | 0.000 | 0.021 | 0.034 | 0.019 | 0.027 |
| 10-20 | 0.000 | 0.094 | 0.111 | 0.061 | 0.026 |
| 20-40 | 0.000 | 0.266 | 0.087 | 0.011 | 0.003 |
| 40-90 | 0.001 | 0.150 | 0.036 | 0.001 | 0.001 |

The grid reaches 50 m above the plate. A short grid miscounts hot gas that has not finished
spreading as moving straight up the axis.

## The dose and your u (item 2)

Each earlier pulse is treated as a ballistic point source of this exhaust; beyond ~1 km the gas is
collisionless. The dose depends only on the density right on the axis, because the next sheet
subtends ~0.02 deg where the dose arrives. Fed your Gaussian, the method reproduces your dose to
2%. The axis value is averaged over a 3 deg cone, the narrowest that has converged.

| closure | flown ~8 m sheet | floor-wide 9.5 m sheet |
|---|---|---|
| 9.86 R_sun, 60 (99.2 km/s) | 83 kJ/m^2, u 5.4, fails | **38 kJ/m^2, u 8.0, 570 K, holds** |
| 9.86 R_sun, X 75 (114.8) | u 6.2, fails | u 9.2, 540 K, holds |
| 4 R_sun, X 90, 5-day gap (119.0) | u 6.6, fails | u 9.8, 506 K, holds |
| 4 R_sun, X 75 (121.8) | u 6.7, fails | u 10.0, 507 K, holds |
| 4 R_sun, 90 (137.1) | u 7.5, fails | u 11.2, 485 K, holds |

"Holds" uses your 400 C rating for 100 um Kapton.

**Fly the sheet as wide as the floor.** It keeps the same impulse per kg, and the thinnest film
that holds is 40-61 um against the ~70 um that 28 kg spreads to. At the 4 R_sun closures this
matters only for Kapton; carbon fabric does not decompose.

## Species (item 3)

| store | at the plate |
|---|---|
| ionization | 4e-6 of the heat still locked: recombined |
| H2 | re-formed, 3.1 decades of margin |
| carbon | bracketed. If atomic, it carries 53 MJ/kg against 18 GJ/kg of kinetic energy in the sheet's frame: ~0.3% of the dose |

## Canting or curving the plate (item 4)

| plate (9.5 m sheet) | dose, 3 deg cone | impulse vs flat skirted |
|---|---|---|
| flat + skirt | 38 kJ/m^2 | 1 |
| flat bare | 63 | 0.85 |
| cone 5 / 10 / 15 / 20 deg | 92 / 119 / 137 / 151 | 0.80 / 0.73 / 0.66 / 0.60 |

The bare and cone doses are not converged in resolution, but every grid ranks the cones worse.

## Kapton at every closure, and the push on the sheet

The pulse's impactor mass spread over the floor-wide sheet, against the thinnest Kapton that holds
in this plume and the face's pressure ceiling at 2 m:

| closure | floor-wide sheet | least that holds | ceiling HY-100 / maraging | impulse vs 8 m sheet |
|---|---|---|---|---|
| 9.86 R_sun, 60 (99.2 km/s) | 70 um | 61 | 143 / 238 | +0.2% |
| 115-122 km/s band | 59-61 | 45-53 | 95-107 / 158-178 | +0.2% |
| 9.86 R_sun, 90 (132.1) | 52 | 48 | 80 / 134 | +0.2% |
| 4 R_sun, 90 (137.1) | 52 | 40 | 75 / 124 | -0.1% |
| 4 R_sun, 150 (207.8) | 33 | 28 | 33 / 54 | -1.2% |

| closure | the plume's push | peak | film slows | arrives |
|---|---|---|---|---|
| 9.86 R_sun headline | 0.54 N s/m^2 | 5.3 Pa | 5.4 m/s | 145 us late |
| 115-137 km/s | 0.27-0.41 | 2.7-4.1 | 3.7-4.8 | 72-113 us |
| 208 km/s | 0.13 | 1.3 | 2.9 | 38 us |

The sheet spans only ~0.02 deg of the plume, so there is no side force or tilt. A heavy node on the
film is held by the film's tension: 100 g dimples it under 10 cm at 24-44 N/m (a ~30-35 rpm
spin, 0.5-0.75 MPa in the film).

## Other ways to protect the sheet

- **A bonded coat or foam does not help.** A 1 um graphite coat lowers the face 3-12 K. A bonded
  20-50 um carbon foam lowers it 15-40 K, because porous carbon conducts well when hot.
- **A carbon veil held off the film by a gap** (MLI-style scrim and ties, ~2-12 g/m^2 of veil or
  CNT sheet) is the fallback at higher speeds or if the plume proves worse, with the Kapton
  carrying the structure.
- **A 50 um glass top layer does not work.** It weighs more than the whole Kapton sheet, and thin
  glass is brittle over a 19 m spun sheet.
- **A 0.1-1 um SiOx-class coat**, the standard atomic-oxygen coat on Kapton, is worth flying
  against particle erosion (estimated ~0.04-0.5 um a pass), not for heat.

## Still open

- The density within 1 deg of the axis is 1.8x the 3 deg value and has not converged.
- Aim scatter, which would smear the axis peak.
- Particle erosion of the face: ~0.04-0.5 um a pass, estimated in
  `docs/kapton_temperature_limits.md`, not modelled.
- The 400 C rating may be very conservative for a sub-second pulse. The decomposition kinetics in
  that note put the limit near 940-1140 K, but they are not yet in a committed model.

## Suggested changes (for the parent to make)

- `sec:charge_pusher_plume`, the sheet: "The baseline sheet is Kapton, as wide as the push plate,
  about 70 um at Parker's depth. It holds in the computed plume at every closure through
  137 km/s, at the same impulse as a narrower sheet. If the plume proves hotter, or at higher
  speed, a carbon veil held off the film by a scrim takes the heat and the Kapton carries the
  structure \cite{Katz_puffsat_impact_sim_2026}." Whether the 4 R_sun headline goes back to
  Kapton is your call; we would keep carbon fabric as its fallback.
- `sec:charge_pusher_plume`, after the sheet: "The plume pushes the sheet straight back and evenly,
  slowing it 3 to 5 m/s, so it arrives up to 0.15 ms late but on course."
- `sec:charge_pusher_plume`: replace the free u with "the impact simulation's exhaust gives an
  effective u of 8 km/s at 99 km/s, if the sheet covers the plate
  \cite{Katz_puffsat_impact_sim_2026}". Print the dose table at u = 8 as the computed case and u = 5
  as the narrow-sheet case.
- State that the sheet is as wide as the plate (ADR-0068), and that canting the plate was checked
  and rejected.
