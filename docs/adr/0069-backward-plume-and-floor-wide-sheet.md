# ADR-0069: The backward plume is a broad lobe, and the floor-wide sheet keeps Kapton alive in it

Status: proposed, 2026-10-19. Result 18 of [`xray_collision_study.md`](../xray_collision_study.md);
code `sweep --plume-2d`, `plume.py` (`make xray-collision-plume`). Answers the parent's S36 (its
`docs/plume_cadence_asks_for_impact_sim.md`, commit 5bde6a0).

## Context

Each sheet of the 2 Hz, 1 g growth push flies face-on through the gas that earlier pulses threw
back up the stream. The parent's Gaussian cloud left the spread u free. Results 13-16 found that
100 um Kapton holds at u >= 8 km/s and fails at u = 5. ADR-0068 then moved the sheet to cover the
whole 9.5 m floor so the kick stops ringing.

## Decision

1. **The exhaust is a broad lobe, mostly 10-40 deg off the axis at 0.2-0.4 w, with 4-5% within
   5 deg.** It is not a hollow cone. Its dose on the next sheet is computed from ballistic point
   sources of the 2-D exhaust, averaged over the 3 deg cone (the narrowest converged in the grid's
   height).
2. **On the floor-wide sheet the equivalent u is 8.0 km/s at the 9.86 R_sun headline** (38 kJ/m^2,
   face 570 K, thinnest Kapton that holds 61 um against ~70 um flown), and 9.2-11.2 km/s at
   114.8-137 km/s. Kapton holds.
3. **The flown ~8 m sheet is worse, at an equivalent u of 5.4 km/s, and Kapton fails there.** So
   ADR-0068's floor-wide sheet is required twice: for the kick and for the plume.
4. **Keep the skirt; no cone.** Cones of 5-20 deg focus a jet onto the axis and lose 20-40% of the
   impulse. The bare plate is worse on the production grid and loses 15%.
5. **Floor-wide Kapton is the baseline at every printed closure through 137 km/s** (Seth,
   2026-10-19). It holds there on HY-100 or maraging at 2 m (52-70 um against 40-61 um needed and
   75-238 um allowed), and at 208 km/s on maraging only. Its impulse is within 0.2% of the 8 m
   sheet's (-1.2% at 208). **The fallback is a carbon veil spaced off the Kapton** (MLI-style
   scrim and ties), with the Kapton carrying the structure. Glass and bonded foams are rejected.
6. **The plume delays the sheet but does not steer it.** It pushes straight back and evenly
   (0.13-0.54 N s/m^2, peak <= 5.3 Pa). The film slows 3-5 m/s and arrives 40-145 us late,
   predictably. A 100 g node dimples the film under 10 cm at 24-44 N/m of tension (0.5-0.75 MPa,
   ~30-35 rpm).
7. **The plume's heat is kinetic.** Recombination adds ~0.3% at most (atomic carbon). The ions and
   H2 recover before the plate.
8. **No thick coat or bonded foam for heat.** They buy 3-40 K. A 0.1-1 um SiOx-class coat is worth
   flying against particle erosion (~0.04-0.5 um a pass, estimated), not for heat.

## Consequences

- The parent can replace its free u with 8 km/s (9.86 R_sun) to 11 km/s (137 km/s) on the
  floor-wide sheet, and should drop the flown 8 m sheet.
- This reverses ADR-0066's "no Kapton window at 137 km/s" for the floor-wide sheet. Whether the
  4 R_sun headline returns from carbon fabric to Kapton is the parent's call; the margins
  (15-30% in thickness) and the open items below argue for keeping carbon as its fallback.
- The thermal margin on 70 um is thin under the 673 K rule (61 um needed). The decomposition
  kinetics in `docs/kapton_temperature_limits.md` suggest the face can run ~940-1140 K for a
  sub-second pulse. Porting them to a committed module is the next upgrade.
- Open: the axis value inside 1 deg (1.8x the 3 deg value, unconverged), aim scatter, particle
  erosion of the face (~0.04-0.5 um a pass, estimated), and resolution of the bare and cone cases.
