# ADR-0064: At 1 Hz the flat plate is set by its absorber, and the sheet's thickness by its face

Status: accepted, 2026-10-18. Result 13 of [`xray_collision_study.md`](../xray_collision_study.md);
code `sweep --sheet-plate-2d`, `cadence_plate.py` (`make xray-collision-cadence-plate`). Answers
the parent's S37.3 and supplies S34 inputs (its `docs/plume_cadence_asks_for_impact_sim.md` and
`docs/plume_cadence_asks_for_aim_repo.md`, commit 5bde6a0).

## Context

The parent slows the 400 km overtake push to 1 Hz at 14.7 MN s (~1 g) so that the 100 km/s Kapton
sheets survive the earlier pulses' gas, and offers a 1 mm sheet as the alternative. Seth asked for
the plate at that force, whether a flat sheet's even load lightens it, and whether thicker sheets
on a smaller plate pay.

## Decision

1. **The plate is sized by the inertia bound** `I^2 / (2 s M a)` at the calcs repo's 100 g cap and
   2 m stroke, with `M` the unit at the end of the push. The planar sheet's even load pushes the
   plate as one body, so the calcs repo's point-blast survival floor (30 kg/m^2 scaled as
   energy^(2/3)) is not used. At 1 Hz the bound is the larger anyway: 48.8 t against 46.5 t.
2. **The sheet's peak is `c(h) sigma v^2`**, fitted from 2-D runs at the flown k = 6.74:
   c = 0.65 / 0.50 / 0.35 at 2 / 4 / 8 m standoff (production grid; a lower value). The peak does
   not depend on the sheet's radius.
3. **The face is solid hard steel**, since a ~49 t plate is 4-43 mm of steel over any 5-20 m
   plate. Its unprestressed allowables (1.5 GPa HY-100 / 4340, 2.5 GPa maraging 300) replace the
   400 MPa sandwich facesheet, and admit 140-250 um sheets at a 2 m standoff.
4. **A small plate keeps a skirt.** Sheet radius + 1.5 m with a 2 m skirt keeps 94-97% of the
   flown plate's impulse at 2 m standoff; bare it keeps 72-82%.
5. **The 1 mm sheet is rejected** at any standoff run (5-9 GPa).

## Consequences

- The plate is not lighter at 1 Hz; it is ~49 t at 2 m stroke, 24 t at 4 m and 10 t at 10 m. A
  smaller plate only makes it thicker.
- The sheet can be several times thicker than the flown 35-62 um. That is what the plume needs.
- Open: the peak grows on refinement (~37% from coarse to production); the face's rise time against
  its 4-15 us ring period; the plume's spread `u` (S36), which sets the sheet's floor.
