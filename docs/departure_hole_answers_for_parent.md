# The visible departure on a second plate with a hole: answers for `Balloon-Pulse-Propulsion`

2026-10-19. **Written to be copied into the parent.** It answers S38 item 4 (the canted law at the new
angle, asked of the calcs repo and passed on here) and gives the calcs repo inputs for items 2 and 5.
Evidence: Result 14 of [`xray_collision_study.md`](xray_collision_study.md),
[ADR-0065](adr/0065-staged-departure-plate-through-a-hole.md), `make xray-collision-departure-hole`
(`departure_hole_law.csv`, `_loss.csv`, `_burn.csv`, `_stage.csv`). The routes' angles, closing
speeds and burns are your S38 table's.

## The design (Seth's, proposed)

The visible sign's stream arrives at 129.9 deg (117.0 at 4 R_sun), from the ship's side.
- **The plates are staged.** The 12.2 t push plate is dropped.
- **The departure plate is a bowl with a hole.** Compact impactors thread a ~1 m hole and stop in a
  compact polyethylene charge at its focus (δ/D 0.5 paraboloid, charge at 1.5 F), as at the Sun.
- **Operations:** the hole slides as the stream turns over the burn, and its rim takes extra spray.
- **1 Hz at ~0.25 g mean.**

## Answers

1. **Impulse lost to the angle.** The lost share is the impactor's own momentum,
   `cos theta m w`.
   - **Per `m w`:** against 90 deg, the ceiling falls 11-26% at k = 30-5, and the solved bowl
     15-35% against its own 89.2 deg runs.
   - **Per kg of impactor:** the visible stream closes faster (84 against 75 km/s), so the loss is
     only 5% at k = 30, 9% at 20, 17% at 10 and 28% at 5.
2. **Effective Isp, 9.86 R_sun, carbon atomic (condensed):**

   | k | canted flat plate, 89.2 deg (flown) | hole bowl, 129.9 deg (visible) |
   |---|---|---|
   | 5 | 1185 s (1258) | 1387 s (1454) |
   | 10 | 758 s (851) | 1093 s (1196) |
   | 20 | 443 s (571) | 781 s (944) |
   | 30 | 290 s (451) | 589 s (806) |

   At 4 R_sun (117 deg, 107.9 km/s) the hole bowl gives 2071 / 1569 / 1149 / 921 s. So the visible
   sign, on the bowl, is *better* than the flown sign was on the flat plate. Use these in place of
   the held-flat `earth_ve` extrapolation.
3. **The hole loses little.** At 130 deg the hole sits 3.2 m off the fireball's axis. There a 1 m
   hole takes 0.5-0.6% of the impulse and ~0.01% of the gas; across four bowls, 0.3-0.9% and under
   0.03%. The Sun plate's on-axis bore takes 0.5-0.6% and 0.05%.
4. **Large k works up to ~30.** Above that, carbon's bonds at 84 km/s eat the heat: 27% of the speed
   at k = 30 if carbon stays atomic. The departing mass per kg of impactor peaks at k ≈ 20-30: ~23
   (carbon atomic) to 30-37 (condensed) at 9.86 R_sun, and 35-50 at 4 R_sun.
5. **The plate is ~10-12 t at k = 20-30** (15 t at k = 10, 21 t at k = 5). It is set by its
   point-blast survival floor, the calcs repo's paraboloid law, at 10-20 kg of impactor per 1 Hz pulse.
   The absorber's inertia bound is only 1.6 t at 0.25 g.
6. **The 0.25 g burn costs ~460 m/s** (11% of 4.21 km/s) over ~1900 s, fixed direction, on the
   10-day orbit. The calcs repo's `finite_burn_loss` model was re-implemented and checked against its
   pins. 1600 s costs 360 m/s, and the flown 8 minutes cost 45. Doubling the mean to 0.5 g roughly
   halves the burn and cuts the loss to ~170 m/s (the 1000 s row). At 1 Hz each pulse then doubles, and
   the survival floor grows as pulse energy^(2/3); at 2 Hz it would not.

## For the calcs repo (S38 items 2 and 5)

- Fly the visible sign on the hole bowl's law (`departure_hole_law.csv`) rather than the held-flat
  canted law.
- Stage the plates: drop the 12.2 t push plate and carry a ~10-12 t departure bowl.
- Charge the 0.25 g burn's ~460 m/s, or pick the mean thrust on the trade.

## Not modelled

- The fireball's sideways drift (the runs are axisymmetric).
- The real compact-into-charge collision at 130 deg (the Sun plate's factor stands in).
- The sliding hole's travel and mechanism, the oblique hole's shape and the collar's spray.
