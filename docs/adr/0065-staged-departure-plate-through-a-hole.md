# ADR-0065: The visible departure stages to a second plate, and the impactor threads a hole

Status: proposed, 2026-10-19 (Seth's design, being explored). Result 14 of
[`xray_collision_study.md`](../xray_collision_study.md); code `sweep --hole-plate-2d`,
`departure_hole.py` (`make xray-collision-departure-hole`). Responds to the parent's S38 item 4
(its `docs/departure_occultation_ask_for_aim_repo.md`, commit 8413f3f).

## Context

The flown mirror sign of the 3S departure puts periapsis behind Earth (the parent's S38). The visible
sign meets the stream at 129.9 deg to the thrust axis (117.0 at 4 R_sun), so the impactor arrives from
the ship's side. Result 12's canted flat plate was fitted only near 101 deg, and it lets the fireball
slide past.

## Decision (proposed)

1. **Two plates.** The 12.2 t push plate (ADR-0064) is dropped before the departure.
2. **The departure plate is a bowl with a hole.** Compact impactors thread a ~1 m hole and stop in a
   compact polyethylene charge at the focus of a δ/D 0.5 paraboloid at 1.5 F, as at the Sun. The
   hole slides as the stream turns, and its rim takes extra ablative spray.
3. **1 Hz at ~0.25 g mean**, for aiming through the hole.
4. **k ≈ 20-30.** It maximizes the departing mass per kg of impactor, ~23 at 9.86 R_sun (carbon
   atomic).

## Consequences

- At 130 deg the impactor's own momentum is lost: -11 to -26% of the ceiling per `m w` (k = 30 to 5),
  but only -5 to -28% per kg of impactor at the visible stream's 84 against 75 km/s.
- Isp 781-1093 s at k = 20-10 (carbon atomic), against 443-758 s for the canted flat plate at
  89 deg. The bowl holds the fireball the flat plate let slide.
- The hole takes 0.3-0.9% of the impulse and under 0.03% of the gas.
- The plate is ~10-12 t at k = 20-30, set by its point-blast survival floor. Inertia is 1.6 t.
- The 0.25 g burn (~1900 s) loses ~460 m/s, 11% of the 4.21 km/s burn. The 8-minute burn lost 45.
- Open: the sideways drift (not run), the real collision at 130 deg, the sliding hole and its spray,
  and the apoapsis turn on the visible sign (S38 item 3, the calcs repo's).
