# ADR-0063: The 3S departure uses a flat canted plate, and k is optimized at both waves for growth per total impactor

Status: accepted, 2026-10-15. Result 12 of [`xray_collision_study.md`](../xray_collision_study.md);
code `euler2d::canted`, `sweep --canted-2d`, `growth.py` (`make xray-collision-growth`). Builds on
[ADR-0062](0062-sun-wave-redirects-earth-wave-breeds.md). Asks G1-G9 in
[`solar_environment_4rsun_ask_for_parent.md`](../solar_environment_4rsun_ask_for_parent.md).

## Context

Seth wants the Earth wave to drive the calcs repo's 3S cycle (Earth -> Jupiter -> Sun -> Earth),
and k chosen at each wave to maximize growth per 3S. Two reasons favour the Sun dive:
- it avoids head-on rockets at the departure;
- the solar boost has energy and the Oberth effect to spare.

The calcs repo's own closure puts the departure's thrust axis at ~101° to the returning stream.
So the impactor sweeps past the plate and cannot rebound off it.

## Decision

1. **A flat canted plate.** It takes no tangential force, so the stream's sideways momentum leaves
   with the gas. A cant that turned the stream into the plate (θ < 90°) was considered and
   rejected. It gains ~0-5%, and any sheet that misses the charge would then strike the plate at
   ~100 km/s.
2. **Priced against the canted ceiling** `m w (cos theta + sqrt(1 + k - sin^2 theta))`, in planar
   2-D strips. A wide sheet strikes a polyethylene charge held 1 m off a 20 m flat plate. The real
   collision is the design value; the point fireball is the ideal it could approach.
3. **The square block is the charge.** Slab, round and half-cylinder charges were tried, and all
   did worse (Seth asked about hemispheres).
4. **Growth counts all impactor.** The optimizer maximizes the Earth-wave impactor returned per kg
   of Earth- *and* Sun-wave impactor spent, over `k_E` in 5-60 and `k_S` in 20-100. The calcs
   repo's closures are copied inputs (commit 868e88a); this repository cannot import that code.
5. **The closures are 9.86 R_sun / 100 km/s (baseline) and 4 R_sun / 150 km/s (upgrade).**

## Consequences

| closure | `k_E` | `k_S` | growth per 3S (square charge; plate dead to released) | doubling |
|---|---|---|---|---|
| 9.86 R_sun | 24 | 30 | 7.0-8.4x | 1.1-1.2 yr |
| 4 R_sun | 54 | 30 | 9.2-10.9x | 1.0 yr |

The Sun stage is priced on the real disk-into-charge collision (ADR-0060 amendment). That holds
`k_S` near 30 in both cycles.

- The canted plate keeps ~1/3 of its ceiling (0.29-0.37), against 0.70-0.84 for head-on plates.
  It is the weak link. The ideal point fireball's 0.43 bounds what a better charge could recover:
  ~+30% growth.
- The optimum k is set by polyethylene's bond toll at the Earth wave and by the real collision's
  coupling at the Sun. If the carbon condenses, or a bond-free charge is carried, growth rises
  10-27%, with `k_E` at the edge of the simulated range.
- **Defect found on the way.** All disk-into-charge runs had started the charge at the disk's
  temperature, not its pressure, giving it spurious heat (+10% at k = 10, +90% at k = 100). Both
  runners now start cold and assert it. Results 10-11's affected rows were rerun (see ADR-0060 and
  ADR-0061, amended 2026-10-15).
