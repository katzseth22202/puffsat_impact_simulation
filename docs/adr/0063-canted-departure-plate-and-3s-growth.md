# ADR-0063: The 3S departure uses a flat canted plate, and k is optimized at both waves for growth per total impactor

Status: accepted, 2026-10-15. Result 12 of [`xray_collision_study.md`](../xray_collision_study.md);
code `euler2d::canted`, `sweep --canted-2d`, `growth.py` (`make xray-collision-growth`). Builds on
[ADR-0062](0062-sun-wave-redirects-earth-wave-breeds.md). Asks G1-G10 in
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
   did worse (Seth asked about hemispheres). *Amended 2026-10-15:* two half-cylinder runs were
   mis-aimed (see below); with the aim fixed, the dome-to-the-plate half-cylinder edges the
   square by 1-6%. The square stays the design charge.
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

## Amendment, 2026-10-15: the impactor stream splits from the ship at Jupiter

Code `mission.jupiter_split` (`make xray-collision-mission` -> `mission_jupiter_split.csv`).

**Decision.** The Sun-wave impactors and the Sun ship leave Earth together and split at Jupiter.
The ship falls into a *prograde* dive and rides the cycle unpowered. The impactors fall into the
mirror-image *retrograde* dive. Mirror-image dives that leave Jupiter together reach periapsis
together, so the two streams meet head-on there.

A flyby rotates the Jupiter-relative excess but never changes its size. So the impactors need a
perijove burn that raises it to the retrograde floor, and that burn gets Jupiter's Oberth boost.
The burn is done by a carrier bus with a storable stage (330 s, stated not derived). It takes all
the disks through the flyby, burns once at perijove, and then releases them. The calcs repo's
floors and arrival excesses are copied as cited inputs (`aim_is_all_you_need` 868e88a,
`opposing_stream_placement`):

| closure | excess at arrival -> retrograde floor | burn at 1.1 / 2 / 5 / 10 R_J | propellant, share of the ship |
|---|---|---|---|
| 9.86 R_sun | 11.70 -> 14.78 km/s | 0.70 / 0.92 / 1.37 / 1.77 km/s | 0.40 / 0.55 / 0.88 / 1.21% |
| 4 R_sun | 12.17 -> 14.16 km/s | 0.45 / 0.59 / 0.88 / 1.14 km/s | 0.25 / 0.34 / 0.52 / 0.70% |

The impactors are `(1 - 1/MR_S)/k_S` ≈ 1/60 of the ship, so the split costs ≲1% of the ship in
propellant. The carrier's dry mass is not counted.

**Rejected alternative.** A separate, faster Earth departure for the impactors (~11.3 km/s excess,
~0.6 km/s more burn than the ship's). That would make the two streams' timing at Jupiter a second
problem to solve.

**Open (ask G10).**
- **Radiation.** A 1.1 R_J perijove is severe for the disks' guidance electronics. Higher
  perijoves are gentler but cost about twice the burn.
- **Bend.** At the higher excess the flyby bends the path less, and only the ship's bend is
  checked (135° available against 77° needed).
- **Timing.** The streams must leave Jupiter as true mirror images.

Bend and timing are the calcs repo's flyby machinery, so they are asked of it, not computed here.

## Amendment, 2026-10-15: the half-cylinders' aim, and a methane-faced charge

Code `euler2d::canted` (`PlateLayer`), `sweep --canted-2d shapes` and `--canted-2d layers`
(`make xray-collision-growth`).

**Defect.** The sheet was always 1 m wide about the charge's centre. Half-cylinders lying on a face
(dome up, dome to the plate) fill only half of that across the stream, so ~41% of the sheet missed
them. The sheet now spans the across-stream extent of what it is aimed at, and full-width charges
are unchanged. Corrected shares at k = 10 / 30 / 60: dome to the plate 0.38 / 0.33 / 0.30 (was
0.28 / 0.24 / 0.22), dome up 0.32 / 0.27 / 0.23 (was 0.23 / 0.19 / 0.17). The square block is
0.36 / 0.32 / 0.29.

**Methane-faced charge.** A dense plate-side layer (solid methane) under polyethylene foam, with
the sheet striking the foam alone. The best is the dome-to-the-plate half-cylinder with 25-50% of
k in a layer 5x the foam's density: 0.397-0.399 at k = 10 and 0.337-0.346 at k = 30, about +10% and
+7% over the square block. A thin layer at 15x does nothing. Both layers are one gamma gas, so
methane's EOS, rebound and bonds are not modelled.

**Not yet decided.** The design charge stays the square block until the growth optimizer is rerun
with the layered charge and methane's rebound and bonds are priced.

