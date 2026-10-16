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
| 9.86 R_sun | 24 | 30 | 7.8-8.0x | 1.1 yr |
| 4 R_sun | 52 | 30 | 10.3-10.6x | 1.0 yr |

*Corrected 2026-10-16* (amendment below): the departure burn is 4.34 / 4.42 km/s and the kept
plate 2.5% of the ship. The first version had 7.0-8.4x and 9.2-10.9x.

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

| closure | excess at arrival -> true mirror excess | burn at 1.1 / 2 / 5 / 6.5 R_J | propellant, share of the ship |
|---|---|---|---|
| 9.86 R_sun | 11.70 -> 15.07 km/s | 0.77 / 1.02 / 1.51 / 1.67 km/s | 0.45 / 0.62 / 0.99 / 1.13% |
| 4 R_sun | 12.17 -> 14.34 km/s | 0.49 / 0.65 / 0.97 / 1.07 km/s | 0.27 / 0.37 / 0.58 / 0.65% |

*Corrected 2026-10-16* (amendment below): the first version burned only to the retrograde floor
(14.78 / 14.16 km/s), and listed 10 R_J, where the flyby cannot turn the disks.

The impactors are `(1 - 1/MR_S)/k_S` ≈ 1/60 of the ship, so the split costs ≲1% of the ship in
propellant. The carrier's dry mass is not counted.

**Rejected alternative.** A separate, faster Earth departure for the impactors (~11.3 km/s excess,
~0.6 km/s more burn than the ship's). That would make the two streams' timing at Jupiter a second
problem to solve.

**Open (ask G10), closed 2026-10-16 by the calcs repo's ADR 0072; see the amendment below.**
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

## Amendment, 2026-10-16: the companions' corrections (calcs ADR 0072, parent E1-E6 and G5-G9)

The round of asks G1-G10 (`docs/solar_environment_4rsun_ask_for_parent.md`) is closed. The
calcs repo answered G1-G4 and G10 in its ADR 0072 (`aim_is_all_you_need` 127162f); the parent
answered E1-E6 and G5-G9 (`docs/solar_environment_answers_from_parent.md`, copied verbatim).

**Corrections taken into code** (`make xray-collision-mission`, `make xray-collision-growth`):
1. **The departure burn is 4.34 km/s (9.86 R_sun) and 4.42 km/s (4 R_sun)**, not 4.15 / 4.23. The
   burn starts at the 10-day orbit's 600 km periapsis, not from parabolic speed at 200 km.
   `mission.departure_dv` now computes it by vis-viva from the departure excess, reproducing the
   calcs repo's figures, so `Closure.departure_dv` is no longer a copied input.
2. **The impactors need the true mirror excess, 15.07 km/s (14.34 at 4 R_sun)**, not the
   retrograde floor. The floor reaches the dive, but not the ship's periapsis at the same instant.
   The flyby turns the disks only below 6.72 R_J (6.90), so the split is listed to 6.5 R_J, where
   the burn is 1.67 km/s, as the calcs repo found. `mission.JUPITER_FLOORS` carries both.
3. **The Sun ship's kept plate is 2.5% of the ship at most** (ADR-0062 amendment), not 0.16.

Net effect: the burn costs ~5% of growth, the plate gives most of the bottom of the range back.
The square charge grows **7.8-8.0x per 3S at 9.86 R_sun** (k_E 24, k_S 30, 21.3 kg departs per kg
of Earth impactor) and **10.3-10.6x at 4 R_sun** (k_E 52, 45.8 kg). The point fireball bounds them
at 10.4-10.6x and 13.4-13.7x.

**Confirmed:** G1 (the closures) and G4 (the 5.3 km/s outbound).

**Not this repository's physics, recorded for the Sun-wave budget** (parent, E2-E3):
- the disks fly edge-on (face normal out of the orbital plane) and flip 90 degrees in the last
  minute; an extreme CME then costs <= 0.6 m/s and ~1 km of drift. Edge-on at 4 R_sun the disk
  runs at ~700 K, not 2400-2900 K;
- guidance electronics >= 100 krad behind ~1 g/cm^2 (a few hundred grams per 1 kg disk), and ~1%
  of passes lost to SEP events at solar maximum;
- dose scales with fluence, r^-0.7 to r^-1.5, not the peak flux's r^-2 to r^-3.
Steady sunlight stays the design heat load at both depths. The shielding mass may change how a
disk strikes its charge if it does not hit like graphite. That is open, not modelled.

**The calcs repo's own ledger** prices the whole cycle with plates, dry mass, the carrier and
launch, counted in launch units: 3.60x per cycle at 9.86 R_sun and 4.92x at 4 R_sun (3.23 and 4.35
per total impactor), against this repository's 7.8-10.6x counted per impactor only. Its optima sit
on this repository's table edges (k_E 5-6, k_S 20); widening them is that repo's ask, not part of
this round.

## Amendment, 2026-10-16: CME density and non-radial flow at depth (parent item 4)

Code `cme_drift.py` (`make xray-collision-cme` -> `cme_{density,tilt,drift,edge_sunlight}.csv`).
Answer: `docs/cme_drift_answers_for_parent.md`, for the parent's 53ca8e8 items 4(a), 4(b) and 5.

- **Density.** Temmer et al. 2021's 29 coronagraph CMEs, scaled inward as r^-3 with 2-6%/R_sun of
  mass still lifting, give a median of 2.7-6.1e5 cm^-3 at 4 R_sun (densest 1.0-2.2e6). The Helios
  in-situ fits, extrapolated, give 1.0e5. The parent's 1e6 / 1e7 assumption is 1.6-10x high
  typical and ~2x high extreme. Peak to mean is ~2: Parker's measured shock compression is 2.03.
- **Tilt.** Inside ejecta it is ~1 degree (1 au). In a sheath it is 23 degrees, which Parker
  measured downstream of the 2022-09-05 shock at 13-15 R_sun.
- **Drift of the edge-on disk.** The CME's own speed is 91-96% of the ram at 1000 km/s. In the
  ejecta even the extreme case costs <= 0.8 m/s and 1.4 km. Only the worst stacked sheath case
  costs ~14 m/s, ~25 km and ~23 g of cold gas at 4 R_sun (2.6 m/s and 12 km at 9.86 R_sun).
- **Edge-on sunlight.** The parent's 11% / 5% of face-on is confirmed (0.108 / 0.048,
  limb-darkened).

Recorded for the Sun-wave budget, not modelled further: a ~1 km tracking requirement and ~25 g of
cold gas per 1 kg disk. Open: the out-of-plane flow inside near-Sun ejecta, and the sheath's
thickness at 4 R_sun (an MHD CME run).
