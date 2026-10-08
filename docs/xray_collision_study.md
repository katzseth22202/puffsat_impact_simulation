# X-ray collision study: the parent's X-ray pusher at 67 km/s

The parent paper's pass-through X-ray pusher (`templateArxiv.tex` §`why_not_xray_drive`,
§`pass_through_pusher`) was written for 618 km/s near-Sun dives, where a carbon-foam collision
reaches about 470 eV. This study asks whether a collision at the rocket's speeds can drive the
same scheme. The impactor arrives at 56 km/s and the ship moves toward it at 11 km/s, so they
close at 67 km/s in the ship's frame. The question is what share of the collision's heat leaves
as light before expansion turns it into bulk motion, and at what photon energy.

This study is separate from `f(v)`, the tamper study, the water and argon plates, and the walled
nozzle. Its Python is `python/puffsat/xray_collision/`, its outputs are under
`data/results/xray_collision/`, and its targets are `make xray-collision`,
`make xray-collision-test` and `make xray-collision-tops`.

## Model

Everything is one-zone, and everything is called through
`fireball.collide(material, impactor_radius, mass_ratio, v_rel, lowering=, radiation=, kappa_scale=, rho0=)`.

- **Kinematics (ship frame).** An impactor of mass `m` hits a resting target of `k m`. Only the
  reduced-mass energy `Q = (1/2) mu v_rel^2` can become heat, which is `k/(1+k)` of the
  impactor's kinetic energy. The rest is centre-of-mass motion that leaves with the debris.
  Every `f_*` share below is a fraction of `Q`. `light_of_ship_ke` is the same quantity as a
  fraction of the impactor's ship-frame kinetic energy.
- **Shocked state.** The same material is on both sides, so the contact moves at `v_rel/2` and
  the first-shocked matter has `u_p = 33.5 km/s`, or `e = u_p^2/2 = 561 MJ/kg`. Its density comes
  from the Rankine-Hugoniot jumps. The fireball starts at that density, at rest in the
  centre-of-mass frame, with the mass-averaged heat `Q / M`.
- **EOS.** Single-element Saha over every stage, with NIST ionization energies and ground-term
  weights (`materials.py`). It is run two ways. `ideal` has no continuum lowering and
  under-ionizes at condensed density. `lowered` subtracts the ion-sphere shift
  `(z+1) e^2/(4 pi eps0 a_i)`. The pair brackets the shocked temperature.
- **Opacity.** TOPS (LANL OPLIB) gray Rosseland and Planck means for pure argon and pure iron,
  over 1e-9 to 20-30 g/cc and 0.5-60 eV (`data/tables/tops/tops_{argon,iron}_fireball_gray.html`).
- **Expansion.** A homologous uniform sphere, with `Rddot = (20 pi/3) R^2 p / M` and
  `dE_int = -p dV - L dt`, integrated with RK4 down to 0.5 eV (the OPLIB floor). Energy closes to
  about 1e-7.
- **Luminosity.** `L = 4 pi R^2 sigma T^4 / (1 + 1/(4/3 kP rho R) + 3/4 kR rho D)`. When thick
  this is the parent's `sigma T^4/(1 + 3/4 tau)` (eq. `cooling_race`). When thin it is volume
  emission `4 kP rho sigma T^4 V`. The spectrum is a blackbody at the zone temperature. That
  overstates its hardness while the zone is opaque, so the photon-band shares (`f_euv` above
  30 eV, `f_xray` above 100 eV) are upper bounds.

## Result 1: two equal 10 kg balls (k = 1) -- `baseline.csv`

| | Argon (liquid, 1395 kg/m^3) | Iron (7874 kg/m^3) |
|---|---|---|
| energy per atom | 232 eV | 325 eV |
| shocked state | 6.5-7.3 t/m^3, kT 27-32 eV, zbar 2.4-3.7 | 32-36 t/m^3, kT 43-52 eV, zbar 2.4-4.0 |
| fireball radius | 9 cm | 5 cm |
| optical depth at start | ~1e6 | ~5e6 |
| parent's cooling race at start | 1.6-3.7e7 | 0.6-1.8e8 |
| **heat radiated** (opacity x1) | **0.8-1.0%** | **3e-7 - 3e-6** |
| above 30 eV / above 100 eV | ~2e-8 / ~1e-8 | ~6e-9 / ~4e-9 |
| luminosity-weighted kT | ~1 eV | ~0.6-0.8 eV |

The verdict holds with the opacity scaled from 0.01x to 10x. Iron gets hotter, as heavier atoms
should, but it radiates 1e3-5e4 times less than argon. It starts denser and smaller, and it stays
opaque all the way down. Argon's 1% is the late glow of cooled gas going transparent near 1 eV.
That is visible and near-UV light, not X-rays.

**At 67 km/s a dense kinetic collision is not an X-ray source.** It is hot only while it is about
a million optical depths thick, and by the time expansion thins it out it has cooled to a few eV.
33.5 km/s gives 5.8 eV per atomic mass unit. Heavier atoms collect more per atom, but they spread
it over more free electrons, and they are more opaque.

## Result 2: a 3:1 target, condensed balls, radius 0.1 um - 10 cm -- `sweep_k3.csv`

A ship-frame target of three impactor masses turns 75% of the impactor's kinetic energy into heat.
The other 25% is centre-of-mass motion, which no geometry can recover. The goal was 80% of that
heat leaving as light. The radii asked for are 1 mm, 1 cm and 10 cm. The sweep runs down to
0.1 um to see whether shrinking ever helps.

| impactor radius | argon `f_rad` | iron `f_rad` | argon / iron kT_emit |
|---|---|---|---|
| 10 cm | 0.55-0.74% | 2e-8 - 4e-7 | 1.0 / 0.8-1.5 eV |
| 1 cm | 0.47-0.63% | 2e-7 - 4e-6 | 1.1 / 0.8-1.5 eV |
| 1 mm | 0.30-0.39% | 2e-6 - 3e-5 | 1.2 / 0.8-1.6 eV |
| 10 um | 0.08-0.11% | 5e-5 - 2e-4 | 3-4 / 4-5 eV |
| 0.1 um | 0.34-0.57% | 0.07-0.19% | 17-20 / 27-34 eV |

(Ranges span the ideal and lowered EOS at opacity x1. The CSV also holds opacity x0.1 and x10.)

**None of the 84 cases reach 80%. The best is 0.74%.** Shrinking a condensed ball does not let its
light out, because the limit is not opacity. Even a fully transparent zone can do no better than a
blackbody, so the race bottoms out at the blackbody ceiling

    t_cool / t_exp  >=  n k T c_s / (2 sigma T^4)

For shocked argon at ~7 t/m^3 and 27 eV this ceiling is ~30-50. At 0.1 um the start race is 50-200,
already near it. Along an adiabat `n/T^3` grows as `rho^-1`, so expansion only raises the ceiling.
At 67 km/s the temperature is pinned at 25-50 eV by the energy per atom, so the only free lever is
the density the heat is deposited at. Argon's ~0.5% shrinks with radius between 10 cm and 10 um
for a related reason. It is the late transparent glow near 1 eV, and the time available for that
glow scales with R.

A lattice of condensed elements therefore cannot beat a condensed ball by surface area alone. Each
element radiates at most at the blackbody ceiling, which is ~30x too slow.

## Result 3: what would reach 80% -- `porous_k3.csv`

The same 3:1 collision, started from foam or spray density `rho0` (ideal Saha, opacity x1). The
pores are assumed to collapse fully onto the plasma EOS.

| rho0 [kg/m^3] | argon `f_rad`, r = 1 mm / 1 cm / 10 cm / 1 m | iron `f_rad`, same radii |
|---|---|---|
| 100 | 0.02 / 0.03 / 0.04 / 0.05 | 0.02 / 0.02 / 0.01 / 0.01 |
| 1 | 0.21 / 0.20 / 0.20 / 0.20 | 0.30 / 0.22 / 0.16 / 0.12 |
| 0.1 | 0.54 / 0.45 / 0.39 / 0.36 | 0.70 / 0.58 / 0.43 / 0.32 |
| 0.01 | 0.78 / **0.82** / 0.74 / 0.63 | 0.79 / **0.91** / **0.84** / 0.68 |
| 0.001 | 0.75 / **0.88** / **0.94** / **0.97** | 0.61 / **0.96** / **0.98** / **0.95** |

80% needs **both** conditions:

- **A low start density.** About 1e-3 to 1e-2 kg/m^3, which is 1e-6 to 1e-5 of the condensed
  density and roughly the density of a gas at a few hundred pascals. This drops the race below 1.
- **An optical depth near one.** The pre-collision column `rho0 r` must be about 1e-4 to
  1e-3 kg/m^2. Thicker, and the light is trapped. Much thinner (argon or iron at 1e-3 kg/m^3 and
  1 mm), and the zone is too dilute to emit before it flies apart.

The passing bodies in this grid hold 4 ug to 4 g. That is a property of the grid, not a limit.
What has to stay small is the column, not the mass. With the column held fixed and the cloud made
larger and thinner, the share only rises:

| column [kg/m^2] | r = 1 cm | 10 cm | 1 m | 10 m | 50 m |
|---|---|---|---|---|---|
| 1e-4, argon / iron | 0.82 / 0.91 | 0.94 / 0.98 | 0.98 / 0.99 | 0.99 / 0.98 | 0.99 / 0.98 |
| 1e-3, argon / iron | 0.45 / 0.58 | 0.74 / 0.84 | 0.97 / 0.95 | 0.99 / 0.99 | 0.99 / 0.99 |

At 1e-3 kg/m^2, a 50 m cloud (2e-5 kg/m^3) holds a 10.5 kg impactor and radiates 99% of its heat.
The cost is softer light: kT_emit falls from ~5 eV to ~3 eV, and the share above 30 eV from ~15%
to 3-6%. That is because more of the heat goes into ionization at low density. So the
requirement is a tenuous cloud, or a sheet of about 1e-3 kg/m^2, not a small piece of foam. Those
densities sit far below any solid foam (the lightest aerogels are ~0.2-1 kg/m^3). The medium is a
gas, vapour or dust cloud, or a lattice whose mean density is that low.

**Even then the light is EUV, not X-ray.** Low density moves heat into ionization, so the zone
starts at only 6-8 eV and emits at kT 4-6 eV (peak photon ~12-17 eV). 10-20% of the heat comes
out above 30 eV, and 4e-6 to 2e-4 above 100 eV. At 67 km/s any radiation-driven pusher is an
EUV pusher. EUV is absorbed in tens of nanometres of almost anything, which suits ablation, but
the radiation case and channel filler of the parent's design would have to be re-chosen for it.

**Open for the lattice study.** A sheet of condensed elements may behave like foam at its mean
density, if each element's plume merges with its neighbours' before it radiates. That
re-thermalization is the mechanism that would make a lattice work, and this one-zone model cannot
see it. Collisionality also has to be checked. At 1e-2 kg/m^3 the two streams may interpenetrate
rather than shock.

## Result 4: the parent's 618 km/s carbon case -- `near_sun.csv`

This is the parent's pass-through pusher as written: a 1 kg carbon-foam projectile meets a 1 kg
carbon-foam target (k = 1) at 618 km/s. The material is carbon (NIST ladder, TOPS pull to 3 keV).
The rows bracket trapped radiation in the EOS (`a T^4` in or out) and scale the opacity by 0.1-10.
The two `>100 eV` and `kT_emit` values are the photosphere and zone bounds.

| foam rho0 | kT0 | tau0 (x1) | f_rad, opacity x0.1 / x1 / x10 | >100 eV (x1) | >1 keV (x1) | kT_emit (x1) |
|---|---|---|---|---|---|---|
| 10 kg/m^3 | 390-470 eV | 0.9-1.4 | 0.70 / **0.59-0.61** / 0.45-0.48 | 0.53-0.57 | 0.25-0.28 | 180-280 eV |
| 30 kg/m^3 | 430-470 eV | 2.9-3.9 | 0.57 / 0.43-0.45 / 0.26-0.29 | 0.39-0.43 | 0.23-0.25 | 180-320 eV |
| 100 kg/m^3 | 450-470 eV | 13-15 | 0.39 / 0.22-0.23 / 0.07-0.08 | 0.20-0.22 | 0.13-0.14 | 160-355 eV |
| graphite, 2200 | 490 eV | ~860 | 0.01 / 0.004 / 0.003 | 0.001 | 0.001 | 17-86 eV |

**The parent's assumption holds for its stated design.** At 10 kg/m^3 the model reproduces the
parent's 470 eV (468 eV) and its assumed 60% radiated (59-61%). More than 90% of that light is
above 100 eV, and about 45% is above 1 keV. These are soft X-rays, as the parent says. The margin
is thin, though. 30 kg/m^3 gives 43-45%, 100 kg/m^3 gives 22-23%, and solid graphite gives 0.4%.
Ten times TOPS's opacity drops even the 10 kg/m^3 case to 45-48%. So "the lightest foam that still
holds together" is a requirement, not a preference. Trapped radiation (`a T^4` is 14-30% of the
matter energy at 10 kg/m^3) moves the result by only 1-2 points.

**The opacity is not what the parent guessed, and the mechanism differs.** TOPS gives fully stripped
carbon at the start state (46 kg/m^3, 468 eV) kR = 0.04 and kP = 0.18 m^2/kg, which is 25-250x below
the parent's 1-10 m^2/kg bracket. The fireball starts about one optical depth thick. Bremsstrahlung
from a thin plasma is weak, so the race is closer than the parent's `sigma T^4/(1 + 3/4 tau)`
suggests. The parent's formula saturates at a full blackbody when thin, and its race of
0.004-0.009 would imply that nearly everything radiates. Most of the light comes as the fireball
cools through ~100-300 eV and K-shell recombination switches on. At 200 eV kP = 9 m^2/kg, which is
where the parent's bracket is right.

**Not checked here.** The shock heats ions first. Electrons are heated by electron-ion equilibration,
estimated at ~50 ns at 470 eV and 1.4e22 cm^-3, against ~1 us for expansion. That is fast enough
to keep the radiated share, but electrons running cooler than ions would soften the spectrum.
Non-LTE emission at these densities and the collisionality of foam on foam at 309 km/s each
(parent §needle_through_fog) are also open.

## Result 5: a lattice target and the mass ratio at 618 km/s -- `near_sun_target.csv`

Suppose the projectile has to be dense to survive its trip, while the ship's target can be a
carbon lattice deployed briefly. `fireball.piston` models a dense 1 kg graphite projectile
snowplowing `k` kg of a target at mean density `rho0`. In a snowplow each swept kilogram absorbs
`u^2/2` at the current speed, so **all of the heat Q lands in the light target**. The projectile
is only compressed by the ram pressure (~4 Mbar at 10 kg/m^3), which this model neglects. Rows
marked `both` are the foam-on-foam collision for comparison.

| k | rho0 | foam on foam: f_rad / light per KE | dense into lattice: f_rad / light per KE | column length (lattice) |
|---|---|---|---|---|
| 1 | 1 | 0.86 / 0.43 | 0.94 / 0.47 | 140 m |
| 1 | 10 | 0.61 / 0.30 | **0.81 / 0.41** | 14 m |
| 1 | 100 | 0.23 / 0.12 | 0.62 / 0.31 | 1.4 m |
| 2 | 10 | 0.55 / 0.37 | 0.71 / 0.47 | 28 m |
| 3 | 10 | 0.46 / 0.35 | 0.60 / 0.45 | 42 m |
| 3 | 1 | 0.81 / 0.61 | 0.86 / 0.64 | 420 m |

The column length is the target length a 4.8 cm graphite ball must sweep: `k m / (rho0 pi r^2)`.

- **A dense projectile helps rather than hurts.** The heat sits in the dilute target alone, which
  runs at 1030 eV instead of 470 eV (k = 1). At 10 kg/m^3 the radiated share rises from 61% to 81%,
  with 60% of the heat above 1 keV.
- **Mean density matters, not surface area.** A lattice behaves as a foam at its mean density once
  the shock merges its struts. Shrinking a 10 kg/m^3 foam fireball from 1 kg to 0.01 kg (more
  surface per mass) only moves 61% to 68%. Going from 10 to 1 kg/m^3 moves it to 86%.
- **3:1 trades target mass for light.** Light per joule of projectile kinetic energy peaks near
  k = 2 (0.47 for the dense-into-10 kg/m^3 case), against 0.41 at k = 1 and 0.45 at k = 3. Light
  per kilogram of ship-carried target falls from 78 GJ (k = 1) to 45 GJ (k = 2) and 29 GJ (k = 3).
  Higher k also cools the fireball (470 eV at k = 3) and thickens it. The parent counts the target
  as spent mass, so k = 1 to 2 is the likely optimum for specific impulse. The pass-through
  specific-impulse optimizer should decide it. This study does not.
- **Solar survival (hand estimate).** At 4 solar radii (3.93 MW/m^2), an isolated gray strut that
  radiates from both sides settles at (S / 2 sigma)^(1/4) ~ 2400 K. Carbon's sublimation is slow
  there, so a lattice should survive far longer than a second as long as its struts can see cold
  sky.

## Limits

One zone, so there is no cool photosphere and no temperature gradient. For k > 1 the shock decays
into the target, and the one-zone mean heat stands in for a graded heating. A collision of two balls
is a pancake, not a sphere. LTE is assumed throughout, including the thin, cold tail. The
lowered EOS is a bracket, not a model.
