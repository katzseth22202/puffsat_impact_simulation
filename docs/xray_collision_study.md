# X-ray collision study: the parent's X-ray pusher at 67 km/s

The parent paper's pass-through X-ray pusher (`templateArxiv.tex` §`why_not_xray_drive`,
§`pass_through_pusher`) was written for 618 km/s near-Sun dives, where a carbon-foam collision
reaches about 470 eV. This study asks whether a collision at the rocket's speeds can drive the
same scheme. The impactor arrives at 56 km/s and the ship moves toward it at 11 km/s, so they
close at 67 km/s in the ship's frame. The question is what share of the collision's heat leaves
as light before expansion turns it into bulk motion, and at what photon energy.

This study is separate from `f(v)`, the tamper study, the water and argon plates, and the walled
nozzle. Its Python is `python/puffsat/xray_collision/`, its outputs are under
`data/results/xray_collision/`. Its targets are `make xray-collision`, `make xray-collision-test`,
`make xray-collision-tops`, `-argon-bags`, `-anvil`, `-ch2-table`, `-rebound`, `-charge-2d`,
`-charge`, `-charge-mix` and `-freeze`. The study grew past its title. Results 7-10 ask what
pushes a ship at a Parker-like periapsis (400 km/s closing).

## Conclusions

1. **At rocket speeds (66-67 km/s) a collision is not an X-ray source** (Results 1-3, 6). Dense
   matter at a few tens of eV stays opaque until it has cooled to ~1 eV. Light carries most of the
   heat only from gas spread to ~1 g/m^3, and a cloud that thin cannot thread the bore.
2. **At 618 km/s the parent's X-ray pusher works** (Results 4-5). A dense disk into a ~10 kg/m^3
   carbon lattice radiates 47% of the projectile's energy as soft X-rays. The light partner of a
   collision takes the heat, and that rule recurs throughout.
3. **At 400 km/s the X-ray route survives only with the density contrast reversed** (Results 7-9,
   ADR-0059). A light foam projectile stops on a dense graphite anvil inside the case. Iron does
   not beat carbon, and a z-pinch adds lossy stages. The design clears the parent's 2650 s bar
   only because the polyethylene propellant rebounds from the plate (`e_eff` 0.63-0.72, computed),
   and its case-wall mass is unresolved.
4. **The simplest near-Sun pusher has no X-rays and no case** (Result 10, ADR-0060).
   - **Design:** a graphite disk threads the bore and stops in a compact charge of near-solid
     polyethylene, 20-100 times its mass. The fireball is too dense to radiate, so all of its heat
     becomes expansion onto the plate.
   - **Mission:** on the parent's own 20 m bowl, a Jupiter fall boosted at 9.86 R_sun crosses
     Earth's orbit at 100 km/s with a mass ratio of 2.0-2.9. A δ/D 0.5 paraboloid lowers that to
     1.7-2.3.
   - **Cost:** projectiles come to 1-3% of the final mass.
   - **Checks:** ionization recombines and polyethylene's bonds cost 1.6-5.4% of velocity. A solid
     charge mixes well enough, and foam does not.
   - **The one hazard** is a back-splash jet of several GPa onto the ring around the bore.
5. **Without a nozzle, a plate recovers less directed energy per joule.** In jet-efficiency terms
   it gets 0.4-0.6, against the 0.8-1 the parent assumes for a magnetic nozzle. With cheap reaction
   mass that is enough: the projectile bill, not Isp, is the figure of merit.

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
is only compressed by the ram pressure, which this model neglects. At 10 kg/m^3 that is ~50 Mbar
at 618 km/s (the Hugoniot; `rho0 v^2` alone is 38 Mbar), falling to ~5 Mbar at 206 km/s. One shock
at 50 Mbar puts ~0.8 GJ/kg into graphite, under 1% of Q. At that pressure the disk flows like a
liquid, though, and rim relief runs 10-20 cm inward during the sweep (see Limits). Rows
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
  with 60% of the heat (73% of the light) above 1 keV.
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

## Result 6: argon gas bags at 66 km/s, 1 moving : 2 at rest -- `argon_bags_k2.csv`

`make xray-collision-argon-bags` (`argon_bags.py`, `multigroup.py`). A 10 kg argon impactor at
66 km/s meets 20 kg of argon at rest in the ship's frame (k = 2). The centre-of-mass motion takes
a third of the 21.8 GJ, and the heat is Q = 14.5 GJ. Both bodies are taken to be at the same
uniform density at the moment of impact. The question is whether more than half of Q leaves as
light, and whether the light is trapped.

**Why a multigroup check.** The gray model cannot see where the light sits relative to argon's
absorption. Cold argon is opaque above its 15.76 eV ionization edge (~1e5 m^2/kg) and clear below
its 11.6 eV resonance lines (< 0.01 m^2/kg). That is the right range for a ~4 eV fireball.
`multigroup.py` re-runs the same one-zone fireball from a new TOPS pull of 120 log groups over
1-300 eV (`data/tables/tops/tops_argon_groups.html`). Each group uses the gray model's thin/thick
bridge. The Planck and Rosseland group means in emission bracket line saturation inside a group.
The cold-target column is a lower bound. It removes everything the whole unshocked target
`rho0 R_target` would absorb at 0.5 eV, and counts it as lost.

| density at impact [kg/m^3] | impactor cloud radius | gray | **multigroup** | after cold target | < 11.6 eV | > 30 eV | half the light out by |
|---|---|---|---|---|---|---|---|
| 16 (intact 10 atm bag) | 0.53 m | 0.09 | -- | -- | -- | -- | -- |
| 1 | 1.3 m | 0.20 | -- | -- | -- | -- | -- |
| 0.1 | 2.9 m | 0.35 | -- | -- | -- | -- | -- |
| 1e-2 | 6.2 m | 0.57 | **0.50-0.52** | 0.40-0.42 | 0.35-0.37 | 0.014 | 0.42 ms |
| 3e-3 | 9.3 m | 0.72 | **0.66-0.68** | 0.50-0.54 | 0.45-0.46 | 0.031 | 0.39 ms |
| 1e-3 | 13 m | 0.88 | **0.83-0.85** | 0.65-0.69 | 0.55-0.57 | 0.042 | 0.25 ms |
| 3e-4 | 20 m | 0.99 | **0.95-0.96** | 0.79-0.81 | 0.61-0.63 | 0.034 | 55 us |
| 1e-4 | 29 m | 0.99 | **0.95-0.97** | 0.71-0.76 | 0.51-0.57 | 0.03-0.05 | 14-20 us |
| 1e-5 | 62 m | 0.99 | **0.91-0.96** | 0.28-0.42 | 0.21-0.29 | 0.04-0.08 | 3-20 us |

Shares are of Q. Denser than 1e-2 kg/m^3 the shocked gas lies above the pull's 0.1 kg/m^3 top,
so those rows are gray only. At 1e-2 the first shocked state (0.12 kg/m^3) sits just above it.

- **More than half of Q radiates once the bags are below about 1e-2 kg/m^3**, or about 10 g/m^3.
  80% needs about 1e-3 kg/m^3, which is one gram per cubic metre. At 1e-3 that is 12 GJ of light,
  56% of the impactor's ship-frame kinetic energy (the ceiling is 2/3). An intact bag radiates 9%.
- **The light is not trapped, because it sits in argon's windows.** The multigroup answer is
  within 2-8 points of the gray one. More than half of Q leaves below 11.6 eV, where cold argon is
  clear. Even if the whole cold target absorbs everything else and keeps it, 65-69% leaves at
  1e-3. Ionizing all 20 kg of target once costs 5.2% of Q, so the hard part should burn through.
- **It is vacuum ultraviolet, not EUV.** The fireball starts at 6-7 eV and emits at ~3-4 eV. At
  1e-3 the light is 40% at 6.2-11.6 eV, 19% at 11.6-15.8 eV, and only 4% above 30 eV.
- **Too thin is worse again.** Below ~1e-4 kg/m^3 the fireball is so thin that it radiates in
  lines above the 15.76 eV edge, which the cold target then absorbs (0.28-0.42 at 1e-5). The best
  densities are 3e-4 to 1e-3 kg/m^3.
- **Bag mass hardly matters at a fixed density.** At 1e-3 the share runs from 0.88-0.89 at
  0.1 kg to 0.76-0.78 at 1000 kg. Density is the lever, not size.

**A lattice of small bags must be burst, not left intact (hand estimates).** Argon at 10 atm and
300 K is 16 kg/m^3, so a 2.5 cm bag holds 0.13 g. At a mean 1e-3 kg/m^3 the bags sit 0.5 m apart
and fill 6e-5 of the volume. A path across the target meets a bag with a probability of the
target column over one bag's areal density, `rho_mean L / (4/3 rho_bag r)` = 0.034/0.27 = 0.13.
So 87% of the impactor's gas passes through the gaps. The 13% that hit a bag meet 16 kg/m^3 gas
and radiate ~9%. The total is about 1% of Q. A lattice of intact elements acts as a continuous
medium only if each element's areal density is well below the column: sub-millimetre bags at
10 atm, or a mist. Bursting the bags works instead. A cubic lattice of Gaussian clouds of width
sigma at spacing d has a density ripple of ±25% at sigma/d = 0.4, ±11% at 0.45 and ±4% at 0.5
(`argon_bags.lattice_ripple`). Room-temperature argon spreads at about 250 m/s per axis, so a
0.5 m lattice needs about 0.9 ms of lead. One 20 kg bag needs ~25 ms and arrives densest at its
core. Two costs come with it. A bag skin is condensed matter that radiates almost nothing here; a
10 um film is ~20% of the argon's mass on a 2.5 cm bag at 10 atm. The expanding argon will also
probably condense into ~50 nm clusters (Hagena's jet-scaling, applied loosely), but at 66 km/s
those are shredded and re-merged within a nanosecond.

**Not usable in the pass-through pusher.** The cloud that radiates is metres to tens of metres
across. It cannot thread the ship's bore, and the radiation case would have to enclose it. A
compact projectile that fits the bore needs hundreds of metres of dilute target (`piston`). See
[`xray_pusher_66kms_handoff.md`](xray_pusher_66kms_handoff.md).

Not modelled beyond the general limits below: the clouds are taken to be uniform; the collision
is taken to be collisional (ion-ion mean free path ~0.2 m at 1e-3 kg/m^3 against a 13 m cloud,
marginal by 1e-5); and LTE is assumed, though electron excitation at 1e-3 kg/m^3 would radiate
the heat in ~0.1 us (a rough coronal estimate, loss rate ~3e-31 W m^3) against a ~1 ms expansion, so emission is not the limit and escape is.

## Result 7: 400 km/s, a carbon or a steel lattice -- `near_sun_400.csv`

The same dense 1 kg graphite projectile and `piston` model as Result 5, at a 400 km/s closing
speed: a 200 km/s periapsis, met head-on by a ship moving at about the same speed. The target is
a carbon lattice or a steel one. Steel is modelled as pure iron (NIST ladder, Saha), with a new
TOPS pull to 3 keV (`data/tables/tops/tops_iron_hot_gray.html`). The older iron pull stops at
60 eV, and these fireballs start at 130-340 eV. The question was whether the heavier atom, with
more energy per atom, gets hotter and radiates more.

Per kilogram of projectile, ship frame, ideal Saha. "Light/KE" is the radiated share of the
projectile's 80 GJ of kinetic energy. The iron range spans opacity 0.1x to 10x.

| k | lattice rho0 | C kT0 | Fe kT0 | **C light/KE** | **Fe light/KE** (0.1x-10x) | C / Fe heat above 1 keV |
|---|---|---|---|---|---|---|
| 1 | 1 | 376 eV | 233 eV | **0.42** | 0.39 (0.44-0.30) | 21% / 13% |
| 2 | 1 | 218 | 189 | **0.49** | 0.46 (0.54-0.32) | 5% / 6% |
| 3 | 1 | 140 | 166 | 0.45 | 0.46 (0.57-0.31) | 0.6% / 3% |
| 2 | 3 | 218 | 207 | 0.35 | 0.32 | 5% / 5% |
| 2 | 10 | 219 | 229 | 0.18 | 0.17 | 3% / 1% |
| 3 | 10 | 145 | 197 | 0.10 | 0.15 | 0.1% / 0.4% |

Continuum lowering moves iron by under 0.005. A 1 kg/m^3 target of 2 kg (k = 2) is a column 290 m
long behind the 4.7 cm graphite ball, or 16 m behind the handoff's 40 cm disk.

- **About half the projectile's kinetic energy can still leave as light, but only from a lattice
  ten times lighter.** The best case is carbon at k = 2 and 1 kg/m^3: 0.49, against 0.47 for the
  618 km/s recommendation at 10 kg/m^3. At 10 kg/m^3 and 400 km/s the share falls to 0.18.
- **The light is softer.** At k = 2 it is emitted near 110 eV, against 320-380 eV at 618 km/s.
  Above 1 keV is 3% of the kinetic energy, against 28%.
- **Steel does not run hotter.** An iron atom gets 4.6 times carbon's energy, but at about
  Ne-like stripping it shares it with ~20 free electrons and pays several keV to free them. At
  k = 1 iron starts at 233 eV against carbon's 376 eV. Iron only matches or beats carbon at k >= 3
  or in denser lattices, where carbon has fallen below its K shell and iron still emits from its
  L shell.
- **Iron's difference is opacity.** It starts 30-1000 optical depths thick where carbon starts
  near one. Its answer therefore depends far more on the opacity (0.54-0.32 at k = 2).
- **Carbon is the choice.** It is level with or better than iron at k <= 2, its result is less
  sensitive to opacity, and it is the material already shown to survive the near-Sun exposure.

The limits of Result 5 carry over. At 1 kg/m^3 the sweep is longer still: ~80 us behind the
40 cm disk and ~1.5 ms behind the ball. So the time-graded heating described under Limits
matters more here, not less.

## Result 8: a light foam stopped on a dense anvil, and its thrust ledger -- `anvil_sweep.csv`, `anvil_ledger.csv`

`make xray-collision-anvil` (`anvil.py`). This reverses Result 5's density contrast. A **light**
projectile, 10 kg of carbon foam at 10 kg/m^3, arrives at 400 km/s and stops on a **dense**
ship-held **anvil**. The anvil sits at the far end of a radiation case, which sends the light into
the parent's channel filler and propellant plate. The settled design:

- **Projectile.** A Parker Solar Probe-style sandwich: a foam cylinder 1 m across and 1.3 m long,
  so it threads the bore, with C-C sheets of 0.5 mm or less on the **sunward side only** (~1.8 kg)
  and **no front cap**. The 10 kg includes the sheets.
- **Anvil.** 50 kg (k = 5) of graphite, 2.9 cm thick, with a thin high-Z liner (tungsten or lead)
  for albedo. It is expended and leaves astern at 67 km/s.
- **Case and plate.** The case re-emits into the channel filler. The filler blows off astern and
  drives the propellant plate onto the pusher plate at no more than 150 km/s.

**The collision is Result 5 in another frame.** A dense anvil swept by light foam is the same
collision as a dense projectile sweeping a light target, so `fireball.piston` is reused with the
roles swapped. The whole heat `Q = (1/2) mu v^2` lands in the light partner. Here that is the
foam, which stagnates on the anvil's face at ~700 eV. Shares below are of the projectile's
800 GJ (80 GJ/kg).

| anvil | heat share | kT0 | light, foam at 10 kg/m^3 (opacity 0.1-10x) | light, 1 kg/m^3 | above 1 keV (10 kg/m^3) | RT burn-through, graphite / iron |
|---|---|---|---|---|---|---|
| 20 kg | 0.67 | 534 eV | 0.42 | 0.57 | 0.22 | 1.6 / 0.8 us |
| 40 kg | 0.80 | 661 eV | 0.55 (0.62-0.48) | 0.71 | 0.34 | 3.1 / 1.7 us |
| **50 kg** | **0.83** | **692 eV** | **0.59 (0.65-0.51)** | **0.74** | **0.37** | **3.9 / 2.1 us** |
| 100 kg | 0.91 | 764 eV | 0.66 (0.73-0.59) | 0.82 | 0.44 | 7.8 / 4.1 us |

- **Dense mass in the projectile is dead weight.** A C-C sheet that strikes the anvil is a
  dense-on-dense collision, and it stays opaque (Result 4: solid graphite radiates 0.4%). With
  1.8 kg of sheets the design's light share falls from 0.59 to **0.50**. A front cap would also
  hit the liner with a ~1 Gbar spike (about `rho u^2`) before the foam arrives, so it is omitted.
- **The anvil must outlast the pulse.** The foam pushes on the anvil at about `1.5 rho0 v^2`
  (24 Mbar) for its whole arrival, 3.2 us (2.6 us for the 8.2 kg of foam inside the sheets). A light
  plasma pushing a heavy plate is Rayleigh-Taylor unstable. Bubbles eat in as `h ~ 0.05 g t^2`,
  so at a fixed areal mass the burn-through time goes as `sqrt(thickness)`. The anvil should
  therefore be **thick and low-density**. A 50 kg graphite anvil holds for 3.9 us, but the same
  mass of iron fails at 2.1 us. This classical mixing-layer rule is good to about a factor of 2.
- **The anvil's material does not change the foam's light.** The heat is in the carbon foam. The
  anvil sets the case's albedo through its liner, which favours high Z: tungsten, lead or gold
  rather than iron. Spaced "blanket" foils do not work. At a ~250 eV drive each foil turns to
  plasma within nanoseconds, and a case re-emits from a thick high-Z wall.

**The case.** In a steady cavity the filler's share is `A_f / (A_f + (1 - alpha) A_wall + A_hole)`
(`case_delivery`). Two geometric facts bear on it. First, the drive lasts microseconds, not
nanoseconds, and a wall's loss fraction falls with drive time at a fixed temperature, so a high
albedo is plausible. It is not computed here. Second, the 1 m projectile needs a 1 m entry hole
through the filler. The arriving foam's own cold tail plugs that hole while it is still arriving
(cold carbon is ~200 m^2/kg at 1 keV against its 13 kg/m^2 column). Light emitted after the tail
has passed through the hole leaves toward the ship. Both are open. The ledger takes the parent's 80% delivery as its
reference and also shows 50% and 65%.

### Thrust and energy efficiency (per kg of projectile)

The parent's plate stage (§`pass_through_pusher`, eq. `pulse_split`) uses the same optimizer as
the parent's table. `best_split` reproduces `tab:pass_through_isp` to within 0.2% (tested). The fixed
spent mass is the anvil (5 kg/kg) plus case walls and liner. The ship's impulse is `J = c p` in
terms of the propellant plate's momentum `p`. The parent takes `c = 0.8`, with no bounce credit.
If the propellant rebounds from the plate, `c = 2f`, where `f` is this repository's fudge factor:
`c = 1.6` at `f = 0.8`. That `f` was measured at 3-16 km/s and has not been tested at the
~100 km/s used here.

Case delivery 80%, impingement cap 150 km/s, walls and liner 1.5 kg/kg (which fills the parent's
6.5 kg/kg case budget):

| | no bounce (c = 0.8) | c = 1.2 | rebound, f = 0.8 (c = 1.6) |
|---|---|---|---|
| light | 40.2 GJ (50% of KE) | same | same |
| into filler and propellant motion | 32.2 GJ (40%) | same | same |
| propellant plate toward the pusher | 3.3 kg at 99 km/s, 16.1 GJ (20%) | same | same |
| filler blow-off astern | 3.3 kg | same | same |
| impulse on the ship | 259 kN s | 389 kN s | 518 kN s |
| impulse / projectile momentum (400 kN s) | 0.65 | 0.97 | 1.30 |
| spent mass | 13.0 kg | 13.0 kg | 13.0 kg |
| **specific impulse** | **2030 s** | **3050 s** | **4060 s** |
| jet efficiency `(1/2) M_spent v_e^2 / KE` | 3.2% | 7.3% | 13% |

With 65% and 50% case delivery, the no-bounce Isp is 1830 s and 1610 s. The cap never binds, because
the propellant arrives at 79-99 km/s.

- **Against the parent's bar of 2650 s, this design clears only with a rebound.** Without one it
  gives 2030 s. With the parent's 618 km/s design at a 5 kg case (2930 s, 14 kg spent) it shares
  the same impulse-to-projectile-momentum ratio, 0.65. It needs no lattice target and runs at
  400 km/s instead of 618.
- **Energy efficiency is low by construction.** The projectile arrives at 400 km/s. The spent mass
  leaves at an effective 20 km/s. The 67 km/s anvil and the blow-off astern carry momentum the
  ship does not want. About 3% of the projectile's kinetic energy ends up as useful
  exhaust-equivalent energy (13% with a full rebound). The parent's 618 km/s design gets 3.0% by
  the same measure.
- **A lighter anvil buys Isp.** With walls held at 1.5 kg/kg and no bounce, Isp is 2400 s at
  20 kg, 2150 s at 40 kg, 2030 s at 50 kg and 1610 s at 100 kg. The extra light from a heavier
  anvil does not pay for its spent mass. Rayleigh-Taylor survival sets the floor: a 40 kg graphite
  anvil (3.1 us) still covers the 2.6 us pulse of the sheeted foam, and 30 kg (2.4 us) does not.
- **Which side should be heavy?** Eq. `pulse_split` is symmetric in the two masses. At a fixed
  total, equal masses give the most momentum. With 6.5 kg/kg already spent on the anvil and case,
  the optimum is **equal filler and propellant**. A light forward layer and a heavy backward one
  (1.6 + 4.9 kg) lose 13% of the impulse. They also send the propellant at 172 km/s, over the
  cap. Adding heavy backward mass instead (3.3 + 9.8 kg) raises the impulse 23% but lowers Isp
  18%. Materials matter for other reasons. The backward layer is the filler's blow-off, and it
  must be low Z to absorb the light, because high Z re-emits it (that is the case wall's job). The
  forward layer is free in mass. A hydrogen-rich propellant such as polyethylene radiates least
  where it stagnates on the plate, and that is what any rebound depends on.

Open: the case-wall mass needed to hold a ~5 Mbar fireball for microseconds (a rough
inertial-tamping estimate is ~20 kg/m^2, or ~100 kg, against the 15 kg the budget leaves); the
case albedo and the entry hole; the rebound at ~100 km/s; sheet survival through a ~10 solar
radii dive; and the time-graded stagnation, which the one-zone model replaces with an average.

## Result 9: the propellant's rebound and the entry hole -- `propellant_rebound.csv`, `anvil_hole.csv`, `anvil_design.csv`

Result 8 left two inputs as assumptions. One is how much of the forward momentum the pusher plate
catches. The parent takes 0.8 with no bounce credit. The other is how much of the light the case
delivers to the filler. The parent assumes 80%. This result computes both for the Result 8
design: a low-Z filler and a polyethylene propellant plate of equal mass, the iron backing dropped.

### Rebound: the f(v) bounce at 60-150 km/s (`make xray-collision-rebound`)

This is the `f(v)` study's coupled bounce (`hydro1d` `CoupledBounce`: a 1-D Lagrangian slug on a
rigid wall, with flux-limited radiation), run with a new polyethylene table. `data/tables/ch2.json` (`make
xray-collision-ch2-table`) is an atomic C/H Saha EOS (`spray_eos.AtomMixture`, carbon ladder
added) with TOPS gray means for C:H = 1:2 (`tops_ch2_gray.html`, 0.5-300 eV). A cold uniform slug of
density `rho` and length `L` strikes the plate at `v`, and `e_eff = J_wall / p_in - 1`. At
100 km/s the stagnated gas runs at about 38 eV.

| `rho`, `L` (column) | 60 km/s | 80 | **100** | 120 | 150 |
|---|---|---|---|---|---|
| 0.1, 1 m (0.1 kg/m^2) | 0.66 | 0.55 | 0.44 | 0.36 | 0.29 |
| 0.3, 1 m (0.3) | 0.71 | 0.64 | 0.51 | 0.42 | 0.34 |
| 1, 1 m (1) | 0.73 | 0.72 | **0.63** | 0.53 | 0.42 |
| 1, 3 m (3) | 0.73 | 0.73 | **0.68** | 0.59 | 0.49 |
| 3, 1 m (3) | 0.74 | 0.74 | **0.72** | 0.66 | 0.61 |
| 10, 1 m (10) | 0.75 | 0.75 | 0.75 | 0.75 | 0.74 |

At 100 km/s, opacity scaled from 0.1x to 10x gives 0.54-0.73 at `rho = 1` and 0.65-0.74 at
`rho = 3`. Two rows (`rho = 10`, `L = 3` m at 120 and 150 km/s) did not converge. They are flagged
in the CSV and left out.

- **The propellant does rebound.** At the ~95 km/s the optimizer chooses, `e_eff` is 0.63-0.72
  for an arriving column of 1-3 kg/m^2. With the parent's 0.8 geometric capture, the plate catches
  `0.8 (1 + e)` = **1.30-1.38** times the forward momentum, against the parent's 0.8. The ceiling
  is about 0.75, which the densest columns reach at every speed.
- **Thin and fast is lossy.** Below ~1 kg/m^3, the stagnated layer is thin enough for its light
  to escape upstream before it re-expands. The rebound falls to 0.44-0.51 at 100 km/s and to
  ~0.3 at 150 km/s. That is a reason to stay well under the impingement cap, not just at it. The
  column is the propellant's mass over the plate's area. 30-40 kg per pulse needs a plate under
  about 30 m^2 to stay at or above 1 kg/m^2.

### The entry hole (`hole_delivery`)

The projectile enters the case through a hole its own diameter wide, and the hole is lost area. Two
facts set how much it costs:

- **Emission tracks arrival.** A stagnated parcel sends out half its light within 3.5 ns and 90%
  within 0.3 us (the one-zone history, `emission_kernel`), while the pulse lasts microseconds.
- **The foam's own tail plugs the hole**, but only while it is still in the hole. The plug helps
  only if the foam is **longer than the case**. A 1 m wide foam is 1.0 m long, so in a 1.5 m case
  its tail has passed the hole before the head stagnates, and all of the light meets an open
  hole.

The model is a time-dependent cavity balance:

    A_filler / (A_filler + (1 - alpha) A_wall + (1 - alpha_c) A_column + A_open)

The foam in flight absorbs light and recycles it at `alpha_c = 0.7` (about its radiated share
once stagnated), bracketed at 0.5-0.9. The wall albedo is 0.95. The table gives the share of the
light delivered to the filler:

| hole = foam diameter | foam length / pulse | case 1.6 x 0.5 m | case 1.6 x 1.5 m | case 2.5 x 0.5 m |
|---|---|---|---|---|
| 1.0 m | 1.0 m / 2.6 us | 0.53 | 0.42 | 0.76 |
| 0.7 m | 2.1 m / 5.3 us | 0.70 | 0.55 | 0.83 |
| **0.5 m** | 4.2 m / 10 us | **0.78** | 0.62 | 0.87 |
| 0.35 m | 8.5 m / 21 us | 0.82 | 0.67 | 0.88 |

**A smaller hole does help.** In a short case it takes delivery from 0.53 to 0.78, which is the
parent's 80%. Most of the gain comes from the plug, which covers all but 16% of the light at
0.5 m. Below ~0.35 m the returns are small. A wider, shorter case helps as much as a smaller hole.

Three costs come with a narrow foam:
- **Its sheets weigh more.** The same volume has `4V/d` of side area, so the sunward C-C is 2.6 kg
  at 0.5 m against 1.5 kg at 1 m, and the light share falls from 0.52 to 0.46.
- **The pulse is longer:** 9.4 us. The 50 kg anvil is then 0.5 m wide and four times as thick per
  area, so it still outlasts the pulse (Rayleigh-Taylor 15.7 us).
- **The stagnated plasma spreads further.** At ~250 km/s it moves ~2.5 m sideways during the
  pulse, so it reaches the case walls early. That makes the case-tamping problem worse.

### The design ledger (per kg of projectile, 400 km/s)

The 50 kg anvil, walls and liner at 1.5 kg/kg (the parent's 6.5 kg/kg cap), a case 1.6 m across
and 0.5 m long, the computed delivery, and a 150 km/s cap that does not bind:

| | light | delivered | rebound | **Isp** | impulse | propellant speed | jet efficiency |
|---|---|---|---|---|---|---|---|
| 1 m hole, no bounce (parent's convention) | 0.52 | 0.54 | none | 1680 s | 215 kN s | 82 km/s | 2.2% |
| 1 m hole, 1 kg/m^2 column | 0.52 | 0.54 | 0.63 | **2740 s** | 350 kN s | 82 km/s | 5.9% |
| 1 m hole, 3 kg/m^2 column | 0.52 | 0.54 | 0.72 | 2890 s | 369 kN s | 82 km/s | 6.5% |
| 0.5 m hole, no bounce | 0.46 | 0.78 | none | 1920 s | 244 kN s | 94 km/s | 2.9% |
| **0.5 m hole, 1 kg/m^2 column** | 0.46 | 0.78 | 0.63 | **3130 s** | 398 kN s | 94 km/s | 7.6% |
| 0.5 m hole, 3 kg/m^2 column | 0.46 | 0.78 | 0.72 | 3300 s | 420 kN s | 94 km/s | 8.5% |

- **With the computed rebound, the design clears the parent's 2650 s bar.** It gives 2740-2890 s
  with the 1 m hole and 3130-3300 s with a 0.5 m one. Without a rebound it fails at 1680-1920 s.
  The rebound is worth a factor of 1.6 in impulse. It is the largest single lever in the ledger,
  and the parent does not count it.
- **The 0.5 m hole is worth about 400 s.** The higher delivery more than pays for the heavier
  sheets.
- **Energy efficiency rises to 6-8.5%**, from 2-3% without a bounce. It stays low for the reason
  given in Result 8: the exhaust leaves at an effective 27-32 km/s, while the projectile arrived
  at 400 km/s.

Open, beyond Result 8's list:
- The plate is rigid and does not ablate. At ~95 km/s the `v/4` rule gives it about 24 kJ per
  N s to absorb. The f(v) study's ablating-wall arm (`--ablating`) is the tool for that.
- The arriving propellant is a uniform cold slug. A real jet is stretched, with a velocity
  gradient, and the drive shock has preheated it.
- The column the plate sees depends on its area, which is the parent's to fix.
- The cavity balance is steady at each instant, with no shadowing by the in-flight foam column.
  Its column albedo is a bracketed estimate, not a calculation.

## Result 10: a cheap polyethylene charge and a point fireball onto the plate -- `charge_ledger.csv`

`make xray-collision-charge-2d`, `make xray-collision-charge` (`charge.py`,
`euler2d::bounce::run_fireball_bounce`). This design uses no radiation and no case. It suits an
architecture where projectiles cost money and ship-carried mass is cheap. A 1 kg graphite disk at
400 km/s threads the bore. It stops in a **charge** astern of the plate: near-solid polyethylene,
k = 20-40 times its mass, with a graphite skin, kept in the ship's shadow until the pulse. The
merged fireball expands onto the parent's 20 m plate. The mission falls from Jupiter to a 9.86
solar radii periapsis (Parker Solar Probe's) and boosts there to cross Earth's orbit at 100 km/s.
That takes **20.8 km/s**, so a mass ratio of 3 needs Isp 1930 s. The fall from Jupiter arrives
3.5 km/s faster than a fall from 1 AU.

**The fireball is dense, so it does not radiate.** On the polyethylene table, a k = 30 charge
starts at 32 eV and about 2e5 optical depths thick, and radiates under 0.01% of its heat. Shocked
solid polyethylene, ~4x 950 kg/m^3, lies above the table's 1284 kg/m^3 ceiling. The one-zone
start is therefore held at the ceiling, which only understates the opacity. Its heat
all becomes expansion: `u_rms` = 85, 71 and 62 km/s at k = 20, 30 and 40. It drifts away from the
plate at `w = v/(1+k)` = 19, 13 and 10 km/s. That holds only if the charge is **compact and in
the disk's path**. The disk couples only to the column it crosses: at k = 20 it needs ~160 kg/m^2
of polyethylene, 17 cm of solid against 5 m of 30 kg/m^3 foam. A long foam column would be a line
source exploding sideways. A loose cloud of foam balls would mostly be missed. The snowplow also
puts half its heat into the first projectile-mass of polyethylene it sweeps. Inside a solid charge,
shocks spread that heat in ~4 us, long before the gas reaches the plate.

**The plate, in 2-D.** The merged mass starts as a uniform hot sphere, 0.5 m in radius, drifting
away. The run is lossless effective-gamma Euler (gamma 1.4; the CH2 table gives 1.3-1.5 over these
states) with an immersed paraboloid and cups. As in the f(v) study, the gas-dynamic rebound
`1 + e_gamma` = 1.749 (the confined plane wave) is swapped for polyethylene's radiative
`1 + e_eff` from Result 9, at the plate speed and arriving column. The impulse inside the 0.25 m
bore is subtracted. It is 0.1-0.5% everywhere. "Capture" is the plate impulse per `M u_rms`. The
energy ceiling is `1 - w/u_rms` = 0.82 at k = 30.

| plate (20 m across), k = 30 | capture | Isp | mass ratio | projectile kg per kg final | peak pressure, 1 kg projectile | wetted area |
|---|---|---|---|---|---|---|
| paraboloid δ/D 0.25, charge at focus (5.0 m) | 0.37 | 2530 s | 2.31 | 0.044 | 95 MPa | 383 m^2 |
| **parent's bowl δ/D 0.30, at focus (4.2 m)** | **0.44** | **2990 s** | **2.03** | **0.034** | 163 MPa | 409 m^2 |
| paraboloid δ/D 0.50, at focus (2.5 m) | 0.58 | 3980 s | 1.70 | 0.023 | 783 MPa | 533 m^2 |
| paraboloid δ/D 0.75, at focus (1.7 m) | 0.64 | 4380 s | 1.62 | 0.021 | 2580 MPa | 713 m^2 |
| bowl 0.30 + 2 m skirt (the parent's cup) | 0.50 | 3400 s | 1.86 | 0.029 | 163 MPa | 535 m^2 |
| bowl 0.30 + 4 m skirt | 0.52 | 3540 s | 1.82 | 0.027 | 162 MPa | 660 m^2 |
| bowl 0.30 + 8 m skirt | 0.53 | 3580 s | 1.81 | 0.027 | 162 MPa | 912 m^2 |
| flat + 4 m skirt, charge 2.7 m off | 0.43 | 2940 s | 2.06 | 0.035 | 465 MPa | 565 m^2 |
| bowl 0.30, charge at 1.5 x focus | 0.37 | 2540 s | 2.30 | 0.043 | 55 MPa | 409 m^2 |
| **paraboloid δ/D 0.50, charge at 1.5 x focus (3.75 m)** | **0.57** | **3860 s** | **1.73** | **0.024** | **280 MPa** | 533 m^2 |

Rebound `e_eff` = 0.60 at a 1 kg projectile. Its fireball lays only ~0.05 kg/m^2 on the bowl,
below the rebound table's 0.1 kg/m^2 edge where it is held, so `e` may be a few hundredths lower.
A 10 kg projectile (0.5 kg/m^2, `e` 0.69) gains 5-6% in Isp.

| leading designs | k = 20 | k = 30 | k = 40 |
|---|---|---|---|
| parent bowl: Isp / mass ratio / projectile per final / peak | 3200 s / 1.94 / 0.047 / 131 MPa | 2990 s / 2.03 / 0.034 / 163 MPa | 2800 s / 2.13 / 0.028 / 184 MPa |
| bowl + 4 m skirt | 3830 s / 1.74 / 0.037 / 132 MPa | 3540 s / 1.82 / 0.027 / 162 MPa | 3290 s / 1.90 / 0.023 / 184 MPa |
| paraboloid 0.50 at 1.5 F | 4210 s / 1.65 / 0.033 / 214 MPa | 3860 s / 1.73 / 0.024 / 280 MPa | 3580 s / 1.81 / 0.020 / 320 MPa |

On a grid twice as fine (k = 30), capture moves by under 1.2% (0.438 to 0.442, 0.519 to 0.525,
0.566 to 0.570). Peak pressure is not converged. It falls 18-30% (163 to 133 MPa, 280 to 195 MPa),
so the production values above are upper ones.

- **The design works, with a wide margin.** Every plate reaches the mass ratio of 3. The
  parent's own bowl, unchanged, gives 2.0. Its 2990 s is 1.4 times the ballistic estimate of
  2080 s for an infinite flat plate at the same k, because the gas is pressure-borne and the bowl
  turns it toward the mouth. The plug study found the same.
- **A deep paraboloid beats a bowl with a skirt, and its advantage is enclosure, not focus.**
  Moving the charge from the δ/D 0.5 focus to 1.5 F costs 3% of capture but cuts the peak
  pressure 2.8 times. A skirt on the 0.30 bowl saturates near 4 m (0.52). The 0.5 paraboloid
  gets 0.57 with 19% less wetted area than bowl plus 4 m skirt. Per square metre, and so per
  kilogram of plate at the parent's survival floor of tens of kg/m^2, it is the best of the set.
  A flat plate with a skirt does no better than the parent's bowl and peaks three times higher.
- **The mass-ratio limit does not bind, so the projectile bill decides k.** Raising k from 20 to
  40 costs 10-15% of Isp and cuts projectiles per kilogram of final mass by 40%. The optimum lies
  above k = 40.
- **Peak pressure caps the projectile at about 1-2.5 kg per pulse on a 20 m plate.** Pressure
  scales with the fireball's mass at fixed geometry. Against the f(v) study's 400 MPa facesheet
  limit, the parent bowl takes about 2.5 kg at k = 30 and the 0.5 paraboloid at 1.5 F about 1.4 kg
  (2 kg on the fine grid). Per pulse that is about 0.9 MN s from a 1 kg projectile, so a 100 t
  final mass at a mass ratio of 1.8 needs ~2700 pulses, about one every 13 s across a ~10 h
  periapsis arc.

### Higher k, the bond toll and the finite arc

The rebound table now runs down to 20 km/s (`--propellant-rebound low`). Polyethylene rebounds
*better* when slow: `e_eff` = 0.68-0.73 at 30-50 km/s and 0.1-1 kg/m^2, and 0.73-0.79 at 20 km/s.
The leading plates were re-run at k = 60, 80 and 100, where `u_rms` is 51, 44 and 40 km/s and the
drift is 6.6, 4.9 and 4.0 km/s. The table gives Isp and its value if no polyethylene bond re-forms
(84 MJ/kg, which the CH2 table, counting atoms only, does not charge). It also gives the mass ratio
for an impulsive boost and for a 10 h burn arc centred on periapsis (`finite_arc_mass_ratio`, a
heliocentric integration). Values are for a 1 kg projectile on the production grid.

| k | parent bowl: Isp (bonds frozen) / mass ratio (10 h arc) / projectile per kg final / peak | δ/D 0.5 paraboloid at 1.5 F: same |
|---|---|---|
| 30 | 2990 (2940) s / 2.03 (2.05) / 0.034 / 163 MPa | 3860 (3800) s / 1.73 (1.74) / 0.024 / 280 MPa |
| 60 | 2450 (2370) s / 2.37 (2.40) / 0.023 / 214 MPa | 3100 (3000) s / 1.98 (2.00) / 0.016 / 349 MPa |
| 80 | 2200 (2100) s / 2.62 (2.65) / 0.020 / 231 MPa | 2760 (2640) s / 2.15 (2.17) / 0.014 / 395 MPa |
| 100 | 2015 (1905) s / 2.86 (2.90) / 0.019 / 239 MPa | 2520 (2380) s / 2.32 (2.34) / 0.013 / 405 MPa |

(The bowl with a 4 m skirt sits between them: 2330 s, mass ratio 2.48 and 0.015 at k = 100.)

- **High k pays.** From k = 30 to 100 the projectile bill falls by 45%. At k = 100 the 0.5
  paraboloid needs 0.013 kg of projectile per kilogram of final mass. Even the parent's unchanged
  bowl meets the mass ratio of 3 (2.86). If no bond re-forms it is at the edge (Isp 1905 s, mass
  ratio 3.04).
- **The polyethylene bond toll is minor.** It costs 3% of exhaust velocity at k = 60 and 5.5% at
  k = 100, in the worst case where no bond re-forms. Ionization holds ~18% of the energy at the
  start. If it froze at once, it would cost ~10% more of the velocity. The fireball starts at
  solid density, where three-body recombination is very fast, so most of it should return.
  Freeze-out is not computed here.
- **The finite arc is minor.** A 10 h arc raises the mass ratio by about 1%. Against impulsive,
  20 h raises it 3-5% and 40 h 9-15%.
- **More projectile mass per pulse does not come from higher k.** At fixed plate geometry the peak
  pressure goes as the projectile's energy and rises slowly with k: 163 to 239 MPa per kilogram
  for the bowl. The levers are geometry:
  - Moving the charge off the focus: the bowl at 1.5 F peaks at 55 MPa/kg at k = 30, so it takes
    ~7 kg per pulse at 400 MPa, for capture 0.37 against 0.44.
  - A larger plate: peak pressure goes as `m/R^3` at a fixed shape, so a plate twice the diameter
    takes eight times the projectile.
  - A stronger facesheet.

  A graphite disk 1 mm thick that fits the 0.5 m bore weighs 0.43 kg (2.3 mm makes 1 kg). It
  needs only ~2.2k kg/m^2 of charge column to stop.

### Checks: freeze-out, mixing, the exit flow and the back-splash

**Freeze-out (`make xray-collision-freeze`, `freeze.py`).** The fireball's expansion was followed
on the CH2 table from solid density to the plate. Its ionization store was relaxed toward Saha
equilibrium at `recombination.py`'s three-body plus radiative rates, and H2 re-formation was raced
with `rates.H_ATOMIC_H`.

| k | ionization at start | still frozen at the plate | slowest race (Da) | gas at the plate | bond toll on velocity |
|---|---|---|---|---|---|
| 20 | 18.9% | 0.03% | 2e5 | 1.36 eV | 1.1% |
| 30 | 19.3% | 0.03% | 2e5 | 1.21 eV | 1.6% |
| 60 | 17.9% | 0.01% | 1e4 | 0.77 eV | 3.3% |
| 100 | 17.6% | ~0% | 1.5 | 0.43 eV | 5.4% |

- **Ionization is not a loss at any k.** The fireball starts at solid density, where ions recombine
  far faster than it expands. At k = 100 the race nears freezing (Da 1.5), but by then almost
  nothing is left to freeze.
- **The bonds are paid in full.** The gas reaches the plate at 0.4-1.4 eV, too hot for H2, and
  stagnation reheats it. No bond re-forms while the gas can still push, so the "bonds frozen"
  Isp in `charge_ledger.csv` is the one to quote.

**Mixing (`make xray-collision-charge-mix`).** The point fireball assumed that the disk's energy
is shared through the charge. The test runs the real collision in 2-D. A 0.5 m disk moving away
from the plate strikes a charge of its own footprint: near-solid (0.16 m at k = 30, 0.54 m at
k = 100, on a grid twice as fine) or 30 kg/m^3 foam (5.1 m). Capture is per `M u_rms` throughout.

| plate | k | point fireball | disk into solid charge (2x coarser grid) | disk into foam charge |
|---|---|---|---|---|
| parent bowl | 30 | 0.439 | **0.473** (0.422) | 0.267 |
| δ/D 0.5 at 1.5 F | 30 | 0.567 | **0.580** (0.553) | 0.481 |
| parent bowl | 100 | 0.512 | **0.642** (0.588) | -- |
| δ/D 0.5 at 1.5 F | 100 | 0.641 | **0.761** (0.730) | -- |

- **A solid charge mixes well enough.** The real collision gives as much impulse as the uniform
  fireball or more: -4% to +25%, and more at k = 100 and on the finer grid. The snowplow puts its
  hottest, first-swept mass on the plate's side, and that mass leaves toward the plate. Result
  10's capture is conservative.
- **Foam does not mix well.** A 5 m column is a line source. It loses 15-39% of the impulse, and
  its exhaust is the least uniform of all (below).

**The exit flow: no nozzle, so uneven speeds.** The plate impulse comes from integrating the wall
pressure of the full 2-D flow, so gas leaving at uneven speeds and angles is already in every
number above. Separately, the gas's mass, momentum and speed moments are tallied as it crosses
the domain's outer faces. Nothing is lost (exit mass 1.00):
- the **direction factor** `P_z / Σ m|u|` is 1 if all gas leaves along the axis;
- the **speed uniformity** `(Σ m|u|)^2 / (M Σ m u^2)` is 1 if it all leaves at one speed;
- their product, `dir^2 x unif`, is what a nozzle designer would call jet efficiency.

| exhaust, k = 30 | direction | uniformity | dir^2 x unif |
|---|---|---|---|
| point fireball, bowl / paraboloid | 0.70 / 0.87 | 0.83 / 0.83 | 0.41 / 0.62 |
| disk into solid, bowl / paraboloid | 0.76 / 0.89 | 0.65 / 0.64 | 0.37 / 0.51 |
| disk into foam, bowl / paraboloid | 0.61 / 0.91 | 0.45 / 0.44 | 0.17 / 0.36 |

An ideal nozzle would have 1 and 1. The parent assumes 0.8-1 for its magnetic nozzle. A plate
recovers less directed energy per joule than that. It still meets the mission (Result 10) because
the charge, not the projectile, supplies the reaction mass. The deep paraboloid is the most
nozzle-like of the plates.

**The back-splash jet.** The real collision throws a narrow, early jet of the first-shocked
charge back along the axis. Its impulse is small (0.3-0.5% falls inside the 0.25 m bore, about
five times the point fireball's), but its pressure is not. Just outside the bore, at r 0.26-0.56 m
and within the first ~5% of the pulse, the plate's peak pressure per kilogram of projectile is
5-12 GPa at k = 30 and ~18-20 GPa at k = 100 on the production grid. On a grid twice as coarse
it is 2.1-7.1 GPa, so the jet is under-resolved and the true peak is at least this. The point
fireball's peaks are 137-301 MPa. **The ring around the bore, out to ~0.6-1 m, therefore needs a
sacrificial or ablative collar.** Possible alternatives are a charge shaped to tamp its own
splash, or a greater standoff. The rest of the plate sees the point fireball's loads.

Open:
- The collar. Its size, its ablation per pulse and its mass are not computed, and the jet's true
  peak needs a resolved, material-strength calculation.
- The energy the bore passes toward the ship (now ~5x the point fireball's share).
- The charge as a uniform equilibrated sphere. The snowplow's hot core and its early X-rays are
  not followed.
- The rigid, non-ablating plate. At ~70 km/s the `v/4` rule asks it to absorb ~18 kJ per N s.
- Peak-pressure convergence.
- The charge's thermal survival in the ship's shadow.

## Limits

**The piston's heat goes in over the sweep, not at once.** This is the largest unmodelled effect
in Result 5. The one-zone k = 2 fireball radiates 50% of Q by 10 ns, 61% by 0.1 us and 65% by
1 us. The dense disk takes ~5 us to sweep its 1.6 m column (`L (1 + k/2) / v`). In reality early
target mass is shocked at 618 km/s to ~2.2 keV and late mass at 206 km/s to ~160 eV. Each
parcel's thin cooling time at its Hugoniot state is 0.2-6 ns, so it radiates almost as it is
shocked, while the ram pressure holds it against the disk. Confinement removes the expansion
loss, so the radiated share probably rises. The spectrum becomes a sweep-weighted mix of
parcel temperatures rather than one cooling zone, and the k = 2 optimum may move. Two related
geometric effects are open. The unswept lattice ahead is optically thick to this light: cold
carbon absorbs ~200 m^2/kg at 1 keV against a column of up to 16 kg/m^2. Forward light
therefore preheats the lattice rather than escaping, so escape is sideways past the disk's
rim. Also, at 5-50 Mbar the 3.6 mm disk flows, and rim relief moves 10-20 cm inward (half its
radius) during the sweep, so it does not stay a flat rigid plate.

One zone, so there is no cool photosphere and no temperature gradient. For k > 1 the shock decays
into the target, and the one-zone mean heat stands in for a graded heating. A collision of two balls
is a pancake, not a sphere. LTE is assumed throughout, including the thin, cold tail. The
lowered EOS is a bracket, not a model.
