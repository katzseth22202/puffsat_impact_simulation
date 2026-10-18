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
     Earth's orbit at 100 km/s with a mass ratio of 2.2-2.3 at k = 20-30, priced on the real
     collision (corrected 2026-10-15). A δ/D 0.5 paraboloid gives 1.9-2.0. Above k ≈ 30 the mass
     ratio climbs steeply while the projectile bill barely falls, so k ≈ 30 is the practical
     choice.
   - **Cost:** projectiles come to ~3-4% of the final mass (0.033 kg/kg on the paraboloid at
     k ≈ 30).
   - **Checks:** ionization recombines and polyethylene's bonds cost 1.6-5.4% of velocity. The
     real disk-into-solid-charge collision gives 13-18% less impulse than the point fireball at
     k = 30, and 32-36% less at k = 100 (corrected 2026-10-15). On that basis mass ratio 2 needs
     k ≈ 19 / 26 / 28 on the bowl, the skirted bowl and the deep paraboloid. Foam is far worse.
   - **The one hazard** is a back-splash jet of several GPa onto the ring around the bore. Sized
     by impulse, its collar is ~0.5-1 t of extra structure on a 100 t ship, out to ~2-3.5 m. The
     plate also needs grading for the smooth fireball's centre load (+3.6-4.1 t in all). The
     collar's erosion, 0.1-9 t over the boost, leaves as reaction mass, and the mass ratio stays
     under 3 even at the ceiling.
5. **The same charge launches the ship from Earth** (Result 11, ADR-0061). A 100 km/s impactor
   overtaking a polyethylene charge held off the plate gives an 11 km/s departure at a mass ratio
   of 1.5 with k ≈ 10-12.
   - **Coupling:** with the merged mass drifting onto the plate it is better than at 400 km/s.
   - **Planar sheets:** a wide impactor sheet against a wide charge sheet changes the impulse
     little, but spreads the load. The plate then takes ~11 kg per pulse instead of ~2-6 kg.
   - **Cost:** the impactor is ~$70-100 per kg of propellant.
   - **The two stages close a loop** (ADR-0062). The Sun wave redirects Sun ships whose whole
     remaining mass becomes the Earth wave. Priced on the real collisions, the impactor grows
     ~9-12x a generation, and the expensive projectile costs ~$1.6-2.3 per kg of payload
     departure.
6. **The Earth wave can also drive the 3S Jovian growth cycle** (Result 12, ADR-0063). The stream
   meets the departure's thrust axis almost sideways (~101°).
   - **The plate:** a flat canted plate keeps only ~1/3 of its ceiling.
   - **Why it still works:** the burn from the parking orbit is only 4.34 km/s.
   - **Growth:** with k ≈ 24 at Earth and ~30 at the Sun, the cycle grows ~8x per 3S at
     9.86 R_sun, counting all impactor. At 4 R_sun it grows ~10-11x.
7. **Without a nozzle, a plate recovers less directed energy per joule.** In jet-efficiency terms
   it gets 0.4-0.6, against the 0.8-1 the parent assumes for a magnetic nozzle. With cheap reaction
   mass that is enough: the projectile bill, not Isp, is the figure of merit.
8. **The 1 Hz growth push needs a ~49 t plate, at any size** (Result 13, ADR-0064). The plate is
   set by the absorber's inertia bound `I^2 / (2 s M a)`, and the 1 g pulse is larger. A small
   flat plate with a 2 m skirt keeps 94-97% of the impulse. The peak is `~0.65 sigma v^2` at 2 m,
   so a heavy solid-steel plate admits 140-250 um sheets; 1 mm puts 5-9 GPa on the face. At
   2 Hz and 1 g the plate is ~12 t and ~100 um sheets hold if the plume spreads at u >= 8 km/s;
   JWST-style coatings do not help, because the heat comes as particles, not light.

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

| plate | k | point fireball | disk into solid charge (2x coarser grid) | real / point | disk into foam charge |
|---|---|---|---|---|---|
| parent bowl | 30 | 0.439 | **0.382** (0.345) | 0.87 | 0.129 |
| δ/D 0.5 at 1.5 F | 30 | 0.567 | **0.465** (0.442) | 0.82 | 0.208 |
| parent bowl | 100 | 0.512 | **0.348** (0.289) | 0.68 | -- |
| δ/D 0.5 at 1.5 F | 100 | 0.641 | **0.410** (0.369) | 0.64 | -- |

*Corrected 2026-10-15.* The first version of this table started the charge at the disk's
temperature rather than its pressure, giving it spurious heat (+28% of the collision's at k = 30,
+90% at k = 100). It read the real collision as 4-25% *above* the point fireball. Both runners now
start cold and assert it (ADR-0060 amendment).

- **The real collision gives less than the point fireball:** 13-18% less at k = 30 and 32-36% less
  at k = 100. The disk leaves its heat strung along the charge's axis, and a heavier charge is a
  longer cylinder. The finer grid raises the real collision 8-20%, so these ratios are if anything
  pessimistic, but the point fireball is **not** conservative.
- **The ledger now prices the real collision.** `charge.mixing_factor` interpolates the ratio in
  ln k (`charge_mixing.csv`). `mission.py` and `growth.py` apply it to the point fireball's `v_e`.
  `charge_ledger.csv` itself stays the point fireball's.
- **Foam is far worse:** 0.30-0.37 of the point fireball. A 5 m column is a line source.

**The exit flow: no nozzle, so uneven speeds.** The plate impulse comes from integrating the wall
pressure of the full 2-D flow, so gas leaving at uneven speeds and angles is already in every
number above. Separately, the gas's mass, momentum and speed moments are tallied as it crosses
the domain's outer faces. Nothing is lost (exit mass 1.00):
- the **direction factor** `P_z / Σ m|u|` is 1 if all gas leaves along the axis;
- the **speed uniformity** `(Σ m|u|)^2 / (M Σ m u^2)` is 1 if it all leaves at one speed;
- their product, `dir^2 x unif`, is what a nozzle designer would call jet efficiency.

| exhaust, k = 30 | direction | uniformity | dir^2 x unif |
|---|---|---|---|
| point fireball, bowl / paraboloid | 0.70 / 0.87 | 0.83 / 0.82 | 0.41 / 0.62 |
| disk into solid, bowl / paraboloid | 0.77 / 0.89 | 0.59 / 0.59 | 0.35 / 0.47 |
| disk into foam, bowl / paraboloid | 0.69 / 0.88 | 0.21 / 0.21 | 0.10 / 0.16 |

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
- The collar's mass and erosion are now priced (next section). The jet's true peak still needs a
  resolved, material-strength calculation.
- The energy the bore passes toward the ship (now ~5x the point fireball's share).
- The charge as a uniform equilibrated sphere. The snowplow's hot core and its early X-rays are
  not followed.
- The rigid, non-ablating plate. At ~70 km/s the `v/4` rule asks it to absorb ~18 kJ per N s.
- Peak-pressure convergence.
- The charge's thermal survival in the ship's shadow.

### The bore collar: its mass, its erosion, and its vapour as propellant -- `collar_ledger.csv`

`make xray-collision-collar` (`collar.py`; `sweep --charge-mix K collar` and `--collar-ablation`).
The disk-into-charge runs above were repeated at k = 30 and 100, on both plates and two grids,
with the plate cut into six rings from the bore's edge to the rim (0.25-0.6, 0.6-1, 1-2, 2-3.5,
3.5-6 and 6-10 m). Each ring records its impulse, the gas arriving on it head-on (column, speed,
energy) and its peak pressure. The rings tile the plate: their impulses add up to the plate's,
less the bore's, to four figures.

**The load is broad, not a needle.** The impulse per area on the inner ring is 8-10 times the
plate's mean. It is still ~5 times the mean at 1-2 m and ~3 times at 2-3.5 m, and it falls to the
mean near 6 m. The point fireball the ledger used is itself centre-heavy, 2-3 times the mean
inside 3.5 m, because the vertex is the nearest and squarest-on part of the bowl. The jet roughly
doubles that again inside 2 m.

**Sizing by impulse, not by peak.** The pulse lasts microseconds, far shorter than any structural
response, so the ring is damaged by its impulse (the impulsive branch of a pressure-impulse
diagram), not by its 1-20 GPa peak. The collar is sized to be *impulse-matched*. Its areal mass is
the plate's 30 kg/m^2 scaled by the ring's impulse over a reference, so it takes the same velocity
kick per pulse as that reference. Two references bracket it:
- **Against the point fireball on the same ring:** this charges the jet's increment alone, for a
  plate already built to the smooth load the ledger used.
- **Against the plate's mean:** this assumes no structural spreading, and also charges the smooth
  fireball's own centre-heavy load.

All the impulses scale with the disk's mass, so neither collar mass does.

Impulse converges where peak pressure does not. On the coarse grid the inner rings' impulse
changes by 10-60%, while their peak changes up to six times.

**Erosion and thrust.** Each ring with a real arriving column (r < 2-3.5 m, at least 1e-4 kg/m^2
per kg of disk) is re-run as the study's 1-D radiation-hydrodynamic bounce on the CH2 table:
- once against a rigid wall;
- once against graphite that boils off at `q_wall / Q*` and is injected into the gas (the f(v)
  study's Rung E, ADR-0014).

The 1-D slug carries the ring's arriving column and speed. Its density is set to give the ring's
2-D peak. Further out the rings are loaded by gas sliding along the wall, not arriving on it, so
their erosion is the plate-wide ablation item.

`Q*` is bracketed between graphite's sublimation (59 MJ/kg) and a third of it. The f(v) study's
vapour curtain is switched off, because it creates energy here. It was calibrated for ablated
mass well below the cloud's, and the inner ring boils off ~10 times its arriving column. Each row
is checked: `(J - p_in)^2 / 2M` must not exceed the slug's kinetic energy.

The 1-D slug arrives cold, but in 2-D the jet carries 3-36 times its kinetic energy as enthalpy.
The erosion is therefore bracketed:
- **Floor:** the cold slug at 59 MJ/kg.
- **Ceiling:** 20 MJ/kg, scaled up by that enthalpy ratio.

**The vapour is credited with only its own push.** That is the 1-D slug's extra impulse per area.
Most of a ring's 2-D impulse is pressure from gas already stagnated against the wall. Scaling it
by the 1-D gain gives the vapour an impossible exhaust speed (~0.6-2.7 Mm/s).

The ledger for a 100 t final mass (production grid; jet-only coarse-grid values in brackets):

| plate, k | collar, jet only | collar, vs the mean | eroded over the boost, floor-ceiling | worst-ring depth | impulse the vapour adds | mass ratio change, floor / ceiling |
|---|---|---|---|---|---|---|
| bowl, 30 | +0.73 t (0.50) | +3.6 t to 6 m | 0.25-9.2 t | 0.10-2.3 m | 0.02-0.4% | +0.00 / +0.12 |
| paraboloid 1.5 F, 30 | +0.53 t (0.26) | +4.0 t | 0.08-5.0 t | 0.03-1.2 m | 0.2-0.4% | -0.00 / +0.06 |
| bowl, 100 | +0.29 t (0.12) | +3.8 t | 0.68-9.2 t | 0.11-1.7 m | 1.3% | -0.03 / +0.11 |
| paraboloid 1.5 F, 100 | +0.21 t (0.05) | +4.1 t | 0.45-5.5 t | 0.12-1.2 m | 0.4-1.2% | -0.00 / +0.06 |

*Corrected 2026-10-15.* The ring loads come from the disk-into-charge runner, which first started
the charge as hot as the disk (ADR-0060 amendment). This is the rerun. The earlier version had
~1 t for the jet and 6-11 t of ceiling erosion. Values in brackets are the coarse grid.

- **The jet itself costs ~0.2-0.7 t of structure**, centre-heavy out to 2-3.5 m. That is under 1%
  of a 100 t ship, and it scales with no disk size. It rises on refinement (0.05-0.5 t coarse), so
  it is not converged: ~0.5-1 t is the honest figure.
- **The smooth fireball's centre load costs more.** The conservative bound, which also grades the
  plate for it, is **+3.6-4.1 t** out to 6 m, of which 2.0-3.1 t is the smooth fireball's own
  load. This is dry mass, so it costs payload, not mass ratio. **The plate should be graded, and
  the jet is a modest addition.**
- **Erosion is a consumable, concentrated on the inner ring.** The floor eats 0.1-0.7 t over the
  boost, a few cm to 12 cm of graphite on the worst ring. The hot-jet ceiling eats 5-9 t, 1.2-2.3 m
  there: a fed sleeve at the ceiling, a thick fixed ring at the floor.
- **The vapour adds little either way.** It leaves at 3-68 km/s depending on the ring's flow, and
  adds 0.02-1.3% to the pulse's impulse. Charged as propellant, the eroded collar moves the mass
  ratio by -0.03 to +0.12, never past 3. (These mass ratios are the point-fireball ledger's; the
  changes are what the collar adds.)

Not modelled here: mechanical erosion (spall and crushing; a bonded, backed collar has no free
rear face to spall), the collar's own heat soak between pulses, the jet's true resolved peak, and
the wall flow's erosion of the rest of the plate.

## Result 11: the Earth-departure overtake pusher -- `overtake_ledger.csv`

`make xray-collision-overtake` (`overtake.py`, `sweep --overtake-2d`, ADR-0061). This is the
departure wave at Earth. An impactor arriving at **100 km/s from behind** the departing ship meets
a polyethylene charge held a few metres off the plate. The boost is **11 km/s** to a high
apoapsis. The target is a propellant fraction of 1/3 (mass ratio 1.5), which needs `v_e` of at
least 27.1 km/s.

**What changes from Result 10.**
- **No bore.** The merged mass drifts *onto* the plate at `v/(1+k)`. The impactor's own momentum
  is delivered as well as the expansion's, `J ≈ m v (1 + η √(1+k))`.
- **The ratio is still scale-free.** `drift/u_rms = 1/√k` at any closing speed, so the 2-D runs
  carry over. 100 km/s enters only through the one-zone fireball, the rebound swap
  (`e_eff` 0.67-0.77 at the 20-45 km/s arrival speed) and the bond toll.
- **The fireball is dark again.** It is ~10^6-10^7 optical depths thick and radiates ~1e-7 of its
  heat.
- **The chemistry is kinder.** Ionization (11-18% of the heat) recombines completely. The gas then
  cools below 2500 K before the plate, and H2 re-forms with Damköhler numbers of 3e2-9e5, two to
  six decades of margin on the H + H + H rate. Only carbon's share of the bonds is in doubt. It
  is bracketed between condensed (no toll) and atomic, a velocity loss of 3% at k = 5, 6% at
  k = 10 and 12% at k = 20.

Three geometries were run on two grids (0.05 m and 0.1 m cells):
- **Point fireball**, k = 5-40, on the parent bowl, the δ/D 0.5 paraboloid at 1.5 F, and a flat
  plate with a 4 m skirt (charge 2.7 m off).
- **Compact disk-into-charge collision**, the point fireball's check.
- **Planar collision**: a wide impactor meeting a wide polyethylene sheet of the same radius
  (4 m or 7 m) 2 m in front of the flat skirted plate.

| plate | geometry | k | `J/(m v)` | `η` | mass ratio (carbon atomic) | max impactor per pulse (400 MPa) |
|---|---|---|---|---|---|---|
| parent bowl | point | 5 / 10 / 20 / 40 | 2.53 / 3.05 / 3.82 / 4.93 | 0.62 at k 10 | 1.27 / **1.49** / 1.94 / 3.23 | 4-10 kg |
| paraboloid 1.5 F | point | 5 / 10 / 20 / 40 | 2.81 / 3.45 / 4.38 / 5.72 | 0.74 at k 10 | 1.24 / **1.42** / 1.79 / 2.76 | 3-6 kg |
| flat + 4 m skirt | point | 10 / 20 | 3.04 / 3.79 | 0.61 | 1.49 / 1.97 | ~2 kg |
| parent bowl | compact | 10 / 20 | 2.88 / 3.24 | 0.57 | 1.53 / 2.17 | 2.0 / 15.7 kg |
| flat + 4 m skirt | 4 m sheet | 10 / 20 | 3.03 / 3.69 | 0.61 | 1.50 / 2.00 | 3.6 / 4.4 kg |
| flat + 4 m skirt | 7 m sheet | 10 / 20 | 3.10 / 3.79 | 0.63 | 1.48 / 1.94 | 10.9 / 13.6 kg |

`η = (J/mv - 1)/√(1+k)` is the share of the expansion's ideal momentum the plate keeps.

**Mass ratio 1.5 needs k ≈ 10-12** (power-law fit over the point rows; the first value has carbon
left atomic, the second has the bonds returned):
- parent bowl: k ≈ 10.3-11.6;
- deep paraboloid at 1.5 F: k ≈ 12.0-14.1;
- flat + skirt: k ≈ 10.1-11.4.

The impactor bill is then 0.042-0.049 kg per kg of final mass. At $1000/kg that is $42-49 per kg
of final mass, or $1000/k ≈ $70-100 per kg of propellant. k = 40 would cost $25/kg of propellant,
but its mass ratio is 2.1-3.3: at 122 MJ/kg the gas is too cool, and its bonds are too large a
share of its heat.

- **The overtake couples better than the departing case.** `η` is 0.61-0.74 here, against
  0.40-0.55 when Result 10's fireball drifted away from the plate.
- **The planar sheets match the point fireball; the compact charge falls short.** On the flat
  skirted plate a 7 m sheet gives 0-2% more impulse than the point fireball. The compact charge
  on the bowl gives 6% less at k = 10 and 15% less at k = 20: it leaves its heat strung along
  the axis, as at 400 km/s. (The first version had the compact charge at 1-4% and the sheets 3-6%
  above it. Both runs had started the charge as hot as the impactor; ADR-0060 amendment.)
- **Planar coupling greatly helps the plate's loading.** The peak scales with the impactor's mass
  per area, so a sheet in front of the flat plate takes 0.071-0.088 kg/m^2 at 400 MPa whatever its
  radius. That is ~11-14 kg on a 7 m sheet against ~2 kg as a compact charge. Inside the
  paraboloid it is lower, 0.02-0.07 kg/m^2.
- **The ledgers carry the compact factor.** `mission.overtake_mixing` applies the measured
  compact-over-point ratio to the breeding cycle's Earth stage.
- **The grid study is conservative.** The finer grid gives 2-4% more impulse for the point
  fireball and 5-8% more for the sheets. The compact and sheet peaks rise on refinement, as at
  400 km/s, so their per-pulse limits are upper values; ~0.05 kg/m^2 is the safer sheet figure.

**Confinement follow-up** (`sweep --overtake-2d confine`, `overtake_2d_confine.jsonl`). A perfect
plate, with every parcel leaving along the axis at one speed, would get `J = m v (1 + √(1+k))`.
That ceiling alone caps mass ratio 1.5 at k ≈ 21. The runs get 68-80% of it (column
`ceiling_share`). To find how much of the gap confinement can close, these were added at
k = 10, 15 and 20:
- skirts of 4 m and 8 m on the δ/D 0.5 paraboloid;
- a deeper δ/D 0.75 paraboloid, with the charge at 1.5 F and 2 F;
- 4 m and 5.5 m sheets inside the δ/D 0.5 paraboloid.

| plate (k = 10) | ceiling share | direction | uniformity | k for mass ratio 1.5 | impactor per pulse |
|---|---|---|---|---|---|
| δ/D 0.5 paraboloid, 1.5 F | 80% | 0.90 | 0.83 | 12.0-14.1 | ~4 kg |
| + 4 m or 8 m skirt | 80% | 0.96-0.97 | 0.76-0.78 | 12.2-13.9 | ~4 kg |
| δ/D 0.75 paraboloid, 2 F | **84%** | 0.96 | 0.83 | **13.0-14.9** | 3.5-4.5 kg |
| δ/D 0.75 paraboloid, 1.5 F | 84% | 0.96 | 0.81 | 13.0-14.9 | ~1.3 kg |
| 4-5.5 m sheet inside the δ/D 0.5 paraboloid | 75-76% | 0.90-0.91 | 0.66-0.70 | 11.1-12.4 | ~2-2.5 kg |

- **Direction is solved; speed spread is the rest.** Deep bowls and skirts send the exhaust out at
  a direction factor of 0.96. What remains is the spread of exit speeds (uniformity ~0.8), which
  confinement does not touch. Skirts straighten the flow but widen the spread, so the impulse is
  unchanged. A sheet in the bowl is slightly worse than a compact charge.
- **The deeper bowl gains ~1 in k but costs plate.** δ/D 0.75 has 712 m^2 of plate against 533 m^2,
  ~5 t more at 30 kg/m^2. It saves ~0.003 kg of impactor per kg of final mass, ~0.3 t on a 100 t
  ship. **The δ/D 0.5 paraboloid at 1.5 F, k ≈ 12, stays the recommendation.** Its charge at 2 F
  rather than 1.5 F matters for the deeper bowl: the impulse is the same, with three times the
  impactor per pulse.
- **The next lever is the speed spread**, for instance a layered charge whose hot face pushes its
  cold bulk like a piston. It is not modelled.
- **The grid study holds.** The finer grid gives 1.5-6% more impulse than the coarse one for every
  case.

**The impactor sheet** (`overtake_sheet.csv`, `overtake.sheet_budget`). At 100 km/s the impactor
carries ~5 GJ/kg, a thousand times any binding energy, so its material does not matter to the
coupling. Seth chose polyimide (Kapton).
- **Mass per area:** the flat-plate sheets' peak gives 0.071-0.088 kg/m^2 at 400 MPa on the
  production grid (column `max_areal`). It is an upper value, since the peak rises on refinement;
  0.05 kg/m^2 is the safer figure. That is 35-62 µm of polyimide or Kevlar, or 31-55 µm of
  carbon-fibre laminate.
  The ship's polyethylene sheet at k = 10 is ~0.5-0.8 mm thick.
- **Survival from 0.5 to 1 AU:** a black sheet facing the Sun and radiating from both faces runs
  at 468 K at 0.5 AU and 331 K at 1 AU.
  - Polyimide is proven as a solar sail (IKAROS, spin-deployed, 7.5 µm) and good to ~400 °C.
  - Carbon fibre needs a dry or polyimide matrix, not epoxy.
  - Kevlar needs a UV-blocking coating.
- **The sheet is a solar sail.** At 0.5 AU it accelerates at 0.23-0.73 mm/s^2, from absorbing to
  perfectly reflecting (0.21-0.41 mm/s^2 at 0.088 kg/m^2). Over an 8.7-day coast from 0.5 to 1 AU
  (at an assumed 100 km/s outward) that is **29,000-100,000 km of drift**. It is predictable only with a known attitude, so the
  sheet spins to stay flat, flies face-on to its direction of travel (so the push mostly shifts
  its arrival time), and is tracked and corrected late.
- **Tolerances are loose.** A 10 cm sag spreads the arrival by 1.0 µs, and a 1° tilt across 7 m by
  2.4 µs, against a ~100 µs pulse.

Open:
- The effective-gamma gas. Polyethylene's real EOS through dissociation, and the soot that carbon
  may condense into, are not followed in 2-D; the toll is bracketed instead.
- The sheet-on-sheet collision's peak (under-resolved), and how a non-ideal sheet behaves: sag,
  holes, a spin-flattened cone.
- The plate's own ablation at 20-45 km/s, and the charge sheet's deployment and survival on the
  ship.

## Result 12: the canted departure plate and the 3S growth cycle -- `growth_cycle.csv`

`make xray-collision-growth` (`growth.py`, `euler2d::canted`, `sweep --canted-2d`, ADR-0063).

**The cycle.** The Earth wave now launches ships onto the calcs repo's 3S Jovian cycle: Earth ->
Jupiter -> dive periapsis -> Earth, 3.28 years. The closures were run in that repo
(`solve_synodic_closure`, `departure_nozzle_ledger`) and are copied as cited inputs:

| closure | departure burn | stream vs thrust axis | stream closing speed | periapsis boost | Sun-wave closing speed |
|---|---|---|---|---|---|
| 9.86 R_sun, 100 km/s return | 4.34 km/s | 100-102° | 103 km/s | 20.45 km/s | 392 km/s |
| 4 R_sun, 150 km/s return | 4.42 km/s | 94-96° | 158 km/s | 35.04 km/s | 617 km/s |

The departure burn is computed here (`mission.departure_dv`, vis-viva from the 10-day orbit's
600 km periapsis); it reproduces the calcs repo's ADR 0072 correction of the 4.15 / 4.23 km/s
first copied.

- **The burn is small.** It starts at the 10-day parking orbit's periapsis, near escape speed, and
  Jupiter only has to stop the ship, not reverse it.
- **The stream arrives almost sideways**, ~10° *away* from the thrust direction. It cannot rebound
  off the plate. The impulse has to come from the merged fireball exploding beside a flat plate as
  it slides past.

**The canted ceiling.** A frictionless flat plate takes no tangential force. Conserving the
stream's tangential momentum and spending all the energy normally gives

`J_n <= m w (cos theta + sqrt(1 + k - sin^2 theta))`,

which is the overtake ceiling at θ = 0.

**The canted plate in planar 2-D** (`euler2d::canted`). This is a strip model: wide sheets, a flat
20 m plate, and the charge 4 m in at a 1-4 m standoff. The share of the canted ceiling, best at the
1 m standoff:

| charge, 1 m across | k = 5 | 10 | 30 | 60 |
|---|---|---|---|---|
| thin slab, 0.3 m along the stream | 0.30 | 0.29 | 0.26 | 0.24 |
| **square block (the design)** | **0.37** | **0.36** | **0.32** | **0.29** |
| round | | 0.34 | 0.30 | 0.27 |
| half-cylinder: flat face to the sheet / dome to the plate / dome up | | 0.31 / **0.38** / 0.32 | 0.28 / **0.33** / 0.27 | 0.25 / **0.30** / 0.23 |
| ideal point fireball (the merged state) | 0.42 | 0.43 | 0.43 | 0.44 |
| head-on reference, point fireball | | 0.53 | 0.49 | |

- **The sideways geometry itself costs ~20%** against head-on in the same model (0.43 against 0.53
  for the point fireball).
- **The rest is how the collision leaves its heat.** A sheet-on-charge strike leaves hot gas
  strung along the stream, so it expands along the plate rather than into it.
- **The half-cylinder with its dome to the plate edges out the square block** (+6% at k = 10,
  +1-3% at k = 30-60, within ~2x the grid difference). The round charge and the half-cylinder
  with its flat face to the sheet do worse: they give the sheet's edges too thin a column. Dome
  up is worst, because its mass sits farther from the plate.
- **Corrected 2026-10-15.** The first shape runs aimed a sheet 1 m wide at the charge's centre,
  so ~41% of it passed beside the two half-cylinders lying on a face. They read 0.23-0.28 at
  k = 10. The sheet now spans what it strikes. The square, slab, round and flat-face runs are
  unchanged, and so is the growth optimizer, which uses the square block.
- **The torque fix is placement.** The impulse centroid lies 0.5-1.3 m downstream of the charge at
  the 1 m standoff, so the ship's centre of mass goes there, or the charge is graded to match.
- **Grids agree to 3%** (5% for the half-cylinders).

**A methane-faced charge** (`--canted-2d layers`, `canted_2d_layers.jsonl`, Seth's question). A
dense layer on the charge's plate side (solid methane, ~520 kg/m^3) holds part of k, under a
polyethylene foam, and the sheet strikes the foam alone. Methane vaporizes at ~0.6 MJ/kg, so it is
a safe home for mass that a solid polyethylene charge might leave unvaporized. In the gamma-gas
model only the density differs between the layers. Share of the canted ceiling at the 1 m standoff:

| charge | layer share / density ratio | k = 10 | k = 30 |
|---|---|---|---|
| square block | none | 0.361 | 0.322 |
| square block | 0.5 / 5 | 0.357 | 0.330 |
| half-cylinder, dome to the plate | none | 0.383 | 0.326 |
| **half-cylinder, dome to the plate** | **0.25-0.5 / 5** | **0.397-0.399** | **0.337-0.346** |
| half-cylinder, dome to the plate | 0.25-0.5 / 15 | 0.368-0.379 | 0.313-0.326 |
| half-cylinder, dome up (methane on the flat face) | 0.5 / 5 | 0.345 | 0.295 |

- **The best layered charge gains ~10% (k = 10) and ~7% (k = 30) over the square block.** That is
  the dome-to-the-plate half-cylinder with a thick methane layer (density ratio 5, i.e. a foam of
  ~100 kg/m^3). About half the gain is the shape and half the layer.
- **A thin, dense layer (ratio 15) does nothing or harm.** The methane has to be a thick cushion,
  not a skin.
- **Methane on the flat face of a dome-up charge loses** to the plain square.
- **The gains are a few times the grid difference, but small.** The coarse grid ranks them the
  same. The growth optimizer has not been rerun with this charge.
- **Not modelled:** methane's own EOS, rebound and bonds (both layers are one gamma gas), and
  whether the foam's fireball mixes into the methane.

**The optimizer.**
- **Earth stage:** the canted share (fitted in ln k), times the canted ceiling at the closure's
  angle and speed. Polyethylene's rebound is swapped in at the fireball's `u_rms`. The bond toll
  comes from `freeze.freeze` at the stream's speed.
- **Sun stage:** the δ/D 0.5 paraboloid's `v_e(k)` from `charge_ledger.csv`, carried from 400 km/s
  to the dive's closing speed with its rebound and bond toll re-evaluated.
- **Objective:** Earth-wave impactor returned per kg of *all* impactor spent (Earth and Sun wave)
  per 3S, maximized over `k_E` (5-60) and `k_S` (20-100).

| closure | Earth charge | `k_E` | `k_S` | Earth / Sun mass ratio | departs per kg of Earth impactor | **growth per 3S** (per total impactor) | per Earth-wave kg | doubling |
|---|---|---|---|---|---|---|---|---|
| 9.86 R_sun | **square** | **24** | **30** | 2.13 / 1.98 | 21.3 kg | **7.8-8.0** | 10.5-10.8 | 1.1 yr |
| 9.86 R_sun | point (ideal) | 34 | 30 | 2.07 / 1.98 | 32.2 kg | 10.4-10.6 | 15.8-16.3 | 1.0 yr |
| 4 R_sun | **square** | **52** | **30** | 2.13 / 2.33 | 45.8 kg | **10.3-10.6** | 19.2-19.7 | 1.0 yr |
| 4 R_sun | point (ideal) | 60* | 40 | 1.78 / 2.59 | 77.2 kg | 13.4-13.7 | 29.1-29.8 | 0.9 yr |

How to read the table:
- The ranges run from the Sun ship's plate kept as dead mass (2.5% of the ship, the calcs repo's
  ADR 0072 sizing; ADR-0062 amendment) to the plate released as
  impactor.
- Polyethylene's carbon is left atomic.
- The Sun stage is always the real disk-into-charge collision (`charge.mixing_factor`), and that
  keeps `k_S` at ~30. Above it, the real collision's coupling falls faster than the projectile
  bill.
- If the carbon condenses before the plate, growth rises 10-27%. Those optima sit on the `k_E` = 60
  edge of the simulated range, and are flagged rather than extrapolated.
- \* = an optimum on the edge of the range.
- The Sun wave is 27-55% of all the impactor spent.

**Conclusions.**
- **The departure plate is the weak link, and the cycle still closes.** The plate keeps only ~1/3
  of its ceiling (head-on plates kept 70-84%), but the burn is only 4.34 km/s.
- **The real optimum is k ≈ 24 at Earth and ~30 at the Sun at 9.86 R_sun**, growing ~8x per 3S
  (doubling in ~1.1 years).
- **4 R_sun is worth ~30% more growth**, ~10-11x, at k ≈ 52 and ~30. Its faster stream raises the
  Earth-side k, but its larger boost needs ~half of all impactor to be Sun wave.
- **Per Earth-wave kilogram the square charge returns 10.5-19.7 kg.** The calcs repo's 4 R_sun
  ledger, with a magnetic nozzle from the same parking orbit, quotes 23.0.
- **The levers left are the charge and the energy.**
  - A charge whose merged state is rounder approaches the point fireball's 0.43 share, worth
    +30% growth.
  - A bond-free charge lifts the `k_E` cap.

**The Jupiter split** (`mission.jupiter_split`, ADR-0063 amendment, `mission_jupiter_split.csv`).
The Sun-wave impactors travel with the ship from Earth and split from it at Jupiter. The ship
falls into a prograde dive and the impactors into its mirror-image retrograde dive, so the two
meet head-on at periapsis. A flyby cannot change the size of the excess, so the impactors' carrier
burns at perijove to raise it from 11.70 to the true mirror excess, 15.07 km/s (9.86 R_sun), with
Jupiter's Oberth boost:
- 0.77 km/s at 1.1 R_J and 1.67 km/s at 6.5 R_J, near the 6.72 R_J above which the flyby cannot
  turn the disks;
- 0.45-1.13% of the ship in propellant at 330 s, because the impactors are only ~1/60 of the ship;
- at 4 R_sun (12.17 -> 14.34 km/s), 0.49-1.07 km/s.

The mirror excess and the turn ceiling are the calcs repo's (ADR 0072, 2026-10-16); the first
version burned only to the retrograde floor (14.78 / 14.16 km/s).

**The round of asks G1-G10 and E1-E6 is closed** (2026-10-16; ADR-0063 amendment). The calcs
repo answered G1-G4 and G10 in its ADR 0072, and the parent E1-E6 and G5-G9
(`docs/solar_environment_answers_from_parent.md`):
- G1 and G4 confirmed; the departure burn and the mirror excess corrected (above);
- G5: the Sun ship's plate is kept as dead mass, 1-2.5% of the ship;
- G6: sheets coast furled, unfurl ~100 s out and aim off the plate until the unfurl is confirmed;
  guidance is 2-10% of sheet mass;
- G7: the disks ride with the ship and the carrier burns 1.67 km/s at 6.5 R_J;
- G8: the disk survives at 4 R_sun (~2430-2890 K face-on, ~700 K edge-on); its guidance needs a
  shade and ~1 g/cm^2 of shielding;
- G9: polyethylene stays the charge for now;
- E1-E6: steady sunlight stays the design heat load at both depths. SEP dose (>= 100 krad parts,
  ~1% of passes lost at solar maximum) and the CME push on the disks (flown edge-on) are the
  design drivers, and neither is this repository's physics.

Left with the calcs repo: the departure charged from the pad (G3, ~6.6x), and its ledger's optima
on this repository's table edges (k_E 5-6, k_S 20).

Not modelled:
- real 3-D sheets and charges (the strips are planar);
- the plate's ablation at 15-30 km/s arrival speeds;
- the canted plate's peak pressure, which sets the impactor per pulse.

## Result 13: the plate for the 1 Hz growth push, and how thick its sheets may be -- `cadence_plate.csv`, `cadence_sheet.csv`

The parent (`sec:charge_pusher_plume`, asks S36-S37 in its
`docs/plume_cadence_asks_for_impact_sim.md`) slows the 400 km overtake push from 4 Hz at 12 MN s
to 1 Hz at 14.7 MN s (~1 g on the 1500 t unit), so that the 100 km/s Kapton sheets survive the
gas of earlier pulses. Seth asked what that does to the plate, whether a flat sheet's even load
lightens it, and whether thicker sheets on a smaller plate pay. ADR-0064;
`make xray-collision-cadence-plate` (`sweep --sheet-plate-2d`, `cadence_plate.py`).

**The 2-D runs** (`sheet_plate_2d.jsonl`, two grids). The sheet and its charge sheet merge at the
flown k = 6.74 and drift onto a flat plate. There are three thick-sheet radii (2.5 / 3.6 / 5 m),
three standoffs (2 / 4 / 8 m) and three plates each (sheet + 1.5 m bare, sheet + 4 m bare,
sheet + 1.5 m with a 2 m skirt). The reference is the flown 7 m sheet on the 20 m plate with its
4 m skirt (J/mv 2.81).
- **The peak is `c(h) sigma v^2`, set by the sheet's mass per area.** It does not depend on the
  sheet's radius: c = 0.646 / 0.496 / 0.349 at 2 / 4 / 8 m on the production grid, with a spread
  under 1% at 2 m. It is a lower value, since the coarse grid gives 0.47 / 0.36 / 0.22.
- **A small flat plate that keeps a skirt keeps the impulse.** At 2 m standoff, sheet + 1.5 m with
  a 2 m skirt keeps 94-97% of the reference. Bare, the same plate keeps 72-82%, and a plate 4 m
  wider still only 77-86%: the bare plate throws the gas sideways (direction 0.57-0.69).
- **Standoff costs impulse.** At 8 m the skirted small plates keep 73-81%.

**The plate's mass is the inertia bound, not survival.** Each pulse kicks the plate, and the
absorber must stop it within its stroke `s` under the cargo's peak `a`:
`m_p >= I^2 / (2 s M a)`, with `M` the unit at the end of the push. This depends on the impulse
per pulse, not on the plate's size or the mean thrust. At a fixed mean, `I` grows as 1/cadence.

| push | sheet per pulse | calcs floor (E^(2/3)) | inertia at 100 g: 2 / 4 / 10 m stroke |
|---|---|---|---|
| 4 Hz, 12 MN s (priced) | 46.4 kg | 40.7 t | 32.4 / 16.2 / 6.5 t |
| 1 Hz, 14.7 MN s (1 g) | 56.8 kg | 46.5 t | **48.8** / 24.4 / 9.8 t |

- Slowing to 1 g lowers the force but raises `I^2` 1.5x, so the plate gets heavier.
- The flat sheet's even load removes bending, which retires the calcs repo's survival floor, a
  point-blast law (30 kg/m^2 scaled as pulse energy^(2/3)). The inertia bound was already larger.
- **A smaller plate is not lighter, only thicker.** 48.8 t over a 5-20 m plate is 4-43 mm of
  steel. The lever is the stroke.
- The 100 g cap is the calcs repo's `G_CAP` for the growth units' cargo. The parent's 16 g case
  (on a 10 m stroke) needs 61 t at the end-of-push mass.

**The sheet's thickness has a floor and a ceiling.**
- **Floor:** the parent's plume dose (full absorption) held inside Kapton's 400 C rating:
  - 1 Hz: 82 / 33 / 14 um at spread u = 5 / 8 / 12 km/s;
  - 4 Hz: 1.1 mm / 435 / 197 um.
- **Ceiling:** the face pressure on a face that holds it:
  - 35-62 um: 0.3-0.6 GPa, the sandwich facesheet;
  - 140 um: 1.3 GPa at 2 m, HY-100 / 4340;
  - 250 um: 2.3 GPa at 2 m (1.2 at 8 m), maraging 300;
  - 0.5 mm: 4.7 / 2.5 GPa at 2 / 8 m, no face (allowables from
    `docs/spray_plate_steel_face_limits.md`, without prestress);
  - 1 mm (the parent's S37.3 option): 9.2 / 5.0 GPa, above every steel and in the range
    `docs/spray_plate_fatigue_and_spall.md` found fatal in fatigue.
- **So at 1 Hz, 140-250 um sheets** give a 2-3x margin over the slow-spread plume. They need a
  solid hard-steel face and a 7-10 m plate with a 2 m skirt, at the same ~49 t (12-18 mm of
  steel). Because the plate must be ~49 t anyway, making it solid steel costs nothing, and that is
  what admits the thicker sheet.

Not settled:
- **The peaks grow ~37% from the coarse grid to the production one.** The grid's sheets are ~0.2 m
  thick against the real millimetres, so the ceilings may be a third lower.
- **The face's rise time is not resolved.** The face's pressure must rise slower than one round
  trip of a stress wave through the steel, 4-15 us for these plates.
- **The plume dose is the parent's model.** Its spread `u` (S36) is not fitted to our gas.

### 2 Hz at 1 g, and a JWST-style coated sheet -- `coated_sheet.csv`, `coated_sheet_dose.csv`

Seth asked whether a thicker sheet with JWST's coatings (an emissive doped-silicon face over
reflecting aluminium or silver, on ~100 um of Kapton) flies through the old plume at 2 Hz and 1 g.
`coated_sheet.py` re-implements the parent's expanding-cloud model (`eq:jovian_plume_density`,
summed over the train). It reproduces the parent's doses: 633/251/113 against its 635/251/114
kJ/m^2 at 4 Hz, and 47/19/8 exactly at 1 Hz. The sheet is 1-D conduction through the Kapton, with
each face radiating.

**The dose at 2 Hz and 1 g is twice the 1 Hz dose: 96 / 38 / 17 kJ/m^2** at u = 5 / 8 / 12 km/s.
Each pulse carries half the gas (220 kg, 7.35 MN s), but the sheet meets the clouds twice as young.
At a fixed mean thrust the dose goes as the cadence. It arrives as a pulse peaking at 375 kW/m^2
about 96 ms before impact (u = 8).

**The coating does nothing useful.**
- **The metal film cannot reflect this heat.** The heat arrives as ~50 eV/amu particles that stop
  in the top ~0.1 um. A metal film reflects photons, and a 50-150 nm film under the silicon sits at
  the face's temperature within nanoseconds.
- **Radiation is too slow.** Even ideal emitters on both faces lower the face's peak by only
  10-50 K. At 600 K a face radiates ~7 kW/m^2, against the 375 kW/m^2 pulse.
- **JWST's own stack is 3-7 K worse than bare Kapton**, because its aluminised back face barely
  radiates.

**100 um at 2 Hz, 1 g** (face peak; Kapton's rating is 673 K):

| u | bare | Si / Al / Kapton | ideal emitters | no radiation |
|---|---|---|---|---|
| 12 km/s | 406 K | 406 K | 403 K | 420 K |
| 8 km/s (reference) | 572 K | 572 K | 568 K | 589 K |
| 5 km/s | 1095 K: all of it past 800 K ~100 ms out | 1094 K | 1086 K | 1139 K |

- **Thickness does not fix the slow spread.** The face heats faster than heat can soak into the
  sheet: at u = 5, a 250 um sheet's face still reaches 985 K, passing 400 C ~102 ms out. The skin
  past 800 K is only 31-38 um, though, and the rest of the sheet stays intact (mean 537-692 K).
- **A sacrificial 25 um shield film ahead of the sheet, the JWST lesson of gaps, buys little.** It
  fails ~208 ms out at u = 5, before the flux peaks, and the main sheet's face then reaches
  940-990 K. At u = 8 it holds the main sheet's face to ~445 K.
- **At 1 Hz, 1 g, 100 um holds at every spread** (619 K at u = 5).

**The plate at 2 Hz and 1 g is a quarter of the 1 Hz plate**, since `I` halves: 12.2 t at 100 g on
a 2 m stroke. With 100 um sheets (28.4 kg per pulse, an 8 m sheet) the 2 m peak is 0.92 GPa,
inside HY-100's allowable, on a 9.5 m plate with a 2 m skirt; the 30 kg/m^2 floor and the
inertia bound then meet at ~12 t (3.9 mm of steel). The lob and booster costs are the 1 Hz push's,
because the mean thrust is the same 1 g.

**So 2 Hz at 1 g works with ~100 um sheets if the plume spreads at u >= 8 km/s**, the parent's
reference, and needs no coating. If u is nearer 5 km/s, every sheet loses its face over the last
~100 ms. Whether a sheet with a decomposed 30-40 um face still works as the impactor is the
parent's S37.2 question, and is not settled here. The spread `u` (S36) therefore decides 2 Hz.
The doses are the parent's full-absorption bound; Kapton's c_p is held at its room value, which
overstates the rise; and the 800 K decomposition threshold is assumed.

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
