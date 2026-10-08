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
`fireball.collide(material, impactor_radius, mass_ratio, v_rel, lowering=, kappa_scale=)`.

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

## Limits

One zone, so there is no cool photosphere and no temperature gradient. A collision of two balls
is a pancake, not a sphere. LTE is assumed throughout, including the thin, cold tail. The
lowered EOS is a bracket, not a model.
