# Pass-through pusher at 66 km/s: argon clouds radiate, but cannot be used

Handoff to `Balloon-Pulse-Propulsion` (parent), §`why_not_xray_drive` and §`pass_through_pusher`.
It is written to be copied there. It answers whether the pass-through pusher could run at the
rocket's own collision speed, before any Sun dive. It replaces the "67 km/s, for completeness"
note in [`xray_pusher_handoff.md`](xray_pusher_handoff.md). The evidence is
[`xray_collision_study.md`](xray_collision_study.md) Result 6 (`make xray-collision-argon-bags`,
`data/results/xray_collision/argon_bags_k2.csv`).

**Status: provisional.** These are one-zone fireballs, with a Saha EOS and a new 120-group TOPS
(LANL OPLIB) pull for argon. The light's escape is followed group by group. Expect factor-of-
order-unity errors in geometry, not order-of-magnitude ones.

## Conclusion for the paper

At 66 km/s, light carries most of a collision's heat only from an optically thin argon cloud of
about a gram per cubic metre, and such a cloud is tens of metres across. A projectile cannot be
spread that thin and still pass through a narrow bore and the radiation case's entry hole. A
compact one that can pass needs a target column hundreds of metres to kilometres long. The
pass-through pusher therefore does not work at the rocket's collision speed. It needs the Sun
dive, where 618 km/s puts about 130 times more energy into each kilogram, before it makes
sense.

## What radiates at 66 km/s

One kilogram of argon at 66 km/s strikes two kilograms at rest in the ship's frame (k = 2). Its
2.18 GJ splits into 0.73 GJ of centre-of-mass motion, which no design recovers, and 1.45 GJ of
heat. Both bodies are taken at one uniform density at the moment of impact.

| argon density at impact | projectile cloud, 1 kg (diameter) | target cloud, 2 kg (diameter) | heat radiated |
|---|---|---|---|
| 16 kg/m^3 (an intact bag at 10 atm) | 0.5 m | 0.6 m | 9% |
| 0.1 kg/m^3 | 2.7 m | 3.4 m | 35% |
| 10 g/m^3 | 5.8 m | 7.3 m | 52-54% |
| 3 g/m^3 | 8.6 m | 11 m | 66-68% |
| **1 g/m^3** | **12 m** | **16 m** | **86-88%** |
| 0.1 g/m^3 | 27 m | 34 m | 93-96% |

The 16 kg/m^3 and 0.1 kg/m^3 rows are gray (10 kg bags; cloud sizes scaled to 1 kg). The rest
use the multigroup model, for a 1 kg projectile except at 3 g/m^3, which is the 10 kg value.

- **More than half the heat leaves as light below about 10 g/m^3, and 80% at about 1 g/m^3.**
  That density is a millionth of liquid argon's, and like a gas at 60 Pa.
- **The light is not trapped.** Cold argon absorbs strongly above its 15.8 eV ionization edge but
  is clear below its 11.6 eV resonance lines. The fireball emits at kT ~3-4 eV, so more than half
  of the heat leaves below 11.6 eV. Even if the unshocked target absorbed every harder photon and
  kept the energy, 65-70% would still leave at 1 g/m^3.
- **It is vacuum ultraviolet, not X-rays.** At 1 g/m^3, 40% of the heat leaves at 6-12 eV, 19% at
  12-16 eV, and 4% above 30 eV. Almost nothing reaches 100 eV.
- **Bag size does not help; density does.** At 1 g/m^3 the share falls only from ~89% to ~77%
  between a 0.1 kg and a 1000 kg projectile.
- **A lattice of small bags works only if each bag is burst first.** Intact 2.5 cm bags at
  10 atm fill 6e-5 of the volume, so 87% of the projectile's gas passes through the gaps, and the
  rest hits 16 kg/m^3 gas. That gives about 1% of the heat as light. Bursting a 0.5 m lattice
  about 1 ms ahead merges the clouds to within ±11%. That works, and it is easier than one big
  bag, which needs ~25 ms of lead and arrives densest at its core.

## Why it cannot be used in the pass-through pusher

The parent's pusher sets the pulse unit 5-10 m astern. The projectile threads the ship's clear
bore and the radiation case's entry hole, then strikes the target inside the case.

1. **A dilute projectile cannot thread the bore.** A 1 kg projectile at 1 g/m^3 is a cloud 12 m
   across. It cannot be launched compact and spread after the bore either. At 66 km/s it crosses
   the 5-10 m from bore to target in 76-152 us. Argon freely expanding into vacuum leaves at under
   1 km/s, so it grows less than 15 cm on the way.
2. **A compact projectile needs a target column hundreds of metres long.** A projectile that fits a
   1 m bore carries at least 1.3 kg/m^2. To sweep k = 2 of target it must cross 2.5 kg/m^2 of it
   (§`needle_through_fog`). Even with all of the heat landing in the target (the snowplow of the
   618 km/s design):

   | target density | target column needed | heat radiated |
   |---|---|---|
   | 1 kg/m^3 | 3 m | 21% |
   | 0.1 kg/m^3 | 25 m | 38% |
   | 10 g/m^3 | 255 m | 62% |
   | 1 g/m^3 | 2.5 km | 93% |

   A column a case can hold radiates about a fifth of the heat.
3. **The radiation case would have to enclose the cloud.** The 618 km/s target is 2 kg at
   10 kg/m^3, a column 1.6 m long behind a 40 cm disk. At 66 km/s the same 2 kg at 1 g/m^3 fills
   about 2000 m^3, a sphere 16 m across. That is about 10^4 times the volume and hundreds of
   times the wall area, for a case the parent already finds dominates specific impulse.
4. **And it pays back about 75 times less.** Per kilogram of projectile the best case gives
   ~1.3 GJ of vacuum ultraviolet, against 90 GJ of soft X-rays at 618 km/s.

| per kg of projectile, k = 2 | 618 km/s, dense disk into 10 kg/m^3 lattice | 66 km/s, argon cloud at 1 g/m^3 |
|---|---|---|
| ship-frame kinetic energy | 191 GJ | 2.18 GJ |
| heat | 127.3 GJ | 1.45 GJ |
| **light** | **90.1 GJ (47% of KE)** | **1.25-1.27 GJ (57-58% of KE)** |
| photon energy | kT 320-380 eV, 94% above 100 eV | kT ~3-4 eV, ~4% above 30 eV |
| target size | 1.6 m column, 40 cm across | 16 m cloud |
| fits through a bore | yes | no |

The share radiated is not the obstacle. At 66 km/s a gram-per-cubic-metre cloud radiates a
larger share of its kinetic energy than the 618 km/s lattice does. The obstacle is the cloud's
size. At 618 km/s, about 7.9 keV per carbon atom against about 0.2 keV per argon atom here, the
heat radiates from a target ten thousand times denser. That target fits through a bore and
inside a case.

## Suggested text for the parent

For §`why_not_xray_drive`, or a short closing paragraph of §`pass_through_pusher`:

> At the rocket's own 66 km/s collision speed, the pass-through pusher does not work. Each argon
> atom then carries only about 0.2 keV, so the fireball is at a few electronvolts. Its light is
> vacuum ultraviolet, and it escapes before expansion takes the heat only if the gas is very thin.
> Argon colliding at 1:2 radiates most of its heat only when spread to about
> \SI{1}{\gram\per\cubic\meter}, and an intact bag radiates under a tenth
> \cite{Katz_xray_collision_2026}. A kilogram at that density is a cloud \SI{12}{\meter}
> across. It cannot pass through the ship's bore or the case's entry hole, and it cannot spread
> after the bore, since it reaches the target within \SI{150}{\micro\second}. A projectile compact
> enough to pass would need hundreds of metres of equally dilute target to radiate half its heat.
> The scheme needs the energy per atom of the 4 solar radii dive, where a target ten thousand times
> denser still radiates.

## Not modelled; open for the parent to flag

- **Uniform clouds.** A real burst bag is densest at its centre. The lattice estimate is a hand
  calculation of overlapping Gaussian clouds.
- **One zone, a sphere, LTE.** As in the rest of the study. Emission is not the limit at these
  densities (a rough coronal estimate radiates the heat in ~0.1 us, against ~1 ms of expansion).
  Escape is the limit, and it is followed group by group.
- **Collisionality.** The ion mean free path is ~0.2 m at 1 g/m^3 against a 12 m cloud, which is
  fine. It becomes marginal near 0.01 g/m^3.
- **Bag skins.** They are condensed matter that radiates almost nothing at this speed. A 10 um
  film is ~20% of the argon's mass on a 2.5 cm bag at 10 atm.
- **Where the light goes.** Collecting vacuum ultraviolet from a 16 m cloud is not designed here.
  Since the pusher fails on geometry first, that design is not needed.
