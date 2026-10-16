# ADR-0062: The Sun wave redirects mass; the Sun ship's whole remainder becomes the Earth wave

Status: accepted, 2026-10-15. Code: `mission.breeding_cycle` (`make xray-collision-mission`,
`mission_breeding.csv`). Builds on [ADR-0060](0060-charge-pusher-point-fireball-2d-and-rebound-swap.md)
(the Sun stage at mass ratio 2) and [ADR-0061](0061-overtake-pusher-charge-off-the-plate.md) (the
Earth departure at mass ratio 1.5). Terms: **Sun wave**, **Sun ship**, **Earth wave** and
**Breeding factor** in `CONTEXT.md`.

## Context

The two charge pushers were priced separately: the Sun boost (400 km/s closing, 20.8 km/s to a
100 km/s Earth crossing) and the Earth departure (100 km/s overtake, 11 km/s). Seth closed the loop.
The ships boosted at the Sun exist only to become the Earth wave's impactors. The expensive Sun-wave
projectile is the "first impactor" that aims the rest.

## Decision

1. **A Sun ship's whole remaining mass is Earth-wave impactor.** It releases it as Kapton sheets
   with guidance before the Earth crossing. Anything it cannot release (its plate, if kept whole)
   is dead mass `f_dead`, bracketed between 0 and the plate's share. The δ/D 0.5 paraboloid at
   30 kg/m^2 is 0.16 of a 100 t final mass.
2. **The cycle is priced per kg of Earth impactor.** With `k_E` from Result 11 and `k_S` from
   `mission.k_for_mass_ratio`, both on the δ/D 0.5 paraboloid:
   - `D = k_E/(MR_E - 1)` departs;
   - a Sun ship of `S = MR_S/(1 - f_dead)` returns one kg of impactor;
   - the breeding factor is `G = D/S`.

   The payload share of each departure is `1 - 1/G`.

| Earth stage's bonds | `f_dead` | `k_E` | departs per kg of impactor | `G` | payload share | Sun-wave impactor per kg of payload | in Earth orbit per kg of payload |
|---|---|---|---|---|---|---|---|
| carbon atomic | 0 | 10.6 | 21.2 kg | **10.6** | 91% | 1.86 g ($1.86) | 1.66 kg |
| carbon atomic | 0.16 | 10.6 | 21.2 kg | **8.9** | 89% | 2.26 g ($2.26) | 1.69 kg |
| bonds returned | 0 | 12.0 | 24.0 kg | 12.0 | 92% | 1.63 g ($1.63) | 1.64 kg |
| bonds returned | 0.16 | 12.0 | 24.0 kg | 10.1 | 90% | 1.97 g ($1.97) | 1.67 kg |

`k_S` = 27.9 in every row. *Corrected 2026-10-15:* both stages are now priced on the real
collision, not the point fireball (ADR-0060 amendment). The first version had `G` 10.1-14.1 at
`k_S` = 51.7 and ~$0.74-1.06 per kg.

## Consequences

- **The expensive impactor costs about two dollars per kg of payload.** One kg of Sun-wave
  projectile sends ~28 kg of Earth impactor back, and that launches ~600 kg at 11 km/s.
- **The loop multiplies impactor, not total mass.** Every departing kg must first be in Earth orbit
  with its polyethylene: ~1.64 kg in orbit per kg of payload departing. The chemical equivalent for
  11 km/s is a mass ratio of ~12.
- **It costs less per kg of impactor than the Jupiter-return scenario appears to, but buys less
  mass per kg.** The calcs repo's departure wave meets argon PuffSats at 45.6-65 km/s at k ≈ 8.5 for
  a 5.3-7.2 km/s departure, with no Sun stage to pay. Here the 100 km/s impactor exists only
  because a Sun ship burned as much polyethylene as the impactor it returned (mass ratio 2), and
  the departure is the harder 11 km/s. That is why the departure mass per kg of impactor looks
  similar.
- **Open:** releasing (or cutting up) the Sun ship's plate as impactor; the generation time
  (~5 years: Earth to Jupiter's distance, then the fall); the Earth wave's guidance mass inside
  `f_dead`.

## Amendment, 2026-10-16: the kept plate is 1-2.5% of the ship, not 0.16

The calcs repo sizes the Sun ship's paraboloid as pulse energy^(2/3) at a fixed peak pressure, with
pulses capped at 2 Hz over a 2-5 h arc. That plate is 6-19 t, **1-2.5% of the ship** (its ADR 0072,
`aim_is_all_you_need` 127162f), and the parent adopts it (G5 in
`docs/solar_environment_answers_from_parent.md`): the plate is kept as dead mass, and cutting it
into sheets stays an option not taken. This repository's 0.16 was the geometric floor of a fixed
δ/D 0.5 paraboloid on a 100 t final mass, without the arc's pulse cap.

`mission.PLATE_DEAD_CALCS` carries the bracket as a cited input, and `mission.dead_fractions()`
now runs released (0), kept at 2.5%, and the old 0.16 as a pessimistic row
(`make xray-collision-mission` -> `mission_breeding.csv`):

| Earth stage's bonds | `f_dead` | `G` | payload share | Sun-wave impactor per kg of payload |
|---|---|---|---|---|
| carbon atomic | 0.025 | **10.4** | 90% | 1.91 g ($1.91) |
| bonds returned | 0.025 | 11.7 | 91% | 1.67 g ($1.67) |

So the growth is ~10-12x a generation, not ~9-12x. The plate item under **Open** is closed.
