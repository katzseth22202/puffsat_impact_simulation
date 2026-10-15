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
| carbon atomic | 0 | 12.0 | 24.1 kg | **12.0** | 92% | 0.88 g ($0.88) | 1.64 kg |
| carbon atomic | 0.16 | 12.0 | 24.1 kg | **10.1** | 90% | 1.06 g ($1.06) | 1.66 kg |
| bonds returned | 0 | 14.1 | 28.2 kg | 14.1 | 93% | 0.74 g ($0.74) | 1.61 kg |
| bonds returned | 0.16 | 14.1 | 28.2 kg | 11.9 | 92% | 0.89 g ($0.89) | 1.64 kg |

`k_S` = 51.7 in every row.

## Consequences

- **The expensive impactor costs about a dollar per kg of payload.** One kg of Sun-wave projectile
  sends ~52 kg of Earth impactor back, and that launches ~1250 kg at 11 km/s.
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
