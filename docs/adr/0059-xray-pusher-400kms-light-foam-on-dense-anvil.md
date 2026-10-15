# ADR-0059: The 400 km/s X-ray pusher stops a light foam on a dense anvil, and prices the plate with a rebound

Status: accepted, 2026-10-15. Results 7-9 of
[`xray_collision_study.md`](../xray_collision_study.md); handoff
[`xray_pusher_400kms_handoff.md`](../xray_pusher_400kms_handoff.md). Extends the study's 618 km/s
work (Results 4-5, [`xray_pusher_handoff.md`](../xray_pusher_handoff.md)).

## Context

Seth asked whether the parent's pass-through X-ray pusher works at a Parker-like periapsis: a
200 km/s periapsis met head-on, closing at 400 km/s, against the parent's 618 km/s. At 400 km/s,
the parent's layout (a dense disk sweeping a light lattice) radiates half its kinetic energy only
from a lattice ten times thinner, so the column it must sweep becomes very long (Result 7).
Several alternatives were raised and settled in a grilling session.

## Decision

1. **Reverse the density contrast.** A light carbon-foam projectile (10 kg, ~10 kg/m^3, a
   Parker-style sandwich with C-C on the sunward side only and no front cap) stops on a dense
   **anvil**: 50 kg of graphite with a high-Z liner, at the far end of the case. In a collision
   the light partner takes the heat, so the foam stagnates at ~700 eV and half the projectile's
   kinetic energy radiates. Graphite rather than iron, because at a fixed areal mass a thick,
   low-density anvil outlasts Rayleigh-Taylor burn-through.
2. **Carbon, not steel, for the radiating mass.** Iron gets more energy per atom but shares it
   among ~20 free electrons and pays keV to free them. It starts cooler than carbon and depends
   more on opacity (Result 7).
3. **Reject a z-pinch or any magnetic conversion.** It adds two lossy stages (kinetic to
   electrical and back) ahead of the same radiating collision. Its self-generated field is about
   0.02 T. Spaced blanket foils are rejected too, because they explode within nanoseconds.
4. **Price the plate stage with the parent's own pulse-split optimizer** (`anvil.best_split`
   reproduces `tab:pass_through_isp` to 0.2%), with **equal filler and propellant masses**. A
   light-forward, heavy-backward stack loses 13% of the impulse and breaks the impingement cap.
5. **Compute the propellant's rebound and credit it.** It comes from the f(v) study's coupled
   bounce on a new polyethylene table: an atomic C/H Saha EOS plus a TOPS pull. It gives
   `e_eff` 0.63-0.72 at ~95 km/s, which the parent does not count. **Compute the case's delivery**
   from a time-dependent cavity balance in which the foam's own tail plugs the entry hole.

## Consequences

- The design reaches 2740-3300 s only with the rebound credit, and fails at 1680-1920 s without
  it. The rebound is therefore a quantity the parent's ledger should adopt or dispute.
- The entry-hole plug works only when the foam is longer than the case. A narrower foam (0.5 m)
  delivers 0.78 against 0.53, but carries heavier sheets and a longer pulse.
- The case-wall tamping mass (~100 kg against a ~15 kg budget) is the open risk. Result 10's
  charge pusher (ADR-0060) needs no case at all.
