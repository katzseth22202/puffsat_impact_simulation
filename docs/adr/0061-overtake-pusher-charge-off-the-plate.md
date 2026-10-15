# ADR-0061: The Earth-departure pusher holds its charge off the plate and lets the merged mass drift onto it

Status: accepted, 2026-10-15. Result 11 of [`xray_collision_study.md`](../xray_collision_study.md);
handoff [`overtake_pusher_handoff.md`](../overtake_pusher_handoff.md). Builds on
[ADR-0060](0060-charge-pusher-point-fireball-2d-and-rebound-swap.md).

## Context

The departure wave at Earth gives a ship 11 km/s toward a high apoapsis. Impactors overtake it at
100 km/s, and the ship carries cheap mass ($1/kg against the impactor's $1000/kg). Seth asked for a
propellant fraction of 1/3 (mass ratio 1.5). He proposed a flat impactor meeting a flat target on
the plate, with the plume pressure-coupled across it, and asked whether k = 40 was reasonable.

## Decision

1. **The charge is held off the plate, not on it.** A solid impactor striking a target layer
   in contact with the plate at 100 km/s makes a shock of order `rho v^2 ~ 10 TPa`. The collision
   happens a few metres out, and the merged mass reaches the plate as a dense, dark gas.
2. **It is priced with Result 10's machinery, with the drift reversed.** The point fireball's drift
   is set toward the plate (`-1/sqrt(k)` of `u_rms`), with no bore. The real compact collision and
   a planar sheet-on-sheet collision check it. The rebound swap and the 2-D runs are unchanged.
3. **The bond toll is bracketed by carbon alone.** At 100 km/s the fireball cools below 2500 K
   before the plate, and `freeze.py` re-forms H2 with two to six decades of margin. Carbon has no
   rate to take, so it is bracketed between condensed and atomic.
4. **The recommended k is 10-12, not 40.** Mass ratio 1.5 is met at k ≈ 10.3-11.6 on the parent
   bowl and k ≈ 12.0-14.1 on the deep paraboloid. At k = 40 the gas has 122 MJ/kg, too little
   against polyethylene's 84 MJ/kg of bonds, and the mass ratio is 2.1-3.3.
5. **Planar sheets are chosen for loading, not impulse.** They add 3-6% to the impulse, but the
   plate's peak then scales with the impactor's mass per area (~0.05-0.07 kg/m^2 at 400 MPa). That
   lets a 7 m sheet carry ~11 kg per pulse.

## Consequences

- The impactor costs $1000/k per kg of propellant, ~$70-100, against ~$20 at the 400 km/s stage,
  because 100 km/s brings 16 times less energy per kg of impactor.
- The impactor sheet is ~35-50 µm of polyimide or carbon fibre. It is a solar sail, so it needs
  spin, attitude control and late guidance.
- Open: polyethylene's real EOS and soot in 2-D, the sheet collision's resolved peak, plate
  ablation at 20-45 km/s, and the charge sheet's deployment.
