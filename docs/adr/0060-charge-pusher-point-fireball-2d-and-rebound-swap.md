# ADR-0060: The charge pusher is priced from a 2-D fireball, with polyethylene's rebound swapped in and the bonds charged

Status: accepted, 2026-10-15. Result 10 of [`xray_collision_study.md`](../xray_collision_study.md);
handoff [`charge_pusher_handoff.md`](../charge_pusher_handoff.md). Builds on
[ADR-0059](0059-xray-pusher-400kms-light-foam-on-dense-anvil.md) (the rebound table and
pulse-split ledger) and on the f(v) study's 2-D/1-D split (ADR-0003).

## Context

Seth wanted a dumber, lower-risk near-Sun pusher. Projectiles cost money and ship-carried mass is
cheap. A graphite disk at 400 km/s threads the bore and detonates a cheap polyethylene **charge**
astern, and the fireball pushes the parent's plate directly, with no radiation case. Several
modelling choices could each move the answer by tens of percent, and some are not obvious.

## Decision

1. **The figure of merit is projectile kg per kg of final mass, subject to a mass ratio of 3 or
   less** (Seth's objective), not Isp. That is why high k (60-100) is recommended despite lower Isp.
2. **The charge is near-solid polyethylene in the disk's path, not foam.** A projectile couples
   only to the column it crosses. A foam column long enough to stop it (5 m at k = 30) is a line
   source that loses 15-39% of the impulse in 2-D. The term is **Charge**, cross-referenced to
   the spray plate's **Plug** (CONTEXT.md).
3. **The plate impulse comes from lossless effective-gamma 2-D Euler** (`euler2d`, gamma 1.4, as
   the plug study; the CH2 table gives 1.3-1.5), using a **point fireball**. This is a hot
   uniform sphere drifting away at `v/(1+k)`, which the real disk-into-charge collision checks.
   The collision gives -4% to +25% more impulse, so the point fireball is conservative. It is
   kept as the production model because it resolves cleanly and is scale-free.
4. **Polyethylene's radiative rebound is swapped in** as `(1 + e_eff)/(1 + e_gamma)`, with
   `e_eff` from the 1-D radiation-hydrodynamic table at the plate speed and arriving column, as
   in the f(v) study. The open bore's impulse is tallied separately and subtracted.
5. **The bond toll is charged in full.** The CH2 table counts atoms only. `freeze.py` shows the
   gas reaches the plate too hot for H2 and that ionization recombines. So the 84 MJ/kg of bonds,
   but not the ionization, is a real loss: 1.6% of velocity at k = 30, 5.4% at k = 100.
6. **Peak pressure is quoted as an upper value for the point fireball, and as a hazard for the
   real collision.** The disk-into-charge back-splash puts a narrow jet of several GPa on the
   ring around the bore. It grows on refinement, so no converged value is claimed, and a
   sacrificial collar is required.

## Consequences

- The parent's unchanged 20 m bowl meets the mission (mass ratio 2.0-2.9). A δ/D 0.5 paraboloid
  is offered as an upgrade, justified by enclosure rather than focus.
- Plate survival away from the bore follows the point fireball (~0.1-0.3 GPa per kg of
  projectile). The bore collar is the open structural item. A charge shaped to tamp its own
  splash is the alternative to study.
- Shocked solid polyethylene lies above the CH2 table's 1284 kg/m^3 ceiling. The one-zone start
  is held there, which only understates the fireball's opacity.

## Amendment (2026-10-15): the collar is sized by impulse, and its erosion is charged as propellant

`collar.py`, `collar_ledger.csv`; the study's "bore collar" section.

- **The collar is sized by impulse, not peak pressure.** The back-splash lasts microseconds, so
  the ring is damaged by its impulse, and impulse converges on refinement where the peak does not.
  The collar's areal mass is impulse-matched: the plate's 30 kg/m^2 scaled by the ring's impulse
  over a reference. The jet's increment (reference: the point fireball on the same ring) gives
  ~1-1.5 t. Grading for the smooth load as well (reference: the plate's mean) gives 3.3-3.8 t.
- **Erosion is charged as propellant, and credited with only its own push.** That push is the 1-D
  ablating slug's extra impulse per area. It is not the 1-D gain applied to the ring's 2-D impulse,
  which is mostly pressure from stagnated gas and would give the vapour an impossible exhaust
  speed. The f(v) study's vapour curtain is off here: it makes energy when the ablated mass
  exceeds the arriving column.
- **Consequence.** At the floor the vapour is better reaction mass than the charge (45-72 against
  20-29 km/s). At the hot-jet ceiling the mass ratio rises to at most 2.93. The open structural
  item is now the inner ring at k = 30, which may have to be a fed sleeve, and grading the plate.
