# ADR-0071: Iapetus water bags push on an ice charge in a graded paraboloid, spread by a gun's ice bullet

Status: proposed, 2026-10-21. Result 20 of [`xray_collision_study.md`](../xray_collision_study.md);
code `sweep --iapetus-2d [spread|bumper|bullet [key]] [coarse]`, `sweep --iapetus-rebound`,
`euler2d::bounce::run_bumper_charge_bounce`, `iapetus_bag.py` (`make xray-collision-iapetus`).
Answers the parent's S50 (its `docs/iapetus_bag_plate_asks_for_impact_sim.md`, commit 5b940a7);
answers doc [`iapetus_bag_plate_answers_for_parent.md`](../iapetus_bag_plate_answers_for_parent.md).

## Context

The parent's Iapetus cycle returns 25 kg bags of water that overtake an Earth-launched ship at
51.1 km/s (60 and 68 km/s as sensitivities) and push it to 10.916 km/s. It carried 10.3x per
kilogram of bag at k = 8.5 and `eta_jet` = 0.775 as a placeholder. A bag bursts from near a point,
unlike the floor-wide Kapton sheet of Results 13-19. Seth (2026-10-20) asked whether the charge
should be argon or water given water's bonds, whether the skirted plate is enough, and what the
plate weighs at 1 g and 2 Hz. In conversation he added the bag's miss, and (2026-10-21) chose a gun
that fires an ice bullet back along the bag's track, met ~30 m before the charge, to spread it.

## Decision

1. **The charge is water ice.** Water's bonds re-form in time (`H + OH + M` wins by Damkohler
   70-5000 in the fireball; the rebound stays above ~5500 K while it pushes). Its cost is
   thermodynamic: it expands like a gamma 1.17-1.27 gas, against argon's 1.30-1.44. That alone
   leaves the two level, ice slightly ahead on its higher rebound (`e_eff` 0.62-0.67 against
   0.45-0.65). What decides it is the face. Argon radiates 6-12% of the arriving energy into the
   wall against water's 0.7-4%, and the face must be replaced in kind. With a graphite face
   carried as charge, ice gives 7.1-8.2x and argon 6.4-7.3x on the paraboloid. The f(v) study's
   frozen-rebound bracket, about 0.9x off the ice columns, is carried as water's pessimistic edge.
2. **The plate is a delta/D 0.5 paraboloid with the charge at two focal lengths.** It gives
   7.7-8.4x for a realistic charge (compact column or spread disk) against 5.7-6.6x for the skirted
   flat plate, sends 79-90% of the exhaust's momentum along the axis, peaks at 0.25-2.2 GPa for a
   25 kg bag, and rocks a third as much per metre of miss. The 0.75 paraboloid adds ~0.5x but rings
   and rocks more.
3. **Carry ~8x, not 10.3x, at k ~ 4-8.5.** `eta_jet` is 0.5-0.6 for a real charge and falls with k
   for a compact column (0.53 at k = 8.5 to 0.32 at k = 20), which is a line source. The point
   fireball does not fall, so the parent's optimum near k = 20 holds only for perfect mixing.
4. **The plate is 25-32 t on the 8 m spring, graded to the bag's load.** One bag gives
   2.9-3.7 MN s, a 590-760 t ship at 1 g and 2 Hz. The centre takes 3-12x the mean impulse per
   area, so an even floor rings with 47-81% of each kick. Graded to the load (5-18 cm of steel at
   the centre, 30 kg/m^2 at the rim), it rings with 2-8%. The 63.4 t plate is twice the need.
5. **The face's inner ~6 m is a renewable graphite liner.** The face takes 0.3-50 MJ/m^2 a pulse
   against ~0.15 MJ/m^2 to melt steel. Ice consumes 1-9 kg of graphite a pulse (2-18 t a push),
   argon 13-49 kg at k = 8.5.
6. **A gun spreads each bag.** A 100 g ice bullet fired back along the track meets the bag ~30 m
   before a 2.5 m charge disk, where the merged cloud has spread to the disk's radius. It is the
   only spreader tested whose miss the plate survives (1.03 GPa on the production grid). Whether it
   shocks the bag's middle quarter or its whole face changes the impulse by 1-4%, so the solid-ice
   fracture question does not decide it. It costs leverage: 6.7x at k = 8.5 on ice against 7.2x
   for the compact charge and 8.35x for a perfectly self-spread bag, and 3.4x at k = 20, so it is
   flown at k <= 8.5.
7. **Not adopted:** a thin layer or mesh as a spreader (under 1% of the bag's column); a laser
   (surface only, MW-class at 2 Hz); rupture alone (~14 m/s, centimetre drops); a wide bumper 4-6 m
   ahead (its miss reaches 3.5-9.6 GPa); a 500 g bullet at 14 m (6 GPa miss); bag-borne H2/O2
   bubbles (Seth: too easy to set off); small spread disks of 0.8-1.6 m (4-21 GPa on the axis).

## Consequences

- The parent's S51 ledger should carry ~8x for a well-spread bag (7.1-8.2x with the face's
  graphite) and ~6.7x with the gun (k = 8.5), k ~ 4-8.5 of water ice, a 25-32 t graded paraboloid, a graphite liner consumed at 2-18 t a push, and the gun with
  ~200 kg of bullets a push.
- The ship carries water ice, not cryogenic argon.
- A bag the gun misses is a compact 33 GJ projectile. The plate is not a shield against it: a bare
  compact bag peaks at 29-40 GPa, holes the plate, and dumps ~25-30 GJ into the nitrogen spring.
- Rocking: grading the floor makes its rim kick 1.8-2.7x harder per unit of impulse offset than an
  even floor's, which the paraboloid's smaller lever (0.28 m per metre of miss) partly offsets.
- Open: the bullet at k below 8.5 and its face heat (the ring tally over-counts on its long runs);
  strength, fracture and atomization of a solid ice bag; mechanical erosion of the liner; the gun's
  design; the 1-D radiative share as the face heat.
