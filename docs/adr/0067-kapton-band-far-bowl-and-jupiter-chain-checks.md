# ADR-0067: Kapton holds at 115-122 km/s where it holds at 99, and the 2.5 F bowl costs ~9%

Status: proposed, 2026-10-19. Result 16 of [`xray_collision_study.md`](../xray_collision_study.md);
code `sweep --band-bowl-2d`, `band_bowl.py` (`make xray-collision-band-bowl`). Answers the parent's
S44 items 1-2 (its `docs/fast_growth_wave_asks_for_impact_sim.md`, commit 5d1379c) and S42 (its
`docs/plate_jupiter_loop_asks.md`).

## Context

Result 15 (ADR-0066) ran Kapton at 99 and 132-208 km/s. Three printed closures sit between, and
one of them (4 R_sun, X 75) would move the parent's headline if it flew Kapton. Result 15 also
moved the departure charge to 2.5 F so that a bare impactor clears the bowl's rim, but `hole_ve`
had been run only at 1.5 F. The calcs repo's Jupiter-only chain (its ADR 0075) leaned on two
extrapolations: the bowl at 160-180 deg, interpolated between runs at 129.9 and 180 deg, and the
push at 54-69 km/s, taken as the 100 km/s law scaled by `w`.

## Decision

1. **Kapton holds in the 100-130 km/s band on the same terms as at 99 km/s.** At the flown 2 m
   standoff, with each closure's own k, the window at u = 8 is 79-128 / 78-112 / 84-109 um on
   HY-100 and 79-214 / 78-187 / 84-182 um on maraging 300 (114.8 / 119.0 / 121.8 km/s). There is
   none at u = 5, as at 99 km/s. Fly ~100 um on a maraging face. HY-100's window would close if
   the peak rose again on refinement as it did from the coarse grid (37%); maraging's would not.
2. **The bowl with its charge at 2.5 F costs 8-10% of `hole_ve`** at k = 10-30 and both visible
   routes (4 R_sun 117 deg at k 10: 15.4 -> 14.0 km/s, 1423 s; mass ratio 1.381 -> 1.427). The
   hole moves out (4.1 -> 5.8 m) and takes less impulse (0.39 -> 0.34%); the bowl's peak falls
   x0.24-0.28. Use 2.5 F only at 117 deg, where Result 15 needs it.
3. **The calcs repo's head-on bowl stands, ~2% high.** Its `cos theta` interpolation is 0-3%
   above runs at 160-176 deg (+2.1-2.3% at k 12.5-14). The hole takes ~0.5% of the impulse.
4. **The calcs repo's push law overstates v_e by 8-20% at 54-69 km/s** (x0.86-0.92 at k 6.5),
   because the bond toll grows to 10-18% as the heat per kg falls with `w^2`. The sheets run
   105-150 K cooler than at 99 km/s at u = 8.

## Consequences

- The parent may price 4 R_sun X 75 on Kapton (its $83 dear row), carrying the same u >= 8
  caveat as the 9.86 R_sun headline.
- The Jupiter-only chain gets slightly worse on both plates. That supports ADR 0075's verdict
  (the Sun dive is needed) and changes no decision.
- Open: the 2-D peak's convergence (unchanged), S36's u, the push plate's S30 / S31.
