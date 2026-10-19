# Kapton at 115-122 km/s, the bowl at 2.5 F, and the Jupiter chain's plates: answers (S44, S42)

Answers 2026-10-19 to the parent's S44 (its `docs/fast_growth_wave_asks_for_impact_sim.md`,
commit 5d1379c) and S42 (its `docs/plate_jupiter_loop_asks.md`). **Written to be copied into
the parent.** Evidence: Result 16 of [`xray_collision_study.md`](xray_collision_study.md) and
[ADR-0067](adr/0067-kapton-band-far-bowl-and-jupiter-chain-checks.md);
`make xray-collision-band-bowl` (`band_*.csv`). **Status: proposed.** S44 item 3's open items
(S36, S30 / S31) are not answered here; S46 is answered in
[`gas_spring_answers_for_parent.md`](gas_spring_answers_for_parent.md).

## Short answers

1. **S44.1: Kapton holds in the 100-130 km/s band wherever it holds at 99 km/s.** At the flown
   2 m standoff and u = 8, all three closures have a window: 84-109 um on HY-100 and 84-182 um
   on maraging 300 at 4 R_sun X 75. None does at u = 5, as at 99 km/s. **4 R_sun X 75 can fly
   ~100 um Kapton on a maraging face**, so the $83 row stands on the same footing as the
   9.86 R_sun headline.
2. **S44.2: the charge at 2.5 F costs 8-10% of `hole_ve`.** At 117 deg and k 10 it is
   15.4 -> 14.0 km/s (1423 s), mass ratio 1.381 -> 1.427. The hole's share falls (0.39 -> 0.34%)
   and no gas lands on it at either height. The rim's peak falls x0.4 and the bowl's peak x0.25.
3. **S42.1: the head-on bowl interpolation is ~2% optimistic.** Runs at 160-176 deg give 0-3%
   less than the calcs repo's `cos theta` interpolation; +2.1-2.3% at k 12.5-14. The hole takes
   ~0.5% of the impulse.
4. **S42.2: the push law overstates v_e by 9-17% at 54-69 km/s** (k 6.5). The runs cover k 6-8;
   speed enters only through the rebound and the bond toll, and the toll grows to 10-18% at
   these speeds. The sheets run 105-150 K cooler than at 99 km/s.

## 1. Kapton in the band (S44.1)

The planar sheet is re-run at each closure's own k on the flown 9.5 m plate: `c(k, 2 m)` =
0.625 / 0.664 / 0.653 at k = 7.45 / 6.20 / 6.53. The thermal floor is the thinnest Kapton whose
face stays under 400 C to the last 10 ms (Result 13's rule). The ceiling is the face's allowable
over `c w^2`.

| closure | closing, k | dose u = 5 / 8 / 12 | 100 um face, u = 8 | window u = 8: HY-100 / maraging | u = 12 | u = 5 |
|---|---|---|---|---|---|---|
| 9.86 R_sun, X 75 | 114.8, 7.45 | 114 / 45 / 20 kJ/m^2 | 634 K | 79-128 / 79-214 um | 29-128 / 29-214 | none |
| 4 R_sun, X 90, 5-day gap | 119.0, 6.20 | 111 / 44 / 20 | 632 K | 78-112 / 78-187 | 28-112 / 28-187 | none |
| 4 R_sun, X 75 | 121.8, 6.53 | 116 / 46 / 21 | 647 K | 84-109 / 84-182 | 29-109 / 29-182 | none |

- The 100 um face at 121.8 km/s peaks 26 K under the rating at u = 8 (647 K against 673 K). At
  99 km/s the margin was 100 K.
- The peak coefficient still grows on refinement (coarse 0.46-0.48 to production 0.63-0.66). If
  it rose another 37%, maraging's ceilings would still be 133-156 um and its windows would stay
  open. HY-100's would fall to 79-93 um and close at 4 R_sun. **Quote maraging 300.**

## 2. The bowl with its charge at 2.5 F (S44.2)

The point fireball is re-run on the δ/D 0.5 paraboloid with its charge at 2.5 F, priced as
Result 14 prices the 1.5 F bowl (rebound swapped in, Sun-plate mixing, carbon atomic).

| route, k | J/(m w), 1.5 -> 2.5 F | v_e | mass ratio | hole: r, impulse share |
|---|---|---|---|---|
| 4 R_sun 117 deg, k 10 | 1.903 -> 1.726 | 15.4 -> 14.0 km/s | 1.381 -> 1.427 | 4.1 -> 5.8 m, 0.39 -> 0.34% |
| 4 R_sun 117 deg, k 20 | 2.861 -> 2.619 | 11.3 -> 10.3 | 1.553 -> 1.618 | 0.39 -> 0.35% |
| 4 R_sun 117 deg, k 30 | 3.600 -> 3.307 | 9.0 -> 8.3 | 1.733 -> 1.819 | 0.39 -> 0.35% |
| 9.86 R_sun 129.9 deg, k 10 | 1.739 -> 1.557 | 10.7 -> 9.6 | 1.547 -> 1.627 | 3.2 -> 4.8 m, 0.54 -> 0.33% |
| 9.86 R_sun 129.9 deg, k 20 | 2.692 -> 2.442 | 7.7 -> 7.0 | 1.841 -> 1.960 | 0.58 -> 0.34% |

- **The law** loses 8-10% at every k, because the charge sits farther from the bowl's floor.
- **The hole** moves out and takes slightly less. No gas lands on it at either height.
- **The rim spray:** the peak on the hole's ring falls x0.34-0.43. Energy arriving per area there
  rises x1.8-1.9 at 117 deg and falls x0.3 at 129.9 deg. That ring is one coarse 3.5-6 m band, so
  this split is not resolved.
- **The bowl** peaks x0.24-0.28 lower anywhere on it.
- Use 2.5 F only where Result 15 needs it, at 117 deg. At 130 deg the bare impactor already leaves
  through the mouth at 1.5 F.

## 3. The bowl nearly head-on (S42.1)

The flown 1.5 F bowl is run at 160, 164, 168 and 176 deg. `J/(m w)` is listed with the calcs
repo's interpolation between 129.9 and 180 deg, and the hole's share:

| theta | k 7 | k 10 | k 12.5 | k 14 | hole (k 12.5): r, impulse, gas |
|---|---|---|---|---|---|
| 160 | 1.131 (interp +3.0%) | 1.494 (+0.1%) | 1.759 (+2.2%) | 1.906 (+2.3%) | 1.30 m, 0.52%, 0.03% |
| 168 | 1.101 (+3.1%) | 1.464 (0.0%) | 1.728 (+2.2%) | 1.874 (+2.3%) | 0.78 m, 0.51%, 0.05% |
| 176 | 1.087 (+3.1%) | 1.448 (0.0%) | 1.712 (+2.1%) | 1.859 (+2.3%) | 0.26 m, 0.51%, 0.04% |

At k 12.5, `hole_ve` is 6.0-6.2, 7.3-7.5 and 8.4-8.7 km/s at 66, 76 and 85 km/s closing. The
bowl keeps 59-65% of its `cos theta + sqrt(1 + k)` ceiling. The coarse grid agrees to 1%.

## 4. The push at 54-69 km/s (S42.2)

| closing | k 6.5: calcs -> re-priced v_e | bond toll | impactor / pulse | 100 um face u = 8 / u = 5 | floor u = 8 / u = 5 |
|---|---|---|---|---|---|
| 54 km/s | 21.3 -> 18.2 (x0.86) | 17.2% | 62 kg | 426 / 654 K | 29 / 91 um |
| 57 | 22.4 -> 19.6 (x0.87) | 15.3% | 58 | 433 / 679 | 31 / 104 |
| 61 | 24.0 -> 21.4 (x0.89) | 13.2% | 53 | 444 / 713 | 34 / 139 |
| 65 | 25.6 -> 23.1 (x0.90) | 11.5% | 49 | 456 / 748 | 36 / 266 |
| 69 | 27.2 -> 24.9 (x0.92) | 10.2% | 45 | 467 / 784 | 39 / 865 |

Over k 6-8 the ratio is x0.83-0.93.

- **The runs cover it.** k 6-8 lies inside Result 11's 5-40 runs, and a gamma-law run has no
  speed scale.
- **Speed enters only through the rebound and the bond toll.** The rebound is e = 0.68-0.72. The
  toll is the bonds' share of the heat, and the heat per kg falls as `w^2`. That is why the linear
  law drifts high as the stream slows.
- **The sheets are milder, as you expected.** Each pulse needs 1.6-2.2x the impactor, but the dose
  falls faster. At u = 5, 100 um holds at the slow end and fails above ~57 km/s.

## Suggested changes (for the parent to make)

- `sec:charge_pusher_plume`, the 100-130 km/s band: replace "nothing between has been run, so
  closures in that band are priced both ways" with: "Kapton holds at 115-122 km/s on the same
  terms as at 99 km/s: about 100 um on a maraging face at a 2 m standoff if the plume spreads at
  u >= 8 km/s, and no thickness at u = 5 \cite{Katz_puffsat_impact_sim_2026}." The 4 R_sun X 75
  row can then be printed on Kapton.
- The departure bowl at 117 deg: if the charge moves to 2.5 F, its v_e falls 8-10%
  (15.4 -> 14.0 km/s at k 10).
- `tab:no_sun_loop` and `sec:why_sun_dive`: the bowl near head-on is ~2% below the
  interpolation, and the push at 54-69 km/s is 9-17% below `push_ve`. Both make the Jupiter-only
  chain slightly worse, so the Sun dive's case only strengthens.
