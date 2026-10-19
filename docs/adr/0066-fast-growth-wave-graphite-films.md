# ADR-0066: Above ~130 km/s the growth push needs dry carbon sheets, and a bare impactor misses the 130 deg bowl

Status: proposed, 2026-10-19. Result 15 of [`xray_collision_study.md`](../xray_collision_study.md);
code `sweep --fast-wave-2d`, `fast_wave.py` (`make xray-collision-fast-wave`). Answers the parent's
S40 (its `docs/fast_growth_wave_asks_for_impact_sim.md`, commit d9c7fa1).

## Context

ADR-0064 adopted the 2 Hz, 1 g push on bare ~100 um Kapton, checked at 100 km/s. The parent's
closures cross 400 km at 99-208 km/s, and the 4 R_sun headline at 137 km/s. The parent's own
estimate says Kapton has no window above ~130 km/s and asks whether graphite works, what it does
to the overtake law and the plate, and what a bare compact impactor does to the departure bowl.

## Decision

1. **The planar sheet's peak is `c(k, h) sigma w^2` at the closure's own k**: c(2 m) = 0.65 at
   k = 6.67-6.74, 0.61 at 7.93 and 0.60 at 8.21 (production grid, still a lower value). One run per
   k serves every speed, since a gamma-law gas has no speed scale.
2. **Kapton is kept only where its window is real.** At 99 km/s (the 9.86 R_sun headline) it is
   61-169 um on HY-100 at 2 m (u = 8). On maraging at 2 m the window is 141-165 um at 132 km/s and
   nil at 137, so it needs a 4 m standoff (143-188 um, 94% of the impulse) or HY-100 at
   8 m (84%). At 208 km/s no thickness up to 2 mm holds at u = 8. At u = 5 Kapton fails at every
   speed.
3. **Kapton is kept where it holds; dry carbon replaces it only where the speed rules it out**
   (Seth, 2026-10-19): the 4 R_sun headline and above, and fast departures that maximize the
   crossing speed. The carbon is flown as a planar sheet as Kapton is, best as a dry woven fabric
   (spread-tow ~64-90 g/m^2, triaxial or stitched against bias shear). The runs use graphite's
   properties; any all-carbon form heats the same per kg/m^2. They hold thermally at every spread run, including u = 5. They take the same
   pressure ceiling per kg/m^2 (26-64 um on HY-100, 44-106 um on maraging at 2 m), and sublimation
   is negligible below ~2600 K.
4. **Compact impactors, of any material, are not used on the push plate.** As a point blast they exceed HY-100
   below 6 m of standoff (2.1-9 GPa at 2.7-4 m), and at 6 m keep 80-83% of the sheet's impulse. The point-blast
   floor puts the plate at 35-49 t against the planar sheet's 12.2 t.
5. **The overtake law is re-priced at speed** for the push: the rebound and bond toll are taken at
   the gas's real arrival speed and heat. It sits 1-4% below the calcs repo's linear `push_ve` at
   99-137 km/s, mostly that law's interpolation in k, and 9% below at 208 km/s, where the hotter gas
   radiates (`e_eff` 0.68 -> 0.51). The impactor's material does not change it (< 0.1%).
6. **The 130 deg departure bowl needs no abort spray**: a bare impactor's line leaves through the
   mouth, 2.5 m inside the rim. At 117 deg (4 R_sun) it strikes the far wall near the rim and
   vaporizes a >= 0.9 m^2 patch (0.14-0.55% of the impulse); moving the charge to 2.5 F clears it.

## Consequences

- The 4 R_sun headline's push needs carbon, or Kapton at a 4 m standoff and maraging face if
  the plume's spread is u >= 8 km/s. Carbon also removes the push's dependence on S36's `u`.
- The calcs repo's `push_ve` should take `fast_wave_law.csv` above 130 km/s (or a factor of 0.97
  at 99-137 km/s and 0.91 at 208 km/s).
- Open: the peak's convergence; the carbon fabric's bias shear and seams at 6-9 m; the bowl's
  capture with the charge at 2.5 F; the vapour blast around a bare hit's hole.
