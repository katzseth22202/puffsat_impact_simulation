# The sheets in the old plume: answers for `Balloon-Pulse-Propulsion` (S37, partly)

Answers 2026-10-18 to the parent's `docs/plume_cadence_asks_for_impact_sim.md` (S36-S37, commit
5bde6a0). **Written to be copied into the parent.** Evidence: Result 13 of
[`xray_collision_study.md`](xray_collision_study.md) and
[ADR-0064](adr/0064-cadence-plate-inertia-bound-and-sheet-thickness.md);
`make xray-collision-cadence-plate` (`cadence_plate.csv`, `cadence_sheet.csv`, `coated_sheet*.csv`).

**Adopted 2026-10-18 (Seth): the push flies at 2 Hz and 1 g mean**, 7.35 MN s per pulse, with bare
~100 um Kapton sheets (28.4 kg per pulse, an ~8 m sheet) on a 12.2 t flat plate with a 2 m skirt
and a solid HY-100 / 4340 steel face. No coating.

**S36 (the backward plume's spread) is not answered yet.** Everything below uses your
expanding-cloud model with `u` as a free input, and `u` decides the cadence.

## S37.3: the 1 mm sheet is rejected

- The planar sheet's peak is `c(h) sigma v^2`, independent of the sheet's radius. From 2-D runs at
  k = 6.74: c = 0.65 / 0.50 / 0.35 at 2 / 4 / 8 m standoff (a lower value; it grows ~37% from the
  coarse grid to the production one).
- A 1 mm sheet (1.42 kg/m^2) puts 9.2 / 7.0 / 5.0 GPa on the face, above every steel allowable.
  0.5 mm gives 4.7 / 3.6 / 2.5 GPa. The usable ceiling is ~140 um on HY-100 / 4340 (1.5 GPa)
  and ~250 um on maraging 300 (2.5 GPa) at a 2 m standoff.
- A small flat plate keeps its impulse only with a skirt: sheet radius + 1.5 m with a 2 m skirt
  keeps 94-97% of the 20 m skirted plate at 2 m standoff, and 72-82% bare.

## S34 inputs: the plate is set by its absorber

On the calcs repo's 100 g cap and 2 m stroke (`m_p >= I^2 / (2 s M a)`, end-of-push mass):

| push | plate |
|---|---|
| 4 Hz, 12 MN s | 32.4 t |
| 1 Hz, 14.7 MN s (1 g) | 48.8 t (24.4 t at 4 m, 9.8 t at 10 m) |
| 2 Hz, 7.35 MN s (1 g) | 12.2 t |

A planar sheet loads the plate evenly, so the point-blast survival floor (30 kg/m^2 scaled as
energy^(2/3), your 41 t and 47 t) does not apply. Its size only sets its thickness.

## S37.1-2: the sheet's heating, at 1 and 2 Hz

Your model, re-implemented, gives your doses (633/251/113 kJ/m^2 at 4 Hz; 47/19/8 at 1 Hz). At
2 Hz and 1 g it gives **96 / 38 / 17 kJ/m^2**, twice the 1 Hz dose, peaking at 375 kW/m^2 ~96 ms
before impact (u = 8). The heat was conducted through the film, still at full absorption:

- **1 Hz, 1 g:** 100 um holds at every spread (face 619 K at u = 5).
- **2 Hz, 1 g:** 100 um holds at u = 8 (face 572 K) and u = 12 (406 K). At u = 5 every film up to
  250 um passes 400 C ~100 ms out; 100 um decomposes through, 150-250 um lose a 31-38 um face.
- **Coatings do not help.** The heat is ~50 eV/amu particles stopping in the top ~0.1 um, not
  light, so a metal reflector under an emissive face sits at the face's temperature. Ideal
  emitters on both faces lower the peak 10-50 K. A 25 um shield film ahead fails before the flux
  peaks.

So 2 Hz at 1 g, with a quarter of the 1 Hz plate and the same lob, is workable with ~100 um sheets
if u >= 8 km/s. At u ~ 5 it depends on S37.2: whether a sheet with a decomposed 30-40 um face
still works as the impactor. Accommodation and sputtering (S37.1) and the strands (S37.4) are not
yet done.
