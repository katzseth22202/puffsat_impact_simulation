# ADR-0058: Hot hydrogen is judged with a thick-gas wall, and the impact is re-solved for each lighter fill

Status: accepted, 2026-10-14. Answers `aim_is_all_you_need`'s
`docs/hydrocode_asks_2026-10-14_hot_hydrogen.md` (H1-H5, its ADR 0063); results in
[`walled_nozzle_hot_hydrogen.md`](../walled_nozzle_hot_hydrogen.md). Extends
[ADR-0056](0056-hydrocode-asks-blowdown-clock-and-fixed-charge.md) and
[ADR-0057](0057-finite-rate-nozzle-particle-count-and-carbon-bracket.md).

## Context

The calcs repo asks whether the bulged hydrogen chamber (106 m^3 at 2.5 kg, 0.84 m^2 throat,
37.0 t dry wrap) should be run at 7000 K or 10 000 K on a smaller charge. Its estimate scaled
today's efficiency by a chemistry ceiling and held the wall, the rate and the wall heat at
today's values. Three things about the hot charge break that:

- **its mass:** the polyethylene rod and plug are 22-27% of it and do not store hydrogen's
  energy, so pure-hydrogen `u(T)` misses the charge;
- **its radiation:** at 7000-10 000 K and ~0.25-0.3 kg/m^3 the gas is optically thick (TOPS),
  so the wall sees `sigma T^4`. Up to now the study's wall heat bracketed it between the thin H-
  floor and a blackbody;
- **its impact:** the spike beside the plug is a line blast nearly independent of the gas, but
  the cone's hit is not. Today's hydrogen fill already lets the remnant reach the cone at 1.19 GPa.

## Decision

1. **The charge is solved, not scaled:** `near_term.solve_charge` at the bulge's volume, rod and
   plug included. `solve_charge` takes the pulse energy, so the ends of the burn (74 / 81 km/s)
   are solved with k following w.
2. **Recombination and impulse are ADR-0057's march, unchanged**, through the AR 300 cone on the
   0.84 m^2 throat. Both carbon models are reported; the exchange model is the quoted row.
3. **Radiation is optically thick from the core:** `sigma T^4 (1 - exp(-tau_R))`, with `tau_R`
   the LANL TOPS gray Rosseland mean across the beam length, pulled for each charge's H:C
   (`data/tables/tops/tops_hydrogen_{7000,10000}_gray.html`). Below TOPS's 5802 K floor the
   in-band H- floor is used. Convection stays the study's Bartz bracket, quoted on an 800 K copper
   face (`bartz_liner`, which now also takes a station area ratio for the throat).
   `blowdown.integrate` accepts any `HeatFlux`.
4. **The impact is re-solved per fill.** `sweep --vessel-blast --hot-hydrogen` reruns the r 3.0 m
   bulge with the rod and plug as material and each lighter charge in it (0.5 cm). `zone_loads`
   scores it at the charge's own QSP.
5. **The wall's heat is judged against the charge's capacity as regen coolant** (77 K → 700 K,
   9.0 MJ/kg), not only against a flux limit. Thick radiation makes the wall heat set by the
   temperature, while the coolant scales with the charge, which shrinks.
6. **The cycle is `cycle.seals` on the lossy blowdown**, at both seal pressures and both fill
   rates. Nothing is chosen for the parent.

## Consequences

- 7000 K passes the wrap with more margin than today (1.77% fling, 29.8 t). 10 000 K fails it
  (3.48% in the taper; 2.24 GPa on the cone). The impact, not the gas state, rules out 10 000 K
  on today's bulge.
- Bonds returned are 60% / 55%, not 94%. Impulse per pulse is 340 / 285 kN s, and per kilogram
  carried it is flat from 7000 to 10 000 K.
- Today's own row moves: 683 kN s finite-rate against B5's 657 (Bray). The calcs repo's
  `HYDROGEN_BULGED_IMPULSE` is owed that correction.
- The parent's ledger, not this repo, decides whether 7000 K's +12% exhaust beats its ~29% thrust
  loss.

## Not done

- A 0.25 cm check of the 7000 K fill (as is owed for today's).
- The boundary layer's own absorption of the core's radiation, and the liner's base temperature.
- Any wall other than copper for the hot rows.
