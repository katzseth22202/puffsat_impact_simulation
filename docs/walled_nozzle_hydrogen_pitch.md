# Pitch under the hydrogen pulse: the coupled surface balance (N18 item 1)

Reply to `aim_is_all_you_need` @ `acfe859`, "N18 as it now stands", item 1. Scope and method:
[ADR-0054](adr/0054-hydrogen-pitch-chemical-ablation-and-blowing.md). Reproduce with
`make walled-nozzle-wall-layers-hydrogen`; tests in `test_walled_nozzle_blowing.py` and
`test_walled_nozzle_wall_layers.py`. This is the agreed first step: **0.2 mm pitch over bare
Cr-Mo steel**, `hydrogen_5500` (`near_term.PAIRINGS[0]`), one pulse from 300 K, both flux edges
and both throats. `g` is taken against a 3900 K wall.

## Result

| throat | edge | chemical (um) | thermal (um) | total (um) | hot band 3x A* (kg) | whole wall 35.6 m^2 (kg) | surface | steel peak |
|---|---|---|---|---|---|---|---|---|
| 0.149 | low | 2.0 | 0 | 2.0 | 0.001 | 0.09 | 2510 K, never at cap | 350 K |
| 0.149 | high | 27.3 | 10.9 | 38.1 | 0.022 | 1.76 | 3900 K for 74 ms | 420 K |
| 0.198 | low | 1.9 | 0 | 1.9 | 0.001 | 0.09 | 2475 K, never at cap | 345 K |
| 0.198 | high | 25.7 | 8.0 | 33.7 | 0.026 | 1.56 | 3900 K for 53 ms | 407 K |

- **The chemical term reproduces the companion's 26 um** at the high edge (25.7-27.3 um). `g0` comes
  out at 0.244-0.446 and 0.306-0.560 kg/m^2/s, matching its table, because it is built from the
  same `h_g / c_p`.
- **The thermal term is small: 8-11 um at the high edge, against methane's 94-120 um.** Two
  reasons. The hydrogen gas is cooler and far less radiative than the methane pulse. And the
  blown film cuts the convective heat to **0.35 of its unblown value** at the cap (`ln(1+B')/B'`
  with B' = 5.3), which more than halves the thermal term (18.7 um with chemistry off).
- **The low edge never reaches the cap.** The surface peaks near 2500 K, where B' is 0.2-0.7, so
  chemistry is 2 um, not the 13 um a 3900 K wall would give. The companion's low-edge estimate
  assumed a wall at 3900 K and is an overestimate.
- **Hot spot (2x load):** 64-70 um total, 0.041-0.049 kg in the hot band, 3.0-3.2 kg for the whole
  wall.
- **Reaction-enthalpy bracket** (graphite -> C2H2, 9.43 MJ/kg C): thermal falls by 3-4 um further;
  chemical moves under 2%.

## Against N18's thresholds

- **The hot-band text stands.** The high edge with a 2x hot spot is 0.05 kg, an order of magnitude
  under the 0.5 kg bar, and 70 um sits well inside a 0.3 mm spray. The 120 um methane thermal bound
  the paper borrowed is about ten times what this solve finds.
- **The whole-wall decision reopens.** Thermal recession is well under methane's 120 um, as N18
  said would reopen it. A whole-wall coat at the high edge costs 1.6-1.8 kg per pulse (2.6-2.9% of
  the 59.9 kg charge; 3.2 kg, 5.4% at a 2x load) against the paper's 4 kg. Whether it pays is now a
  question about the hotter-charge gain, not the pitch.
- Item 4 (strip threshold) is unblocked, but not done.

## Limits

Single pulse from 300 K, adiabatic. The substrate is Cr-Mo steel, not GRCop-84 or A-286, and at
350-420 K peak the steel is not what limits it; chasing the other substrates' properties is not
justified by these numbers. B' is the pure-hydrogen upper bound, with the cold tail extrapolated.
Item 2 (Lewis number) is not applied: if `St_m` sits tens of percent off `St_h`, the chemical term
moves by the same fraction.
