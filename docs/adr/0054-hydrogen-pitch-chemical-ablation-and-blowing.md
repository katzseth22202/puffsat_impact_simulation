# ADR-0054: Chemical ablation and blowing in the coated-wall solve (N18 item 1)

Status: accepted, 2026-10-05. Extends [ADR-0053](0053-coated-wall-conduction-under-the-near-term-pulse.md).

## Context

`aim_is_all_you_need` @ `acfe859` revised N18: the chemical loss of a pitch coat in hot hydrogen is
13-26 um per pulse, not 5-90, and the unsolved term is the pitch's *thermal* recession under
hydrogen. Item 1 asks for the A5 wall solver run on `hydrogen_5500` with the chemical term
coupled in. The agreed order is bare Cr-Mo steel first: the repository holds no GRCop-84 or A-286
properties, and inventing them would break the provenance rule.

## Decision

1. **Mass flux and heat flux share one blowing parameter.**
   `mdot = g ln(1 + B')` and `q_conv = [ln(1 + B')/B'] h (T_gas - T_s)`, with
   `g = h_g / c_p` taken on the *same* Bartz edge as `h` (frozen `c_p` low, equilibrium high), so
   the `c_p` bracket cancels. `near_term.liner_heat_capacities` is the one place that builds them.
2. **B' is the companion's full-species table**, copied into `blowing.py` with its provenance.
   `near_term.carbon_blowing_parameter` is not used: it lacks C4H2, C6H2 and C2H4 and reads
   `ln(1+B')` ~15% low at 3900 K.
3. **B' is lagged one Picard pass** (five passes when chemistry is on). It is steep in `T_s`, and
   the step-refinement test holds the result to 2%.
4. **Chemistry runs at every surface temperature**, unlike the thermal removal, which exists only
   once the surface is pinned at the cap. Chemical and thermal recession are reported separately.
5. `Coating.reacts` marks carbon surfaces; `GasHistory.g0 = 0` switches the chemistry off, which
   leaves the methane study's numbers unchanged.

## Not modelled, and what that does to the answer

- **Reaction enthalpy** is omitted by default. `chem_heat = CHEM_HEAT_C2H2` (9.43 MJ/kg C, the
  graphite -> C2H2 formation enthalpy at 298 K) is run as a bracket; it lowers the thermal term
  further and moves nothing else. It is not a solved enthalpy of the equilibrium mixture.
- **The grid does not move** as material is removed (as in ADR-0053).
- **B' below 2500 K** is an Arrhenius extrapolation of the two coolest table columns.
- **Real edge composition** (carbon already in the gas) is not applied; the table is an upper
  bound on B'.
- **No Lewis-number correction** (item 2).
- **Substrate** is Cr-Mo steel only.

## Consequences

- `make walled-nozzle-wall-layers-hydrogen` writes
  `data/results/walled_nozzle/near_term/wall_layers_hydrogen.csv`.
- Results and the reply to the companion: `docs/walled_nozzle_hydrogen_pitch.md`.
