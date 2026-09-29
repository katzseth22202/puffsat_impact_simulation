# The near-term wall: 1-D conduction through a coating on steel, driven by the solved pulse

Status: **accepted**, 2026-09-29. Follow-up to asks A4-A5 of
`todos/impact_sim_asks_methane_7000K.md` (paper `sec:methane_7000_near_term`). Extends
[ADR-0052](0052-near-term-chamber-feeds-janaf-and-carbon-branches.md).

## Context

The paper puts 0.2 mm of pitch over Cr-Mo steel for methane at 7000 K. It rejects graphite
because graphite "conducts too well". The design question raised in review was whether ~1 mm of
*sprayed* graphite over bare steel, cooled between pulses by spraying the next charge over it,
could replace both the pitch and any coolant plumbing. ADR-0052's chemistry answered half of it:
methane's boundary layer is carbon-saturated and deposits carbon below ~3700 K, so a carbon coat
is chemically stable. The thermal half needs the heat to be followed into the wall.

## Decision 1: one dimension, backward Euler, a node on the interface

`walled_nozzle/wall_layers.py` solves transient conduction through coating + steel:

- node-based finite volume, backward Euler, one Picard pass per step for the temperature-dependent
  properties;
- uniform nodes through the coating, geometric growth to 30 mm of steel, adiabatic back face;
- perfect contact at the interface. That is conservative for the steel: a real contact
  resistance would insulate it.

One dimension is sound because the heat penetrates millimetres in a 250 ms cycle, against a
1.68 m radius of curvature. It is verified against the semi-infinite convective (erfc) solution
to 1%.

## Decision 2: the gas side is coupled to the live surface temperature

Convection is Bartz (`near_term.bartz_liner`'s bracket), scaled over the pulse as `p^0.8`, times
`T_gas - T_s`. Radiation is grey exchange, `eps_w eps_g sigma (T_gas^4 - T_s^4)`: an opaque gas at
the high edge, the H- in-band floor at the low edge. The chamber history is `near_term.blowdown`.
Coupling to `T_s` matters for a coating that runs at 3000-4800 K: its own temperature cuts the
convective drive, and it re-radiates into an opaque gas. A fixed-flux boundary would overstate
the heat in.

## Decision 3: graphite's conductivity is bracketed, and the surface is capped at the carbon melt

Through-thickness conductivity varies 10x between a porous sprayed coat (~3 W/m/K) and dense
isotropic graphite (~40 W/m/K). That range *is* the answer, so it is bracketed (3, 10 and 40)
rather than picked. `c_p(T)` is the NASA-Glenn fit already used for the graphite ceiling. At
496 bar the chamber is above carbon's triple point (107 +- 2 atm, Haaland 1976), so the surface is
capped at the melt (~4800 K). Excess heat melts graphite at the heat of fusion (~9 MJ/kg), which
is the *low* removal heat, so the melted thickness is an upper edge. The paper's pitch model (k
0.2, surface held at 3900 K, 59.3 MJ/kg) is carried unchanged as the cross-check. Under this flux
it reproduces the paper's steel temperatures and pitch loss.

## Decision 4: the verdict turns on heat retained, not only on peak steel temperature

A coating can hold the steel cool through one pulse and still fail the cycle. The heat it keeps
has to be removed before the next pulse, and in the paper's scheme the next charge removes it.
The module therefore reports heat absorbed and heat carried off by removal. The result: porous
graphite (k ~ 3, 1-2 mm) protects the steel through a pulse, but it retains 110-243 MJ per pulse.
Pitch retains 28-50 MJ. Warming the whole 68.8 kg methane charge to 500 K absorbs ~100 MJ. **Pitch
stays; graphite alone would need coolant channels.**

## Consequences

- `make walled-nozzle-wall-layers` produces `data/results/walled_nozzle/near_term/wall_layers.csv`.
- The pitch **consumption** still rests on the paper's unverified 59.3 MJ/kg removal heat,
  borrowed from graphite sublimation, which is suppressed at this pressure. The retained heat
  does not depend on that number.
- Not modelled: the between-pulse spray quench (Leidenfrost), pitch pyrolysis-gas blowing,
  melt-layer runoff, and contact resistance.
