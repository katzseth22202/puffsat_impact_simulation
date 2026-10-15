# ADR-0056: The cycle's clock is the real-gas 0-D blowdown; volumes are compared at fixed charge

Status: accepted, 2026-10-12. Answers `aim_is_all_you_need`'s hydrocode asks B1-B6
([`walled_nozzle_hydrocode_asks.md`](../walled_nozzle_hydrocode_asks.md)). Extends
[ADR-0050](0050-walled-thermal-nozzle-methane-chamber.md) and
[ADR-0052](0052-near-term-chamber-feeds-janaf-and-carbon-branches.md).

## Context

The parent times its chamber cycle (ADR 0053 there) with a closed-form adiabatic blowdown and
asks this repository to confirm or replace it with transient solves. It also asks how the wall
loads, the cycle and the impulse change when the 212 m^3 chamber shrinks to 106 or 53 m^3 "at the
same rod, charge, temperature and throat". At a fixed pulse energy, the charge and the temperature
cannot both be held: the specific energy is the same, but a denser gas dissociates less.

The obvious transient solve is the 2-D chamber (`euler2d::vessel`) run through the whole
blowdown. Two things stop that:
- its ghost-cell wall is not conservative over tens of milliseconds;
- its transmissive throat edge feeds reversed lip flow back in.

## Decision

1. **The blowdown clock is 0-D on the real equilibrium gas** (`walled_nozzle.blowdown`). Its
   choked flux is the local `gamma_s` form, calibrated against the exact real-gas sonic flux on the
   isentrope (within 1%). The wall's heat loss is carried as `near_term`'s low/high bracket.
2. **The 2-D chamber gives only the first milliseconds** (the blast's head start through the
   throat). Its wall loses 3-7% of the mass by 20 ms in the rocket chambers and a third of it by
   30 ms in the sphere, and its throat outflow stalls.
   - `sweep --vessel-blowdown` records the leak and the outflow, so that a future conservative
     (cut-cell) wall can be checked against them.
   - Whole-blowdown discharge and sloshing are out of scope until that wall exists.
3. **A throat that discharges to vacuum uses `Bc::Outflow`:** transmissive while gas leaves,
   reflecting while it returns. Existing blast runs keep `Bc::Transmissive`, so they reproduce
   unchanged.
4. **Volumes are compared at fixed charge,** the 212 m^3 chamber's 108.4 kg of methane per 5 kg
   rod (54.22 kg per 2.5 kg). The temperature follows: 7,000 / 7,457 / 7,976 K at 212 / 106 /
   53 m^3. The parent's rows are fixed-charge too. `volume_trade.py`'s fixed-temperature trade
   stays as it was.
5. **The 5 kg chamber is the 2.5 kg chamber's state at half the volume.** Its clock, wall area and
   throat are the real ones. Wall histories are read at 5 kg by cube-root scaling and scaled to
   each chamber's own fixed-charge QSP (`wall_zones.History.qsp_real`).
6. **Net impulse per pulse carries the 1% gate and the H + H + M freeze,** with the head-on debit:
   `I = m_rod ((1 - gate)(1 + k + P) u_e - w)`. The chamber-blast note's 906 / 971 kN s had
   neither; 861 / 892 is the comparable figure.
7. **The cycle's gate is the throat film, not the fill pressure.** The residue the film traps
   joins the next charge, and the film must hold the mixed charge, bounded between the residue's
   thermal energy and all of it (`walled_nozzle.cycle`).

## Consequences

- The parent's blowdown estimate is confirmed at its fast end (221 / 281 ms against 215-260 /
  278-336). Its cycle gains a residue and film problem it did not have.
- The cushioned wall (parent ADR 0054) is rejected on the solved loads (`walled_nozzle.cushion`):
  the tiles bend at GPa.
- B3's unbulged chambers fail the dry wrap, so the bulge stays.

## Amendment, 2026-10-13: the throat doubles to 2.66 m^2

The 1.33 m^2 throat was sized to fit 2 Hz, not to recover chemistry. A wider throat recovers
more: a longer nozzle at the same area ratio gives H + H + M more time. Doubling it leaves
150-200 ms of the 2 Hz period, against 61 ms (or less than nothing with a polymer film sealed at
0.3 bar). The costs are 2.6 kg of throat film per shot and a ~32 m exit at area ratio 300; nozzle
size is not binding. The wall check passes: solved in 2-D with the rod and plug, the dry wrap is
46.2 t (against 52.4 t) with a 0.99 GPa hit on the shorter cone. The fully compact shape (159 m^3,
20° cone) fails the wrap at the taper, so the chamber keeps its 12° cone and only the throat
changes (`walled_nozzle_hydrocode_asks.md` §9).
The hydrogen chamber keeps its 0.84 m^2 throat: doubled, its dry wrap fails at the taper (2.69%
flung strain against 2.4%), and its cold residue already lets it cycle at ~2.4 Hz.

## Amendment, 2026-10-15: the wrap is sized for the liner's hoop swing too

The wall-zone score sized the dry wrap for strain only and never grew a crack in the 10 mm Cr-Mo
liner. Each pulse swings it from its autofrettage rest (-0.75 Y) up by `E * strain`. `wall_zones`
now integrates the shell through 50 ms of ring-down and rainflow-counts the liner's hoop stress.
It grows a 1 mm surface crack by Paris on the CC2938 hydrogen curve, and stiffens any zone whose
crack fails before 4,200 pulses. Every strain-sized wrap failed that test (621-1,340 pulses).
Sized for both, the hydrogen chamber's wrap is **43.8 t** (was 37.0) and the chosen methane
chamber's (throat x2) is **55.2 t** (was 46.2). The methane increment is owed only if hydrogen
reaches its bore: on the air curve its liners already pass. Flung strain improves slightly
(hydrogen 2.20% -> 2.16%). `walled_nozzle_hydrocode_asks.md` §9a.
