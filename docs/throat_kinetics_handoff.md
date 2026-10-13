# S26 answered: recombination rates in the bulged chamber's nozzle

Handoff to `Balloon-Pulse-Propulsion`, answering its
`docs/throat_kinetics_asks_for_impact_sim.md` (register item S26). Written to be copied there.
The evidence is [`walled_nozzle_throat_kinetics.md`](walled_nozzle_throat_kinetics.md); the method
is [ADR-0057](adr/0057-finite-rate-nozzle-particle-count-and-carbon-bracket.md).

**Status: provisional.** Quasi-1D (one-dimensional, area-varying) steady flow through a straight
cone. Adiabatic and inviscid. Hydrogen recombines at its measured rate. Carbon has no measured
rate, so it is bracketed between two limits.

**Revised 2026-10-13 (second issue).** The first issue had no finite-rate row for the 2.66 m^2
throat at A/A* 100, which is the row the calcs ledger flies, and it tied the calcs ledger to the
0.669 ceiling. Since its ADR 0055 that ledger takes eta from `SURVIVABLE_IMPULSE`. Both are fixed
below, and the ledger has been re-run on the new figures.

## Conclusions for the paper

1. **The 906 / 971 kN s in `tab:one_chamber` are full equilibrium**, for carbon and hydrogen alike,
   and gross of the 1% gate. With the gate they are 893 / 958.
2. **Replace them with 885 / 941 kN s on the doubled 2.66 m^2 throat** (A/A* 100 / 300), the
   throat the calcs ledger flies since its ADR 0055. On the old 1.33 m^2 throat they are 883 / 938. These assume
   carbon trades atoms freely through two-body reactions. If carbon is instead frozen at the
   chamber state, the floor is 810 / 849 (808 / 845 on 1.33 m^2).
3. **These supersede the 861 / 892 (868 / 899 on 2.66 m^2) sent on 2026-10-12.** Those charged a sudden freeze of
   hydrogen recombination. That bound tripled the hydrogen loss, because the gas keeps
   recombining after the freeze point.
4. **The wider throat is safe on chemistry.** At the same area ratio it gains 0.2-0.3%. Option (a),
   2.66 m^2 on the old 22.5 m exit (A/A* 150), loses 3.5% because it expands less, not because
   of chemistry. Option (b), the 31.9 m exit, is the one to print.
5. **The equilibrium "bond energy recovered" overstates the finite-rate value.** At A/A* 300,
   0.51 of the bond store comes back in equilibrium, 0.44-0.45 at finite rate with free carbon
   exchange, and 0.31-0.32 with carbon frozen.
6. **Wall heat is not in these figures.** The 2-12% wall-heat debit of the 2026-10-12 handoff
   still applies on top.
7. **In the calcs ledger, A/A* 100 loses its lead.** The gain over 868 / 899 is uneven, +2% at
   A/A* 100 and +5% at 300, because the longer nozzle had more to lose to the sudden freeze.
   Re-run with these figures (24 kg pitch), doubling is 3.36 yr at A/A* 100 and 3.33 yr at 300,
   against ADR 0055's 3.57 / 3.86. The two ratios are now level. The carbon-frozen floor gives
   4.66 / 4.94 yr. Details: the calcs repo's `docs/companion_replies_2026-10-13.md`.

## The asks

**1. Which limit is 0.539?** Equilibrium in every species. The ask's "factor of 40 slower"
concern is real but small. The hydrogen freeze costs 1.3-1.5% of thrust at A/A* 300 with free
carbon exchange, and 3.2-3.5% with carbon frozen.

**2. The finite-rate expansion** (whole blowdown, 15° cone):

| nozzle | eta: equilibrium / finite / carbon frozen + finite / all frozen | bond returned: eq / finite / carbon frozen | net impulse [kN s]: eq / finite / carbon-frozen floor |
|---|---|---|---|
| 1.33 m^2, 22.5 m exit, A/A* 300 | 0.530 / 0.525 / 0.459 / 0.289 | 0.51 / 0.44 / 0.31 | 958 / **938** / 845 |
| (a) 2.66 m^2, 22.5 m exit, A/A* 150 | 0.499 / 0.497 / 0.443 / 0.289 | 0.47 / 0.44 / 0.32 | 919 / **908** / 827 |
| (b) 2.66 m^2, 31.9 m exit, A/A* 300 | 0.530 / 0.526 / 0.461 / 0.289 | 0.51 / 0.45 / 0.32 | 958 / **941** / 849 |
| 2.66 m^2, 18.4 m exit, A/A* 100 | 0.480 / 0.479 / 0.430 / 0.289 | 0.45 / 0.43 / 0.31 | 893 / **885** / 810 |
| 1.33 m^2, 13.0 m exit, A/A* 100 | 0.480 / 0.478 / 0.429 / 0.289 | 0.45 / 0.42 / 0.31 | 893 / **883** / 808 |

Here eta is the share of the rod's kinetic energy leaving as directed exhaust, over the whole
blowdown. The impulses debit the vacuum thrust, which counts exit pressure, rather than eta,
which does not.

**The whole blowdown matters.** At the peak chamber state alone the hydrogen freeze costs only
0.15-0.2% of thrust. Later parcels leave a thinner chamber and freeze sooner, so the whole pulse
loses 7-9 times as much.

**3. The sphere** (20 m^3, 3.4 kg/m^3, 496 bar, A/A* 300) loses 0.6% of thrust to the hydrogen
freeze, against 1.5% in the bulged chamber. Its gas falls out of equilibrium at A/A* ~100, against
~39. So density is the cause. The penalty stays small because most of the hydrogen store has
returned before either nozzle freezes.

**4. The wall angle.** The contour is a straight 15° cone from the throat, 39.6 m long on today's
nozzle and 56 m on option (b). A 10° cone is 52% longer and recovers only 0.1-0.4% of thrust.
**It is not worth a row.** The rate's own uncertainty (half a decade) moves the answer 2-4 times
as much.

## Carbon

No one has measured the rate of carbon atoms joining with hydrogen in a three-body collision. Only
one measurement exists for two carbon atoms joining (Slack 1976, argon bath). Neither limit can be
promoted to a result. The equilibrium carbon store is half released just past the throat, at
A/A* ~2.5. Keeping up there needs a coefficient about 100 times smaller than Slack's measurement.
For the carbon-frozen floor to hold, carbon would have to join about 100 times more slowly than the
one carbon reaction anyone has measured. This favours the free-exchange figures in item 2 but does
not prove them.

## What changes in the paper

- `tab:one_chamber`: on the doubled 2.66 m^2 throat, 885 / 941 kN s (A/A* 100 / 300), with the
  carbon-frozen floor (810 / 849) in a note. Option (a) is 908.
- `tab:nozzle_area_ratio`: the "bond energy recovered" column is the equilibrium figure. Add the
  finite-rate value (0.44-0.45 at A/A* 300) beside it, or label the column as equilibrium.
- Wherever the paper prints ADR 0055's doubling times (3.57 / 3.86 yr) or says A/A* 100 leads,
  it should take the re-run ledger once the calcs repo adopts these impulses: 3.36 / 3.33 yr,
  the two ratios level.
