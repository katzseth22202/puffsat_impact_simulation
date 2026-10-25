# ADR-0073: A snow or water charge couples like solid where its column stays short; argon wins only at low speed

Status: proposed, 2026-10-25. Result 22 of [`xray_collision_study.md`](../xray_collision_study.md);
code `sweep --snow-charge [coarse] [sun|departure|push] | confined`, `sweep --iapetus-rebound
fast|slow`, the graded column of `euler2d::bounce::DiskCharge` (`ChargeLayer`,
`init_disk_charge_grid`), and `snow_charge.py` (`make xray-collision-snow-charge`). Answers the
parent's S60 (its `docs/snow_charge_asks_for_impact_sim.md`, commits dd4dcf8 and 860b570); answers
doc [`snow_charge_answers_for_parent.md`](../snow_charge_answers_for_parent.md). Builds on
ADR-0060 (the real collision), ADR-0065 (the departure bowl), ADR-0066 (the push sheet) and
ADR-0071 (water and argon races).

## Context

Seth (2026-10-24) asked what happens if the polyethylene charge on all three plates (the growth
push, the departure bowl, the Sun bowl) is water snow, froth, ice, liquid water or argon snow in a
thin bag. His reasons are debris (lost ice sublimes) and launch (water is dense and is made into
snow in flight). His premise is that a dense impactor couples well into a light charge whose pores
collapse. The parent expected a sweet spot: porous enough to heat, dense enough to stay a compact
source. In the session Seth added the 100 km/s Earth runs and the 3S cycle's 44-55 km/s Jupiter
return, and stated **k ≈ 10 at 44 km/s** as that cycle's design point.

## Decision

1. **Bulk density enters only as the column's length.** The gamma gas has no material, and the
   impactor couples only to the column in its footprint, so a lighter charge is a longer column.
   On the production grid (0.025 m on the bowls):
   - Sun bowl (k = 20, 0.25 m disk): within 2% of solid down to ~580 kg/m^3, within 5% to ~250.
     Snow at 300-600 costs 2-3%.
   - Departure bowl (k = 12.4, 0.4 m disk, on the axis): within 2% to ~590, 5% to ~430. Snow at
     600 costs 2%, froth at 500 3%, snow at 300 9%.
   - Push sheet (k = 6.4): the sheet's grid thickness from 0.1 to 0.8 m moves the impulse by 1.9%.
     The real 0.5-2.5 mm sheet's density does not matter.
   - Ice, liquid water and solid polyethylene are the same to within the grid (0-0.5%).
2. **Pores do not change coupling at the plates' own speeds.** In the snowplow the last kilogram
   is struck at `w / (1 + k)`, and its shock heat is exact. It is 5.5-57x the charge's vaporization
   energy (argon 76x or more), and 3.5x after the largest cold-compression share a linear-Hugoniot
   bracket allows (35%, polyethylene at the departure). Pore collapse can add at most that share
   to heat already in surplus. The 2-D gas already treats all dissipated energy as heat, which is
   the porous limit. A P-alpha kernel is therefore not needed for these plates. **At 44 km/s and
   k = 20 this fails for water** (its last kilogram gets 0.7x), which is one more reason for k ≈ 10
   there (2.6x).
3. **Material is a few percent at the plates' speeds, except argon at the Sun.** Water's and argon's
   bonds and ions do not freeze (no stranded share; Damköhler 1.5-1e8, the lowest on water's minor
   ion store at the departure). Over polyethylene, per kg of impactor, with polyethylene's carbon
   bracketed frozen / condensed:
   - push: water 1.01 / 0.98, argon 0.96 / 0.93;
   - departure: water 1.06 / 0.95, argon 1.02 / 0.91;
   - Sun: water 1.02 / 1.01, argon 0.85 / 0.84 (argon's rebound falls to 0.27 at 83 km/s: it
     radiates its stagnation heat into the face).
   **Water snow at 600 kg/m^3 is within -7% to +4% of solid polyethylene on every plate.**
4. **A graded charge does not raise impulse; it can cut the peak.** Ice at 917 under snow at 1/3 or
   1/5 of it, 25-75% dense: never above uniform solid (best -1% Sun, -3% departure, -1% push), and
   within about +/-2% of a uniform column of the same length. Dense at the entry is slightly better
   for impulse. **Light at the entry (the parent's rule) cuts the on-axis back-splash peak 4-10x on
   the departure and 2-5x on the Sun** for 1-6% impulse. Note the parent's orientation is inverted
   for the Sun bowl: there, as on the departure, the impactor enters from the plate side. Peaks are
   not converged (Result 10's 18-30% per grid halving), so this is a lever to study, not a load.
5. **The bag matters only on the push.** 0.05 kg/m^2 of film on both faces is 22% of the push's
   0.45 kg/m^2 charge (4% at 10 um), and it is polyethylene debris. On the bowls it is 0.05-0.2%.
   The gas a sealed bag holds is under 2e-4 of the charge.
6. **A sealed argon bag holds no 22% vapour.** Flashing liquid argon gives 18% vapour at the triple
   point and 25% on cooling to 60 K, but at 60 K the frost holds 904 Pa (0.072 kg/m^3), so 22% as
   vapour needs 3.0 m^3 per kg of charge. The bag is a dense charge only after it radiates the
   flash heat (42 kJ/kg): ~31 days for the Sun charge, ~64 for the departure's, shaded at 70 K.
7. **At low speed argon is the clean charge, and at k ≈ 10 water ties it.** At 44 km/s and k = 10
   (the 3S design point) water is 1.01x argon at a common gamma, 0.95x with the sheet's gamma and
   1.01x with the bowl's; polyethylene 0.71-0.73x (its carbon toll). Water's bonds stay broken
   through the push because stagnation re-heats the gas, not because they freeze: a softness, priced
   in its gamma and rebound. Argon radiates 10-100x more into the face at 6-12 km/s.
8. **Debris favours ice, but the push's film and any intact bag stay.** Sunlit at 1 AU an ice
   fragment of 1 cm lasts 1-3 days, a whole 0.4 m bag 49-120 days, snow at 300 a third as long;
   near the Sun, minutes to hours. A struck column vaporizes (decision 2). The film is 20 m^2/kg.

## Consequences

- aim can price each charge from `snow_charge_ledger.csv` (and `snow_charge_low_speed.csv` for
  44-100 km/s), and the bond toll as water's softness rather than the parent's frozen 19-30%.
- The departure's density cost is larger than the Sun's: its solid column is already as long as
  the impactor is wide. A wider impactor or a shorter charge (Result 21's lens) is the lever.
- Argon is not worth its cryostat at the Sun bowl (-15%) or the push (-4 to -7%), and is level at
  the departure; at 44-55 km/s it is the only charge with no softness, but water ties it at k = 10.

## Not modelled

- Strength, fracture and the bag's skin; how snow is made and held in flight.
- The departure's 130 deg entry (axisymmetric, on the axis, as Result 14).
- Peak convergence; solid argon on the Sun bowl (2.4 cells, held at solid's).
- A porous EOS: unnecessary at the plates' speeds by decision 2, and open at 44 km/s, k = 20.
- The low-speed panel's plate: it uses a 15 kg impactor at 3.75 m and both gamma brackets, since
  the 3S plate at 44 km/s is not named.
