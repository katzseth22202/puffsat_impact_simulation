# Hot hydrogen (asks H1-H5 of 2026-10-14): answers for the calcs repo

Handoff to `aim_is_all_you_need`, answering its `docs/hydrocode_asks_2026-10-14_hot_hydrogen.md`
(ADR 0063). It is written to be copied there. The evidence is
[`walled_nozzle_hot_hydrogen.md`](walled_nozzle_hot_hydrogen.md); the method is
[ADR-0058](adr/0058-hot-hydrogen-thick-gas-and-lighter-fills.md).

**Status: provisional.**
- The impact is a gamma-law 2-D rod-and-plug solve at 0.5 cm cells.
- Recombination is the finite-rate quasi-1D march of S26.
- Wall heat comes from the 0-D blowdown with TOPS-thick radiation and a Bartz bracket.
- Mechanism times (films 30 ms, lock 20 ms, fill 0.5 t/s) are B2's assumptions.

**This repo has not run the chain.** The rows below are physics inputs for `departure_levers`.

## Conclusions

1. **7000 K helps modestly. 10 000 K does not.**
   - 7000 K raises the effective exhaust 12% (10.3 → 11.5 km/s in your convention, impulse per kg
     of charge plus plug).
   - It halves the impulse per pulse (683 → 340 kN s), and the faster cycle does not make it up:
     thrust is ~29% below today's, like for like.
   - 10 000 K gains no exhaust speed over 7000 K, loses more thrust, and fails the dry wrap.
2. **The bonds come back at 60% (7000 K) and 55% (10 000 K)**, not 94%. On your ADR 0063 table
   that is between the 70% and 50% rows. Take eta from the impulses below
   (`efficiency_for_impulse`), not from the ceiling formula.
3. **The impact rules out 10 000 K on today's bulge.**
   - The 18.7 kg fill flings the wrap to 3.48% in the taper and hits the cone at 2.24 GPa.
   - The 24.1 kg (7000 K) fill passes at 1.77%, more margin than today's 2.20%, on a 29.8 t wrap
     (today 37.0 t).
4. **The wall heat is the open problem at 7000 K.** The gas is optically thick, so the wall takes
   3-5% of each pulse, mostly as radiation. That is 1.0-1.5x what the hydrogen can carry away as
   regen coolant (today 0.04-0.6x). The copper face rises 383 K per pulse.
5. **Correction for today's chamber:** its finite-rate impulse is **683 kN s**, not 657 (B5 used a
   sudden freeze). It is 661 with the high wall-heat edge. Update `HYDROGEN_BULGED_IMPULSE`.

## The rows for `departure_levers`

Bulged chamber, 106 m^3, 2.5 kg rod, 0.84 m^2 throat, AR 300, P = 1.76. Net impulse per pulse,
gate charged, head-on debit taken, finite-rate recombination, free carbon exchange.

| | k (75 km/s) | I at 74 / 75 / 81 km/s [kN s] | I with wall heat (75) | carbon-frozen floor (75) | bonds returned | wrap |
|---|---|---|---|---|---|---|
| today, 4884 K | 23.97 | - / 683 / - | 661 | 681 | 82% | 37.0 t, passes (2.20%) |
| **7000 K** | **9.65** | **328 / 340 / 416** | **327** | 336 | **60%** | **29.8 t, passes (1.77%)** |
| **10 000 K** | **7.47** | **275 / 285 / 345** | **264** | 281 | **55%** | **fails (3.48%)** |

k at 74 / 81 km/s: 9.34 / 11.64 (7000 K) and 7.25 / 8.90 (10 000 K). Your `slug_ratio_at` (constant
`u`) reproduces these to 1% at 7000 K and 3% at 10 000 K (9.17 against 8.90 at 81 km/s).

Intermediate temperatures, same convention (wall not solved in 2-D):

| T | k | I (75) | I with wall heat |
|---|---|---|---|
| 5500 K | 16.13 | 504 | 485 |
| 6000 K | 12.65 | 417 | 402 |
| 8000 K | 8.53 | 311 | 296 |
| 9000 K | 7.92 | 296 | 279 |

**Pulse rate** (blowdown with wall heat, films 30 ms, lock 20 ms):

| | sealed at 0.3 bar: 0.5 / 1 t/s fill | sealed at fill pressure: 0.5 / 1 t/s |
|---|---|---|
| today | 2.27 / 2.63 Hz | 2.85 / 3.43 Hz |
| 7000 K | **3.28** / 3.57 Hz | 3.77 / 4.14 Hz |
| 10 000 K | 3.61 / 3.87 Hz | 3.99 / 4.31 Hz |

- Today's 0.3 bar row is the 2.2-2.4 Hz you fly, so **the like-for-like hot rates are 3.3 Hz
  (7000 K) and 3.6 Hz (10 000 K)**.
- The mixed charge stays at ≤290 K everywhere, so a polymer film seals behind the hotter residue.
- The throat film per shot is 0.39 kg of aramid at 7000 K, against 0.56 kg today.

**Wall, if you carry it:** dry wrap 29.8 t at 7000 K, against 37.0 t. Tanks per kg of hydrogen are
unchanged. The extension is unchanged (same throat and area ratio).

## H3: the wall heat (per pulse)

| | wall share of pulse | of which radiation | heat per kg of charge (coolant capacity 9.0 MJ/kg) | copper face rise | burn-average flux at the 0.3 bar rate |
|---|---|---|---|---|---|
| today | 0.4-4.4% | ~0 | 0.4-5.1 MJ/kg | 202 K | 5.6 MW/m^2 |
| 7000 K | 3.0-4.7% | 1.5 MJ/m^2 of 1.7-2.6 | 9-14 MJ/kg | 383 K | 8.7 MW/m^2 |
| 10 000 K | 7.6-8.1% | 4.1 of 4.3-4.6 | 29-31 MJ/kg | 811 K (copper melts at 1358 K) | 16.6 MW/m^2 |

- **Optically thick:** Rosseland depth 7 / 65 across the chamber (LANL TOPS). The wall sees
  ~`sigma T^4` from the core: 136 / 567 MW/m^2 at peak.
- **Regen cooling cannot carry it at 7000 K without help** (extra coolant, a radiator, or a hotter
  outlet). At 10 000 K it is short 3.4x, and copper is not a liner.
- The throat's convection (0.1-1.7 GW/m^2 on the bracket) is the same order today and hot. It is
  open for every hydrogen row.

## H5: the sweet spot

| | I per kg carried (charge + tank, plug, film) | thrust at its 0.3 bar rate |
|---|---|---|
| today | 8.51 km/s | 1.50 MN |
| 6000 K | 9.30 | 1.18 |
| **7000 K** | **9.57** | **1.07** |
| 8000 K | 9.65 | 1.02 |
| 9000 K | 9.66 | 0.99 |
| 10 000 K | 9.63 | 0.95 |

Per kilogram carried, the curve is flat past 7000 K while thrust keeps falling, so 7000 K is the
hottest worth flying. Whether 6000-7000 K beats today depends on how the chain values exhaust
speed against thrust (finite-burn loss, burn time), which is your ledger's to run.

## Still open

- A 0.25 cm check of the 7000 K fill's wrap.
- A heat path for the 0-33% of wall heat the charge cannot absorb at 7000 K.
- The liner's base temperature and the boundary layer's own absorption.
