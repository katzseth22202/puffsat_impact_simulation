# Walled nozzle: finite-rate recombination through the bulged chamber's nozzles (S26)

Status: analysis note, 2026-10-13. Answers the paper repo's S26
([`throat_kinetics_asks_for_impact_sim.md`](https://github.com/katzseth22202/Balloon-Pulse-Propulsion/blob/main/docs/throat_kinetics_asks_for_impact_sim.md)).
Method and limits are [ADR-0057](adr/0057-finite-rate-nozzle-particle-count-and-carbon-bracket.md).
Reproduce with `make walled-nozzle-throat-kinetics` (~10 min, 8 cores). Outputs are in
`data/results/walled_nozzle/near_term/finite_rate/` (`pulses.csv`, `bray.csv`,
`carbon_threshold.csv`). The summary for the paper repo is
[`throat_kinetics_handoff.md`](throat_kinetics_handoff.md).

## 1. Bottom line

| ask | answer |
|---|---|
| **1. Which limit is 971 kN s?** | **Full shifting equilibrium, carbon and hydrogen both**, gross of the 1% gate (958 with it). This repo said so on 2026-10-12 (`hydrocode_asks_handoff.md` item 4). Its 892 replaced it with a Bray sudden freeze of H. |
| **2. Finite rate** | The hydrogen freeze costs **1.3-1.5% of thrust** at A/A* 300 if carbon exchanges freely, and **3.2-3.5%** if carbon is frozen. That gives **938 / 941 kN s** on the 1.33 / 2.66 m^2 throats, with a carbon-frozen floor of **845 / 849**. B6's Bray debit (892 / 898) overstated the hydrogen loss threefold in impulse (66 kN s against 20) and tenfold in `eta`. |
| **3. Sphere** | It lands nearer equilibrium (0.6% thrust lost against 1.5%), so density is the cause. The effect is about 25 times smaller than the factor of 40 in recombination time suggested. |
| **4. Wall angle** | A straight 15° cone from the throat. 10° buys 0.1-0.4% of thrust for 52% more nozzle length. **It does not earn a row.** A rate wrong by half a decade moves the result 2-3 times as much. |
| **Carbon** | No rate exists, so carbon stays bracketed. The equilibrium carbon store is half released by A/A* ~2.5. Keeping up there needs a three-body coefficient of ~5e-47 m^6/s, **2.1 decades below** the one measurement (Slack's C + C + Ar). The carbon-frozen floor therefore needs carbon about 100 times slower than anything measured. |
| **Wider throat** | Kinetics do not oppose it. At A/A* 300 the 2.66 m^2 throat is +0.3% (941 against 938). Option (a), A/A* 150 on the 22.5 m exit, loses 3.5% to the shorter expansion, not to chemistry. |

## 2. Method in brief

- **One kinetic variable: particles per kilogram.** Bimolecular exchanges conserve it and three-body steps change it. The only sink is `H + H + M`, at the within-study efficiencies of Jacobs, Giedt & Cohen (`freeze.mixture_coefficient`). The reverse rate comes from detailed balance.
- **Carbon is bracketed by two composition laws:**
  - `frozen_carbon`: carbon frozen at the peak, and H alone moves. This is the floor.
  - `exchange_carbon`: constrained equilibrium on the particle count (Keck). Carbon re-forms molecules freely, but only from particles `H + H + M` has removed. This is the upper realistic case.
  Full equilibrium is the ceiling above both.
- **Quasi-1D march** from A/A* 1.05, stepping in `ln rho`, with the area from continuity and the clock from the cone. The gas is in equilibrium to there: the Damkoehler number at A/A* 1.05 is 900-1300 in the bulged chamber and 6000-7500 in the sphere.
- **Every parcel of the blowdown** is marched (14 parcels), not only the peak (§5).
- **Impulse is debited by the vacuum-thrust ratio** `u + p/(rho u)`, applied to B6's shifting impulse (gate included, adiabatic). `eta` counts kinetic energy only, as the study's convention does, and is reported too. The two disagree in one place. At a fixed area ratio a gas that freezes early keeps a higher gamma, so it can leave slightly *faster* than the shifting gas, but colder and at lower pressure. The thrust ratio is ordered in the rate; the speed is not (tested).
- **Tests:** a rate 1e6 times the measured one reproduces each branch's isentrope to 0.3%. A zero rate freezes the particle count exactly. Bond return and thrust rise monotonically with the rate. Step halving moves the result under 0.2% (`test_walled_nozzle_finite_rate.py`).

## 3. Ask 2: the three nozzles (whole blowdown, 15°, measured rate)

`eta` is on the blowdown (B6's convention). "Bond returned" is the mass-weighted share of each
parcel's chamber bond store returned by the exit; the chamber holds 87.9 MJ/kg of a 115 MJ/kg
charge. "Frozen" means carbon frozen at the peak and the particle count frozen at the throat.

| nozzle | `eta` eq | exch. + finite H | FC + shifting H | FC + finite H | frozen |
|---|---|---|---|---|---|
| 1.33 m^2, AR 100 | 0.480 | 0.478 | 0.435 | 0.429 | 0.289 |
| **1.33 m^2, AR 300** (22.5 m) | **0.530** | **0.525** | 0.478 | **0.459** | 0.289 |
| **(a) 2.66 m^2, AR 150** (22.5 m) | 0.499 | 0.497 | 0.452 | 0.443 | 0.289 |
| **(b) 2.66 m^2, AR 300** (31.9 m) | 0.530 | 0.526 | 0.478 | 0.461 | 0.289 |

| nozzle | bond returned: eq / exch. / FC / frozen | H freeze, thrust lost: exch. / FC | net impulse [kN s]: eq / exch. / FC / frozen | B6 Bray |
|---|---|---|---|---|
| 1.33 m^2, AR 100 | 0.45 / 0.42 / 0.31 / 0.06 | 0.8% / 2.0% | 893 / **883** / 808 / 597 | 861 |
| **1.33 m^2, AR 300** | 0.51 / 0.44 / 0.31 / 0.06 | 1.5% / 3.5% | 958 / **938** / 845 / 605 | 892 |
| **(a) 2.66 m^2, AR 150** | 0.47 / 0.44 / 0.32 / 0.06 | 0.8% / 2.2% | 919 / **908** / 827 / 600 | 880 |
| **(b) 2.66 m^2, AR 300** | 0.51 / 0.45 / 0.32 / 0.06 | 1.3% / 3.2% | 958 / **941** / 849 / 605 | 898 |

- **Carbon matters more than hydrogen's rate.** Freezing carbon costs 0.05 of `eta` at equilibrium hydrogen. The hydrogen rate then costs 0.005 (exchange) or 0.02 (frozen carbon).
- **Why B6's Bray figure was so pessimistic.** Bray freezes H at the first station with `Da < 1` (A/A* 39 on the 1.33 m^2 throat) and allows no recombination after it. The finite-rate gas keeps recombining past that point. With exchange carbon its exit `n_H^2/n_H2` is still ~9x the equilibrium value, but most of the store is back.
- **Frozen-carbon H freezes earlier** (A/A* 24-28 against 39-47). Carbon frozen as free atoms cannot take up hydrogen, so more H stays free to strand.
- **(a) against (b).** (a) loses the least to chemistry, but at A/A* 150 it expands less: 908 against 941.
- **Wall heat is not in these rows.** B6's bracket (2-12% of impulse at 212 m^3) still applies on top. The doubled throat empties sooner and loses less (B6: 840 against 789 at the high edge).

## 4. Ask 3: the sphere as calibration (A/A* 300)

| | Da at A/A* 1.05 | first `Da < 1` (exch.) | `eta` eq / exch. / FC | H freeze, thrust lost: exch. / FC | bond eq / exch. / FC |
|---|---|---|---|---|---|
| sphere, 0.198 m^2 (15.3 m cone) | 7,470 | A/A* 108 | 0.564 / 0.562 / 0.502 | **0.6%** / 1.6% | 0.55 / 0.51 / 0.34 |
| sphere, 0.149 m^2 (13.3 m cone) | 6,480 | A/A* 99 | 0.564 / 0.562 / 0.502 | 0.7% / 1.7% | 0.55 / 0.51 / 0.34 |
| bulged, 1.33 m^2 (39.6 m cone) | 933 | A/A* 39 | 0.530 / 0.525 / 0.459 | **1.5%** / 3.5% | 0.51 / 0.44 / 0.31 |

- **Density is the cause, as S26 suspected.** The sphere is 6.3x denser. Its throat Da is ~8x larger and its gas freezes at ~2.7x the area ratio. The hydrogen loss is 2.4x smaller.
- **The size is small in both.** Recombination is 40x slower at the same station, but the bulged nozzle is 2.6x longer and its throat Da is still ~900. Most of the H store returns before A/A* 40. Only the tail strands.

## 5. Why the whole blowdown, not just the peak

The ask names the peak state. At the peak, finite-rate H costs almost nothing:

| peak parcel only, AR 300 | thrust lost: exch. / FC | bond returned (exch.): shifting / finite | (FC): shifting / finite |
|---|---|---|---|
| 1.33 m^2 | 0.20% / 0.83% | 0.514 / 0.495 | 0.421 / 0.372 |
| 2.66 m^2 | 0.15% / 0.69% | 0.514 / 0.500 | 0.421 / 0.378 |
| sphere 0.198 m^2 | 0.05% / 0.33% | 0.558 / 0.550 | 0.432 / 0.399 |

Later parcels leave an emptier chamber and run through the same nozzle at lower density, so
they freeze sooner. The whole-pulse loss (1.3-1.5%) is 7-9x the peak's. A peak-only answer would
have understated it by that factor.

## 6. Ask 4: wall angle, and the rate's uncertainty

Thrust ratio to shifting equilibrium (bond returned in parentheses):

| nozzle | model | 10° | 15° | 20° | k / 3.16 | k x 3.16 |
|---|---|---|---|---|---|---|
| 1.33 m^2, AR 300 | exch. | 0.9876 (0.451) | 0.9852 (0.442) | 0.9833 (0.435) | 0.9769 (0.414) | 0.9910 (0.465) |
| | FC | 0.9688 (0.320) | 0.9647 (0.312) | 0.9613 (0.306) | 0.9506 (0.288) | 0.9753 (0.334) |
| (a) 2.66 m^2, AR 150 | exch. | 0.9931 (0.443) | 0.9917 (0.437) | 0.9904 (0.432) | 0.9861 (0.415) | 0.9951 (0.452) |
| | FC | 0.9810 (0.324) | 0.9780 (0.317) | 0.9755 (0.311) | 0.9676 (0.294) | 0.9854 (0.336) |
| (b) 2.66 m^2, AR 300 | exch. | 0.9893 (0.458) | 0.9872 (0.449) | 0.9854 (0.443) | 0.9797 (0.423) | 0.9923 (0.470) |
| | FC | 0.9720 (0.327) | 0.9681 (0.319) | 0.9651 (0.313) | 0.9551 (0.295) | 0.9780 (0.340) |

- **The contour is a straight 15° cone from the throat**, `x = r* (sqrt(A/A*) - 1) / tan(15°)`: 39.6 m (1.33 m^2, AR 300), 38.6 m (a), 56.0 m (b). It has no throat rounding and no bell.
- **10° lengthens the nozzle 52%** (60 m and 85 m at AR 300) and recovers **0.1-0.4%** of thrust, ~1-4 kN s. **Not worth a row.** A bell contour is no lever here either, since it opens faster than a cone.
- **Half a decade in the rate** (Baulch's stated uncertainty on `H + H + H2`) moves the result 0.3-0.8% (exch.) and 0.7-1.4% (FC), 2-4x the wall-angle lever. The rate uses Jacobs' `T^-1` fit, extrapolated above its 4,700 K window. Hurle's direct 2,500-7,000 K measurement is flat in temperature, which suggests the extrapolation errs slow (`rates.H_H2_HOT`).

## 7. Carbon: how fast it would have to be

Carbon has no evaluated rate (ADR-0050), so the floor and the upper case stay limits. This is the
Bray threshold on the equilibrium expansion: the three-body coefficient that keeps carbon up to
the station where half its free atoms have associated (`required_carbon_rate`). Partners are all
free atoms (C + H), and third bodies are all neutral heavies.

| nozzle | half-release at | T | free C: throat / half / exit | k needed [m^6/s] | Slack C + C + Ar | margin |
|---|---|---|---|---|---|---|
| 1.33 m^2, AR 300 | A/A* 2.6 | 4,936 K | 0.60 / 0.31 / 0.01 | 5.1e-47 | 6.1e-45 | 2.1 decades |
| (a) 2.66 m^2, AR 150 | A/A* 2.5 | 4,948 K | 0.60 / 0.31 / 0.02 | 3.5e-47 | 6.1e-45 | 2.2 |
| (b) 2.66 m^2, AR 300 | A/A* 2.6 | 4,936 K | 0.60 / 0.31 / 0.01 | 3.6e-47 | 6.1e-45 | 2.2 |
| sphere 0.198 m^2 | A/A* 2.4 | 5,182 K | 0.27 / 0.14 / 0.00 | 5.0e-48 | 5.7e-45 | 3.1 |

- **Carbon releases its store just past the throat**, at A/A* ~2.5 where the gas is still dense. Keeping up there takes a coefficient ~100x below both Slack's measurement and `H + H + M` at the same station (8.3e-45).
- **Not a verdict.** Slack's value is one bath gas (Ar), one channel (C + C), a 5,000-6,000 K window (4,936 K sits 1% below it), and ±55%. No `C + H + M` record exists. The real carbon path runs to clusters and soot, not diatomics. But the 2.1-decade margin is ten times Slack's stated uncertainty. **For the carbon-frozen floor to hold, carbon association would have to be ~100x slower than the one carbon association anyone has measured.** That leans the answer towards the exchange row, though it does not prove it.

## 8. Limits

- **Only `H + H + M` removes particles.** Unmeasured carbon three-body channels would lift both models towards equilibrium, and none would lower them. The exchange model also makes every bimolecular step instantaneous, which no rate supports. It is a limit.
- **The throat is equilibrium up to A/A* 1.05** (Da ~900 there in the bulged chamber).
- **Adiabatic and inviscid:** no wall heat, boundary layer, divergence loss or throat rounding. B6's wall-heat bracket still applies.
- **The 0-D chamber clock is not re-run.** The parcels are taken on the chamber's isentrope, as in `near_term.efficiencies`.
- **An EOS seam:** `eos_methane.janaf_correction` holds `delta` below 1000 K but drops its slope, so `e(T)` steps by ~2e-4 of `h0` there. Frozen parcels cross 1000 K; the march brackets the step (ADR-0057). Equilibrium results never see it at these exits.
