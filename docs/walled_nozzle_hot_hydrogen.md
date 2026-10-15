# Walled nozzle: the bulged hydrogen chamber at 7000 K and 10 000 K (asks H1-H5)

Status: analysis note, 2026-10-14. Answers `aim_is_all_you_need`'s
[`docs/hydrocode_asks_2026-10-14_hot_hydrogen.md`](https://github.com/katzseth22202/aim_is_all_you_need/blob/main/docs/hydrocode_asks_2026-10-14_hot_hydrogen.md)
(its ADR 0063). Method and limits are [ADR-0058](adr/0058-hot-hydrogen-thick-gas-and-lighter-fills.md).
Reproduce with `make walled-nozzle-hot-hydrogen` (the 2-D rod-and-plug solves take ~1 h on 8
cores, then ~25 min). Outputs are in `data/results/walled_nozzle/near_term/hot_hydrogen/`. The
summary for the calcs repo is [`hot_hydrogen_handoff.md`](hot_hydrogen_handoff.md).

The chamber is B5's bulged hydrogen chamber (`walled_nozzle_hydrocode_asks.md` §6a): 106 m^3 at the
2.5 kg rod, the r 3.0 m bulge, a 0.84 m^2 throat, AR 300 on a 15° cone, polyethylene rod and 4.4
kg plug. Only the charge changes: it is sized to hold each temperature.

## 1. Bottom line

| ask | answer |
|---|---|
| **H1** peak state | **k = 9.65 at 7000 K (24.1 kg) and 7.47 at 10 000 K (18.7 kg)**, against the parent's 9.5 and 5.5. The rod and plug (6.9 kg of polyethylene) are 22-27% of the hot chamber's mass and do not store hydrogen's ~330 MJ/kg, so the 10 000 K charge is a third heavier than estimated. Peak pressure stays near today's: **127 / 155 bar** against 142. The gas is ideal (coupling 0.08-0.17) and **optically thick** (Rosseland depth 7 / 65 across the chamber, TOPS). |
| **H2** recombination | AR 300 returns **60% of the bonds at 7000 K and 55% at 10 000 K** at the measured rate (82% today). That is not the 94% the parent's best row took from the dense 5500 K sphere. Even full equilibrium returns only 73% / 65%. **Net impulse 340 / 285 kN s per pulse** (683 today, finite rate); 328-416 / 275-345 over 74-81 km/s. Carbon is a minor term here (≤1%). |
| **H3** wall | **The impact favours 7000 K; the heat does not.** The 2-D rod-and-plug solve with the 24 kg fill passes the dry wrap at **29.8 t, fling 1.77%** (today 37.0 t, 2.20%). **The 18.7 kg fill fails it: 3.48% in the taper, a 2.24 GPa hit on the cone.** Radiation is close to `sigma T^4` (136 / 567 MW/m^2 at peak), and the wall takes 3.0-4.7% / 7.6-8.1% of the pulse. **That is 1.0-1.5x / 3.2-3.4x what the charge can carry away as regen coolant** (today 0.04-0.57x). The copper face rises 383 / 811 K per pulse. |
| **H4** rate | **3.3 Hz at 7000 K and 3.6 Hz at 10 000 K**, sealed at 0.3 bar with a 0.5 t/s fill (today 2.3 Hz). Sealed at fill pressure, **3.8 / 4.0 Hz**; with a 1 t/s fill, up to 4.1 / 4.3 Hz. The residue is hotter (3200-3600 K) but there is less of it: **the mix stays at 100-290 K, so a polymer film still seals.** |
| **H5** sweet spot | Impulse per kilogram carried rises 12% from today to 7000 K and is then **flat to 10 000 K** (9.57 / 9.66 / 9.63 km/s at 7000 / 9000 / 10 000 K). **Thrust falls monotonically**: 1.50 MN today, 1.07 MN at 7000 K, 0.95 MN at 10 000 K at their own rates. **7000 K is the hottest worth flying.** Hotter gains no exhaust speed, loses thrust, doubles the wall heat and fails the wrap. |

**So 7000 K helps only modestly, and 10 000 K does not help.** At 7000 K the effective exhaust
rises 12% (the parent's convention: 10.3 to 11.5 km/s), but only 60% of the bonds come back. The
faster cycle recovers only part of the halved impulse per pulse: thrust is ~29% below today's,
like for like. Whether that nets a faster doubling is the calcs repo's chain to run.

**Correction owed for today's chamber:** the finite-rate impulse of today's bulged chamber (0.84
m^2, AR 300) is **683 kN s**, not 657 (B5's Bray sudden freeze; the same overstatement S26 found
for methane). It is 661 with the high wall-heat edge.

## 2. H1: the peak state

`near_term.solve_charge` at 106 m^3: the charge (plus the polyethylene rod and plug) that the
pulse energy heats to the temperature. `eos_methane`'s partition functions, Saha for H and C, no
fitted polynomial extrapolated.

| | charge | k | rho | p | u | bonds | H free | C free | ionised | Gamma_eff | H:C | tau_R (TOPS) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| today (4884 K) | 59.9 kg | 24.0 | 0.630 | 142 bar | 105 MJ/kg | 44% | 0.23 | 0.09 | 2e-7 | 1.215 | 123 | (below TOPS) |
| 5500 K | 40.3 | 16.1 | 0.446 | 127 | 149 | 57% | 0.43 | 0.36 | 2e-6 | 1.192 | 83 | |
| 6000 K | 31.6 | 12.7 | 0.363 | 123 | 183 | 63% | 0.59 | 0.62 | 8e-6 | 1.186 | 66 | 2.1 |
| **7000 K** | **24.1** | **9.65** | 0.293 | **127** | 227 | 68% | 0.82 | 0.86 | 6e-5 | 1.191 | 51 | **7.1** |
| 8000 K | 21.3 | 8.53 | 0.266 | 136 | 249 | 67% | 0.91 | 0.93 | 3e-4 | 1.205 | 45 | 16 |
| 9000 K | 19.8 | 7.92 | 0.252 | 146 | 263 | 66% | 0.95 | 0.96 | 8e-4 | 1.220 | 42 | 33 |
| **10 000 K** | **18.7** | **7.47** | 0.241 | **155** | 275 | 64% | 0.97 | 0.98 | 2e-3 | 1.234 | 40 | **65** |

Over the burn (k follows w, as the parent holds it): 7000 K is 23.4 / 29.1 kg (k 9.34 / 11.64) at
123 / 150 bar for 74 / 81 km/s. 10 000 K is 18.1 / 22.2 kg (k 7.25 / 8.90) at 151 / 182 bar.

- **Why k is not 5.5 at 10 000 K.** The parent held `(k + 1 + P) u(T) = w^2/2` with pure
  hydrogen's `u(T)`. Pure hydrogen at 0.57 kg/m^3 does hold 245 / 329 MJ/kg here (Cantera's 242 /
  340-390, less its 300 K zero). But the rod and plug are ~6 kg of carbon, and carbon stores only
  ~70 MJ/kg even when fully atomised. Above ~8000 K the hydrogen is fully dissociated, so its
  energy grows only as `3/2 kT` per atom. **The charge stops shrinking:** 24.1 kg at 7000 K, 18.7
  kg at 10 000 K.
- **The pressure question is settled:** the smaller charge offsets the temperature, as the ask
  guessed. 127-155 bar against 142.
- **Real gas:** the electron-ion coupling `e^2/(4 pi eps0 a kT)` is 0.08 / 0.17, so the ideal-gas
  equilibrium EOS stands.
- **Optically thick.** TOPS gray Rosseland means for the two H:C mixtures (LANL OPLIB, pulled
  2026-10-14; `make walled-nozzle-hot-hydrogen-tops`) give 7 m^2/kg at 7000 K and ~50 m^2/kg at
  10 000 K at these densities. That is depth 7 / 65 across the 3.4 m beam length. The in-band H-
  floor, the only opacity this study had for this gas, gives 1.9 / 21: thick too, but less so.

## 3. H2: recombination and impulse

`finite_rate.pulse` (ADR-0057) unchanged: H + H + M at the measured rate, every parcel of the
blowdown, quasi-1D 15° cone from A/A* 1.05, carbon bracketed by its two models. Impulse is
`m_rod ((1 - 0.01)(1 + k + P) u - w)`, with `u` the shifting speed times the finite-rate
vacuum-thrust ratio.

| | eta: eq / finite / frozen | bond returned: eq / **finite** / frozen | first Da < 1 (peak) | **I [kN s]**: eq / **finite** / carbon frozen |
|---|---|---|---|---|
| today | 0.865 / 0.839 / 0.691 | 1.00 / **0.82** / 0.21 | A/A* 14 | 701 / **683** / 681 |
| 5500 K | 0.780 / 0.757 / 0.576 | 0.93 / 0.71 / 0.15 | 22 | 523 / 504 / 500 |
| 6000 K | 0.727 / 0.711 / 0.521 | 0.83 / 0.65 / 0.14 | 27 | 432 / 417 / 412 |
| **7000 K** | 0.683 / 0.672 / 0.478 | 0.73 / **0.60** / 0.12 | 32 | 351 / **340** / 336 |
| 8000 K | 0.670 / 0.660 / 0.467 | 0.69 / 0.58 / 0.11 | 34 | 321 / 311 / 307 |
| 9000 K | 0.666 / 0.656 / 0.467 | 0.67 / 0.56 / 0.11 | 35 | 305 / 296 / 292 |
| **10 000 K** | 0.664 / 0.654 / 0.469 | 0.65 / **0.55** / 0.10 | 35 | 294 / **285** / 281 |

| | 74 km/s | 75 km/s | 81 km/s | rate / 3.16 | rate x 3.16 |
|---|---|---|---|---|---|
| 7000 K, I [kN s] (bond) | 328 (0.60) | 340 (0.60) | 416 (0.62) | 334 (0.56) | 344 (0.64) |
| 10 000 K | 275 (0.55) | 285 (0.55) | 345 (0.57) | 280 (0.51) | 288 (0.58) |

- **The 94% does not transfer.** It was the 5500 K sphere at 818 bar (3 kg/m^3). The hot bulge
  is ten times thinner, and its gas holds 2-3x the energy per kilogram. **Even shifting
  equilibrium returns only 73% / 65%**: the exit, at A/A* 300, is still hot enough that hydrogen
  there is partly dissociated. Finite rate costs another 13 / 10 points.
- **The freeze is not sudden.** The throat is in equilibrium (Da ~800-1100 at A/A* 1.05) and the
  race is lost at A/A* 32-35, later than today's 14. The hot store is so large that the stranded
  tail is still large.
- **On the parent's own scale** the returned share is 60% at 7000 K, between its 70% and 50% rows,
  and 55% at 10 000 K. Its eta should come from these impulses, not from the ceiling formula.
- **Carbon is a minor term.** H:C is 40-123, and freezing carbon costs ≤1.2% of impulse.
- **Impulse per pulse halves at 7000 K** (683 → 340) and falls 2.4x at 10 000 K. That is
  `sqrt(2 m E)` with a smaller expelled mass and a lower eta.

## 4. H3: the wall

### 4a. The impact: dry wrap under the lighter fills

`sweep --vessel-blast --hot-hydrogen`: the r 3.0 m bulge solved in 2-D with the rod and plug as
material (0.5 cm, 1.2 ms, gamma 1.2: real `Gamma_eff` 1.19 / 1.23), filled with 18.5 / 13.1 kg of
model gas (the charge less the 10 kg model plug's 5.6 kg excess, as the 59.9 kg run). Scored by
`zone_loads` at each charge's fixed-charge QSP; the limit is Kevlar's 2.4% flung hoop strain.

| fill | QSP | head | bulge ramp | plateau | **taper** | **cone** | dry wrap |
|---|---|---|---|---|---|---|---|
| today, 59.9 kg | 142 bar | 0.55 GPa, 1.17% | 0.44, 1.91% | 0.48, 1.77% | 1.04, **2.20%** | 1.19, 1.08% | **37.0 t, passes** |
| **7000 K, 24.1 kg** | 127 bar | 0.60, 1.05% | 0.14, 0.84% | 0.40, 1.39% | 0.81, **1.77%** | 1.27, 1.48% | **29.8 t, passes** |
| **10 000 K, 18.7 kg** | 155 bar | 0.59, 1.14% | 0.15, 0.82% | 0.31, 1.68% | 1.33, **3.48%** | **2.24**, **2.51%** | **36.3 t, fails** |

(peak pressure, worst-station flung strain)

**Revised 2026-10-15:** these wraps are sized for strain only. Sized also so that the liner's
hoop-swing crack lasts 4,200 pulses on the hydrogen curve
([`walled_nozzle_hydrocode_asks.md`](walled_nozzle_hydrocode_asks.md) §9a), they are **43.8 t**
today, **34.2 t** at 7000 K (fling 1.73%) and 41.4 t at 10 000 K (still failing on fling,
3.38%). The ranking does not change.

- **7000 K has more margin than today** (2.4/1.77 = 1.36x against 1.09x) on a wrap 7 t lighter,
  sized to a lower QSP. The lighter fill weakens the reflections that load the bulge ramp and
  taper.
- **10 000 K fails by 1.45x**, more than the grid's 15-25% per halving. Two things compound:
  - the lightest fill lets the blast remnant reach the cone at 2.24 GPa (today 1.19), so the cone
    flings too;
  - the 22% higher QSP scales every history up.
- **Solid steel in the cone still cracks in hydrogen within ~70 pulses** at 7000 K (80 today).
  The dry wrap is the design for exactly that reason.
- **Not captured:** the 2-D solve is gamma-law and adiabatic. The real thick gas radiates during
  the 1.2 ms, which would soften the spike a little; it cannot rescue 1.45x.

### 4b. The heat

The 0-D blowdown (`blowdown.integrate`) with the thick-gas loss: radiation `sigma T^4 (1 -
exp(-tau_R))` from the core (the in-band H- floor below TOPS's 5802 K), plus Bartz convection on
an 800 K copper face, bracketed frozen to equilibrium `c_p` (`near_term.bartz_liner`).

| | peak radiation | liner convection | throat convection | per pulse, radiation + convection | wall share | heat / kg charge (coolant: 9.0 MJ/kg) | Cu rise per pulse | impulse debit |
|---|---|---|---|---|---|---|---|---|
| today | ~0 MW/m^2 | 5-64 | 125-1,460 | 0.0 + 0.2-2.5 MJ/m^2 | 0.4-4.4% | 0.4-5.1 MJ/kg | 202 K | 2.5% |
| 5500 K | 6 | 6-74 | 137-1,690 | 0.0 + 0.2-2.5 | 0.4-4.5% | 7.9 | 222 | 2.6% |
| 6000 K | 65 | 7-68 | 151-1,550 | 0.4 + 0.2-2.1 | 1.1-4.3% | 9.6 | 238 | 2.4% |
| **7000 K** | **136** | 8-46 | 190-1,050 | **1.5** + 0.2-1.1 | 3.0-4.7% | **9-14** | **383** | 2.5% |
| 8000 K | 232 | 10-34 | 234-780 | 2.4 + 0.2-0.7 | 4.6-5.5% | 18 | 502 | 2.9% |
| 9000 K | 372 | 12-32 | 279-720 | 3.2 + 0.2-0.5 | 6.2-6.7% | 24 | 638 | 3.5% |
| **10 000 K** | **567** | 14-34 | 323-770 | **4.1** + 0.2-0.5 | 7.6-8.1% | **29-31** | **811** | 4.3% |

- **Radiation changes the wall's problem.** Today's wall heat is convection and is bracketed
  10x wide. At 7000 K and above it is radiation, and it is not bracketed: the gas is thick, so
  the wall sees the core's `sigma T^4`.
- **The charge cannot cool the wall.** Taken as regen coolant from the 77 K fill to a 700 K
  outlet, hydrogen absorbs 9.0 MJ/kg.
  - Today it absorbs the wall heat with margin.
  - At 7000 K it falls 1.0-1.5x short (1.08 GW of wall heat at 3.3 Hz, against 0.71 GW).
  - At 10 000 K it falls 3.2-3.4x short (2.07 against 0.61 GW).
  The rest needs another heat path: extra coolant mass, a radiator, or a hotter outlet than
  copper allows.
- **The copper face.** Duhamel on the solved flux history gives a 383 K surface rise per pulse at
  7000 K, on top of a regen-cooled base. On a ~500 K base that is ~880 K: survivable for GRCop-84,
  marginal for OFHC. At 10 000 K the rise is 811 K, so the face reaches ~1,300 K every pulse,
  close to copper's 1,358 K melt. **Copper is not a 10 000 K liner.**
- **Over the burn** (4500 pulses at the H4 rate): 8.7 MW/m^2 burn-average and 12 GJ/m^2 at
  7000 K, against 5.6 MW/m^2 and 11 GJ/m^2 today; 16.6 MW/m^2 and 21 GJ/m^2 at 10 000 K. The burn
  takes 23 min at 7000 K (33 min today).
- **The throat is not the discriminator.** Its convective flux, 0.1-1.7 GW/m^2 on the bracket, is
  the same order today and hot. It is set by `p^0.8`, and the pressures are alike. It is an open
  item for every hydrogen row, not a new one.
- **The impulse debit stays small** (2.5% / 4.3%), because the heat leaves after most of the
  impulse has.

## 5. H4: the cycle

`cycle.seals` on the blowdown with the thick-gas wall heat: seal, films (30 ms), fill (0.5 t/s
at 2.5 kg, `cycle.INJECTOR`, or 1 t/s), lock (20 ms).

| | fill pressure | seal at fill pressure | rate (0.5 / 1 t/s) | seal at 0.3 bar | **rate (0.5 / 1 t/s)** | mix at 0.3 bar | film |
|---|---|---|---|---|---|---|---|
| today | 1.79 bar | 181 ms | 2.85 / 3.43 Hz | 270 ms | **2.27** / 2.63 Hz | 90-91 K, 2.1 bar | 0.56 kg |
| 5500 K | 1.21 | 176 | 3.26 / 3.75 | 245 | 2.66 / 2.98 | 98-103 K | 0.44 |
| 6000 K | 0.95 | 172 | 3.50 / 3.94 | 227 | 2.94 / 3.24 | 102-121 K | 0.41 |
| **7000 K** | 0.72 | 167 | **3.77** / 4.14 | 206 | **3.28** / 3.57 | 107-153 K, 1.0-1.5 bar | 0.39 |
| 8000 K | 0.64 | 165 | 3.88 / 4.23 | 198 | 3.44 / 3.72 | 109-170 K | 0.39 |
| 9000 K | 0.59 | 164 | 3.94 / 4.28 | 193 | 3.54 / 3.80 | 110-180 K | 0.38 |
| **10 000 K** | 0.56 | 163 | **3.99** / 4.31 | 190 | **3.61** / 3.87 | 112-189 K, 0.8-1.4 bar | 0.38 |

- **Today's row reproduces B2's 2.2-2.4 Hz** at a 0.3 bar seal, which is the rate the parent flies.
- **The rate rises only 1.4-1.6x while the impulse per pulse falls 2.0-2.4x.** The fill shrinks
  2.5-3x (120 → 48 → 37 ms), but the blowdown to a 0.3 bar seal shortens only 24-30% (270 → 206 →
  190 ms), since `V/A*` is fixed. The films and lock are a fixed 50 ms.
- **The film still seals behind the hotter residue.** The residue is hotter (3200-3600 K against
  1750-2700 K) but smaller (0.2-0.5 kg against 0.5-1.9 kg), and it joins a cold fill. The high
  bound of the mix is 150-290 K at either seal pressure, far below Kapton's 673 K.
- **The fill pressure is now near 0.3 bar.** At 0.56-0.72 bar, sealing at fill pressure is
  barely earlier than sealing at 0.3 bar. Which seal pressure counts as "weak enough" for laying a
  film over the throat jet is the open mechanism question (B2). The table gives both.
- **4 Hz needs the fast fill:** 4.1 Hz at 7000 K and 4.3 Hz at 10 000 K, sealed at fill pressure
  with 1 t/s.

## 6. H5: the sweet spot

Exchange carbon, finite rate, wall heat debited. Rate: sealed at fill pressure (the film is
cold-safe everywhere), 0.5 t/s. Carried per pulse: charge with its 0.205 kg/kg tank, the plug and
the throat film (the rod is counted against the PuffSats, as the parent does).

| | I per pulse | I / (charge + plug) | **I / carried** | thrust at the sealed-at-fill rate | thrust, both at 0.3 bar |
|---|---|---|---|---|---|
| today | 661 kN s | 10.27 km/s | 8.51 km/s | 1.88 MN | 1.50 MN |
| 5500 K | 485 | 10.85 | 9.02 | 1.58 | 1.29 |
| 6000 K | 402 | 11.17 | 9.30 | 1.41 | 1.18 |
| **7000 K** | **327** | **11.45** | **9.57** | **1.23** | **1.07** |
| 8000 K | 296 | 11.52 | 9.65 | 1.15 | 1.02 |
| 9000 K | 279 | 11.51 | **9.66** | 1.10 | 0.99 |
| 10 000 K | 264 | 11.46 | 9.63 | 1.05 | 0.95 |

- **Impulse per kilogram carried peaks near 9000 K, but the curve is flat past 7000 K** (+0.9%
  from 7000 to 9000 K). Thrust falls 7% over the same span, the wall heat rises 50%, and 10 000 K
  fails the wrap. **7000 K is the hottest temperature worth carrying.** Whether a cooler one (6000
  K: -3% per kg, +10% thrust) flies the departure better is the chain's call.
- **What the ask's estimate got right and wrong.**
  - Right: the shape (hotter thins the charge until lost thrust outweighs exhaust) and the
    pressure offset.
  - Wrong at 7000 K: 11.5 km/s against its 12.1-13.6, because the bonds return 60%, not 70-94%.
  - Wrong at 10 000 K: 11.5 km/s against 11.8-15.5. The charge is a third heavier and returns
    only 55%.

## 7. Limits

- **2-D wall:** gamma-law, adiabatic, one grid (0.5 cm). The bulge's peaks rose 15-25% per halving
  (chamber-blast note §3c). The 7000 K pass (1.36x) and the 10 000 K fail (1.45x) both lie outside
  that, but a 0.25 cm check of the 7000 K fill is owed, as it is for today's hydrogen fill.
- **Opacity:** TOPS is atomic only (no H2 bands) and starts at 5802 K. The 7000 K core is 82%
  dissociated, so this is fine there. The cooled boundary layer's own absorption is not modelled:
  it would trim the radiative flux a little, not by the factor the regen balance needs.
- **Convection** on the copper face is the study's Bartz bracket (`near_term.bartz_liner`) on an
  assumed 800 K face, and the liner's base temperature is not solved.
- **Coolant capacity** counts only the charge (77 K → 700 K). Liquid storage at 20 K adds ~1
  MJ/kg of headroom (boil-off and warming to 77 K).
- **Mechanism times** (films, lock, fill rate) are B2's assumptions.
- The finite-rate march and the blowdown carry ADR-0057's and ADR-0056's limits unchanged.
