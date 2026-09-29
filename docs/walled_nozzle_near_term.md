# The near-term chamber, solved: methane at 7000 K for the 2.5 kg rod

Answers asks **A1-A6** of `todos/impact_sim_asks_methane_7000K.md`, which were
written against `Balloon-Pulse-Propulsion` @ `fe04471` (`sec:methane_7000_near_term`,
`tab:wall_pairings`, `tab:wall_pairing_doubling`, `sec:leak_radiation_trade`).
The decisions are recorded in
[ADR-0052](adr/0052-near-term-chamber-feeds-janaf-and-carbon-branches.md).
Reproduce with `make walled-nozzle-near-term`, which writes
`data/results/walled_nozzle/near_term/{chambers,efficiency,freeze,wall_heat,leak_radiation}.csv`.
Test with `make walled-nozzle-test`.

**Geometry is the paper's, not ADR-0016's.** A 2.5 kg polyethylene rod arrives
at 75 km/s (7.03 GJ) into a 20 m^3 sphere (r = 1.684 m, 35.6 m^2 of wall). The
charge is solved at the chamber temperature. The 4.4 kg polyethylene plug and
the rod are in the feed as CH2. The throats are 0.149-0.198 m^2, on a 15 degree
cone.

## Headline

1. **Methane at 7000 K delivers eta = 0.36-0.39 at A/A* = 14, measured the way the
   paper defines eta** (the rod's kinetic energy that leaves as directed exhaust,
   over the whole blowdown). It gives 0.41-0.44 at A/A* = 30 and 0.47-0.51 at 100.
   The low edge is the frozen-carbon floor and the high edge is equilibrium.
   Wall heat takes another 2-9% of that. The borrowed 53% is not reached at any
   area ratio up to 100.
2. **The paper's decision rule fires.** Methane at 7000 K is at or below ~40% on the
   frozen-carbon branch at A/A* <= 30, and hydrogen at 5500 K on the same geometry
   gives 0.64 (A/A* = 14) to 0.72 (30). The rule was "reopen hydrogen at 5500 K on
   copper".
3. **Carbon is not the binding constraint in this nozzle; exit temperature is.**
   Even in full equilibrium, 34-51% of the pulse leaves the exit still held in
   broken bonds, because the exits sit at 3400-4800 K. Freezing carbon costs only
   2-4.5 points of eta, and allowing soot to condense adds nothing out to
   A/A* = 100. The "39% locked in carbon" argument has the right size but the
   wrong mechanism.
4. **The EOS needed a JANAF correction at exactly this temperature.** Harmonic C3
   is ~30x over-bound at 6000 K. On JANAF partition functions the 7000 K charge is
   **68.8 kg, k = 27.5, 496 bar**. The paper has 62 kg, k = 24.7 and 493 bar, so
   the pressure agrees and the charge does not.

## eta: two conventions, and which one the paper means

The companion's conversion (`eta_peak`) expands the *peak* chamber state from its
stagnation enthalpy `h0 = u + p/rho`, so it credits every kilogram with the flow
work `p/rho`: +14% of `u` for methane at 496 bar, and +24% for hydrogen at 818
bar. Hydrogen's `eta_peak` reaches 0.99 at A/A* = 100 for exactly this reason. A
sealed vessel emptying adiabatically can hand out only `U0` (`integral h drho =
rho0 u0` along the isentrope; tested to 2e-4). `eta_blowdown` follows each parcel
from the chamber state at the moment it leaves. **That is the paper's `eta`**, and
every eta below is blowdown unless marked. Isp uses the impulse-averaged exhaust
speed. That sits within 0.5% of `w sqrt(eta/(k+1+P))` (eq:eta_isp), so the
paper's formula can be used with these etas.

## A1: methane at 7000 K (JANAF thermo)

| | paper | solved |
|---|---|---|
| charge / k | 62 kg / 24.7 | **68.8 kg / 27.5** (H:C = 3.79 with the PE) |
| peak pressure | 493 bar | **496 bar** |
| u | -- | 92.9 MJ/kg |
| bond share | 73% | **69%** |
| H dissociated / C free | -- | 58% / 33% |
| C3 share of store | -- | 1.7% (16% on harmonic thermo) |

| A/A* | eta_peak (companion) | **eta_blowdown** | v_impulse | Isp_true | **Isp_eff** (P = 1.76) | exit T, p |
|---|---|---|---|---|---|---|
| 4 | 0.338 | **0.282** | 7220 m/s | 736 s | **500 s** | 4811 K, 22 bar |
| 14 | 0.461 | **0.389** | 8471 | 864 | **632** | 4125 K, 4.3 bar |
| 30 | 0.521 | **0.441** | 9023 | 920 | **690** | 3803 K, 1.7 bar |
| 100 | 0.600 | **0.510** | 9712 | 990 | **763** | 3389 K, 0.40 bar |

Pitch in the exhaust (paper eq:eta_isp, pitch joining P, eta held fixed) takes
6 s off Isp_eff at 2.0 kg of pitch and 15 s at 5.4 kg, at every area ratio.

## A2: the carbon bracket (methane 7000 K, JANAF, eta_blowdown)

| A/A* | (i) frozen carbon | (ii) equilibrium | (iii) graphite ceiling | net of wall heat, (i)-(ii) |
|---|---|---|---|---|
| 4 | 0.266 | 0.282 | 0.282 | 0.242-0.277 |
| 14 | **0.360** | 0.389 | 0.389 | 0.327-0.381 |
| 30 | **0.406** | 0.441 | 0.441 | 0.369-0.432 |
| 100 | **0.466** | 0.510 | 0.511 | 0.423-0.500 |

Isp_eff on branch (i): 477 / 598 / 651 / 716 s. The three branches are defined
this way:

- **(i)** freezes every carbon-bearing species at the peak chamber composition
  per kilogram, for the whole pulse, and lets hydrogen re-equilibrate. It is the
  floor.
- **(ii)** is gas-phase equilibrium throughout.
- **(iii)** is (ii) plus graphite wherever monatomic carbon is supersaturated.
  Graphite here is NASA-Glenn C(gr) against the EOS's own C(g); it matches
  JANAF's log Kf for C(g) to 0.01 dex at 3000-5000 K. **(iii) is a thermodynamic
  ceiling**: soot at zero kinetic cost. It is not a nucleation result.

Why they sit so close: the gas never gets cold enough for carbon to matter. At
the exit, the chemical share of `u0` is 0.51 / 0.43 / 0.39 / 0.34 (equilibrium,
at A/A* = 4 / 14 / 30 / 100) and 0.55-0.42 frozen. Condensation first appears at
A/A* = 100 (20% of the carbon, at 3445 K and 0.4 bar), at the very end of the
nozzle, where its heat has no expansion left to become speed. **70% therefore
needs exit temperatures below ~2500 K, i.e. area ratios well beyond 100**, not a
faster carbon rate. Methane at 10 kK on this same geometry tells the same story: 0.30 at
A/A* = 4, 0.40 at 14 and 0.52 at 100.

## A3: hydrogen at 5500 K with the polyethylene carbon

Charge **59.9 kg, k = 24.0, 818 bar** (paper: 59 kg, 822 bar); H:C = 123 with
the PE; 18% of the hydrogen dissociated.

| A/A* | eta_peak | **eta_blowdown** | Isp_eff | frozen-carbon eta | H + H + M stranded |
|---|---|---|---|---|---|
| 4 | 0.563 | **0.470** | 752 s | 0.468 | 0 |
| 14 | 0.767 | **0.642** | 929 s | 0.639 | 0 |
| 30 | 0.865 | **0.723** | 1004 s | 0.719 | 0.7-0.8% of u0 |
| 100 | 0.989 | **0.820** | 1087 s | 0.815 | 2.1-2.3% of u0 |

The polyethylene carbon costs hydrogen nothing: the three branches agree to 0.5%.
**Hydrogen's `H + H + M` does freeze in this short cone**: at A/A* = 23-24 on
Jacobs, Giedt & Cohen, a rate that states no uncertainty and is extrapolated
above 4700 K. Beyond that station the stranded store comes off eta as in the last
column. Methane never strands more than 0.04% on the same channel, because its
atomic hydrogen is gone before the density falls.

## A4: wall heat at 7000 K methane (JANAF, equilibrium blowdown)

Throats 0.149 / 0.198 m^2; the pitch surface is held at 3900 K.

| | 0.149 m^2 | 0.198 m^2 | paper |
|---|---|---|---|
| mass e-folding time | 61 ms | 46 ms | 56 / 42 ms (its `M/(0.65 rho a A*)`) |
| 90% emptied | 151 ms | 113 ms | -- |
| radiation-weighted time `integral (T/T0)^4 dt` | **119 ms** | **89 ms** | 15-53 ms bracket |
| convection-weighted time | 44 ms | 33 ms | |
| blackbody at peak | 136 MW/m^2 | 136 | 136 |
| H- in-band optical depth floor, beam length | **>= 12** | >= 12 | "opaque" |
| floor goes thin at | 48 ms | 36 ms | |
| Bartz convection, liner | 4-43 MW/m^2 | 5-55 | 39-74 |
| **fluence per pulse** | **5.2-18.1 MJ/m^2** | **3.9-14.0** | 4.4-11.1 |
| share of pulse to the wall | 2.6-9.2% | 2.0-7.1% | |

- **Opaque: yes**, and the gas does not need carbon's help. H- bound-free alone
  gives an in-band optical depth >= 12 across the 2.25 m mean beam length, in the
  0.25-1.55 um band that holds 90% of the flux. That is a floor: H- free-free, C-
  (affinity 1.26 eV) and the C2, C3 and CH bands are left out.
- **The blowdown history is the paper's biggest miss.** Recombination holds the
  chamber hot while it empties, so `T^4` decays much more slowly than the mass.
  The radiation-weighted time is 1.6-2.1x the e-folding time, beyond the paper's
  upper bracket (`t_empty`). The high-edge fluence is therefore 60% above the
  paper's.
- **Convection is smaller than the paper's scaling** because Bartz is driven by
  `T_gas - T_wall` = 3100 K against a 3900 K pitch surface, which the `T^1.8 p^0.8`
  scaling from a 10 kK graphite-wall cell does not carry. The bracket is wide
  (frozen versus equilibrium c_p, Pr 0.5-0.8, cross-section 2e-20 to 1e-19 m^2).
- Hydrogen at 5500 K is **not** certainly opaque. The H- floor gives tau >= 0.4,
  since at 5500 K without carbon's electrons the H- population is small, so the
  paper's 52 MW/m^2 is an upper edge and the true load could be much less. The
  10 kK methane cell takes 11-30 MJ/m^2 (5.6-15% of the pulse).

## A5: pitch ablation (light: literature and thermochemistry, no ablation model)

- **Carbon melts rather than sublimes at this pressure.** The graphite-liquid-vapour
  triple point is at **107 +- 2 atm** (Haaland, *Carbon* 14 (1976) 357). The
  chamber is at 496 bar. JANAF's total carbon vapour pressure (C1 + C2 + C3)
  reaches 1 atm at ~4080 K and only ~25 bar near 4800 K, where carbon melts. A
  char surface in this chamber therefore rises past the paper's 3900 K (the 1 atm
  figure) toward melting. The triple-point temperature was not determinable in
  Haaland's measurement; ~4800 K is the commonly quoted value and is not checked
  here.
- **The boundary layer is almost carbon-saturated at the surface temperature.**
  Methane gas brought to 3900 K at 496 bar in equilibrium has free carbon at 0.86
  of saturation over graphite, and it is supersaturated below ~3700 K (2.3x at
  3000 K, where 72% of its carbon would condense). Sublimation into this gas has
  almost no chemical driving force, and cooler patches would *plate* carbon. The
  graphite 59.3 MJ/kg sublimation route the paper borrows is therefore suppressed.
- **Consequences for the paper's numbers, direction only.** A hotter surface
  (~4800 K) cuts the Bartz driving temperature by ~30% and re-radiates 30 rather
  than 13 MW/m^2. Removal by melt runoff costs far less than sublimation (fusion
  ~9 MJ/kg plus ~10 MJ/kg sensible to the melt, against 59), which *raises* the mass lost if
  the melt is sheared away. Pyrolysis gas (coal-tar pitch typically cokes at
  ~50-60% yield, a handbook figure not checked here) blows the boundary layer and blocks some convection. **Net: the 2.0-5.4
  kg per pulse cannot be confirmed or bounded from here**; it needs arc-jet data
  for pitch at hundreds of bar. No such data turned up in a literature search.

## A6: leak against wall load, re-solved at 7000 K (light: the paper's cost model)

The paper's model is used, with the gas state and the wall load **solved** at each
volume: X = 2f for leak share `f = A_entrance / A*` with the 10.9 cm port, and
Y = r/2 for wall share `r` (radiation *and* convection over the real blowdown).
V/A* is held at 101-134 m.

| V [m^3] | p [bar] | X | Y (low-high) | X + Y | H- floor tau across R/4 |
|---|---|---|---|---|---|
| 20 | 496 | 0.094-0.125 | 0.010-0.046 | 0.104-0.171 | 2.3 |
| 30 | 329 | 0.063-0.084 | 0.011-0.055 | 0.074-0.139 | 1.7 |
| 45 | 219 | 0.042-0.056 | 0.012-0.067 | **0.054-0.123** | 1.2 |
| 70 | 140 | 0.027-0.036 | 0.013-0.084 | **0.040-0.120** | 0.77 |
| 100 | 98 | 0.019-0.025 | 0.013-0.100 | 0.032-0.125 | 0.54 |
| 150 | 65 | 0.013-0.017 | 0.012-0.123 | 0.025-0.140 | 0.35 |

**The 20 m^3 chamber is leak-dominated at 7000 K too, and should grow to ~45-70
m^3 (140-220 bar).** The high edge of the total is flat at 0.09-0.12 across
45-100 m^3. The H- opacity floor across a quarter radius crosses one near 55 m^3,
the paper's opacity edge. That is a floor: the real edge lies at a larger volume.
Going from 20 to 45 m^3 cuts the total thrust cost from 10-17% to 5-12%.

## What changes in the paper

| paper item | was | now |
|---|---|---|
| tab:wall_pairings, CH4 7000 K charge / k / p | 62 kg / 24.7 / 493 bar | 68.8 kg / 27.5 / 496 bar |
| bond share, CH4 7000 K | 73% | 69% |
| emptying time | 40-53 ms | e-fold 46-61 ms; radiation-weighted 89-119 ms |
| convection | 39-74 MW/m^2 | 4-55 MW/m^2 (pitch at 3900 K) |
| heat per pulse | 4.4-11.1 MJ/m^2 | 3.9-18.1 MJ/m^2 |
| tab:wall_pairing_doubling, CH4 7000 K eta | 41 / 53 / 69% borrowed | 0.36-0.39 (A/A* 14), 0.41-0.44 (30), 0.47-0.51 (100) |
| tab:wall_pairing_doubling, H2 5500 K eta | 41 / 53 / 69% borrowed | 0.64 (14), 0.72 (30), 0.80 (100, less 2% frozen H) |
| sec:leak_radiation_trade at 7000 K | not re-solved | 45-70 m^3, 140-220 bar |

Doubling times need the paper-side synodic-lock harness and are not computed
here.

## Limits

- The 7000 K chamber is 1000 K past JANAF's tables; the partition correction is
  extrapolated there.
- The expansion is quasi-1D, steady per parcel, and adiabatic. Wall heat is
  reported beside eta, not folded into the isentrope. There is no divergence or
  boundary-layer loss.
- Carbon association has no rate, so the branches bracket it.
- The graphite ceiling ignores melting and liquid carbon.
- The H- opacity is a floor over part of the spectrum.
- A5 and A6 are light by agreement: no ablation model, and the paper's own
  thrust-cost model for A6.

## Addendum: area ratios to 1000, and where recovery stops paying

`AREA_RATIOS` now runs 4-1000 (RL10B-2 flies 285:1; high-expansion studies run 500-1000:1). On
these throats 300:1 is a 7.5-8.7 m exit and 1000:1 is 13.8-15.9 m. `freeze.csv` now carries
`eta_net`: equilibrium eta_blowdown less the atomic-H store stranded at the Bray freeze point.

| A/A* | CH4 7000 K eta_net (frozen-C floor) | CH4 bond store recovered | H2 5500 K eta_net | H2 bond store recovered |
|---|---|---|---|---|
| 14 | 0.389 (0.360) | 37% | 0.642 | 85% |
| 30 | 0.441 (0.406) | 43% | 0.716 | 91% |
| 100 | 0.511 (0.466) | 50% | 0.798 | 91% |
| 300 | 0.539 (0.511) | 52% | 0.859 | 94% |
| 500 | 0.550 (0.529) | 53% | 0.880 | 94% |
| 1000 | 0.563 (0.551) | 54% | 0.902 | 94% |

**Chemical recovery stops paying where H + H + M freezes**: at A/A* ~ 100 for methane and
~ 23 for hydrogen, on these 0.15-0.2 m^2 throats. Past that, only thermal energy converts.
- **Hydrogen's chemistry is essentially done by A/A* ~ 30** (91% of its bond store back). It
  keeps gaining eta to 1000 because 64% of its energy is thermal.
- **Methane gets only half its bond store back.** The carbon-bearing half stays bound in C2H2,
  C3 and atoms at exit temperatures of 2800-3400 K. Beyond A/A* = 300 it gains about 0.01 in eta
  per step. A/A* ~ 300 is methane's knee; hydrogen's knee is at A/A* ~ 500-1000.

These efficiencies **exclude the entrance leak** (2f = 9-13% of thrust with an open 10.9 cm port;
see A6). A port door closing within ~3 ms of the rod's passage cuts that to under 1% (the
open-port choked flux is 84-101 kg/s against a 67-76 kg pulse).

## Addendum: the coated steel wall (methane 7000 K)

`make walled-nozzle-wall-layers` runs a 1-D transient solve through coating + Cr-Mo steel on the
solved blowdown flux. It is verified against the convective erfc solution, and the paper's own
pitch model under this flux gives steel 395-441 K and 30-120 um/pulse (paper: 450-490 K,
43-117 um).

- **1 mm of graphite protects the steel for one pulse only if it is a porous spray**:
  k ~ 3 W/m/K gives steel 490-587 K, and 2 mm gives ~325 K. Dense graphite (k = 40) melts the
  steel at the high edge.
- **Graphite retains 3.1-6.8 MJ/m^2 per pulse (110-243 MJ over the wall); pitch retains
  0.8-1.4 MJ/m^2 (28-50 MJ).** The 68.8 kg methane charge warmed to 500 K absorbs ~100 MJ.
  Only pitch's heat fits the between-pulse charge spray, so a graphite-only wall needs coolant
  channels.

## Conclusions to bring back to the paper

Consolidated 2026-09-29. The working log with every number is in the (gitignored)
`todos/impact_sim_asks_methane_7000K.md`, sections 6-7. Items marked *(estimate)* are hand or
scratch-script analyses, not `make` targets.

**Performance** (eta_net = blowdown eta less hydrogen freezing; Isp_eff from eq:eta_isp,
P = 1.76, open-port leak excluded):

| A/A* | CH4 7000 K eta_net | CH4 Isp_eff | H2 5500 K eta_net | H2 Isp_eff |
|---|---|---|---|---|
| 14 | 0.36-0.39 | 602-635 s | 0.64 | 934 s |
| 100 | 0.47-0.51 | 720-766 s | 0.80 | 1075 s |
| **300** | **0.49-0.54** | **740-793 s** | **0.86** | **1126 s** |
| 1000 | 0.50-0.56 | 754-817 s | 0.90 | 1162 s |

1. **Replace the borrowed 41/53/69%** with these. Methane's recombination freezes at A/A* ~ 100
   with half its bond energy back; hydrogen's chemistry is done by ~30. **Design point
   A/A* ~ 300** (7.5-8.7 m exit, a C-C extension in the RL10B-2 family).
2. **The decision rule fires.** Methane at 7000 K is at or below ~40% on the frozen-carbon floor
   for A/A* <= 30. Hydrogen at 5500 K leads by ~0.2 yr of doubling on the paper's own
   tank-inclusive table *(estimate, interpolated)*. Methane's case is handling, not performance.
3. **Neither a colder methane chamber, water, nor a bigger lower-pressure chamber fixes methane's
   eta** *(scratch scans)*.
   - Colder methane raises impulse per rod (1.7x at 5000 K) but lowers Isp_eff (to ~490 s at
     A/A* = 14).
   - Water matches methane's eta and loses ~120 s of Isp.
   - Lower pressure loses eta: high pressure helps recombination both kinetically and
     thermodynamically.
4. **Port door** *(estimate)*: the open 10.9 cm port costs 9.4-12.5% of thrust (~13-16% of
   Isp_eff). A flow-closing petal valve shut by the chamber's own pressure, reopened by springs,
   closing within ~3 ms with a gap of <= 1 mm, cuts that to ~1%. It is essential for methane
   (without it CH4 ties the reactor). The load once shut is 390-640 kN. With the door, there is no
   reason to grow the chamber.
5. **Plug stand-off** *(estimate)*: when the rod hits the plug, `v_cm / v_expansion =
   sqrt(m_rod / m_plug)`. With the paper's heavier plug, the fireball's rear drifts *toward* the
   port (~9 km/s). **Hold the plug near the chamber centre (~1.68 m from the port)**, or the port
   and door take GPa-class loads every pulse. The ~0.2 GPa transient at the wall (~4x the static
   design pressure) needs a shell check.
6. **Walls.**
   - *Methane:* keep pitch on steel (ADR-0053). Graphite is chemically compatible, but it retains
     more heat per pulse than the whole charge can absorb.
   - *Hydrogen:* a bare copper skin survives the pulse (surface <= 1110 K) *(estimate: Duhamel on
     the solved flux)*; never a graphite coat.
   - Carbon melts rather than sublimes at 496 bar (triple point 107 atm).
7. **Table corrections** (CH4 7000 K):
   - charge 68.8 kg, k = 27.5, 496 bar, bond share 69%;
   - convection 4-55 MW/m^2; heat per pulse 3.9-18.1 MJ/m^2;
   - emptying e-fold 46-61 ms, radiation-weighted time 89-119 ms (the paper's 25-33 ms and
     16-21 ms are both wrong).
