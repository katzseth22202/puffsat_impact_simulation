# Walled nozzle: the blast inside the near-term chamber, and what fixes it

Status: analysis note, 2026-10-07. Follows [`walled_nozzle_chamber_fatigue.md`](walled_nozzle_chamber_fatigue.md),
which found that a reflected blast spike, not the mean pressure, decides whether the chamber wall
survives. This note solves that spike. The chamber is the parent's near-term one: the 2.5 kg rod
at 75 km/s (7.03 GJ) in 20 m^3, methane at 7000 K (496 bar) and hydrogen at 5500 K (818 bar).

Reproduce:
- `make walled-nozzle-chamber-blast`: 1-D spherical runs, `euler2d::sphere`.
- `make walled-nozzle-vessel-blast`: 2-D axisymmetric chamber shapes, `euler2d::vessel`.
- `PYTHONPATH=python uv run python -m puffsat.walled_nozzle.volume_trade`.

## 1. Bottom line

1. **The paper's plug sends a 5-6 GPa hammer into the nozzle end on every pulse.** The 4.4 kg
   foam plug of `sec:material_chamber_plug` stops the rod but leaves 36% of its energy, 2.5 GJ,
   as the merged 6.9 kg moving on at 27 km/s. Almost nothing stands in its path: the gas on the
   axis is ~0.2 kg. In the 2-D chamber it strikes the convergent nose at **5.6-5.9 GPa** and the
   throat at **~2.8 GPa**. Steel spalls at ~2.9 GPa. The paper's own estimate, "about 0.2 GPa
   lasting microseconds", counts only the spherical part of the blast.
2. **A longer, lighter plug removes it.** Take ~10 kg of plug spread over 1.75-2.5 m of the
   axis: the paper's 4.4 kg plus ~5.6 kg, 8% of the methane charge, as foam or frozen methane.
   The rod's energy then becomes heat along the axis and no fast remnant survives. The nose
   drops to 0.4-0.8 GPa, and every wall region sits at 0.4-0.8 GPa at 1 cm cells, about half the
   point-blast level and roughly uniform. Mass matters less than length: 20 and 40 kg plugs do
   no better, and short heavy plugs load the port end.
3. **Shape the chamber like a rocket chamber, curved at both ends.** A 2:1 elliptical (domed)
   port head and an elliptical convergence give the lowest loads. A flat port face doubles the
   spike in its corner (2.8 against 1.4 GPa at 1 cm). A straight cone focuses the blast
   mid-length (1.6-1.8 GPa at 2 cm, against 0.7-0.8 on a cylinder).
4. **No carbon overwrap survives at 20 m^3.** The spike puts 0.36-0.6 of itself through the
   wall as tension. The resin holds 64 MPa across the fibres. Even the long-plug loads are about
   5-10x over that.
   - More steel ahead of the carbon does not shield it.
   - A spike longer than a layer's transit passes straight through, so the 5 mm aluminium layer
     does nothing at these pulse lengths.
   - Aluminium-alloy and magnesium backings spall earlier than steel.
   - The working wall at 20 m^3 is steel: multilayer (shrink-fitted) Cr-Mo for methane, about
     35-47 t, with ~2-3x margin on spall under the long plug.
5. **A bigger chamber trades Isp for a gentler spike.** The spike falls as 1/V at fixed pulse
   energy. Holding V/A* (so blowdown stays ~55 ms) and taking the nozzle's freeze debit, for
   methane:
   - 60-80 m^3 costs ~4-5% Isp and brings the long-plug spike to ~0.3-0.4 GPa. That is the
     range for a through-thickness-reinforced or dry Vectran wrap.
   - ~120-160 m^3 costs ~6-7% Isp and reaches ~0.16-0.21 GPa, the plain overwrap's threshold.
   - Hydrogen pays about twice the Isp for the same volume, because recombination freezes more
     at lower pressure.
6. **The liquid-share spray does not soften the wall shocks.** That is the paper's stated intent
   for it (`sec:carbon_overwrap`).
   - Swept-up mist is methane, and so is the gas around it. Both pay the same dissociation toll,
     with Γ_eff 1.13-1.27.
   - Drops of 0.5-1 mm shatter and match the gas in ~3 µs, while the front moves ~3 cm.
   - Inert extra mass at the wall raised the spike 10-20%.
   - The liquid is better spent as plug mass on the axis.
7. **Hydrogen at this pulse size is still unresolved.** Its centred-blast spike, ~2.95 GPa, is
   at steel's spall strength. A heavy axial plug either does not fit as frozen hydrogen or costs
   too much exhaust speed as polyethylene.

## 2. The solvers

- **`euler2d::sphere`** is 1-D spherical MUSCL-Hancock/HLLC, conservative in mass and energy,
  with a well-balanced pressure source. Its tests: gas at rest stays at rest; the closed sphere
  conserves mass and energy; the blast radius matches Sedov (ξ0 = 1.0328 at γ 1.4) within 1%.
- **`euler2d::vessel`** is a closed chamber wall `r_w(z)` (domed head, cylinder, cone-to-elliptical
  convergence, open throat) as a ghost-cell immersed boundary in the 2-D kernel. Its tests: gas
  at rest stays at rest in a curved chamber; a closed-cylinder blast holds its mass within 2%.
- **Gas.** The gas is γ-law at 1.2. Equilibrium methane has Γ_eff = 1 + p/(ρe) of 1.13-1.27 over
  the swept states (`eos_methane`), and 1.4 is the bracket. Each history is read against its
  own quasi-static pressure (QSP) and scaled to the solved real-gas 496 / 818 bar.
- **Deposition.** The rod's 7.03 GJ (rod + plug 6.9 kg) goes in as heat and, where stated,
  axial momentum, with energy and momentum exact on every grid.

## 3. Results

**Centred blast, 1-D sphere, grid series** (methane, γ 1.2, model units, QSP 70.7 MPa):

| cell | 2 mm | 1 mm | 0.5 mm | 0.25 mm | extrapolated |
|---|---|---|---|---|---|
| wall spike | 1.72 GPa | 2.02 | 2.23 | 2.36 | ~2.57 (~36x QSP) |

That is ~1.8 GPa for methane and ~2.95 GPa for hydrogen at the real QSP. The breathing overshoot
D is 2.0-2.3 here, against the paper's 1.7, because the reflections refocus at the centre.

**Chamber shapes, 2-D, 2 cm cells** (methane, γ 1.2, peak in GPa):

| contour (20 m^3) | port-face corner | body | convergence |
|---|---|---|---|
| cylinder + 45° cone | 1.97 | 0.71-0.80 | 0.37-0.76 |
| cylinder + elliptical nose | 1.99 | 0.68-0.81 | 0.33-0.71 |
| curved cone | 1.11 | 0.27-0.47 | 0.83-1.00 toward the throat |
| straight cone | 1.12 | n/a | **1.55-1.82 mid-length (focusing)** |

**Domed head + elliptical nose, grid series** (methane, point blast; peak with median):

| region | 2 cm | 1 cm | 0.5 cm |
|---|---|---|---|
| port face | 0.81 | 1.37 | 2.00 |
| head | 1.01 | 1.41 | 1.93 |
| cylinder | 0.82 (0.72) | 1.18 (1.02) | 1.68 (1.33) |
| nose | 0.60 | 0.69 | 1.03 |

The 2-D peaks rise ~40% per halving of the cell. Read 1 cm values as roughly half to two-thirds
of converged. **Ratios between cases are robust; absolute levels are not.**

**Deposition geometry, domed chamber, 1 cm** (methane, peak GPa; median in brackets):

| deposition | port | head | cylinder | nose | throat | worst |
|---|---|---|---|---|---|---|
| point blast | 1.37 | 1.41 | 1.18 (1.02) | 0.69 | 0.87 | 1.41 |
| **paper's plug + 27 km/s merged plug** | 0.93 | 1.22 | 0.74 (0.24) | **5.56** | **2.81** | **5.56** |
| + 15% eroded along the entry track | 1.42 | 0.67 | 0.54 | **5.86** | **2.92** | 5.86 |
| long column, all heat (ideal) | 0.85 | 0.55 | 0.83 | 0.98 | 0.29 | 0.98 |
| 2 / 5 kg port-side disk | 1.63 / 1.40 | 1.36 / 1.41 | 0.74 | 5.56 | 2.81 | 5.56 |
| 3 / 10 kg shell around the plug | 1.43 / 0.91 | 1.70 / 2.03 | 0.71-0.74 | 5.54 / 5.25 | 2.8 | 5.5 |

Sacrificial polyethylene placed beside or around the plug does nothing for the nose. Only mass
on the merged plug's path matters.

**The long plug, domed chamber, 1 cm** (methane). The plug's mass includes the paper's 4.4 kg,
and the extra comes out of the charge. The rod accretes it slice by slice: each slice releases
`p²/2 (1/m_before - 1/m_after)` as heat where it sits, and the merged body leaves the far end at
`p/(m_rod + m_plug)`.

| plug | port | head | cylinder | nose | throat | worst |
|---|---|---|---|---|---|---|
| 10 kg, 1.75 m | 0.43 | 0.39 | 0.56 (0.51) | 0.77 | 0.44 | 0.77 |
| **10 kg, 2.5 m, back-weighted** | 0.49 | 0.52 | 0.56 (0.46) | **0.50** | 0.61 | **0.61** |
| 20 kg, 1.75 m | 0.65 | 0.48 | 0.62 (0.55) | 0.72 | 0.41 | 0.72 |
| 40 kg, 1.75 m | 0.69 | 0.56 | 0.66 (0.57) | 0.65 | 0.34 | 0.69 |
| 40 kg, 2.5 m | 1.04 | 1.11 | 0.65 (0.56) | 0.75 | 0.50 | 1.11 |

18 cases in `vessel_blast_long_plug.jsonl`. One, 20 kg over 1 m back-weighted, put 1.82 GPa at
the throat, probably where this model places the whole merged body. Treat it as unconfirmed.

**Grid series for the best long plug** (10 kg, 2.5 m, back-weighted; 2 / 1 / 0.5 cm): running at commit time, `make walled-nozzle-vessel-blast` (`--long-plug-fine`). Until it lands, read the 1 cm values above as ~1.5-2x below converged, as in the point-blast series.

**Volume trade** (`volume_trade.py`). Throat grown with volume at V/A* = 134 m. The Isp change
is `sqrt(eta_net ratio)`, where eta_net is the equilibrium blowdown efficiency less the stranded
H-atom store at the Bray freeze, at A/A* = 300.

| methane 7000 K | 20 m^3 | 40 | 60 | 80 | 120 | 160 | 240 |
|---|---|---|---|---|---|---|---|
| chamber pressure (bar) | 496 | 246 | 164 | 122 | 81 | 61 | 41 |
| eta_net | 0.538 | 0.511 | 0.496 | 0.485 | 0.470 | 0.460 | 0.446 |
| Isp vs 20 m^3 | 1 | -2.6% | -4.0% | -5.0% | -6.5% | -7.5% | -9.0% |
| long-plug spike (GPa, est. converged) | ~1.25 | 0.62 | 0.42 | 0.31 | 0.21 | 0.16 | 0.10 |
| exit diameter at A/A* 300 (m) | 7.5 | 10.7 | 13.1 | 15.1 | 18.5 | 21.3 | 26.1 |
| wall area factor | 1 | 1.6 | 2.1 | 2.5 | 3.3 | 4.0 | 5.2 |

| hydrogen 5500 K | 20 m^3 | 40 | 60 | 80 | 120 | 160 | 240 |
|---|---|---|---|---|---|---|---|
| eta_net | 0.858 | 0.783 | 0.723 | 0.681 | 0.631 | 0.599 | 0.557 |
| Isp vs 20 m^3 | 1 | -4.5% | -8.2% | -10.9% | -14.2% | -16.4% | -19.4% |
| centred spike (GPa) | 2.95 | 1.48 | 0.98 | 0.74 | 0.49 | 0.37 | 0.25 |

The chamber temperature does not move with volume: it is set by energy per kilogram. The
pressure does, and the loss is the recombination race in the nozzle. Blowdown stays 54-60 ms
(methane) and 39-44 ms (hydrogen), well inside the 250 ms pulse period.

**Wall materials** (wave model, the solved spike, `wall_waves.py`):
- **Thicker steel, less carbon** (methane, 25 -> 100 mm steel): carbon tension stays 6-7x the
  resin. The steel's own radial tension rises from 0.9 to 1.5 GPa. A no-overwrap wall of ~125 mm
  (35 t at the paper's D = 1.7, ~47 t at the solved D) keeps ~2x on spall.
- **Interlayers** (5-50 mm aluminium, titanium, elastomer, syntactic foam between steel and
  carbon): none lowers carbon tension for spikes of 10 µs or longer. Soft layers cut it 3x only
  for ~1 µs spikes, and raise the steel's tension.
- **Estimates, not sourced or simulated:**
  - Through-thickness reinforcement (z-pins, stitching, 3-D weaves) bridges delaminations; it
    does not stop the resin cracking each pulse. Carrying 0.3-1.3 GPa would take roughly 15-60%
    of the fibre running through the thickness.
  - A dry (unimpregnated) Vectran or aramid wrap cannot delaminate. Instead the reflected wave
    flings its outer layers outward at ~2p/Z (~400 m/s per GPa), which strains the hoop
    ~1-2% at 0.2-0.3 GPa and ~6% at 1 GPa, against ~3-3.6% to break.
  - Aluminium (5083 class) and magnesium alloys spall at roughly 1-1.5 and 0.6-1 GPa. Magnesium
    also forms a hydride in hot hydrogen and loses strength above ~150-250 °C.

## 4. Limits

- **γ-law gas.** Real-gas equilibrium enters only through Γ_eff and the QSP scaling. A
  real-EOS run is owed.
- **Grid.** 2-D spikes are not grid-converged; the sphere's are. 2-D levels are probably 1.5-2x
  higher than the 1 cm tables show.
- **Deposition.** The long-plug model places the hot accreted mass where it was struck. A
  penetration model (rod erosion and crater flow) would refine where along the axis the heat
  goes.
- **Wall response.** The per-station ring model in `vessel_blast.py` is not credible against
  these localized impulses; it goes far past yield. Hoop fatigue under the long plug needs a
  shell model or the global breathing mode.
- **Hydrogen.** Not run with a long plug.
- **Volume trade.** Spikes are scaled as 1/V from the 20 m^3 runs, not re-solved. Nozzle
  extension mass grows with the exit diameter squared.

## 5. What the paper should change (owed to `Balloon-Pulse-Propulsion`)

- **`sec:material_chamber_plug`, the stand-off paragraph** ("the 0.2 GPa at the wall lasts
  microseconds ... the shell's response to it has not been checked").
  - Replace it: the plug as sized leaves 36% of the rod's energy (2.5 GJ) in the merged 6.9 kg
    at 27 km/s. Solved in 2-D, that strikes the convergent nose at 5.6-5.9 GPa and the throat at
    ~2.8 GPa on every pulse, about twice steel's ~2.9 GPa spall strength.
  - The fix: a longer, lighter plug. ~10 kg (the 4.4 kg foam plus ~5.6 kg of frozen methane or
    foam, 8% of the charge) over 1.75-2.5 m of the axis, denser at the back. It turns the rod's
    energy into heat along the axis, so the nose falls below ~1 GPa and the wall loads become
    roughly uniform at about half the centred-blast level.
  - The penetration-depth rule sizes the minimum plug. The plug should be longer than that
    minimum, so the merged body keeps colliding after the rod is consumed.
- **`sec:carbon_overwrap` and `tab:layered_wall_mass`.**
  - At 20 m^3 the reflected blast spike (centred ~36x the quasi-static pressure; ~1-1.5 GPa
    under the long plug) exceeds the overwrap's through-thickness strength (64 MPa resin) by
    about 5-20x everywhere. The aluminium layer and more steel do not change that.
  - The near-term wall should be steel: multilayer shrink-fitted Cr-Mo, about 35-47 t for
    methane. Otherwise the chamber must grow:
    - 60-80 m^3 (~4-5% Isp) for a reinforced or dry-Vectran wrap;
    - ~120-160 m^3 (~6-7% Isp) for a plain overwrap.
  - State that the "softens the local shocks" role given to the liquid share does not hold. The
    drops equilibrate in microseconds and pay the same toll as the gas. The liquid serves better
    as axial plug mass.
- **Chamber geometry.**
  - Draw the chamber as a rocket chamber: a 2:1 elliptical port head and an elliptical
    convergence, not a sphere or a straight cone. A flat port face doubles the corner spike; a
    straight cone focuses mid-length.
  - `fig:rod_port_plug` should show the longer plug.
- **`sec:steel_chamber_service`, `eq:ramp_dlf`.** The blast's reverberation gives D = 2.0-2.3 in
  a centred sphere, not 1.7. Size with the solved D, or with the long plug's lower, spread load.
- **Hydrogen at the 2.5 kg rod.** No wall under ~50 t survives the centred spike (~2.95 GPa, at
  steel's spall strength). Say so. The levers are a larger chamber (~8-11% Isp at 60-80 m^3) or a
  smaller rod per pulse.

## 6. Sources

- [KT07] J. R. Kamm and F. X. Timmes, "On efficient generation of numerically robust Sedov
  solutions", LA-UR-07-2849 (2007): α = 0.851 (γ 1.4, spherical).
- [PE87] M. Pilch and C. A. Erdman, Int. J. Multiphase Flow 13, 741 (1987): drop breakup time
  of 5-6 d/u √(ρ_l/ρ_g). Cited by the parent as `pilch1987_breakup`.
- Parent `templateArxiv.tex`: `sec:material_chamber_plug`, `sec:carbon_overwrap`,
  `sec:batch_filled_steel`, `tab:drop_survival`, `sec:steel_chamber_service`.
- Spall and fatigue data: as in
  [`walled_nozzle_chamber_fatigue.md`](walled_nozzle_chamber_fatigue.md) and
  [`spray_plate_fatigue_and_spall.md`](spray_plate_fatigue_and_spall.md).
