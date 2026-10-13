# Walled nozzle: blowdown, the 2 Hz cycle, wall forces and a cushioned wall

Status: analysis note, 2026-10-12. Answers `aim_is_all_you_need`'s
[`docs/hydrocode_asks_2026-10-12.md`](https://github.com/katzseth22202/aim_is_all_you_need/blob/main/docs/hydrocode_asks_2026-10-12.md)
(asks B1-B6, from its ADRs 0049, 0053 and 0054). Builds on
[`walled_nozzle_chamber_blast.md`](walled_nozzle_chamber_blast.md): the survivable chamber is its
r 3.0 m bulge with a dry Kevlar wrap, scaled to the 5 kg rod (212 m^3, 1.33 m^2 throat). Scope and
conventions are [ADR-0056](adr/0056-hydrocode-asks-blowdown-clock-and-fixed-charge.md).

Reproduce (`make` targets): `walled-nozzle-blowdown` (B1, B6), `walled-nozzle-cycle` (B2, B5 rate),
`walled-nozzle-zone-loads` (B3, B5 wall), `walled-nozzle-cushion` (B4), `walled-nozzle-blowdown-2d`
(the 2-D check). Outputs are under `data/results/walled_nozzle/near_term/` (`blowdown/`,
`zone_loads.csv`, `cushion.csv`).

The summary for the parent calcs repo is [`hydrocode_asks_handoff.md`](hydrocode_asks_handoff.md).

## 1. Bottom line

| ask | answer |
|---|---|
| **B1** blowdown | The parent's estimate holds. Real-gas, the 212 m^3 chamber is 95% empty at **221 ms** and at fill pressure at **281 ms** (parent: 215-260 and 278-336). The gas stays in equilibrium (2.8-3.6 decades of margin). Wall heat takes 3-14% of the pulse and moves the times by under 10%. |
| **B2** 2 Hz | **Closes at 212 m^3 with ~60 ms to spare, if the throat film is laid at fill pressure.** It is not that simple: the film traps ~3 kg of 4,100 K residue, which heats the next charge to 600-1,100 K at 1.6-2.9 bar. A polymer film rated to 400 °C must instead wait for ~0.3 bar (389 ms), and the cycle falls to 1.8 Hz. The fastest each chamber can run is in §3. |
| **B3** wall forces | Without the bulge the plug-band spike grows **5-8x**: 0.44 -> 2.3 GPa and 9 -> 46 kPa s at 106 m^3, **3.5 GPa** and 77 kPa s at 53 m^3 (at steel's ~2.9 GPa spall strength). The dry wrap's hoop overshoot `D` rises from 1.8-3.6 in the band to 8-10. **No unbulged chamber's dry wrap survives.** |
| **B4** cushion | **No: none of 390 design-chamber cases lasts 4,200 pulses (the best, ~190).** The solved spike is **9-18 kPa s** (bulged) to **46-77 kPa s** (unbulged), not the 2-5 the estimate assumed, and it is uneven over a few centimetres. A tile on a soft gas foundation therefore bends at **0.9-8 GPa** (10-50 mm steel, 0.1-0.5 m tiles). The cushion does cut the backwall's `D` in the unbulged 106 m^3 chamber (9.9 -> 2.9), but needs 30 cm of cushion and 15-18 cm of travel. |
| **B5** hydrogen | **The r 3.0 m bulge holds hydrogen too:** its dry wrap passes at **37 t** (2.5 kg rod), on a thin fling margin (2.20% against 2.4%). Hydrogen brings a new 1.2 GPa hit in the cone. At the parent's 60 kg charge the bulge costs **6-8% impulse per rod** (624 / 658 against 666 / 697 kN s), mostly through the nozzle freeze. It needs a **0.75-1.0 m^2 throat** (17-20 m exit at area ratio 300) to cycle at **2.1-3.1 Hz**. **It does not keep 4 Hz.** The sphere keeps 4 Hz only with a ~1 t/s fill. |
| **B6** impulse | Smaller is better, and more so with wall heat: **861 / 901 / 925 kN s** (AR 100) and **892 / 933 / 972** (AR 300) at 212 / 106 / 53 m^3. With the high wall-heat edge, 761 / 848 / 902. **The parent's 906 / 971 kN s for 212 m^3 charged neither the gate nor the nozzle freeze.** Corrected: 861 / 892. |
| **Bulge + gas bags?** | **Keep the bulge and the dry wrap; do not add cushions.** In the bulged chamber the wrap already passes, and the best cushion's backwall `D` (2.5) is no better. The tiles still fail in bending. In a smaller chamber the cushion is what fails. Volume is worth 5-18% impulse per rod, but stand-off is the only thing that cuts the spike. The gain needs a smaller chamber that still stands the wall off from the plug (§8). |

**Decided 2026-10-13: double the throat to 2.66 m^2** (`V/A*` ~80 m, the 5 kg chamber), §9.
It is the timing lever, and it helps recombination rather than costing it. Its wall passes: the
dry wrap is 46.2 t against 52.4 t, with a 0.99 GPa hit on the shorter cone (§9). The combined
compact shape (159 m^3) fails the wrap at the taper and is not adopted.

**Corrections owed to the parent:**
1. `tab:one_chamber`'s 906 / 971 kN s are gross of the 1% gate and the H + H + M freeze. With both,
   861 / 892 kN s (§7).
2. ADR 0054's cushion rests on a 2-5 kPa s spike. The solved one is 2-15x that, and uneven (§5).
3. ADR 0053's fill-pressure gate traps a hot residue. The films hold 1.4-2.6x the cold fill
   pressure, and a 1.33 m^2 throat film is 1.5-5 kg per shot, not grams (§3).

## 2. B1: the blowdown

`blowdown.py` solves the 0-D blowdown on the real equilibrium gas (`eos_methane`, JANAF), with the
wall's heat loss:

    dM/dt = -A* F*(rho, e),     dU/dt = -A* F*(rho, e) h - A_w q.

- **Choked flux:** `F* = sqrt(gamma_s rho p) (2/(gamma_s+1))^(...)` on the local equilibrium
  `gamma_s = rho a^2 / p`. This is within 0.1-0.9% of the exact real-gas sonic flux on the
  isentrope (`near_term.nozzle_exit`), and is calibrated to it.
- **Wall heat:** `near_term`'s bracket along the trajectory. The low edge is the H- radiation floor
  plus frozen-`c_p` Bartz convection; the high edge is an opaque blackbody plus equilibrium-`c_p`
  Bartz.
- **Charge:** fixed at the 212 m^3 chamber's 108.4 kg (ADR-0056). The smaller chambers therefore
  run hotter and denser, as the parent's rows do.

| chamber (5 kg rod; H2 2.5 kg) | T0 | p0 | `gamma_s` | e-fold | 95% out: adiabatic / low / high | to fill pressure | parent |
|---|---|---|---|---|---|---|---|
| CH4 212 m^3 bulged | 7000 K | 92 bar | 1.32 | 65 ms | **221** / 227 / 266 ms | **281** / 284 / 261 ms (1.11 bar) | 215-260 / 278-336 |
| CH4 106 m^3 unbulged | 7457 K | 193 bar | 1.33 | 32 ms | 109 / 111 / 118 | 139 / 140 / 134 (2.23 bar) | 111 / 142 |
| CH4 53 m^3 unbulged | 7976 K | 405 bar | 1.34 | 16 ms | 53 / 54 / 55 | 69 / 69 / 68 (4.45 bar) | 58 / 73 |
| H2 20 m^3 sphere | 5500 K | 818 bar | 1.22 | 33 ms | 110 / 110 / 112 | 143 / 143 / 141 (9.51 bar) | 105 / 135 |

**Wall heat, as a share of the pulse:** 3.3-13.6% (212 m^3), 2.4-7.3% (106), 1.3-3.2% (53),
0.1-1.8% (sphere).
- At the low edge it changes the times by ~1-3%.
- At the high edge a cooled gas flows more slowly at the same pressure, so 95% out comes later
  (+20% at 212 m^3). Its pressure falls faster, so it reaches fill pressure sooner (-7%).
- The heat it takes is debited in B6.

**Recombination:** the chamber gas stays in equilibrium. At fill pressure, H + H + M (Baulch's
`H + H + H2`, half a decade stated) recombines **2.8-3.6 decades** faster than the chamber empties.
So the parent's effective-gamma form is the right physics, and the equilibrium clock lands at its
fast end (it used 4.48 km/s).

**The residue at fill pressure** (adiabatic): ~2.5% of the mass in every methane chamber.

| chamber | residue | T | mole fractions | energy |
|---|---|---|---|---|
| 212 m^3 | 3.1 kg | 4,100 K | H 0.70, H2 0.16, C3 0.04, C 0.04, C2 0.03 | 236 MJ (1.7% of the pulse) |
| 106 m^3 | 3.1 kg | 4,255 K | H 0.68, H2 0.18 | 227 MJ |
| 53 m^3 | 3.0 kg | 4,430 K | H 0.65, H2 0.20 | 220 MJ |
| H2 sphere | 1.6 kg | 3,060 K | H2 0.93, H 0.06 | 66 MJ |

With the high-edge wall heat, the 212 m^3 residue is 6.4 kg at 3,270 K: colder, but more of it.

**The 2-D check** (`sweep --vessel-blowdown`, `blowdown.early_outflow`).
- **The first milliseconds:** the blast throws gas out ahead of the uniform 0-D clock. In 2 ms
  (2.5 kg time; x1.26 at 5 kg) it ejects 8% / 23% / 21% of the mass (212 / 106 / 53 m^3), against
  4% / 9% / 16% for the 0-D clock. That head start is worth a few milliseconds (<2% of the
  blowdown).
- **Not past ~10 ms.** The ghost-cell wall is not conservative over tens of milliseconds:
  - It loses 3-7% of the mass by 20 ms in the rocket chambers, and **a third of it by 30 ms in
    the sphere**.
  - The throat's outflow stalls.
  - A transmissive throat edge also fed reversed lip flow back in, so a 300 ms run filled instead
    of emptying. `Bc::Outflow` (outflow only) fixes that, but not the wall.
  - **The 2-D solver is therefore not used for the whole blowdown.** A conservative cut-cell wall
    is owed before it can be.

## 3. B2: does 2 Hz close?

The cycle: blowdown, films over both openings, plug, fill, lock (`cycle.py`). Mechanism times are
assumptions, stated as such:
- films: 30 ms to lay both;
- fill: 1,000 kg/s per 5 kg of rod, so ~108 ms for 108 kg (the parent sizes 800 kg/s);
- lock: 20 ms (parent).

**The throat film traps whatever is left, and that is hot.** A film cannot be laid across a
sonic jet, so the seal waits for the throat's flow to fade. Whatever is still inside then joins the
next charge. The table gives the mixed charge on two bounds:
- **low:** only the residue's thermal energy reaches the fill;
- **high:** all of it, once its atomic H recombines into the cold gas.

The film must hold the high bound: `t = p r/(2 sigma)` over a hemispherical cap of the throat's
radius (the parent's rule). Kapton yields at 69 MPa; aramid plus a 12.7 µm Kapton liner works at
250 MPa.

| chamber | seal at | seal time | mixed charge | film per shot (aramid / Kapton) | period | fastest rate |
|---|---|---|---|---|---|---|
| CH4 212 m^3 | 1.11 bar (fill) | 281 ms | 600-1,080 K, 1.6-2.9 bar | 1.5 / 5.3 kg | 439 ms | 2.3 Hz |
| | 0.3 bar | 389 ms | 476-660 K, 1.3-1.8 bar | 0.9 / 3.1 kg | 548 ms | 1.8 Hz |
| | 0.1 bar | 488 ms | 440-514 K, 1.2-1.4 bar | 0.7 / 2.4 kg | 647 ms | 1.5 Hz |
| CH4 106 m^3 | 2.23 bar (fill) | 139 ms | 600-1,060 K, 3.3-5.8 bar | 2.9 / 10.3 kg | 297 ms | 3.4 Hz |
| | 0.3 bar | 223 ms | 450-550 K, 2.4-2.9 bar | 1.5 / 5.2 kg | 381 ms | 2.6 Hz |
| CH4 53 m^3 | 4.45 bar (fill) | 69 ms | 610-1,050 K, 6.6-11.4 bar | 5.8 / 20 kg | 227 ms | 4.4 Hz |
| | 0.3 bar | 126 ms | 436-488 K, 4.6-5.2 bar | 2.6 / 9.2 kg | 285 ms | 3.5 Hz |
| H2 20 m^3 sphere | 9.51 bar (fill) | 143 ms | 161-179 K, 20-23 bar | 0.66 / 2.3 kg | 312 ms (500 kg/s) | 3.2 Hz |
| | 0.3 bar | 287 ms | 79-80 K, 9.8 bar | 0.29 / 1.0 kg | 457 ms | 2.2 Hz |

**Margin and what limits the rate:**
- **212 m^3 at 2 Hz:** 61 ms to spare if the film goes on at fill pressure. The blowdown is 64%
  of the cycle.
- **That film then sees 600-1,080 K gas** and holds 1.6-2.9 bar, against the parent's 400 °C seat
  rating and 1.1 bar.
- **Keeping a polymer film under its rating means sealing at ~0.3 bar or below**, and that is
  1.5-1.8 Hz. Two ways to keep 2 Hz:
  - a film that takes ~1,100 K: a metal foil, whose fragments then join the exhaust;
  - a cold purge that flushes the residue out of the still-open throat. It costs ~2-3% of the
    charge per 10 ms at fill pressure (the open-throat leak below).
- **Smaller chambers are limited by the fill and the films, not the blowdown.**
  - 53 m^3 seals in 69 ms, but its 5.8-20 kg film per shot is 5-18% of the charge.
  - The film's mass goes as `p r^3`. The 1.33 m^2 throat's radius is 5.8x that of the parent's
    0.04 m^2 example (~190x in `r^3`), which is why its grams have become kilograms.

**Overlaps:**
- **The fill cannot overlap an open throat.** At the mixed pressure an open throat loses 3-13% of
  the charge per 10 ms (212 -> 53 m^3).
- **The port film and the insert's trim can run during the blowdown tail**, since the port is
  ~1.3% of the throat.
- **The plug can go in during the fill.** It does not need to wait for "the last moment".
  - Fifty milliseconds in the 4,000 K, ~1 bar residue (h ~ 200-1,000 W/m^2 K) costs its
    polyethylene skin 8-40 g, or ≤0.2 kg of the 10 kg of frozen methane.
  - 120 ms in the mixed charge costs ≤60 g.
  - The constraint is the move: 1.5 m rest to rest takes 250 ms at 10 g, 111 ms at 50 g and 78 ms
    at 100 g. So the plug needs ≥50 g to fit within the fill.

**Fastest rates** (seal at fill pressure; the film then needs a ~1,100 K rating):
- 212 m^3: ~2.3 Hz.
- 106 m^3: ~3.4 Hz; ~4 Hz with a 2 t/s fill.
- 53 m^3: ~4.4 Hz.
- H2 sphere: 3.2 Hz at 500 kg/s, and **4 Hz only with a fill ≥ ~1 t/s**: its 250 ms leaves 57 ms
  for 60 kg.

## 4. B3: wall forces per pulse, zone by zone

`zone_loads.py`, on 0.5 cm rod-and-plug solves (`--asks-unbulged` adds 106 and 53 m^3), read at
5 kg and scaled to each chamber's fixed-charge QSP. The 80 m^3 long cone is the existing 40 m^3
run. Zones are the bulge's z-bands. In the unbulged chambers they are bands of the cylinder (and,
at 53 m^3, the cone) at the same z.

| chamber (QSP) | zone | peak | spike impulse | width | D, dry wrap | D, steel |
|---|---|---|---|---|---|---|
| **212 m^3 bulged** (92 bar) | head | 0.51 GPa | 16.6 kPa s | 33 µs | 6.0 | 4.7 |
| | bulge ramp | 0.32 | 12.6 | 25 | 2.3 | 2.2 |
| | plateau | 0.44 | 8.9 | 13 | 1.8 | 1.7 |
| | taper | 0.88 | 17.8 | 5 | 3.6 | 1.9 |
| | cone | 0.15 | 10.1 | 45 | 2.8 | 2.4 |
| **106 m^3 unbulged** (193 bar) | head | 0.90 | 25.1 | 83 | 3.6 | 2.7 |
| | z 1.1-1.8 m | 0.32 | 5.3 | 13 | 1.8 | 1.4 |
| | **z 1.8-3.5 m (beside the plug)** | **2.29** | **45.5** | 8 | **9.9** | 5.0 |
| | z 3.5-6.0 m | 1.61 | 35.2 | 13 | 5.6 | 2.7 |
| | cone | 0.65 | 17.3 | 25 | 4.6 | 3.6 |
| **80 m^3 unbulged** (261 bar) | beside the plug | 2.30 | 49.8 | 8 | 8.1 | 4.0 |
| | cone | 0.79 | 53.2 | 28 | 5.3 | 4.3 |
| **53 m^3 unbulged** (405 bar) | beside the plug | **3.51** | **77.0** | 8 | 7.6 | 2.9 |
| | z 3.5-6.0 m (cone) | 2.70 | 32.9 | 35 | 4.0 | 2.1 |
| | cone near the throat | 2.30 | 22.1 | 8 | 4.0 | 3.3 |

5 kg lengths: z-bands x1.26. The 0.25 cm bulged history gives 0.44-0.67 GPa, 11-17 kPa s and
D 1.8-5.7, and its dry wrap passes at 47.7 t (52.4 t at 0.5 cm).

- **The spike is a line blast beside the plug, `p ~ (E/L)/R^2`.** Unbulged, the wall is back at the
  plug's 1.76 m stand-off whatever the volume, so the band takes 2.3 GPa at both 106 and 80 m^3.
- **Below ~80 m^3 the band also crowds the cone.** At 53 m^3 the plug ends where the cone starts
  (z 3.1 m at 5 kg), and the spike rises to 3.5 GPa and 77 kPa s.
- **The overshoot grows faster than the spike.** D in the band rises 1.8 -> 8-10 because the
  impulse grows 5-8x, while the QSP the wall is sized for only doubles or quadruples.
- **The quasi-static pressure** ends the 1.2 ms solve at 0.5-0.9x QSP beside the plug and up to
  2.2x in the cone, where the gas piles up. It then decays on the B1 clock. The 2-D blowdown's
  first milliseconds show the wall's pressure ringing between 0.2x and 1.9x the mean at 2 ms,
  settling to ±10% by ~20 ms (`zone_quasi_static.csv`).
- **No unbulged chamber's dry wrap passes.** Its flung layers reach 3.3% (106 m^3), 3.8% (80) and
  4.5% (53) beside the plug, against Kevlar's 2.4%.
- **Solid Cr-Mo cracks in hydrogen at once beside the plug** (0 pulses on CC2938) in every
  unbulged chamber. In the bulged taper it lasts 700 pulses (0.5 cm) or 1,400 (0.25 cm), which
  is why the bulged design is the dry wrap.

## 5. B4: a cushioned wall

`cushion.py`. Each station's tile is one degree of freedom on an adiabatic methane cushion with
inner stops and linear orifice damping, as ADR 0054's model, but driven by the **solved** station
history (then the B1 decay). The sweep:
- tile 10 / 25 / 50 mm steel;
- cushion 6 / 8 / 15 / 30 cm;
- precharge 0.64 / 0.85 of the QSP (the parent's 60-80 bar on 94);
- damping 0 / 1e5 / 1e6 Pa s/m;
- tile pitch 0.1 / 0.2 / 0.5 m.

**None of the 390 design-chamber cases (78 designs x 5 chambers) survives.** The best per chamber
(50 mm tile, 30 cm at 0.85 QSP, c 1e5, 0.5 m tiles):

| chamber | tile travel | backwall peak / QSP | backwall D | tile bending | tile spall | energy taken |
|---|---|---|---|---|---|---|
| 212 m^3 bulged | 16 cm | 2.5 | 2.5 | 1.8 GPa | 0.27x | 0.20% |
| 106 m^3 unbulged | 15 cm | 2.1 | 2.9 | 2.2 GPa | 0.83x | 0.38% |
| 80 m^3 unbulged | 28 cm (bottoms at 15 cm) | 26 | 31 | 1.9 GPa | | 0.90% |
| 53 m^3 unbulged | 18 cm | 2.7 | 3.6 | 4.4 GPa | 0.95x | 0.40% |
| H2 106 m^3 bulged (2.5 kg) | 16 cm | 2.3 | 2.4 | 1.9 GPa | 0.42x | 0.30% |

- **The spike is not 2-5 kPa s.** It is 9-18 kPa s on the bulged wall and 46-77 kPa s unbulged
  (§4). A 10 mm tile (78 kg/m^2) leaves the wall at 100-1,000 m/s. The 6-8 cm cushions bottom out,
  and the backwall then takes 15-50x the QSP.
- **A soft cushion stiffens as it is squeezed.** `p_k ∝ (d/(d-x))^1.3`, so a tile that travels
  most of its gap loads the backwall harder than no cushion would. Only 30 cm cushions keep the
  backwall near 2-3x QSP.
- **The binding failure is the tile's own bending.** The solved spike is uneven over a few
  centimetres. A tile on a gas foundation must carry that unevenness in bending, where the dry
  wrap's liner lies on a ~100x stiffer foundation.
  - Free-free plate modes with the cushion as a foundation give 1.8-8 GPa for 10-50 mm tiles.
  - Smaller tiles help only partly, since the spike's structure is finer still. At 0.1 m the
    50 mm tile bends at 0.9 GPa (bulged CH4, ~190 pulses), 1.1 GPa (H2) and 1.8-3.1 GPa
    (unbulged). At that size the tile is only twice its thickness, and plate theory is crude.
  - Bending falls only as 1/thickness. Staying under ~400 MPa takes ~110 mm tiles at 0.1 m
    pitch, or 150-300 mm at 0.5 m. That is 0.9-2.3 t/m^2 over 140-200 m^2, ~130-460 t of tile,
    against the 48 t dry wrap.
  - Plate theory overstates bending when the load varies over less than ~2x the thickness, and
    the spike's centimetre structure is not grid-converged. Even so, the local stress is the
    spike's own 0.9-3.5 GPa.
- **What the cushion does do:** in the unbulged 106 m^3 chamber it cuts the backwall's D from the
  bare wrap's 9.9 to 2.9, and it costs only 0.2-0.9% of the pulse energy (work done on the tiles).
- **What the seams would see** (for the leak discussion):
  - **Pressure across the seal:** up to the spike itself, 0.9 GPa (bulged) to 2.3-3.5 GPa
    (unbulged), for microseconds; then the chamber-to-cushion difference.
  - **Hoop gap:** `x w/R` = 30-55 mm for 0.5 m tiles at 15-18 cm travel; 8-14 mm for 0.1 m tiles.
  - **Step between neighbours along the meridian:** 65-170 mm.
  - **Sliding speed between neighbours:** 30-290 m/s.
  - **Stop impacts on the return:** 30-260 m/s.
- **So a cushion cannot let a 53-106 m^3 chamber survive 4,200 pulses, and it cannot open 4 Hz.**

## 6. B5: the hydrogen chamber

**The wall** (§6a). Real hydrogen's `Gamma_eff` (1.19-1.23 over 20-106 m^3) is within the
solves' gamma 1.2, so the r 3.0 m bulge was solved with the hydrogen charge in it
(`--asks-hydrogen`, 0.5 cm).

**Its cost** (fixed charge: the parent's 60 kg, k = 24):

| chamber | T0 | p0 | throat | to fill pressure | rate | I, AR 100 / 300 | exit, AR 300 |
|---|---|---|---|---|---|---|---|
| sphere 20 m^3 | 5500 K | 818 bar | 0.198 m^2 | 143 ms | 3.2 Hz (4 Hz at ~1 t/s) | 666 / 697 kN s | 8.7 m |
| bulge 106 m^3 | 4884 K | 142 bar | 0.198 | 776 ms | 1.1 Hz | 612 / 647 | 8.7 m |
| | | | 0.50 | 307 ms | 2.1 Hz | 619 / 653 | 13.8 m |
| | | | 0.75 | 205 ms | 2.7 Hz | 622 / 656 | 16.9 m |
| | | | 1.00 | 154 ms | 3.1 Hz | 624 / 658 | 19.5 m |

- **Hydrogen loses more to the bulge than methane does, and mostly in the nozzle.** The bulge
  costs 6.3% / 5.6% impulse (adiabatic, A* 1.0) and 7.8% / 7.0% with the high wall-heat edge.
  - In equilibrium the bulge costs hydrogen little: `eta_blowdown` 0.820 -> 0.798 (AR 100).
  - The freeze costs it more. At 5x lower density, H + H + M cannot keep up, so the stranded
    H-atom store takes 0.08 of `eta`, against the sphere's 0.02.
- **The bulged chamber needs a wide throat to cycle.** On the sphere's 0.198 m^2 throat it takes
  776 ms. At 0.75-1.0 m^2 it reaches 2.7-3.1 Hz.
  - The exit at area ratio 300 then grows to 17-20 m, against 8.7 m.
  - A wider throat also helps the freeze: the debit falls from 0.10 of `eta` (0.198 m^2) to
    0.08 (1.0 m^2).
- **"A smaller rod per pulse" does not help the wall.** By cube-root scaling, a geometrically
  similar chamber sees the same spike. A smaller rod with the same bulge-to-energy ratio is the
  same chamber, smaller: the same pressure and efficiency, only faster on a fixed throat
  (`walled_nozzle_chamber_blast.md` §3i).
- **Verdict:** a bulged hydrogen chamber runs at **~2-3 Hz, not 4**, at -6-8% impulse per rod.
  The ledger's 4 Hz belongs to the unprotected 20 m^3 sphere, whose wall is unresolved at ~2.95 GPa.

### 6a. The hydrogen-filled bulge, solved

2.5 kg rod, 106 m^3, QSP 142 bar (fixed 60 kg charge, 4,884 K):

| zone | peak | spike impulse | D, dry wrap | fling | D, steel | steel in H2 |
|---|---|---|---|---|---|---|
| head | 0.55 GPa | 15.3 kPa s | 4.7 | 1.17% | 3.5 | 2,930 pulses |
| bulge ramp | 0.44 | 14.8 | 1.7 | 1.91% | 1.2 | 18,900 |
| plateau | 0.48 | 16.0 | 2.3 | 1.77% | 1.6 | 20,200 |
| taper | 1.04 | 20.5 | 3.0 | **2.20%** | 1.7 | 226 |
| **cone** | **1.19** | 12.8 | 4.5 | 1.08% | 3.4 | **80** |

- **Beside the plug the line blast is nearly independent of the gas.** In model units the
  hydrogen fill (54 kg) gives 0.7-0.9x the methane fill's (63 kg) peaks at the head, bulge and
  taper. The real hydrogen spike reads ~1.5x methane's only because hydrogen turns more of the
  pulse into pressure: its QSP is 142 bar against 92, and each history is scaled to its own
  QSP. D, which compares the spike with that QSP, is about the same.
- **The cone is the exception.** With less gas to stop it, the remnant carries on into the cone
  and strikes it at 1.19 GPa, against methane's 0.15-0.67.
- **The dry wrap passes at 37.0 t** (sized for the higher QSP). Its margin is fling in the taper,
  2.20% against 2.4% (1.09x), against methane's 1.43x on the same rod. Solid steel cracks in 80
  pulses at the cone.
- **So the bulge and wrap make the hydrogen chamber survive,** with less margin than methane and
  a cone to watch. The 0.25 cm check that confirmed methane's margin is owed for hydrogen.

## 7. B6: impulse per pulse against volume

Net impulse per 5 kg pulse:

    I = m_rod ((1 - 0.01)(1 + k + P) u_e - w),

which includes the head-on debit and the 1% gate. `u_e` is the blowdown's impulse-averaged exhaust
speed less the H + H + M freeze (`eta_net`), on the 1.33 m^2 throat, at the fixed 108 kg charge.
The wall-heat columns debit each parcel's exhaust speed by `sqrt(h0 / h0_adiabatic)`.

| chamber | `eta_blowdown` | `eta_net` | I, adiabatic | I, wall heat low-high | vs 212 m^3 |
|---|---|---|---|---|---|
| AR 100: 212 m^3 | 0.480 | 0.456 | **861 kN s** | 840-761 | |
| 106 m^3 | 0.492 | 0.486 | **901** | 885-848 | +4.6% (+5-11%) |
| 53 m^3 | 0.505 | 0.505 | **925** | 916-902 | +7.4% (+9-18%) |
| AR 300: 212 m^3 | 0.530 | 0.479 | **892** | 871-789 | |
| 106 m^3 | 0.543 | 0.511 | **933** | 917-878 | +4.6% (+5-11%) |
| 53 m^3 | 0.556 | 0.542 | **972** | 963-948 | +9.0% (+11-20%) |

- **A smaller chamber wins three ways:**
  - it is denser, so less of the charge is dissociated;
  - its gas is hotter at the same charge (7,457 / 7,976 K);
  - it has less wall area for a shorter time, so less wall heat.
  The parent's "~7% less impulse per rod than 40 m^3" is the first effect alone.
- **The parent's 906 / 971 kN s are 861 / 892 with the gate and the freeze.** §3h of the chamber-blast
  note reproduces them exactly without both. With the wall heat, 761-840 / 789-871.
- Not counted: pitch in the exhaust (6-15 s of Isp per 2-5 kg, `walled_nozzle_near_term.md`) and
  the throat film's mass (§3).

## 8. Keep the bulge, or bulge and gas bags everywhere?

**Keep the bulge with its dry wrap, and do not add gas bags.**
- **Where the bulge already holds, cushions add nothing.** The dry wrap passes at 48-52 t with D
  1.8-3.6 in the band. The best cushion's backwall D is 2.5, and its tiles still fail in bending
  (0.9-1.8 GPa, ≤190 pulses).
- **Where the bulge is removed, cushions do not save the wall.** They move the failure from the
  backwall to the tile and its seams (§5). The spike's impulse grows 5-8x without the bulge, and
  only stand-off cuts it.
- **The efficiency case for less volume is real, and larger than the parent thought:**
  - +4.6% impulse per rod at 106 m^3 and +7-9% at 53 m^3, adiabatic. With wall heat, +5-20%.
  - Faster cycling: 3.4-4.4 Hz possible against 2.3 Hz.
  - Against that, larger throat films per shot (§3).
- **The cushion's own efficiency cost is small** (0.2-0.9% of the pulse taken by the tiles), and
  it is not the reason to reject it. The reasons are:
  - the tiles' bending;
  - the seams, which see the spike across them and open 1-5 cm;
  - its mass: 50 mm steel tiles are ~390 kg/m^2, ~75 t over the 212 m^3 chamber's 198 m^2 of
    wall before any backwall, against the 48 t dry wrap.
- **The way to the volume savings is a smaller chamber that keeps the stand-off**, not a cushion.
  The r 3.0 m bulge holds 106 m^3 (2.5 kg); its long 12° cone and cylinder hold ~40. Candidates,
  each owed a rod-and-plug solve:
  - a shorter cone behind the bulge. The cone took a late 0.67 GPa hit at 0.25 cm, so shortening
    it must be checked;
  - a bulge only as long as the plug band;
  - a lighter plug that stops sooner, to shorten the band.
  Each could trim 10-30% of the volume, worth ~2-4% impulse per rod.

## 9. Decision: double the throat (2026-10-13)

The 1.33 m^2 throat (`V/A*` 160 m, chamber-blast note §3g) was a timing choice: the narrowest
throat whose blowdown fitted 2 Hz. It was not chosen to recover chemical energy. A wider throat
recovers more, because at a fixed area ratio the nozzle is longer, so H + H + M has longer to keep
up:

| | throat | `eta_net` (AR 300) | freeze debit |
|---|---|---|---|
| CH4 212 m^3 | 1.33 m^2 | 0.479 | 0.051 |
| | **2.66 m^2** | 0.484 | ~0.046 |
| H2 bulge (2.5 kg) | 0.198 m^2 | 0.764 | 0.10 |
| | 1.0 m^2 | 0.784 | 0.08 |

**The cycle at 2 Hz** (methane 212 m^3; the times scale exactly as `V/A*`, since the gas state
depends only on energy per volume):

| throat | period, seal at fill pressure | margin | period, seal at 0.3 bar (polymer film) | margin |
|---|---|---|---|---|
| 1.33 m^2 | 439 ms | 61 ms | 548 ms | **-48 ms** |
| 2.0 m^2 | ~345 ms | ~155 ms | ~417 ms | ~83 ms |
| **2.66 m^2** | 297 ms | **203 ms** | 351 ms | **149 ms** |
| 2.66 m^2 + compact shape (159 m^3; fails the wrap, below) | 263 ms | 237 ms | 313 ms | 187 ms |

- **61 ms was too thin.** The blowdown clock alone is uncertain by ±30-45 ms (the wall-heat
  bracket and a discharge coefficient below 1). Film laying (30 ms) and the fill (1 t/s against
  the parent's 800 kg/s) are assumptions. With a polymer film the margin was negative.
- **The ×2 throat leaves 150-200 ms,** with a polymer film sealed at 0.3 bar.
- **Its costs:**
  - **Throat film:** 2.6 kg of aramid per shot sealed at 0.3 bar, against 0.9 kg (`p r^3`).
    About 2.4% of the charge leaves as film.
  - **The nozzle:** ~32 m exit at area ratio 300, or ~18 m at 100 (against 22.5 / 13.0 m). Not a
    binding constraint (Seth, 2026-10-13).
  - **Impulse per pulse:** 868 / 899 kN s (AR 100 / 300), +0.8% against the 1.33 m^2 throat.
    With the high wall-heat edge, 811 / 840 (+7%), since the gas leaves before it loses its heat.
- **The wall passes** (solved in 2-D with the rod and plug at 0.5 cm, `sweep --vessel-blast
  --compact`, scored by `zone_loads`). Methane, 5 kg, the dry wrap's limit is 2.4% flung hoop
  strain:

| chamber | V | throat | dry wrap | worst fling | peak (zone) |
|---|---|---|---|---|---|
| 212 m^3 bulged (baseline) | 212 m^3 | 1.33 m^2 | 52.4 t, passes | 1.94% (taper) | 0.88 GPa (taper) |
| **throat x2** (12° cone, 1 m shorter) | 210 m^3 | 2.66 m^2 | **46.2 t, passes** | 1.76% (cone) | 0.99 GPa (cone) |
| 2 m taper | 192 m^3 | 1.33 m^2 | 55.8 t, passes | 1.84% (taper) | 0.71 GPa (taper) |
| cylinder 1.5 m shorter | 190 m^3 | 1.33 m^2 | 51.4 t, passes | 1.79% (head) | 0.79 GPa (taper) |
| 30° cone | 196 m^3 | 1.33 m^2 | 48.2 t, **fails** | 2.66% (taper) | 0.71 GPa (taper) |
| all combined (20° cone) | 159 m^3 | 2.66 m^2 | 47.8 t, **fails** | 2.81% (taper) | 1.30 GPa (taper) |

- **So the doubled throat stands.** Its wrap is 6 t lighter at the same QSP (93 against 92 bar),
  mostly because the cone is 1 m shorter. The
  cost is a 0.99 GPa hit on the cone where the baseline had 0.15, much like the hydrogen bulge's
  1.19 GPa (§6a); it stays under spall (0.32 of it).
- **The combined compact shape does not.** Steepening the cone (30° alone, or 20° with the other
  cuts) focuses the reflection back onto the taper, which flings past 2.4%. The 263 / 313 ms row
  of the cycle table above is therefore not available as solved. A wide throat with only the
  milder cuts (2 m taper, shorter cylinder; each passes alone) is the untried next step, worth
  ~10% volume and ~20 ms.
- **Hydrogen does not take the doubled throat.** At 2.5 kg and 1.68 m^2 (`--compact-hydrogen`)
  the dry wrap fails at the taper (2.69% fling, 35.3 t), where the 0.84 m^2 bulge passes at 2.20%
  (§6a); the combined shape fails worse (3.79%, QSP 192 bar). Hydrogen does not need it either:
  its mix is cold (90-190 K), so a polymer film sealed at 0.3 bar works, and the 0.84 m^2 throat
  already cycles at ~2.4 Hz that way (3.0 Hz sealed at fill pressure). Its freeze gain from a
  wider throat (0.764 -> 0.784 at 1.0 m^2) is not worth a failing wrap. The doubling is a methane
  decision.

## 10. Limits

- **Gamma-law 2-D solves.** Real gas enters through the QSP scaling only (as in the chamber-blast
  note). The 0-D blowdown and the impulses are on the real equilibrium gas.
- **Grid.** The 0.5 cm peaks rose 15-25% per halving (the chamber-blast note, §3c). The unbulged
  runs are 0.5 cm only.
- **The 2-D solver's wall is not conservative over tens of milliseconds** (§2): no whole-blowdown
  discharge coefficient, no sloshing past ~10 ms.
- **Wall heat is a bracket.** The high edge (opaque blackbody) is not resolved.
- **The cushion model is one degree of freedom per station.** The tile bending is a modal estimate,
  linear damping stands in for an orifice, and there is no heat into the cushion gas.
- **Mechanism times are assumed:** films 30 ms, fill 1 t/s per 5 kg, lock 20 ms. The film seat's
  own temperature is not modelled. A seat that sees the throat's heat flux recovers slowly unless
  it is a thin, regen-cooled metal ring.
- **Mixed-charge temperature:** the true value lies between the bounds. The residue's H recombines
  quickly; its carbon does not return to methane.
