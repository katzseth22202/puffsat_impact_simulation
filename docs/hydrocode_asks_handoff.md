# Hydrocode asks of 2026-10-12: answers for the calcs repo

Handoff to `aim_is_all_you_need`, answering its `docs/hydrocode_asks_2026-10-12.md` (B1-B6). It is
written to be copied there. The evidence and methods are in
[`walled_nozzle_hydrocode_asks.md`](walled_nozzle_hydrocode_asks.md); the decisions are
[ADR-0056](adr/0056-hydrocode-asks-blowdown-clock-and-fixed-charge.md) and its 2026-10-13
amendment.

**Status: provisional.** The wall loads come from gamma-law 2-D rod-and-plug solves at 0.5 cm
cells, scaled to each chamber's real-gas pressure; peaks rose 15-25% per grid halving. The clock
and impulses come from a 0-D blowdown on the real equilibrium gas. Mechanism times (films 30 ms,
fill 1 t/s per 5 kg of rod, lock 20 ms) are assumptions.

## Conclusions for the paper

1. **Keep the bulge and the dry Kevlar wrap. Do not add gas cushions.** Without the bulge the spike
   beside the plug grows 5-8x and no dry wrap survives. A cushion's tiles bend at GPa, and none of
   390 cases lasts the 4,200 pulses.
2. **Double the methane throat to 2.66 m^2** (5 kg rod; `V/A*` ~80 m). That makes 2 Hz comfortable
   (150-200 ms to spare), recovers slightly more chemistry, and the wall still passes.
3. **The hydrogen chamber keeps its 0.84 m^2 throat** (2.5 kg rod, 106 m^3 bulge). Doubling it
   fails the wrap, and it already cycles at 2-3 Hz. It does not reach 4 Hz.
4. **Correct `tab:one_chamber`:** 906 / 971 kN s become **861 / 892 kN s** once the 1% gate and
   the H + H + M nozzle freeze are charged (868 / 899 on the doubled throat).
   **Superseded 2026-10-13 by S26** ([`throat_kinetics_handoff.md`](throat_kinetics_handoff.md)):
   the sudden freeze tripled the hydrogen loss; the finite-rate figures are 883 / 938 kN s
   (941 on the doubled throat at A/A* 300).
5. **The throat film is a real cost.** It traps a hot residue that heats the next charge, and it
   weighs kilograms per shot, not grams.
6. **The dry wraps are heavier once the liner's hoop swing is scored (2026-10-15).** Each pulse
   swings the 10 mm Cr-Mo liner up from its -0.75 Y rest, and on the hydrogen crack curve the
   strain-sized wraps failed in 600-1,300 pulses. Sized for 4,200: **hydrogen 43.8 t** (was 37.0),
   **methane throat x2 55.2 t** (was 46.2; owed only if hydrogen reaches its bore, since on air it
   passes already). Fling margins are unchanged or slightly better.

## B1: the blowdown clock

**The parent's closed-form estimate holds, at its fast end.** The 212 m^3 methane chamber
(1.33 m^2 throat) is 95% empty at **221 ms** and at fill pressure at **281 ms** (parent: 215-260
and 278-336 ms).

| chamber | 95% out | to fill pressure | parent |
|---|---|---|---|
| CH4 212 m^3 bulged | 221 ms | 281 ms | 215-260 / 278-336 |
| CH4 106 m^3 unbulged | 109 | 139 | 111 / 142 |
| CH4 53 m^3 unbulged | 53 | 69 | 58 / 73 |
| H2 20 m^3 sphere | 110 | 143 | 105 / 135 |
| **CH4 212 m^3, 2.66 m^2 throat** | ~110 | **139** | |

- **Times scale exactly as `V/A*`** for the same gas state, so doubling the throat halves them.
- **The chamber stays in equilibrium:** H + H + M recombines 2.8-3.6 decades faster than the
  chamber empties. The effective-gamma form is the right physics.
- **Wall heat takes 3-14% of the pulse** (212 m^3) and moves the times by under 10%.
- **A residue is left at fill pressure:** ~3.1 kg at ~4,100 K holding ~236 MJ (mostly atomic H)
  in every methane chamber; 1.6 kg at 3,060 K in the hydrogen sphere.
- The whole blowdown cannot yet be solved in 2-D: the 2-D wall leaks mass over tens of
  milliseconds. The 2-D solve shows only the first milliseconds' head start, which is worth a few
  milliseconds.

## B2: does 2 Hz close?

**On the 1.33 m^2 throat: barely, and not with a polymer film.** On the 2.66 m^2 throat: yes.

- **The film traps the residue.** A film cannot be laid across a sonic jet, so the seal waits for
  the throat's flow to fade. What is left joins the next charge.
  - Sealed at fill pressure, the mixed charge is **600-1,080 K at 1.6-2.9 bar**, against ADR 0053's
    400 °C seat rating and 1.1 bar.
  - Keeping a polymer film under its rating means sealing at ~0.3 bar (476-660 K).
- **Film mass goes as `p r^3`.** At the 1.33 m^2 throat it is 0.9-1.5 kg of aramid per shot (or
  3-5 kg of Kapton), not grams; at 2.66 m^2 sealed at 0.3 bar, **2.6 kg of aramid** (~2.4% of the
  charge).
- **The fill cannot overlap an open throat:** it would lose 3-13% of the charge per 10 ms.
- **The plug can go in during the fill.** It loses ≤60 g of its polyethylene skin, but needs
  ≥50 g of acceleration to move 1.5 m in time.

| methane 212 m^3 | period, sealed at fill pressure | margin at 2 Hz | period, sealed at 0.3 bar | margin |
|---|---|---|---|---|
| 1.33 m^2 throat | 439 ms | 61 ms | 548 ms | **-48 ms** |
| **2.66 m^2 throat** | 297 ms | **203 ms** | 351 ms | **149 ms** |

- **61 ms was too thin:** the blowdown clock alone is uncertain by ±30-45 ms.
- **Smaller chambers are limited by their films and the fill, not the blowdown:** 53 m^3 could run
  4.4 Hz, but its film is 5.8-20 kg per shot.

## B3: wall forces by zone, bulged against unbulged

**The spike beside the plug is a line blast, `p ~ (E/L)/R^2`; stand-off is the only lever.**

| chamber (5 kg) | QSP | worst zone | peak | spike impulse | dry wrap |
|---|---|---|---|---|---|
| 212 m^3 bulged | 92 bar | taper | 0.88 GPa | 17.8 kPa s | **52.4 t, passes** (1.94% fling) |
| 106 m^3 unbulged | 193 bar | beside the plug | 2.29 GPa | 45.5 kPa s | fails (3.3%) |
| 80 m^3 unbulged | 261 bar | beside the plug | 2.30 GPa | 49.8 kPa s | fails (3.8%) |
| 53 m^3 unbulged | 405 bar | beside the plug | 3.51 GPa | 77.0 kPa s | fails (4.5%) |

- **"Fling"** is the hoop strain of the wrap's outer layers flung outward by the spike; Kevlar's
  limit is 2.4%. It, not the peak pressure, is what the dry wrap fails on.
- **Unbulged, the wall sits at the plug's 1.76 m stand-off whatever the volume**, so 106 and
  80 m^3 see the same 2.3 GPa. Below ~80 m^3 the plug band also crowds the cone.
- **Solid Cr-Mo steel cracks in hydrogen beside the plug at once** in every unbulged chamber, and
  in 700 pulses in the bulged taper. That is why the bulged design is the dry wrap.

## B4: a gas-cushioned wall (parent ADR 0054)

**Rejected.** None of 390 design-chamber cases (10-50 mm steel tiles, 6-30 cm cushions,
0.1-0.5 m tiles, 5 chambers) lasts 4,200 pulses; the best lasts ~190.

- **The solved spike is 9-18 kPa s (bulged) to 46-77 kPa s (unbulged)**, not the 2-5 kPa s ADR
  0054 assumed, and it is uneven over a few centimetres.
- **The tiles fail in bending,** at 0.9-8 GPa. A tile on a soft gas foundation must carry the
  spike's unevenness itself. Staying under ~400 MPa takes ~130-460 t of tile, against the 48 t
  dry wrap.
- **Thin cushions bottom out** and the backwall then takes 15-50x the QSP. Only 30 cm cushions
  keep it near 2-3x.
- **What it does do:** in the unbulged 106 m^3 chamber it cuts the backwall's overshoot from 9.9
  to 2.9 times the QSP, for 0.2-0.9% of the pulse energy.
- **Its seams would see** the spike itself across them (0.9-3.5 GPa for microseconds), 1-5 cm hoop
  gaps, and 30-290 m/s sliding between neighbours.

## B5: hydrogen

**The r 3.0 m bulge holds hydrogen, on a thinner margin than methane.**

- **2.5 kg rod, 106 m^3, 0.84 m^2 throat:** dry wrap **37.0 t, passes** at 2.20% fling against
  2.4% (methane: 1.94%). **Revised 2026-10-15: 43.8 t** (2.16% fling), so that the liner's
  hoop-swing crack lasts 4,200 pulses on the hydrogen curve.
  Its 0.25 cm grid check is still running.
- **A new hit on the cone,** 1.19 GPa (methane 0.15-0.67): with less gas to stop it, the blast
  remnant carries on into the cone.
- **Impulse per rod (60 kg charge):** 624 / 658 kN s (AR 100 / 300) against the 20 m^3 sphere's
  666 / 697: the bulge costs **6-8%**, mostly in the nozzle freeze at 5x lower density.
- **Cycle rate: ~2-3 Hz, not 4.** At 0.75-1.0 m^2 it reaches 2.7-3.1 Hz sealed at fill pressure.
  Its mixed charge is cold (90-190 K), so a polymer film sealed at 0.3 bar works (~2.4 Hz). The
  ledger's 4 Hz belongs to the unprotected 20 m^3 sphere, and needs a ≥1 t/s fill even there.
- **Hydrogen does not take a doubled throat:** at 1.68 m^2 the wrap fails at the taper (2.69%).
- **A smaller rod per pulse does not help the wall:** a geometrically similar chamber sees the
  same spike.

## B6: impulse per pulse against volume

**Smaller is better, and more so with wall heat.** Net impulse per 5 kg pulse at the fixed
108.4 kg charge, with the 1% gate and the nozzle freeze:

| chamber | AR 100 | AR 300 | with the high wall-heat edge (AR 300) |
|---|---|---|---|
| 212 m^3, 1.33 m^2 throat | 861 kN s | 892 | 789 |
| **212 m^3, 2.66 m^2 throat** | **868** | **899** | **840** |
| 106 m^3 | 901 | 933 | 878 |
| 53 m^3 | 925 | 972 | 948 |

- **Fixed charge means the temperature follows:** 7,000 / 7,457 / 7,976 K at 212 / 106 / 53 m^3.
- **A smaller chamber wins three ways:** less dissociation, hotter gas, and less wall heat. The
  parent's "~7% less than 40 m^3" is the first alone.
- **The wider throat gains +0.8%,** or +7% at the high wall-heat edge, since the gas leaves before
  it loses its heat. Its recombination efficiency rises slightly (`eta_net` 0.479 -> 0.484): a
  longer nozzle at the same area ratio gives H + H + M more time.
- **But the smaller chambers fail the wall (B3),** so this gain needs a smaller chamber that keeps
  the stand-off.

## Smaller bulged chambers (methane, 2.66 m^2 throat considered)

Each solved in 2-D with the rod and plug (5 kg):

| shape | volume | dry wrap | worst fling |
|---|---|---|---|
| baseline bulge | 212 m^3 | 52.4 t, passes | 1.94% |
| **throat x2** (cone 1 m shorter) | 210 m^3 | **46.2 t, passes** | 1.76% |
| 2 m taper | 192 m^3 | 55.8 t, passes | 1.84% |
| cylinder 1.5 m shorter | 190 m^3 | 51.4 t, passes | 1.79% |
| 30° cone | 196 m^3 | fails | 2.66% |
| all combined (20° cone) | 159 m^3 | fails | 2.81% |

- **Steep cones fail:** they send the reflection back onto the taper.
- **The doubled throat's cost on the wall** is a 0.99 GPa hit on its shorter cone (0.32 of spall).
- **Untried:** the doubled throat with only the 2 m taper and the shorter cylinder, ~10% less
  volume (~+2% impulse, ~20 ms faster).

## Corrections owed in the parent

1. **`tab:one_chamber`:** 906 / 971 kN s -> 861 / 892 kN s (1.33 m^2 throat) or 868 / 899
   (2.66 m^2), with the 1% gate and the H + H + M freeze charged.
2. **ADR 0053 (cycle):** the fill-pressure gate traps a ~4,100 K residue. The film then holds
   600-1,080 K at 1.6-2.9 bar, and weighs kilograms per shot. Widen the throat to 2.66 m^2 and
   seal at ~0.3 bar.
3. **ADR 0054 (cushion):** reject it. Its 2-5 kPa s spike is 9-77 kPa s solved, and the tiles bend
   at GPa.
4. **The hydrogen ledger's 4 Hz:** a bulged hydrogen chamber runs at ~2-3 Hz and costs 6-8%
   impulse per rod against the sphere.

## Still open

- **A conservative 2-D wall** (cut cells), so that the whole blowdown, its discharge coefficient
  and its sloshing can be solved rather than taken from 0-D.
- **A 0.25 cm check of the hydrogen bulge,** whose fling margin is only 1.09x.
- **The mild compact shape** with the doubled throat (above).
- **The seam design** for gas-backed pieces (not started).
