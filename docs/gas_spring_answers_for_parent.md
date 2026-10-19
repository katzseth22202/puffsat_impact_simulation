# The gas spring's loss, and the plate's ring-down: answers (S46)

Answers 2026-10-19 to the parent's S46 (its `docs/fast_growth_wave_asks_for_impact_sim.md`,
commit 5d1379c). **Written to be copied into the parent.** Evidence: Result 17 of
[`xray_collision_study.md`](xray_collision_study.md) and
[ADR-0068](adr/0068-gas-spring-loss-and-floor-wide-sheet.md); `make xray-collision-gas-spring`
(`gas_spring_*.csv`). **Status: proposed.**

## Short answers

1. **Your spring law holds when solved nonlinearly.** At 63 t, 8 m and `theta = 2 pi / 3` the
   adiabatic column sweeps 7.99 m and peaks at 1.22x the mean. But its gas column is **18.6 m
   tall** (nitrogen) at any pressure, `V/A = gamma Ma / k` at this stroke.
2. **The spring loses 0.14-0.37% of its 107 MJ swing a cycle, not 1%.** That is 0.16-0.40 MJ a
   cycle, or 0.3-0.8 GJ over 2000 pulses. Thermal hysteresis in the gas is the largest part.
3. **The gas: nitrogen.** Its bags warm +65-74 K over the push with no cooling, inside Vectran's
   ~150 C. Helium saves 2.4 t of gas but warms +140-152 K and must radiate. The bag is ~1.5 t of
   fibre, 2.9 t of nitrogen and 0.1-1.3 t of bladder, at any pressure.
4. **ADR-0064's face and skirt results carry over**, since both are gas-dynamic.
5. **But the plate's bending does not.** Under the flown ~8 m sheet on a 9.5 m floor, ~38% of
   each kick's energy rings in bending. That is 267 MJ a pulse on an evenly spread 63 t plate,
   and the spring cannot take it out. **Fly the sheet as wide as the sprung floor, and hang the
   skirt on the craft's frame.** The kick is then rigid: 0.1 MJ rings, at ~0.3 GPa in a 0.3 m
   sandwich floor.

## 1. The cycle

The orbit is shot for an adiabatic column on the plate-craft reduced mass, kicked by `J/m` = 117
m/s each pulse.

| plate | theta | stroke | gas column | peak / mean | pressure swing |
|---|---|---|---|---|---|
| 63 t | pi/2 | 7.7 m | 33.1 m | 1.11 | x1.37 |
| **63 t** | **2 pi/3** | **8.0** | **18.6** | **1.22** | **x1.75** |
| 63 t | 5 pi/6 | 8.4 | 11.9 | 1.37 | x2.43 |
| 63 t | pi | 8.9 | 8.3 | 1.61 | x3.63 |
| 100 t | 2 pi/3 | 5.0 | 11.7 | 1.22 | x1.75 |

The paper's linear stroke (8.05 m) is within 1% of the solve at 2 pi/3. The column is the new
fact. An 18.6 m column will not fit behind a 9.5 m plate as one bag. The options are a stiffer
spring (8.3 m at `theta = pi`, at a 1.61x peak), a heavier plate (11.7 m at 100 t), or a rigid
reservoir beside the bags. The reservoir adds a duct loss that has not been computed.

## 2. The loss per cycle and the bags' heat

The loss has three channels:
- **Thermal hysteresis:** the gas heats as it compresses and loses some of that heat to the wall
  before it re-expands. The gas is integrated cycle by cycle against a 300 K wall, on the largest
  of three heat-transfer coefficients: the laminar oscillating boundary layer, turbulent
  convection at the gas's breathing speed, and Kornhauser & Smith's measured gas-spring
  correlation. The last is extrapolated far past its data, and it is the largest.
- **Fibre hysteresis:** the hoop fibres cycle with the pressure, at an assumed loss factor of
  1-5%.
- **Wall friction** of the breathing flow.

| bag | per cycle | of the swing | push | rise: bare / radiating | gas | fibre | bladder |
|---|---|---|---|---|---|---|---|
| nitrogen, whole floor, 52 kPa | 0.18 MJ | 0.16% | 0.35 GJ | +65 / +38 K | 2.9 t | 1.5 t | 1.3 t |
| nitrogen, 1 MPa, 4 cells | 0.16 | 0.14% | 0.31 | +74 / +62 | 2.9 | 1.5 | 0.13 |
| helium, whole floor | 0.40 | 0.37% | 0.80 | +152 / +73 | 0.5 | 1.7 | 1.5 |
| helium, 1 MPa | 0.27 | 0.25% | 0.53 | +140 / +106 | 0.5 | 1.7 | 0.15 |
| hydrogen, whole floor | 0.24 | 0.22% | 0.48 | +89 / +50 | 0.2 | 1.5 | 1.3 |
| argon, whole floor | 0.27 | 0.25% | 0.54 | +102 / +51 | 4.9 | 1.75 | 1.5 |

- The rise is the bag's lumped fibre, bladder and gas over 1000 s. "Radiating" lets the bundle's
  outside see a 250 K sink at emissivity 0.85.
- Light gases lose more, because they carry heat to the wall faster and store less of it.
- The warming gas raises its own pressure (x1.22 for nitrogen). A slow trim must vent ~18% of the
  gas or let the column ride ~4 m longer.
- Not included: heat from the plate's own back face. If the plate runs hot, it, not the
  hysteresis, will set the bags' temperature.

## 3. Do ADR-0064's face and skirt results carry over?

- **The face: yes.** The peak `0.65 sigma v^2` is gas against a rigid wall. A plate swinging at
  117 m/s (or 600) is at rest to gas closing at 33-100 km/s. The thicker face rings in thickness
  every 3-5 us, shorter than before, so it follows the load at least as well.
- **The skirt's capture: yes.**
- **The plate's bending: no. ADR-0064 never checked it.** The 2-D runs were re-done with the
  impulse tallied every 0.25 m of radius. The flown ~8 m sheet lands 1.4x the mean per area at the
  centre and 0.4-0.8x near the floor's edge. The skirt takes no axial kick at all. A free
  axisymmetric plate (Rayleigh-Ritz) kicked by that profile over the pulse's ~61 us gives:

| plate | sheet | rings per pulse | peak bending |
|---|---|---|---|
| 12.2 t solid, skirt on floor (ADR-0064) | 8.0 m | 38% = 1.4 GJ | ~90 GPa |
| 63 t, even, skirt on floor (solid or 0.3 m sandwich) | 8.0 m | 38% = 267 MJ | 12-13 GPa |
| 63 t sandwich 0.3 m, skirt on the frame | 8.0 m | 12% = 61 MJ | 6.2 GPa |
| same | 9.0 m (-5%) | 1.5% = 6.5 MJ | ~2 GPa |
| same | **9.5 m (= floor)** | **0.1 MJ** | **~0.3 GPa** |
| same | 10.0 m (+5%) | 0 | 0.13 GPa, but 6% of the impulse lost |

The lossless spring carries only the rigid mode. A sheet that lands where the floor's mass is
not leaves a third of each kick ringing in the plate, at stresses no steel takes. The
fix is geometric:
- **A sheet as wide as the floor** lands evenly to 1% and keeps the same impulse per kg
  (2-D: 2.842 against 2.828 `J/(m v)`). At 9.86 R_sun the pulse's 28 kg over the 9.5 m floor is
  ~70 um of Kapton. That holds at u = 8 (634 K face; the floor is 61 um). Where a thicker sheet
  is needed, shrink the floor to fit: at 121.8 km/s, 84 um needs a floor of ~7.9 m or less.
- **The skirt on the craft's frame**, as in the argon plate's P10, so it is not dead mass on the
  floor's rim.
- **The footprint held to a few percent, erring wide.** 5% narrow rings 6.5 MJ at ~2 GPa. 5% wide
  loses 6% and strikes the skirt's lip at 2.5 `sigma v^2`.
- Bags spread under the whole floor bend it not at all. A 5-8 m ring of 1 MPa cells puts
  15-51 MPa in a 0.3-1 m sandwich.

Not modelled: aim error off the axis (it rocks the plate), core shear in the sandwich (the plate
model overstates short-wave stress), and how the bags fold over their swing.

## Suggested changes (for the parent to make)

- `sec:charge_pusher_plume`, after the spring paragraph: "Solved nonlinearly, the spring holds
  the paper's law: 63 t on an 8 m stroke peaks at 1.22 times the mean. Its gas column is 18.6 m
  tall at any pressure, and its bags lose 0.14-0.37% of the plate's swing a cycle. That is
  0.3-0.8 GJ over the push, warming nitrogen bags by about 70 K \cite{Katz_puffsat_impact_sim_2026}."
- Same section: "The spring carries only the plate's rigid motion. The sheet must therefore land
  as the floor's mass lies: as wide as the sprung floor, with the skirt held on the craft's frame.
  Under a sheet narrower than the floor, about a third of each kick rings in the plate
  \cite{Katz_puffsat_impact_sim_2026}."
- The list of unpriced items: drop "Neither is its bags' loss and heat". Add the 18.6 m column
  (or a stiffer spring), and that the floor's radius now sets the sheet's thickness.
