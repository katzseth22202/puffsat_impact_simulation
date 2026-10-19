# ADR-0068: The gas spring loses ~0.2% a cycle, and the sheet must cover the sprung floor

Status: proposed, 2026-10-19. Result 17 of [`xray_collision_study.md`](../xray_collision_study.md);
code `sweep --plate-profile-2d`, `gas_spring.py` (`make xray-collision-gas-spring`). Answers the
parent's S46 (its `docs/fast_growth_wave_asks_for_impact_sim.md`, commit 5d1379c). Supersedes
ADR-0064's plate rule (its Decision 1).

## Context

ADR-0064 sized the push plate by `I^2 / (2 s M a)`, a plate stopped within each pulse. The parent
found that this plate must shed 2.2 GJ a pulse. It now flies the lossless gas spring of its
`sec:pusher_plate_damping` instead: 63 t on an 8 m stroke, `m s` = 507 t m at 7.35 MN s and
2 Hz. It asks two things. What are the spring's loss and the bags' heat over 2000 cycles? And do
ADR-0064's face and skirt results carry over to the heavier plate?

## Decision

1. **The paper's spring law holds when solved nonlinearly.** An adiabatic column at 63 t and
   `theta = 2 pi / 3` sweeps 7.99 m (8.05 m linear), peaks at 1.22x the mean, and swings its
   pressure x1.75. Its gas column is 18.6 m tall whatever the pressure. A stiffer spring
   (`theta -> pi`) shortens the column to 8.3 m at a 1.61x peak.
2. **The loss is 0.14-0.37% of the 107 MJ swing a cycle**, not 1%: 0.3-0.8 GJ over the push.
   Thermal hysteresis is the largest channel, on the largest of three bracketing heat-transfer
   coefficients. Fibre hysteresis (assumed loss factor 1-5%) and wall friction are smaller.
3. **Nitrogen is the baseline gas.** Its bags warm +65-74 K over 2000 pulses with no cooling,
   inside Vectran's ~150 C. Helium saves 2.4 t of gas but warms +140-152 K bare and needs its
   bundle to radiate. The bag is ~1.5 t of fibre, 2.9 t of nitrogen and 0.1-1.3 t of bladder at
   any pressure. A slow trim (vent ~18%, or ~4 m more column) holds the mean as the gas warms.
4. **ADR-0064's face and skirt results carry over unchanged**, since both are gas-dynamic against
   a rigid wall.
5. **The kick must be rigid, so the sheet covers the sprung floor and the skirt rides the
   frame.** Under the flown ~8 m sheet on the 9.5 m floor, ~38% of the kick's energy rings in
   bending. That is 1.4 GJ a pulse on ADR-0064's 12.2 t plate and 267 MJ on an even 63 t plate,
   at stresses no steel takes, and the spring cannot remove it. A 9.5 m sheet lands evenly to 1%
   with no loss of impulse. At ~70 um of Kapton it holds at u = 8 (634 K). On a 0.3 m sandwich
   floor with the skirt on the craft's frame it rings with 0.1 MJ at ~0.3 GPa. The sheet's
   footprint must hold to within a few percent and err wide.

## Consequences

- The plate is 63 t of sprung floor (a ~0.3 m steel sandwich, ~14 mm faces), plus a skirt on the
  frame and an 18.6 m nitrogen column. The parent's S45 repricing should carry all three.
- Kapton thickness is then set by the floor's area, not chosen: 7.35 MN s / (k v_e) over the
  floor. At 121.8 km/s the 84 um floor needs a plate radius of ~7.9 m or less.
- Open: the plate's back-face heat on the bags, non-axisymmetric aim error (rocking modes), core
  shear, the bags' folding, and a reservoir's duct loss.
