# ADR-0074: The 1 cm Sun disk couples; fly it on a charge ~1.3x wider

Status: proposed, 2026-10-27. **Amended by [ADR-0075](0075-thin-disk-on-a-x1.2-charge.md)
(2026-10-27):** the runs resolve a charge's radius only to whole cells. This ADR's x1.3 ran as
x1.25 / x1.25 / x1.33 and its x1.5 as x1.5 / x1.375 / x1.53, so "converged to 2%" and "x1.5 fell
8%" compare different widths. The fine grid alone still puts x1.3 on top for the 2.5 kg plate.
The thin disk's "1.08x at x1.3" was x1.25; at x1.3 it is 1.047x and its best is x1.2 (1.074x). Result 23 of [`xray_collision_study.md`](../xray_collision_study.md);
code `sweep --thick-disk [coarse|fine] [thin-wide]`, `euler2d::bounce::DiskCharge::charge_radius`, and
`thick_disk.py` (`make xray-collision-thick-disk`). Answers
[`thick_disk_impact_request.md`](../thick_disk_impact_request.md); answers doc
[`thick_disk_answers.md`](../thick_disk_answers.md). Builds on ADR-0060 (the real collision),
ADR-0068 (the lossless spring), ADR-0071 (the ballistic aim tracer) and ADR-0073 (density is column
length; the Sun bowl's disk enters through the bore).

## Context

The Sun boost flies 1 kg graphite disks, 40 cm across and 3.6 mm thick (8 kg/m^2), head-on at
~392 km/s into k = 20 polyethylene charges in a bowl. The request proposes a 1 cm molded plate at
22 kg/m^2 (cheaper, tougher, more dust-tolerant), ideally 2.5 kg so that it keeps the parent's
40 cm footprint and per-pulse energy, flown at 1 Hz. It asks whether a thick disk still couples or
tunnels, what the scaled pulse does to the 5 t bowl, how far off centre it may strike, whether a
cracked disk matters, and whether an alumina coat matters.

## Decision

1. **A thick disk does not tunnel.** A charge as wide as its disk holds k disk masses in its
   column whatever the disk's thickness. The column is `k sigma / rho` long: 0.46 m for every 1 cm
   disk. Its cost is sideways leak along a longer sweep, set by `L / R`. On a charge as wide as the
   disk, the 2.5 kg plate gives 0.94x the thin disk's impulse per kg (3.20 m/s a pulse); a 1 kg
   plate at 1 cm gives 0.86x. Both are converged to 2% on the 0.0125 m grid.
2. **Fly the charge ~1.3x the disk's radius.** The same k spread wider is a shorter column. At
   x1.3 it couples 1.11x the thin baseline (1.18x over its own matched charge; -2% at the last
   refinement). x1.5 led on the coarser grids but fell 8% on the finest, to 1.05x, and x1.15 sits at
   1.05x. At 2x the disk drills a core and coupling falls to 0.86. The gain belongs to the long
   column, not to wide charges in general: the thin disk gains only 1.08x at x1.3 and nothing at
   x1.5, so the plate on x1.3 beats the thin disk on its best wider charge by 3-5% per kg.
3. **The aim tolerance is the bore's clearance.** On a charge as wide as the disk, f falls 5% at
   1.5-2.8 cm off centre. On the x1.3 charge through the 0.25 m bore, the disk is caught whole to
   6.0 cm. Widening to x1.5 and a 0.285 m bore buys 9.5 cm for ~5% of coupling; that is the aim
   repo's trade. The side impulse at the bore's edge is 0.8-1.3% of the axial.
4. **The 5 t bowl is not enough.** Bowl mass on the lossless spring scales as `J / f`. The 2.5 kg
   pulse at 1 Hz is 2.0-2.3 MN s: 34-43 t at an 8 m stroke, or a 54-68 m stroke at 5 t. The
   paper's own 1 kg at 2 Hz already needs 7-8 t at 8 m.
5. **A cracked disk and an alumina coat do not matter** on a charge of x1.15 or wider. Halves a
   few cm apart along the axis couple as one; across the crack, a x1.15 charge catches a 5 cm gap
   whole. 100 um of alumina costs <= 1e-4.

## Code

`DiskCharge` takes `charge_radius` (0 keeps the disk's own; every earlier call passes 0). Its
initializer lowers the ambient's pressure only when the fixed `1e-6 p_cold` would fail the cold
check, which a narrow disk on a fine grid can; every start that passed before is unchanged, and the
run reproduces Result 22's baseline exactly. `iapetus_bag.trace` takes `away` for a burst drifting
off the plate (default unchanged).

## Consequences

- The disk becomes a 1 cm, 2.5 kg molded-graphite plate if the parent accepts a wider charge and
  bore and a heavier bowl. C/C drops out.
- The aim repo should take the bore's clearance (6-9.5 cm) as the 2.5 kg disk's catch radius, not
  the overlap curve of a charge as wide as the disk.
- The bowl's mass, or the cadence, must be re-priced: 2.5 kg at 2 Hz would halve the bowl and
  double the mean force.

## Not decided

None of these bears on the decision:

- **Whether x1.3 and x1.5 settle further.** x1.3 moved -2% at the last refinement and may land at
  ~1.08-1.10x, still well above a matched charge's 0.94x. If x1.5 falls below 1.05x, only the price
  of its extra aim rises; that is the aim repo's trade.
- **The 2 kg disk, the x2 charge and the thin disk on wide charges at 0.0125 m.** None is the
  proposed design. x2's loss (0.86x) is far outside any drift seen, and the thin disk feeds only the
  side comparison (the plate beats it by 3-5% per kg), not the reason for the plate.
- **Off-centre and split strikes in 3-D.** On x1.3 the catch is geometric: the whole disk lands on
  charge until its rim meets the bore. 3-D would sharpen the aim repo's number, not the choice.
- **Absolute peak loads.** The single-cell peak at the bore's lip moves 20-50% between grids for
  every disk. On x1.3 the plate's peak is 0.84-0.87x the thin disk's on the same grid, on both
  fine grids, so the liner sees no worse a peak than now. The bowl is sized by impulse on the
  spring, not by peak.

The bowl's verdict depends only on the ~2 MN s pulse, so it holds whatever these give.

The thin disk on wider charges is `sweep --thick-disk thin-wide [coarse]`, and the 0.0125 m grid is
`sweep --thick-disk fine` (7 h, run by hand; 2026-10-27). `thick_disk_convergence.csv` lists every
case on every grid; the bowl and offset tables use the finest grid present.
