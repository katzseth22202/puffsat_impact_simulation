# ADR-0074: The 1 cm Sun disk couples; fly it on a charge 1.3-1.5x wider, through a matching bore

Status: proposed, 2026-10-27. Result 23 of [`xray_collision_study.md`](../xray_collision_study.md);
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
   disk, the 2.5 kg plate gives 0.95x the thin disk's impulse per kg (3.23 m/s a pulse); a 1 kg
   plate at 1 cm gives ~0.85x (not converged).
2. **Fly the charge 1.3-1.5x the disk's radius.** The same k spread wider is a shorter column. It
   couples 1.13-1.14x the thin baseline (grids agree to 1-2%) and catches an off-centre disk whole
   until the disk's rim reaches the bore. Widen the bore to the charge's radius; at x1.5 that is
   0.285 m and the run already pays for it. At 2x the disk drills a core and coupling falls to 0.86.
   The gain belongs to the long column, not to wide charges in general: the thin disk gains only
   1.08x at x1.3 and nothing at x1.5, so the plate on x1.3-1.5 beats the thin disk on its best
   wider charge by ~6% per kg.
3. **The aim tolerance is the bore's clearance.** On a charge as wide as the disk, f falls 5% at
   1.5-2.8 cm off centre. On the x1.3 charge through the 0.25 m bore, the disk is caught whole to
   6.0 cm; on x1.5 through 0.285 m, to 9.5 cm. The side impulse there is 0.8-1.3% of the axial.
4. **The 5 t bowl is not enough.** Bowl mass on the lossless spring scales as `J / f`. The 2.5 kg
   pulse at 1 Hz is 1.9-2.3 MN s: 33-41 t at an 8 m stroke, or a 52-66 m stroke at 5 t. The
   paper's own 1 kg at 2 Hz already needs ~7 t at 8 m.
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

Convergence of the narrow 1 kg plate and the x1.15 charge (a 0.0125 m grid is `--thick-disk fine`,
running); the thin disk at x1.5 (a 3-cell column); off-centre and split strikes in 3-D; peak loads.

The thin disk on wider charges is `sweep --thick-disk thin-wide [coarse]` (2026-10-27), and
`thick_disk_convergence.csv` lists every case on every grid.
