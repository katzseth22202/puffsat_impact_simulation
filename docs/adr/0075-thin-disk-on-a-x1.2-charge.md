# ADR-0075: The thin 1 kg Sun disk flies on a x1.2 charge; widths are read as the grid resolves them

Status: proposed, 2026-10-27. Result 24 of [`xray_collision_study.md`](../xray_collision_study.md);
code `sweep --thick-disk thin-wide [coarse|fine]` and `thick_disk.py` (`on_grid`, `width_rows`,
`side_rows`, `speed_rows`; `make xray-collision-thick-disk`). Answers the parent's S62
([`thin_disk_wide_charge_asks_for_impact_sim.md`](../thin_disk_wide_charge_asks_for_impact_sim.md));
answers doc [`thin_disk_wide_charge_answers_for_parent.md`](../thin_disk_wide_charge_answers_for_parent.md).
Follows ADR-0074 and corrects how its convergence is read.

## Context

aim's ADR 0086 priced the 1 kg, 8 kg/m^2 Sun disk on a x1.3 charge at $129 / $155. It used Result
23's 1.08x from the 0.025 m grid, an assumed 5 cm bore clearance and the thick disk's side-impulse
slope. The parent held the figure until the thin disk was run at 0.0125 m with its clearance and
side impulse traced, at the best width, and at 4 R_sun.

## Decision

1. **A run resolves a charge's radius only to whole cells.** `DiskCharge` fills the cells whose
   centres lie inside the radius and keeps the mass, so a run resolves `round(R / dr) dr`. Every
   width is quoted at the width resolved (`on_grid`, `thin_disk_width.csv`). Convergence is judged
   only between grids that resolved the same width.
2. **The thin disk's best charge is x1.2.** On the 0.0125 m grid it couples 1.074x its matched
   charge. The neighbours are 1.064x at x1.125, 1.047x at x1.31 and 0.896x at x1.5. Like-for-like
   refinements move it -0.8% (x1.125) and +1.0% (x1.5). Result 23's "x1.3, 1.08x" was x1.25 on
   both grids.
3. **x1.2 keeps the 0.25 m bore.** The clearance is 5.05 cm and the disk is caught whole to
   4.0 cm. At the bore's edge it keeps 1.057x or better.
4. **Side impulse is the bowl's: 0.136 per metre of offset on any wide charge.** The tracer gives
   this for both disks and is linear to 3% from 3 cm out. On a matched charge it is half that, since
   the burst sits at the lens's centroid.
5. **4 R_sun changes no ratio.** The 2-D gas has no speed scale. At 616.6 km/s the pulse and side
   impulse scale 1.57x and the peak 2.47x. The bowl on the 8 m spring at 2 Hz is 12.2-13.0 t
   (x1.2). Polyethylene's frozen toll shrinks (velocity factor 0.9925 to 0.9969) and is the same on
   every width.

## Consequences

- aim reruns with 1.074x, 5.05 cm clearance (whole catch to 4.0 cm) and 0.136 / m at both depths.
- ADR-0074's choice of x1.3 for the 2.5 kg plate stands. On the fine grid x1.33 (1.177) beats x1.13
  (1.113) and x1.53 (1.116). Its "converged to 2%" and "x1.5 fell 8%" compared different resolved
  widths, so read them as the fine grid's values.

## Not decided

- **Off-centre strikes in 3-D.** The catch is geometric to 4.0 cm. Past that, the lone-lens bound
  understates the impulse.
- **Absolute peaks.** They are unconverged, as in ADR-0074. Only the v^2 scaling between speeds is
  quoted.
- **Radiation from the dense charge at 617 km/s.** It is outside the gamma gas, as at 392 km/s.
