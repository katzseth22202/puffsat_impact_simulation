# The thin 1 kg disk on a wider charge, and 4 R_sun: answers to S62 (for the parent and aim)

Answered 2026-10-27 in `puffsat_impact_simulation`. **Written to be copied out.** It answers
[`thin_disk_wide_charge_asks_for_impact_sim.md`](thin_disk_wide_charge_asks_for_impact_sim.md)
(a copy of the parent's S62). It is Result 24 of [`xray_collision_study.md`](xray_collision_study.md)
and [ADR-0075](adr/0075-thin-disk-on-a-x1.2-charge.md). Every number below comes from
`python/puffsat/xray_collision/thick_disk.py` (`make xray-collision-thick-disk`). Its tables are
`data/results/xray_collision/thin_disk_{width,offset,side,speed}.csv` and `thick_disk_convergence.csv`.
The new 2-D runs are `sweep --thick-disk thin-wide [coarse|fine]`: the 1 kg, 8 kg/m^2 disk on
x1.15, x1.2, x1.3 and x1.5 charges at 0.05, 0.025 and 0.0125 m.

## The short answer

**aim's inputs nearly hold, but the thin disk's best charge is x1.2, not x1.3.**

- On the 0.0125 m grid the thin disk on x1.2 couples **1.074x** its matched charge. aim assumed
  1.08x. On x1.3 it couples **1.047x**: the 1.08x quoted for x1.3 was really a x1.25 charge (see
  the premise correction below).
- x1.2 keeps the 0.25 m bore, so its clearance is **5.05 cm**, the 5 cm aim assumed. The disk is
  caught whole to 4.0 cm. At the bore's edge it loses at most 1.5%: 1.057x, lone-lens bound.
- The side impulse is linear in offset at **0.136 per metre**, the 13.5% per metre aim borrowed
  from the thick disk. The slope belongs to the bowl, not the disk. At the bore's edge it is 0.67%
  of the axial, 5.9 kN s.
- **At 4 R_sun every ratio and every clearance is unchanged.** The 2-D gas has no speed scale.
  The pulse grows 1.57x, the peak 2.47x, the side impulse 1.57x, and the bowl on the 8 m spring
  12.2-13.0 t.

So aim's rerun should use 1.074x and a 5.05 cm clearance with whole catch to 4.0 cm, not 1.08x and
5 cm. Our reading is that the $155 survives within about a dollar, but that is aim's number to
rerun, not ours.

## A premise correction: the grid resolves a charge only to whole cells

`DiskCharge` fills every cell whose centre lies inside the charge's radius and rescales the
density to keep the mass. The bore's bookkeeping does the same. So a run resolves
`round(R / dr) dr`, not `R` (`thick_disk.on_grid`). Nominal widths a cell apart therefore ran as one
width on the coarser grids:

| thin disk, nominal | 0.05 m | 0.025 m | 0.0125 m |
|---|---|---|---|
| x1.15 | x1.25: 1.082 | x1.125: 1.073 | x1.125: 1.064 |
| x1.2 | x1.25: 1.082 | x1.25: 1.078 | x1.19: **1.074** |
| x1.3 | x1.25: 1.082 | x1.25: 1.078 | x1.31: 1.047 |
| x1.5 | x1.5: 1.041 | x1.5: 0.887 | x1.5: 0.896 |

Each entry is the width the run resolved and its J over the matched charge on the same grid
(`thin_disk_width.csv`). Result 23 quoted the thin disk at x1.3 as "1.08x, grids agree to 0.4%".
Both of those grids ran x1.25. The finest grid, the only one that resolves x1.3, gives 1.047.

Like-for-like pairs converge. x1.125 moved -0.8% from 0.025 to 0.0125 m, and x1.5 moved +1.0%.
The finest grid's curve rises from 1.064 at x1.125 to 1.074 at x1.19, then falls to 1.047 at x1.31
and 0.896 at x1.5. The optimum is x1.2, at a column `L / R` of 0.49.

**This also touches Result 23.** The 2.5 kg plate's x1.3 ran at x1.25, x1.25 and x1.33 on the
three grids, and its x1.5 at x1.5, x1.375 and x1.53. Part of "x1.5 fell 8% on the finest grid" is a
change of width, not of grid. The fine grid still puts x1.3 on top for the plate: 1.113 at x1.13,
1.177 at x1.33, 1.116 at x1.53. So ADR-0074's choice stands. Its "converged to 2%" compared x1.25
with x1.33; read it as the fine grid's value.

## Item by item

### 1. Thin disk, x1.3 charge, finest grid

At 0.0125 m (resolved x1.31) the thin disk couples **1.047x** its matched charge, J/mv 2.237
against 2.138. x1.2 is better, at 1.074x (J/mv 2.296). Taking the matched disk's 1.36 m/s a pulse,
those are 1.42 and 1.46 m/s.

### 2. Bore clearance and side impulse

The bore's clearance is `max(0.25 m, R_c) - R`, with `R` = 0.1995 m. The disk is caught whole while
the offset stays below `R_c - R`.

| charge | R_c | bore | caught whole to | bore clearance | J over matched at the clearance | side impulse there |
|---|---|---|---|---|---|---|
| x1.00 | 0.200 m | 0.25 m | 0 | 5.05 cm | 0.84-0.91 (lens) | 0.36%, 2.5 kN s |
| x1.15 | 0.229 m | 0.25 m | 3.0 cm | 5.05 cm | 1.02 | 0.67%, 5.7 kN s |
| **x1.2** | 0.239 m | 0.25 m | 4.0 cm | **5.05 cm** | **1.057** | **0.67%, 5.9 kN s** |
| x1.3 | 0.259 m | 0.259 m | 6.0 cm | 5.98 cm | 1.047 | 0.81%, 7.1 kN s |
| x1.5 | 0.299 m | 0.299 m | 10.0 cm | 9.97 cm | 0.896 | 1.38%, 10.4 kN s |

Side impulse against offset is traced directly at each offset (`thin_disk_side.csv`, 392 km/s):

| offset | 1 cm | 2 cm | 3 cm | 4 cm | 5 cm | slope |
|---|---|---|---|---|---|---|
| wide charge, lateral / axial | 0.13% | 0.30% | 0.42% | 0.54% | 0.66% | **0.136 / m** |
| x1.2, kN s | 1.14 | 2.68 | 3.78 | 4.87 | 5.88 | |
| matched charge, lateral / axial | 0.07% | 0.13% | 0.21% | 0.30% | 0.36% | 0.072 / m |

- **It is linear.** From 3 cm out every point lies within 3% of the line through the origin. At 1-2
  cm the tracer's noise floor (5e-5 of the axial for a centred burst) gives -6% to +10%.
- On a wide charge the burst sits at the disk's offset. On a matched charge it sits at the struck
  lens's centroid, half the offset, which halves the slope.
- The slope is the same on every wide charge and for both disks (0.136-0.139 / m), as aim assumed.

### 3. The best width

**x1.2.** It couples 1.074x through the unchanged 0.25 m bore, with 5.05 cm clearance and whole
catch to 4.0 cm. x1.3 trades 2.6% of coupling for 0.9 cm more clearance. x1.5 loses 10% on a
7.5 cm column, the disk drilling a core, and it buys 10 cm of clearance through a 0.30 m bore. It
is not worth flying.

### 4. 4 R_sun

The closing speed is 616.6 km/s, twice the 308.3 km/s perihelion speed (`mission.periapsis`). The
2-D runs are in units of the closing speed and the gamma gas has no speed scale. So every coupling
ratio, catch radius, bore clearance and lateral-over-axial fraction above carries over unchanged
(`thin_disk_speed.csv`):

| 1 kg thin disk | 9.86 R_sun | 4 R_sun |
|---|---|---|
| closing | 392 km/s | 616.6 km/s |
| energy | 76.8 GJ | 190.1 GJ |
| matched pulse | 0.838 MN s | 1.318 MN s |
| x1.2 pulse (1.074x) | 0.900 MN s | 1.416 MN s |
| x1.2 side impulse at 5.05 cm | 5.9 kN s | 9.3 kN s |
| peak, over 392 km/s | 1 | 2.47 |
| bowl on the 8 m spring at 2 Hz, x1.2 | 7.8-8.2 t | 12.2-13.0 t |
| polyethylene's frozen toll, velocity factor | 0.9925 | 0.9969 |

- The one speed-dependent loss we price, polyethylene's frozen ions and carbon bonds (`freeze`),
  is smaller at 4 R_sun and the same on every charge width, so it moves no ratio.
- Radiation from the dense charge is not in the 2-D gas at either speed, as in Results 22-23.

## Caveats

- **Peaks are not converged** (Results 22-23). Only the v^2 ratio between speeds is quoted.
- **Off centre is geometry and the tracer.** The axisymmetric code cannot run an offset strike.
  On x1.2 the struck share is exact to 4.0 cm. Past that it is the lone-lens bound, which
  understates the impulse.
- **The design width is nominal.** The finest grid resolved x1.2 as x1.19 and x1.3 as x1.31.
  Clearances are quoted at the nominal radii.

## For the parent's text

Suggested replacement for the sentence in the harder-boost paragraph ("A charge 1.3 times the
disk's radius catches the disk whole and couples more. On inputs borrowed from the thick disk's run
it breaks even at \$129 and \$155, but the thin disk on that charge has not been simulated closely
enough to adopt the figure."):

> A charge 1.2 times the disk's radius couples 7\% more and catches the disk whole to
> \SI{4}{\centi\meter} off center, through the same \SI{0.25}{\meter} bore
> \cite{Katz_puffsat_impact_sim_2026}. [aim's rerun with 1.074x, 5.05 cm clearance and 0.136 per
> metre: break-even here.]

aim should rerun `make disk-train-aim` with these inputs: 1.074x, a 5.05 cm clearance, whole
catch to 4.0 cm with the lone-lens loss to 1.5% at the edge, and 0.136 per metre. Repeat at
616.6 km/s for the $85 / $89 restatement.
