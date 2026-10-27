# A thicker, heavier Sun disk: answers (for the requesting repo and the parent)

Answered 2026-10-27 in `puffsat_impact_simulation`. **Written to be copied out.** It answers
[`thick_disk_impact_request.md`](thick_disk_impact_request.md) (a copy of the request). It is
Result 23 of [`xray_collision_study.md`](xray_collision_study.md) and
[ADR-0074](adr/0074-thick-sun-disk-on-a-wider-charge.md). Every number below comes from
`python/puffsat/xray_collision/thick_disk.py` (`make xray-collision-thick-disk`), whose tables are
`data/results/xray_collision/thick_disk_*.csv`. Grids: 0.05, 0.025 and 0.0125 m. The finest runs
the seven cases the answers rest on, and the bowl and offset numbers use it; other cases are at
0.025 m.

## The short answer

**The 1 cm disk couples; fly the 2.5 kg plate on a charge ~1.3x its radius, and re-price the
bowl.**

- A thick disk does not tunnel. It loses a little to sideways leak along its longer column: 5% at
  2.5 kg on a charge as wide as itself, giving **3.20 m/s a pulse, not 3.4**.
- Spreading the same charge wider shortens the column. At 1.3x the disk's radius it couples
  **1.11x the thin baseline (3.77 m/s)**, converged to 2%. It also catches the disk whole out to
  the 0.25 m bore, **6 cm** off centre. A x1.5 charge through a 0.285 m bore catches to **9.5 cm**
  but couples ~5% less (1.05x).
- **The 5 t bowl cannot take it.** At 1 Hz the pulse is 2.0-2.3 MN s, which needs a 34-43 t bowl
  on an 8 m spring stroke.

## Item by item

### 1. Sweep or punch-through

The snowplow argument carries over to the bowl charge unchanged. A charge as wide as its disk holds
k disk masses in its column whatever the disk's thickness: the column is `k sigma / rho` long,
0.17 m at 8 kg/m^2 and 0.46 m for every 1 cm disk. The 2.75x longer column is exactly the 2.75x
target column the argument asks for, and the charge supplies it automatically. What a longer column
costs is the shocked gas leaking sideways during a longer sweep; the parameter is the column's
length over its radius.

At the same 1 kg (impulse J/(m v); the point fireball is 2.394):

| disk | radius | column L/R | over 8 kg/m^2: 0.05 m | 0.025 m | 0.0125 m |
|---|---|---|---|---|---|
| 8 kg/m^2 | 0.20 m | 0.84 | 1 | 1 | 1 |
| 22 kg/m^2 | 0.12 m | 3.85 | 0.71 | 0.85 | 0.86 |

The thin disk's own J/mv is 1.95, 2.03 and 2.14 on the three grids (the point fireball's is
2.38-2.41), so the ratio is taken grid by grid. The 1 kg plate at 1 cm settles at 0.86 (+1.8% at
the last refinement): 14% worse. That is one more reason to grow the 1 cm plate to 2.5 kg.

### 2. Scaled pulse

At 1 cm, charges as wide as the disk at k = 20:

| disk | grid | over baseline | ship dv | peak (not converged) |
|---|---|---|---|---|
| 1 kg | 0.0125 m | 0.86 | 1.18 m/s | 2.5 GPa |
| 2 kg | 0.025 m | 0.93 | 2.52 m/s | 5.3 GPa |
| 2.5 kg | 0.0125 m | 0.94 | 3.20 m/s | 9.1 GPa |

The 2.5 kg disk carries 192 GJ, the parent's 191 GJ. The fireball's heat is k/(1+k) of it,
183 GJ. Ship dv scales the paper's 1.36 m/s per kg (9.8 km/s from 7.2 t) by the impulse ratio.

**Peak load.** The absolute peaks above are not converged and likely never will be: each is the
highest single-cell pressure near the bore's lip, and it moves 20-50% between grids for every disk,
the thin one included (11.4 GPa at 0.025 m, 7.5 at 0.0125 m). The comparison that holds is against
today's thin disk on the same grid. **On the recommended x1.3 charge, the 2.5 kg plate peaks at
0.84x the thin disk's peak at 0.025 m and 0.87x at 0.0125 m**, at the same place (0.26-0.29 m, the
bore's lip). So the bowl's liner sees no worse a peak than it does now, for 2.8x the impulse. Other
cases' ratios swing between grids (the matched 2.5 kg plate goes 0.71x to 1.20x), so only x1.3's is
quoted.

**The bowl.** On the lossless spring,
`m s = J tan(theta/4) / (2 theta f)` (our Result 17), at theta 2pi/3-5pi/6:

| pulse train | J | stroke at 5 t | bowl at an 8 m stroke |
|---|---|---|---|
| 1 kg at 2 Hz (the paper) | 0.84 MN s | 11.6-12.3 m | 7.2-7.7 t |
| 2.5 kg at 1 Hz | 1.97 MN s | 54-58 m | 34-36 t |
| 2.5 kg at 1 Hz, charge x1.3 | 2.32 MN s | 64-68 m | 40-43 t |

Impulses are from the 0.0125 m grid; a 5 t bowl takes a 390 m/s kick from each 2.5 kg pulse.

Bowl mass goes as `J / f`. The 2.5 kg disk at 1 Hz is ~4.7x the paper's pulse train; at 2 Hz it
would be half that, at twice the mean force. Even the paper's 5 t is short at an 8 m stroke.

### 3. Off-centre tolerance

The 2-D code is axisymmetric, so this is geometry plus the ballistic tracer. The part of the disk
that misses the charge flies on, away from the plate, and is lost. The struck share is bracketed
between the struck lens alone and the lens joined by the unstruck charge. The off-axis fireball's
side impulse is traced on the paraboloid.

On a charge as wide as the 2.5 kg disk (R_c = 0.19 m):

| f drop | offset | over R_c | side impulse |
|---|---|---|---|
| 5% | 1.5-2.8 cm | 0.08-0.15 | 0.1% of axial, 2 kN s |
| 20% | 6.0-10.6 cm | 0.32-0.56 | 0.4%, 7 kN s |

The disk clips the 0.25 m bore at 6.0 cm whatever the charge, so the 20% line is the bore's edge.

**A wider charge.** The same k spread over a larger radius, the disk centred (2-D runs; 0.0125 m
except x2):

| charge | bore | centred, over the matched charge | over the thin baseline | caught whole to | side impulse there |
|---|---|---|---|---|---|
| x1.15 | 0.25 m | 1.11 | 1.05 | 6.0 cm | 0.8%, 16 kN s |
| x1.3 | 0.25 m | **1.18** | **1.11** | 6.0 cm | 0.8%, 19 kN s |
| x1.5 | 0.285 m | 1.12 | 1.05 | 9.5 cm | 1.3%, 29 kN s |
| x2 (0.025 m) | 0.38 m | 0.91 | 0.86 | -- | -- |

x1.5 led on the two coarser grids (1.14-1.15) and fell 8% on the finest; x1.3 moved 2%. So x1.3 is
the converged best, and x1.5 trades ~5% of coupling for 3.5 cm more aim.

**Does the thin disk gain as much?** No. Spread to x1.3 (0.26 m) it gains 1.08x over its own
matched charge (grids agree to 0.4%). At x1.5 its column is only 7.5 cm and it falls to 0.89 (not
converged; read it as no gain). The wider charge pays most where the matched column is long, which
is the 1 cm plate: 1.18x at x1.3. Per kg, the 2.5 kg plate on x1.3 beats the thin disk on its
best wider charge by 3-5% (1.11-1.13 against 1.08; the thin case was not run at 0.0125 m).

There the tolerance is the bore's clearance, not the overlap. For the aim repo, the guidance
requirement is the disk's rim clearing the bore: 6 cm through 0.25 m on x1.3, or 9.5 cm through
0.285 m on x1.5 at ~5% less coupling, at about 1% side impulse.

### 4. A split disk

- **Halves a few cm apart along the axis** arrive within 1% of the column's sweep time, so they
  couple as one body.
- **Halves pushed apart across the crack** overhang a charge as wide as the disk at their rims:
  `2 g / (pi R)` of the disk misses, 3.3% at 1 cm and 10% at 3 cm (impulse 0.90-0.95x at 3 cm).
- **A x1.15 charge catches a 5 cm gap whole**, and x1.3 catches 10 cm. On the recommended wider
  charge, the pulse does not care.

### 5. Alumina coat

100 um of alumina is 1.8% of the 22 kg/m^2 disk and 0.085% of the pulse's mass (2e-4 of its atoms).
Bounds:
- **Energy:** if its atomization and first ionization were all stranded, that is <= 2e-5 of the
  fireball's heat.
- **Radiation:** argon's 15% loss on the Sun bowl, scaled by alumina's atom share with a 3x
  allowance per ion, is <= 1e-4.

It does not change f measurably.

## Caveats

- **Convergence** (`thick_disk_convergence.csv`, change at the last refinement): 2.5 kg matched
  -0.9%, 1 kg plate +1.8%, x1.3 -2.0%, x1.15 -2.8%, x1.5 -8.0%. x1.3 may settle a point or two
  lower (~1.08-1.10x), still well above a matched charge's 0.94x. x1.5 may fall further, which only
  raises the price of its extra aim. The 2 kg disk, the x2 charge and the thin disk on wide charges
  have only 0.05 and 0.025 m runs; none of them is the proposed design, and x2's loss (0.86x) is
  far outside any drift seen.
- **The gas has no strength.** The 1 cm plate is a 4-cell gas slab of the right mass; at 392 km/s
  that is the right limit for the collision, not for launch or dust.
- **Peaks are not converged** (Result 22).
