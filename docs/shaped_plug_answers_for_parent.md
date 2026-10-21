# The shaped plug and the gun backup: answers to S52 (for Balloon-Pulse-Propulsion)

Answered 2026-10-21 in `puffsat_impact_simulation`. **Written to be copied into the parent repo.**
It answers the parent's S52 (`docs/shaped_plug_asks_for_impact_sim.md` at `ec5995d`) and follows
our S50 answers (Result 20, ADR-0071). It is Result 21 of
[`xray_collision_study.md`](xray_collision_study.md) and
[ADR-0072](adr/0072-shaped-plug-deep-lens-and-geometric-catch.md). Every number below comes from
`python/puffsat/xray_collision/shaped_plug.py` (`make xray-collision-shaped-plug`), whose tables
are `data/results/xray_collision/shaped_plug_*.csv`.

## The short answer

Numbered as the parent's items.

1. **The lens sweep.** The parent's 1x-deep lens (0.38 m, our share 1.0) fails. The bag punches
   through at every k, and its pancake jets onto the vertex at 2.8-2.9 GPa, over maraging's 2.3.
   **Make the plug twice the bag's depth, about 0.73 m.** At k = 4-15 that recovers the point
   fireball's leverage (8.81x at k = 8.5, 7.65x at the frozen edge, against the compact column's
   7.20x) at a face peak of 0.69-0.90 GPa. Keep it at 2 F; the dish does not help.
2. **Punch-through.** The thinnest plug that stops the bag is 1.5x its depth at k = 3-15. At 1x
   the peak is 2.7-2.9 GPa.
3. **The safe overlap.** Our 2-D code cannot credit a partial overlap: its axisymmetric bound
   shrinks with every refinement of the grid (to ~8% of the bag's face on the finest). So the
   catch radius is the plug's radius less the bag's: 0.04 m at k = 4, 0.14 m at k = 8.5, 0.24 m
   at k = 15. The gun fires whenever the predicted offset at the plug exceeds it. At k = 8.5 the
   backup rate is 5% for a track error at the plug of 5.6 cm and 1% for 4.5 cm. **Result 20's
   0.1-0.15 m was what the gun's aim needs at 6 km, not a track we know.**
4. **The backup pulse.** The plug is spent: the spread cloud knocks it off station at ~98 m/s.
   With the plug in the cloud's middle the pulse gives 2.4-2.6 MN s, about two thirds of a
   normal pulse. At k = 8.5 the plate holds if the bullet meets the bag ~46 m out (1.80 GPa).
   At 31 m it reaches 2.27 GPa with no margin. At k = 4 no tested meeting holds (2.41-3.04 GPa).
5. **The trigger window.** For a bare miss, confirmed: ~6 km out at Result 20's implied 230 m/s
   (31 m meeting). With the lens in the cloud, at k = 8.5, the bullet must meet the bag 46 m out,
   so the gun fires 9.4 km out at 230 m/s. The meeting distance is set by how far the cloud must
   spread, not by the gun. At 1 km/s the gun fires 2.2 km out (0.8 km at 3 km/s), inside the
   plug's last 100 ms of steering, so it sees the steered plug. A 50 g bullet does not help.
6. **The Earth bag plate.** It is 27-39 t graded, with 57-70 mm of steel at the centre, rings with
   0.3-4.7%, and its liner takes 0.4-18 kg of graphite a pulse.
7. **Lower k.** The lens pays only above k ~ 4. Below it the compact column is already short, and
   lens, column and point fireball agree to within 2%. At k = 1-2 no plug depth reliably holds
   2.3 GPa, so no catch can be named there.

## How the runs were done

The flown 9.5 m floor as Result 20's delta/D 0.5 paraboloid, the plug's front face at 2 F
(4.75 m), the 2-D effective-gamma 1.4 code, and Result 20's water corrections (expansion gamma,
rebound, the frozen-rebound edge) applied to every run. The bag is 25 kg of water. The plug is
ice (917 kg/m^3), a cylinder `share` times the bag's 0.363 m depth, its radius fixed by its mass
`k x 25 kg`.

Grids: production dr = 0.05 m, with coarse (0.1 m) checks of every run and fine (0.025 m) checks
of the overlap at k = 8.5. Impulse moves 1-3% between grids. **The face peak does not converge
everywhere.** Halving the cell raised the centred 2x lens's peak from 0.69 to 0.85 GPa (+23%).
The thin lenses' axial jets read 2.2-2.8 GPa on the coarse grid and 7-8 GPa on production. Every
verdict below is read on the production grid, and a meeting distance sets the trigger only if
its peak holds 2.3 GPa with a 1.2x margin.

## 1. The lens sweep

Leverage is kilograms pushed to 10.916 km/s per kilogram of bag (frozen-rebound edge in
brackets). Face peak in GPa beneath it.

| k | compact column | lens 1x (0.36 m) | 1.5x | **2x** | 3x | point fireball |
|---|---|---|---|---|---|---|
| 1 | 7.39 (6.69), 7.05 | 7.34 (6.65), 2.07 | 7.36, 2.80 | 7.01, 5.59 | 7.06, 5.10 | 7.46 (6.75), 1.30 |
| 2 | 7.87 (7.09), 1.84 | 7.87 (7.09), 1.84 | 7.86, 1.29 | 7.62, 2.47 | 7.57, 2.71 | 7.89 (7.11), 0.77 |
| 3 | 8.08 (7.25), 1.12 | 8.21 (7.37), 2.71 | 8.11 (7.28), 1.12 | 8.08 (7.25), 1.12 | 7.44, 1.79 | 8.21 (7.37), 0.59 |
| 4 | 8.12 (7.26), 0.42 | 8.41 (7.52), 2.82 | 8.42 (7.53), 1.74 | **8.18 (7.32), 0.90** | 8.07 (7.21), 0.60 | 8.42 (7.53), 0.50 |
| 8.5 | 7.20 (6.21), 0.35 | 8.16 (7.07), 2.82 | 8.83 (7.67), 1.80 | **8.81 (7.65), 0.69** | 8.37 (7.26), 0.44 | 8.83 (7.67), 0.33 |
| 15 | 4.52 (3.87), 0.33 | 6.38 (5.54), 2.90 | 7.72 (6.76), 2.03 | **8.39 (7.37), 0.69** | 8.12 (7.13), 0.38 | 8.37 (7.35), 0.26 |

Plug radius at 2x: 0.22 m at k = 4, 0.32 m at k = 8.5, 0.42 m at k = 15.

- **Why the 1x lens fails.** A plug no deeper than the bag is flattened by it into a thin hot
  pancake that blows out along the axis, and the jet lands on the vertex. A 2x plug stops the
  bag inside itself, the burst is near a point, and the peak falls 3-4x.
- **2x against 1.5x.** At k = 4-8.5 the 1.5x lens gives the same leverage at a peak of
  1.7-1.8 GPa. With the refinement's +23% that is 2.1-2.2 GPa, too near the limit to fly. At
  k = 15 the 1.5x lens loses 8%. The 2x lens costs at most 3% (k = 4) and holds a 2.5x margin.
- **`eta_jet`** for the 2x lens is 0.63 (k = 4), 0.66 (8.5) and 0.63 (15), against the compact
  column's 0.63, 0.53 and 0.37 and the point fireball's 0.66, 0.66 and 0.63.
- **Exit direction** (share of the exhaust's momentum along the axis) is 0.89-0.91 for every
  plug that holds.
- **Ringing** of the graded floor is 4.7%, 1.1% and 0.3% at k = 4, 8.5 and 15 (an even floor
  rings with 46-57%). Item 6 has the plates.

**Height.** For the 2x lens:

| k | 1.5 F | 2 F | 2.5 F |
|---|---|---|---|
| 4 | 8.17 (7.30), 1.60 GPa | 8.18 (7.32), 0.90 GPa | 8.10 (7.24), 0.46 GPa |
| 8.5 | 8.94 (7.77), 1.39 GPa | 8.81 (7.65), 0.69 GPa | 8.56 (7.42), 0.38 GPa |

Keep 2 F. 1.5 F gains 1.5% at k = 8.5 for twice the peak. 2.5 F halves the peak for 1-3%, which
is the place to go if the peak margin is wanted more than the leverage. The 1x lens fails at
every height (2.5-7.4 GPa).

**The dish** (a concave front face recessed half its depth at the axis) does not help. At 2x it
gives 8.43 (7.54) at k = 4 at 1.0 GPa, a small gain on the flat lens, but 8.35 (7.24) at
k = 8.5 at 1.6 GPa, a loss. At 1x it peaks at 5.9-7.5 GPa: the cup funnels the jet.

## 2. Punch-through

The thinnest plug whose peak holds 2.3 GPa, with every deeper plug also holding, is **1.5x the
bag's depth at k = 3-15**. Just past it, at 1x, the peak is 2.7-2.9 GPa. Thinner still, 0.25-0.7x,
it is 3.8-8.1 GPa.

| k | thinnest stopped | its radius | peak there | peak at 1x | geometric catch at 1.5x |
|---|---|---|---|---|---|
| 3 | 1.5x | 0.22 m | 1.12 GPa | 2.71 GPa | 0.04 m |
| 4 | 1.5x | 0.25 m | 1.74 GPa | 2.82 GPa | 0.07 m |
| 8.5 | 1.5x | 0.37 m | 1.80 GPa | 2.82 GPa | 0.19 m |
| 15 | 1.5x | 0.49 m | 2.03 GPa | 2.90 GPa | 0.31 m |

At k = 1-2 no depth reliably holds. The peak wanders between 1.3 and 5.6 GPa from one depth to
the next as a light plug fails to confine the bag, so neither a widest plug nor a catch radius
can be named there.

The 1.5x plug is the widest that stops the bag, but it holds only on the production grid with
~25% to spare. **We recommend the 2x plug** and size the catch from it. If the parent wants the
extra catch, the 1.5x plug gives 0.05 m more at about the same leverage, at the price of the
margin.

## 3. The safe overlap

The 2-D code is axisymmetric, so an off-centre strike cannot be run. We bracketed it two ways
for the 2x lens at k = 4 and 8.5:
- an **open core**: the plug's middle removed (its density kept), so the bag's middle passes and
  is shocked from all round;
- a **narrow plug**: a plug thinner than the bag, so the bag's rim passes, shocked from inside.

The share of the bag's face that may pass before the peak reaches 2.3 GPa, the lesser bracket:

| k = 8.5, 2x lens | coarse (0.1 m) | production (0.05 m) | fine (0.025 m) |
|---|---|---|---|
| open core | 56% | 19% | 8% |
| narrow plug | 75% | 40% | 27% |
| catch radius | 0.32 m | 0.21 m | 0.18 m |

At k = 4 the production grid gives 23% and 39% (catch 0.11 m).

**It does not converge.** In both brackets the part of the bag that misses the plug is focused
onto the axis by the axisymmetry, as in a shaped charge, and its jet's peak climbs as the cell
shrinks: a 25% open core reads 1.0, 3.2 and 19 GPa on the three grids. A real lateral miss
makes a crescent, which has no axis to collapse onto. So the 2-D brackets are a **lower bound**
on the safe share, and a falling one. They do not size the trigger.

**What to carry.** The geometric catch, with no part of the bag allowed past (the plug's radius
less the bag's 0.18 m), needs no overlap credit. The real crescent's allowance is a 3-D question
(see "What we did not do").

### The catch radius and the backup rate

The gun fires when the track predicts the bag's centre will pass the plug's centre by more than
the catch radius. For a Gaussian track error `sigma` per axis at the plug, the backup rate is
`exp(-c^2 / 2 sigma^2)`. For the 2x lens:

| k | plug radius | catch c | sigma for 5% | sigma for 1% | rate at 3 cm | at 5 cm | at 10 cm |
|---|---|---|---|---|---|---|---|
| 3 | 0.19 m | 0.01 m | 0.3 cm | 0.3 cm | 0.97 | 0.99 | 1.0 |
| 4 | 0.22 m | 0.04 m | 1.5 cm | 1.2 cm | 0.46 | 0.76 | 0.93 |
| 8.5 | 0.32 m | 0.14 m | 5.6 cm | 4.5 cm | 2.8e-5 | 0.023 | 0.39 |
| 15 | 0.42 m | 0.24 m | 9.9 cm | 8.0 cm | 7e-15 | 8e-6 | 0.053 |

(`shaped_plug_catch_geometric.csv`. At k = 1-2 the plug catches nothing: no depth holds.)

The 2-D overlap credit would add at most ~0.04 m at k = 8.5 (0.18 m on the fine grid), 5% at a
sigma of 7.2 cm.

**This is the place where k matters most.** The plug's catch grows with its radius, as
`sqrt(k)`, while the bag's radius is fixed. At k = 4 the plug is barely wider than the bag, so
the gun fires on almost every bag unless the track at the plug is known to ~1 cm. At k = 8.5 the
parent's 5% needs ~5.6 cm. At k = 15, ~10 cm.

**The track at the plug.** Result 20's "~0.1-0.15 m of track knowledge at 6 km" was the
accuracy the gun's lateral aim needs for the bullet to hit the bag. It was not a measured or
assumed track error, and it is not the error at the plug. What sets the backup rate is the error
left after the plug's last ~100 ms of steering, which the parent's track and plug actuation
budgets must supply. We have not sized either.

## 4. The backup pulse and the plug's fate

The gun has fired, and the bullet's cloud reaches the plug. Because the gun fires only on a
predicted miss of a few centimetres past the catch radius, while the cloud is ~2.5 m wide, the
plug sits near the cloud's middle. The runs centre it there, the axisymmetric worst case. The
100 g bullet onto the 2x lens, production grid (coarse in brackets):

| k | bullet meets at | impulse | leverage | face peak | holds 2.3 GPa / 1.2x margin |
|---|---|---|---|---|---|
| 4 | 30.6 m | 2.54 MN s | 5.97x | 2.41 GPa (1.23) | no |
| 4 | 45.9 m | 2.61 MN s | 6.17x | 3.04 GPa (1.01) | no |
| 8.5 | 30.6 m | 2.44 MN s | 4.00x | 2.27 GPa (1.34) | yes / no |
| 8.5 | 45.9 m | 2.50 MN s | 4.16x | 1.80 GPa (1.10) | yes / yes |
| none (Result 20's miss) | 30.6 m | 2.22 MN s | | 1.03 GPa | yes / yes |

- **The plug in the cloud raises the peak.** It is a dense obstacle on the axis that the spread
  cloud jets past. A small plug does it worst: at k = 4 meeting farther out makes the peak higher,
  not lower. The coarse grid hides this entirely (1.0-1.3 GPa everywhere).
- **At k = 8.5, meet at ~46 m.** That holds with margin, and the gun's lead is in item 5. At
  k = 4 no meeting we ran holds, so a backup at low k needs the plug clear of the cloud's middle,
  which the plug's steering cannot do in the gun's lead time. This is a further cost of low k.
- **The plug is spent.** Its 2x lens intercepts 0.8-2.9% of the cloud: 2.5 MJ/kg, under ice's
  ~3.3 MJ/kg to vaporize, but a kick of ~98 m/s, far beyond what its steering can return before
  the next bag half a second later. Count the plug as spent on every backup pulse (it is
  consumed on a normal pulse too). Each backup costs a plug of `k x 25 kg` for a pulse of about
  half the normal leverage (4.16x against 8.81x at k = 8.5).
- **The liner.** A backup pulse sends 8.5 MJ/m^2 into the inner ring at k = 8.5 (Result 20's bare
  miss: 53 MJ/m^2), and erodes up to 28-84 kg of graphite. That is an upper bound, since on these
  long runs the ring tally over-counts the arriving energy. Item 6 carries it a push.

## 5. The trigger window

The gun fires early enough for its bullet to meet the bag far enough up-track that the merged
cloud has spread to the plate. That distance is the bullet's, not the gun's. The shortest
meeting distance whose miss holds 2.3 GPa with the 1.2x margin, on the production grid:

| bullet | meets at | miss's peak | nearer meetings |
|---|---|---|---|
| 100 g | 30.6 m | 1.03 GPa | 23.0 m: 2.10 GPa (no margin); 15.3 m: 3.06 GPa |
| 50 g | 32.5 m | 0.99 GPa | 43.3 m: 1.45 GPa |
| 100 g onto the k = 8.5 lens | 45.9 m | 1.80 GPa | 30.6 m: 2.27 GPa (no margin) |
| 100 g onto the k = 4 lens | none | | 30.6 m: 2.41 GPa; 45.9 m: 3.04 GPa |

**The lens case sets the trigger.** The plug is in the cloud on every backup (item 4), so at
k = 8.5 the bullet must meet the bag ~46 m out, not 31 m.

The gun's lead `t = path / v_bullet + D / w`, the path from the 0.5 bowl's rim to the meeting
point on the axis (27.6 m for 100 g):

| muzzle | fires before arrival | bag's distance then | after the plug's last 100 ms? |
|---|---|---|---|
| 230 m/s (Result 20's ~6 km) | 120 ms | 6.15 km | no |
| 500 m/s | 56 ms | 2.85 km | yes |
| 1000 m/s | 28 ms | 1.44 km | yes |
| 2000 m/s | 14 ms | 0.73 km | yes |
| 3000 m/s | 10 ms | 0.50 km | yes |

With the lens in the cloud at k = 8.5 (45.9 m meeting, a 42.3 m path):

| muzzle | fires before arrival | bag's distance then | after the plug's last 100 ms? |
|---|---|---|---|
| 230 m/s | 185 ms | 9.43 km | no |
| 500 m/s | 85 ms | 4.36 km | yes |
| 1000 m/s | 43 ms | 2.21 km | yes |
| 2000 m/s | 22 ms | 1.13 km | yes |
| 3000 m/s | 15 ms | 0.77 km | yes |

**Confirmed for a bare miss: ~6 km at Result 20's implied 230 m/s. For the lens case it is
9.4 km.** A **faster** bullet fires later, and above
~500 m/s it fires after the plug has finished steering, so the decision sees the steered plug and
the tightest track. A ~1 km/s light-gas or coil gun is the natural choice. A **lighter** bullet
does not help: 50 g needs a slightly longer meeting distance than 100 g (32.5 m against
30.6 m), and its peak is not monotone in the distance. The bullet's own speed (up to 3 km/s on
51.1 km/s) adds up to ~12% to the collision's energy, which the runs (a bullet at rest) do not
include. It would spread the cloud slightly faster.

## 6. The Earth bag plate for the lens

Result 20's spring law and grading, re-run on the 2x lens's load (2 F, ice, the graphite liner
as `erosion_case` charges it, 2000 pulses a push):

| | k = 4 | k = 8.5 | k = 15 | compact, k = 8.5 (Result 20) |
|---|---|---|---|---|
| impulse a bag | 3.09 MN s | 3.86 MN s | 4.49 MN s | 3.35 MN s |
| ship at 1 g, 2 Hz | 629 t | 788 t | 916 t | 683 t |
| graded plate | 26.6 t | 33.3 t | 38.7 t | 28.8 t |
| centre over mean impulse | 4.5x | 3.8x | 3.9x | 3.4x |
| centre steel | 70 mm | 57 mm | 69 mm | 53 mm |
| ringing, graded (even floor) | 4.7% (56.5%) | 1.1% (48.6%) | 0.3% (46.2%) | 2.2% (46.7%) |
| inner ring's arriving energy | 41 MJ/m^2 | 182 MJ/m^2 | 197 MJ/m^2 | 48 MJ/m^2 |
| graphite a pulse | 6.0-17.6 kg | 1.4-4.2 kg | 0.4-1.1 kg | 3.2-9.4 kg |
| graphite a push | 12-35 t | 2.9-8.4 t | 0.8-2.3 t | 6.3-18.7 t |
| leverage, graphite charged | 7.88x | 8.74x | 8.37x | 7.06x |

The plate is **27-39 t**, a little heavier than Result 20's 25-32 t because the lens drives a
larger impulse. The liner's range is graphite's heat of ablation, 20-59 MJ/kg. At k >= 8.5 the
lens consumes less of it than the compact column does.

**On backup pulses**, at the 5% rate, the liner takes up to 2.8-8.4 t more a push at k = 8.5
(up to 28-84 kg a backup pulse, item 4). At k = 4 the backup pulse does not hold the plate.

## 7. Lower k

k = 1, 2 and 3 are in the tables of items 1 and 2, with the compact column and the point fireball
at each.

- **The three shapes converge below k ~ 4.** At k = 2 the compact column is 1.2x the bag's depth,
  already near the lens. Lens, column and point fireball give 7.86-7.89x (7.09-7.11). At k = 3
  they give 8.08-8.21x (7.25-7.37).
- **The leverage curve.** Best lens or column that holds 2.3 GPa: 7.34x (6.65) at k = 1, 7.87x
  (7.09) at 2, 8.11x (7.28) at 3, 8.42x (7.53) at 4, 8.83x (7.67) at 8.5 and 8.39x (7.37) at 15.
  The frozen edge peaks near k = 8.5 but is within 4% from k = 4 to 15.
- **The plate below k ~ 3.** At k = 1 only the 1x lens holds (2.07 GPa, with no margin), and at
  k = 2 the 1x and 1.5x lenses hold (1.84 and 1.29 GPa) while the 2x and 3x do not. The light
  plug cannot confine the bag, and its peak moves between the axis and the rim from one depth to
  the next.
- **Catch is the stronger cost of low k.** A plug of `k x 25 kg` at 2x depth is no wider than the
  bag until k ~ 3, so the gun fires on nearly every bag (item 3). For aim's dollars-per-kilogram
  choice of k, the backup rate and its bullets cost as much as the plug's mass does.

## What we suggest the parent change

- The plug is ice, a cylinder **twice the bag's depth (~0.73 m)** at 2 F. Its radius is
  0.22 m at k = 4, 0.32 m at k = 8.5 and 0.42 m at k = 15.
- **Leverage band for the ledger:** 8.18x (7.32) at k = 4, 8.81x (7.65) at k = 8.5, 8.39x (7.37)
  at k = 15, with graphite charged 7.88x, 8.74x and 8.37x. That is the point fireball's to within
  3%, so the parent's "optimistic ceiling" is now the simulated base. The compact column stays
  as the fallback.
- **Replace the 5% placeholder** with the backup rate from item 4's table at the parent's chosen
  k and its track error at the plug. At k = 8.5 that is 5% at 5.6 cm. **k below ~8.5 makes the
  backup fire often.** That should enter aim's dollars-per-kilogram choice of k alongside the
  plug's mass.
- **A ~1 km/s gun**, firing ~2.2 km out (bullet met 46 m out at k = 8.5) after the plug has
  steered. **Count the plug as spent** on each backup pulse, which gives about half a normal
  pulse's leverage (4.16x at k = 8.5).
- **The plate:** 27-39 t graded, 57-70 mm of steel at the centre, and a graphite liner of 1-35 t a
  push depending on k.

## What we did not do

- **A 3-D off-centre strike.** It is the one run that would credit a partial overlap, and so widen
  the catch. The 2-D bound shows the open part of the bag can focus into a jet, so the crescent's
  allowance needs its own geometry.
- **The plug's steering and the track at the plug.** The backup rate depends on them through
  sigma. We give the rate as a function of sigma.
- **Material strength.** The plug and the bag are fluids here, which is right at 51 km/s for the
  collision. The bag's ice skin, if frozen, and the plug's fracture before the bag arrives are not
  modelled.
- **The bullet's own speed** on the cloud's spread (item 5).
