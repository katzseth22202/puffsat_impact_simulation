# Iapetus water bags on the overtake plate: answers (S50)

Answers 2026-10-21 to the parent's S50 (its `docs/iapetus_bag_plate_asks_for_impact_sim.md`,
commit 5b940a7). **Written to be copied into the parent.** Evidence: Result 20 of
[`xray_collision_study.md`](xray_collision_study.md) and
[ADR-0071](adr/0071-iapetus-bags-ice-charge-paraboloid-and-bullet.md); `make xray-collision-iapetus`
(`iapetus_*.csv`, module `iapetus_bag.py`). **Status: proposed.**

Seth's three questions came with the ask (2026-10-20): should the charge be argon or water, given
that water spends energy on bonds; is the skirted plate enough, or must it be a paraboloid; and what
does the plate weigh at 1 g and 2 Hz. In conversation (2026-10-20/21) he added the bag's miss and
chose a gun that fires an ice bullet back along the bag's track to spread it.

All numbers are for a 25 kg bag at the reference 51.1 km/s unless stated, on the production grid.

## Short answers

1. **Use a water-ice charge, not argon.** Water's bonds do re-form in time; they are not the
   problem. Argon's problem is that it radiates 3-9x more of its heat into the plate's face
   (6-12% of the arriving energy against 0.7-4% for water), and the face has to pay for that in
   eroded mass every pulse. With a graphite face carried as charge, ice gives **7.1-8.2x** and argon
   6.4-7.3x on the paraboloid; at k = 8.5 ice erodes 3-9 kg of graphite a pulse, argon 13-49 kg.
   Ice also needs no cryogenics. Argon wins only if water's chemistry froze on the rebound, which the freeze races
   say it does not.
2. **Fly a paraboloid, not the skirted flat plate.** The delta/D 0.5 paraboloid with the charge at
   two focal lengths gives **7.7-8.4x** for a realistic charge, the skirted flat plate 5.7-6.6x.
   The deeper 0.75 paraboloid adds ~0.5x but rings and rocks more.
3. **The leverage is ~8x, not 10.3x.** `eta_jet` is 0.5-0.6 for a real charge, not 0.775. It
   **does fall with k** for a compact charge (0.53 at k = 8.5 to 0.32 at k = 20 for ice), so the
   optimum is low: **k ~ 4-8.5** for ice, not ~20. A charge spread as a disk holds `eta` better
   (0.62 to 0.56), but then the face heat sets k.
4. **The plate is 25-32 t on the 8 m spring, but it must be graded.** One bag gives 2.9-3.7 MN s,
   which holds a 590-760 t ship at 1 g and 2 Hz. The 63 t plate is twice what that needs (or flies
   a 4 m stroke). The bag's load is centre-heavy, and an even floor puts **47-81% of each kick into
   ringing**. Grade the floor to the load (5-18 cm of steel at the centre, thin at the rim) and the
   ringing falls to 2-8%.
5. **The face needs a renewable layer at its centre.** Bare maraging steel would lose millimetres
   a pulse. The face peak is 0.3-1.75 GPa, inside maraging's 2.3 GPa at 51.1 km/s; at 68 km/s the
   bag must shrink to ~19-24 kg to stay inside it.
6. **Spread the bag with a gun.** A 100 g ice bullet fired back along the bag's track meets it
   ~30 m before the charge. Its miss is the only one the plate survives: **1.03 GPa** on the
   production grid. It needs nothing on the bag. It costs leverage: **6.7x** at k = 8.5 on ice,
   against 7.2x for the compact charge and 8.35x for a perfectly self-spread bag, and 3.4x at
   k = 20, so fly it at k <= 8.5 (lower k was not run with the bullet).

## 1. Arrival state (item 1)

| state | size | what it couples to |
|---|---|---|
| intact liquid | 0.363 m across, 242 kg/m^2 | a compact charge in its path (a column 1.05 m long in solid argon or 1.84 m in ice at k = 8.5) |
| intact ice | 0.373 m, 228 kg/m^2 | the same |
| released ahead, 2.5 mm fragments | covers itself out to **1.6 m** | a charge disk of the same radius; wider, the fragments punch through it |
| released ahead, 1 mm | covers out to 2.5 m | the same |
| released ahead, 10-100 um spray | covers the whole floor | a floor-wide disk |

A fragment of radius `a` meets `3 k sigma / (4 rho a)` of charge per its own mass, which is `k`
times the cloud's optical depth. Fragments that do not cover their cloud's area pass through the
charge, each taking only its own track. At 51 km/s a bag carries 1.3 GJ/kg, a thousand times its
binding and latent heats, so **partly frozen or liquid makes no difference to the coupling.** Only
the geometry does. (A black sphere sits at 278 K at 1 AU, so whether the sunlight melt finishes is a
thermal question this repo does not settle.)

**What does not spread the bag** (hand reasoning; this repo has no fracture model):
- A thin layer, mesh or film, cannot. It is under 1% of the bag's column and shocks only a sliver.
- A laser cannot. It ablates the surface, pushes the bag 4-20 m/s per MJ, and breaks it into
  centimetre spall at best; 1 MJ at ~10 GW per bag is MW-class at 2 Hz.
- Rupturing the bag alone cannot. Its own ~100 kPa drives it apart at ~14 m/s, into centimetre
  drops.
- Hydrogen-oxygen bubbles made from the bag's own water would, but Seth rejected them as too easy
  to set off. Dissolved CO2 with warm water is the safe, Saturn-made alternative if a bag-borne
  spreader is wanted.

## 2. The k sweep (item 2)

The 2-D plate impulse is an effective-gamma Euler run (`sweep --iapetus-2d`, 132 runs, two grids).
Each material's own rebound replaces the lossless one at the arrival speed (`sweep
--iapetus-rebound`: water `e_eff` 0.62-0.67, argon 0.45-0.65). Its gas is matched by the share of
heat it turns into expansion before the plate, on its equilibrium table: **water expands like
gamma 1.17-1.27, argon like 1.30-1.44**, and a gamma bracket of the 2-D run (1.17-1.5) carries that
into the impulse (water 0.3-7% below the gamma = 1.4 run, argon within 1.5%).

**Water's bonds re-form.** In the fireball, `H + OH + M` beats the expansion by 2-3 decades
(Damkohler 70-5000, both edges of the H- and OH-limited bracket). On the rebound the gas stays above
~5500 K while it still pushes, so it is not trying to re-form water yet. What does cost water is
that, in equilibrium, it still holds 14-29% of its heat in bonds at the plate; that is the soft
gamma. Argon's ions recombine with Damkohler 1e5-1e7. The f(v) study's frozen-rebound bracket
(`e_eff` 0.15-0.2 lower) is carried as water's pessimistic edge.

Leverage, kilograms pushed to 10.916 km/s per kilogram of bag (the parent's overtake law with the
measured `B = J/(m w)`), at the best k, charge only:

| plate | point fireball (perfect mixing) | compact charge: argon / ice | 2.5 m disk: argon / ice |
|---|---|---|---|
| flat + 2 m skirt, 2.7 m off | 6.7 / 7.0x | 6.5 / 5.9x | - |
| flat + 2 m skirt, 6 m off | 5.7 / 6.4x | 5.7 / 6.3x | 6.2 / 6.6x |
| bowl 0.30 at focus | 7.4 / 7.8x | 6.8 / 7.3x | 7.4 / 7.9x |
| **paraboloid 0.50 at 2 F** | 8.7 / 8.9x | **7.7 / 8.1x** (k 6.3 / 4.0*) | **8.0 / 8.4x** (k 9.8 / 6.4) |
| paraboloid 0.75 at 2 F | 9.9 / 9.6x | 8.2 / 8.6x | - |

\* k = 4 is the lowest run; for the compact ice column the optimum may lie below it. With water's
frozen-rebound edge the ice columns lose ~0.9x. At k = 8.5 on the 2 F paraboloid the
closing speed moves the ice disk from 8.35x (51.1 km/s) to 10.7x (60) and 12.9x (68).

**Does `eta_jet` fall with k? Yes, for a compact charge.** The charge is a column in the bag's own
footprint, and a longer column is a line source. On the 2 F paraboloid, ice: 0.53 (k = 8.5), 0.32
(k = 20); argon 0.56 and 0.41. A charge spread as a disk holds it: 0.62 and 0.56 (ice). The point
fireball, the ideal, does not fall: argon rises from 0.51 (k = 2) to 0.67 (k = 30), ice holds at
0.62-0.66.

**The charge**: water ice, a column in the bag's path or a disk matched to the spread bag, held at
two focal lengths (4.75 m above the vertex).

## 3. Plate shape (item 3)

| plate | best realistic leverage | exhaust direction | face peak, 25 kg (k 4-20) | rocking lever per m of miss |
|---|---|---|---|---|
| flat + 2 m skirt | 5.7-6.6x | 0.20-0.84 | 0.1-2.0 GPa | 0.99 m |
| bowl 0.30 at focus | 6.8-7.9x | 0.64-0.82 | 0.2-2.3 GPa | 0.66 m |
| **paraboloid 0.50 at 2 F** | **7.7-8.4x** | **0.79-0.90** | 0.25-2.2 GPa | **0.28 m** |
| paraboloid 0.75 at 2 F | 8.2-8.6x | 0.94-0.96 | 0.6-2.7 GPa | 0.75 m |

The direction factor is the momentum leaving along the axis; the rest goes out past the rim. The
flat plate's low end (0.2) is the compact column at 6 m off and k = 20-30, which throws its gas
sideways. The
small spread disks (0.8-1.6 m) peak at 4-21 GPa on the axis on the production grid and are out; a
2.5 m disk peaks at 1.5-1.75 GPa. The paraboloid's heat lands inside 3.5-6 m of the axis.

## 4. Plate mass (item 4)

The lossless spring (sim Result 17): `m s = J tan(theta/4) / (2 theta f)`, 507 t m at 7.35 MN s.

| 2 F paraboloid, ice | J at the push's start / end | ship at 1 g, 2 Hz | plate at 8 m | even floor rings | graded floor rings | graded centre |
|---|---|---|---|---|---|---|
| compact, k = 4 | 3.07 / 2.41 MN s | 625 t | 26.4 t | 55% | 4.7% | 57 mm steel |
| compact, k = 8.5 | 3.35 / 2.63 MN s | 683 t | 28.8 t | 47% | 2.2% | 53 mm |
| 2.5 m disk, k = 8.5 | 3.72 / 2.92 MN s | 758 t | 32.0 t | 80% | 3.7% | 175 mm |

- **The 63.4 t plate is enough, twice over.** At these pulses it would fly a ~4 m stroke.
- **It must be graded, not even.** The centre takes 3-12x the floor's mean impulse per area. An
  even floor leaves 47-81% of each kick ringing, worse than the 38% that ruled out Result 17's
  narrow sheet. Graded to the load, with the parent's 30 kg/m^2 at the rim, 2-8% remains. The
  impulse falls 21% over the push as the closing speed drops, so the grading suits the middle of
  the push; Result 19's k trim holds J as before.
- **Grading costs rocking.** A centre-heavy floor has less moment of inertia, so a given impulse
  offset kicks its rim 1.8-2.7x harder than an even floor's (0.75-1.13 per metre against 0.42).

## 5. Aim tolerance (item 5)

The 2-D code is axisymmetric, so an off-axis burst is a ballistic ray trace: an isotropic point
source drifting onto the plate, each ray's normal momentum returned. Per metre of lateral miss:

| plate | impulse tilt | lateral impulse / axial | lever about the plate's centre of mass | rocking share per cm (even floor) |
|---|---|---|---|---|
| flat + 2 m skirt, 2.7 m off | 2.6 deg | 4.5% | 0.99 m | 0.41% |
| bowl 0.30 | 5.9 deg | 10% | 0.66 m | 0.28% |
| **paraboloid 0.50 at 2 F** | **8.6 deg** | **15%** | **0.28 m** | **0.12%** |
| paraboloid 0.75 at 2 F | 12.5 deg | 22% | 0.75 m | 0.32% |

**A centred burst does not tilt any plate.** Off the axis, the near wall is closer, steeper and
struck more squarely, so it takes more of the push and its push has a sideways part. The
paraboloid's sideways push acts high on its wall, above its centre of mass, and partly cancels the
offset axial impulse's moment, so it rocks a third as much as the flat plate. On a graded floor
multiply the rocking share by ~2.6. At S48's 10 cm aim the paraboloid tilts the impulse ~0.9 deg.

## 6. Face survival (item 6)

The 1-D rebound runs put 0.7-3.8% (water) and 6-12% (argon) of the arriving energy into the wall
as radiation, falling with k as the gas arrives cooler. Against a semi-infinite steel face
(~0.15 MJ/m^2 to melt its surface in a 0.2-0.3 ms pulse), the inner rings take 0.3-50 MJ/m^2 a
pulse. **Bare maraging steel would lose millimetres a pulse, metres a push.**

| 2 F paraboloid | face heat inside 0.6 m | water spray to boil it off | or graphite (20-59 MJ/kg) | leverage, charge only / graphite carried |
|---|---|---|---|---|
| ice compact, k = 8.5 | 1.0 MJ/m^2 | 72 kg/pulse | 3-9 kg | 7.20 / 7.06x |
| ice 2.5 m disk, k = 8.5 | 7.0 | 70 | 3-9 | 8.35 / 8.21x |
| ice 2.5 m disk, k = 15 | 2.6 | 23 | 1-3 | 7.47 / 7.43x |
| argon compact, k = 8.5 | 4.8 | 304 | 13-40 | 7.59 / 7.01x |
| argon 2.5 m disk, k = 8.5 | 34.7 | 381 | 17-49 | 8.01 / 7.28x |

Water spray is a poor ablator (2.6 MJ/kg); carried as charge it would cut argon to 2.7-5.7x. A
graphite liner over the inner 3.5-6 m is the practical face, consumed at 2-18 t a push for ice and
26-100 t for argon. The eroded graphite leaves as propellant. Mechanical erosion and spall are not
modelled.

## 7. The bag's miss, and the gun

Seth's concern: a bag that misses the charge reaches the plate as it came. Face peaks with no
charge in the way (coarse grid; peaks rise 20-45% on refinement):

| what reaches the plate | peak, flat / paraboloid |
|---|---|
| the compact bag, nothing in front | **29 / 40 GPa**: the plate is holed, and ~25-30 GJ goes into the nitrogen spring (hand estimate: ~10 MJ per kg of its gas, ~50x its design pressure) |
| bag hits a 0.1-0.2 bumper 4-6 m ahead | 3.5-9.6 GPa |
| bag spread itself to 1.6 m | 7.7-10.8 GPa |
| bag spread itself to 2.5 m | 3.1-4.6 GPa |
| **100 g bullet ~30 m ahead** | **1.03 GPa (paraboloid, production grid)**; 1.3 GPa for a bullet spread across the bag's face (coarse) |

**The gun** (Seth, 2026-10-21): fire a 100 g ice bullet back along the bag's track when the bag is
~6 km away. They meet ~30 m before the charge, which is where the merged cloud, expanding at
`sqrt(5/3) sqrt(k_b) w / (1 + k_b)` (~3-4 km/s at its edge), reaches the charge disk's 2.5 m.
- Fired along the track, timing only moves the meeting point (+-10 ms is +-2.5 m). Only the lateral
  aim matters: ~0.1-0.15 m of track knowledge at 6 km and ~5 mrad of pointing. From the rim, at
  ~17 deg to the track, the bag's along-track position must be known to ~+-100 m.
- 100 g at 51 km/s carries ~130 MJ, ~1000x what spreading needs. Whether that energy reaches the
  whole ice ball (a solid shatters; the gas model has no strength) barely matters: a bullet
  striking only the bag's middle quarter and one spread across its face give plate impulses within
  1-4% (coarse grid).
- Leverage with the bullet, production grid, 2 F paraboloid: J/mv 2.67 at k = 8.5 and 2.78 at
  k = 20 (lossless), against 2.99 and 3.09 for the compact charge and 3.10 and 3.98 for the ideal
  disk. On ice that is **6.7x at k = 8.5** (compact 7.2x, ideal disk 8.35x) and 3.4x at k = 20.
  The cloud the bullet makes is moving sideways and spread along the axis when it reaches the
  disk, and that motion is partly lost to the push. The bullet's face heat is not read: on these
  long runs the ring tally counts gas recirculating in the bowl (it sums to 1.35-1.43 of the bag's
  energy), so it says nothing.
- Earth-borne cost: ~200 kg of bullets a 2000-pulse push, plus the gun. A 500 g bullet meeting
  ~14 m out gives a few percent more push but its miss reaches 6 GPa.

## Limits

- Effective-gamma 2-D gas, with each material's rebound and expansion matched in. Neither material's
  real EOS is followed in 2-D.
- No strength, fracture or atomization. Every spreading mechanism is a gas cloud in the runs.
- Peaks are not converged (median +20-45% from the coarse grid to the production one); impulses
  are (+2-5%).
- The bumper runs and most bullet runs are on the coarse grid. The bullet's three deciding runs
  are on the production grid; the fourth, the wide bullet's miss, went unstable there (J/mv 174,
  exhaust moving backward) and is dropped by the energy ceiling.
- The face heat is the ring tally of gas moving toward the wall. It counts recirculating gas, so
  it is an upper estimate; it reads 0.15-0.93 of the bag's energy on the runs used.
- The ray trace for aim is ballistic and single-gas.
- The face heat is the 1-D slug's radiative share applied to the 2-D arriving energy.

## Suggested changes for the parent

- `tab:mass_scenarios` and S51: carry the Iapetus push at **~8x** per kilogram of bag (7.1-8.3x
  with the face's graphite charged), not 10.3x, at **k ~ 4-8.5 of water ice** on a delta/D 0.5
  paraboloid at 2 F. The 60 and 68 km/s arrivals give ~10.7x and ~12.9x at k = 8.5.
- The push plate: 25-32 t on the 8 m spring for a 600-760 t ship at 1 g and 2 Hz, its floor graded
  to the bag's centre-heavy load, with a graphite liner over the inner ~6 m consumed at 2-18 t a
  push.
- Spreading and safety: a rim-mounted gun firing 100 g ice bullets to meet each bag ~30 m out,
  at ~6.7x (k = 8.5) rather than ~8x, in exchange for a miss the plate survives.
  "A bag the gun missed is a compact 33 GJ projectile; the plate is not a shield against it."
