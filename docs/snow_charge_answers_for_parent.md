# Water or argon snow as the charge: answers to S60 (for Balloon-Pulse-Propulsion)

Answered 2026-10-25 in `puffsat_impact_simulation`. **Written to be copied into the parent repo.**
It answers the parent's S60 (`docs/snow_charge_asks_for_impact_sim.md` at `860b570`, item 7 added
there) and folds in Seth's follow-ups of the same session: the 100 km/s Earth runs, the 3S cycle's
44-55 km/s Jupiter return, and his design point of **k ≈ 10 at 44 km/s**. It is Result 22 of
[`xray_collision_study.md`](xray_collision_study.md) and
[ADR-0073](adr/0073-snow-charges-density-is-column-length.md). Every number below comes from
`python/puffsat/xray_collision/snow_charge.py` (`make xray-collision-snow-charge`), whose tables are
`data/results/xray_collision/snow_charge_*.csv`.

## The short answer

**Snow works where its column stays short, water couples like polyethylene, pores do not matter at
your plates' speeds, and argon is worth it only at low speed.**

- **The premise needs one correction.** Good coupling does not come from density contrast or from
  pores. At k = 20, k/(1+k) = 95% of the impactor's energy becomes heat whatever the charge is
  made of. What a lighter charge changes is the length of the column the impactor crosses, and a
  longer column is a worse source on the plate.
- **Water snow at 600 kg/m^3 is within -7% to +4% of solid polyethylene on every plate.** Ice and
  liquid water are the same as solid polyethylene to within the grid.
- **Your Sun-bowl grading rule is inverted.** On the Sun bowl the impactor enters through the bore,
  from the plate side, as on the departure bowl. So "light where it enters" puts the dense end
  toward the mouth on both bowls, not toward the floor.

## Item by item

### 1. Density on each plate

Impulse over the same plate with solid polyethylene (950 kg/m^3), production grid:

| bulk density [kg/m^3] | Sun bowl, k = 20 | departure bowl, k = 12.4 | growth push, k = 6.4 |
|---|---|---|---|
| 1000 (liquid water) | 1.000 | 1.003 | 1 |
| 917 (ice) | 1.000 | 1.000 | 1 |
| 600 (snow) | 0.981 | 0.981 | 1 |
| 500 (froth) | 0.975 | 0.966 | 1 |
| 300 (snow) | 0.969 | 0.909 | 1 |
| 100 | 0.862 | 0.621 | |
| 30 (the old foam) | 0.533 | -- | |
| **within 2% down to** | **~580** | **~590** | any |
| within 5% down to | ~250 | ~430 | |

- **Breakpoints.** The Sun bowl holds within 2% to ~580 kg/m^3 and within 5% to ~250. The departure
  holds within 2% to ~590 and within 5% to ~430.
- **Why the departure is more sensitive.** Its solid column (0.39 m) is already about as long as
  the 0.4 m impactor is wide. The Sun's (0.11 m) is less than half the 0.25 m disk's radius, so it
  has room to triple.
- **The push does not care.** Sweeping the sheet's grid thickness over 0.1-0.8 m moves the impulse
  by 1.9%. The real sheet is 0.5-2.5 mm thick at any of these densities.
- **The real collision against the point fireball** is 0.88 at k = 20 on the Sun bowl and 0.95 at
  k = 12.4 on the departure, measured on solid polyethylene. Our ledger had held 0.82-0.87 down
  from k = 30.
- **Snow also cuts the on-axis back-splash peak.** At the departure the peak is 0.56x solid's at
  600 kg/m^3 and 0.25x at 300; at the Sun, 0.57x at 300. Peaks are not converged on our grids (they
  move 18-30% per halving), so take this as a trend.

### 2. Material at fixed density, and porosity

**Is Seth's "like polyethylene foam" right?** Yes, for the hydrodynamics. Water and polyethylene at
the same bulk density make the same column, and the 2-D gas cannot tell them apart. They differ only
in chemistry and rebound, raced on each material's own equilibrium table at each plate's heat. Over
polyethylene per kg of impactor, with polyethylene's carbon either frozen atomic or condensed:

| plate (heat per kg) | water | argon | why |
|---|---|---|---|
| push (573 MJ/kg) | 1.01 / 0.98 | 0.96 / 0.93 | argon's gamma 1.33 and rebound 0.61 |
| departure (244 MJ/kg) | 1.06 / 0.95 | 1.02 / 0.91 | water's bonds hold 13% at the plate, none stranded |
| Sun (3.5 GJ/kg) | 1.02 / 1.01 | **0.85 / 0.84** | argon's rebound falls to 0.27: it radiates into the face |

**Bond toll.** The parent's 19-30% figure takes every bond as frozen. Water's bonds are not frozen:
no share strands, since H + OH + M beats the expansion at every plate. So water's real price is
softness. While it is still dissociated near the plate it expands with an effective gamma of
1.27-1.40, and that is already inside the numbers above. Polyethylene's carbon is bracketed, as in
the parent's ledger.

**Porosity.**
- **The test.** In the snowplow the last kilogram is struck at `w / (1 + k)`, and its shock heat
  `V^2/2` is exact by conservation whatever the EOS. Pore collapse can add at most what a solid
  holds as cold compression. That share is bounded with linear Hugoniots and a Mie-Grüneisen release
  (with `s` +/- 0.2 and Gamma 0.5-2), at that slowest strike only, since the hotter states are past
  the fits.

| plate | last strike | heat over vaporization | cold share (max) |
|---|---|---|---|
| push | 13.4 km/s | 25x (PE), 29x (water) | 0.10-0.21 |
| departure | 6.3 km/s, inside the fits | 5.5x (PE), 6.4x (water) | 0.19-0.35 |
| Sun | 18.7 km/s | 48x (PE), 57x (water) | 0.07-0.16 |
| argon, every plate | | 76x or more | |

- **The verdict.** Every kilogram of the column vaporizes, porous or not. Pore heating adds to
  heat already in surplus, so it changes neither the vaporized share nor the impulse.
- **The 2-D gas is already the porous limit.** It counts all dissipated energy as heat, so a
  P-alpha model is not needed at these plates.
- **Where this breaks.** At 44 km/s and k = 20, water's last kilogram gets only 0.7x its
  vaporization energy. See the low-speed section.

### 3. Argon

- **No bonds, and its ions do not freeze.** The frozen share is about 0.1% of the heat, with
  Damköhler numbers of 3e3-1e7, against 39-42% of the heat that starts in ionization. So its frozen
  toll is near zero.
- **Its cost is the plate, not the toll.** Its effective gamma is 1.30-1.33, softer than you might
  expect from a monatomic gas, because the ionization store soaks up heat as it expands. Its
  rebound is weak because hot argon radiates into the face: 0.59-0.61 at the push and departure
  speeds, 0.27 at the Sun bowl's 83 km/s.
- **Net:** level at the departure, 4-7% behind on the push, 15% behind at the Sun.
- **Your premise on argon snow needs a correction.** A sealed bag cannot hold 22% of its argon as
  vapour at 60 K:
  - Flashing liquid argon gives 18% vapour at the triple point and 25% on cooling to 60 K.
  - At 60 K the frost holds only 904 Pa (0.072 kg/m^3), so 22% as vapour would need 3.0 m^3 per kg
    of charge.
  - The vapour must deposit as frost, which means the bag must radiate the flash heat (42 kJ/kg).
    Shaded at 70 K that takes ~31 days for the 20 kg Sun charge and ~64 days for the 190 kg
    departure charge, longer than the 10-day parking orbit.
  - Until then the bag is a near-vapour line source.

### 4. The bag

**On the bowls it is negligible:** 0.05-0.2% of the charge.

**On the push it is not.** Film of 0.05 kg/m^2 on both faces of a 0.45 kg/m^2 charge is 22% of it.
That film is polyethylene debris and pays polyethylene's toll. At 10 um it is 4%.

**The gas held is negligible:** water vapour or argon at the frost's pressure is under 2e-4 of the
charge.

### 5. The growth push's geometry

- **Thickness does not matter.** The 1.9% spread over 0.1-0.8 m of grid thickness means a 1-2.5 mm
  pillow couples as the 0.5 mm sheet does.
- **What does matter is the film, not the shape.** A different charge shape gains nothing. The
  cheaper options are thinner film or no bag at all: a loose snow sheet made in real time just
  ahead of the pulse, with no film to lose.

### 6. Debris

- **A struck column vaporizes** (item 2), so a hit leaves no charge debris.
- **What a miss leaves: the intact bag, or the uncrossed part of a grazed one.** Sunlit lifetimes
  (fast-rotating, isothermal; Earth's shadow ignored):

| piece | 1 AU, bright / dark ice | 9.86 R_sun |
|---|---|---|
| 1 mm | 0.12-0.30 d | 17-32 s |
| 1 cm | 1.2-3.0 d | 3-5 min |
| 10 cm | 12-30 d | 0.5-0.9 h |
| 0.4 m (the departure's whole bag) | 49-120 d | 2-4 h |
| snow at 300, 0.4 m, dark | 16 d | 37 min |
| argon frost, 0.4 m | 8 d | 25 min |

- **The film does not sublime.** At 20 m^2/kg (against 0.005 m^2/kg for a 0.2 m solid
  polyethylene block) solar pressure and drag act on it far harder. Less film is the debris lever.

### 7. A graded charge

Ice at 917 kg/m^3 under snow at 306 or 183, with the dense part holding 25-75% of k, at either end
(on the Sun bowl the 25% dense layer is too thin for the grid and was not run):

| | best impulse over solid | over a uniform column of the same length | on-axis peak, dense at entry / light at entry |
|---|---|---|---|
| departure | 0.969 (3x step, 75% dense, dense at entry) | 0.94-1.06 | ~1.0-1.4x / **0.09-0.25x** |
| Sun | 0.989 (3x step, 50% dense, dense at entry) | 0.97-1.02 | ~0.9-1.6x / **0.20-0.51x** |
| push | 0.989 (dense at the plate) / 1.000 (dense at entry) | | |

- **Grading never beats a uniform solid charge.** Against a uniform column of the same length it
  is within a few percent either way. The gamma gas cannot see pore heating, but item 2 says there
  is none to see at these speeds.
- **The orientation.** Dense at the entry is slightly better for impulse. Light at the entry (your
  rule) cuts the on-axis peak several times over for 1-6% impulse at a 3x step.
- **The departure's mouth-facing orientation is right, for the peak.** The same holds on the Sun
  bowl once the orientation is corrected (dense toward the mouth).
- **Caveats.** Peaks are not converged, and on the departure the real hole is off the axis. This
  is a lever for the collar and hole loads to study, not a load.

## Low speeds: the 100 km/s Earth runs and the 3S Jupiter return

Seth asked whether this holds at 100, 55 and 44 km/s.
- **The hydrodynamics carries over unchanged.** The 2-D gas has no speed scale.
- **The chemistry scales as v^2 / k.** Heat per kilogram falls while the bond energies do not.
- **The set-up.** A 15 kg impactor at 3.75 m. Both gamma brackets are used, since the 3S plate at
  44 km/s is not named. Water and argon rebound runs were added at 6-12 km/s.

Impulse against argon at the same k:

| closing | k | water: common gamma / sheet / bowl | polyethylene | water's last-strike margin |
|---|---|---|---|---|
| 44 km/s | 5 | 1.10 / 1.06 / 1.10 | 0.98 | 8.8x |
| 44 km/s | **10** | **1.01 / 0.95 / 1.01** | **0.71** | **2.6x** |
| 44 km/s | 20 | 0.99 / 0.90 / 0.99 | 0 (carbon takes all of v_e) | 0.7x |
| 55 km/s | 10 | 1.04 / 1.00 / 1.04 | 0.89 | 4.1x |
| 100 km/s | 10 | 1.06 / 1.07 / 1.06 | 1.03 | 13.5x |
| 100 km/s | 20 | 1.04 / 1.03 / 1.04 | 0.97 | 3.7x |

- **At Seth's k ≈ 10 and 44 km/s, water and argon tie within about 5%.** Water's rebound (0.57)
  matches argon's (0.55). Water's bonds stay broken through the push because stagnation re-heats
  the gas, not because they freeze. At the end of the rebound push 19-48% of the stagnation energy
  is still in broken bonds at k = 5-20 (39% at k = 10), against 4-16% at 100 km/s. That is a
  softness, priced in its gamma and rebound.
- **Polyethylene loses ~29%** there to its carbon toll.
- **Argon radiates 10-100x more into the face at 6-12 km/s,** a thermal point in water's favour.
- **Argon is the only clean charge at k = 20 and 44 km/s,** but k ≈ 10 makes that moot.
- **At 100 km/s all three are within ~4% at k = 10-20.**

## Not modelled

- Strength, fracture and the bag's skin; how snow is made and kept in flight.
- The departure's 130 deg entry (run on the axis, as Result 14).
- Peak convergence; solid argon on the Sun bowl (too short for the grid, held at solid's).
- Pore physics at 44 km/s and k = 20, where it could matter.
- The 3S plate's own geometry at 44-55 km/s.
