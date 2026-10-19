# The growth wave above 100 km/s: answers for `Balloon-Pulse-Propulsion` (S40)

Answers 2026-10-19 to the parent's `docs/fast_growth_wave_asks_for_impact_sim.md` (S40, commit
d9c7fa1). **Written to be copied into the parent.** Evidence: Result 15 of
[`xray_collision_study.md`](xray_collision_study.md) and
[ADR-0066](adr/0066-fast-growth-wave-graphite-films.md); `make xray-collision-fast-wave`
(`fast_wave_*.csv`). **Status: proposed.** Kapton stays wherever its window is real; dry carbon
fabric replaces it only where the speed rules Kapton out (Seth, 2026-10-19).

## Short answers

1. **Your ballpark is right**, to within a few kelvin and ~7% in pressure.
2. **At 137 km/s Kapton's window is nil at your 2 m standoff.** It opens to 143-188 um on maraging
   at 4 m, if the plume's spread is u >= 8 km/s. At 208 km/s there is none. At u = 5 there is none
   at any speed, 100 km/s included.
3. **Where the speed rules Kapton out, fly dry carbon instead, still as a planar sheet**: an
   all-carbon sheet of 0.06-0.23 kg/m^2, best a dry woven carbon fabric (section 3). Kapton is
   kept at 99 km/s, where it holds. The 4 R_sun headline (137 km/s) needs carbon at the flown 2 m
   standoff, as do faster returns and fast interplanetary departures that maximize the crossing
   speed. Compact impactors of any material do not work: as a point blast they need a 35-49 t
   plate and lose ~1/5 of the impulse. The overtake law does not depend on the impactor's
   material, but it is 1-4% below your `push_ve` at 99-137 km/s and 9% below at 208.
4. **A bare compact impactor never touches the 130 deg bowl.** It leaves through the mouth. At
   117 deg (4 R_sun) it vaporizes a ~1 m^2 patch near the rim (0.14-0.55% of the impulse); moving
   the charge to 2.5 F avoids that.
5. S36 stays open. A carbon sheet removes the push's dependence on it.

## 1. Result 13 at the real closing speeds

Your plume model with each row's speed, k and v_e, conducted through bare Kapton (both faces 0.8).
The plate side now uses the 2-D peak coefficient at each row's own k: `c` = 0.65 / 0.61 / 0.60 at
k = 6.67 / 7.93 / 8.21. That is why our 208 km/s peak is 3.7 GPa against your 4.0.

| closing | dose, u = 5 / 8 / 12 | 100 um, u = 8 | 250 um, u = 8 | 100 um on the face, 2 m |
|---|---|---|---|---|
| 99.2 | 97 / 38 / 17 kJ/m^2 | 573 K | 528 K | 0.89 GPa |
| 132.1 | 134 / 53 / 24 | 708 K | 648 K | 1.52 GPa |
| 137.1 | 131 / 52 / 23 | 708 K | 648 K | 1.73 GPa |
| 207.8 | 212 / 83 / 37 | 1045 K, 52 um lost | 950 K, 25 um lost | 3.71 GPa |

## 2. The Kapton window at 137 km/s

Thermal floor: the thinnest film whose face stays under 400 C to the last 10 ms. Pressure ceiling:
the face's allowable over `c w^2`. At u = 8:

| face | 2 m | 4 m (94% of the impulse) | 8 m (84%) |
|---|---|---|---|
| HY-100 / 4340 (1.5 GPa) | none (143 > 87 um) | none (143 > 113) | 143-161 um |
| maraging 300 (2.5 GPa) | 143-144 um | **143-188 um** | 143-268 um |
| maraging, prestressed (3.5 GPa) | 143-202 um | 143-263 um | 143-375 um |

- At 132 km/s the floor is 141 um and maraging at 2 m gives 141-165 um. At u = 12 the windows are
  wide (34-87 um on HY-100 at 2 m).
- The pressure ceilings are upper values: the 2-D peak still grows ~38% from the coarse grid to
  the production one. So the 2 m maraging window at 137 km/s is nil, and the 4 m one is thin.
- **No face opens one by layering or coating.** The pulse is gas stagnation pressure, which
  reaches the steel through any layer (`docs/spray_plate_steel_face_limits.md` in the sim). The
  levers are standoff (lower `c`), a stronger face, or a material that runs hotter: graphite.

## 3. Carbon growth impactors, where Kapton fails

**Geometry: a film.** It is flown like the Kapton sheet, at the same areal-mass ceiling. With
graphite's real specific heat (710 J/kg/K at 300 K, ~1800 at 1000 K) and radiation from both faces
at 0.85, every film run holds. The face at the ceiling thickness, at the worst spread (u = 5):

| closing | HY-100, 2 m | maraging 300, 2 m | maraging 300, 4 m |
|---|---|---|---|
| 132.1 km/s | <= 64 um (6.9 m radius), 956 K | <= 106 um (5.4 m), 751 K | <= 139 um, 669 K |
| 137.1 km/s | <= 56 um (7.5 m), 1011 K | <= 93 um (5.8 m), 790 K | <= 121 um, 702 K |
| 207.8 km/s | <= 26 um (8.9 m), 1972 K | <= 44 um (6.9 m), 1533 K | <= 57 um, 1326 K |

- Sublimation is under 0.3 nm in every case run, the 2580 K face of a 10 um film at 208 km/s
  included. The C atoms are taken from thermochemistry, x10 for C2/C3.
- Your lumped check is reproduced: 1542 K against your 1591 for 100 um at 208 km/s, u = 5. With
  the real c_p and conduction it is 976 K.
- Your 1 mm plate runs hotter than lumped (575 K, not 435), because heat crosses graphite's planes
  slowly (~5 W/m/K). It is ruled out anyway: at 208 km/s it puts 57 GPa on the face.
- At 208 km/s on HY-100 the 26 um film is 8.9 m in radius, so the plate grows to ~10.4 m. On
  maraging it stays near the flown 9.5 m.

**The material: dry woven carbon fabric.** These runs use graphite's properties, but every
all-carbon form heats the same per kg/m^2: the same specific heat, emissivity and sublimation, and
nothing to decompose. So the choice is set by handling. A graphite film is brittle; a dry
(resin-free) carbon fabric drapes, rolls and folds.
- Spread-tow fabrics are sold at ~64-90 g/m^2, inside every ceiling at 132-137 km/s (2-3 layers fit
  on maraging). At 208 km/s, 64 g/m^2 fits maraging at 2 m and HY-100 at 4 m.
- In-plane, a 0/90 weave has almost no shear stiffness on the bias, and a spinning sheet's radial
  and hoop stresses differ, so it would trellis and wrinkle. A triaxial weave (0/+-60, flown in
  silicone-matrix reflectors) or light stitching fixes it. The flatness tolerance is loose (a 10 cm
  sag spreads the arrival by ~1 us against a ~100 us pulse; Result 11).
- Fabric comes ~1-1.3 m wide, so the 6-9 m sheet is seamed; its mass per area must stay even on the
  metre scale. Avoid pressed creases: carbon filaments break at a bend radius of ~0.2 mm.
- Ceramic yarns (alumina, silica, SiC) are no better here: they are more brittle (except silica),
  weaker per kg, poorer emitters when hot, and softer at 1900-2400 K. Their usual advantage,
  surviving hot air, does not apply in vacuum.

**The overtake law, graphite onto polyethylene** (Result 11's flat plate, carbon atomic, re-priced
at the closing speed):

| closing | k | your `push_ve` | re-priced, graphite |
|---|---|---|---|
| 99.2 | 7.07 | 36.7 km/s | 35.5 km/s |
| 132.1 | 7.93 | 44.7 | 44.1 |
| 137.1 | 6.67 | 53.0 | 51.3 |
| 207.8 | 8.21 | 68.4 | 62.3 |

- **The bond toll does not choose the material.** The toll is the charge's polyethylene; the
  impactor is only `1/(1 + k)` of the gas. Graphite's 59.7 MJ/kg against Kapton's ~53 MJ/kg moves
  `v_e` by under 0.1%. Charging the impactor's bonds at all (never done before) costs 0.2-0.8%.
- **The law is not linear in w.** Above 130 km/s the toll's share falls (1.3% at 208 against 5% at
  100), but the gas arrives hot enough to radiate. Polyethylene's rebound falls from 0.68 to 0.51,
  and `v_e` is 6% below linear at 208 km/s.
- At 99-137 km/s most of the gap is your law's interpolation, linear in log k between the k = 5
  and 10 rows, which overstates by ~2.5% near k = 7. The rows are in `fast_wave_law.csv`.

**Compact impactors, any material: the point-blast plate** (a 0.25 m disk into a compact polyethylene charge, on the
flown 9.5 m plate with its 2 m skirt; 14-29 kg per 7.35 MN s pulse):

| closing | sheet at 2 m | compact at 4 m | compact at 6 m | point-blast floor (energy^(2/3)) |
|---|---|---|---|---|
| 132.1 km/s | 1.5 GPa | 3.1 GPa | 1.1 GPa, 81% of the impulse | 34.8 t |
| 137.1 km/s | 1.7 GPa | 3.9 GPa | 1.6 GPa, 83% | 37.1 t |
| 207.8 km/s | 3.7 GPa | 3.1 GPa | 1.8 GPa, 80% | 49.4 t |

- The compact peak is on the axis. Its hottest ring takes 1.4-7.5x the plate's mean impulse per
  area; the sheet's takes 1.7x.
- At 2.7 m it is 5-9 GPa. So a compact impactor costs ~1/5 of the impulse and 3-4x the plate.
  This comes from the shape, not the material: the impactor never touches the plate, only the
  merged gas does.
  ADR-0064's 12.2 t holds only for a planar sheet.

## 4. A bare compact impactor on the bowl (S30/S31 extended)

**The geometry decides first.** With no charge at 1.5 F, the impactor's line runs through the
hole and the charge's place. It meets the far wall only if the stream is shallower than the line
through the rim: 122.0 deg for the δ/D 0.5 bowl with the charge at 1.5 F.

- **9.86 R_sun, 129.9 deg: no hit.** The impactor leaves through the mouth 2.5 m inside the rim and
  continues on its hyperbolic path. No abort spray is needed. Keep the stream above 122 deg
  through the burn.
- **4 R_sun, 117.0 deg: a hit**, at r = 9.18 m, 1.6 m below the rim, 34 deg to the wall, 61 km/s
  normal.
  - *Peak:* 4.4-6.6 TPa planar (a 9.2 kg, 8.3 mm graphite disk on 3.2 mm of steel; linear
    Hugoniots, extrapolated).
  - *Spall:* none, because the patch does not spall, it vaporizes. The merger dissipates 386 MJ per
    kg of disk and patch, 55x what vaporizes iron.
  - *The hole:* at least the disk's 0.89 m^2 footprint. Its loss is 0.14% of the impulse, 0.55% at
    twice the diameter.
  - *Debris:* vapour and droplets leave along the stream at ~91 km/s, away from the ship.
  - *Later bare hits:* they pass through the same hole. As the stream turns over the burn, the line
    slides and would cut a slot.
- **The fix at 117 deg is geometric:** with the charge at 2.5 F the critical angle is 110.6 deg,
  and the bare impactor leaves 2.6 m inside the rim. At 2.0 F it clears by only 0.2 m and the
  0.4 m disk clips the lip. The bowl's capture at 2.5 F has not been run.
- **The buffer spray (S31's question):** to bring the struck wall's peak down to a normal pulse's
  (0.18-0.25 GPa), the methane must sit 6.7-8.4 m off the wall. That is 18-63 kg at k = 2-5,
  sprayed 0.07-0.4 s ahead at 20-100 m/s. That is a spare charge near the axis, not a face buffer.
  Moving the charge is cheaper.

The push plate's own S30 (a bare Kapton or graphite film on the flat plate) is not answered here.

## Suggested changes to the paper (for the parent to make)

- Where the paper says the faster returns "probably need graphite sheets": *Kapton serves the push
  up to about 100 km/s. Above about 130 km/s, which includes the 4 R_sun loop and faster departures
  that trade growth for crossing speed, the sheet is a dry woven carbon fabric of 60-230 g/m^2,
  flown flat like the Kapton. Carbon's heat capacity triples by 1000 K and it radiates hard, so
  even at the plume's slowest spread it peaks at 750-2000 K and sublimes less than a nanometre.
  Kapton at 137 km/s holds only on a maraging face 4 m off, and only if the spread is above
  ~8 km/s.*
- The table of sheet heating: use section 1's rows.
- `push_ve` above 130 km/s: use `fast_wave_law.csv` (62.3 km/s at 208, not 68.4).
- The departure abort: *At the visible sign's 130 deg a bare impactor leaves through the bowl's
  mouth without touching it. At 4 R_sun's 117 deg the charge is held at 2.5 focal lengths so that
  it does too.*
- Register: S40 answered (items 1-4); S36 still open; the push plate's S30 and S31 still open.
