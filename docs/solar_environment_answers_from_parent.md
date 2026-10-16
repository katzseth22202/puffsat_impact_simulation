> **Copied verbatim** from `katzseth22202/Balloon-Pulse-Propulsion` commit `53ca8e8`
> (`docs/solar_environment_answers_for_impact_sim.md`), 2026-10-16; first copied at `1938f94`. It
> answers this repository's `docs/solar_environment_4rsun_ask_for_parent.md`. What this repository
> did with it is under "Answered" in that ask and in ADR-0063's amendments of 2026-10-16. Its items
> 4-5 (the CME asks) are answered in `docs/cme_drift_answers_for_parent.md`.

# The near-Sun environment at 4 and 9.86 R_sun: answers to E1-E6, G5-G8 (and G9)

Answered 2026-10-16 in `Balloon-Pulse-Propulsion` (after paper commit `86b6671`), for
`puffsat_impact_simulation` @ `a6742a2`, `docs/solar_environment_4rsun_ask_for_parent.md`
(copied here as `docs/solar_environment_ask_from_impact_sim.md`). **Written to be copied
verbatim back into the simulation repository.**

Ground rules: neither document is the source of truth. Numbers marked **(lit.)** come from the
cited source. Numbers marked **(est.)** are the parent's own scaling, shown so they can be checked.
They are not results. Several items have no published value at 4 R_sun at all, because nothing
has flown there, and they say so rather than inventing one.

Distances used: 9.86 R_sun = 0.0459 AU (Parker's flown perihelion); 4 R_sun = 0.0186 AU.

## The answer, in one table

| item | 9.86 R_sun | 4 R_sun | verdict |
|---|---|---|---|
| steady sunlight | 0.65 MW/m^2 (est., 1/r^2) | 3.93 MW/m^2 (est.) | **design heat load at both depths** |
| CME heat | ~1e-3 of sunlight (est.) | 2e-4 typical, 3e-2 extreme (est.) | negligible beside sunlight |
| CME push on the ship's plate | < 1 N typical (est.) | 0.5 N typical, ~30 N extreme on 300 m^2 (est.) | negligible for a ~600 t ship |
| **CME push on a 1 kg, 0.5 m disk** | flown edge-on: ~16 m typical, ~1 km extreme (est.) | flown edge-on: ~17 m typical, ~1 km extreme (est.) | **a small guidance load if the disk flies edge-on; tens of m/s if face-on** |
| flare heat (X17-class) | ~3e-4 of steady (lit. at 1 AU, scaled) | <= ~5e-4 of steady (est.) | negligible |
| SEP dose, one large event | x9 to x100 the 1 AU fluence (lit. scaling) | x16 to x400 (lit. scaling, extrapolated) | **the open hazard; sets electronics hardening** |
| dust | measured flux, rising inward (lit.) | inside or at the edge of the dust-depletion zone (lit.) | open in number; head-on retrograde impacts are fast |
| charging | Parker flew **positive** inside ~25 R_sun (lit.) | not measured | not a design driver for plates or stowed charges |

**One-line verdict:** steady sunlight stays the design heat load at both depths. Two items
overtake it as design drivers, and neither is heat:
1. **SEP dose** sets the hardening of every disk's and ship's guidance electronics.
2. **CME dynamic pressure on the disks** sets the disks' attitude rule. Flown edge-on (face normal out
   of the orbital plane, flipped 90 degrees to face the ship in the last minute), the correction is
   under 0.6 m/s even in an extreme CME, at roughly a 10% chance per pass at solar maximum (est.).
   Face-on to the radial flow it would be tens of m/s.

---

## E1. Coronal structure

- **Closed field sits below 4 R_sun almost always.** Helmet-streamer cusps are generally observed
  below 2-2.5 R_sun, near the inner edge of LASCO C2's field (lit., the streamer-cusp literature
  surveyed for this answer; the PFSS source surface is conventionally 2.5 R_sun). Some
  current-sheet source-surface models put the cusp near 2.5 R_sun with the source surface farther
  out. So a 4 R_sun pass flies above the closed loops and through open field, in a streamer
  **stalk** rather than its dome. At solar maximum, streamers rise radially from their sources
  instead of bending toward the equator (lit., tomography), so the streamer belt is spread over
  all latitudes. A pass near the ecliptic will cross it.
- **Current-sheet crossings.** Parker crosses the heliospheric current sheet near perihelion
  (encounters 8 and 11 imaged it closely; Liewer et al. 2023, ApJ). The plasma sheet around it has
  densities about **twice** the surroundings, with strong transverse structure (lit., Lavraud et
  al., first orbit; WISPR "fan of rays"). Streamer-belt wind is high-beta (Huang et al. 2023,
  ApJS).
- **What it means here:** a few-hour pass may cross the current sheet once or more, with a density
  contrast of ~2x and a field reversal. Neither is a heat or push driver: twice a tenuous wind is
  still tenuous. It matters only as part of the disks' plasma drag budget (see E2).
- **Not found:** cusp-height statistics at solar maximum specifically, and a crossing rate per
  hour at 4 R_sun. Both need the eclipse/coronagraph literature or an MHD model run (PSI's
  coronal models), which this answer did not do.

## E2. CME encounters

- **Rate at solar maximum:** LASCO counts reach ~4-5 CMEs a day in active rotations (lit., LASCO
  catalog studies; rates track the cycle in amplitude and phase). Catalogs disagree on totals.
- **Speeds:** typical CMEs run a few hundred to ~1000 km/s. The fastest of Cycle 25's rising phase,
  which Parker flew through on 2022-09-05 at ~13 R_sun, moved at ~2500-2700 km/s (lit., Patel et al.
  2023, ApJL, doi 10.3847/2041-8213/acf2f0; JHU APL release). Parker survived it.
- **Chance per pass (est.):** at ~5 a day, with a typical angular width covering ~5-10% of the sky
  from the Sun, a given direction is swept ~0.35 times a day. A CME body takes ~1 h to pass a point
  near 4 R_sun. So a 5-hour pass has roughly a **10%** chance of meeting one at solar maximum, and
  much less at minimum. This is an order of magnitude, not a statistic.
- **Density (est., not found in the literature searched):** taken as ~1e6 cm^-3 typical and
  ~1e7 cm^-3 extreme inside a CME near 4 R_sun. These are assumptions. Verify them against Parker's
  in-situ moments and coronagraph mass estimates before quoting.
- **Heat:** `rho v^3 / 2` is ~8e2 W/m^2 typical (2e-4 of sunlight at 4 R_sun) and ~1.3e5 W/m^2
  extreme (3e-2). **The simulation's "negligible" stands** for heat.
- **Push on the ship:** `rho v^2` is 1.7e-3 Pa typical and 0.1 Pa extreme. On a 300 m^2 plate that
  is 0.5 N and 30 N, or <0.1 m/s and ~0.5 m/s on a ~600 t ship over half an hour. **Negligible
  stands** for the ship.
- **Push on the disks, which the ask did not cover.** A 0.5 m graphite disk of ~1 kg is 2.8 mm
  thick. Its face has ~0.2 m^2 per kg, about 2000 times the ship's area per kilogram, but its edge
  has 139 times less than its face. **Attitude decides the answer.** Near perihelion the disk's
  velocity is tangential, while sunlight and a CME stream radially. Over a CME passage of ~1 h at
  4 R_sun and ~2.5 h at 9.86 R_sun (est.):

  | attitude | 4 R_sun typical | 4 R_sun extreme | 9.86 R_sun typical | 9.86 R_sun extreme |
  |---|---|---|---|---|
  | face square to the radial flow | 1.2 m/s, 2.1 km | 74 m/s, 133 km | 0.5 m/s, 2.1 km | 30 m/s, 133 km |
  | strike attitude (normal along velocity) | 0.12 m/s, 220 m | 1.7 m/s, 3 km | 0.02 m/s, 100 m | 0.4 m/s, 1.8 km |
  | **normal out of the orbital plane** | **0.01 m/s, 17 m** | **0.5 m/s, 1 km** | **0.004 m/s, 16 m** | **0.2 m/s, 1 km** |

  **Recommendation (author's proposal, 2026-10-16):** fly the disks with the face normal out of the
  orbital plane, so sunlight, the CME and the disk's own ram all meet the edge. Flip 90 degrees to
  face the ship in the last minute. Even face-on to an extreme CME for those final 60 s costs only
  1.2 m/s and 37 m. A few m/s of cold gas covers every case. The remaining requirement is
  **sensing**: tracking a 1 km drift against the ship, not paying for it.
  - **Depth barely matters.** At 9.86 R_sun the push is ~6-15x weaker (thinner plasma, slower
    disk), but a CME takes ~2.5x longer to pass, and drift grows as the square of exposure. The
    drift is within ~2x at both depths, so this is not a 4 R_sun-only item. Parker met its
    2022-09-05 CME at ~13 R_sun.
  - **Heat bonus:** edge-on at 4 R_sun the sunlit rim takes ~5.6 kW and the disk radiates it at
    **~700 K**, against 2430-2890 K face-on. That eases G8.

## E3. Energetic-particle dose (the one most likely to matter)

- **Radial scaling (lit.):** transport modelling of 8-35 MeV protons inside 1 AU gives peak flux
  falling as r^-1.3 to r^-3 and **fluence as r^-0.7 to r^-1.5**, depending on energy, the shock
  and injection (Lario et al. 2007, as reviewed in a 2024 A&A Solar Orbiter paper; Lario et al.
  2006 on 77 Helios/IMP-8 events). Dose follows fluence. **So the simulation's 1/r^2 to 1/r^3 is
  pessimistic for dose.** That scaling describes peak flux.
- **Extrapolated to the two depths (est.).** These are fluence multipliers on the same event at
  1 AU. They extrapolate a law fitted at 0.3-1 AU to 0.02-0.05 AU, so treat them as a bracket:

  | | r^-0.7 | r^-1.5 | r^-2 (simulation's floor) |
  |---|---|---|---|
  | 9.86 R_sun | x9 | x100 | x480 |
  | 4 R_sun | x16 | x390 | x2900 |
  | 4 against 9.86 | x1.9 | x3.9 | x6.1 |

- **The 1 AU anchor (lit.).** For the August 1972 event, the largest dose event of the space age,
  NASA's published bone-marrow doses are 1.6 Gy behind ~0.4 g/cm^2 (a space suit), 1.3 Gy behind
  ~1 g/cm^2 of aluminium (a pressure vessel) and 0.5 Gy behind ~5 g/cm^2 (an equipment room). Skin
  doses inside a spacecraft could reach ~2 Gy, and 15-20 Gy behind spacesuit-like shielding
  (lit., National Academies 2006 workshop report summarising the earlier studies). 1 Gy is 100 rad.
  October 1989 is the de facto electronics worst case (CREME96; 99%-confidence in the JPL model's
  >10, >30 and >60 MeV channels), and it is softer than 1972 at medium energies. No tabulated
  krad(Si) behind 1, 3 and 10 mm Al was found for it. Run CREME96 or SPENVIS on it.
- **What that gives (est.).** Take ~2 krad behind ~1 mm Al (from the thin-shield skin dose) and
  ~0.15 krad behind ~3.7 mm (1 g/cm^2) for one 1972-class event at 1 AU. Scaling gives:
  - **9.86 R_sun:** ~20-200 krad behind 1 mm, ~1-15 krad behind 1 g/cm^2.
  - **4 R_sun:** ~30-800 krad behind 1 mm, ~2-60 krad behind 1 g/cm^2.
  - These are whole-event fluences. A few-hour pass sees only part of an event at 1 AU, but near
    the Sun the event is shorter and sharper, so a pass can catch most of it. Do not discount
    for pass length until this is modelled.
- **Against parts:** Parker's EPI-Hi electronics were designed to >100 krad, with no latch-up and
  no processor crash from single-event upset (lit., ISIS PDR material, 2013). Rad-hard parts
  tolerate 100 krad to 1 Mrad. Commercial parts are often 10-30 krad.
- **Rate (est.):** large SEP events occur a few times a year at solar maximum. The chance that one
  coincides with a given few-hour pass is of order 1%. A Sun-wave disk sees one perihelion, so
  most disks will never see a large event, but a fleet of thousands will.
- **Verdict:** behind ~1 g/cm^2 a rad-hard design survives one worst-case event at either depth
  on the r^-0.7 to r^-1.5 laws. Behind 1 mm it does not reliably survive at 4 R_sun. **What the
  simulation should carry:** guidance electronics at >= 100 krad rad-hard, ~1 g/cm^2 equivalent
  shielding on the disks' guidance (a few hundred grams per disk), and a fleet-level attrition
  of ~1% of passes lost to SEP events at solar maximum. The 4 R_sun dose is ~2-4 times the
  9.86 R_sun dose on the evidenced laws.
- **Owed:** a CREME96/SPENVIS run of Oct 1989 and Aug 1972 behind 1, 3 and 10 mm Al at 1 AU,
  scaled by the table above. That is a one-afternoon job for whoever has the tool.

## E4. Flares

- **Anchor (lit.):** the X17 flare of 2003-10-28, the first definite flare in the total solar
  irradiance record, raised TSI by **270 ppm** (Woods et al. 2004, GRL, doi 10.1029/2004GL019571).
  About half its radiated energy was below 190 nm. Its visible peak was 0.028% of TSI.
- **At depth (est.):** flare light and steady sunlight both fall as 1/r^2 to first order, so the
  ratio carries over. At 4 R_sun a flare under the sub-spacecraft point is ~3 R_sun away rather
  than 4, which is worth (4/3)^2 = 1.8x. So an X10-X17 flare is **<= ~5e-4 of the steady load** for
  minutes to tens of minutes. Hard X-rays are a small part of that energy.
- **Verdict:** negligible as heat. The flare's particle burst is E3, not E4.

## E5. Dust

- **Depletion zone (lit.):** WISPR saw the F-corona's brightness fall at short elongations on
  Parker's first two passes, which Howard et al. (2019, Nature, doi 10.1038/s41586-019-1807-x)
  read as the first sign of the long-predicted dust-free zone. Press coverage puts the thinning
  at ~7 million miles (~16 R_sun) and less dust still at ~4 million miles (~9 R_sun). Those are
  secondary figures. The sublimation zone for silicate grains is conventionally put at ~4-5 R_sun.
  So **4 R_sun sits at or inside the depletion zone, and should have less dust than 9.86 R_sun.**
- **Flux (lit.):** Parker's FIELDS antennas count impacts by voltage spikes (Page et al. 2020, ApJS
  246, 51). The impactors are mostly beta-meteoroids escaping on hyperbolic orbits, with a total
  outflow of 0.5-10 t/s (Szalay et al. 2020, ApJS 246, 27). A two-population model fits the first
  six orbits, and the impact rate rises as perihelion shrinks. No per-square-metre rate at
  perihelion was extracted for this answer.
- **What matters here (est.):** the Sun-wave disks are on **retrograde** orbits, so prograde dust
  meets them at roughly twice perihelion speed, 400-600 km/s. A 1 um grain (~4e-15 kg) at 500 km/s
  carries ~5e-4 J. That pits a plate and does nothing to its impulse. It could degrade a disk's
  guidance optics if they face forward. The ship is prograde and meets dust more slowly.
- **Owed:** Parker's measured impact rate at its 9.86 R_sun perihelia (the PSP dust database of
  Malaspina et al. 2023), scaled to a 0.2 m^2 disk face and a 300 m^2 plate.

## E6. Charging and plasma interaction

- **Measured (lit.):** Parker's floating potential was **positive for most of each encounter inside
  ~25 R_sun**, directly contradicting the pre-launch prediction of large negative charging
  ("Testing spacecraft charging predictions as Parker Solar Probe approaches the Sun"; the FIELDS
  charging Parts I and II). Photoelectrons dominate. An electrostatic barrier forms at the sunlit
  shield and returns photoelectrons, charging some shaded surfaces negative. A bias range of
  -5 to +5 V relative to the spacecraft was recommended at the first perihelion.
- **Deep dielectric charging** comes from penetrating SEP electrons and so is part of E3, not a
  separate environment.
- **Stowed Kapton sheets and polyethylene charges** in the ship's shadow are dielectrics, so they
  will charge differentially. Near the Sun the plasma is dense and its Debye length short, which
  limits that. Volts to tens of volts would not move a stowed sheet or a charge. A kilovolt-class
  discharge would matter only to adjacent electronics.
- **Verdict:** not a design driver for plates, charges or sheets. Ground the guidance packages
  and keep them off floating dielectric surfaces, as Parker does.

---

## G5-G9 (the parent's share of the growth-cycle asks)

- **G5. The Sun ship's plate.** Kept as dead mass. ADR 0072 of the calcs repo shows it is
  6-19 t (1-2.5% of the ship), not 16%, because pulses capped at 2 Hz over a 2-5 h arc are small.
  Cutting it into sheets stays an option not taken.
- **G6. The Earth-wave sheet's guidance.** The calcs repo brackets the guidance mass at 2-10% of
  sheet mass, and the bracket moves the dear-seed break-even by under $1/kg (ADR 0072). So the
  parent does not size it further for the economics. Its design rules, now in the paper
  (`sec:charge_pusher_abort`):
  - sheets coast **furled** and unfurl ~100 s out, which also shrinks the solar-sail drift on the
    coast;
  - sheets **aim off the plate until their unfurl is confirmed**, so a failed unfurl (a compact,
    penetrating roll) misses by construction;
  - sheets carry enough cold gas to **divert tens of metres with >10 s of warning**, which is the
    stream-steered abort.
  The paper's four coasting co-flyers (`sec:jupiter_coflyers`) are the tracking architecture for
  each returning stream.
- **G7. How the disks reach their orbit.** Closed by your G10 and the calcs repo. They ride out
  with the ship, and a carrier burns 1.67 km/s at a 6.5 R_J perijove (about 1% of the ship) to
  put them on the true mirror dive at 15.07 km/s.
- **G8. The disk at 4 R_sun.**
  - **Heat (est.):** a flat graphite disk face-on at 4 R_sun reaches ~2430 K radiating from both
    faces and ~2890 K from one. That is below graphite's sublimation. Flown edge-on, as E2
    recommends, it runs at ~700 K.
    The disk survives. Its **guidance package does not survive bare**, and needs a shade. The
    paper's pass-through disk already carries a sacrificial foam brow for this. The ship's shade
    is held at 0.5% Parker-scaled, and the calcs repo's "shade scales with sunlight" sensitivity
    (x6) costs 4.92 -> 4.70 growth and $67 -> $75 dear.
  - **Dose:** E3. About 1 g/cm^2 equivalent over the electronics.
  - **Dust:** E5, the forward-facing optics.
  - **CME push:** E2, a correction budget of ~1 m/s typical.
  - **Peak plate pressure per pulse.** The calcs repo sizes the Sun plate as pulse energy^(2/3)
    at a fixed peak pressure, with pulses capped at 2 Hz. At 4 R_sun the bowl grows to 18.7 t
    (against 5.9 t at 9.86 R_sun) on a 2.0 h arc, so the pressure limit is already charged.
- **G9. The charge material.** Polyethylene for now (author, 2026-10-16). At the flown k ~ 6 the
  bonds take ~12% of the sheet's energy, and returning every bond moves the dear break-even by
  ~$2. A bond-free argon charge would recover that toll at the cost of a cryostat.

## New items for this repository

1. **Fly the Sun-wave disks edge-on and flip at the last minute** (E2): the face normal out of the
   orbital plane, rotated 90 degrees to face the ship in the final minute. A few m/s of cold gas then
   covers even an extreme CME (<= 0.6 m/s at either depth). Size the tracking for a ~1 km drift.
   The parent's estimate uses assumed CME densities and purely radial flow, so check both first
   (item 4).
2. **Carry ~1% of Sun-wave passes lost to SEP events at solar maximum**, and >= 100 krad
   electronics behind ~1 g/cm^2, in the disk mass and the attrition (E3).
3. **Use fluence scaling, not peak-flux scaling, for dose**: r^-0.7 to r^-1.5.
4. **Two CME data asks that decide item 1** (added 2026-10-16, after the edge-on recommendation).
   The edge-on rule assumes the CME flows exactly radially. Radial flow from the Sun does lie in the
   disk's orbital plane, but the edge has 1/139 of the face area, so the rule holds only while the
   flow stays within ~0.4 degrees of that plane. Non-radial motion (deflection, sheath flow along
   the shock, flux-rope expansion) tilts the flow onto the face, and the projected area goes as
   sin(tilt) + 1/139 of the face. Parent scratch estimate (est.), on this doc's assumed densities:

   | out-of-plane flow, % of CME speed | tilt | 4 R_sun extreme, 1 h | 9.86 R_sun typical, 2.5 h |
   |---|---|---|---|
   | 0 (item 1) | 0 deg | 0.5 m/s, 1 km | 0.004 m/s, 16 m |
   | 5 | 2.8 deg | 4 m/s, 7.5 km | - |
   | 10 | 5.7 deg | 8 m/s, 14 km | - |
   | 20 | 11 deg | 15 m/s, 27 km | 0.10 m/s, 450 m (1.1e-5 N mean on 1 kg) |

   So a few degrees of tilt brings back most of the face-on push, and the extreme case no longer
   fits "a few m/s of cold gas". Against that, a mass check suggests the assumed densities run
   ~6-10x high: a typical 1.6e12 kg CME in a 50-degree shell r/2 thick at 4 R_sun is ~1.5e5 cm^-3,
   and a heavy 1e13 kg one ~1e6 cm^-3 (est.; it spreads the mass evenly, and real cores and sheaths
   are denser). The two asks:
   - **(a) CME density at depth.** Parker's measured proton density (SWEAP moments) inside the
     2022-09-05 CME at ~13 R_sun, and coronagraph mass-and-volume densities for CMEs at 4-10 R_sun,
     typical and 99th percentile. Peak-to-mean ratio across the sheath and core too, since the
     guidance sees the peak.
   - **(b) Non-radial flow.** A bound on the out-of-plane velocity component inside CMEs near the
     Sun, as a fraction of the bulk speed (Parker in-situ flow angles, or MHD CME runs). This sets
     the tilt in the table above.
   Together these decide between a ~1 km tracking requirement and a tens-of-m/s correction budget
   for roughly one pass in ten at solar maximum.
5. **Correction to item 1's sunlight figure.** The Sun is not a point source at depth. Its angular
   radius is 14.5 degrees at 4 R_sun and 5.8 degrees at 9.86 R_sun, so an edge-on disk still sees
   the upper and lower halves of the Sun on its two faces. The pressure across the face cancels by
   symmetry, but the outward push does not. Averaged over a uniform solar disk, edge-on sunlight is
   ~11% of face-on at 4 R_sun and ~5% at 9.86 R_sun, not 1/139. A 1% modeling error then misses by
   ~80 m at 4 R_sun (not ~5 m), still a steady, plannable term.
6. **Sun-wave co-flyers: how the disks' drift is measured** (Seth's proposal, 2026-10-16; parent
   scratch numbers, est.). The paper gives the Sun wave no tracking architecture yet. The four
   co-flyers of `sec:jupiter_coflyers` serve the Earth-return streams only.
   - **Layout.** Four co-flyers ride the disks' retrograde dive, ~500,000 km off the train. The
     train is ~3.4e6 km long (4.85 h at 197 km/s at 9.86 R_sun). Put them at the corners: each
     offset both radially and out of plane, two near a quarter of the way along the train and two
     near three quarters. The Sun ship ranges along-track as a fifth station.
   - **Geometry (dilution of precision, two-way ranging).** Corners give a worst case of 2.0 in the
     radial and normal components; two out-of-plane plus two radial give 3.7. Losing one corner
     co-flyer still leaves 3.7, so four carries one spare. Corners also keep the Sun at least 45
     degrees off every line of sight. A purely radial-outward co-flyer would look at its nearest
     disks straight against the Sun.
   - **Cost of the offsets (two-body, 9.86 R_sun dive from 5.2 AU, aphelion speed 1.73 km/s).**
     Along-track is free (timing). Radial, a perihelion of 9.14 or 10.6 R_sun, is ~60 m/s at
     aphelion. Out of plane is ~1.3 km/s, kicked near 2.6 AU inbound (1.7 km/s at 1 AU, 3.1 km/s
     at 0.3 AU). A plane change at Jupiter does not help: every orbit leaving its aphelion with
     purely tangential velocity keeps the same perihelion point, whatever its plane. A
     differently targeted Jupiter flyby might supply the turn for free. That is unchecked against
     the flyby's turning limit, and is ask (c) below. The inner co-flyers at 9.14 R_sun see ~16%
     more sunlight.
   - **What sets the precision: the network's shape.** Ranges measure distances only, so they fix
     the shape of the network (ship, co-flyers, disks) up to a rigid rotation and translation. A
     hit or miss depends only on the disk's state relative to the ship inside that shape, which a
     rigid motion leaves unchanged. Cross-ranging between all nodes fixes the shape to
     centimetres, at the same precision as ranging to the disks. Orientation enters only through
     the Sun's differential gravity on the predicted closing path: a 1 microradian error on a disk
     1e6 km (~43 min) out mispredicts its miss by ~4 mm. So the floor is centimetres, not the
     ~40 m an assumed 20 m co-flyer position error once suggested. Lunar laser ranging is
     millimetre-class for the same reason it works here: the geometry it needs is known, there
     because its stations are tied into the ground reference frame, here because the network
     ranges itself.
   - **Light-path terms.** The Shapiro delay between two nodes 2e6 km apart at ~10 R_sun adds ~860 m
     of path. General relativity predicts it to ~2 cm (PPN gamma to ~2e-5, Cassini), so it is
     modelled, not budgeted. Coronal plasma delay at optical wavelengths is nanometres.
   - **Co-flyer accelerometers.** They only bridge the gaps between range fixes, not the whole arc.
     A bias b grows into ½·b·t² of position error. Across the whole 4.85 h arc with nothing else,
     1 cm would take ~1e-10 m/s^2. Re-measured by ranging every ~100 s, it takes ~1e-6 m/s^2,
     because each fix wipes out what built up since the last. Delta-DOR is not what saves this
     (it is ~150 m at 1 AU, and fixes only the network's position, which a hit does not need).
     The accelerometer may even be optional. A string of range fixes already measures how the
     shape accelerates. It earns its place by smoothing over gaps (a disk's transponder pointed
     at another co-flyer), and by telling whether a co-flyer or a disk was pushed, which ranges
     alone cannot separate.
     GOCE-class is the wrong spec. Its full range of 6.5e-6 m/s^2 worked only because GOCE flew
     drag-free, and sunlight alone pushes a co-flyer with 1 m^2 per 50 kg at ~4e-5 m/s^2 at
     9.86 R_sun. Needed: a range of ~1e-3 to 1e-2 m/s^2, and a bias holding ~1 cm over ~100 s
     between fixes, i.e. ~1e-6 m/s^2. That is the same electrostatic technology set for wider
     range. A ~13 kg, ~20 W simplified gravitational reference sensor (Conklin et al., arXiv
     2107.08545) bounds the mass.
   - **On each disk:** a pointed laser transponder, stepped between co-flyers. 1 W peak in a 1 mrad
     beam into a 30 cm co-flyer telescope returns ~1e5 photons/s at 2e6 km. The paper's wide,
     milliwatt-average beacon aimed at the ship does not reach the co-flyers. Cold gas for the
     drift: ~25 g covers the worst tilted extreme CME (15 m/s at 60 s specific impulse), under 1 g
     the typical tilted case at 9.86 R_sun (0.1-0.3 m/s).
   - **Calibration disks (Seth, 2026-10-16).** The Sun-wave form of the paper's sacrificial lead
     PuffSat. The first one or two disks aim to miss the ship by tens of metres. The ship's
     terminal tracker sees the pass almost perfectly (3 microradians at 100 m is under 1 mm), and
     the gap between the network's prediction and the measured pass gives the error all disks
     share: transponder delay biases, the ship's place in the shape, sunlight-model error, and a
     CME push common to the train. Correcting for it fixes that shared error for every disk behind.
     Errors that differ from disk to disk stay each disk's own job. Aiming off first follows the
     Earth sheets' rule (aim off until unfurl is confirmed), so a gross network error shows up as
     a measured miss, not an off-centre strike. After the first hits every disk gives a residual
     anyway, through the terminal tracker in its last seconds. But the shared error can change
     mid-arc, so add a calibration miss after any long pause in pulses and after a CME passes.
     Cost: even a dozen of ~8,100 disks is 0.15% of the wave. A missed disk is harmless, since
     nothing orbits at that depth, and it coasts back out on its own orbit.
   - **Split of the correction.** A CME pushes neighbouring disks at the same attitude alike, so
     much of it is common-mode and can go to the ship. An extreme CME's ~27 km at 4 R_sun is too
     much for it: moving 1500 t that far in an hour takes ~6 kN. The disks cancel most of it
     themselves, and the ship trims what is left. CME density also varies along the train, which
     needs per-disk thrust anyway.
   - **Asks.**
     - **(c) aim repo:** can a differently targeted Jupiter flyby put a co-flyer ~500,000 km out
       of the disks' plane at perihelion, inside the flyby's turning limit, and at what cost
       otherwise?
     - **(d) control sim:** the disk-to-charge landing tolerance for the Sun wave (the paper says
       only "within meters"), and an error budget for the co-flyer network against it: shape
       noise, accelerometer gaps between fixes, and the CME push between fixes.

## Sources

- Banks, Backus & de Groh, NASA/TM-2008-215490 (for the paper's debris sentence; not an E-item).
- Howard et al. 2019, Nature 576, 232, doi 10.1038/s41586-019-1807-x.
- Szalay et al. 2020, ApJS 246, 27; Page et al. 2020, ApJS 246, 51.
- Patel et al. 2023, ApJL, doi 10.3847/2041-8213/acf2f0; JHU APL release 2023-09-15.
- Woods et al. 2004, GRL, doi 10.1029/2004GL019571.
- Lario et al. 2006 (Helios/IMP-8, 77 events); Lario et al. 2007 transport modelling, as reviewed
  in A&A 2024 (doi 10.1051/0004-6361/202452591).
- National Academies 2006, *Space Radiation Hazards and the Vision for Space Exploration*, ch. 3;
  NASA technical report on August 1972 BFO doses (one of NTRS 19880016043 / 19910021747; the
  search summary did not say which, so confirm before citing).
- Parker FIELDS charging: "Testing spacecraft charging predictions as Parker Solar Probe approaches
  the Sun"; FIELDS instrument charging Parts I-II (ESS Open Archive, 10.1002/essoar.10504357.1,
  10.1002/essoar.10504358.1).
- EPI-Hi PDR material (SwRI ISIS wiki, 2013): >100 krad design goal.
- Huang et al. 2023, ApJS (high-beta streamer-belt wind); Liewer et al. 2023, ApJ (WISPR imaging
  of the current-sheet plasma sheet); Lavraud et al. (first-orbit current sheet).
