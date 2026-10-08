# Pass-through X-ray pusher: recommended target and per-projectile energy ledger

Handoff to `Balloon-Pulse-Propulsion` (parent), §`pass_through_pusher` and
§`why_not_xray_drive`. It is written to be copied there. The evidence is this repository's
X-ray collision study ([`xray_collision_study.md`](xray_collision_study.md), Results 4-5;
`make xray-collision`; `data/results/xray_collision/near_sun.csv` and `near_sun_target.csv`).

**Status: provisional.** Every number comes from a one-zone fireball. It uses a Saha EOS with
NIST ionization energies, TOPS (LANL OPLIB) gray opacities for pure carbon, Rankine-Hugoniot
starting states, and homologous expansion with a thin/thick luminosity bridge. Expect factor-of-
order-unity errors in geometry and spectrum, not order-of-magnitude ones.

**Corrections after the parent's review.**
- The ram pressure on the projectile at 10 kg/m^3 is ~50 Mbar at 618 km/s, falling to ~5 Mbar at
  206 km/s, not ~4 Mbar (that was the 1 kg/m^3 value). It heats the graphite by under 1% of the
  heat, so the ledger stands. The disk flows at these pressures (see "Not modelled").
- Commit 1bc3143's message says "60% of it above 1 keV". It means 60% of the heat (57 of
  95.5 GJ, dense into 10 kg/m^3 at k = 1), which is 73% of the light.
- The one-zone model radiates most of the heat within 0.1 us, while the disk takes ~5 us to sweep
  its column. This time-graded heating is now the first item under "Not modelled".

## Recommendation

1. **Projectile: dense carbon, flattened.** Use graphite or carbon-carbon rather than foam, so it
   survives the dive, and make it a disk rather than a ball. One kilogram as a 40 cm disk 3.6 mm
   thick carries ~8 kg/m^2 instead of a 4.8 cm ball's ~140 kg/m^2.
2. **Target: a ship-held carbon lattice at a mean density of ~10 kg/m^3** (3-10 kg/m^3 is the
   useful range). It is deployed on the projectile's line just before each pulse. An isolated
   strut in sunlight at 4 solar radii (3.93 MW/m^2) settles near (S/2 sigma)^(1/4) ~ 2400 K, so the
   lattice survives the exposure. That is a hand estimate.
3. **Mass ratio: 2 kg of target per kg of projectile (k = 2)**, not 1:1 and not 3:1.

## Per-kilogram-of-projectile ledger at 618 km/s

Ship frame, with the target at rest. Every row has the same projectile kinetic energy:
**191 GJ per kg**. "Light" is the radiated part of the collision heat. The photon bands give the
photosphere-to-zone bounds. "Fireball KE" is the heat that did not radiate and became isotropic
expansion in the centre-of-mass frame. "CM KE" is the bulk motion of all the debris astern, which
no design recovers.

| case | target per kg | heat | **light** | > 100 eV | > 1 keV | emitted kT | fireball KE | CM KE (debris speed) |
|---|---|---|---|---|---|---|---|---|
| parent today: foam on foam, 10 kg/m^3, k = 1 | 1 kg | 95.5 GJ | 57.8 GJ (30% of KE) | 52-55 GJ | 27 GJ | 210-280 eV | 37.6 GJ | 95.5 GJ (2 kg, 309 km/s) |
| dense into 10 kg/m^3 lattice, k = 1 | 1 kg | 95.5 GJ | 77.6 GJ (41%) | 75-76 GJ | 57 GJ | 520-570 eV | 17.9 GJ | 95.5 GJ (2 kg, 309 km/s) |
| **dense into 10 kg/m^3 lattice, k = 2 (recommended)** | **2 kg** | **127.3 GJ** | **90.1 GJ (47%)** | **84-87 GJ** | **54 GJ** | **320-380 eV** | **37.2 GJ** | **63.7 GJ (3 kg, 206 km/s)** |
| dense into 3 kg/m^3 lattice, k = 2 | 2 kg | 127.3 GJ | 104.0 GJ (54%) | 96-100 GJ | 55 GJ | 310-345 eV | 23.3 GJ | 63.7 GJ (3 kg, 206 km/s) |
| dense into 10 kg/m^3 lattice, k = 3 | 3 kg | 143.2 GJ | 85.8 GJ (45%) | 77-81 GJ | 40 GJ | 210-280 eV | 57.4 GJ | 47.7 GJ (4 kg, 155 km/s) |

The first row checks the parent's current text. The one-zone model gets 468 eV (the parent says
about 470 eV) and 61% of the heat radiated (the parent assumes 60%). Within this model, the
parent's assumption is right for its own design.

For the recommended case, the energy splits as
`191 GJ = 90.1 light + 37.2 fireball KE + 63.7 CM KE`, which is 47% + 19% + 33%.

## Why a lattice does better

More light gets out, but not mainly because the struts add surface. When the shock reaches the
lattice it vaporizes the struts. At 10 kg/m^3 with ~10 um fibres the struts are about 130 um
apart, and their plasmas merge in well under a nanosecond. By the time anything radiates, the
lattice is a uniform plasma at its mean density. Its strut surfaces no longer exist. The lattice
wins for two other reasons.

1. **Low mean density lets light leave from the whole volume, not only its skin.** A fireball's
   light escapes only from within about one optical depth of its edge. Below that depth the light
   is reabsorbed, and the gas spends the heat on expansion instead. A 10 kg/m^3 lattice shocks to
   ~50 kg/m^3 at ~1 keV, where fully stripped carbon has kR ~ 0.04 and kP ~ 0.2 m^2/kg (TOPS). The
   whole fireball is then about one optical depth across or less. Every part of it is effectively
   at the surface, and it radiates as it cools through the K-shell recombination range
   (100-300 eV, kP ~ 9 m^2/kg at 200 eV). Solid graphite in the same collision is ~860 optical
   depths thick and radiates 0.4%.
2. **A dense projectile puts all the heat into the lattice.** A dense body sweeping a light
   medium is a snowplow. Each swept kilogram absorbs `u^2/2` at the current speed, so the whole
   collision heat lands in the target and almost none in the projectile. The projectile is only
   compressed by the ram pressure (~50 Mbar at the start, ~5 Mbar at the end). The dilute target alone holds the heat, so it starts at
   ~1030 eV (k = 1) instead of the 470 eV of foam on foam, which shares the heat with the
   projectile's own mass. Hotter and thinner both mean faster emission.

The checks that separate these effects from surface area: at a fixed 10 kg/m^3, shrinking the
fireball from 1 kg to 0.01 kg (100x more surface per kilogram) only moves the radiated share from
61% to 68%. Lowering the density from 10 to 1 kg/m^3 at fixed mass moves it to 86%. The lattice is
the practical way to build a 3-10 kg/m^3 carbon target that survives sunlight. Its geometry
matters only through its mean density.

At the rocket speed of 67 km/s the same study found the opposite limit. There, even unlimited
surface cannot help a dense collision, because a perfect blackbody surface still radiates ~30x too
slowly. See the 67 km/s note below.

## Why 2:1

- **More of the kinetic energy becomes heat:** k/(1+k), so 50%, 67%, 75% for k = 1, 2, 3.
- **Each kilogram of target gets less of it.** The fireball runs cooler (1030, 660, 470 eV start
  for k = 1, 2, 3 in the lattice case) and thicker, so it radiates a smaller share of that heat.
- **Light per kilogram of projectile peaks near k = 2:** 77.6, 90.1, 85.8 GJ at 10 kg/m^3.
- **Light per kilogram of target the ship carries falls:** 78, 45, 29 GJ. The target is spent
  mass in the parent's specific-impulse bookkeeping. So k = 2 is the most light per projectile
  without the steep loss per target kilogram of k = 3. The parent's optimizer for
  `tab:pass_through_isp` should make the final call between k = 1 and k = 2.
- **The debris carries less away astern:** the CM share falls from 50% to 33%.

## Geometry the parent should state

For the projectile to sweep rather than punch through, the target column along its path must carry
about the projectile's own areal density times `k` (§`needle_through_fog` rule). The target column
length is `k * sigma_projectile / rho_target`:

| projectile (1 kg graphite) | areal density | column at 10 kg/m^3, k = 2 | at 3 kg/m^3, k = 2 |
|---|---|---|---|
| ball, r = 4.8 cm | ~140 kg/m^2 | 28 m | 94 m |
| disk, 40 cm x 3.6 mm | ~8 kg/m^2 | 1.6 m | 5.3 m |

The radiation case has to enclose this column, and the parent's table shows that case mass
dominates the specific impulse. That is the reason for the disk, and the reason not to go below
~3 kg/m^3. The one-zone model is insensitive to the column's shape and treats the target fireball
as a sphere, which slightly understates escape.

## Suggested edits to the parent text

- **§pass_through_pusher, opacity.** Replace "The opacity of carbon at 470 eV is a number this
  paper has not computed" and the 1-10 m^2/kg bracket. TOPS gives fully stripped carbon at the
  fireball's start (46 kg/m^3, 468 eV) kR = 0.04 and kP = 0.18 m^2/kg. That puts the 10 kg/m^3
  fireball at about one optical depth, not 18-183. The 1-10 bracket applies only once the gas has
  cooled to ~200 eV.
- **§why_not_xray_drive, eq. `cooling_race`.** `F = sigma T^4 / (1 + 3/4 tau)` saturates at a full
  blackbody when the gas is thin, so it overstates emission from a transparent fireball. Thin
  emission is `4 kP rho sigma T^4 V`. With the formula as written, the 10 kg/m^3 case gives a race
  of ~0.004, which suggests nearly total radiation. The integrated answer is 61%. The verdicts
  stand, but the formula needs a thin-limit factor or a caveat.
- **§pass_through_pusher, density claims.** The trend holds and can now be quantified as radiated
  shares of 61%, 45% and 23% at 10, 30 and 100 kg/m^3 (foam on foam, TOPS opacity). Solid graphite
  gives 0.4%.
- **§pass_through_pusher and `tab:pass_through_isp`, the design and its inputs.** Adopt the dense
  disk projectile into a carbon lattice at ~10 kg/m^3 with k = 2. For the table, the "60% of the
  fireball's heat leaves as X-rays" becomes **71% of 127 GJ, or 90 GJ per kg of projectile**
  (94% of it above 100 eV), with 2 kg of target spent per kg of projectile. Re-optimize the filler
  and propellant split for the new heat and target mass. The case-delivery (80%) and plate-catch
  (80%) assumptions are untouched by this work.
- **The "lightest foam" sentence.** It can now say why: low mean density makes the whole fireball
  a radiating surface. A lattice is how to build that target.

## Not modelled; open for the parent to flag

- **Electron-ion equilibration.** The shock heats ions first. It takes an estimated ~50 ns at
  470 eV and 1.4e22 cm^-3, against ~1 us of expansion, so the radiated share should survive.
  Electrons running cooler than ions would soften the spectrum.
- **Time-graded heating under confinement (the largest open effect).** The one-zone k = 2 fireball
  radiates 61% of Q within 0.1 us and 65% by 1 us. The disk takes ~5 us to sweep its 1.6 m
  column. Early target mass is shocked at 618 km/s to ~2.2 keV and late mass at 206 km/s to
  ~160 eV. Each parcel's thin cooling time is 0.2-6 ns, so it radiates almost as it is shocked,
  while the ram pressure holds it against the disk. The radiated share probably rises because
  confinement blocks the expansion loss. The spectrum (the 54 GJ above 1 keV) and the k = 2
  optimum could shift. The ledger above does not capture this.
- **Where the light escapes.** The unswept lattice ahead is optically thick to keV photons
  (cold carbon ~200 m^2/kg at 1 keV, column up to 16 kg/m^2). Forward light preheats it rather
  than leaving, and the disk blocks the rear, so escape is mainly sideways past the rim. The
  projectile also absorbs part of the light while the hot layer is pressed against it.
- **The disk under ram pressure.** The pressure is ~50 Mbar at 618 km/s and ~5 Mbar at 206 km/s.
  One shock at 50 Mbar puts ~0.8 GJ/kg into graphite, under 1% of the heat, so the energy ledger
  stands. At these pressures the 3.6 mm disk flows like a liquid, though. Rim relief moves
  10-20 cm inward during the sweep, half the disk's radius, so it does not stay a flat rigid
  plate.
- **Non-LTE emission and two-temperature effects** at these densities.
- **The disk's thermal and structural survival,** and the lattice's survival when its struts
  shade each other, are hand estimates only.

## 67 km/s (rocket case), for completeness

The same scheme at the rocket's 56 + 11 = 67 km/s closing speed is **not viable as an X-ray
source**.

- **Condensed argon or iron**, at 1:1 or 3:1 and from 0.1 um to 10 cm radius, radiates at most
  ~1% of the heat, at ~1 eV.
- **The limit is the blackbody ceiling** `n k T c_s / 2 sigma T^4` ~ 30-50 at shocked condensed
  density. Opacity is not the problem, so surface area cannot fix it.
- **Diluted clouds of ~1e-5 kg/m^3** with a column of ~1e-3 kg/m^2 do radiate 80-99%, but as
  EUV at kT 3-6 eV.

The 618 km/s case escapes this ceiling because the energy per atom is about 100 times larger
(5.9 keV per carbon atom against 0.23 keV per argon atom).

The full 66 km/s verdict, with argon bags at 1:2, a multigroup check that the light is not
trapped, and why a gram-per-cubic-metre cloud cannot be used in the pass-through geometry, is in
[`xray_pusher_66kms_handoff.md`](xray_pusher_66kms_handoff.md).
