# Spray plate: spall and fatigue of the maraging floor over 3,000 pulses

Status: analysis note, 2026-10-06. Serves [ADR-0055](adr/0055-argon-arm-of-the-spray-plate.md)
and [`argon_plate_owed_to_companions.md`](argon_plate_owed_to_companions.md) P12. Reproduce with
`make water-plate-plate-fatigue`. The code is `python/puffsat/water_plate/plate_fatigue.py`.

## 1. Bottom line

**The merged pulse is the design, and it is what keeps the plate alive.**

- **Merged (the design intent).** The 2-D bowl at 12 MN·s peaks at 0.89 GPa and rises over ~36
  µs. That is slower than a stress wave's 10.6 µs round trip through the 30 mm floor, so the floor
  never goes into tension. Notched fatigue damage over 3,000 pulses is 0.03, a life of ~10^5
  pulses. A 0.5 mm flaw does not grow. **It survives with large margin.**
- **Unmerged for 3,000 pulses: dead, in every case.** The stratified pulse peaks at 2.5-4.1 GPa
  and rises in 0.5 µs. It reflects off the back face as 0.9-3.6 GPa of tension. That is
  0.22-0.87 of the measured spall strength, so no single pulse spalls, but every case fails
  in fatigue:

  | unmerged case | peak face | peak tension | pulses to notched-fatigue failure | pulses to critical crack, 300 / 350 |
  |---|---|---|---|---|
  | argon 45.58 km/s | 2.5 GPa | 0.9-2.1 GPa | 0.7-23 | 4-1,053 / 0-393 |
  | argon 65.13 km/s | 3.6 GPa | 1.6-3.1 GPa | 0.1-3 | 0-200 / 0 |
  | water 45.58 km/s | 3.7 GPa | 1.8-3.2 GPa | 0.04-1.8 | 0-115 / 0 |
  | water 65.13 km/s | 4.1 GPa | 2.0-3.6 GPa | 0.03-1.2 | 0-68 / 0 |

  The ranges span floor thickness (30-36 mm) and damping (Q 30-1,000). Here 0 means a 0.5 mm flaw
  is already past critical: fast fracture on the first pulse. The 3.6-4.1 GPa faces also exceed
  the 2.5 GPa allowable and maraging 300's ~3.1-3.6 GPa HEL, so the face yields every pulse.
- **Fault tolerance.** An occasional failed merge is survivable only at the gentle end: argon at
  45.58 km/s, on maraging 300, with damping. Even there the budget is tens to hundreds of
  unmerged pulses over the plate's life, not thousands. Water at 65 km/s can break a flawed plate
  on its first unmerged pulse. **A failed merge must be detected and the push stopped within a
  few pulses.** Merge quality is a structural requirement, not just an efficiency one.
- **Maraging 300, not 350.** Its `K_Ic` is double (66.5 against 33 MPa m^0.5), so its critical
  crack is four times larger. In the unmerged cases 350 tolerates essentially nothing.

## 2. The requirement this sets: rise time at least one round trip

The 2-D grid (12.5 cm cells) smears the merged front, so its 36 µs rise is partly numerical.
Imposing a sharper front on the same merged pulse (30 mm floor) gives:

| front rise | Q 100: tension, damage, crack pulses (300) | Q 1,000: tension, damage, crack pulses (300) |
|---|---|---|
| 0.5 µs | 0.36 GPa, 0.65, 3,256 | 0.68 GPa, 21, 95 |
| 2 µs | 0.30 GPa, 0.52, 5,506 | 0.58 GPa, 12, 133 |
| 5 µs | 0.20 GPa, 0.36, 15,834 | 0.38 GPa, 4.2, 400 |
| **10 µs** | **0.02 GPa, 0.04, >10^5** | **0.04 GPa, 0.04, ~10^5** |
| 20 µs | 0.01 GPa, 0.03, >10^5 | 0.04 GPa, 0.03, >10^5 |

Damage is notched, tensile cycles only, over 3,000 pulses.

**The face pressure must rise over at least one thickness-mode period, `2L/c` ≈ 10.6 µs for 30
mm.** A ramp that long cancels the floor's ringing, independent of damping. A faster front makes
survival depend on damping, which nothing in the design supplies yet. Steel alone has Q of order
10^3, and a 150 µm pitch film adds little. Two levers follow:

- **Keep the front slow.** The merged gas and the spray near the face act as a cushion. Resolving
  the real merged front needs a finer 2-D run near the plate, and this is the next check.
- **Add damping, or remove the free back.** If the floor's back face is backed by a lossy layer
  matched in impedance, the transmitted wave leaves instead of reflecting as tension. This is a
  wave absorber, not a toughness layer, and its bond line would see the cycling itself.

## 3. Method

- **Loads.**
  - Merged: the d/D 0.30, 2 m skirt bowl from `make water-plate-argon-2d-shape`, hardest-hit
    floor band, scaled to 12 MN·s as `skirt_hoop.py` does.
  - Unmerged: new 1-D real-EOS runs at the 150 t design point (k = 8.52, 92 / 65 kg PuffSat,
    4 m cloud, 4 m pulse, bare face, `make water-plate-face-history`).
- **Stress wave.** The model is linear elastic, in uniaxial strain, with `c_L` = 5,654 m/s. The
  front carries the traction `sigma = -p(t)` and the back is free. Damping is a quality factor Q
  on the thickness mode, swept over 30, 100 and 1,000. The model is checked against the closed
  forms: a held load settles to the rigid-body gradient, and a short pulse reflects as tension P.
- **Spall.** Peak tension against 4.1-5.4 GPa. That is `0.5 ρ c_b W` from the pull-back velocities
  of an 18Ni-9Co-5Mo maraging (σ0.2 2.0 GPa, the 300 chemistry) [M17].
- **Stress-life.** ASTM E1049 rainflow on `sigma_x` at 13 depths. Miner sums against NI76's
  lower-band constant-life diagrams for 18Ni1700/1900, smooth (Fig. 14) and notched Kt 2.4-2.5
  (Fig. 15) [NI76], read by eye to about ±20 MPa. A compressive mean is taken as zero. There is
  no endurance limit: the 10^5-10^7 slope is extended. Both choices are conservative.
- **Damage tolerance.** An embedded penny crack parallel to the face (the spall plane) starts at
  0.5 mm radius, about an ultrasonic call-out. It grows by `da/dN = 1.36e-10 ΔK^2.25` (m/cycle,
  MPa m^0.5), the martensitic-steel upper bound [BR], on tensile cycles only, with no threshold.
  It fails when `K = (2/π) σ_max √(πa)` reaches `K_Ic`: 66.5 (300) or 33 (350) MPa m^0.5 [NI76
  Table 12].

## 4. Limits

- Elastic only. Above the HEL the face yields and the ringing is overstated, but those cases fail
  on peak pressure anyway.
- **Not modelled:**
  - In-plane bending. The taper is meant to remove it.
  - Thermal stress at the face.
  - The joint to the gas spring behind the floor, which is treated as a free surface. That is
    the worst case for spall.
- The merged load is effective-gamma 2-D, with no radiation or real EOS. It merges in that model;
  merging with real chemistry is unshown (handoff, step 2c).
- The Paris constants are a generic martensitic bound and were not re-fetched. No maraging-specific
  ΔK_th is used. NI76 notes that maraging fatigue is sensitive to environment and surface
  condition, and recommends vacuum-melted, shot-peened parts.

## 5. Sources

- [NI76] Nickel Development Institute, "18 per cent nickel maraging steels: engineering
  properties", Publication 4419 (1976), Table 12, Figs. 14-15. Read via
  https://web.archive.org/web/2022id_/https://nickelinstitute.org/media/1598/18_nickelmaragingsteel_engineeringproperties_4419_.pdf
- [M17] Yu. I. Mescheryakov, N. I. Zhigacheva, A. K. Divakov, Yu. A. Petrov, "Shock-induced
  structural instabilities and spall-strength of maraging steels", Materials Physics and Mechanics
  32, 152-164 (2017), Table 5. https://mpm.spbstu.ru/userfiles/files/MPM232_07_mescheryakov.pdf
- [BR] J. M. Barsom and S. T. Rolfe, *Fracture and Fatigue Control in Structures*: martensitic
  steels, `da/dN = 0.66e-8 ΔK^2.25` (in/cycle, ksi √in). Not re-fetched here.
- ASTM E1049, rainflow counting. Implementation checked against its worked example.
