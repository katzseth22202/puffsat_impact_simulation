# CME density at depth, non-radial CME flow, and the edge-on disk's drift: answers to item 4 (and a check of item 5)

From `katzseth22202/puffsat_impact_simulation`, 2026-10-16. It answers items 4(a), 4(b) and 5 of the
parent's `docs/solar_environment_answers_for_impact_sim.md` (`Balloon-Pulse-Propulsion` 53ca8e8;
our verbatim copy is `docs/solar_environment_answers_from_parent.md`). It is written to be copied
into the parent as `docs/cme_drift_answers_from_impact_sim.md`. Asks (c) (the aim repo) and (d)
(the control sim) are not this repository's.

Every number below comes from `python/puffsat/xray_collision/cme_drift.py` (`make
xray-collision-cme` -> `data/results/xray_collision/cme_{density,tilt,drift,edge_sunlight}.csv`).
The literature inputs are cited in the module, and the scaling, drift and sunlight integral are
computed there. **(lit.)** marks a published value and **(scaled)** this repository's inward
extrapolation of one.

## The answer, in four lines

1. **Density (4a).** The parent's typical 1e6 cm^-3 at 4 R_sun is 1.6-4x high on coronagraph
   masses and ~10x high on in-situ fits. Its extreme 1e7 cm^-3 is ~2x high against the densest of
   29 measured CMEs times Parker's measured shock compression. In pressure the parent is ~4x high
   typical (1.7e-3 against 4.6e-4 Pa) and ~2x high extreme (0.10 against 0.047 Pa).
2. **Peak to mean: ~2.** Parker measured a shock compression of 2.03, and the 1 au sheath-to-ejecta
   ratio of the same 29 events has a median of 2.25. A strong shock bounds it at 4.
3. **Tilt (4b).** Inside the ejecta the flow is near radial: ~1 degree (1 au, median). In the sheath it is not. Parker measured
   **23 degrees** out of its orbital plane downstream of the 2022-09-05 shock, at 13-15 R_sun on
   the flank. The largest of 13 years at 1 au is 18.4 degrees, also in a sheath.
4. **Verdict on item 1.** Edge-on works in the ejecta: <= 0.8 m/s and 1.4 km even in the extreme
   case. Only a tilted sheath held on the disk for the whole hour breaks it: **~14 m/s, ~25 km and
   ~23 g of cold gas at 4 R_sun** in the worst stacked case, and ~2.6 m/s and ~12 km at 9.86 R_sun.
   A typical CME at 1000 km/s costs <= 0.3 m/s and <= 0.55 km at any measured tilt. So it is a
   ~1 km tracking requirement, plus a ~25 g cold-gas reserve per disk for the rare extreme sheath.
   The parent's ~25 g figure stands.

## 4(a). CME density at depth

**Sources.**
- **Coronagraph (lit.):** Temmer et al. 2021, JGR Space Physics 126, e2020JA028380, Table A1. It
  covers 29 CME-ICME pairs of 2008-2014, with masses deprojected from STEREO COR2 at ~15 R_sun
  (1.6e12 - 1.8e13 kg). Densities are quoted at 21.5 R_sun over the GCS flux-rope volume, assuming
  self-similar expansion and uniform plasma. The paper gives a mean of (2.2-17.6)e3 cm^-3 over
  15-30 R_sun.
- **Inward scaling (scaled):** volume as r^3 (their own assumption), times the share of the mass
  already lifted. Bein et al. 2013 (ApJ 768, 31) find the mass still grows by 2-6% of its final
  value per R_sun, so the mass at r is bracketed as `1 - f (15 - r)` with f = 0.02 and 0.06.
- **In situ, extrapolated (lit. fits, scaled):** Temmer & Bothmer 2022 (A&A 665, A70) fit Helios
  0.3-1 au ICMEs as ejecta `7.1 R^-2.4` and sheath `22.3 R^-1.7` cm^-3 (R in au). The authors note
  that the coronagraph densities sit above this extrapolation.
- **Parker, measured (lit.):** Romeo et al. 2023 (ApJ 954, 168), on the 2022-09-05 CME. Electron
  density was 1962 -> 3989 cm^-3 across the shock at ~15 R_sun, ~1e4 cm^-3 in the closed-field
  interval, and below 1e3 in most of the rest. Parker crossed the **flank and wake**, not the core.
  SWEAP's SPAN-I misses part of the ion distribution there, so the paper uses quasi-thermal-noise
  electron densities.

**Densities (cm^-3).** "Densest" is the densest of the 29 events, before shock compression.

| estimate | 4 R_sun median | 4 R_sun densest | 9.86 R_sun median | 9.86 R_sun densest |
|---|---|---|---|---|
| coronagraph, mass grows 2%/R_sun | 6.1e5 | 2.2e6 | 4.7e4 | 1.7e5 |
| coronagraph, mass grows 6%/R_sun | 2.7e5 | 9.7e5 | 3.6e4 | 1.3e5 |
| Helios ejecta fit, extrapolated | 1.0e5 | - | 1.2e4 | - |
| Helios sheath fit, extrapolated | 2.0e4 | - | 4.2e3 | - |
| ambient wind fit, for scale | 3.0e4 | - | 4.5e3 | - |
| **parent's assumption** | **1e6** | **1e7** | - | - |

- **Check at Parker's distance.** At 13.3 R_sun the coronagraph scaling gives a median of
  1.9-2.1e4 and a densest event of 7.5e4, and the Helios ejecta fit gives 5.7e3. Parker measured
  4e3 in the sheath and 1e4 at the closed-field peak. Its flank crossing therefore sits near the
  in-situ fit, below the coronagraph median. That fits a flank pass and does not test the core.
  So the coronagraph row is the high (core) estimate and the Helios row the low one.
- **The parent's mass check** (~1.5e5 cm^-3 typical, "6-10x high") agrees with the in-situ
  extrapolation and the 6% mass-growth bracket. It is more pessimistic than the coronagraph median at 2%.
  Temmer et al.'s sample is weighted to large, well-observed CMEs (median mass 6.4e12 kg, against
  the parent's typical 1.6e12 kg).
- **Sheaths are thin at 4 R_sun.** The Helios sheath fit meets the ambient wind at ~13 R_sun,
  which Temmer & Bothmer read as the average distance where sheaths start to form. A fast CME
  drives its shock lower; the 2022-09-05 shock was already formed at 15 R_sun.

**Peak to mean.** The guidance sees the peak:
- Parker's measured shock compression is **2.03** (Romeo et al. Table 2).
- Temmer et al.'s 29 events give a sheath/ejecta ratio at 1 au with a median of **2.25** (range
  0.6-8.3), computed from their table.
- A strong shock bounds the sheath at 4x the gas ahead of it.
- The extreme rows below use 2.0 on the densest event.

## 4(b). Non-radial flow

| source | where | tilt out of the plane | out-of-plane share of the speed |
|---|---|---|---|
| radial (parent's item 1) | - | 0 deg | 0 |
| Al-Haddad et al. 2022, STEREO-A PLASTIC 2007-2019 (lit.) | 1 au, inside ejecta, median normal flow 7 km/s | 0.9 deg | 0.016 |
| same, the largest CME deflection of 13 years (lit.) | 1 au, a sheath | 18.4 deg | 0.32 |
| **Romeo et al. 2023, Parker shock downstream (lit.)** | **13-15 R_sun, flank sheath** | **23.4 deg** | **0.40** |

- The Parker row comes from v_RTN = (455, -358, 251) km/s downstream (their Table 2). The shock
  normal itself is 23 deg out of the RT plane and 51 deg off radial within it. The in-plane
  component (-358 km/s) meets the edge and costs nothing. Only v_N tilts the flow onto the face.
- At 1 au only ~2% of all solar-wind points lie more than 5 degrees north or south (Al-Haddad).
  The large angles are sheath flows pushed around the ejecta.
- **Self-similar expansion is radial.** A cone expanding self-similarly moves every parcel along
  its own radius from the Sun. "Flux-rope expansion" adds non-radial flow only beyond
  self-similarity: deflection, sheath draping or rotation.
- **Not found:** an out-of-plane flow statistic *inside the ejecta* near the Sun. The 0.9 degree
  row is at 1 au and may underrate the near-Sun ejecta, where deflection happens (below ~10 R_sun).
  That is the remaining gap.

## The drift, and how much of it the CME's own speed causes

The parent's disk is used throughout: 1 kg, 0.5 m, 2.8 mm graphite, face normal out of the orbital
plane, at periapsis speed (307 km/s at 4 R_sun, 196 km/s at 9.86 R_sun). It is held in the flow
for the parent's whole passage (1 h at 4 R_sun, 2.5 h at 9.86 R_sun). The disk is taken as absorbing, so
`F = rho v_rel^2 (A_face |sin tilt_rel| + D t cos tilt_rel)`.

**The CME's speed is most of the push.** The ram is `rho v_rel^2` with
`v_rel^2 = v_CME^2 + v_disk^2`.
- At 1000 km/s the CME's own flow is **91%** of it at 4 R_sun and **96%** at 9.86 R_sun.
- The drift goes as the square of the CME's speed at fixed density: 500 -> 1000 -> 2500 km/s gives
  ~1 : 3.5 : 21 at 4 R_sun (the disk's speed softens the low end).
- The face term goes as `rho v_rel v_CME sin(phi)`. The disk's own motion dilutes the tilt but
  raises the ram, so it barely changes the face push.

**4 R_sun, 1 h (dv m/s / drift).**

| case | 0 deg | 0.9 deg (ejecta) | 18.4 deg | 23.4 deg (Parker sheath) |
|---|---|---|---|---|
| in situ (Helios ejecta), 680 km/s | 0.000 / 1 m | 0.001 / 3 m | 0.02 / 35 m | 0.025 / 44 m |
| coronagraph median, 500 km/s | 0.002 / 3 m | 0.005 / 9 m | 0.07 / 123 m | 0.09 / 154 m |
| coronagraph median, 680 km/s (median) | 0.003 / 5 m | 0.009 / 16 m | 0.12 / 212 m | 0.15 / 266 m |
| coronagraph median, 1000 km/s | 0.006 / 10 m | 0.018 / 32 m | 0.24 / 437 m | 0.30 / 548 m |
| coronagraph median, 2500 km/s | 0.03 / 58 m | 0.10 / 188 m | 1.5 / 2.6 km | 1.8 / 3.3 km |
| **densest x 2.0, 2500 km/s** | **0.24 / 437 m** | **0.78 / 1.4 km** | **10.9 / 19.6 km** | **13.7 / 24.6 km (23 g gas)** |

**9.86 R_sun, 2.5 h.**

| case | 0 deg | 0.9 deg | 18.4 deg | 23.4 deg |
|---|---|---|---|---|
| coronagraph median, 680 km/s | 0.000 / 2 m | 0.002 / 7 m | 0.02 / 96 m | 0.03 / 121 m |
| coronagraph median, 1000 km/s | 0.001 / 5 m | 0.003 / 15 m | 0.045 / 204 m | 0.06 / 256 m |
| **densest x 2.0, 2500 km/s** | **0.05 / 208 m** | **0.15 / 670 m** | **2.1 / 9.4 km** | **2.6 / 11.7 km (4 g gas)** |

- **Against the parent's table:** its 20% row at 4 R_sun (15 m/s, 27 km) is close to our 23.4
  degree extreme (13.7 m/s, 24.6 km). Its pressure was ~2x high, but its tilt was half Parker's.
  Its 9.86 R_sun typical tilted row (0.10 m/s, 450 m) is ~2-4x high against ours.
- **The extreme row stacks four pessimisms:**
  - the densest of 29 CMEs;
  - Parker's shock compression;
  - the 2022-09-05 speed;
  - a fully tilted sheath for the whole hour, at a depth where sheaths are thin.
  It is a bound, not an expectation.

## 5. Edge-on sunlight (check)

Integrated over the finite solar disk (absorbing faces, plus the edge's own `D t`):

| | Sun's angular radius | uniform disk | limb-darkened (u = 0.6) | small-angle `4 alpha/(3 pi)` |
|---|---|---|---|---|
| 4 R_sun | 14.5 deg | 0.115 | **0.108** | 0.107 |
| 9.86 R_sun | 5.8 deg | 0.050 | **0.048** | 0.043 |

**The parent's ~11% and ~5% stand**, and so does its ~80 m miss per 1% modelling error at 4 R_sun.
Limb darkening lowers the share by ~6%, because it weights the Sun's centre.

## What the parent might take

- Item 1's edge-on rule stands for the ejecta and for every typical CME. Item 4's "tens of m/s"
  is a bound only for an extreme CME's sheath. Carry ~25 g of cold gas per 1 kg disk (2.5%) for
  it, or accept losing the rare disk caught in one.
- CME densities for the paper: ~1e5 to 6e5 cm^-3 typical and ~1e6 to 5e6 cm^-3 peak in an extreme
  event at 4 R_sun; ~1e4 to 5e4 typical at 9.86 R_sun.
- **Still open:** the out-of-plane flow inside near-Sun ejecta, and the thickness of the sheath at
  4 R_sun. Both need an MHD CME run (e.g. EUHFORIA or PSI). Parker's later CME passes might also
  give them.

## Sources

- Temmer, M., et al. 2021, *Deriving CME Density From Remote Sensing Data and Comparison to
  In-Situ Measurements*, JGR Space Physics 126, e2020JA028380 (arXiv 2011.06880). Table A1.
- Temmer, M., & Bothmer, V. 2022, *Characteristics and evolution of sheath and leading edge
  structures of interplanetary coronal mass ejections in the inner heliosphere based on Helios and
  Parker Solar Probe observations*, A&A 665, A70. Fig. 4 fits.
- Romeo, O. M., et al. 2023, *Near-Sun In Situ and Remote-sensing Observations of a Coronal Mass
  Ejection and its Effect on the Heliospheric Current Sheet*, ApJ 954, 168,
  doi 10.3847/1538-4357/ace62e. Table 2 and Section 2.
- Al-Haddad, N., et al. 2022, *Investigating the Cross Sections of Coronal Mass Ejections through
  the Study of Nonradial Flows with STEREO/PLASTIC*, ApJ 927, 68, doi 10.3847/1538-4357/ac32e1
  (arXiv 2110.10682).
- Bein, B. M., et al. 2013, *The Height Evolution of the "True" Coronal Mass Ejection Mass*, ApJ
  768, 31.
