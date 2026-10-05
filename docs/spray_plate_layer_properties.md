# Spray plate: properties of a thin metal layer between the film and the maraging

Status: research note, 2026-10-05. Serves [ADR-0055](adr/0055-argon-arm-of-the-spray-plate.md)
decision 9 and [`spray_plate_steel_face_limits.md`](spray_plate_steel_face_limits.md). Question:
does a thin copper-alloy or superalloy layer under the ablative film help a maraging 300 plate?
Every value carries a source tag. Values not traced to a primary source say "not verified".
Values marked "derived" are arithmetic on cited values. Temperatures in K unless marked °C.

## 1. Bottom line

- Effusivity `e = sqrt(k ρ c_p)` relative to maraging 300: GRCop-84 **3.5x at 300 K, 3.3x at
  700 K**; CuCrZr **3.7-3.9x at 300 K** (700 K not verified); Inconel 718 **0.67-0.73x and
  0.82-0.85x**; Haynes 230 **0.64x and 0.83x**. Surface ΔT under a flux pulse scales as `1/e`.
- A copper layer cuts the face temperature rise by about 70%, but only if it is thicker than the
  ~0.3 mm diffusion depth of a 1 ms pulse. Its yield (0.1-0.55 GPa) is far below the 2.5 GPa
  face allowable, so it goes plastic on every pulse.
- A superalloy layer runs hotter than bare maraging (+18 to +55% ΔT), but it tolerates 650-700 °C
  against maraging's 480 °C aging limit. Inconel 718 is the only candidate here whose
  uniaxial-strain elastic range spans a 2.5 GPa pulse (derived; Section 7).
- No primary data were found on ratcheting of thin Cu or Ni-alloy layers bonded to steel under
  repeated GPa compression. The nearest evidence is rocket-liner thermal ratcheting (Section 6).

## 2. Maraging 300 (18Ni1900), the baseline

| property | 300 K | 500 K | 700 K | 900 K | 1100 K | source |
|---|---|---|---|---|---|---|
| k (W/m·K) | 21 (20 °C) | 26-27 (200-300 °C) | 28 (400-480 °C) | not verified | not verified | [NI76] T19 |
| c_p (J/kg·K) | 460 | 460 (single value) | 460 | not verified | not verified | [NI76] T19 |
| c_p, M350 fit (J/kg·K) | 519 | 553 | 577 | 592 | 597 | [K23] eq. 5 (derived) |
| ρ (kg/m³) | 8000 | | | | | [NI76] T19 |
| 0.2% proof / RT value | 1.00 (20 °C) | 0.90 (200 °C), 0.87 (300 °C) | 0.82 (400 °C), 0.73 (480 °C) | none | none | [NI76] T14 |

- RT 0.2% proof stress 1790-2070 MPa, so about 1600-1860 MPa at 200 °C and 1470-1700 MPa at 400 °C
  (derived from [NI76] T7 and T14). E 190 GPa, ν 0.30 [NI76] T7. E(T) only in a figure: not verified.
- α mean: 9.9 (20-100 °C), 10.6 (20-300 °C), 11.3 (20-480 °C) x10⁻⁶/K [NI76] T19.
- Service limits: aged 3 h at 480 °C; at that temperature over-ageing is slight even after 200 h
  [NI76]. Table 16 (200 h exposures, 18Ni1700, room-temperature proof stress) reads 1730 MPa
  as-maraged, 1730 after 370 °C, then 1140, 830 and 280 MPa after 480, 540 and 650 °C. The column
  alignment in the scanned table is uncertain at 480 °C. DSC on maraging 350 shows Ni3(Ti,Mo)
  precipitation near 550 °C and austenite reversion from 697.6 °C [K23]. Solidus: not verified.
- [K23] c_p runs 13% above [NI76] at 300 K. The [NI76] single value is used for ratios below.

## 3. GRCop-84 (Cu-8Cr-4Nb at.%)

| property | 300 K | 500 K | 700 K | 900 K | 1100 K | source |
|---|---|---|---|---|---|---|
| k (W/m·K) | 286 | 300 | 304 | 298 | 281 | [E00] eq. 15, high-Fe lots, 296-1173 K |
| c_p (J/kg·K) | 382 | 413 | 424 | 436 | 470 | [E00] eq. 12 |
| ρ (kg/m³) | 8620 ± 80 | | | | | [E02] 2.2.1 |
| YS tension, as-received (MPa) | 196 (25 °C) | 172 (200 °C) | 140 (400 °C) | 87 (600 °C) | 27 (800 °C) | [D07] T4 |
| YS tension, after 935 °C braze (MPa) | 187 | 158 | 126 | 83 | 28 | [D07] T4 |
| YS compression, as-received (MPa) | 312 (23 °C) | 167 (350 °C) | 118 (500 °C) | 67 (650 °C) | 24 (800 °C) | [D07] T7 |
| α instantaneous (10⁻⁶/K) | 14.8 | 17.1 | 18.8 | 20.2 | 21.4 | derived from [E00] eq. 22 |

- k: low-Fe powder raised room-temperature k by 7% [E05], and [E05] quotes 305-320 W/m·K over
  the liner range. Use 286 as the conservative fit.
- α: [E00] eq. 22 is `ΔL/L (%) = -0.3287 + 2.265e-4 T^1.285`. The exponent is garbled in the scan.
  1.285 is the value that reproduces [D07]'s 0.80% at 480 °C, and it gives zero near 293 K.
- E and ν: the handbook sections are empty [E02]. Not verified.
- Limits: solidus 1080 ± 10 °C, no phase change below it [E02]. Usable to about 700 °C [E05].
  A 935 °C braze costs at most 30 MPa [E05]. Warm rolling above 500 °C degrades properties [E02].

## 4. CuCrZr (Cu-1Cr-0.1Zr, C18150 / ITER grade)

| property | 300 K | 500 K | 700 K | 900 K | 1100 K | source |
|---|---|---|---|---|---|---|
| k (W/m·K) | 314-335 | 350-360 (to 573 K) | not verified | not verified | not verified | [LZ12] T1; [P20] LFA |
| c_p, pure Cu proxy (J/kg·K) | 385 | 409 | 425 | 440 | 464 | [NIST] Shomate, Cu(s) |
| ρ (kg/m³) | 8900 | | | | | [LZ12] T1 |
| YS, hard-drawn as-received (MPa) | 550 (25 °C) | 452 (200 °C) | 283 (500 °C) | 113 (650 °C) | 20 (800 °C) | [D07] T4 |
| YS, after 935 °C braze (MPa) | 106 | 94 | 66 (500 °C) | 36 (650 °C) | 19 | [D07] T4 |
| YS, SA + aged, MPH average fit (MPa) | 288 | 249 | 191 | 115 | outside fit | [Z19] eq. 3 (derived) |
| α instantaneous (10⁻⁶/K) | 15.9 | 17.8 | 19.8 | 21.8 | not verified | derived from [D07] T12 |

- c_p is pure copper's (the alloy is >98.8% Cu). It is a proxy, not a CuCrZr measurement.
- [Z19] shows the ITER k(T) curves (its refs 42, 43) only as a figure. The ITER Materials
  Properties Handbook itself and the Pintsuk et al. 2010 interlaboratory data were not retrieved.
- E: 123 GPa at room temperature [LZ12]. [Z19] eq. 1 fits tensile-curve moduli, falling to about
  81 GPa at 427 °C. Sonic data do not exist [Z19], so E(T) is not verified. ν: missing [Z19].
- α mean 20-480 °C about 18.1 x10⁻⁶/K (derived from [D07] T12).
- Limits: melting point 1075 °C [LZ12]. Aged near 460 °C for 3 h ([D07] citing ASM) or 475 °C
  for 2 h (the [Z19] batch). It over-ages above 500 °C [D07]. A 935 °C braze cuts strength by
  about 80% [D07]. Strength recovers only if the bond cycle is followed by quench and re-age.

## 5. Superalloys

**Inconel 718 (precipitation hardened)**

| property | 300 K | 500 K | 700 K | 900 K | 1100 K | source |
|---|---|---|---|---|---|---|
| k, LFA (W/m·K) | 9.94 | 13.24 | 16.61 | 20.1 (eq. 9 end) | 22.7 | [A19] T2 |
| k, from resistivity, aged (W/m·K) | 11.4 | 14.8 | 17.9 | 21.0 | 24.3 | [SMC] T5 (°F interp.) |
| c_p (J/kg·K) | 425 | 468 | 510 | 592 (eq. 2) | 635 | [A19] T2, eq. 2 |
| ρ (kg/m³) | 8170-8230 | | | | | [A19]; [H718] |
| YS, plate aged (MPa) | 1157 (RT) | 1027 (427 °C) | 1006 (538 °C) | 965 (649 °C) | 723 (760 °C), 236 (871 °C) | [H718] |

- k(900 K) is the end of [A19] eq. 9 (298-800 K) and is slightly extrapolated. [A19] reports
  phase-transition peaks starting at 823 K and 993 K on heating, so c_p near 900 K is effective.
- E: 200 GPa at 21 °C, 179 at 400 °C, 159 at 700 °C [H718]. ν 0.294 (21 °C), 0.271 (427 °C),
  0.283 (649 °C), 0.321 (816 °C) [SMC] T4.
- α mean from 25 °C: 12.8 (to 100 °C), 14.1 (400 °C), 15.5 (700 °C) x10⁻⁶/K [H718].
- Limits: melting 1260-1336 °C [SMC]. Aged at 718 °C then 621 °C [SMC]. Use to 704 °C (1300 °F)
  [SMC]. Yield halves between 649 and 760 °C [H718].

**Haynes 230 (solid solution)**

| property | 300 K | 500 K | 700 K | 900 K | 1100 K | source |
|---|---|---|---|---|---|---|
| k (W/m·K) | 8.9 | 12.4-14.4 (200-300 °C) | 16.4 (400 °C) | 20.4-22.4 (600-700 °C) | 24.4-26.4 (800-900 °C) | [H230] |
| c_p (J/kg·K) | 397 | 435-448 | 465 | 486-574 | 595-609 | [H230] metric column |
| ρ (kg/m³) | 8970 | | | | | [H230] |
| YS, plate (MPa) | 383 (21 °C) | none | 263 (538 °C) | 267 (649 °C), 260 (760 °C) | 234 (871 °C), 116 (982 °C) | [H230] |

- The [H230] c_p table has errors: "5595" at 800 °C, and the metric and °F columns disagree
  between 500 and 650 °C (473 vs 561 J/kg·K). Treat c_p above 500 °C as not verified.
- E: 209 GPa RT, 186 GPa (400 °C), 168 GPa (700 °C). ν 0.31-0.34 over 20-700 °C [H230].
- α mean from 25 °C: 11.8 (to 100 °C), 13.2 (400 °C), 14.7 (700 °C) x10⁻⁶/K [H230].
- Limits: melting 1301-1371 °C; solution treated at 1177-1246 °C; no aging. No sigma or mu
  phase after 16,000 h at 649-871 °C [H230].

## 6. Effusivity and pulse heating

Derived from the tables. Constant flux `q` for time `t` raises the surface by
`ΔT = 2 q sqrt(t/π) / e`. Diffusion depth `sqrt(κ t)` at 1 ms is also shown.

| material | e at 300 K (W·s^½/m²·K) | ratio | e at 700 K | ratio | depth at 1 ms |
|---|---|---|---|---|---|
| maraging 300 | 8,790 | 1.00 | 10,150 | 1.00 | 76 µm |
| GRCop-84 | 30,700 | 3.49 | 33,400 | 3.29 | 295 µm |
| CuCrZr | 32,800-34,600 | 3.73-3.94 | not verified | - | ~300 µm |
| Inconel 718 | 5,880-6,390 | 0.67-0.73 | 8,320-8,660 | 0.82-0.85 | 56 µm |
| Haynes 230 | 5,630 | 0.64 | 8,420 | 0.83 | 50 µm |

With [K23] c_p for maraging, its e rises 6-12% and every ratio falls by that much.
A layer thinner than its diffusion depth behaves partly like the substrate beneath it.

## 7. Cyclic behaviour, bonding, and the layer under GPa pulses

Published evidence found:
- OFHC copper rocket liners failed by "ratcheting-enhanced creep rupture", with progressive wall
  thinning each cycle. Thinning reached 90% at the failure site, and lives were 39-150 cycles at
  2.6-2.9% strain range. Amzirc failed by ordinary low-cycle fatigue [HKP76].
- Porowski et al. model this as plastic ratcheting of the hot wall. For NARloy-Z the thinning
  dies away as the net ratchet increment vanishes [P82].
- CuCrZr: "there have been till now no available research on the ratcheting behavior" [Z19].
- HIP diffusion bonding of CuCr1Zr to 316L steel (950 °C, 105 MPa): joints broke in the CuCr1Zr
  at 360 ± 2 MPa, not at the interface. Interface conductivity matched the weaker bulk [P20].
- GRCop-84 has been diffusion bonded, brazed, inertia welded and EB welded [E05]. VPS NiAl/Ni
  coatings on GRCop-84 cool down in compression and were not expected to spall [N05].
- Explosively welded Cu/steel: interface shear of 222-260 MPa appears in publisher abstracts
  (Coatings 13, 600, 2023; Indian Weld. J. 2018). Full texts were blocked: not verified.
- Cold-sprayed Cu on steel, and any Ni-alloy layer on steel under cyclic compression: no primary
  data found. Note that ASTM C633 pull tests are capped by the epoxy (not verified here).

Derived estimates (formula from [BF91] via the steel note; static Y, so rate uplift is ignored):
uniaxial-strain HEL `= Y (1-ν)/(1-2ν)`, and repeated pulses shake down only if `P <= 2 HEL`.

| layer at room temperature | Y (MPa) | ν | HEL (GPa) | 2 HEL (GPa) | vs a 2.5 GPa pulse |
|---|---|---|---|---|---|
| GRCop-84 | 196 | 0.34 (not verified) | 0.40 | 0.8 | yields both ways every pulse |
| CuCrZr, hard drawn / brazed | 550 / 106 | 0.34 (not verified) | 1.13 / 0.22 | 2.3 / 0.4 | yields every pulse |
| Haynes 230 | 383 | 0.31 | 0.70 | 1.4 | yields every pulse |
| Inconel 718, aged | 1157 | 0.294 | 1.98 | 4.0 | shakes down; elastic after pulse 1 |
| maraging 300 | 1790-2070 | 0.30 | 3.1-3.6 | 6.2-7.2 | elastic |

- A bonded layer is held in-plane by the steel, so its flow is reversed plasticity in uniaxial
  strain, not free ratcheting. That is low-cycle fatigue, 1,500 cycles per push.
- The bulk heat-up adds an in-plane misfit of `Δα ΔT`. Over a 300 K rise this is about 0.2-0.3%
  for copper on maraging (Δα ≈ 7-9 x10⁻⁶/K) and about 0.1% for the superalloys (derived). Copper
  yields at ~0.1-0.4% strain, so the misfit alone yields it once per heat-up.
- Pulse pressure does not fall through any layer (steel note, Section 1). The maraging under a
  layer still sees the full `P`.

## 8. Comparison

- **Thermally, copper helps.** GRCop-84 or CuCrZr has 3.3-3.9x the effusivity of maraging, so a
  layer at least ~0.3 mm thick cuts the per-pulse face rise to about 30% of bare steel. It also
  spreads hot spots sideways. GRCop-84 is preferred because it keeps strength after a bonding
  cycle and to 700 °C. CuCrZr loses ~80% of its strength after a 935 °C bond and over-ages above
  500 °C. Its usual ~460-475 °C age is close to maraging's 480 °C age, so a shared final age may
  be possible (inference, not verified). The cost is that copper yields on every 2.5 GPa pulse.
  Over 1,500 pulses that is a low-cycle-fatigue and thinning risk, as in the rocket-liner data.
- **Structurally, Inconel 718 helps.** It is the only layer that stays elastic at 2.5 GPa after
  shakedown. It keeps 965 MPa yield at 649 °C, against maraging's 480 °C aging limit. It is a worse
  heat sink, though: the face runs about 1.2-1.6x hotter than bare maraging for the same flux. It
  helps only if the extra temperature stays in the Inconel and the maraging below stays cooler
  than 480 °C. Haynes 230 has the same thermal penalty but yields at 0.7 GPa, so it brings no
  structural gain at these pressures.
- For the solver: use the Section 2-5 rows with a source tag. Mark CuCrZr k above 573 K, all
  maraging values above 753 K, and GRCop-84 E and ν as gaps.

## Sources

- [NI76] Nickel Development Institute, "18 per cent nickel maraging steels: engineering properties", Publication 4419, Tables 7, 14, 16, 19. Read via https://web.archive.org/web/2022id_/https://nickelinstitute.org/media/1598/18_nickelmaragingsteel_engineeringproperties_4419_.pdf
- [K23] P. Koniorczyk et al., "Experimental studies of thermophysical properties and microstructure of X37CrMoV5-1 hot-work tool steel and maraging 350 steel", Materials 16, 1206 (2023). https://doi.org/10.3390/ma16031206 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC9920995/
- [E00] D. L. Ellis, D. J. Keller, M. Nathal, "Thermophysical properties of GRCop-84", NASA/CR-2000-210055. https://ntrs.nasa.gov/api/citations/20000064095/downloads/20000064095.pdf
- [E02] D. L. Ellis, "Aerospace Structural Materials Handbook supplement: GRCop-84" (2002). https://ntrs.nasa.gov/api/citations/20020070630/downloads/20020070630.pdf
- [E05] D. L. Ellis, "GRCop-84: a high-temperature copper alloy for high-heat-flux applications", NASA/TM-2005-213566. https://ntrs.nasa.gov/api/citations/20050123582/downloads/20050123582.pdf
- [D07] H. C. de Groh III, D. L. Ellis, W. S. Loewenthal, "Comparison of GRCop-84 to other high thermal conductive Cu alloys", NASA/TM-2007-214663, Tables 4, 7, 12, 13. https://ntrs.nasa.gov/api/citations/20070017311/downloads/20070017311.pdf
- [LZ12] M. Li, S. J. Zinkle, "Physical and mechanical properties of copper and copper alloys", Comprehensive Nuclear Materials 4.20 (2012), Table 1. https://digital.library.unt.edu/ark:/67531/metadc829254/m2/1/high_res_d/1036615.pdf
- [Z19] K. Zhang et al., "Development of the material property handbook and database of CuCrZr", Fusion Eng. Des. 144, 148 (2019). https://scientific-publications.ukaea.uk/wp-content/uploads/DEVELOPMENT-OF-THE-MATERIAL-PROPERTY-HANDBOOK-AND-DATABASE-OF-CUCRZR.PDF
- [P20] S. Pianese et al., "HIP assisted diffusion bonding between CuCr1Zr and AISI 316L ... SPS internal beam dump at CERN", arXiv:2011.07942. https://arxiv.org/abs/2011.07942
- [NIST] NIST Chemistry WebBook, copper solid phase heat capacity (Shomate, Chase 1998 JANAF). https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440508&Mask=2
- [SMC] Special Metals, "INCONEL alloy 718", SMC-045 (2007), Tables 2, 4, 5, 19. https://www.specialmetals.com/documents/technical-bulletins/inconel/inconel-alloy-718.pdf
- [H718] Haynes International, "HAYNES 718 alloy" brochure. https://haynesintl.com/wp-content/uploads/2023/09/718-brochure.pdf
- [A19] A. Sh. Agazhanov et al., "Thermophysical properties of Inconel 718 alloy", J. Phys.: Conf. Ser. 1382, 012175 (2019), Table 2. https://iopscience.iop.org/article/10.1088/1742-6596/1382/1/012175/pdf
- [H230] Haynes International, "HAYNES 230 alloy" brochure. https://haynesintl.com/wp-content/uploads/2023/09/230-brochure.pdf
- [HKP76] N. P. Hannum, H. J. Kasper, A. J. Pavli, "Experimental and theoretical investigation of fatigue life in reusable rocket thrust chambers", NASA TM X-73413 (1976). https://ntrs.nasa.gov/api/citations/19760019182/downloads/19760019182.pdf
- [P82] J. S. Porowski et al., "Development of a simplified procedure for thrust chamber life prediction", NASA CR-165585 (1982). https://ntrs.nasa.gov/api/citations/19820013379/downloads/19820013379.pdf
- [N05] NASA Glenn, "Coating development for GRCop-84 liners for reusable launch vehicles aided by modeling studies" (2005). https://ntrs.nasa.gov/api/citations/20050192259/downloads/20050192259.pdf
- [BF91] Ballard et al. 1991, as cited in `spray_plate_steel_face_limits.md`.

## Blocked sources

No host returned the sandbox's "Approval required" or "Blocked by network policy" responses.
Publisher-side walls met:
- mdpi.com: Akamai "Access Denied" 403 (explosive-weld Cu/steel paper); Wayback was rate-limited.
- sciencedirect.com: 403 bot wall (CuCrZr divertor heat-sink paper).
- juser.fz-juelich.de: JavaScript challenge (Pintsuk 2010 CuCrZr interlaboratory record).
- nickelinstitute.org: Cloudflare challenge (read via Wayback). web.archive.org gave transient 429s.
- MMPDS / MIL-HDBK-5: no open copy sought further (everyspec.com returned 403 in the earlier note).
