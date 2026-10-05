# Spray plate: steel face limits, prestress, and alternative face materials

Status: research note, 2026-10-05. Serves [ADR-0055](adr/0055-argon-arm-of-the-spray-plate.md) and
[`argon_plate_handoff.md`](argon_plate_handoff.md). Load case: 0.3-7 GPa sharp-fronted pulses of
0.1-1 ms, ~1,500 per push at 4 Hz, on a tapered steel plate of ~4 cm mean thickness and 20 m
loaded diameter. Anything not traced to a primary source is marked "not verified".

## 1. Bottom line

- The prestress mechanics are **correct**. Under uniaxial strain with an equibiaxial in-plane
  compression `p` (0 <= p <= Y), the face stays elastic up to `(Y + p)(1-ν)/(1-2ν)`. The ceiling
  is 2 HEL, reached at `p = Y`. Ballard et al. (1991) give the same 2 HEL bound from the
  laser-peening side (Section 3).
- The **delivery** of that prestress is the problem, not the mechanics. A rim band cannot do it:
  a 20 m, 4 cm plate buckles under ~1-4 MPa of membrane compression, about 1000x too little.
  A self-equilibrated face/backing pair works only for a face layer, and costs the backing margin,
  because a flat-topped pulse longer than the 7 µs transit loads the whole thickness to `P`.
- No face layer lowers the pressure the steel behind it sees. For a pulse much longer than the
  layer transit times, every layer equilibrates to about `P`. SiC tiles or aramid change the
  thermal and erosion picture, not the steel's pressure limit.
- Recommended working allowable for a prestressed hard-steel face: **3.5 GPa** (maraging 300 or
  350, `p ≈ 0.5 Y`, face ΔT <= 300 K). Without prestress: **2.5 GPa** for maraging, **1.5 GPa**
  for HY-100 or 4340 at Rc 38-45. The unmerged 4.7-7 GPa pulses fail every option here.

## 2. Steel data: quasi-static Y, HEL, dynamic yield

`F(ν) = (1-ν)/(1-2ν)`: 1.54 at ν = 0.26, 1.67 at 0.286, 1.69 at 0.29, 1.75 at 0.30.

| steel | quasi-static Y (GPa) | measured HEL (GPa) | condition | source |
|---|---|---|---|---|
| HY-100 | 0.69 min (spec; not verified) | **1.89 ± 0.01** (4 shots, 1.87-1.90) | symmetric impact 0.198-0.394 km/s, 3.5-6.2 GPa peak; ν = 0.286, K = 163.85 GPa, G = 81.79 GPa | Thomas et al. 2018 [T18] |
| Vascomax 350 (18Ni2400) | 2.39 typical, ν 0.26, E 191-199 GPa | **4.8 ± 2.0** (6.4 mm sample) | gas gun / explosive, 180-1000 kbar series | Gust and Royce 1970 [GR70]; Y from [NI76] Table 7 |
| maraging 300 (18Ni1900) | 1.79-2.07, ν 0.30, E 190 GPa | not found | none | [NI76] Table 7 |
| 4340, Rc 35 | 1.03 (150 ksi), ν 0.30 | steels in general "5 to 40 kbar depending on hardness" | explosive lens / gas gun | Weirick 1994 [W94] |
| 4340, Rc 30 / 38 / 45 / 50 | 0.73 / 0.97 / 1.25 / 1.49 | not found per temper | fit `Y = exp(0.0355 Rc + 5.5312)` MPa to ASM data | Banerjee [B05] |
| AF9628 (tempered martensitic low-alloy) | (abstract: dynamic Y "only slightly higher" than quasi-static) | **2.35** | symmetric gun impacts | Neel et al. 2019 [N19], abstract only |
| HY-80 | 0.55 min (spec; not verified) | not found | none | none |
| 2.25Cr-1Mo, 4130 | ~0.2-0.3 for 2.25Cr-1Mo (not verified) | not found | none | none |

Dynamic yield and rate effects:
- HY-100: `Y_HEL = 1.89 / 1.668 = 1.13 GPa`, about 1.6x the 0.69 GPa spec minimum. This is
  the usual rate uplift of a lower-strength bcc steel at ~10^5-10^6 /s.
- 4340: the Johnson-Cook rate term used by Banerjee [B05] is `C = 0.014` (ε̇0 = 1 /s). That gives
  1.16x at 10^5 /s and 1.19x at 10^6 /s, applied to all tempers. The original Johnson-Cook
  constants (Rc 30, A = 792 MPa) are from Johnson and Cook 1983 [JC83], not verified here.
- Maraging and AF9628: hard martensites show little rate uplift. [N19] says so for AF9628. The
  Vascomax value implies `Y_HEL = 4.8 / 1.542 = 3.1 ± 1.3 GPa`, consistent with 2.39 static
  within its large error bar.
- Correction to the brief: Rohde 1969 [R69] is on **iron** (76-573 K, ~10^5 /s), not 4340.
  The LASL Shock Hugoniot Data compilation (Marsh 1980) was not retrieved.
- HEL in steels drops with sample thickness (precursor decay). Gust and Royce note thick samples
  give lower values [GR70]. A 4 cm plate sits at the low end of any thin-sample value.

## 3. The prestress relation, checked

Take x normal to the face, compression negative. Uniaxial strain `ε_x = e`: `σ_x = (λ+2μ)e` and
`σ_y = σ_z = -p + λe`. Von Mises with `σ_y = σ_z` reduces to `|σ_x - σ_y| = Y`. At e = 0 the
difference is `+p`. It falls as `2μe` and reaches `-Y` when `2μ|e| = Y + p`. Then
`|σ_x| = (λ+2μ)(Y+p)/(2μ) = (Y+p)(1-ν)/(1-2ν)`. The prestressed state itself needs `p <= Y`
(equibiaxial von Mises). So the limit runs from 1 HEL (p = 0) to 2 HEL (p = Y). The relation in
the brief holds exactly for a perfectly plastic von Mises solid.

Primary cross-check: Ballard, Fournier, Fabbro and Frelat [BF91] write the HEL as
`(1 + λ/2μ) σ_Y`, which equals `Y(1-ν)/(1-2ν)`. Their single-pulse plastic strain is zero below
HEL, rises linearly to `-2(1+λ/2μ)σ_Y/(3λ+2μ)` at `P = 2 HEL`, and saturates there. The optimum
peening pressure is `2 HEL`, and the surface residual stress tends to `-σ_Y` for a large spot.
That is the same physics read backwards: one shock at 2 HEL builds `p = Y`, after which a second
identical shock is elastic. Their 35CD4 steel (50 HRC, σ_Y = 1250 MPa) measured about
-400 MPa residual over a 5 mm spot, matching their model.

Achievable elastic limits (static Y; multiply by 1.1-1.6 for rate uplift where Section 2 allows):

| steel | Y (GPa) | HEL (GPa) | p = 0.5 Y | p = 0.9 Y |
|---|---|---|---|---|
| maraging 350 | 2.39 | 3.7 calc, 4.8 ± 2.0 meas. | 5.5 | 7.0 |
| maraging 300 | 1.79-2.07 | 3.1-3.6 calc | 4.7-5.4 | 6.0-6.9 |
| 4340 Rc 50 | 1.49 | 2.5 calc | 3.8 | 4.8 |
| 4340 Rc 45 | 1.25 | 2.1 calc | 3.2 | 4.0 |
| 4340 Rc 38 | 0.97 | 1.6 calc | 2.5 | 3.1 |
| HY-100 | 1.13 (from HEL) | 1.89 meas. | 2.8 | 3.6 |
| 2.25Cr-1Mo | ~0.2-0.3 (not verified) | ~0.4-0.5 calc | <= 0.8 | <= 1.0 |

Caveats on delivering `p`:
1. **Rim band.** Net membrane compression on a thin disk buckles it. For a circular plate under
   uniform radial compression, `N_cr = k D / a²` with k ≈ 14.7 (clamped) or 4.2 (simply
   supported, ν = 0.3) (Timoshenko and Gere, Theory of Elastic Stability; coefficients not
   verified against the text). With `D = E h³/12(1-ν²) = 1.17 MN·m`, `a = 10 m`, `h = 4 cm`:
   σ_cr ≈ 1.2-4.3 MPa. A rim band can supply at most a few MPa, against ~1 GPa needed.
2. **Face layer on a tensioned backing.** Net force is zero, so there is no global buckling.
   But the backing carries `+q = p t_face / t_back` and sees the same `P` on the first transit.
   Its limit drops to `(Y - q) F(ν)`. A 5 mm face at `p = 0.8 Y` on 35 mm puts q = 0.11 Y in the
   backing. The eccentric layer also curls the plate unless balanced.
3. **Bonded internal tendons.** A member prestressed by bonded tendons is not pushed into
   buckling by its own prestress (a standard prestressed-concrete result; not verified here).
   Concentric fibre tendons could hold the whole steel section in compression. Arrangement and
   wave reflections at the tendon layers are not analysed.
4. **Self-generated prestress does not rescue a thin plate.** In [BF91] the plastic layer is
   restrained by elastic material around and below it. Here the whole 4 cm section yields in
   uniaxial strain, so pulse 1 makes the footprint grow laterally rather than leaving residual
   compression. The surrounding plate resists that growth only by membrane compression, which
   buckles at MPa (item 1). This is an inference, not verified.
5. **Uniaxial strain lasts only so long.** [BF91] require `τ << r0 √(ρ(λ+2μ)/(4μ(λ+μ)))`.
   For a 5 m footprint that is of order 1 ms. So 0.1-0.3 ms pulses are uniaxial; 1 ms pulses
   partly relax toward biaxial bulk yield at the footprint edge.

## 4. Repeated pulses above the limit: ratcheting and fatigue

- Per-pulse plastic strain above HEL in uniaxial strain is `|ε_p| = 2(P - HEL)/(3K)` [BF91].
  For maraging 300 (HEL ≈ 3.4 GPa, K ≈ 164 GPa [T18]) at P = 4.7 GPa this is 0.5% through-
  thickness strain per pulse. If nothing shakes it down (Section 3, item 4), 1,500 pulses give
  an accumulated strain of about 8 (800%). Any `P > limit` design is therefore excluded, not merely penalised.
- Shakedown bound: the elastic range of the deviator is 2Y, so repeated uniaxial-strain pulses
  can shake down only for `P <= 2 HEL`, and only if `p = Y` is reacted structurally. Above
  2 HEL, reverse yield on unloading gives alternating plasticity every pulse. This is the
  [BF91] saturation read as a cyclic limit.
- Fatigue under compression-dominant cycling: the [BF91] notched 35CD4 tests show compressive
  residual stress raises fatigue strength (680 MPa reference, 720 shot-peened at -650 MPa,
  900 laser-shocked at -550 MPa). That is a benefit only while the face stays elastic.
- Heat checking (in-plane compressive yield of a hot skin, then tension on cooling) is the
  likely crack driver for pulses below HEL. No primary multi-pulse data for steels in this
  regime were found. Orion pusher reports (GA-5009 series) were not located on OSTI or NTRS:
  **not verified**.

## 5. Thermal interaction

In-plane thermal stress in a heated skin, fully restrained: `σ_th = E α ΔT / (1-ν)`.
Diffusion depth over 0.5 ms is `√(κt)` ≈ 50 µm (maraging, k = 21 W/m·K, ρc = 8.0 x 0.46 MJ/m³·K
[NI76]). So the thermal stress lives in a ~0.05-0.1 mm skin.

| material | E (GPa) | α (10^-6/K) | ν | σ_th at ΔT 100 / 300 K (MPa) | basis |
|---|---|---|---|---|---|
| maraging 300 | 190 | 9.9 (20-100 C) to 11.3 (20-480 C) | 0.30 | 277 / 830 | [NI76] Tables 7 and 19 |
| maraging 350 | 191-199 | 11.4 (20-480 C) | 0.26 | ~290 / ~880 | [NI76] |
| HY-100 | 210 (from K, G of [T18]) | ~12 (not verified) | 0.286 | 353 / 1060 | [T18] + not verified |
| 4340 | ~205 (not verified) | ~12.3 (not verified) | 0.29 | 355 / 1065 | not verified |
| SiC | 434 (from ρ, c_L, c_S of [G93]) | ~4 (not verified) | 0.16 | 207 / 620 | [G93] + not verified |

The budget rule: the skin's in-plane compression is `p + σ_th`, and it must stay `<= Y(T)`.
- While `p + σ_th <= Y`, heating does not lower the normal-stress limit. It raises it, since the
  limit is `(Y + p + σ_th) F(ν)`.
- When `p + σ_th > Y`, the skin yields in compression during heating. On cooling it relaxes
  to less compression or to tension. The prestress is lost and every pulse is a thermal cycle.
- So `p_max = Y(T) - σ_th`. For maraging 300 (Y 1.79 GPa) at ΔT = 300 K, `p <= 0.96 GPa ≈ 0.5 Y`.
  For HY-100 or 4340 at ΔT = 300 K, σ_th is 0.85-0.94 of Y, leaving almost nothing.
- **Bulk temperature.** ADR-0055 budgets a ~500 K bulk rise. Maraging is aged at 480 C [NI76],
  and the Rc 45 4340 temper is 425 C ([B05], citing Chi et al. 1989). A plate reaching
  ~500 C would over-age or over-temper and soften. That rise budget is incompatible with these
  tempers (inference; softening rates not verified).

## 6. Alternative face materials

**Aramid (Kevlar 49) / epoxy.** Jones 1983 [J83], UD Kevlar 49/epoxy:
- Longitudinal tension 1,800 MPa at Vf 0.53. Longitudinal compression **230 MPa** (yield, Vf 0.60).
- Transverse (through-thickness proxy) compression **140 MPa** (Vf 0.48). In-plane shear 73 MPa.
- In [J83]'s words, the "poor compressive strength of Kevlar composites renders them unsuitable for
  any applications where significant compressive or bending stresses exist".
- Confined (uniaxial-strain) through-thickness behaviour at GPa was not found: not verified.
- Thermal: decomposition 427-482 C in air; continuous use to ~200 C; strength losses
  noted above ~150-177 C (DuPont Kevlar Technical Guide, quoted from a search snippet; full text
  blocked; **not verified**). That is far below the steel's 500 K bulk-rise budget.
- Verdict: the expectation is **confirmed**. Aramid is a tension member, such as a rim band or
  bonded tendons that prestress the steel (Section 3, item 3). It is not a GPa face. Under an
  ablative it would still have to stay below ~150-200 C over 1,500 pulses. Its low impedance
  (ρ ≈ 1.4 g/cm³) gives it no shielding, since a long pulse equilibrates every layer to `P`.

**SiC (or B4C) tiles on steel.**
- HEL: Grady 1993 [G93] gives dynamic yield 12.5 GPa for an Eagle Picher SiC (ρ 3177 kg/m³,
  c_L 12.06, c_S 7.67 km/s). That converts to **HEL ≈ 15.5 GPa**. B4C gives 13.7-15.1 GPa
  dynamic yield, or HEL ≈ 17-19 GPa. Gust and Royce [GR70] measured 9.1 ± 1.5 GPa for a lower
  grade (Carborundum KT) SiC and 13.7-16.2 GPa for B4C, falling with thickness.
- Spall: SiC-B **0.90 GPa** at 1.6 GPa impact, peaking at **1.3 GPa** at 3.75 GPa, then
  0.82-0.90 GPa by 12 GPa. Other SiCs run 0.5-1.85 GPa. Spall strength declines once shocked
  past ~5.8 GPa [DB01]. Grady [G93] gives 0.2-0.5 GPa as typical for ceramics (OCR-garbled
  text) and near-zero spall after shocks above HEL in alumina.
- Thermal shock: no primary ΔT_c was retrieved. The Kingery estimate `R = σ_f(1-ν)/(Eα)` with
  σ_f ≈ 0.4 GPa (not verified) gives ~190 K. Face heating puts the hot skin in compression,
  which SiC tolerates well. The risk is the tensile subsurface and cool-down, and the 4 vs
  10-12 x 10^-6/K CTE mismatch to steel over a 500 K bulk rise.
- Impedance check of the coordinator's claim, using [G93] and [T18] data:
  Z_SiC = 3177 x 12,060 = 3.83e7 Pa·s/m. Z_steel elastic = 7840 x 5,900 = 4.63e7, giving
  `R = +0.094`. **Compression reflects as compression, confirmed**: the reverse of ADR-0011's SiC-on-Ti
  (R ≈ -0.15). Caveat: once the steel goes plastic, its shock impedance is about ρU_s =
  7840 x 4,760 = 3.73e7 [T18], so `R ≈ -0.01`, essentially zero. The interface stays benign.
- Verdict: SiC survives the pressure itself to ~10 GPa. But the steel behind still sees `P`.
  SiC on steel is a thermal and erosion face, not a fix for the steel's HEL.

## 7. Material comparison

| row | steel, plain | steel, prestressed | aramid/epoxy | SiC tile on steel |
|---|---|---|---|---|
| compressive / HEL limit | 1.9 (HY-100) to ~3.4-4.8 GPa (maraging) | up to 2 HEL; ~4.7-5.5 GPa maraging at p = 0.5 Y | 0.23 GPa long., 0.14 GPa transverse (unconfined) | SiC HEL 9-15.5 GPa; steel backing still limits |
| ductility | ductile; ratchets above HEL | same, above raised limit | fibre kinking; brittle-ish matrix | brittle; spall 0.5-1.3 GPa |
| thermal limit | temper/aging 425-480 C | same, and σ_th consumes `p` | ~150-200 C use, 427-482 C decomposition (not verified) | >1000 C; ΔT_c ~200 K estimate; CTE mismatch |
| cyclic behaviour | shakedown impossible in thin plate if P > HEL | elastic if P <= limit and `p + σ_th <= Y` | creep and stress rupture in tension (not verified) | spall loss after >HEL shocks [G93, DB01] |
| suitable role | backing / structure, P <= ~1.5-2.5 GPa | face, P <= ~3.5 GPa | rim band or bonded tendons (tension only) | thermal/erosion face over steel |

## 8. Recommendation

| steel | Y (GPa) | HEL (GPa) | prestressed elastic limit (GPa) | basis |
|---|---|---|---|---|
| maraging 350 (18Ni2400) | 2.39 | 3.7 calc; 4.8 ± 2.0 meas. | 5.5 (p = 0.5 Y) to 7.0 (0.9 Y) | [NI76], [GR70], Section 3 |
| maraging 300 (18Ni1900) | 1.79-2.07 | 3.1-3.6 calc | 4.7-5.4 (0.5 Y) to 6.0-6.9 (0.9 Y) | [NI76], Section 3 |
| 4340 Rc 45-50 | 1.25-1.49 | 2.1-2.5 calc (+16-19% rate) | 3.2-3.8 (0.5 Y) to 4.0-4.8 (0.9 Y) | [B05], [W94] |
| HY-100 | 0.69 spec; 1.13 dynamic | 1.89 meas. | 2.8 (0.5 Y) to 3.6 (0.9 Y) | [T18] |
| HY-80 | 0.55 spec (not verified) | ~1.5 est. (not verified) | ~2.2-2.8 (not verified) | scaled from HY-100 |
| 2.25Cr-1Mo | ~0.2-0.3 (not verified) | ~0.4-0.5 calc | <= 1.0 | unsuitable |

**Working allowable, prestressed hard-steel face: 3.5 GPa peak normal stress.** Basis: maraging
300 at its lower-bound Y (1.79 GPa); `p = 0.5 Y`, the most the ΔT = 300 K thermal budget allows
(Section 5); limit 4.7 GPa; knocked down by 0.75. The knock-down covers HEL scatter (±40% on the
only measured maraging value), thickness decay, and the unverified delivery of `p`. Without
structurally reacted prestress, use **2.5 GPa** (maraging) or **1.5 GPa** (HY-100, 4340 Rc 38-45).
All of these hold only while the bulk stays below the aging or tempering temperature. That means
well under ADR-0055's 500 K rise.

Against the handoff's loads: the deep-cloud cases (0.34-0.77 GPa) and step 2b's 0.64-1.8 GPa
pass on plain steel. Step 2b's 2.7 GPa needs prestressed maraging. The unmerged 4.7-7 GPa pulses
exceed every allowable here and would ratchet ~0.5% per pulse.

## Sources

- [T18] S. A. Thomas, M. C. Hawkins, M. K. Matthes, G. T. Gray III, R. S. Hixson, "Dynamic strength properties and alpha-phase shock Hugoniot of iron and steel", J. Appl. Phys. 123, 175902 (2018). https://doi.org/10.1063/1.5019484 ; open copy https://www.osti.gov/servlets/purl/1735933
- [GR70] W. H. Gust, E. B. Royce, "Dynamic yield strengths of light armor materials", UCRL-50901 (1970). https://doi.org/10.2172/4119679 ; https://www.osti.gov/servlets/purl/4119679
- [NI76] Nickel Development Institute, "18 per cent nickel maraging steels: engineering properties", Publication 4419 (1976), Tables 7 and 19. https://nickelinstitute.org/media/1598/18_nickelmaragingsteel_engineeringproperties_4419_.pdf (read via https://web.archive.org/web/2022id_/https://nickelinstitute.org/media/1598/18_nickelmaragingsteel_engineeringproperties_4419_.pdf)
- [W94] L. J. Weirick, "Plane shock generator explosive lens: shock characterization of 4340 and PH13-8Mo steels, C360 brass and PZT 65/35", SAND93-3919 (1994). https://www.osti.gov/servlets/purl/10142257
- [B05] B. Banerjee, "The Mechanical Threshold Stress model for various tempers of AISI 4340 steel", arXiv:cond-mat/0510330. https://arxiv.org/abs/cond-mat/0510330 (secondary; used for its temper/yield fit and JC rate constant).
- [JC83] G. R. Johnson, W. H. Cook, Proc. 7th Int. Symp. Ballistics, The Hague (1983) 541-547. Not retrieved; not verified.
- [N19] C. Neel et al., "Shock and spall in the low-alloy steel AF9628", J. Dyn. Behav. Mater. 6, 64 (2020). https://doi.org/10.1007/s40870-019-00228-5 (abstract only, via https://ui.adsabs.harvard.edu/abs/2019JBDM....6...64N/abstract).
- [R69] R. W. Rohde, "Dynamic yield behavior of shock-loaded iron from 76 to 573 K", Acta Metall. 17, 353 (1969). https://doi.org/10.1016/0001-6160(69)90075-3 (cited for scope only).
- [BF91] P. Ballard, J. Fournier, R. Fabbro, J. Frelat, "Residual stresses induced by laser-shocks", J. Phys. IV 01 (C3), C3-487 (1991). https://doi.org/10.1051/jp4:1991369 ; https://hal.science/jpa-00250513v1
- [G93] D. E. Grady, "Dynamic properties of ceramic materials", SAND93-0610 (1993/1995), Table 3. https://www.osti.gov/servlets/purl/10187138
- [DB01] D. P. Dandekar, P. T. Bartkowski, "Tensile strengths of silicon carbide (SiC) under shock loading", ARL-TR-2430 (2001). https://www.govinfo.gov/content/pkg/GOVPUB-D101-PURL-gpo22081/pdf/GOVPUB-D101-PURL-gpo22081.pdf
- [J83] M. L. C. Jones, "The basic ply properties of a Kevlar 49/epoxy resin composite system", RARDE Report 11/83 (1983), DTIC AD-A136614, Table 16. https://web.archive.org/web/20250904214118/https://apps.dtic.mil/sti/tr/pdf/ADA136614.pdf

## Blocked sources

No host returned the sandbox's "Approval required" or "Blocked by network policy" responses.
These publisher-side walls were met, and open copies or the Wayback Machine were used where possible:
- apps.dtic.mil: "Service unavailable, the request is blocked" (site CDN). ADA136614 read via
  Wayback. ADA428178 (Dandekar 2004 SiC spall) had no snapshot; ARL-TR-2430 used instead.
- link.springer.com: login redirect (AF9628 paper; abstract used).
- everyspec.com: HTTP 403 (MIL-HDBK-5J; steel E and α for 4340/HY-100 therefore not verified).
- nickelinstitute.org: Cloudflare challenge (read via Wayback).
- ptacts.uspto.gov: CloudFront 403 (DuPont Kevlar guide; thermal figures not verified).
- srd.nist.gov/JPCRD/jpcrd522.pdf returned an unrelated paper; SiC α and σ_f not verified.
