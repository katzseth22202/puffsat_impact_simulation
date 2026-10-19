# How hot can the Kapton face run for a sub-second pulse?

Research note, 2026-10-19. The question: what raises the short-term temperature tolerance of the
100 um Kapton HN growth-wave sheet's plume face, and by how much.

The model is `python/puffsat/xray_collision/coated_sheet.py`:

- `T_RATING = 673.15 K`: the 400 C continuous rating. The face must stay under it until the last
  10 ms (`LAST`).
- `T_FAIL = 800 K`: "rapid decomposition", marked in the code as assumed and not verified.
- `KAPTON_RHO 1420`, `KAPTON_CP 1090`, `KAPTON_K 0.12`. Both faces have an emissivity of 0.8.

**Short answer.** Both limits come from long-exposure data. Neither is a sub-second decomposition
limit.

- Kapton HN's own kinetics come from an ESA TGA study. They put the face that loses 1 um of mass
  in a 0.1-0.5 s pulse at about 880-1100 C (1150-1370 K). Losing 0.1 um takes about 760-970 C.
  Both are far above 673 K and 800 K.
- That number is an extrapolation. The kinetics were fitted at 550-610 C and run here up to
  1000 C. Each factor of 10 in the true rate moves the threshold by about 100 K.
- The real limits below about 900 C are not decomposition. They are softening above Tg
  (360-410 C), particle erosion of the top ~0.1-1 um, and outgassing.

Labels: **[S]** means sourced (citation given). **[E]** means my own estimate or arithmetic.
The scripts behind every [E] number are listed at the end. They live in the session scratchpad
and are not committed modules. To put any of these numbers in an ADR, port them into a committed
module and Make target first.

---

## 1. Decomposition kinetics of PMDA-ODA (Kapton HN) and what a sub-second pulse allows

### 1.1 The sourced kinetics

The best primary source found is an ESA/ESTEC TGA study of Kapton HN and Upilex S, run under inert
(nitrogen) purge at 1, 2, 5 and 10 C/min. It was checked against isothermal TGA and 350 C
high-vacuum ageing.

> Heltzel, S. & Semprimoschnig, C. O. A. (2003). *A detailed study on the thermal endurance of
> Kapton HN and Upilex S.* Proc. 9th Int. Symp. Materials in a Space Environment, ESA SP-540,
> pp. 179-184. http://esmat.esa.int/Materials_News/SP-540/papers/s4_heltz.pdf

| quantity (10 C/min unless stated) | Kapton HN | Upilex S | |
|---|---|---|---|
| T5 (5 % mass loss) | 608 C | 631 C | [S] |
| Tmax (peak rate) | 631 C | 662 C | [S] |
| residue after full decomposition | 54.6 % | 62.4 % | [S] |
| visible onset of decomposition | ~400 C | ~500 C | [S] |
| E, ASTM E1641 (first order, at 5 % conversion) | 250 kJ/mol | 312 kJ/mol | [S] |
| log10 A [A in min^-1], ASTM E1641 | 13.5 | 16.8 | [S] |
| E, model-free kinetics at 5 % conversion | 290 kJ/mol | 300 kJ/mol | [S] |
| "60-min half-life temperature" (50 % conversion in 1 h) | 574 C (measured 51 %) | 598 C (measured 48 %) | [S] |
| 350 C high-vacuum ageing, 2674 h / 4452 h | 1.1 % mass loss | 0.0 % | [S] |

In that paper, "conversion" means the fraction of the eventual volatile loss. For Kapton that is
45.4 % of the mass, so 5 % conversion is 2.3 % of the mass.

**Check [E].** Take the ASTM set, `k = A exp(-E/RT)` with `A = 10^13.5 / 60 = 5.27e11 s^-1`.
It gives 51.8 % conversion after 60 min at 574 C. The paper measured 51 %. So the set reproduces
its own isothermal check.

Other kinetic sources agree on E but are weaker:

- An isoconversional (FWO / KAS / Starink) study of a Kapton-type PI membrane gave mean E of
  284.6 / 286.1 / 286.5 kJ/mol, best fit second or third order (F2/F3). [S] Li, Q. et al. (2024),
  *The thermal behavior and pyrolysis mechanism of a polyimide gas separation membrane*, RSC
  Sustainability. https://doi.org/10.1039/d4su00489b (abstract and search snippet only; the full
  text was blocked).
- Kapton in an oxidative purge was fitted with a five-step model. The initial competitive steps
  have E = 189 kJ/mol. [S] Hondred, P. R., Bowler, N., Kessler, M. R. (2013), *Electrothermal
  lifetime prediction of polyimide wire insulation with application to aircraft*, J. Appl. Polym.
  Sci. 130, 1639. https://doi.org/10.1002/app.39304. That figure is for air and does not apply in
  vacuum.
- Lua, A. C. & Su, J. (2006), *Isothermal and non-isothermal pyrolysis kinetics of Kapton
  polyimide*, Polym. Degrad. Stab. 91, 144-153,
  https://doi.org/10.1016/j.polymdegradstab.2005.04.021, is the standard Kapton-in-N2 reference.
  It is paywalled; I could not read its numbers, so none are quoted.
- DuPont's summary of properties (H-38492) shows isothermal weight-loss curves in He and air at
  400-600 C (its Fig. 10). The digitised text does not carry the values. DuPont does state that
  life in helium is "at least an order of magnitude greater than in air". [S]
  https://www.electro-wind.com/web-files/DuPont/Literature/Kapton-summary-of-properties.pdf

In hard vacuum there is no oxygen, so the inert-gas kinetics are the right ones. No vacuum-specific
Arrhenius set was found. The ESA 350 C high-vacuum ageing matches the TGA model (Tables 5-6 of the
paper), which supports using the inert set in vacuum. [S]

### 1.2 Isothermal arithmetic at the face [E]

This is first-order conversion of the volatile fraction, `alpha = 1 - exp(-k t)`, holding a
temperature for the whole pulse. It is the upper bound on damage at that temperature.

| T | k [1/s] (ASTM) | time to 1 % conv. | alpha after 0.1 s | alpha after 0.5 s | same, Upilex S, 0.5 s |
|---|---|---|---|---|---|
| 400 C | 2.1e-8 | 5.5 days | 2e-9 | 1e-8 | 3e-10 |
| 500 C | 6.8e-6 | 25 min | 7e-7 | 3.4e-6 | 4.4e-7 |
| 600 C | 5.8e-4 | 17 s | 6e-5 | 2.9e-4 | 1.1e-4 |
| 700 C | 0.020 | 0.5 s | 0.002 | 0.010 | 0.009 |
| 800 C | 0.36 | 28 ms | 0.035 | 0.16 | 0.29 |
| 900 C | 3.9 | 2.6 ms | 0.32 | 0.86 | ~1 |
| 1000 C | 29 | 0.3 ms | 0.95 | 1 | 1 |

What the table shows:

- At the model's `T_FAIL` of 800 K (527 C), a 0.5 s hold converts about 1e-5 of the volatiles.
  That is not "rapid decomposition".
- Rapid means alpha of order 0.1 inside the pulse. That happens near 800 C for 0.5 s and near
  850-900 C for 0.1 s.
- A third set uses E = 290 kJ/mol, with A anchored to the ASTM rate at 580 C (A = 1.5e14 s^-1).
  This is my estimate, since the paper does not print A for the model-free fit. At 800-900 C it
  runs 2-5 times faster than the ASTM set, which shifts the thresholds down by about 30-60 K.

### 1.3 Depth lost through a real temperature profile [E]

A face held at T only matters through the thin layer near the peak. With a gradient `G = q/k` into
the film, the reacting layer is about `(R T^2/E) / G` thick. At 1000 C and E = 250 kJ/mol that is
about 54 K. At `q` = 1 MW/m^2 the gradient is 8.3 K/um, which gives a reacting layer of about 6 um.

A rising face also spends only part of the pulse near its peak. So I ran the case directly:

- 1-D backward-Euler conduction through 100 um of Kapton, with the same `rho`, `cp`, `k` and
  emissivities as `coated_sheet.py`, starting from 300 K.
- A flux pulse of length `t_p`, either flat or a linear ramp, tuned so the face peaks at the
  stated temperature at the end of the pulse.
- First-order Arrhenius conversion integrated in every cell.

"Mass lost" below is `0.454 x integral(alpha dx)`. It is the thickness-equivalent of mass actually
gone; the rest stays as char.

**Face peak temperature at which the mass lost reaches 0.1 / 1 / 10 um:**

| kinetics | pulse | 0.1 um | 1 um | 10 um |
|---|---|---|---|---|
| Kapton HN, ASTM (E 250) | flat 0.01 s | 1082 C | 1254 C | > 1500 C |
| | flat 0.1 s | 912 C | 1027 C | 1256 C |
| | flat 0.5 s | 789 C | 876 C | 984 C |
| | ramp 0.1 s | 972 C | 1103 C | 1468 C |
| | ramp 0.5 s | 846 C | 945 C | 1080 C |
| Kapton HN, E 290 anchored at 580 C (estimate) | flat 0.01 s | 996 C | 1131 C | > 1500 C |
| | flat 0.1 s | 864 C | 956 C | 1156 C |
| | flat 0.5 s | 764 C | 836 C | 926 C |
| | ramp 0.1 s | 912 C | 1016 C | 1343 C |
| | ramp 0.5 s | 811 C | 892 C | 1005 C |
| Upilex S, ASTM (E 312) | flat 0.01 s | 1006 C | 1135 C | > 1500 C |
| | flat 0.1 s | 881 C | 969 C | 1167 C |
| | flat 0.5 s | 784 C | 853 C | 940 C |
| | ramp 0.1 s | 927 C | 1025 C | 1353 C |
| | ramp 0.5 s | 830 C | 907 C | 1017 C |

The doses these runs need to reach their peaks are close to the sheet's own 40-100 kJ/m^2. For
example, a flat 0.1 s pulse peaking at 900 C carries 99 kJ/m^2, and a flat 0.5 s pulse peaking at
800 C carries 127 kJ/m^2. The model uses room-temperature cp; the real cp rises with temperature,
so real faces would peak lower for the same dose.

**What this means for the model's constants [E]:**

- Thermochemically, a sub-second face can go about 450-700 K above the 400 C rating before it
  loses a micrometre. 800 K is not a decomposition threshold on these time scales. Taken at face
  value, the kinetics put "1 um lost" at about 1150 K (0.5 s flat) to 1300-1370 K (0.1 s).
- **The extrapolation is the weak point.** The fit used data at 550-610 C, and the table runs it
  to 1000 C. At 1000-1200 K, `(R T^2/E) ln 10` is about 75-110 K, so each decade of error in the
  rate moves every threshold by that much.
- **One fast-heating datum points the unsafe way.** Kuper, Brannon & Brannon (1993), *Threshold
  behavior in polyimide photoablation: single-shot rate measurements and surface-temperature
  modeling*, Appl. Phys. A 56, 43-50, https://doi.org/10.1007/BF00351902, put the modelled
  surface temperature at the first detectable excimer ablation of Kapton at about 850 C for
  248/308/351 nm. That is with ~tens-of-ns pulses (abstract via search snippet). [S]
- The ESA kinetics would give almost no conversion at 850 C in nanoseconds. So either
  photochemistry assists (likely at UV), or the modelled temperature is off, or the
  high-temperature rate exceeds the extrapolation. None of these is resolved here.
- A defensible engineering limit is therefore the 0.1 um column minus one decade: about 670-870 C
  (940-1140 K) for 0.1-0.5 s pulses. That is still 270-470 K above `T_RATING`. **[E]**
- Mechanics, not chemistry, may bind first. Kapton's second-order transition (taken as Tg) is at
  360-410 C. [S] DuPont H-38492, op. cit. Above it the face layer softens; it does not melt. A
  softened 0.1-10 um skin on a sheet that will be used as an impactor probably does not matter,
  but if flatness or the sheet's strength in its last 10 ms matters, Tg is the real limit, not
  decomposition.
- **Outgassing / blistering.** Water uptake is 1.3-2.8 % (DuPont, op. cit.; Upilex S 1.4 %). [S]
  A film fast-heated straight from the ground can blister. After a long vacuum coast it is dry.
  **[E]**

### 1.4 A caution about microsecond laser data

Kappes et al. (2014), *A study of photothermal laser ablation of various polymers on microsecond
time scales*, SpringerPlus 3, 489, https://doi.org/10.1186/2193-1801-3-489, found that a 15 us
pulse ablated "their polyimide" at about 570 K, close to its TGA onset. [S]

That polyimide was a low-molecular-weight (Mn 746) 6FDA/hydroxy-biphenyl oligomer with much
monomer. It is not Kapton, so the result does not transfer.

---

## 2. Higher-temperature film chemistries

| film | chemistry | T5 / onset | Tg | density | ~25-125 um film available? | source |
|---|---|---|---|---|---|---|
| Kapton HN | PMDA-ODA | T5 608 C, Tmax 631 C (N2, 10 C/min); residue 54.6 % | 360-410 C | 1.42 | yes, 7.5-125 um | [S] ESA SP-540; DuPont H-38492 |
| PMDA-ODA (lab) | PMDA-ODA | T5 574 C | - | - | - | [S] Wu et al. 2020, below |
| Upilex S | BPDA-PPD | T5 631 C, Tmax 662 C; residue 62.4 %; E 312 kJ/mol | no Tg reported by UBE; "used in processes just below 500 C" | 1.47 | yes: 12.5, 25, 50, 75, 125 um | [S] ESA SP-540; [UBE Upilex-S catalogue](https://ube.com/upilex/catalog/pdf/upilex_s_e.pdf) |
| Kapton HPP-ST | PMDA-ODA, treated surface | same balance as HN per DuPont; no separate TGA found | not stated | 1.42 | yes | [S] https://www.dupont.com/products/kapton-hpp-st.html |
| Apical AH/NP (Kaneka) | PMDA-ODA type | no decomposition data found; Kaneka quotes -269 to 400 C | Kaneka says none | ~1.42 | yes | [S] https://kaneka.com/kaneka-americas/products/apical |
| PBO (Zylon), lab-cast film | PBO | decomposition reached 670 C in N2 | none (rigid rod) | ~1.54 (fibre) | **no commercial film found**; Toyobo sells fibre (decomposition 650 C quoted) | [S] Fukumaru et al. 2012, Macromolecules, https://doi.org/10.1021/ma3006526 |
| PBI (lab films) | polybenzimidazole | T5 533-535 C (N2) | 281-319 C (TMA) in patent films; ~410-427 C for bulk Celazole in vendor tables | ~1.3 | Goodfellow lists PBI film; thin film is not a mainstream product | [S] US patents 11548983 / 11655339 (search snippets); https://www.goodfellow.com/usa/polybenzimidazole-film-group |
| POSS-PMDA-ODA (side-chain, 30 wt %) | PMDA-ODA + DABA-POSS | T5 **512 C** (lower than the 574 C of the neat film) | - | - | lab | [S] Wu, H. et al. (2020), Polymers 12, 2865, https://doi.org/10.3390/polym12122865 |
| phosphine-oxide / fluorinated (6FDA) PIs | various | generally **lower** onset than PMDA-ODA; 6FDA PIs in these papers have Tg 230-280 C | - | - | lab | [S] Butnaru et al. 2018, JAAP, https://doi.org/10.1016/j.jaap.2018.06.015 (title/abstract only); Ghaemy et al. 2013, https://doi.org/10.1080/15685551.2013.840473 |

Notes:

- **Upilex S is the only drop-in.** It is a real, space-qualified film in the right thicknesses
  (125S has an area factor of 5.4 m^2/kg). It is better than Kapton HN by 23-31 K in T5/Tmax, and
  its higher activation energy (312 against 250 kJ/mol) makes its rate climb faster above the fit
  range. The Arrhenius lines cross near 700 C. Above that, Upilex S decomposes *faster* than
  Kapton HN. Its sub-second thresholds in section 1.3 come out 5-60 K *lower* than Kapton's at the
  0.1-1 um levels, and up to ~115 K lower at 10 um.
  So on the sub-second scale it buys nothing chemically. Its through-thickness conductivity is 0.29 W/m/K against 0.12 for
  Kapton [S] (UBE catalogue). That lowers the face for a given flux, which helps more than its
  chemistry. **[E]**
- Fluorination and phosphorus raise flame retardance (char) and transparency, not the onset. They
  are the wrong direction here.
- PBO and PBI films are not commercial products at 25-125 um that I could find. PBO film reaches
  670 C onset, which is about 40-60 K over Kapton.
- Costs: no primary-source prices were found. None are quoted.

---

## 3. Fillers and additives

| additive | loading | measured shift | atmosphere | source |
|---|---|---|---|---|
| graphene oxide | 0.1 wt % | T5 to ~587 C, T10 to ~603 C (~+2.4 % / +1.7 %, about +14 K) | not stated | [S] Chinese journal *Insulating Materials* 2021, https://castjournals.cast.org.cn/joweb/jycl/CN/PDF/10.16790/j.cnki.1009-9239.im.2021.11.018 (search snippet) |
| graphene + clay | 12.5 vol % graphene + 0.45-9 vol % clay | decomposition temperature up to **+50 C**; char decomposition in air up to +80 C | N2 and air | [S] Akinyi & Iroh (2023), Polymers 15, 299, https://doi.org/10.3390/polym15020299 |
| nanographene | up to 50 wt % | ~43 % lower decomposition **rate**, higher char | - | [S] Akinyi et al. (2021), Minerals 11, 168, https://doi.org/10.3390/min11020168 |
| organoclay (alkylammonium) | 3-8 wt % | the modifier itself degrades during imidization (darkening) | air / N2 | [S] Delozier et al. (2002), NASA, NTRS 20020023603 |
| synthetic clay, carbon nanofibre | - | "increases the onset temperature of matrix degradation", mostly oxidative | air | [S] Miller, S. G. (2007), NASA, NTRS 20070031911 |
| silica (sol-gel) | 5 % | weight loss "shifted to higher temperature" (BTDA copolyimide) | - | [S] Karatas et al., Marmara (search snippet); Rehman, Schmidt & Ahmad (2006), J. Macromol. Sci. A, https://doi.org/10.1080/10601320600602613 (PMDA-ODA hybrids, no number in abstract) |
| POSS | 30 wt % | T5 **-62 K** | - | [S] Wu et al. (2020), op. cit. |
| boron nitride | up to 30 wt % | decomposition 460-536 C reported; no matched neat control in the abstracts | - | [S] search snippets only; no shift can be quoted |

Reading: fillers move inert-atmosphere TGA onset by roughly 0 to +50 K. The larger quoted gains
are oxidative (barrier) effects, which do not exist in vacuum. Heavy graphene loading cuts the
*rate* by about 2x. On the kinetics above, a factor of 2 in rate is worth about 30 K. **No filler
plausibly gives >100 K.** **[E]**

What a filler does change is through-thickness conductivity. BN or graphite at tens of percent
raises it, which lowers the face temperature for a given flux. That lever is physical, not
chemical, and would need its own conduction run. **[E]**

---

## 4. Surface treatments and overlayers (thin, bonded)

- **SiOx and POSS-Kapton (CORIN XLS POSS, main-chain POSS-Kapton).** These are atomic-oxygen
  defences. AO turns the surface into silica. [S] Tomczak et al. (2009), AIP Conf. Proc.,
  https://doi.org/10.1063/1.3076863: 3.5 wt % POSS-PI eroded ~2 um over 3.9 years on MISSE-1,
  against ~32 um in 4 months for unfilled PI. [S] Wu et al. (2020): 30 wt % POSS gives an erosion
  yield of 1.1e-25 cm^3/atom, against 3.0e-24 for PMDA-ODA.
  - Their benefit is chemical (oxidation). The plume's O/C/N/H at ~50 eV/amu is not 5 eV LEO AO,
    and the thermal limit is set by heat, so **these do not raise the thermal limit**. POSS
    copolymers *lower* T5 (section 2).
- **Bonded thin coats (SiOx, Al2O3, metal, 1 um graphite).** A 50 nm-1 um coat conducts across
  itself in microseconds. It sits at the face temperature, stores almost nothing, and cannot
  reflect heat that arrives as particles. This repo already checked a 1 um graphite coat: it
  lowers the face 3-12 K. An Al coat warms the face 4-39 K. **[S, this repo]**
  `docs/fast_growth_wave_answers_for_parent.md`. No thin bonded coat moves the thermal limit.
- **Where a thin inorganic coat could matter is erosion, not heat.** See section 5.2. Inorganic
  overlayers sputter far less than polymers. A 0.1-1 um SiOx or carbon coat could stop the
  sub-micron particle erosion of the face. **[E]**

---

## 5. Char, and particle bombardment

### 5.1 Does a char layer protect the film in vacuum?

Kapton chars rather than melts: 54.6 % residue in inert gas [S] (ESA). Char made from Kapton HN
between 500 and 1200 C has a thermal conductivity of 0.2-2 W/m/K, rising with carbonisation
temperature. [S] Kondratenko et al. (2022), *Nanoscale thermal conductivity of Kapton-derived
carbonaceous materials*, J. Appl. Phys., https://doi.org/10.1063/5.0074407. That is 2-15 times
Kapton's 0.12. Laser-converted Kapton (LIG) survives >900 C in Ar, against ~550 C for the PI.
[S] Lin et al. (2014), Nat. Commun. 5, 5714, https://doi.org/10.1038/ncomms6714.

So a char skin is not an insulator. It conducts *better* than the film beneath it, and so spreads
heat inward. It does survive to much higher temperatures, and a black char face (emissivity likely
~0.9 against Kapton's ~0.8 in the model) radiates a little more.

At the 1-10 um depths reached in section 1.3, the char adds almost no heat capacity and only
slightly raises emission. The pyrolysis endotherm takes some heat, but no enthalpy was found in a
primary source. **The char is self-limiting, not protective.** Once formed, it does not keep the
Kapton beneath any cooler. **[E]**

Ablation-style cooling (blowing) needs mass loss of order the absorbed dose divided by the heat of
ablation. For 100 kJ/m^2 and a few MJ/kg that is ~0.03 kg/m^2, about 20 um of Kapton. **[E]** So
ablation is only "protective" at the cost of the 10 um+ losses the sheet is trying to avoid.

### 5.2 Particle bombardment at ~50 eV/amu

At 50 eV/amu the projectiles carry H 50 eV, C 600 eV, N 700 eV and O 800 eV. **[E]**

- **Range.** SRIM-class ranges for these in a polymer are a few nm, not 0.1 um. **[E]** Over
  0.1-0.5 s the thermal penetration depth is `sqrt(alpha t)` ≈ 90-200 um, so the deposition depth
  is irrelevant to conduction.
- **Measured / modelled erosion.**
  - Di Sarcina & Scaglione (2009), *The thermal spike model of polyimide erosion by low energy ion
    bombardment*, J. Optoelectron. Adv. Mater. 11, 1157-1162, bombarded Kapton with End-Hall Ar
    and O beams. Their SRIM-plus-thermal-spike model, which matched the measured erosion, gives up
    to ~2 atoms/ion for Ar and up to ~35 atoms/ion for O, at <=150 eV. [S]
    https://joam.inoe.ro/articles/the-thermal-spike-model-of-polyimide-erosion-by-low-energy-ion-bombardment/fulltext
  - Below 500 eV the PI "maintains its integrity" under Ar+. [S] (search snippet; secondary review).
  - LEO AO (5 eV) erodes Kapton at 3.0e-24 cm^3/atom. [S] Wu et al. 2020; NASA AO literature.
- **Arithmetic [E].** The dose is 40-100 kJ/m^2 at 50 eV/amu (4.8e9 J/kg), which is 8-21 mg/m^2.
  If mostly heavy atoms (~14 amu), that is about 3.6-8.9e16 atoms/cm^2. Kapton's atomic volume is
  1.15e-23 cm^3. At 10-50 atoms removed per heavy projectile the erosion is **~0.04-0.5 um per
  pass**. The AO-style chemical yield alone would give only a few nm.
- **This is a floor that no temperature limit removes.** It is comparable to the "lose 1 um"
  criterion. Inorganic layers sputter at ~0.1-1 atom/ion (order-of-magnitude, **[E]**; no primary
  source checked here). So a 0.1-1 um SiOx or carbon coat is the one overlayer with a real job:
  against erosion, not heat.

---

## 6. Porous low-density carbon overlayers (10-50 um) as an insulating radiator

The idea: a porous carbon skin bonded on the plume face runs hot (1500-2500 K) and radiates. Its
low through-thickness conductivity holds the Kapton beneath cool.

### 6.1 What it would have to do [E]

The layer's diffusion time `L^2/alpha` is milliseconds for a 10-50 um aerogel. Against a 0.1-0.5 s
pulse it is therefore quasi-steady. The flux splits into what the outer face radiates,
`eps sigma Ts^4`, and what it conducts, `(Ts - Ti) k / L`.

To radiate everything, a 0.9-emissivity face must sit at 1407 K for 200 kW/m^2, or 2104 K for
1 MW/m^2. The share conducted into the Kapton (face Ti held at 700 K or 1100 K):

| q | layer | k = 0.003 W/m/K | k = 0.03 | k = 0.3 |
|---|---|---|---|---|
| 200 kW/m^2, Ti 700 K | 20 um | 40 % (Ts 1237 K) | 89 % | 93 % |
| | 50 um | 19 % | 78 % | 93 % |
| 1 MW/m^2, Ti 700 K | 20 um | 19 % (Ts 1994 K) | 87 % | 98 % |
| | 50 um | 8 % | 59 % | 97 % |
| 1 MW/m^2, Ti 1100 K | 20 um | 14 % | 70 % | 91 % |
| | 50 um | 6 % | 44 % | 88 % |

**The requirement is an effective conductivity of a few mW/m/K at 1500-2000 K, in a 10-50 um
layer.** Radiation inside a porous layer adds `k_r = 16 sigma T^3 / (3 beta)`:

| extinction coefficient beta | optical depth over 20 um | k_r at 1000 K | at 1500 K | at 2000 K |
|---|---|---|---|---|
| 1e4 /m | 0.2 | 0.030 | 0.10 | 0.24 |
| 1e5 /m | 2 | 0.003 | 0.010 | 0.024 |
| 1e6 /m | 20 | 0.0003 | 0.001 | 0.0024 |

- Only a layer with beta of order 1e6 /m (optically thick within ~1 um) and solid conduction of a
  few mW/m/K helps meaningfully. A low-density (~0.1 g/cm^3) carbon at a specific extinction of
  ~100-1000 m^2/kg gives beta ~1e4-1e5 /m (estimate, not measured here).
- A layer that thin is also nearly transparent in the IR, so the Kapton sees the hot layer's
  radiation directly.
- **Measured high-temperature carbon aerogels are 100x too conductive.** Wiener's thesis
  extrapolates a lowest ~0.17 W/m/K at 2500 C. A Netzsch poster reports 0.23-0.48 W/m/K at 1500 C
  (in Ar). A C-fibre/carbon-aerogel composite gives 0.69 W/m/K at 1800 C. [S] (below)
- From these numbers, **a bonded 10-50 um porous carbon skin still passes about 60-95 % of the
  flux into the Kapton.** A spaced, free-standing shield does better, because its "conductivity"
  is only the radiative exchange across a gap. That is `two_film()` in `coated_sheet.py`. **[E]**

### 6.2 What exists

| material | density | areal mass at ~20 um | through-thickness k (vacuum) | emissivity | flexibility / foldability | adhesion to PI | can it go on Kapton below ~400 C? | maturity |
|---|---|---|---|---|---|---|---|---|
| **Polyimide aerogel thin film** (NASA Glenn; Meador et al.) | 0.14-0.3 g/cm^3 | 3-6 g/m^2 [E] | as low as 14 mW/m/K at RT (in air; lower in vacuum); no high-T data | not a carbon; IR-semi-transparent | rolls and folds back on itself without cracking (films cast roll-to-roll; ~0.5 mm "nominal") | same chemistry; could be cast on Kapton | **yes** (cast, chemically imidised, supercritically dried) | TRL ~4-5; licensed commercially (Blueshift AeroZero) |
| **Carbon aerogel** (RF- or PI-derived) | 0.1-0.7 g/cm^3 | 2-14 g/m^2 [E] | 0.17-0.7 W/m/K at 1500-2500 C (in Ar, extrapolated); low-T vacuum values much lower | high (black), not measured in these sources | brittle monolith; thin flexible films not found as a product | must be bonded | **no**: needs 800-2500 C pyrolysis, so it must be made separately and bonded; adhesives limit service to <~400 C | lab / bulk-insulation |
| **PI aerogel carbonised in place** | - | - | - | - | - | - | **no**: carbonising needs >600-800 C, which also destroys the Kapton | concept |
| **Laser-induced graphene (LIG)**: porous carbon made *from the Kapton face itself* | porous foam (density not given) | from the Kapton's own mass | not measured in source (char made in bulk: 0.2-2 W/m/K) | black | as flexible as the film it is on | integral (it *is* the film) | **yes**: CO2 laser in air; local >2500 C, bulk stays cool; thickness up to ~20-35 um vs laser power | lab, widely reproduced |
| **Vertically aligned CNT forest** (super-growth, Vantablack-type) | ~0.03-0.08 g/cm^3 (8 vol % forests) | ~0.6-1.6 g/m^2 [E] | **high along the tubes**: 15 W/m/K (8 vol %), 37-127 W/m/K reported; ~4 W/m/K across | 0.98-0.99 (0.2-200 um) | forests crush, cannot be touched or folded | must be transferred; weak | **grown at ~700-800 C**: no. Vantablack S-VIS spray version cures at 100-150 C | grown: flight-proven on metal; spray: commercial |
| **Drawn CNT sheet / Miralon-type** | ~1.5 mg/cm^3 aerogel sheet (single drawn sheet), densifies to ~0.5 g/cm^3 | single sheets ~µg/cm^2; stacks per need | in-plane 50 W/m/K; normal to plane ~1 W/m/K; stacked with spacers 0.01-0.5 W/m/K | high | very flexible | lamination or adhesive | **yes** (room-temperature lamination) | commercial sheets |
| **Acktar black (Magic/Fractal Black)** | inorganic (metal-based), not carbon | 3-14 um coat | not given | >0.97 (3-10 um), >0.88 (3-30 um) for Fractal Black | on 25-75 um Kapton (MaxiBlack) | qualified on PI | **yes**: vacuum deposition, sold on Kapton; rated -269 to +350/450 C | flight-qualified |

Sources for section 6:

- Meador, M. A. B. et al. (2012), *Mechanically strong, flexible polyimide aerogels cross-linked
  with aromatic triamine*, ACS Appl. Mater. Interfaces 4, 536.
  https://doi.org/10.1021/am2014635. Densities from 0.14 g/cm^3; films fold back on themselves
  without cracking; Tg 270-340 C; decomposition onset 460-610 C. [S]
- NASA Glenn PI aerogel technology summary: 0.1-0.3 g/cm^3, as low as 14 mW/m/K at RT.
  https://ntrs.nasa.gov/citations/20120011884 ; flexibility only in thin (~0.5 mm) films:
  https://ntrs.nasa.gov/api/citations/20210011477/downloads/PIU%20Aerogels%20-%204-13-2021%20STRIVES.pdf
  [S] (search snippets)
- Wiener, M. (2009), dissertation, Univ. Wurzburg, carbon aerogels for >1000 C insulation; lowest
  ~0.17 W/m/K extrapolated at 2500 C.
  https://opus.bibliothek.uni-wuerzburg.de/files/3569/Disse_Wiener_2009_Druckversion_final.pdf
  [S] (snippet; the PDF would not download here)
- Zhang, X., Netzsch poster: carbon aerogel 0.23 (with GO) / 0.48 W/m/K at 1500 C in Ar.
  https://hausmesse-at.netzsch.com/_Resources/Persistent/f/4/e/5/f4e5bf793f38820e9b6c5c41cb33d60e57d98c15/Xiaomeng%20Zhang.pdf
  [S]
- C-fibre-reinforced carbon aerogel, 0.69 W/m/K at 1800 C, 0.68 g/cm^3. RSC Adv. 2022,
  https://pubs.rsc.org/en/content/articlelanding/2022/ra/d2ra00276k [S]
- Lin, J. et al. (2014), *Laser-induced porous graphene films from commercial polymers*, Nat.
  Commun. 5, 5714, https://doi.org/10.1038/ncomms6714. LIG made on Kapton by CO2 laser in air;
  thickness rises with laser power (figure 3b, up to ~20-35 um); LIG stable >900 C in Ar against
  ~550 C for PI. [S]
- Mizuno, K. et al. (2009), *A black body absorber from vertically aligned single-walled carbon
  nanotubes*, PNAS 106, 6044, https://doi.org/10.1073/pnas.0900155106. Near-perfect absorption
  0.2-200 um; 0.98-0.99 from spectral emissivity (follow-up article by the same group). [S]
- Hata, K. et al. (2004), *Water-assisted highly efficient synthesis of impurity-free
  single-walled carbon nanotubes*, Science 306, 1362, https://doi.org/10.1126/science.1104962.
  CVD super-growth forests to 2.5 mm. [S] Growth temperature (~750 C) is from general knowledge of
  that process. **[E]**
- Ivanov, I. et al. (2006), VANTA longitudinal conductivity 15.3 W/m/K at 8 vol %, via search
  snippet. https://impact.ornl.gov/en/publications/fast-and-highly-anisotropic-thermal-transport-through-vertically-/
  [S] Sci. Rep. srep21014: 127 W/m/K along the tubes, ~4 W/m/K out of plane.
  https://www.nature.com/articles/srep21014 [S]
- Surrey NanoSystems, Vantablack S-VIS: the substrate must withstand 100-150 C.
  https://www.coatingsworld.com/exclusives/the-worlds-blackest-material-is-now-in-spray-form [S]
- CNT sheet conductivity: in-plane ~100 W/m/K and normal ~1 W/m/K; stacked insulators 0.01-0.5
  W/m/K. US patent 8722171, *Nanotube-based insulators* (Nanocomp).
  https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/8722171 [S] Drawn MWNT sheets,
  50 W/m/K along the draw:
  https://somwritinglab.utdallas.edu/research/nanotech/publications/documents/ThermalconductivityofcarbonMWNTsheets-blackbodyradiationandquenchingofphononmodesinnanotubeb.pdf
  [S] Drawn-sheet density ~1.5 mg/cm^3 is from Zhang, Atkinson & Baughman (2005), Science 309,
  1215 (from memory, not re-checked). **[E]**
- Acktar: https://acktar.com/coated-black-foil/ (MaxiBlack on 25/50/75 um Kapton);
  https://acktar.com/?p=433 (Fractal Black 5-14 um, emissivity); [satsearch Fractal Black on PI
  film](https://satsearch.co/products/acktar-fractal-black-diffusive-polyimide-film). [S]

### 6.3 Reading [E]

- **The vertically aligned forest is the wrong geometry.** The tubes conduct along their length,
  normal to the film, at 4-127 W/m/K. That is a heat pipe into the Kapton, not an insulator. Its
  near-perfect emissivity is real, but a 20 um forest at, say, 15 W/m/K has `L/k` ~1e-6 m^2K/W, so
  the Kapton face sits within a kelvin of the forest top. It also cannot be grown on Kapton.
- **Polyimide aerogel** is the only low-k porous layer that can be made on Kapton below 400 C, and
  it folds. It is a polyimide, though, with onset 460-610 C, so it decomposes long before reaching
  1500 K. It cannot be the radiator. It could be a spacer *beneath* a separate radiating skin, but
  at 20 um its R is only ~1e-3 m^2K/W even at 14 mW/m/K.
- **Carbon aerogels** survive the temperature but need furnace pyrolysis. They are brittle and,
  at 1500-2500 C, about 100x too conductive for 10-50 um to matter.
- **Stacked drawn CNT sheets** are flexible, laminate cold, and reach 0.01-0.5 W/m/K normal to the
  plane in stacks. The low end (0.01) at 50 um would conduct ~40-60 % of the flux (table 6.1).
  This is the best candidate on paper, but there is no high-temperature through-plane data and the
  sheets are optically thin individually.
- **LIG** is the cheapest way to get a black porous carbon face, because it is made out of the
  sheet itself. It is a char, though, and conducts like one (section 5.1).
- **Acktar-type coats** are thin, flight-qualified and very emissive, but they are 3-14 um dense
  inorganic films, not insulators. They act like the thin coats in section 4.

The overall verdict: no bonded 10-50 um porous carbon layer from the literature reaches the
few-mW/m/K effective conductivity at 1500-2000 K needed to keep most of a 0.2-1 MW/m^2 pulse out
of the Kapton. The physics favours a *gap* (a free-standing shield) over a bonded porous skin.

### 6.4 The thinnest commercial carbon textiles

These are candidates for a free-standing carbon shield or spacer, or a bonded skin.

| product | form | lightest areal weight | thickness | binder / matrix and its limit | tear / handling | source |
|---|---|---|---|---|---|---|
| NTPT ThinPreg (North Thin Ply Technology) | spread-tow UD carbon **prepreg** tape | 15 g/m^2 fibre | not stated in the source; ~0.015-0.02 mm at that weight **[E]** | high-Tg epoxy (ThinPreg 120EPHTg-402); organic, so it fails long before carbon | UD tape: splits easily along the fibre | [S] https://www.compositesworld.com/products/unidirectional-carbon-fiber-prepreg-tapes-just-15-gsm ; https://www.unmannedsystemstechnology.com/2015/07/ntpt-introduces-ultra-thin-15gsm-prepreg-tape-for-uavs/ |
| TeXtreme (Oxeon) spread-tow woven | dry woven fabric | 60 g/m^2 (HR40 HM, reseller); 76 g/m^2 IM ±45° | ~0.12 mm as a 35 % Vf laminate (R&G calculation) | small binder, ~10 g/m^2 added (R&G) | woven, so it tears less readily than tape | [S] https://compositeenvisions.com/?p=1332 ; https://shop1.r-g.de/en/art/190276-45 . No Oxeon datasheet was reached; lighter custom grades may exist |
| Chomarat C-Ply, Hexcel | thin-ply NCF / fabric | not found | - | - | - | no datasheet found |
| conventional woven 1k | plain weave | 90 g/m^2 | ~0.1 mm per layer | none (dry) | - | [S] https://easycomposites.co.uk/90g-plain-weave-1k-carbon-fibre-cloth |
| Optiveil 20332A (TFP / James Cropper) | wet-laid nonwoven carbon veil, 12 mm × 7 um fibres | **2 g/m^2** | **0.03 mm** at 10 kPa | cross-linked styrene acrylic. The TDS gives no temperature limit. Acrylics depolymerise at roughly 300-400 C **[E]**, so above that the veil is loose fibres with no binder | very weak: tensile 4 N/15 mm along the machine direction, 0.8 N/15 mm across it | [S] https://advancedmaterials.jamescropper.com/wp-content/uploads/sites/2/2026/06/JC-OPTIVEIL-20332A-Carbon-Veil-TDS.pdf |
| Hollingsworth & Vose carbon veils | nonwoven | not found | - | - | - | no datasheet reached |
| Miralon sheet (Nanocomp / Huntsman) | CNT nonwoven sheet | 12 g/m^2 reference grade; 6-30 g/m^2 custom (8-12, 10-15, 15-20 standard) | 25 um (Huntsman TDS) or 20 um (reseller); density 0.65 g/cm^3 | none (pure CNT) | flexible, tough | [S] Huntsman MIRALON Sheet TDS, https://huntsman.com/docs/Documents/MIRALON%20Sheet%20Technical%20Datasheet.pdf (via search snippet; direct download 404'd); https://www.perigeedirect.com/products/miralon-sheets-carbon-graphene-nanotubes-s-t01avb-12-100812-100816 |
| Lintec drawn CNT sheet | drawn aligned CNT sheet | not found in a datasheet | - | none | - | no source reached |

Reading **[E]**:

- The lightest carbon you can buy as a sheet is the 2 g/m^2, 30 um Optiveil veil. Its binder burns
  off at a few hundred C, and across the machine direction it has almost no strength.
- Miralon at 6-12 g/m^2 and 20-25 um is the lightest *strong* all-carbon sheet. It has no binder,
  and its 0.65 g/cm^3 density is far from aerogel-like.
- Woven spread-tow fabric starts at about 60 g/m^2. That is about 40 % of the 100 um Kapton
  sheet's own 142 g/m^2, a large mass toll for a shield.

---

## 7. Ranked summary

**Plausibly moves the face limit by >100 K for a sub-second pulse:**

1. **Recognising the limit is kinetic, not the 400 C rating** (section 1). By the ESA inert-gas
   kinetics, a 0.1-0.5 s face loses <0.1 um below about 760-970 C, so the limit is about
   1040-1245 K. Even allowing one decade of rate error it is about 940-1140 K, against the model's
   673 K. This is +300 K or more, with no material change. **[S kinetics, E extrapolation]**
2. **A spaced sacrificial shield** (already `two_film()`): the only geometry that keeps most of the
   heat off, because its "conductivity" is radiative exchange across a gap. **[E]**

**Moves it by tens of K, not >100 K:**

3. Upilex S (BPDA-PPD): +23-31 K in slow TGA T5/Tmax, but 5-60 K *worse* on the sub-second
   thresholds because its kinetics overtake Kapton's above ~700 C. Its 2.4x higher conductivity
   lowers the face for a given flux, which is the only real gain. That gain is not computed here.
   **[S/E]**
4. Fillers (graphene, clay, silica): 0 to +50 K in inert TGA, up to 2x lower rate. **[S]**
5. PBO film: about +40-60 K onset, but not commercial as film. **[S]**

**Does not help, or hurts:**

6. Thin bonded coats (SiOx, Al, 1 um graphite, Acktar black): -39 to +12 K. They neither block
   nor store the heat. **[S, this repo]**
7. POSS-PI (-62 K T5), fluorinated and phosphorus PIs (lower onset). **[S]**
8. Char: conducts better than the film, so it is self-limiting rather than protective. **[S/E]**
9. Bonded porous carbon (CNT forests, carbon aerogel, LIG): passes 60-95 % of the flux at the
   conductivities measured. CNT forests conduct *along* the tubes and are the worst choice.
   **[S/E]**

**Unknown:**

- Kapton's decomposition rate at 800-1100 C under ms-s heating; no direct measurement was found.
  The UV-ablation datum (~850 C at ns) disagrees with the extrapolation.
- Erosion of Kapton by 0.6-0.8 keV C/N/O. The estimate is 0.04-0.5 um per pass; its yield could
  be several times higher.
- The face's mechanical state above Tg in the impactor's last milliseconds.
- High-temperature through-plane conductivity of stacked CNT sheets.

---

## Calculation scripts (session scratchpad, not committed)

- `kin.py`: the Arrhenius sets and the isothermal table (section 1.2).
- `depth.py`, `thresh.py`: 1-D conduction plus conversion, and the 0.1 / 1 / 10 um face thresholds
  (section 1.3).
- `overlayer.py`: the quasi-steady overlayer split and radiative conductivity (section 6.1).

They were run under the project's `.venv` Python (numpy, scipy). Per the repo rule that numbers in
docs/ADRs come from committed modules, port these to `python/puffsat/xray_collision/` with a Make
target before any of these figures enter an ADR or a handoff.
