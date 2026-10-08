# Hydrogen and the wall steels: maraging on the plate, Cr-Mo in the chambers

Status: research note, 2026-10-07. The question: is maraging steel damaged by hydrogen? Does an
ablative coat protect it? Is there a hydrogen-tolerant steel nearly as strong? It covers the spray
plate ([`spray_plate_fatigue_and_spall.md`](spray_plate_fatigue_and_spall.md); maraging 300) and
the walled chambers ([`walled_nozzle_chamber_blast.md`](walled_nozzle_chamber_blast.md); Cr-Mo).

**Hydrogen reaches both walls.**
- **Plate:** the water PuffSat dissociates at 14-60 kK, so atomic H and H2 strike the face every
  pulse. The argon arm's PuffSat is water too.
- **Hydrogen chamber:** it is hydrogen throughout, during and between pulses.
- **Methane chamber:** 58% of the hydrogen is dissociated at 7000 K (`chambers.csv`), so atomic H
  reaches its wall during every pulse.

## 1. Maraging steel is among the worst steels in hydrogen

NASA's compilation [L16], Table 3 gives the ratio of notched tensile strength in hydrogen to that
in helium. It was measured at 24 °C and 69 MPa H2, and reads below 0.49 as "extreme". The rows
that matter here:

| steel | notched strength ratio in H2 | rating | yield (MPa) |
|---|---|---|---|
| **18Ni-250 maraging** | **0.12** | extreme | 1709 |
| AerMet 100 (peak aged) | 0.15 | extreme | |
| 17-4 PH | 0.18 | extreme | 1075 |
| 4140 (Q&T) | 0.25-0.40 | extreme | 1233 |
| 4340 | 0.35 | extreme | 1302 |

Maraging keeps about an eighth of its notched strength. Susceptibility rises with strength, and
maraging is about as strong as steel gets. Gaseous-hydrogen crack growth in 18Ni maraging was
measured from -6 to +100 °C [GW]. Growth fades above a critical temperature inside that range,
and the critical temperature rises with hydrogen pressure.

## 2. The hydrogen-tolerant alternatives are about half as strong

Rated negligible or small [L16], Tables 3-4:

| alloy | notched strength ratio | rating | yield (MPa) | note |
|---|---|---|---|---|
| **NASA-HR1** (Fe-Ni-Co superalloy) | 0.98 (34.5 MPa H2) | negligible | **944** (UTS 1323) | developed at MSFC for hydrogen engines |
| Inconel 718, solution treated at 1900 °F | 0.92 | small | 1075 | the common 1750 °F treatment is "extreme" (0.46-0.53) |
| Waspaloy (powder metallurgy) | 0.95 | small | 1199 | |
| MERL 76 (powder metallurgy) | 0.96 | small | 1034 | |
| A-286, solution treated | 0.97 | negligible | 847 | **aged A-286 is "severe" (0.51)** |
| Nitronic 50 (22-13-5) | 1.00 | negligible | 586 | |
| 316 | 1.00 | negligible | 441 | |
| Al alloys (2024, 6061, 7075...), Cu, GRCop-84 | ~1.0 | negligible | 120-450 | |

None approaches maraging's 1.7-2.0 GPa yield. The best are ~0.95-1.2 GPa, at roughly
half the Hugoniot elastic limit (`Y (1-nu)/(1-2nu)`): ~1.6-2.1 GPa against 3.1-3.6. Their spall
strengths are not sourced here. Low-alloy steel held under ~950 MPa tensile strength (the
parent's Cr-Mo) is not on the negligible list. It is qualified for hydrogen by fracture mechanics
on the ASME CC2938 crack-growth curves, which is what `chamber_fatigue.py` already does.

## 3. Is the ablative coat a barrier? Not one we can count on

- **Pitch is not designed as a barrier.** It is a hydrocarbon that chars and gives off H2 and
  light hydrocarbons where it is hot. Its cool layer next to the steel is the part that stays
  intact. No hydrogen-permeation data for coal-tar pitch or its char were found.
- **Hydrogen reaches the steel quickly at these temperatures.** The plate steel runs at 400-507 K
  and the chamber steel at 350-420 K, where hydrogen diffuses quickly in ferritic and martensitic
  steels. Per pulse it enters only the near-surface (tens to hundreds of µm over a push), but
  that is where cracks start.
- **Temperature helps maraging, but cannot be relied on.** The plate face sits above the tested
  range in which maraging's gaseous crack growth fades. But the pulse's hydrogen is atomic and at
  GPa-level pressure, the critical temperature rises with pressure, and the first pulses of every
  push meet a cold plate.

## 4. Recommendations

**Spray plate (maraging 300, ADR-0055):**
1. **Keep maraging 300 for the face, but put a designed hydrogen barrier between it and the
   pitch.** Aluminium is rated negligible [L16] and its oxide is a known strong permeation
   barrier. A few tens of µm of sprayed aluminium, or an aluminide or alumina coating, runs well
   below aluminium's limits at the plate's 400-507 K. Its permeation reduction under pulsed plasma
   exposure needs a test.
2. **Or the hydrogen-tolerant fallback face: NASA-HR1 or 1900 °F-treated Inconel 718**, yield
   ~0.95-1.1 GPa. Its face allowable falls to roughly 1.2-1.5 GPa against maraging's 2.5. The
   merged pulse (0.89 GPa, `spray_plate_fatigue_and_spall.md`) still fits. The failed-merge pulse
   (2.5-4.1 GPa) does not, so fault tolerance gets worse.
3. Either way, the paper's P12 should say maraging is "extreme" in hydrogen, and that the face
   needs a barrier or a different alloy.

**Walled chambers. Superseded 2026-10-08 by `walled_nozzle_chamber_blast.md` §1 and §3d.** The
design is now a 10 mm Cr-Mo shell under a dry Kevlar wrap on a bulged chamber.
- Its thin shell sees ≤0.03x spall, so hydrogen barely matters to it. The alumina film stays as
  defence in depth.
- Solid Cr-Mo at the plug band cracks on the hydrogen curve within 30-660 pulses.
- A **maraging 300 band behind the alumina barrier** is the fallback. It is scored on the air
  curve, so it holds only if the barrier works.
- The steel under the pitch also meets H2 from the pitch's own pyrolysis, and it runs at
  ~350-420 K, where embrittlement is strongest.

The original assessment, kept as the record:

**Walled chambers (Cr-Mo, as the paper has it; maraging was never recommended there):**
1. **Keep the parent's Q&T Cr-Mo held under 950 MPa UTS,** qualified on CC2938. The multilayer
   steel wall of the chamber note is this steel.
   - Its spall strength, ~2.9 GPa (38KhN3MFA, a Q&T Cr-Ni-Mo-V steel [M17]), covers the long-plug
     spike.
   - Maraging would add little: the chamber's binding limit is spall and resin, not the face's
     compression limit, and maraging would bring the "extreme" hydrogen penalty.
2. **Hydrogen chamber:** the GRCop-84 liner (negligible [L16]) carries most of the hydrogen
   exposure. Its effectiveness as a permeation barrier for the Cr-Mo behind it is unquantified.
3. **Methane chamber:** pitch alone faces atomic H during each pulse. Add the same aluminium or
   alumina barrier on the steel under the pitch, or run the steel on the hydrogen crack-growth
   curve (the bound already given in `walled_nozzle_chamber_fatigue.md`).
   - The multilayer wall's inner shell can be a hydrogen-tolerant alloy as a sacrificial liner.
     A-286 (solution treated) or Nitronic 50 are negligible-rated with 0.6-0.85 GPa yield, under
     the Cr-Mo shells.

## 4b. The barrier: a thin, dense ceramic film (estimates, not yet sourced)

**Recommended:** a sub-micron dense **alumina** film on the steel, under the ablative.
- **Why thin works.** A film cracks only if the elastic energy it stores exceeds its fracture
  energy, and that energy scales with thickness. The channel-cracking estimate
  `h_c ~ Gamma E_f / (2 sigma^2)` (alumina E ~380 GPa, assumed Gamma ~30 J/m^2) gives:
  - crack-free to ~70 µm at the face's 0.06% thermal-mismatch strain;
  - ~3 µm at 0.3%;
  - ~1 µm at 0.5%.
- **Why the plate face suits it.** The face is compressed through its thickness without in-plane
  stretch, and it is always in compression. Reflected tension appears at the plate's back, not
  its face. The film is far thinner than the shock and simply rides with the steel.
- **The real limit is pinholes, not cracking.** Answers: atomic layer deposition (pinhole-free at
  tens of nm), ceramic/metal nanolaminates, or an aluminized-and-oxidized surface that regrows
  its oxide.
- **Process fit to the steel.**
  - **Maraging:** it must stay below its ~480 °C aging temperature after heat treatment, so
    sputtered (PVD) or ALD alumina, applied after the final age. A 20 m face would be coated in
    segments.
  - **Cr-Mo chambers:** tempered hotter, so pack or slurry aluminizing then oxidation (~900-1100
    °C) can be sequenced before the final temper.
- **Cost:** order 10^2-10^3 $/m^2. On ~315 m^2 of plate face that is ~$0.03-0.3M, a few percent of
  the 150 t of maraging alone.
- **Alumina film, not an aluminium layer, on the plate face and chamber walls.** Armor-grade
  aluminium (5083/7039) has a shock-compression limit of roughly 0.4-0.6 GPa (estimate), below even
  the merged 0.9 GPa pulse. It would yield every pulse, ratchet thinner, and re-crack its own oxide
  skin, which is where its hydrogen resistance comes from. 7xxx alloys also soften above
  ~120-150 °C. **Aluminium cladding or aluminized steel suits the low-pressure parts:** the skirt
  (5.5-9.5 MPa) and hydrogen-wetted lines.
- **Rejected:**
  - *Tungsten:* brittle when cold, and it blisters under intense low-energy hydrogen plasma. Its
    fusion role is erosion and heat, not permeation.
  - *CMCs:* microcracked and porous by design, so they would need a dense seal coat anyway.
  - *Silicone:* among the most gas-permeable solids, used industrially as gas membranes. It can
    be the ablator, not the barrier.
  - *Copper or aluminium layers on the plate:* they yield every 2.5 GPa pulse.
- **Risks:**
  - Where the pitch burns through, the film meets ~3,900 K plasma and alumina melts at ~2,070
    °C, so it depends on the respray.
  - There is no data for the permeation reduction under thousands of shock pulses with atomic
    hydrogen plasma. **That is a test, and the paper should list it as one.**

## 5. Sources

- [L16] J. A. Lee, "Hydrogen Embrittlement", NASA/TM-2016-218602 (2016), Tables 1-4.
  https://ntrs.nasa.gov/citations/20160005654
- [GW] R. P. Gangloff and R. P. Wei, gaseous hydrogen crack growth in 18Ni maraging steels,
  Metall. Trans. A (1977); temperature range and critical-temperature behaviour as summarized at
  https://www.science.gov/topicpages/c/crack+growth+kinetics.html (critical value not retrieved).
- [M17] Mescheryakov et al., Mater. Phys. Mech. 32, 152 (2017): spall strength of Q&T and
  maraging steels.
- Parent `templateArxiv.tex`, `sec:steel_chamber_service`, citing San Marchi et al. 2017 and
  Lee 2016.
