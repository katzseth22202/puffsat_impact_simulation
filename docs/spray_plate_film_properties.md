# Spray-plate film: carbon-loaded oil against pitch

Research note for the spray plate's sacrificial film (ADR-0055, `docs/argon_plate_handoff.md`).
It asks whether a carbon-loaded oil (~20% carbon by mass, the child paper's film and Orion's
"antiablation oil") should replace the paper's pitch in `walled_nozzle/wall_layers.py`.
That model takes five numbers: k, ρ, c, a surface cap temperature, and a removal heat booked
per kilogram removed at the cap (sensible heat up to the cap is carried by c, not by it).

Bottom line. In the repo's 1-D model a lower cap runs the steel much cooler, as expected.
But the removal heat is 20-60x below pitch's per unit volume, so the same radiant energy eats
0.7-5 mm of oil per argon pulse. That is 5-35x the 150 µm film. Oil wins only if vapor and
soot shielding cut the energy reaching the surface by about that factor. Orion's own design
relied on exactly that, but the size of the cut for this plate is not yet known.

## Sources (tags)

- [GA-5009-III] General Atomic, *Nuclear Pulse Space Vehicle Study, Vol. III: Conceptual
  Vehicle Designs and Operational Systems*, GA-5009, NASA contract NAS8-11053, 19 Sep 1964.
  Open copy: projectrho.com/public_html/rocket/supplement/GA-5009vIII.pdf (page image checked).
- [GA-6307] E. A. Day and J. C. Nance, *Nuclear Pulse Propulsion (Orion) Technical Status
  Summary and Ground-Oriented Development Plan*, GA-6307 (AD0361710), c. 1965. Open copy:
  documents.theblackvault.com/documents/space/AD0361710.pdf.
- [AIAA-2000] Schmidt, Bonometti, Morton, *Nuclear Pulse Propulsion: Orion and Beyond*,
  AIAA 2000-3856. NTRS 20000096503.
- [KF96] Shin-Etsu, *KF-96 Silicone Fluid* technical catalogue (kf96_e.pdf), property table
  and sections 10-12.
- [ROYCO] Royco 602 (MIL-PRF-87252 PAO coolant) product data sheet, lubritec.com and
  bigcommerce TDS SGP55352. Vendor data, ASTM D2766 for c.
- [XU2011] J. Xu, B. Yang, B. Hammouda, *Thermal conductivity and viscosity of self-assembled
  alcohol/polyalphaolefin nanoemulsion fluids*, Nanoscale Res. Lett. 2011 (NIST pub 908311).
- [THERMTEST] Thermtest application note, carbon black (LITX CB) in ethylene glycol, transient
  hot wire, 20-30 °C. Vendor measurement, not peer reviewed.
- [WALTERS] R. N. Walters, PhD thesis (UCLan), Tables 6-7, which reproduce Stoliarov and
  Walters, Polym. Degrad. Stab. 93 (2008) 422-427. DSC heats of gasification.
- [CAMINO] G. Camino et al., *Polydimethylsiloxane thermal degradation*, Parts 1-2,
  Polymer 42 (2001) 2395 and 43 (2002) 2011. Abstract level only.
- [JANAF] NIST-JANAF table C-002, carbon (reference state, graphite).
- [NISTWB] NIST Chemistry WebBook: n-dodecane ΔfH°liquid -352.1 kJ/mol (Prosen and Rossini
  1945); ethylene ΔfH°gas 52.47 kJ/mol (Chase 1998).
- [JP7] Supercritical pyrolysis of a JP-7-type fuel, J. Anal. Appl. Pyrolysis 2025
  (ScienceDirect S0165237025004656). Abstract level only.
- [NASA-MET] NASA Ames, *Ablation and Heating During Atmospheric Entry and Its Effect on
  Airburst Risk*, NTRS 20180002971 (laser and arc-jet chondrite ablation).
- [LAUB] B. Laub and D. Curry, *Tutorial on Ablative TPS*, NASA Ames 2004, NTRS 20210022472.
- [RUN] Scratch run of this repo's `wall_layers.run` through `plate_thermal`'s pulse and
  maraging substrate, this session. Not committed.

## 1. Orion's ablative

- What it was. GA calls it "a layer of carbonaceous material such as oil (~6 mil thick)",
  applied to the pusher's lower face between pulses [GA-5009-III p. 4, page image read].
  6 mil = 152 µm. This is the child paper's ~150 µm.
- That figure is the applied film, not a measured ablation per pulse. No per-pulse ablated
  depth was found in the open volumes. Status: **not verified** as a measured loss.
- GA-6307 calls it "a relatively thin layer of hydrogeneous material (e.g., oil)". It keeps
  the pusher's rise "within a very small fraction of a degree per pulse", and "any ablation
  is limited to the layer of oil" [GA-6307 §2.3.3].
- "A few thousandths of an inch of oil will keep the metal of the pusher at temperatures
  below 250 °F even after ... a thousand pulses" [GA-6307 §2.7]. GA-5009-III gives the same
  sentence with 600 °F and "a few thousand pulses" [GA-5009-III, ground-test section].
- The ablating layer is "approximately a few thousandths of an inch thick" [GA-6307 §1.2].
- The test plasma stagnated at 150,000 °F (7-8 eV), "high enough that radiation is the
  dominant mode of energy transfer" [GA-6307 §2.7.1]. Stagnation lasted ~100 µs at
  ~100,000 K [GA-5009-III] or ~120,000 K [GA-6307].
- GA built "scaling laws ... which relate the amount of ablation ... to properties of the
  ablative material" and of the gas, via the 1-D code SPUTTER, checked against HE plasma
  shots [GA-5009-III p. 24, GA-6307 §1.3]. The laws themselves are in classified volumes
  (AFWL TDR-64-xx). **Not verified.**
- Graphite loading. Only the later review says "graphite-based grease" [AIAA-2000 p. 5].
  The GA documents read here say "oil" or "carbonaceous material" and give no carbon
  fraction. Orion's 20% carbon is **not verified**. The 20% is the child paper's choice.
- Nance (IEEE Trans. Nucl. Sci. NS-12, 1965) and Dyson (*Death of a Project*, Science 149,
  1965) were not read. Both are paywalled and no open copy was found.

## 2. Thermophysics of candidate oils

| Base fluid | k [W/m/K] | ρ [kg/m³] | c [J/kg/K] | Source |
|---|---|---|---|---|
| PDMS silicone (KF-96, ≥100 cSt) | 0.16 at 25 °C | 965-975 | 1500 at 25 °C | [KF96] |
| PAO (MIL-PRF-87252, Royco 602) | 0.146 at -18 °C, 0.141 at 38 °C, 0.119 at 260 °C | 806 at 0 °C, 739 at 100 °C, 677 at 190 °C | 2050 at -18 °C, 2260 at 38 °C, 2640 at 149 °C, 3010 at 260 °C | [ROYCO] |
| PAO (lab measurement) | 0.143 at room T | - | - | [XU2011] |
| Mineral oil | ~0.13 | ~850-880 | ~1700-1900 | not verified |
| PFPE | ~0.08-0.09 | ~1900 | ~1000 | not verified |

Royco values were converted from cal/(h·cm·°C) and cal/(g·°C) at 1 cal = 4.186 J.

Carbon filler:

- Graphite c = 8.517 J/mol/K = 709 J/kg/K at 298 K, 1217 at 500 K, 1799 at 1000 K [JANAF].
- Carbon-black true density 1800-2100 kg/m³. **Not verified** (no datasheet fetched).
- At 20 wt% in a 850 kg/m³ oil, that is about 10 vol% carbon.

Mixture rules at 20 wt% carbon:

- ρ: 1/ρ = 0.8/ρ_oil + 0.2/ρ_C gives 956 (PAO), 1075 (PDMS), ~1900 (PFPE).
- c: mass-weighted, 0.8·2260 + 0.2·709 = 1950 for PAO at 38 °C. It is 2900 at 260 °C.
  For PDMS it is 1340 at 25 °C.
- k: Maxwell for 10 vol% highly conducting spheres gives k/k0 = (1+2φ)/(1-φ) = 1.33.
  The one measured carbon-black dispersion found gives +38% at 11 wt%
  (0.260 to 0.360 W/m/K in ethylene glycol) [THERMTEST]. Extrapolated linearly, that is
  about +70% at 20 wt%. A percolated black can go higher. That is **not verified** for oils.
- So a 20% carbon PAO film has k ≈ 0.19-0.25 W/m/K. That is the same as pitch's 0.2.
  The thermal diffusivity is also close: oil 0.2/(956·1950) = 1.1e-7 m²/s, pitch 0.96e-7.
- The bulk conduction numbers therefore do not separate oil from pitch. The cap and the
  removal heat do.

## 3. The cap: boil, crack or char?

- Boiling does not apply at the plate. The face sees 0.4-0.9 GPa (handoff ladder), and Orion's
  peak was ~100,000 psi [GA-6307]. Hydrocarbon critical pressures are ~1-2 MPa. The film is
  supercritical, so mass leaves by pyrolysis (cracking), not boiling. Critical pressures are
  textbook values, **not verified** here.
- Slow-heating decomposition peaks (10 K/min DSC): polyethylene 478 °C (751 K),
  polypropylene 447 °C [WALTERS Table 6]. PDMS depolymerises at 400-650 °C (670-920 K) to
  cyclic siloxanes and leaves no solid residue in inert gas below 500-600 °C [CAMINO].
  KF-96 oxidises above 200 °C and autoignites near 450 °C [KF96 §10-12]. At high heating
  rates it gives methane as well [CAMINO].
- Fast heating raises the cracking temperature. Under 10^9-10^10 W/m², the surface recedes at
  v = q/(ρQ) ≈ 0.3-4 m/s. A layer then spends α/v² ≈ 10 ns-1 µs near the surface. Cracking a
  C-C chain (activation ~250 kJ/mol, A ~10^15-10^16 s^-1) in that time needs about
  1200-1600 K. The derivation is mine. The kinetic constants are typical and **not verified**.
- Recommended cap: 1300 K, bracket 900 K (slow-cracking floor) to 1600 K. The user's
  600-1000 K is the low-pressure, slow-heating value. It is too low for this load.
- Silicone is no better as a cap. It cracks in the same window, and its residue is silica,
  not carbon, if oxygen is present. Its oxidation chemistry does not apply in argon or water
  vapor. Its 1.5 kJ/kg/K c is lower. No reason to prefer it was found.
- The carbon's fate. Carbon black in a liquid has no binder. PAO and mineral oil leave no
  char, and PDMS leaves none in inert gas [CAMINO]. So no coherent char layer forms behind
  the cracking front. The particles leave with the cracked gas. Their mass-absorption
  (~7.5 m²/g, handoff) heats a particle at ~10^12-10^13 K/s, so they sublime in the vapor
  layer within nanoseconds. Their heat is spent screening the surface. In this repo that is
  κ_vapor (Rung E), not the wall model's surface removal heat. This is reasoning, **not
  verified** by a high-flux test of a carbon-loaded oil.
- The opposite bracket is a loose soot layer left on the face. It would behave like pitch
  (cap ~3900 K) for as long as it stays. It is carried here as the high case.

## 4. Removal heat

- Polyethylene, DSC: melt 0.218 MJ/kg plus decomposition 0.92 MJ/kg = 1.14 MJ/kg beyond
  sensible heat. The total heat of gasification from 25 °C is 2.51 ± 0.25 MJ/kg. PP is 2.54
  [WALTERS Table 7]. PE is the closest polymer model for a PAO or mineral oil.
- JP-7-type fuel: chemical endotherm 1.19 MJ/kg at 665 °C and 40 bar. Total heat sink
  3.29 MJ/kg including sensible heat [JP7].
- Full cracking ceiling: C12H26(l) → 6 C2H4 + H2 absorbs 6·52.47 + 352.1 = 667 kJ/mol =
  3.92 MJ/kg at 298 K [NISTWB]. Faster heating pushes the product toward this.
- So the oil removal heat beyond sensible heat is 1.1-3.9 MJ/kg, central ~2.
- Per kg of film with 20% carbon carried off cold: 0.8 × (1.1-3.9) = **0.9-3.1, central
  1.6 MJ/kg**. If the carbon instead sublimes at the surface: + 0.2 × 59.3 = **~13.5 MJ/kg**.
- High-flux check. Laser ablation of chondrite at 5-16 kW/cm² gave an effective heat of
  ablation rising from the canonical 8 to 18 MJ/kg. NASA attributes the rise to "radiation
  blockage from ablation products" [NASA-MET slides 8-9]. Q* is a correlation parameter,
  not a material property [LAUB BL-12]. At 10^5-10^6 W/cm², blockage should dominate more.
  No measured Q* for an oil or polymer at 10^9-10^10 W/m² was found. **Not verified.**
- Pitch's own 59.3 MJ/kg is generous. It is graphite's sublimation heat. Real coal-tar pitch
  is roughly half volatile, and that half cracks at polymer-like 1-3 MJ/kg (coking yield
  **not verified**). An honest pitch might be ~30 MJ/kg.

## 5. Recommended set beside pitch

| Quantity | Pitch (current) | Oil, central | Oil low | Oil high | Basis |
|---|---|---|---|---|---|
| k [W/m/K] | 0.2 | 0.20 | 0.14 | 0.30 | §2, [ROYCO][XU2011][THERMTEST], Maxwell |
| ρ [kg/m³] | 1300 | 960 | 900 | 1100 | §2 mixture rule (PFPE ~1900 excluded) |
| c [J/kg/K] | 1600 | 1800 | 1350 | 2600 | §2 mixture rule, [ROYCO][KF96][JANAF] |
| cap [K] | 3900 | 1300 | 900 | 1600 | §3, fast-cracking estimate, not verified |
| removal heat [MJ/kg] | 59.3 | 1.6 | 0.9 | 3.1 | §4, [WALTERS][JP7][NISTWB] |
| soot-layer case | - | cap 3900 K, removal 13.5 MJ/kg | | | §3-4, carbon left on the face |

Use PAO as the base fluid. It has the highest c, is clean and is a qualified coolant [ROYCO].
The 20% carbon belongs in κ_vapor for the shielding run (handoff Step 1b, κ 1500).

## 6. What the model says (scratch run)

Setup: 1,500 pulses, 4 Hz, cold respray only, 1 ms top-hat, maraging 300. The film is 150 µm
for argon and 30 µm for water, and face energies come from `plate_thermal` [RUN]. Pitch
matches the committed `plate_thermal.csv`. Results:

| Film (cap, removal) | Argon steel peak | Argon removed/pulse | Water steel peak | Water removed/pulse |
|---|---|---|---|---|
| Pitch (3900 K, 59.3) | 507 K | 139-162 µm | 991-993 K | 5-30 µm |
| Oil (1300 K, 3.0) | 352 K | 3.6-4.2 mm | 484 K | 0.15-0.79 mm |
| Oil (1300 K, 15) | 352 K | 0.72-0.84 mm | 484 K | 30-157 µm |
| Oil (900 K, 2.5) | 315 K | 4.8-5.6 mm | 392 K | 0.20-1.05 mm |
| Oil soot layer (3900 K, 15) | 534 K | 0.72-0.83 mm | 984-985 K | 26-154 µm |

- The steel peak depends only on the cap and the film's k, ρ and c. Changing the removal heat
  from 3 to 15 MJ/kg left it unchanged to 1 K. The steel heats by soak from the hot film
  after the pulse. A 1300 K surface stores about a third of the heat a 3900 K one does.
- Removal scales as 1/(ρ × removal heat). At the central 1.6 MJ/kg, argon would remove about
  6.8 mm per pulse.
- The solver's fixed grid keeps the film whole while it ablates (ADR-0053). Every oil row
  removes more than the film holds, so its steel numbers assume a film that would be gone.
  In reality the film burns through mid-pulse and the rest of the pulse hits bare steel.

## 7. Expected verdict

The hypothesis is confirmed for heat and refuted for film. With the film intact, a 900-1600 K
cap conducts much less into the steel than pitch at 3900 K: 315-352 K against 507 K for
argon, and 392-484 K against ~990 K for water. The water arm would then need no cooling spray.
But per cubic metre, oil absorbs only 1-3 GJ while it is removed, against pitch's 77 GJ (book
value) or ~39 GJ (if half volatile). Under the unshielded surface flux the 1-D model applies,
the same pulse consumes 5-40x more oil than pitch. That is 0.7-7 mm of oil per argon pulse,
which is impossible. Oil is the better film only if the vapor and soot curtain blocks all but
~3-15% of the radiant energy before it reaches the surface. Orion's GA team believed a thin
hydrogenous film did this. Their own scaling laws tie ablation to the ablator's properties,
and the open volumes do not give them. The decision therefore rests on the Rung E shielded run
with the oil set (cap 1300 K, removal 1.6 MJ/kg, κ_vapor 1500-10^4). Switch to oil if that
run's oil consumption stays below about twice the film thickness per pulse. Until then, pitch
stays the stand-in.

## Blocked and unreached sources

- No host returned the sandbox's "Approval required" or "Blocked by network policy" responses.
- apps.dtic.mil (ADA406019, AFRL endothermic heat sink of hydrocarbon fuels) returned the
  DTIC bot wall ("The request is blocked"). Its Wayback copy returned HTTP 429.
- Nance 1965 (IEEE TNS) and Dyson 1965 (Science) are paywalled. No open copy was found.
- GA's ablation scaling laws and SPUTTER results (AFWL TDR-64-xx, RTD-TDR-63-xxxx) are cited
  in GA-6307 as classified. No open copy was found.
