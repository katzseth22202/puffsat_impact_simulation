# Facesheet atomic/displacement damage is out of scope: per-atom energy is sub-threshold

A reviewer will ask whether hypervelocity gas atoms crystal-damage the SiC facesheet. They do not,
anywhere in the envelope, for a quantitative reason — so displacement/penetration/sputtering damage
is explicitly out of scope, and the survivability frontier is thermal/structural-limited only.

**Per-atom energy is eV-scale, at or below the SiC displacement threshold.** The damaging species is
the heaviest abundant gas atom — oxygen (16 amu). Its kinetic energy `½mv²`:

| v (km/s) | O atom KE |
|---|---|
| 3.2 | 0.85 eV |
| 11 | 10 eV |
| 16 | 21 eV |

(Hydrogen is ~16× lower; 1.3 eV even at 16 km/s.) The SiC displacement threshold is `E_d ≈ 20–35 eV`
(C ~20, Si ~35). Elastic transfer O→C/Si is near-maximal (`γ = 4m₁m₂/(m₁+m₂)² ≈ 0.93–0.98`), so even
at 16 km/s the transferred energy (~20 eV) only marginally reaches the carbon threshold and never the
silicon one; below ~13 km/s it is sub-threshold outright. Gas atoms cannot crystal-damage the SiC.

**Contrast with Project Orion.** Orion's tungsten jets (184 amu) at ~100 km/s carry ~9.5 keV per atom
— ~10³× the PuffSat per-atom energy, far above any displacement/sputter threshold. Orion was genuinely
in a sputtering/displacement regime; PuffSat is not. Both lower atomic weight (H, O vs W) and lower
velocity contribute, the `v²` term dominating, so `m·v²` lands ~3 orders of magnitude below Orion on
the benign side of the threshold.

**A few-micron ablator stops the atoms with enormous margin.** A ~20 eV atom has a stopping range
< 1 nm; even a 21 keV atom (O at the out-of-scope near-Sun ~500 km/s) ranges only ~tens of nm. An atom
would need ~MeV (oxygen at ~thousands of km/s) to penetrate microns of ablator. Penetration is never
the limiting factor in or beyond the envelope; the binding constraint on the few-micron ablator is
**thermal burn-through per pulse** (ablation depth < layer thickness) — the ablation-per-pulse
requirement reported at the frontier — not atomic penetration.

**The real SiC threats are chemical and thermal, both already handled (§5).** Oxidation by dissociated
atomic O (→ the SiC choice and SiO₂ passivation) and transient thermal shock (→ the pressure-limit
table). The ablator's role is vapor/thermal shielding plus a sacrificial oxidation surface, not
atom-stopping.

## Consequence

Crystal displacement, sputtering, and penetration damage are not modeled. The survivability frontier
is limited by peak facesheet pressure (structural) and ablation burn-through (thermal/logistical) only.
