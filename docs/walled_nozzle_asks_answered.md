# The walled-nozzle asks N9 and N10 — first return

Answers to `katzseth22202/Balloon-Pulse-Propulsion` @ `6faf171`,
`docs/nozzle_asks_for_impact_sim.md` (N9–N11 raised 2026-09-08 against that repo's ADR-0016).
Worked in `katzseth22202/puffsat_impact_simulation` on 2026-09-09; scope is
[ADR-0050](adr/0050-walled-thermal-nozzle-methane-chamber.md).

**This document is written to be copied into the paper repository and worked there.** Every
number in [What the paper should change](#what-the-paper-should-change) is stated in full,
because an agent in the paper repository cannot run the code that produced it. Every figure came
from a run; the [provenance table](#provenance) says which.

---

## Status

| item | verdict | does the paper change? |
| --- | --- | --- |
| **N10.4b** (chamber dissociation equilibrium) | **Answered, and it reverses the premise.** Hydrogen is 93–98% dissociated at the flown 10 kK, not 49–76%. The store is 95–99% charged. | **yes** (W1) — and it gives back most of the Isp the 2026-09-09 correction took |
| **N10, Project 242** | **Answered, and it settles the tension.** 2700 s is *above* the frozen ceiling of 2085–2229 s. The arithmetic is right; the prose is wrong. | **yes** (W2) |
| **N9.0** (sealed vessel) | **Answered, and the geometry question dissolves.** Column length cancels exactly; the verdict is a bore/throat area ratio. But the sound speed used was hydrogen's. | **yes** (W3, W4) |
| **N10.1–3** (freeze stations, `H+H+M`, `N+N+M`) | **Not started.** Needs literature rate coefficients this repository will not invent. | no |
| **N10.5** (carbon nucleation) | **Not started.** | no |
| **N9.1–7** (contact station, wall fluence, throat carbon, convective flux) | **Not started.** | no |

**Four items fall out, W1–W4.** W1 is the one that moves a headline number. **W2 is the strongest
result here**, because it is a published case reproduced and a stated contradiction resolved in
the direction the design needs.

**If you read only two:** W1 (the correction that cost methane 40 s should be largely unwound)
and W2 (Project 242's own number requires the recombination its prose denies).

**One of the asks was itself mistaken** and one of ADR-0016's inputs is wrong; both are in
[What the ask got wrong](#what-the-ask-got-wrong), because one of them would have sent N9 item 0
at a target that is not there.

---

## What the paper should change

### W1. The chamber's hydrogen is essentially fully dissociated, and the 2026-09-09 correction overshot

**Locate:** ADR-0016, "with the dissociation equilibrium solved rather than assumed" and
"Hydrogen runs 55 to 76% dissociated across these cases"; and N10 item 4b, "only 49 to 76%
dissociated".

**Now:** the chemical store is taken as partly charged, which moved methane from 1132 s to
1059–1095 s.

**Should be:** at the flown 10,000 K the store is charged to 95–99%. Solved on a from-scratch
nine-species equilibrium EOS (`CH4, CH3, CH2, CH, C2H2, C2H, C2, C3, H2`, atomic `C` and `H`,
`H+`, the `C+..C6+` Saha ladder, `e-`) with NIST-JANAF 0 K thermochemistry:

| chamber | `rho` | `p` [bar] | `u` [MJ/kg] | H dissociated | C free | **store charged** | `k` |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 200 m³ | 2.570 | 642 | 136.8 | 0.934 | 0.902 | **0.947** | 19.56 |
| 400 m³ | 1.257 | 321 | 139.9 | 0.966 | 0.961 | **0.976** | 19.11 |
| 673 m³ | 0.740 | 190 | 141.2 | 0.980 | 0.980 | **0.986** | 18.92 |

All at 10,000 K. **`k` is an output, not an input**: it is the fixed point of ADR-0016's own
identity `(1+k) u = w²/2`, and at 673 m³ it returns `k = 18.92` against the ask's stated 18.96 —
which is the check that both repositories are using the same energy ledger.

**Why the paper-side estimate came out low.** `H2 ⇌ 2H` at 10,000 K has `D0/kT = 5.2`, so the
Boltzmann factor barely suppresses it: the equilibrium constant `n_H²/n_H2` is `1.22e28 m⁻³`.
Half dissociation needs `n(H nuclei)` equal to that constant, which for methane is **81.5 kg/m³
and ~14,000 bar** — thirty times the densest chamber on the table, and the full solve there
returns 0.458. Le Chatelier pushes in the direction the paper argues; at 600 bar it is nowhere
near strong enough to do what was booked.

**What to do with the Isp.** On ADR-0016's own exhaust-speed identity the impulse per kilogram of
launched slug goes as `sqrt(1+k)/k`, so the ratio survives every normalisation the paper's
effective-Isp column applies. Against the ask's own `k = 18.96` geometry:

| chamber | 8,000 K | 10,000 K | 12,000 K |
| --- | ---: | ---: | ---: |
| 200 m³ | 0.878 | 0.984 | 1.035 |
| 400 m³ | 0.918 | 0.996 | 1.043 |
| 673 m³ | 0.940 | **1.001** | 1.049 |

At the flown point the three chambers are within 1.8% of each other and within 1.6% of the
ask's geometry. **The 200 → 673 m³ gain is 1.8%, not the +3.4% ADR-0016 books.**

**Where the effective Isp lands depends which of ADR-0016's own rows you anchor on, and the two
do not agree** — which is itself worth knowing, because it says the effective-Isp column is not a
pure function of `k`:

| anchor | its `k` | 200 m³ | 400 m³ | 673 m³ |
| --- | ---: | ---: | ---: | ---: |
| pre-correction, full atomisation, 1132 s | 18.67 | 1105 s | 1118 s | 1124 s |
| post-correction, 1080 s at 400 m³ | 21.59 | 1137 s | 1151 s | 1157 s |

The first is the like-for-like comparison — both cases have the store fully charged, and it
differs from 1132 s only because the solved `u` is 139.9 MJ/kg rather than the assumed 143. **Take
1105–1124 s**, and treat the second row as evidence that something other than `sqrt(1+k)/k` is in
the paper's column and should be identified before either number is printed.

Either way the direction is the same: **up from the current 1059–1095 s**, by 2–6%. Subject
throughout to N10.1–3 and N10.5, which are not answered here and which decide whether the charged
store is *returned* at all.

**Two things this does not license.** The 8,000 K row is a real 6–12% penalty, so the argument
for 10,000 K over a cooler chamber is now stronger than the argument against a hotter one. And a
charge is not a return: `store charged = 0.976` is a ceiling on what recombination could give
back, not a claim that it does.

### W2. Project 242's 2700 s cannot be reached frozen — its arithmetic is right and its prose is wrong

**Locate:** ADR-0016, "Their prose says molecular recombination 'is not' fast. Their arithmetic
says otherwise"; N10, "Those two statements are in tension and the solver should say which is
right."

**Should be:** stated as settled. Equilibrium hydrogen at 10,000 K, solved on the same frame:

| `p` | dissociated | `h` [MJ/kg] | `h` frozen | Isp equilibrium | **Isp frozen** | store return 2700 s needs |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 bar | 0.9999 | 453.1 | 238.8 | 3070 s | **2229 s** | 0.52 |
| 3 bar | 0.9997 | 439.3 | 225.0 | 3022 s | **2163 s** | 0.59 |
| 10 bar | 0.9988 | 430.5 | 216.5 | 2992 s | **2122 s** | 0.63 |
| 100 bar | 0.9884 | 420.9 | 209.1 | 2959 s | **2085 s** | 0.67 |

Both columns are *ceilings on a perfect nozzle* — everything converted, nothing left in the
exhaust — which is what makes the comparison decisive in one direction. **2700 s is 24–29% above
the frozen ceiling at every pressure in "a few bar".** No nozzle, however good, reaches it
without recombination. It sits comfortably inside the equilibrium ceiling of ~3000 s, so the
number is not absurd either: it requires **52–63% of the held dissociation store to return**,
which is where the ask's own 67% estimate lands.

The ask's supporting figures are right in shape and 5–10% low: it carries 206 MJ/kg for
sensible-only hydrogen against a solved 217–239, and 421 for full recovery against 421–453. The
351 MJ/kg that 2700 s demands is exact.

**Why this matters more than a validation.** ADR-0016 rests its whole case on running 139× the
number density Rubbia does. That argument needs Rubbia's own case to *have* recombination in it,
otherwise the density scaling starts from zero. It does.

### W3. Item 0's answer does not depend on the column length, so it cannot pick the geometry

**Locate:** ADR-0016, "Until it returns, 400 m³ over 14 m is the conservative middle"; N9 item 0,
"it changes the geometry ... Run 200 m³ over 7.1 m, 400 m³ over 14 m and 673 m³ over 23.8 m".

**Should be:** the equilibration count is `V / (L A* f)`, and `V = pi r² L` at a fixed bore, so

> **n_eq = A_bore / (f A*)** — the column length cancels exactly.

(`f = (2/(gamma+1))^((gamma+1)/(2(gamma-1)))` is the standard choked-flow coefficient, 0.65 at
`gamma = 1.2`.) The run confirms it: 23.8, 24.1 and 23.9 crossings at a 2 m² throat for the
three volumes; 6.8, 6.9 and 6.8 at 7 m². **Item 0 is a throat question, not a length question**,
and the short column is admissible on exactly the same terms as the long one. The geometry
choice returns to the thermal and Isp arguments — which W1 shows are worth under 2% between
200 and 673 m³ at 10,000 K.

The verdict itself holds but with less margin than booked: `n_eq ≈ 24` at a 2 m² throat,
`≈ 6.8` at 7 m². Not ~50.

### W4. The 11 km/s sound speed is hydrogen's, not methane's

**Locate:** ADR-0016, "at an 11 km/s sound speed a 7 m column equilibrates some fifty times".

**Now:** 11 km/s. **Should be:** 6.0–6.3 km/s.

Dissociated methane at 10,000 K averages **3.33 amu per particle**, not 1: each `CH4` gives one
carbon and four hydrogens, so `16.043/5 = 3.2`. Atomic hydrogen at 10,000 K is
`sqrt(5/3 kT/m_H)` = **11.7 km/s**, which is where the number came from. The solved equilibrium
sound speed is 6036 m/s at 200 m³, 6200 at 400 m³ and 6253 at 673 m³. The equilibration count is
linear in it, which is half of why the ~50 became ~7–24.

**A separate consequence, and it protects N9 items 1–3.** The impactor crosses the column at
75 km/s against a 6 km/s sound speed, so it completes **0.08 of one acoustic crossing** while it
is inside. Equilibration sets `k` — it happens over the blowdown, long after — but it does
**not** soften the wall's arrival transient. The cone problem `sec:needle_through_fog` describes
is untouched by the sealed vessel, and N9 items 1–3 remain exactly as posed.

---

## Smaller corrections, no argument attached

- **The 1.4 bar pre-charge belongs to the 200 m³ chamber alone.** The same methane mass sits at
  0.70 bar in 400 m³ and 0.41 bar in 673 m³ at 111 K. ADR-0016 quotes 1.4 bar without a volume.
- **Atomisation energy: 102.35 MJ/kg, not 103.7; formation enthalpy 4.15 MJ/kg, not 4.66.** The
  gap is the reference state — ADR-0016 used 298 K heats of formation and this side uses 0 K,
  matching the partition functions. Not a disagreement; state whichever, but state which.
- **The 52% / 43% split is confirmed exactly.** Methane returns 52.6% of atomisation as `H2` and
  43.3% as condensing carbon, with the remaining 4.2% the formation enthalpy that never returns.
  Derived here from heats of formation alone, so it is an independent check of both sides.
- **The sensible fraction is confirmed.** ADR-0016's 0.273 against a solved 0.277–0.287 at
  10,000 K, so the argument against a hotter chamber stands on a number this repository now owns.
- **Chamber pressures are confirmed.** 642 bar at 200 m³ against ADR-0016's 636; 190 bar at
  673 m³ against N10's 206.

## What the ask got wrong

**N9 item 0 is posed as a geometry question and it is not one.** "It changes the geometry ... Run
200 m³ over 7.1 m, 400 m³ over 14 m and 673 m³ over 23.8 m" expects the three to separate. They
cannot: the length divides out of the equilibration count analytically (W3). Running the three
was still worth doing — it is how the cancellation was found — but the item cannot deliver what
it was asked to deliver, and the geometry has to be chosen on the thermal and Isp arguments
instead.

**N10 item 4b's premise is the thing that needed checking, and it does not survive.** The item is
written as "report the equilibrium composition, we already know it is 49–76%". The composition is
reported and it is 93–98% (W1).

## Deferred, with cost

- **N10 items 1–3 (freeze stations and frozen fractions for `H+H+M` and `N+N+M`).** Not started,
  deliberately. `recombination.py` states that it "computes no rate coefficients of its own" and
  every constant in it is a named literature value with a stated uncertainty. Supplying `H+H+M`
  and `N+N+M` coefficients from recall would break that discipline on the exact numbers the
  1080 s / 793 s fork turns on. **What is needed to unblock it:** evaluated low-pressure
  three-body coefficients and their temperature exponents for `H + H + M` and `N + N + M`, with
  third-body efficiencies, from Baulch et al. or Park's high-temperature set. Everything else —
  the Bray criterion, the expansion history, the Damköhler comparison, the decades-of-margin
  reporting — already exists and needs only a no-field variant of `nozzle_history`.
  **Cost of the delay:** W1 raises the *charged* store; only this decides how much comes back.
- **N10 item 5 (carbon nucleation).** Not started. It is not a three-body reaction and needs
  classical-nucleation machinery that does not exist here. It is 43% of the store.
- **N9 items 1–7.** Not started. W4 confirms they are still live: the sealed vessel does not
  dissolve the arrival transient.
- **N11.** Not started.

## Provenance

| number | produced by | command |
| --- | --- | --- |
| chamber composition, `k`, `u`, pressures, Isp ratios, item 0 timescales | `puffsat.walled_nozzle.chamber` | `make walled-nozzle-chamber` |
| Project 242 bracket, pure-hydrogen rung | `puffsat.walled_nozzle.hydrogen` | `make walled-nozzle-hydrogen` |
| the equilibrium EOS underneath both | `puffsat.eos_methane`, `puffsat.walled_nozzle.hydrogen` | `make walled-nozzle-test` |

Committed artifact: `data/results/walled_nozzle/chamber.csv`.

## Known weaknesses of this answer

Stated rather than buried, in the order they would bite:

1. **`C3`'s bend.** `C3` is treated as a harmonic rotor-oscillator and its `nu2` bend at 63 cm⁻¹
   is quasilinear, so the harmonic partition function (≈90 per mode at 8 kK) certainly overstates
   it. `chamber.csv` carries `c3_store_exposure`, the *exact upper bound* on what this costs — all
   of the store `C3` is withholding. At the flown 10,000 K it is **0.1–1.2%**; at 8,000 K it is
   **3.4–9.9%**. The 10,000 K conclusions are safe; the 8,000 K row is soft in the direction of
   understating the charge.
2. **Ground electronic terms only**, as in `eos_water`. `C2`'s low-lying `a 3Pi_u` at 716 cm⁻¹ is
   not summed, understating `C2`, which never exceeds a few percent of the carbon here.
3. **Ideal mixture.** `p = n k T` with no excluded volume and no Coulomb correction. At the
   densest chamber the electron coupling parameter is `Gamma ≈ 0.25` and the excluded volume is
   ~1%, so this is a percent-level effect at 200 m³ and less elsewhere — the same assumption
   `eos_water` makes, and the same place it would first fail.
4. **Condensed carbon is absent by construction** (ADR-0050). Below ~4,000 K the true equilibrium
   has graphite in it and this EOS does not, so the cold end of the expansion is not yet a
   physical answer. That is N10 item 5 and it is not answered here.
