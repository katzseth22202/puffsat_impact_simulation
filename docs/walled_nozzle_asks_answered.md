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
| **N10.1–3** (freeze stations, `H+H+M`, `N+N+M`) | **Answered, and the fork closes on the good side.** The gas never freezes in this nozzle — `Da` stays above the threshold everywhere, with 0.4–2.2 decades of margin. The walled nozzle converts **39–41% of the wall-cap energy into directed kinetic energy** at the ask's own 7 m² throat (methane; no ammonia Isp is produced — see W7). That is *not* the paper's effective-Isp column and must not be compared with 1080 s directly. | **yes** (W5, W6, W7) |
| **N10.3** (the water and hydrogen rungs) | **Answered.** The ask's own `k = 37.70` and `k = 7.99` are reproduced as *outputs* (39.55, 7.77). Neither water nor ammonia beats methane; hydrogen nearly doubles it. | **yes** (W8) |
| **N10.5** (carbon nucleation) | **Reframed, and much smaller than booked.** 84% of the carbon store returns as gas-phase acetylene with no nucleation at all; only 16% of it (5.2% of the energy budget) is genuinely hostage to soot. | **yes** (W9) |
| **N9.1–7** (contact station, wall fluence, throat carbon, convective flux) | **Not started.** | no |

**Nine items fall out, W1–W9.** W1 is the one that moves a headline number. **W2 is the strongest
result here**, because it is a published case reproduced and a stated contradiction resolved in
the direction the design needs. **W5 settles the 1080 s / 793 s fork** the ask calls high priority,
and **W6 says the remaining lever is the throat, not the chemistry**.

**If you read only three:** W1 (the correction that cost methane 40 s should be largely unwound),
W2 (Project 242's own number requires the recombination its prose denies), and W5 (the walled
nozzle does not freeze, so the *equilibrium* branch applies and the frozen one does not — read its
note on Isp conventions before quoting any number from it).

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
throughout to N10.5, which is not answered here. N10.1–3 *is* now answered (W5) and it confirms
the charged store is genuinely returned rather than frozen out — but note that W5's exit Isp is
computed directly from the expansion rather than from `sqrt(1+k)/k`, and lands at 1074–1076 s at
the ask's own throat — but on a *different normalisation*, so **the two must not be compared
directly**; W5 sets out why. The `sqrt(1+k)/k` ratio above is the safe way to move a number in
this column, and it is what this section uses.

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

**And at the wide end of the ask's own throat range the verdict does not hold at all.** The
criterion this study applies is `n_eq >= 10` — a column that turns over ten times has no memory of
which parcels the cone touched first. The run gives `n_eq ≈ 24` at a 2 m² throat, `11.9` at 4 m²
and **`6.8` at 7 m², which fails it**. Not ~50 anywhere, and not passing at all at 7 m².

So item 0's answer is conditional on the throat: **the sealed vessel sets `k` only for throats of
about 4 m² or narrower.** At 7 m² the chamber empties before it has equilibrated and `k` is once
again set by what the cone sweeps — the very coupling problem the walled chamber was adopted to
escape. This turns out to point the same way as W6's Isp argument, which is the strongest reason
in this document to prefer a narrow throat.

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

### W5. The plume does not freeze in this nozzle, so 1080 s is the right branch and 793 s is not

**Locate:** N10, "For the flown methane slug it is the fork between 1,080 s and 793 s", and its
items 1–3; ADR-0016's hand-estimated margins of ×3400 for `H+H+M` and ×403 for `N+N+M`.

**Now:** the fork is open, and the equilibrium branch is an assumption.

**Should be:** the fork is closed in favour of equilibrium. Running Bray's criterion station by
station along the de Laval expansion, on the literature `H + H + M` coefficients tabulated
below, the Damköhler number never falls to the freezing band anywhere in the bore:

**All of the following is methane** (`eos_methane`, ADR-0016's 25 kg at 75 km/s into CH₄).
Ammonia is not run here at all; W7 says why. **The `Isp` column is a *total*, unnormalised
exhaust-velocity Isp and is not the paper's effective-Isp column** — see the note under the table
before comparing it with anything.

| chamber | throat | `A/A*` | exit `T` [K] | min `Da` | verdict | store returned | **exit Isp [s]** | margin [decades] |
| ---: | ---: | ---: | ---: | ---: | :--- | ---: | ---: | ---: |
| 200 m³ | 7.0 m² | 4.04 | 5584 | 1461 | equilibrium | 0.262 | **1074.3** | 2.16 |
| 400 m³ | 7.0 m² | 4.04 | 5425 | 815 | equilibrium | 0.243 | **1076.3** | 1.91 |
| 673 m³ | 7.0 m² | 4.04 | 5296 | 343 | equilibrium | 0.226 | **1073.8** | 1.53 |
| 200 m³ | 4.0 m² | 7.07 | 5067 | 732 | equilibrium | 0.313 | 1153.6 | 1.86 |
| 400 m³ | 4.0 m² | 7.07 | 4931 | 341 | equilibrium | 0.294 | 1155.5 | 1.53 |
| 673 m³ | 4.0 m² | 7.07 | 4822 | 134 | equilibrium | 0.278 | 1152.7 | 1.13 |
| 200 m³ | 2.0 m² | 14.14 | 4561 | 216 | equilibrium | 0.367 | **1229.7** | 1.33 |
| 400 m³ | 2.0 m² | 14.14 | 4443 | 65.2 | equilibrium | 0.348 | **1232.0** | 0.81 |
| 673 m³ | 2.0 m² | 14.14 | 4349 | 25.8 | equilibrium | 0.332 | **1229.1** | 0.41 |

**`Isp` here is the nozzle's own exhaust velocity over `g0`, and it is NOT the paper's
effective-Isp column.** The distinction matters and an earlier draft of this document got it
wrong, so it is spelled out:

- **What is quoted above is *total*** — exit flow speed divided by `g0`, counting every kilogram
  expelled (slug *and* vaporised impactor), with **no** divergence or geometry loss, no drift
  term and no vessel mass.
- **ADR-0016's 1080 s is *effective*** — impulse per kilogram of *launched slug*, which is a
  factor `(1+k)/k = 1.051` larger for the same exhaust speed, and it additionally carries the
  launch-ledger normalisations this repository does not own (`eta_geom`, the drift term of
  `eq:reflection_baseline`, vessel mass). Those are worth a further **factor 0.609 on velocity**:
  ADR-0016's own identity `w/sqrt(1+k)` gives 16 541 m/s, which on its slug convention would be
  1773 s, and the column says 1080 s.

Put on the paper's slug convention our exhaust speeds give 1129–1133 s at 7 m², but that number
still lacks the 0.609 of normalisation, so **it must not be read as reproducing 1080 s either.**
Any near-agreement between an unnormalised total and a normalised effective is two errors
cancelling, and this document should not trade on it. `chamber.isp_scaling` refuses to emit an
absolute Isp for exactly this reason and that discipline applies here too.

**The convention-free statement of the same result is the energy conversion fraction**, which is
what this study actually owns:

| throat | wall-cap energy appearing as directed KE |
| ---: | ---: |
| 7.0 m² | **0.393–0.406** |
| 4.0 m² | 0.452–0.468 |
| 2.0 m² | **0.514–0.532** |

At the ask's own throat the walled nozzle turns about **40% of `u` into directed kinetic energy**,
and that figure carries no normalisation of anyone's. It is the number to compare against, and
the ratios in W6 are safe because a constant normalisation cancels out of them.

**What the fork verdict does and does not need.** The 1080 s / 793 s choice is a question about
whether the composition freezes, and that is settled by the Damköhler columns above without any
Isp convention entering. **The equilibrium branch is self-consistent everywhere in the bore, so
793 s does not apply** — that conclusion stands. What this study cannot yet do is quote the
paper's own 1080 s back at it, because that would need the paper's normalisations. **The clean
like-for-like test would be the equilibrium-to-frozen *ratio*** — ADR-0016's own 1080/793 = 1.36 —
run on a frozen-composition methane EOS, which `eos_water` has (`pressure_energy_frozen`) and
`eos_methane` does not yet. That is the first thing to build next on this rung.

**The margin is the part to keep.** `margin` is how many decades the rate coefficient could be
wrong, in the pessimistic direction, before the verdict changes. The stated uncertainty on the
evaluated coefficient is 0.5 decades (a factor of 3.16), so **the flown 200 m³ / 7 m² case, at
2.16 decades, survives its own rate uncertainty roughly four times over.**

**One corner does not, and should be said out loud: 673 m³ at a 2 m² throat has 0.41 decades,
which is *less* than the coefficient's own uncertainty.** That configuration — the largest, and
therefore thinnest, chamber pushed to the largest area ratio — is the one place in this study
where the verdict genuinely rests on the rate constant being right. It is also the corner W6
recommends moving toward, so the two results have to be read together: narrow the throat on the
200 m³ chamber, where the margin is 1.33 decades, rather than on the 673 m³ one.

ADR-0016's hand-estimated ×3400 (3.5 decades) is too generous even for the best case here, and
its ×403 is the wrong reaction for methane entirely — but both erred in the safe direction.

**Where it would eventually freeze.** Pushed past the bore as a diagnostic, the 400 and 673 m³
cases do freeze — at `A/A*` of 188 and 95, at 3300–3500 K and around 0.1 bar, with about half the
store still held. So there is a real ceiling: no amount of nozzle recovers all of it. The 200 m³
case does not freeze out to `A/A* = 400` at all, being the densest.

### W6. The binding constraint is the expansion ratio, not the chemistry — narrow the throat

**Locate:** N10 item 4, "Sensitivity to throat area, since that is the knob that sets how long the
gas stays dense."

**Now:** the item is posed as a chemistry-preservation knob — hold density up so recombination
keeps up.

**Should be:** that is the right knob for the wrong reason, and it points the other way. The
chemistry is not the thing in short supply: even the *worst* case in W5 has 0.4 decades of margin,
and the flown one has 2.16. What limits the return is that the nozzle **does not expand far
enough** — at `A/A* = 4.04` the gas leaves at 5300–5600 K still holding 70–76% of its store, not
because it froze but because at 5500 K equilibrium itself still holds the bonds broken.

So narrowing the throat helps, and substantially:

| throat | `A/A*` | exit `T` [K] | exit Isp [s] — *total*, methane, unnormalised (see W5) | vs the 7 m² baseline | `n_eq` (N9.0) | sets `k`? | blowdown, 200 m³ |
| ---: | ---: | ---: | ---: | ---: | ---: | :--- | ---: |
| 7.0 m² | 4.04 | 5296–5584 | 1073.8–1076.3 | — | 6.8 | **no** | 8.0 ms |
| 4.0 m² | 7.07 | 4822–5067 | 1152.7–1155.5 | **+79 s (+7.4%)** | 11.9 | yes | 14.0 ms |
| 2.0 m² | 14.14 | 4349–4561 | 1229.1–1232.0 | **+155 s (+14.4%)** | 23.8 | yes | 28.0 ms |

**Why not go the other way and widen it?** Because the bore is fixed at 28.3 m², so
`A/A* = A_bore / A*` and the throat is the *only* expansion-ratio knob available inside a 3 m
bore. Widening shrinks the expansion, and a de Laval nozzle converts enthalpy to directed speed
only in proportion to how far it opens. What a wider throat buys is chemical margin — the one
quantity already in surplus by two decades. It is paying, in the currency that is limiting, for
more of the commodity that is not.

**Three independent arguments converge on the same recommendation**, which is what makes it worth
acting on: the Isp gain above; N9 item 0's equilibration count, which *fails* at 7 m² and passes
at 4 m² (W3); and the freeze margin, which stays comfortable at 4 m² in every chamber.

**The cost is blowdown time, and it is not costed here.** Choked mass flow is proportional to `A*`,
so the pulse stretches from 8.0 ms to 28.0 ms at 200 m³ as the throat goes 7 → 2 m². That is 3.5×
longer for the wall to absorb the same pulse, at the same time as the throat itself is passing the
same power through a third of the area. **Both are N9 items 1–7 and neither is answered**, so the
Isp gain quoted here is an upper bound on what is actually collectable.

**This is the largest single lever found in either ask**, and it costs no new physics — only a
smaller throat and whatever that does to the wall loading, which is N9's problem. The chemistry
still keeps up at 2 m² in every chamber — but read W5's caveat before acting on the bottom row:
at 2 m² the 673 m³ chamber's margin falls to 0.41 decades, below the rate coefficient's own
uncertainty, while the 200 m³ chamber still holds 1.33. **Narrow the throat on the small chamber,
not the large one**, and note that the 2 m² exits at 4349–4561 K walk into the temperature range
where this EOS's omission of condensed carbon starts to matter (weakness 4).

**A second result falls out of the same table: chamber volume barely moves the specific impulse.**
At fixed throat the three chambers land within 0.3% of each other (1073.8 / 1076.3 / 1074.3 at
7 m²). A larger chamber charges more of the store (W1) but is thinner and returns less of it by
the exit plane, and the two effects cancel almost exactly. Taken with W3 — where the sealed-vessel
count also turned out not to separate the geometries — **N9's volume choice should be made on
thermal and structural grounds, because the Isp argument does not distinguish them.**

### W7. The ammonia rung rests on a rate that is uncertain by 1.3 decades, and carbon has no rate at all

**Locate:** N10, "Nitrogen is about eight times slower, which is why it freezes in arcjets at 0.1
to 1 bar", and item 2's 68.9 MJ/kg for ammonia; item 5's carbon condensation.

**Now:** the 8× is quoted as a settled ratio.

**Should be:** quoted with its spread, because the nitrogen channel is the *least* well determined
number in this ask. Comparing atomic third bodies at 6000 K:

| source | `k(N+N+N)` [m⁶/s] | `k(H+H+H)/k(N+N+N)` |
| --- | ---: | ---: |
| Byron 1966 (shock tube, 6000–9000 K) | 1.43e-44 | **0.6** — nitrogen *faster* than hydrogen |
| Notey, Jo & Panesi 2025 (ab initio master equation) | 7.24e-46 | **12.7** |

**The two disagree by 1.3 decades.** The ask's 8× sits inside that band but only near the ab
initio end; Byron's measurement contradicts the claim outright. The statement is defensible, not
robust, and should not be carried without the spread — especially as ammonia's whole rung depends
on it. **An ammonia specific impulse is not produced here**: that needs an `eos_ammonia` this
repository does not have, and building one on a rate this uncertain would be misleading precision.

**Carbon is a different and harder situation, and the literature search is itself the finding.**
The evaluated databases hold exactly one measurement of `C + C + M -> C2 + M` (one third body, a
1000 K window, ±55%) and **nothing at all** for `C + H + M -> CH + M`, against six independent
evaluations and three shock-tube studies for hydrogen. An association channel gets measured when
it controls something measurable, and carbon's do not — because carbon in a cooling carbon-rich
gas does not form diatomics, it forms clusters and then soot. **So the 43% of methane's
atomisation locked in carbon cannot be settled by a rate coefficient the way the 52% held as H₂
just was.** N10 item 5 needs classical nucleation theory, and no amount of further rate-hunting
substitutes for it.

### W8. The propellant ladder, solved: neither water nor ammonia beats methane, and hydrogen doubles it

**Locate:** N10 item 3, "The same for a water slug at `k = 37.70` and for pure hydrogen at
`k = 7.99`, so the ladder in ADR-0016 rests on solved chemistry rather than on the equilibrium
assumption."

**Now:** the ladder's slug ratios are quoted as inputs.

**Should be:** they are outputs, and they come out close to the quoted values. The chamber
identity mentions no species —

> `rho = (1+k) m_p / V`,  `u = e(rho, T)`,  `k = (w²/2)/u − 1`

— so pointing it at a different EOS gives the whole ladder. Run at 10 000 K, 75 km/s, through the
ask's own 7 m² throat, with each fluid on **its own** equilibrium EOS so the conversion fraction
is solved rather than borrowed:

`Isp total` is the exhaust speed over `g0`; `Isp effective` is that same speed multiplied by
`(1+k)/k` — so the effective column is *derived* from the total, not independently computed.
Neither carries the paper's launch-ledger normalisations (W5).

| fluid | `m̄` [amu] | `k` | `u` [MJ/kg] | exit `T` [K] | conversion | `u_e` [m/s] | Isp total | **Isp effective** |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **hydrogen** | 1.01 | **7.77** | 320.8 | 5168 | 0.532 | 18 482 | 1885 s | **2127 s** |
| **methane** | 3.21 | 19.56 | 136.8 | 5584 | 0.406 | 10 536 | 1074 s | **1129 s** |
| ammonia\* | 4.26 | 28.54 | 95.2 | — | *assumed* | ~8 800 | 833–944 s | **862–977 s** |
| **water** | 6.00 | **39.55** | 69.4 | 5200 | 0.415 | 7 588 | 774 s | **793 s** |

\* **Estimate, not a solved rung, and the estimate is load-bearing.** There is no `eos_ammonia`
— nitrogen appears in no EOS in this repository — so ammonia's `u` is assembled from the ask's own
68.9 MJ/kg atomisation plus an exactly computed translational term, and **its conversion fraction
is assumed to be methane-like (0.35–0.45) rather than solved.**

**That assumption is probably wrong in ammonia's favour, so the methane-over-ammonia ordering is
NOT established.** The three solved rungs each returned a different share of their store, and the
reason methane converts only 0.406 is W9: at the exit its carbon is parked in C₃ with just 30% of
atomisation handed back. Ammonia has no C₃ trap — its store returns as **N₂ (941 kJ/mol, among
the strongest bonds in chemistry) and H₂**, two simple diatomics — so a materially higher
conversion is physically plausible.

**The crossover is computable even though the rung is not:** ammonia draws level with methane at a
conversion of **0.601**, against methane's solved 0.406 — a factor of 1.48. That is a large ask,
but it is not absurd: **hydrogen reaches 0.532 at this same area ratio**, so conversions above
0.5 are demonstrably attainable here. Until an `eos_ammonia` exists, read this table as

> **hydrogen > methane > water, solved** — and ammonia unplaced, somewhere between water and
> possibly above methane.

The `1/sqrt(m̄)` scaling below cannot settle it either, because that law implicitly assumes a
common conversion fraction and this is exactly the case where that fails.

**Two of these are checks rather than results, and they pass.** Nothing in the fixed point was
told the ask's numbers, and it returns **`k = 39.55` against the stated 37.70 for water** and
**`k = 7.77` against 7.99 for hydrogen** — both within 5%. The ladder is therefore resting on
solved chemistry, which is what item 3 asked for.

**The ordering is one variable, and it is not negotiable by nozzle design.** Effective Isp tracks
`1/sqrt(m̄)`, the mean mass of a particle once atomised, to within 10% across a factor of six:

| | H₂ | CH₄ | NH₃ | H₂O |
| --- | ---: | ---: | ---: | ---: |
| atomises to | 2 H | C + 4H | N + 3H | O + 2H |
| `1/sqrt(m̄)` ÷ methane | 1.78 | 1.00 | 0.87 | 0.73 |
| **solved Isp ÷ methane** | **1.88** | 1.00 | 0.86 | **0.70** |

A light fluid stores more per kilogram, because at 10 kK most of `u` is `3/2 kT` **per particle**
and a kilogram of a light gas is more particles. Storing more per kilogram means less mass is
needed to absorb the pulse, so `k` falls — and then `k` enters a *second* time through the
effective convention `(1+k)/k`, which is 1.13 for hydrogen and 1.03 for water. Low `k` is rewarded
twice. This is the molecular-weight law of any rocket, reappearing in a chamber heated by impact
rather than by combustion.

**So the carbon exposure cannot be dodged by switching to a carbon-free fluid without paying for
it** — with ammonia the open question above. Water is carbon-free and 30% worse. Only hydrogen is
both carbon-free and better, by 88% — and its store returns through `H + H + M`, the one channel
W5 has established with margin to spare. **ADR-0016 treats hydrogen as a validation rung; on
these numbers it is the strongest candidate in the table**, and its real cost is storage density,
which is a vehicle problem that appears nowhere in these asks.

### W9. Most of the carbon store returns as acetylene, not as soot — which moves N10 item 5 from a nucleation problem to an expansion-ratio problem

**Locate:** N10 item 5, "Methane returns 52% of its atomisation as H2 and 43% as condensing
carbon", and its framing of the 43% as recoverable only by nucleation.

**Now:** the whole 43% is treated as hostage to soot formation, with ADR-0016's 880 s floor if
carbon does not condense.

**Should be:** the great majority of that 43% comes back **in the gas phase**, before any particle
forms, because acetylene is nearly as strongly bound per atom as methane was. The stoichiometry is
exact and needs no rate constant at all:

> **2 CH₄ → C₂H₂ + 3 H₂** returns **89.0%** of full atomisation.

against 52.6% for the H₂ channel alone and 95.9% for H₂ plus fully condensed carbon. All three
follow from 0 K atomisation energies this repository already carries (`CH4` 1642.0, `C₂H₂` 1625.7,
`H2` 432.1 kJ/mol).

**So the carbon store splits, and only the small part needs nucleation:**

| carbon store, 43.3% of atomisation | share | how it returns |
| --- | ---: | --- |
| recovered by reaching acetylene | **36.3 points (84%)** | gas-phase chemistry, species already in `eos_methane` |
| requires actual condensation | **7.0 points (16%)** | nucleation — N10 item 5 proper |

That is **5.2% of the chamber's energy budget** genuinely hostage to soot, not the 32.4% implied
by treating the whole carbon store as condensation-limited. **And the two-phase lag exposure
shrinks with it**: carbon that leaves as C₂H₂ is gas, so it stays momentum-coupled and there are
no particles to lag.

**The equilibrium path is C → C₃ → C₂H₂, and it is temperature-gated.** Solved on `eos_methane`
at 1 kg/m³, the fraction of carbon nuclei by species:

| `T` [K] | atoms | C₃ | C₂H₂ | store still held |
| ---: | ---: | ---: | ---: | ---: |
| 6000 | 0.16 | 0.67 | 0.03 | 0.599 |
| 5000 | 0.03 | 0.67 | 0.16 | 0.397 |
| 4000 | 0.00 | 0.30 | **0.63** | 0.206 |
| 3500 | 0.00 | 0.08 | **0.89** | 0.144 |
| 3000 | 0.00 | 0.01 | **0.98** | **0.118** |

`held = 0.118` is the 89% return, arrived at independently. **C₃ is the intermediate**, which is
worth stating because it is also this EOS's largest admitted weakness (weakness 1) — the species
that carries the carbon between 4000 and 6000 K is the one whose partition function is least
trustworthy.

**And the flown nozzle does not get there — it parks the carbon in C₃.** At the exit plane:

| throat | exit `T` [K] | exit `rho` | free atoms | **C₃** | C₂H₂ | store held |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 7.0 m² | 5584 | 0.207 | 0.250 | **0.658** | 0.009 | 0.698 |
| 4.0 m² | 5067 | 0.110 | 0.160 | **0.750** | 0.015 | 0.650 |
| 2.0 m² | 4561 | 0.052 | 0.090 | **0.819** | 0.024 | 0.599 |

That is better than free atoms — C₃ is strongly bound and has already banked about 62% per carbon
of what condensation would give (439.9 against 711.2 kJ per mole of C) — but it is a long way from
done. **The gap that matters is not the soot gap.** From the exit state to acetylene equilibrium
is **59 points of atomisation**; from acetylene to fully condensed carbon is **7**. The dominant
missing energy is simply unfinished gas-phase recombination, and it is unfinished because the
nozzle stopped too hot, not because any nucleation failed.

**This is the same finding as W6 arriving from the carbon side: the binding constraint is
expansion ratio.** The acetylene energy is real, large, gas-phase and needs no nucleation
theory — and this geometry leaves essentially all of it on the table because it stops 1000–1500 K
too hot. Reaching ~3500 K collects it; reaching it *at density* (≥0.2 kg/m³) collects it before
the chemistry can freeze, which W5's deep-expansion diagnostic puts at 3300–3500 K.

**What this does not settle.** W5's Damköhler analysis is run on `H + H + M` and therefore prices
the H₂ channel only. **The acetylene path has its own kinetics and they are not tested here** — but
they are a far better prospect than the three-body carbon association W7 could find no data for,
because C₂H₂ forms through *bimolecular* hydrocarbon reactions (`C₂H + H₂`, `C₂H₂ + H` and
relatives) that the combustion literature characterises heavily. **That is the next thing to
build**, and it is ordinary chemical kinetics rather than nucleation theory.

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

- **N10 item 5 (carbon nucleation).** Not started, and W7 now shows it *cannot* be started as
  a rate problem: the evaluated literature has one `C + C + M` measurement and no `C + H + M`
  at all. It needs classical-nucleation machinery that does not exist here. It is 43% of the
  store, and it is the last thing standing between this study and a complete methane answer.
- **An ammonia specific impulse (N10 items 2–3 for `NH3`).** The `N + N + M` rates are in
  hand, but the rung needs an `eos_ammonia` alongside `eos_methane`, and W7's 1.3-decade
  spread on the nitrogen channel would dominate the answer. **What is needed to unblock it:**
  a resolution of Byron 1966 against the 2025 ab initio result, or an ammonia number reported
  as a band.
- **N9 items 1–7.** Not started. W4 confirms they are still live: the sealed vessel does not
  dissolve the arrival transient.
- **N11.** Not started.

## The rate literature this rests on

**Stated in full because the answer is only as good as these, and because an agent in the paper
repository cannot run the code that consumes them.** Every coefficient below was read from the
NIST Chemical Kinetics Database (`kinetics.nist.gov`), which reproduces each evaluation's own
fit rather than refitting it; the squib is NIST's record identifier, so any line can be traced
back. Nothing here was computed on this side. The form throughout is

    k(T) = A (T / 298 K)^n   in cm^6 molecule^-2 s^-1;   1 cm^6 = 1e-12 m^6.

They live in `python/puffsat/walled_nozzle/rates.py`, one named constant per row.

### `H + H + M -> H2 + M` — methane's recoverable store

| M | `A` | `n` | range [K] | unc. | source | squib |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| H₂ | 8.85e-33 | −0.60 | 100–5000 | ×3.16 | Baulch et al. 1992, *J. Phys. Chem. Ref. Data* **21**, 411 | `1992BAU/COB411-429:216` |
| H₂ | 9.04e-33 | −0.60 | 50–5000 | ×2.51 | Cohen & Westberg 1983, *JPCRD* **12**, 531 | `1983COH/WES531:42` |
| Ar | 6.48e-33 | −1.00 | 77–5000 | ×2.00 | Cohen & Westberg 1983 | `1983COH/WES531:44` |
| N₂ | 9.1e-33 | −1.30 | 77–2000 | ×2.00 | Tsang & Hampson 1986, *JPCRD* **15**, 1087 | `1986TSA/HAM1087:253` |
| **H₂** | **4.83e-33** | **0** | **2500–7000** | — | **Hurle, Jones & Rosenfeld 1969**, Proc. 8th Int. Shock Tube Symp. 253 | `1969HUR/JON253-276:1` |
| Ar | 1.68e-32 | 0 | 2500–7000 | — | Hurle et al. 1969 | `1969HUR/JON253-276:2` |

Two points about this table matter more than any single row.

**The channel is genuinely known.** Baulch and Cohen & Westberg evaluated `M = H₂`
*independently* and landed within 2% of each other with the same exponent. That agreement is a
stronger statement than either stated uncertainty bar, and it is what a regression test in
`test_walled_nozzle_freeze.py` now pins.

**Most of the nozzle is measured, and the extrapolation that remains is conservative.** Hurle et
al. measured 2500–7000 K directly, which covers the expansion from just below the throat all the
way to the exit; only the throat itself, near 8500 K, sits above their ceiling, and the atomic-H
efficiency is read from a 2900–4700 K experiment. So this is *not* an extrapolation-free result —
see weakness 6. What makes it safe is the **direction**: Hurle found the rate
*temperature-independent* across their window, so the steep negative exponents in the cold
evaluations are low-temperature behaviour that does not continue upward. Extrapolating Baulch's
`T^−0.6` to 5000 K gives 1.6e-45 m⁶/s where Hurle measures 4.8e-45 — **using the cold fits above
their window understates the rate by about 3×**. Every margin quoted in W5 is therefore a floor,
not an estimate.

### Third-body efficiency — the find that changes the answer

Three independent shock tubes each measured `H + H + M` with M = Ar, H₂ **and atomic H** over
2800–5330 K, all fitted with the same `n = −1`. Ratios taken *within* one experiment cancel that
experiment's systematics:

| study | `k_H / k_Ar` | `k_H2 / k_Ar` | squibs |
| --- | ---: | ---: | --- |
| Jacobs, Giedt & Cohen 1967, *J. Chem. Phys.* **47**, 54 | 20.0 | 2.5 | `1967JAC/GIE54-57:1,2,3` |
| Rink 1962, *J. Chem. Phys.* **36**, 262 | 6.7 | 2.0 | `1962RIN262-265:1,2` |
| Patch 1962, *J. Chem. Phys.* **36**, 1919 | 66.6 | 10.0 | `1962PAT1919-1924:1,2,3` |

**Atomic hydrogen is 7–67× better at stabilising the collision than argon, and 2–10× better than
H₂**, and all three studies agree on the ordering. This is not a detail. W1 established that at
the flown 10 kK the charge is 93–98% dissociated — so the third body in this nozzle *is* atomic
hydrogen. A margin computed on H₂, or on a generic collider, understates the recombination rate
by roughly an order of magnitude. The paper-side hand estimates of ×3400 and ×403 do not say
which third body they assumed; if it was H₂ or an inert, they are pessimistic.

### `N + N + M -> N2 + M` — the ammonia rung only

| M | `A` | `n` | range [K] | source | squib |
| --- | ---: | ---: | ---: | --- | --- |
| N₂ | 4.16e-33 | −0.50 | 6000–9000 | Byron 1966, *J. Chem. Phys.* **44**, 1378 | `1966BYR1378-1388:4` |
| N | 1.29e-30 | −1.50 | 6000–9000 | Byron 1966 | `1966BYR1378-1388:5` |
| Ar | 1.60e-33 | −0.50 | 6000–9000 | Byron 1966 | `1966BYR1378-1388:6` |

Byron's window *straddles* the chamber temperature, so nitrogen needs no upward extrapolation at
all. Against it we set a modern ab initio result:

> Notey, Jo & Panesi, "Master equation study of three-body recombination of nitrogen and oxygen
> in non-equilibrium hypersonic flows", *J. Chem. Phys.* **163**, 194311 (2025); arXiv:2506.17452.

That paper solves a state-to-state master equation on ab initio potential energy surfaces for
**atoms held at 10 000 K and then plunged into a cold bath** — which is exactly what a nozzle
does to a chamber. Its quasi-steady `N + N + N` values are 9.125e-46 m⁶/s at 2500 K and
7.60e-46 at 5000 K. **It sits 1.3 decades below Byron.** That disagreement is the honest width of
the nitrogen channel, and it is why the ammonia rung is reported below as a band rather than a
number.

### Carbon — searched, and the answer is that there is nothing to pick

| reaction | evaluated records | best available |
| --- | ---: | --- |
| `C + C + M -> C2 + M` | 1 | Slack 1976, *J. Chem. Phys.* **64**, 228: 5.46e-31 (T/298)^−1.60, M = Ar, 5000–6000 K, **±55%** |
| `C + H + M -> CH + M` | **0** | nothing |

Hydrogen's channel has six independent evaluations and three shock-tube studies spanning
2500–7000 K. Carbon's has one measurement, one third body, a 1000 K window and a 55% error bar —
and the `C + H + M` channel is absent from the evaluated literature entirely.

**That asymmetry is a result, not a failed search.** An association channel gets measured when it
controls something measurable, and the carbon ones do not, because carbon in a cooling
carbon-rich gas does not go to diatomics — it goes to clusters, and then to soot. So the 43% of
methane's atomisation locked in carbon **cannot be settled by a rate coefficient at all**, the
way the 52% held as H₂ can. This is independent support for N10 item 5 being a nucleation problem
needing different machinery, rather than one more row in the table above.

### One citation checked and rejected

Bourdon & Vervisch, *Phys. Rev. E* **54**, 1888 (1996), "Three-body recombination rate of atomic
nitrogen in low-pressure plasma flows", is the obvious hit for `N + N + M` and **is not that
reaction**. It is electron–ion collisional-radiative recombination (`N⁺ + e + e`), i.e. the
ionisation store, which `recombination.py` already carries from Zel'dovich & Raizer. Recorded
here so the next reader does not spend the search twice.

## Provenance

| number | produced by | command |
| --- | --- | --- |
| chamber composition, `k`, `u`, pressures, Isp ratios, item 0 timescales | `puffsat.walled_nozzle.chamber` | `make walled-nozzle-chamber` |
| Project 242 bracket, pure-hydrogen rung | `puffsat.walled_nozzle.hydrogen` | `make walled-nozzle-hydrogen` |
| freeze stations, Damköhler margins, throat sensitivity, exit Isp, the `N+N+M` comparison | `puffsat.walled_nozzle.freeze` | `make walled-nozzle-freeze` |
| the rate coefficients themselves, one named constant per published value | `puffsat.walled_nozzle.rates` | `make walled-nozzle-test` |
| the equilibrium EOS underneath all of them | `puffsat.eos_methane`, `puffsat.walled_nozzle.hydrogen` | `make walled-nozzle-test` |

Committed artifacts: `data/results/walled_nozzle/chamber.csv` and `data/results/walled_nozzle/freeze.csv`.

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
   physical answer. That is N10 item 5 and it is not answered here. **This bites W5's 2 m² row
   hardest**, whose exit sits at 4349–4561 K, and the W6 recommendation to narrow the throat walks
   further into it — the +155 s is real for the hydrogen store and silent about the carbon one.
5. **The nozzle's clock rests on an assumed shape.** Temperature against area ratio is solved, but
   turning it into a time needs the area profile along the bore, which ADR-0016 does not give; a
   linear opening over the chamber's own length is assumed, exactly as `expansion.cooling_history`
   does for the magnetic nozzle. `Da` scales linearly in that length, so a nozzle half as long
   halves every margin in W5 — the flown case survives it, the 673 m³ / 2 m² corner does not.
6. **The `H + H + M` rate is extrapolated in temperature at the hot end**, though conservatively.
   The evaluations stop at 5000 K and the throat starts near 8000 K. Hurle et al.'s direct
   2500–7000 K measurement is *faster* than the extrapolated cold fits by about 3×, so the
   extrapolation understates the recombination rate and the margins in W5 are floors, not
   estimates. Being wrong in the safe direction is still being wrong.
7. **Third-body efficiencies are read across from Ar/H₂/H shock-tube ratios.** The hydrocarbon and
   carbon-bearing colliders are charged at argon's efficiency, the slowest of the three measured,
   because nothing better exists. At the flown state they are a small minority of the third bodies
   (the gas is 93–98% dissociated) so this is a minor exposure — but it grows at the cold end,
   which is where the freeze question lives.
