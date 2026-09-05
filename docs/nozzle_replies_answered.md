# The replies R1–R15, answered

Answers to `katzseth22202/Balloon-Pulse-Propulsion` @ `b6dba6f`,
`docs/nozzle_replies_to_impact_sim.md` (raised 2026-09-04, replying to this repository's
`docs/nozzle_asks_answered.md` @ `49a5e49`). Worked in
`katzseth22202/puffsat_impact_simulation` on 2026-09-05.

**This document is written to be copied into the paper repository and worked there.** Section
["What the paper should change"](#what-the-paper-should-change) is self-contained: every number
needed to make every edit is stated in it, because an agent in the paper repository cannot run the
code that produced them. Items continue the **P-sequence** (P11 onward) so the paper keeps one
namespace across both batches. Every figure here came from a run on 2026-09-05, not from recall;
the [provenance table](#provenance) says which.

**The headline is a concession.** R1 is right that `jet.py` modelled the wrong regime, and P2's
`eta_geom` = 0.70–0.88 is withdrawn. What replaces it is **not** R1's 0.48–0.64 either. The
corrected number is **0.58–0.86 as flown**, and the reason both previous estimates missed is the
same: neither had the divergence traced against a real field solve.

---

## Status

| # | verdict | does the paper change? |
| --- | --- | --- |
| **R1** | **Conceded, and it retires `jet.py`.** `Kn` = 2.5e-7 at the most collisionless station anywhere; a parcel collides ~1e6 times in transit. But the replacement number is **higher** than R1's, not lower. | **yes** (P11, P12) |
| **R2** | **Answered on the real Biot–Savart solve.** The paraxial probe is good to ~4 points and errs *generous*, not harsh. Inside the winding there is no divergence at all. | **yes** (P16) |
| **R3** | Recorded, nothing owed. The conclusion survives the framework change. | no |
| **R4** | **Premise withdrawn.** Its factor of 2–3 is `mu`-framework arithmetic that R1's own amendment retires. Still worth ~3% on `f_d`. | **yes** (P12) |
| **R5** | **Answered geometrically.** 17.0% / 14.0% / 0% against R5's 12.9% / 11.6% / 0%. The wall-interaction half stays deferred. | **yes** (P16) |
| **R6** | Conceded. P5 is a required edit at two sites, not a clarification. | **yes**, paper-side already done |
| **R7** | **Done**, and with a counter-suggestion the paper should take. | **yes** (P22, small) |
| **R8** | **Declined, with the physical reason.** Nothing is needed for detachment; a wall cannot act where the pointing loss happens; and R8's own thermal gate puts the bell past the station where the field lines have already turned. | **yes** (P14) |
| **R9** | **Answered, and both of the paper's front numbers are confirmed from our EOS.** The cone construction is validated. | **yes** (P18) |
| **R10** | **P8's constraint withdrawn.** The loss cone is a collisionless object and does not apply; the correct criterion is choking, and it binds at the chamber. | **yes** (P17) |
| **R11** | **Answered, and it clears everywhere above ~56 km/s.** Leg 2 (75 km/s) clears with a 10.9 m extension and on one branch with none. Leg 1's cold tail cannot clear at any length. | **yes** (P14) |
| **R12** | **Answered — and `phi` is the wrong quantity to multiply by.** 0.99 flat on the frozen branch, 0.04–0.85 on the equilibrium one, and anti-correlated with the total radiated. The passive-structure claim does not fail. | **yes** (P21) |
| **R13** | **Confirmed as the binding item**, quantitatively, and then **deferred with a cost.** | **yes** (P15) |
| **R14** | **Do not adopt 672.9 m³.** Our pair was never derived as a pair, and there is a four-way consistent solution neither repo found. | **yes** (P20) |
| **R15** | **The liner**, on both questions — and **yes**, P9's 5 T is a bore-referenced standoff number, so the larger version is back on the table. | **yes** (P19) |

**Fifteen items fall out, P11–P21 and P23–P26**, plus one small one (P22). P23–P26 were raised
by Seth mid-batch and are not replies to any R-item.

**If you read only three:** **P12** (the corrected `eta_geom`, which is the number R1 asked for and
is better than R1's), **P15** (where the chain actually fails, which is a much smaller place than
R13 gets built), and **P19** (the standoff wall, which unblocks ADR-0013 and is a one-line answer).

**Two of the replies rest on arithmetic their own amendments retire** — R4 entirely and R10's
threshold — and those are in [What the replies got wrong](#what-the-replies-got-wrong).

**Two things we need back before the last question can be closed**, both owned by
`aim_is_all_you_need` rather than by either of us:

1. **What fraction of leg 1's pulses land below ~56 km/s** (the closing-speed schedule). That is
   the difference between the chemistry shortfall being a footnote and being the design driver.
   See P15.
2. **Which figure of merit the growth tables maximise**, and whether `k` may vary per pulse.
   Momentum per projectile and momentum per kg expended want `k` moved in opposite directions,
   and P24 finds the two legs pushed apart by a drift term that changes sign between them.
   See P23 and P24.
3. **Which efficiency multiplies which term** in `sqrt(1+k) x eta_jet ± 1`. Our `eta_geom` is an
   expansion-frame quantity and the `sqrt(1+k)` bound includes the drift's energy; the drift is
   turned by a magnetic wall rather than by the flare. See P24.

---

## What the paper should change

### P11. R1 is right: the expansion is collisional, and `jet.py` modelled the wrong regime

**This is a withdrawal of our own P2. Make it before anything else in this document.**

`jet.py`'s conversion `alpha_exit = 1 - (1 - alpha_0)/(A/A*)` is a guiding-centre argument. It
needs two things, and its own docstring listed both under "what this is not":

1. **Magnetization** — the field varies slowly over a gyroradius. `adiabaticity_parameter` tested
   this and it passes.
2. **Collisionlessness** — a particle keeps its own `mu` long enough for the field to act on it.
   **Nothing tested this**, and it fails.

Tested now, on our own solved cooling history rather than on a bag average, and on the **longest**
of the three heavy-particle collision channels — neutral hard-sphere, which is the reading most
favourable to keeping `mu`:

| leg | branch | longest mean free path | `Kn` vs the 6 m bore | collisions per parcel in transit |
| --- | --- | ---: | ---: | ---: |
| 45.58 | equilibrium | 1.43 µm | 1.19e-7 | 8.0e6 |
| 45.58 | frozen | **2.96 µm** | **2.46e-7** | 2.0e6 |
| 56.53 | equilibrium | 1.50 µm | 1.25e-7 | 9.7e6 |
| 75.00 | equilibrium | 1.96 µm | 1.63e-7 | 7.4e6 |
| 75.00 | frozen | 1.23 µm | 1.02e-7 | 6.1e6 |

The most collisionless station anywhere is `Kn` = **2.5e-7**, against the `Kn` ~ 0.1 a
guiding-centre treatment would want — short by a factor of 4e5. **A parcel collides between two
and ten million times crossing the nozzle**, and `mu` does not survive the first one.

**So P2's `eta_geom` = 0.70–0.88 is withdrawn**, and so is `jet.py`'s docstring claim that "the
deeper the flare the more complete the conversion". R1's replacement framing is right: at
`beta` = 0.013–0.073 the field exceeds plasma pressure 15–75× and *is* the wall, so this is a de
Laval nozzle with magnetic walls. That is the paper's own "walled by field rather than fenced by
it", and it needs no hardware change.

**Reproduce:** `make analysis-continuum`.

---

### P12. The corrected `eta_geom` is 0.58–0.86 as flown — above P2's floor and above R1's estimate

**This is the number R1 asks for, and it is the item that moves the paper's chain.**

In a collisional exhaust a parcel has a directed speed `u` along its streamline plus isotropic
thermal motion about it. Thermal motion averages to zero in the mean and contributes fully to the
mean square, so

    eta_geom = <cos theta> / sqrt(1 + <v_th^2>/u^2),    <v_th^2> = 3 k T / m_bar

which is exactly R1's `<cos theta>/sqrt(1 + 3/(gamma M^2))` with `R_sp T = c_s^2/gamma` substituted
out. **We evaluate it without `gamma`**, from the mean *heavy-particle* mass out of the equilibrium
composition, because the mixture's effective `gamma` varies through the expansion. (Verified
identical to R1's form on a `gamma`-law gas, to machine precision, by test.)

**Two terms, and the paper needs to know which one is doing the work.**

**The thermal term is not the problem.** At the flown exit it is **0.86 to 0.97** across both legs
and both branches — already above the 0.775 the growth tables sweep to. The plume is well
converted by the time it leaves; there is very little thermal store left to recover.

**The divergence term is the whole of it, and it is entirely a downstream quantity.** Traced
through the real Biot–Savart field of the current winding rather than a paraxial expansion,
`<cos theta>` **inside** the magnet is 0.999 — a solenoid bore is axial and nothing is lost in it.
The loss is all past the last coil, where the vacuum field opens toward its own return path:

| distance past the exit plane | `<cos theta>`, capped 12 T | `<cos theta>`, flown 20 T | tubes still running forward |
| ---: | ---: | ---: | ---: |
| 0 (exit plane) | 0.9985 | 0.9968 | 100% |
| 0.5 exit radii | 0.9960 | 0.9889 | 100% |
| 1.0 exit radii | 0.9543 | 0.9132 | 100% / 75% |
| **1.44 exit radii (P3 low)** | **0.8166** | **0.7372** | 62% / 50% |
| **2.00 exit radii (P3 high)** | **0.6890** | **0.5682** | 50% |

R1's paraxial probe gives 0.703 and 0.632 at those two stations. **The real solve is more generous
at the near edge and close at the far edge**, so R1's own note that "0.48–0.64 is more likely
generous than harsh" is the wrong way round.

**The last column is a result in its own right.** By 1.4 exit radii, 38–50% of the traced tubes
have *turned* — the field line curves back toward the return flux. A plume cannot follow a field
line back, so **detachment has to happen inside about 1.5 exit radii or it does not happen at
all**. That is an independent constraint the paper does not have, and it is the thing a staged
nozzle (R8) or an extension (R11) is really buying.

**Put together, at the flown geometry, per leg and branch:**

| leg | branch | profile | `A/A*` | `M` | `M_A` | detach | `<cos theta>` | thermal | **`eta_geom`** | `eta_jet` (water) |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 45.58 | equilibrium | flown 20 T | 4.00 | 3.233 | 0.353 | 2.00 | 0.572 | 0.909 | **0.520** | 0.380 |
| 45.58 | frozen | flown 20 T | 4.00 | 3.440 | 0.355 | 2.00 | 0.570 | 0.972 | **0.554** | 0.405 |
| 75.00 | equilibrium | flown 20 T | 4.00 | 2.726 | 0.576 | 1.45 | 0.735 | 0.892 | **0.656** | 0.597 |
| 75.00 | frozen | flown 20 T | 4.00 | 3.441 | 0.578 | 1.44 | 0.737 | 0.956 | **0.704** | 0.641 |
| 45.58 | equilibrium | **capped 12 T** | 2.40 | 2.421 | 0.439 | 1.73 | 0.669 | 0.869 | **0.582** | 0.425 |
| 45.58 | frozen | **capped 12 T** | 2.40 | 2.669 | 0.444 | 1.72 | 0.656 | 0.914 | **0.599** | 0.438 |
| 75.00 | equilibrium | **capped 12 T** | 2.40 | 2.291 | 0.708 | 1.26 | 0.922 | 0.863 | **0.795** | 0.724 |
| 75.00 | frozen | **capped 12 T** | 2.40 | 2.669 | 0.723 | 1.24 | 0.925 | 0.928 | **0.859** | 0.782 |

**`eta_geom` = 0.58–0.86 on the winding the paper now flies.** Every case clears
`sec:mass_interest`'s forward-thrust floor of `1/sqrt(1+k)` = 0.324 by a wide margin, and the hot
leg reaches the swept 0.775 without any extension at all.

**What to write.** P2's mechanism claim was right for the wrong reason and its number is
withdrawn. The correct statement is that the plume leaves the magnet **well converted and badly
aimed**: the thermal term is 0.86–0.97 and the divergence term is 0.57–0.93, so the chain's
remaining nozzle loss is a *pointing* loss, not a conversion loss. R1's structural criticism
stands and its arithmetic does not.

**Reproduce:** `make analysis-nozzle-extension`, `make analysis-nozzle-fluxtube`.

---

### P13. ADR-0012's 12 T cap is a gain, not a compromise — and its own table understates it by 30×

**The paper is currently pricing this backwards, in its own favour but for the wrong reason.**

`A/A*` is the field ratio by flux conservation. Capping the peak at 12 T while P9 holds the exit
at 5 T cuts the magnet's own area ratio **from 4.00 to 2.40** — a 40% cut in the expansion. The
obvious expectation is that it costs `eta_geom`, and ADR-0012's table says it costs less than
0.005.

**It costs the thermal term and buys back more on the divergence term.** From the table in P12:

| | flown 20 T | capped 12 T | change |
| --- | ---: | ---: | ---: |
| exit `M` (45.58, equilibrium) | 3.233 | 2.421 | worse |
| thermal term | 0.909 | 0.869 | **−0.040** |
| exit `M_A` | 0.353 | 0.439 | better |
| detachment station | 2.00 radii | 1.73 radii | better |
| `<cos theta>` | 0.572 | 0.669 | **+0.097** |
| **`eta_geom`** | **0.520** | **0.582** | **+0.062** |

The mechanism: a smaller area ratio leaves the plume **denser** at the exit at the same 5 T, so
`v_A = B/sqrt(mu0 rho)` is lower and `M_A` is *higher*. A higher `M_A` detaches sooner, and sooner
means less fanning.

**Across all four flown cases `eta_geom` moves by +0.06 to +0.15, in the paper's favour**, against
ADR-0012's stated "less than 0.005". The two frameworks disagree because in the `mu` picture
`alpha` depends on the field *ratio*, which the cap preserves at the exit; in the continuum picture
the cap moves the detachment surface.

**What to write.** ADR-0012's decision is right and its justification is incomplete. The cap is not
a trade of `eta_geom` for field energy — it improves both. Its `eta_geom` column should be
restated, and the sentence "why `eta_geom` does not pay for this" becomes "why `eta_geom` is paid
*by* this".

**Reproduce:** `make analysis-nozzle-extension`.

---

### P14. R11's extension works, the short one is enough, and it clears on one leg only

**R11 asks four things. All four are answered, and the fourth dissolves.**

#### 1. Where does the `M_A` = 1 crossing actually land?

Run on our solver rather than the 1-D isentropic area–Mach relation, on the capped profile, both
legs, both branches. `M_A` at the end of each prescribed extension:

| extension | add-on | total `L` | `A/A*` | `B_exit` | `M_A` at its exit | crossing inside? |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| none | — | 23.8 m | 2.40 | 5.00 T | 0.44–0.72 | no |
| **hot-favourable** | 10.9 m | 34.6 m | 6.31 | 1.90 T | **0.75–1.24** | **leg 2 yes, leg 1's cold end no** |
| hot-pessimistic | 12.3 m | 36.0 m | 6.96 | 1.72 T | 0.80–1.31 | hot yes, cold no |
| cold-favourable | 23.9 m | 47.5 m | 13.44 | 0.89 T | 1.12–1.88 | **both** |
| cold-pessimistic | 30.3 m | 53.8 m | 17.88 | 0.67 T | 1.30–2.19 | both |

**R11's central claim holds: an extension does move the crossing inside the flare.** The lengths
are close to R11's own — its 10.9 m case reaches `M_A` = 1.23–1.24 at 75 km/s, so the crossing
lands comfortably inside. Its cold cases are needed only if the crossing must land inside at
45.58 km/s, which needs about 24 m.

#### 2. Is the collision front contained past the bag?

**The question dissolves, and R11 should drop it rather than answer it.** R11 applies P9's
`r_front/r_bore = beta^{1/(2 gamma)}` past `z` = 23.8 m and gets `beta` = 1.38 at the extension
exit. But P9's formula is about a *collision* front standing against a field, and **there is no
collision front past the bag**: the snowplow terminates where the mist ends. The plow runs from
the chamber to the exit, arrives having merged with all 213 kg, and what leaves the magnet is the
merged plume, not a front.

What the extension actually sees is the plume, whose `beta` we solve directly: **0.039–0.268 at
the 10.9 m extension's exit**, two orders below the 1.38 R11 feared. No containment problem
exists there, and the ~1 m of extra winding clearance R11 books for it can come out of the
conductor column.

#### 3. `eta_geom` through the extension

| leg | branch | none | 10.9 m | 30.3 m |
| --- | --- | ---: | ---: | ---: |
| 45.58 | equilibrium | 0.582 | 0.862 | 0.930 |
| 45.58 | frozen | 0.599 | 0.917 | 0.976 |
| 75.00 | equilibrium | 0.795 | 0.894 | 0.918 |
| 75.00 | frozen | 0.859 | 0.960 | 0.978 |

Against R1's 0.48–0.64. **Most of the gain is bought by the short extension**: going from 10.9 m
to 30.3 m adds 0.02–0.07 while nearly tripling the add-on. R11's cold-pessimistic 30.3 m case is
not worth building.

#### 4. Narrowing `beta` does not narrow the design, and does not need to

R11 hopes that narrowing `beta` will narrow the 3× spread between an 11 m add-on and a 30 m one.
**It runs the other way.** `beta = p/(B^2/2 mu0)`, and as the flare opens `B` falls faster than
`p`, so the spread *widens* with extension length: 4.6× with no extension, 6.9× at 10.9 m, 10.7×
at 30.3 m.

**But the decision does not turn on it**, which is the useful answer. See the verdict below.

#### The verdict, and it is a split one

| leg | branch | no extension | 10.9 m extension | vs the swept 0.775 |
| --- | --- | ---: | ---: | --- |
| **leg 1, cold end** (45.58) | equilibrium | 0.425 | 0.630 | **misses at every length** |
| **leg 1, cold end** (45.58) | frozen | 0.438 | 0.671 | **misses at every length** |
| **leg 2** (75.00) | equilibrium | 0.724 | **0.813** | clears with the short extension |
| **leg 2** (75.00) | frozen | **0.782** | **0.873** | clears with *no* extension |

**Leg 2 clears, and the short extension is enough** — 10.9 m of add-on at under 2 T, which
is cheap per metre. On the frozen branch it clears without any extension at all, which is new.

**Leg 1's cold end misses at every length**, and no nozzle can fix it: `eta_chem` = 0.731 caps it
below 0.775 before geometry enters. See **P15**, which locates exactly where that starts — and it
is *not* the whole leg.

**What to write.** `sec:jet_efficiency`'s original claim that detachment happens inside the nozzle
is **restorable above about 56 km/s with a 10.9 m magnetic extension**. P3's correction becomes a
design change rather than a retraction, exactly as R11 hoped. The 30.3 m case should be dropped.

#### And this retires R8. **Declined, with the physical reason.**

R11 predicted it would, and it does — for R11's reason and for two stronger ones it did not have.

**First, separate the two jobs, because R8 conflates them.**

- **Detachment needs no structure at all.** Free expansion crosses `M_A` = 1 on its own at
  1.24–1.73 exit radii. That is P3, unchanged. Nothing — bell, flare, or otherwise — is *required*
  for the plume to let go. What P3 killed was the claim that it lets go *inside* the nozzle, which
  is a wording fix and not a hardware problem.
- **Pointing is the whole job**, and **it finishes before detachment.** `<cos theta>` is 0.999
  inside the magnet, 0.997 at the exit, and 0.57–0.93 by the time `M_A` reaches 1. All of the
  damage is done in the sub-Alfvénic stretch.

**Second, a wall cannot act in that stretch, for two independent reasons.**

1. **The plume is glued to field lines there.** `beta` = 0.033–0.154 means the field is 6–30×
   over-strength. To be steered by a wall, plasma has to *cross* field lines, which at that `beta`
   it cannot. This is R8's own stated difficulty, now quantified.
2. **The wall would sit outside the plume's magnetic envelope entirely.** Traced on R8's own
   geometry, the plume's **outermost flux tube reaches 6.34 m and turns** — 5% past the 6.03 m
   exit — against R8's *smallest* proposed bell radius of **7.0 m**. The bell is not a surface the
   plume slides along; it is a surface the plume would have to break out through.

**Third, R8's own thermal gate pushes the bell the wrong way.** Its table has graphite *not*
surviving at the magnet exit (16 224 K on the equilibrium branch) and surviving only from 7.0 m
outward, once the gas has expanded and cooled. So the wall is *required* to sit past exactly the
station where the field lines have already turned and the fanning has already happened. **The gate
that lets the wall survive is the gate that makes it useless.**

**The squeeze.** A bell is either close in — where it must fight a 6–30× over-strength field, and
where its own thermal table forbids it — or far out, where the plume is free but has already
fanned to 0.62–0.74. Both ends fail, and the ablation work in `detachment.py` establishes only
that a wall would *survive*, which turns out not to be the question.

**The flare wins by working with the coupling instead of against it.** The plume is glued to field
lines, so bend the field lines. A magnetic extension holds it inside a controlled 15° cone to
`M_A` = 1 and releases it at `<cos theta>` = **0.983**, with no interface, no wall erosion over
eleven cycles, and no MHD question outstanding.

**A correction we owe on our own side.** `detachment.py`'s ablation docstring said "the gas
dynamics were never the obstacle: every exit is supersonic, and a diverging wall on supersonic
flow is exactly what a de Laval nozzle is. The obstacle is heat." **That was the wrong Mach
number**, and it is the sentence that made the bell look attractive to both repositories. Being
supersonic (`M` = 2.4–3.4) is not being super-Alfvénic (`M_A` = 0.44–0.72). Corrected in place.

**What we are and are not claiming.** This is geometry plus our own `beta` and `M_A`, not an MHD
solve — we have not proved a bell impossible. But R8 asked what would settle whether the paper can
carry a staged nozzle, and the answer is that it should not try: **the structure is not needed for
detachment and cannot act where the pointing loss occurs.** R8 comes off the critical path.

#### What the paper should say, and it should say it

**Seth's call, 2026-09-05: this belongs in the paper**, not only in a reply document. A physical
bell is the obvious thing for a reader to propose — a magnetic nozzle that sprays looks like it
wants a wall — and the paper should foreclose it explicitly rather than leave the reader to
re-derive R8. Suggested substance, for `sec:jet_efficiency` or `sec:two_leg_nozzle`:

> A physical bell downstream of the winding is the natural remedy for the divergence loss, and it
> does not work here. It is not needed for detachment: the plume crosses `M_A` = 1 on its own
> 1.2 to 1.7 exit radii past the last coil. And it cannot supply the exhaust direction it would be
> built for, because the loss happens *before* that crossing, while the plume is still
> sub-Alfvénic and tied to the field. At `beta` = 0.03 to 0.15 the field exceeds plasma pressure
> by six to thirty times, so plasma cannot cross field lines to reach a wall; the plume's own
> bounding flux tube turns back at 6.3 m, inside the 7.0 m a graphite bell would have to start at
> to survive the exit temperature. The condition that lets the wall survive is the condition that
> places it beyond the region where it could act. The remedy is to bend the field lines rather
> than to obstruct them, which is what the magnetic extension does.

**Three points to keep if the wording changes**, because each closes a different escape:

1. **Detachment needs nothing.** Otherwise a reader assumes the bell is load-bearing.
2. **The loss is upstream of detachment.** Otherwise a reader assumes a bell placed after
   detachment could still fix the aim — it cannot, the aim is already set.
3. **The thermal gate and the geometry gate point the same way.** Otherwise a reader proposes a
   more refractory wall material, which does not help: the problem is where the wall must sit, not
   what it is made of.

**Reproduce:** `make analysis-nozzle-extension`.

---

### P15. R13 is right — and the failure is confined to below ~56 km/s, not to "the cold leg"

**This is the item that decides what to build next, and the second half of it narrows the problem
a long way.**

#### The proof R13 wanted

The chain is `eta_jet = eta_geom x eta_chem`, and `eq:eta_chem` charges 45.58 km/s **0.731**. So
even a perfect nozzle gives

    eta_jet(45.58 km/s, water)  <=  1.000 x 0.731  =  0.731  <  0.775

**That misses the paper's own swept target before any nozzle exists.** Our extension sweep reaches
`eta_geom` = 0.976 there — within 2.4% of perfect — and still lands at `eta_jet` = 0.713. R13's
claim was an argument; it is now arithmetic.

#### But "the cold leg" is the wrong unit, and the correct one is much smaller

**Naming first, because it changes the conclusion.** 45.58 km/s and 75 km/s are not two ends of
one range — they are **different legs**, as this repository's own `CONTEXT.md` insists and as the
N-batch's correction 3 already told us:

- **Leg 1**, the overtake / growth push. Closing speed **falls** through it, so it *sweeps* a range.
  45.58 km/s is its cold **end**, and specifically the cold end of the **slowest three-synodic
  cadence** — the extreme corner of the whole envelope, not a typical operating point.
- **Leg 2**, the head-on departure burn at ~75 km/s. A single point, and **never affected**.

Computed across the range, with `k` = 8.52 and water:

| closing speed | `eta_chem` | `eta_jet` at `eta_geom` = 0.93 | verdict |
| ---: | ---: | ---: | --- |
| 45.58 | 0.730 | 0.679 | misses |
| 50 | 0.782 | 0.728 | misses |
| 55 | 0.824 | 0.767 | misses, barely |
| **56.5** | 0.835 | **0.776** | **clears** |
| 65 | 0.878 | 0.816 | clears |
| **70** | 0.896 | **0.833** | **clears** |
| **75** | 0.910 | **0.846** | **clears** |

**The crossover, by how good the nozzle is:**

| nozzle | clears 0.775 above |
| --- | ---: |
| no extension | 71.8 km/s |
| **10.9 m extension, equilibrium branch** | **56.3 km/s** |
| best extension, frozen branch | **51.2 km/s** |
| a hypothetical perfect nozzle | **49.3 km/s** |

**49.3 km/s is the hard floor**: below it, water chemistry alone fails 0.775 whatever the nozzle
does. Between 49 and 56 km/s the nozzle decides. **Above 56 km/s there is no problem at all**, and
leg 2 clears by a wide margin.

**This reframes what the extension buys.** It does not merely raise a number from 0.68 to 0.83 —
**it moves the crossover speed from 71.8 km/s down to about 56**, converting roughly 16 km/s of
the flight envelope from failing to passing. That is a better justification for 10.9 m of magnet
than an efficiency delta.

#### What is now owed, and it is not ours

The remaining question is no longer "does the chain work". It is **"how many of leg 1's pulses
land below ~56 km/s?"** That is the closing-speed **schedule**, which this repository explicitly
does not own (`CONTEXT.md`: "This repository owns neither — both are `aim`'s"). Until `aim` states
the fraction of leg 1's pulses in the cold tail, nobody can say what this costs the mission.

**Ask for that fraction.** If it is small, the whole item is a footnote and argon is optional. If
it is large, argon is the design.

#### What the paper should say

The paper's own sentence is right and should be promoted from an aside to the conclusion of
`sec:jet_efficiency` — *"Both of argon's gains therefore land on the leg where water is
weakest"* — with the qualification that it is **leg 1's cold tail below about 56 km/s**, not
leg 1 as a whole and not leg 2 at all. R1's amended box is likewise right that "the shortfall this
item reports is a chemistry problem, not a nozzle one", and the qualification is the same.

#### What we are not delivering, and what it would take

R13's four sub-asks need an argon EOS — an `eos_argon` beside `eos_water` with an Ar I–III Saha
ladder — and then the whole expansion / detachment / fireball chain re-run on it. That is a new
species, not a new calculation, so it is out of this batch's scope by construction. Costed in
[Deferred, with cost](#deferred-with-cost). **We still recommend it be the next thing built** —
but note that P23 finds argon and the low-`k` lever spend the same currency, and that sub-ask 3
(the radiated share, which goes as `T^4` and multiplies every liner figure in R12) could take the
gain back on either.

**A caution on sub-ask 1.** R13 asks for "`eta_chem` for argon, stated". The recombination
argument in the paper's `CONTEXT.md` is sound as far as it goes — three-body electron-ion
recombination at `n_e^2` with `alpha ~ T_e^-4.5` is nanoseconds at 0.32 kg/m³ against a ~100 µs
expansion — **but it is evaluated at the bag density, and the plume is 13× thinner by the exit**.
Three-body recombination goes as `n_e^2`, so that is a 170× slowdown, and this repository's own
`fireball.py` result on water is that recombination *freezes* below 0.01 kg/m³. Argon will almost
certainly still win that race, but `eta_chem` = 1.000 should not be printed until it has been run.

---

### P16. What misses the winding, on the real solve — and the two repositories are modelling different launches

R2 asks for the flux-tube accounting redone against a real field solve. Done, on the Biot–Savart
field of a 48-coil winding, with tubes traced by RK4 on `dr/dz = B_r/B_z` and truncated at the
winding contour. The tracer's acceptance test is the analytic invariant `B r^2 = const`, which it
holds to under 5% along a tube.

| winding | this repo | R5's paraxial figure |
| --- | ---: | ---: |
| straight 3.5 m, flown 20 T | **17.0%** | 12.9% |
| straight 3.5 m, capped 12 T | **14.0%** | 11.6% |
| **flared, capped (ADR-0011 as amended)** | **0.0%** | 0% |

**R5's conclusion is confirmed and its number is a little generous.** The real tubes fan ~4 points
more than the paraxial estimate, in the direction R5's own status note guessed at but assigned the
wrong sign. At 17.0% the straight winding puts **36 kg/pulse** on the liner against the 4.9 kg
booked, a factor of **7.4**, not 2.3%. The flare zeroes it, which is ADR-0011's whole point and
survives the framework change intact.

**A modelling difference worth naming, because it is worth a factor of three.** We first launched
every tube at the chamber plane — a steady nozzle with one inlet — and got 40–60% missing. That is
wrong, and R5's station-weighted picture is right: **the plume is not fed through an inlet**, it is
a 23.8 m column that already fills the bore when it starts to leave, so a parcel at `z` = 20 m
never sees the chamber field and rides a tube with only `B(20)/B_exit` left to fall. The tables
above use the column launch. The chamber launch is reported alongside in the artifact as the
pessimistic bound.

**But that same weighting does not transfer to R4** — see
[What the replies got wrong](#what-the-replies-got-wrong).

**Reproduce:** `make analysis-nozzle-fluxtube`.

---

### P17. P8's ">= 36 coils" is withdrawn. The loss cone is a collisionless object

R10 disputes P8's pass criterion — existence of a minimum rather than its depth — and proposes
`R > 1/sin^2(theta)` with `sin^2(theta)` from P1's `alpha`. **R10 is right that existence is the
wrong test, and the correct reason retires the criterion rather than sharpening it.**

**Mirror trapping is `mu` conservation seen from the other side.** A particle turns around because
its own `mu` is fixed, and the loss cone is the set of pitch angles for which that never happens.
P11 shows there is no `mu`. The fluid statement is cleaner: `(J x B) . B = 0` identically, so
**the Lorentz force has no component along a field line**. In a collisional plasma the pressure is
a scalar, the parallel momentum equation is `rho Du/Dt = -dp/ds` with no magnetic term, and a
`|B|` variation exerts no force along the flow at all. The mirror force of kinetic theory,
`-(p_perp - p_par) d ln B/ds`, is proportional to an anisotropy that collisions destroy on the
mean-free-path timescale — about 1e-10 s here, against a 2 ms transit.

**So there is no trap, at any coil count.**

**What can still stop the flow, and it is a real constraint with a sharp threshold.** A flux tube's
cross-section is `A ~ 1/|B|`, so a local `|B|` **maximum** is a local area **minimum** — a throat.
Supersonic flow through a contraction decelerates and can **choke**. Isentropic flow at Mach `M`
survives a contraction up to `A/A*(M)`, and the ripple's contraction ratio is exactly the mirror
ratio `field.py` already reports. So the test is **`R < A/A*(M_local)`**.

**The binding station moves from the winding to the chamber**, because `A/A*` goes to 1 at the
sonic point:

| local `M` | margin `A/A*` at `gamma` = 1.25 |
| ---: | ---: |
| 1.1 | 1.009 |
| 1.2 | 1.033 |
| 2.0 | 1.825 |
| 3.4 | 9.674 |

Judged at each minimum's own station, with `M(z)` from the solved cooling history:

| criterion | straight, flown | flared, capped |
| --- | ---: | ---: |
| **choking (correct)** | **>= 12 coils** | **>= 18 coils** |
| R10's loss cone at `R` = 1.096 | >= 24 | >= 24 |
| P8's existence test | >= 36 | >= 72 |

**P8's requirement drops by a factor of two to four**, and R10's own proposal — while much closer
than P8 — is still about twice as strict as the physics demands.

**R10's caveat about the flat shelf is correct and is the reason the flared winding needs more
coils than the straight one.** Every binding minimum lands at `z` = 0.6–2.8 m, on the 12 T shelf,
where there is no background gradient to swamp ripple and where the flow has almost no Mach margin
to spare. That is exactly the interaction R10 predicted, arrived at from the other direction.

**Reproduce:** `make analysis-nozzle-residence`.

---

### P18. R9: both of the paper's front numbers check out, and the cone construction is validated

**And the answer to "9 T or 12 T" is neither — it is 11.0 T, or 8.3 T if the sound speed holds.**

`sec:needle_through_fog` carries two figures cited to this repository — a 94 600 K shocked layer
and a 21.1 km/s spreading speed — that had never actually been reproduced here. They are both the
**strong-shock piston solution** on `eos_water`, and they are exact. A piston driven into cold gas
at `v` leaves the shocked gas moving at the piston speed, so `e_shocked = v^2/2` directly. At
45.58 km/s and fourfold compression of the flown bag that is 1038.8 MJ/kg, which inverts to
**94 632 K**, whose sound speed is **21.12 km/s**. Against 94 600 K and 21.1 km/s: **0.03% and
0.1%**, with `c_s` varying only 20.8–21.3 km/s across a 2×–16× compression bracket.

**The cone construction survives, for a reason ADR-0012 does not give.** ADR-0012 draws a straight
cone at the entry angle. Integrating the front instead moves two things that are *not* independent:
the projectile decelerates by `k+1`, which would open the cone faster — but `c_exp` is the shocked
layer's sound speed and the shock is driven by that same `v`, so both halves of `c_exp/v` fall
together. **Freezing `c_exp` while letting `v` fall is an artifact, not a correction**, and it moves
contact by up to 1.3 m in the wrong direction.

What survives is the part that does not scale: dissociation and ionisation absorb energy at
thresholds, so `T` rises more slowly than `v^2` and `c_exp/v` **falls** — 0.463 at 45.58 km/s to
0.29 near 10 km/s. The front opens slightly more slowly than the cone and touches **later**:

| | cone (ADR-0012) | frozen `c_exp` (artifact) | closed (correct) |
| --- | ---: | ---: | ---: |
| 45.58, sound speed, bag bore | 6.11 m | 5.29 m | **6.28 m** |
| 45.58, 1.9x, bag bore | 3.22 m | 2.95 m | **3.27 m** |

**The straight cone is right to within 0.3 m**, so ADR-0012's construction stands as arithmetic.

**The cap, with R15's liner as the wall** (see P19):

| spreading speed | wall | contact | **peak field demanded** |
| --- | --- | ---: | ---: |
| sound speed | bag bore 3.02 m | 6.22 m | 8.85 T |
| **sound speed** | **liner 3.50 m** | **7.29 m** | **8.25 T** |
| 1.9x (the bracket that forced 12 T) | bag bore 3.02 m | 3.27 m | 11.75 T |
| **1.9x** | **liner 3.50 m** | **3.83 m** | **10.95 T** |

**ADR-0013's proposed 11 T is confirmed almost exactly** (10.95 T), and if the sound-speed spread
is the right one the shelf is **8.25 T**, better than the 9 T R9 hoped for. ADR-0012's 12 T is
1.05 T conservative.

**And R9's own flag is worth keeping in the paper.** Taking the sound speed *is* conservative for
the coupling argument and anti-conservative for the field cap. That asymmetry should be stated
wherever the 21.1 km/s is used, because the same assumption is protective in one place and
exposing in the other.

**Reproduce:** `make analysis-nozzle-front`.

---

### P19. R15: the wall is the liner, `k` survives, and P9's 5 T *is* a bore-referenced standoff number

**All three of R15's questions, answered. Question 1 is the one that unblocks ADR-0013.**

**1. The standoff requirement is written against the liner.** The mist is inside the 3.0 m bag and
the clearance gap between bag and liner is **vacuum**. So a front that expands past the bag bore
sweeps no *less* mass — it already spans the bag's full cross-section, which is all there is to
sweep — and no *more*, because there is nothing out there. **Coupling is indifferent to the
excursion.** The only thing the field must prevent is plasma reaching the graphite. So ADR-0012's
own sentence ("the field is there to keep plasma off the liner") is the operative one, its table's
3.02 m is the stale one, and **ADR-0013 should be accepted.**

**2. `k` survives, and the excursion helps twice.** Beyond the swept-mass argument above,
expanding from 3.00 m to 3.50 m dilutes the front's pressure by `(3.5/3.0)^(-2 gamma)` = **0.60**
at `gamma` = 5/3, so the standoff demand at the liner is lower than at the bore for two independent
reasons. Our integration confirms the swept mass is identical to nine digits whether the wall is
placed at 3.02 m or 3.50 m, because the swept area is capped by the bag either way.

*The one caveat, and it is ours not yours:* the lateral excursion into the gap is momentum going
sideways, and against a hard wall that momentum would be reflected. The gap is 0.5 m against a 3 m
radius, so the volume involved is 1.36× and the effect is second-order. This is an argument, not a
resolved run; the resolved version is the rod-resolved plow that R4 asks for on other grounds.

**3. Yes — and this reopens the larger version you set aside.** The graded profile 20/12/9/5 T is
derived by setting the snowplow's pressure against `eq:bore_from_length`'s bore area **at every
station**, so every station of it, the 5 T exit included, is a standoff number evaluated at
3.02 m. R15 says "if P9's 5 T is in fact a standoff number evaluated at a 3.02 m bore, say so, and
the larger version comes back on the table." **It is, and it does.**

**And in the corrected framework the larger version is better than merely cheaper.** R15 notes
that letting the whole profile follow the flared bore gives 8.44 T peak and 2.87 T exit, with
`A/A*` improving from 2.40 to 2.94. P13 shows why that improvement is not a rounding: a larger
`A/A*` raises the exit Mach number *and* — because it also lowers the exit field — has to be
checked against `M_A`. We have not run that profile; it needs a new field grading, not a new
model. **We recommend the paper ask for it as a numbered item** rather than adopting it.

---

### P20. R14: do not adopt 672.9 m³. Our pair was never derived as a pair

**This is a correction to an adoption already made, so it is worth doing before it propagates.**

R14 adopts 23.8 m / 672.9 m³ on the strength of three of the paper's round figures matching it.
R14 also asks, correctly, "if your 23.8 m came from something we should be matching more carefully
than three round numbers, say so." **It did, and you should.**

**Provenance of our two numbers, which come from different places:**

- `expansion.FIELD_LENGTH = 23.8` came from the **paper's own** 660 m³ at the **paper's own**
  aspect ratio 4: `8 pi r^3 = 660` gives `r` = 2.9722 m and `l = 8r` = **23.777 m**.
- `expansion.CHAMBER_RADIUS = 3.0` came separately from the paper's quoted "3 m bore" and its
  28 m² cross-section.
- `BAG_RHO = 0.323` = 213/660, i.e. on the **660 m³** basis.

So our own radius and our own density sit on different bases, by 2% in area — a defect on our
side, now recorded. **We never asserted (23.8 m, 3.00 m) as a pair, and its 672.9 m³ contradicts
the 660 m³ that `PV = nR_gT` fixes.** The volume is the physics; the bore and the aspect ratio are
consequences of it. Adopting 672.9 inverts that.

**The four-way consistent solution, which neither repository found:**

| | paper as stated | paper as quoted | R14's adoption | **this** |
| --- | ---: | ---: | ---: | ---: |
| standoff volume | 660 m³ | — | 672.9 m³ | **660 m³** |
| length | 23 m | — | 23.8 m | **23.78 m** |
| bore | 3.022 m | "3.0 m" | 3.000 m | **2.972 m** |
| cross-section | 28.70 m² | "28 m²" | 28.27 m² | **27.75 m²** |
| aspect `l/2r` | 3.81 | "4" | 3.97 | **4.000** |

**All three of the paper's quoted round figures still round correctly, and the 660 m³ constraint
is preserved.** R14's adoption keeps a bore figure that was itself a rounding and breaks the one
number that is a physics requirement.

**What this costs to change on our side:** `CHAMBER_RADIUS` would move 3.000 → 2.972 m, a 1.9%
change in area. We have **not** made that change, because it would invalidate every artifact this
document cites and the effect on `eta_geom` is below the width of every bracket here. It is
recorded as owed, and it should be made in the same pass as whichever geometry the paper settles
on. **Tell us which and we will re-run on it.**

---

### P21. R12: `phi` is delivered — and it is the wrong quantity to multiply by

**The numbers R12 asked for, one per leg and branch.** Radiated power per unit mass is `e/t_rad`,
integrated over a parcel's own nozzle transit for the numerator and continued down the free jet on
the same isentrope (`fireball.py`'s clock, 45 deg half-angle) for the denominator.

| leg | branch | **`phi`** | transit | jet tail | in-nozzle radiated | of the plume's internal energy |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 45.58 | equilibrium | **0.853** | 2.72 ms | 10.69 ms | 3.55 MJ/kg | 4.2% |
| 45.58 | frozen | **0.993** | 2.84 ms | 10.12 ms | 0.99 MJ/kg | 1.2% |
| 56.53 | equilibrium | **0.734** | 2.22 ms | 7.88 ms | 12.03 MJ/kg | 9.1% |
| 56.53 | frozen | **0.991** | 2.30 ms | 8.20 ms | 1.77 MJ/kg | 1.3% |
| 65.00 | equilibrium | **0.211** | 1.95 ms | 6.63 ms | 6.10 MJ/kg | 3.5% |
| 65.00 | frozen | **0.989** | 2.02 ms | 10.31 ms | 2.19 MJ/kg | 1.2% |
| 75.00 | equilibrium | **0.041** | 1.69 ms | 5.69 ms | 4.35 MJ/kg | 1.9% |
| 75.00 | frozen | **0.985** | 1.74 ms | 8.90 ms | 2.65 MJ/kg | 1.1% |

**The frozen branch is `phi` ~ 0.99 flat. The equilibrium branch runs 0.04 to 0.85** and falls
steeply with closing speed. The mechanism is the one the ADR-0026 bracket is about: a frozen plume
has no recombination store to hand back, so it is at 3000-5300 K by the exit and has done nearly
all its radiating inside. An equilibrium plume is buffered at 11 700-16 200 K at the exit and goes
on radiating hard, for a long way, into nothing.

#### The important part: do not use `phi` as a multiplier

R12's construction is `liner load = radiated share x pulse energy x phi`. **That is a product of
two strongly anti-correlated numbers, and it amplifies a variation that is not in the physical
quantity.** The cases with the smallest `phi` are exactly the cases that radiate most in total — a
hot equilibrium plume radiates 106 MJ/kg overall and only 4% of it inside; a cold frozen plume
radiates 1.0 MJ/kg and 99% of it inside.

**What the liner actually sees is the product, and the product is nearly flat:**

    in-nozzle radiated energy = 1.0 to 12.0 MJ/kg across all eight cases,
                                and 1.0 to 4.4 MJ/kg for seven of them

against a `phi` that varies by a factor of **24**. So the fourth column of the table above, not
`phi`, is the number the paper should carry. **Give the liner the energy, not the fraction.**

#### Converted into the paper's own currency

At 238 kg of merged slug (213 kg of bag plus the 25 kg projectile) against a 62.9 GJ pulse — *the
arithmetic is stated so it can be redone against whatever pulse energy the leg actually has, which
this repository does not own*:

| branch | in-nozzle radiated share of pulse energy |
| --- | ---: |
| **frozen** | **0.38% – 1.00%** |
| **equilibrium** | **1.34% – 4.55%** |

**This brackets `tab:bag_sizing`'s own 1.0–3.6%** — computed independently, from a solved
expansion and a TOPS opacity table, and it lands on the paper's own range. It also sits well
inside P5's Jupiter gate of **7.95%**, in every case, on both branches.

**And it carries a double-count warning.** If `tab:bag_sizing`'s 1.0–3.6% is *already* an
in-nozzle figure — which its agreement with our in-nozzle number suggests — then multiplying it by
`phi` again subtracts the same downstream radiation twice. **Check what that column is a share of
before applying anything to it.**

#### What it does to the shield stack

Taking R12's own scaling (3.32 GW at a 3.6% share, and `T ~ load^{1/4}`), our worst case is the
56.53 km/s equilibrium leg at 4.55%:

| | R12's `phi` = 1.0, 3.6% | **our worst case (4.55%)** | **our best case (0.38%)** |
| --- | ---: | ---: | ---: |
| liner load | 3.32 GW | **4.20 GW** | **0.35 GW** |
| against the 8.28 GW graphite ceiling | 0.40x | **0.51x** | **0.04x** |
| liner equilibrium temperature | 3105 K | **~3290 K** | **~1770 K** |

**The passive-structure claim does not fail anywhere.** R12 identifies `phi` = 1.0 on the
equilibrium branch as "the only case in this whole exercise where the passive-structure claim
actually fails", putting the liner at 4300 K. **That case does not arise**: the equilibrium
branch's `phi` is 0.04–0.85, not 1.0, and its largest in-nozzle share is 4.55%, which puts the
liner near 3290 K against graphite's 3900 K.

The shield count is the paper's model rather than ours, but on these numbers it is
**branch-dependent and small**: essentially none on the frozen branch, and R12's own table would
want one to three on the equilibrium branch.

#### Two caveats, stated because they cut in opposite directions

1. **The denominator has not converged on the 75 km/s equilibrium case.** Its last decade of area
   ratio contributes 24.8% of the total, so carrying the free jet further would push `phi` below
   0.041. Every other case is converged to under 7%. This direction is *safe* for the liner.
2. **The perturbative treatment fails on that same case's denominator.** Integrating `e/t_rad`
   against an adiabat implies it radiates 45% of its internal energy down the jet, which it cannot
   do while staying on the adiabat, so its true `phi` is *higher* than 0.041 — the unsafe
   direction. **Both caveats attach to the same case and they fight.** The numerator is solidly
   perturbative everywhere (1.1–9.1% of internal energy), so the in-nozzle energy we are
   recommending the paper carry is unaffected by either.

**Reproduce:** `make analysis-nozzle-phi`.

---

### P22. R7: renamed, and the paper should go further than it proposes

`expansion.THROAT_RADIUS` is now **`CHAMBER_RADIUS`**, and the term is in `CONTEXT.md`'s glossary.
That is R7's first suggestion and we take it.

*A note on how, because it is the kind of thing that bites later:* the old name was **removed**
rather than aliased. An alias would still have read correctly while silently no longer
*controlling* anything — `audit/freeze_scaling_check.py` rebinds that constant to scale the
geometry, and against an alias the rebinding would have been a no-op that changed the audit's
answer without failing. A missing name fails loudly instead.

**We do not recommend R7's second suggestion**, that both repositories adopt "throat" for the 5 T
end. R1's own amendment — correctly — reframes this device as a **de Laval nozzle** with magnetic
walls. In that vocabulary "throat" is a term of art meaning the sonic station where `A/A*` = 1,
and in this nozzle **the sonic station is at the chamber**, at the 20 T (now 12 T) end. The 5 T end
is a Mach 2.4–3.4 exit. Standardising on "throat = exit" would have the paper call a Mach 3 station
"the throat" in a section that has just introduced de Laval language, which will mislead every
reader with fluid-dynamics training.

**Suggested convention: `chamber` for the strong-field end, `exit` for the weak-field end, and
reserve `throat` for the sonic station** (which coincides with the chamber). That is unambiguous in
both vocabularies and costs the paper one search-and-replace of "throat" to "exit" in
`sec:two_leg_nozzle`.

---

### P23. Lowering `k` does relieve the dissociation toll — and the paper already sits on the momentum optimum

**Raised by Seth, 2026-09-05, and it is the obvious lever nobody had priced.** If the bond bill is
a fixed 50.9 MJ/kg, spread the same collision energy over *less* mass and the bill is a smaller
slice of it. `toll.py`'s closed form says so directly:

    eta_chem = sqrt(1 - 2 phi E_B (1 + k) / w^2)

`k` enters only through `(1 + k)`, so **the lever is real and it is strong.** At 45.58 km/s:

| `k` | `eta_chem` | `eta_jet` at 0.93 | vs 0.775 | `T_stag` |
| ---: | ---: | ---: | --- | ---: |
| 2.0 | 0.924 | 0.859 | clears | 24 515 K |
| 4.0 | 0.869 | 0.808 | clears | 20 726 K |
| **5.2** | **0.833** | **0.775** | **clears, exactly** | ~19 500 K |
| 6.0 | 0.810 | 0.754 | misses | 18 017 K |
| **8.52 (flown)** | **0.730** | **0.679** | **misses** | 15 147 K |

**Dropping `k` from 8.52 to about 5.2 clears the target at the coldest point in the envelope.**
That is the honest headline, and it is a cheaper-sounding fix than an argon programme.

#### Three costs, and the first one is decisive

**1. The paper is already within 0.3% of the momentum optimum.** `toll.py`'s own bound for gross
exhaust momentum is `m w sqrt(1+k) x eta_chem`. Substituting the closed form, that is
`m w sqrt((1+k) - 2 phi E_B (1+k)^2/w^2)`, which has a maximum at

    (1 + k)* = w^2 / (4 phi E_B)     ->   k* = 9.20 at 45.58 km/s

**The paper flies 8.52.** Momentum at `k` = 8.52 is 2.2529 in those units against 2.2531 at the
optimum — the design is sitting on the peak. Moving to `k` = 5.2 costs about **9% of the delivered
momentum**.

So the lever improves `eta_jet` by shrinking what `eta_jet` is a ratio *of*. **You clear 0.775 by
making the target easier, not by delivering more push.** That is worth saying plainly, because a
table of `eta_jet` values cannot show it.

**2. The margin over the *physical* floor gets worse, not better.** `sec:mass_interest`'s
forward-thrust floor is `1/sqrt(1+k)`, which **rises** as `k` falls: 0.324 at `k` = 8.52 against
0.401 at `k` = 5.2. So the ratio of `eta_jet` to the floor — which is the threshold that actually
means something physically — goes from **2.10× down to 1.93×** while the headline number improves.
The 0.775 is a swept design point; the floor is physics.

**3. It costs radiation, because the whole mechanism is "make the plasma hotter".** Run on the
solved expansion at each `k`'s own stagnation temperature:

| `k` | `T_stag` | in-nozzle radiated, equilibrium | frozen |
| ---: | ---: | ---: | ---: |
| 8.52 | 15 147 K | 4.31 MJ/kg (1.00×) | 1.08 MJ/kg (1.00×) |
| 6.00 | 18 017 K | **12.71 MJ/kg (2.95×)** | 1.56 MJ/kg (1.45×) |
| 4.00 | 20 726 K | 8.44 MJ/kg (1.96×) | 1.97 MJ/kg (1.83×) |
| 2.00 | 24 515 K | 4.91 MJ/kg (1.14×) | 2.45 MJ/kg (2.27×) |

**1.5× to 3× more radiation onto the liner.** At `k` = 6 the equilibrium branch reaches
12.71 MJ/kg, which is the same load as P21's current worst case — a liner near 3340 K against
graphite's 3900 K. Affordable, but it spends the margin P21 just established. (The equilibrium
column is non-monotone because the radiated share peaks at the opacity crossover, which the
plume's temperature walks across as `k` moves; the frozen column is clean and monotone.)

#### The structural point, and it is the useful one

**Argon and low-`k` are the same trade.** Both work by putting more energy on fewer particles;
both raise the plume temperature; both pay in `T^4` radiation onto the liner. **There is one
currency here and both remaining levers spend it.** So they are not independent options to be
stacked — pricing them together is the right way to look at it, and R13's sub-ask 3 is the gate
on both.

#### The one thing we cannot settle, and it is the crux

Two figures of merit point in **opposite directions**, and which one governs is `aim`'s growth
model, not ours:

| figure of merit | scales as | wants |
| --- | --- | --- |
| momentum per **projectile thrown** (if cadence, aiming or tracking is the scarce thing) | `sqrt(1+k) x eta_chem` | **`k` = 9.2** — the flown value |
| momentum per **kg of total mass expended** (if propellant mass is scarce) | `eta_chem / sqrt(1+k)` | **`k` as small as possible** |
| `eta_jet` as the growth tables charge it | `eta_geom x eta_chem` | `k` small |

Per kg expended, `k` = 2 delivers **2.25×** what `k` = 8.52 does. Per projectile, it delivers
**0.71×**. **We are not going to adjudicate that** — the mass-interest architecture exists
precisely because propellant mass is meant to be cheap, which argues one way, and the cadence and
aiming constraints argue the other.

**And note the optimum is speed-dependent**, `k* = w^2/(4 phi E_B) - 1`, so it is **9.2 at
45.58 km/s and 26.6 at 75 km/s**. A single flown `k` across the envelope is a compromise, and if
`k` is a per-pulse dial — it is just bag mass over projectile mass — then the cold tail and the
departure burn want it moved in *opposite* directions.

**What we need back.** Which figure of merit the growth tables are actually maximising, and
whether `k` is allowed to vary per pulse. With those two answers this becomes a one-line
optimisation; without them it is a genuine fork.

---

### P24. The drift's bounce is real, already modelled, and worth a third of the overtake's impulse

**Raised by Seth, 2026-09-05, as a check on whether the chain accounts for two competing effects
of the merged centre of mass. It does — and stating them explicitly changes how close to 0.775
the design actually is.**

#### The two effects, and where each already lives

| effect | where it is in this repository | size at `k` = 8.52 |
| --- | --- | ---: |
| **The drift steals energy from the fireball.** Energy in bulk centre-of-mass motion is not available to the expansion. | `nozzle_ledger.drift_fraction(k)` = `f_d` = `1/(1+k)`; `reflection_baseline` uses `u = sqrt(1-f_d) v_g` for the thermal part | **10.5% of the pulse** |
| **The drift can be bounced.** On the overtake it points prograde — the wrong way — so reversing it is impulse the nozzle gets for free. | the sign flip in `nozzle_ledger.passthrough_floor` (`-sqrt(f_d)` on the overtake), so `nozzle_work = base + sqrt(f_d)`; and `toll.py`'s untolled `±1` | **+1.000 in units of `m w`** |

**Both scale with `f_d = 1/(1+k)`, so both grow as `k` falls** — which is exactly why they matter
to P23's lever and had to be checked before it.

#### The drift is transferred perfectly, not merely efficiently

This is the part worth stating in the paper. The drift momentum is `M V = m w` **by momentum
conservation alone**. It never passes through the nozzle as energy, so `eta_chem` cannot charge it
and `eta_geom` cannot misaim it. It enters the impulse twice and neither entry is tolled: once as
the incoming momentum `p_in` (the `±1`), and once as energy inside the `sqrt(1+k)` bound. The only
thing required of the hardware is that the plume be **turned at all** — not turned well.

#### Impulse per projectile, in units of `m w`

Ship frame, impulse = `p_in - p_out`, with `p_in` = `+m w` on the overtake (prograde, the wrong
way) and `-m w` head-on:

| `k` | overtake @ 45.58: jet | + drift | **total** | drift's share | head-on @ 75: jet | − drift | **total** |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2.00 | 1.600 | +1.0 | **2.600** | 38% | 1.684 | −1.0 | **0.684** |
| 4.00 | 1.943 | +1.0 | **2.943** | 34% | 2.132 | −1.0 | **1.132** |
| 5.20 | 2.077 | +1.0 | **3.077** | 32% | 2.346 | −1.0 | **1.346** |
| **8.52 (flown)** | 2.253 | +1.0 | **3.253** | **31%** | 2.807 | −1.0 | **1.807** |
| 12.00 | 2.171 | +1.0 | **3.171** | 32% | 3.153 | −1.0 | **2.153** |

#### What this does to "how far short are we"

**Nearly a third of the overtake's impulse never passes through the nozzle.** So the shortfall in
`eta_jet` is *not* the shortfall in delivered impulse:

| | at 45.58 km/s, `eta_geom` = 0.93 |
| --- | ---: |
| `eta_jet` reached | 0.679 against a 0.775 target — **12.4% short** |
| overtake impulse reached | 3.095 against 3.391 — **8.7% short** |

**The efficiency number overstates the miss by about half again.** P15's crossover analysis is
still the right way to size the problem, but when the paper says "misses 0.775" it should also say
what that costs in impulse, because on the leg where the miss happens they are different numbers.

#### And it changes P23's verdict on lowering `k`

Seth's question was whether the bounce compensates enough to avoid lowering `k` as far. **On the
overtake leg it substantially does; on the head-on leg it makes things worse.**

| `k` move | overtake impulse | (jet term alone) | **head-on impulse** |
| --- | ---: | ---: | ---: |
| 8.52 → 6.00 | 0.967× | 0.952× | **0.815×** |
| 8.52 → 5.20 | **0.946×** | 0.922× | **0.745×** |
| 8.52 → 4.00 | 0.905× | 0.862× | **0.627×** |
| 8.52 → 2.00 | 0.799× | 0.710× | **0.379×** |

The `+1` cushions about a third of the overtake's loss — dropping to `k` = 5.2 costs 5.4% rather
than 7.8%. **But the `−1` on the head-on leg amplifies it instead**, because shrinking the jet term
subtracts against a fixed debit: the same move costs **25.5%** there.

**So a single flown `k` cannot serve both legs.** Cutting `k` to clear 0.775 on leg 1's cold tail
would cost leg 2 a quarter of its impulse, and leg 2 has no problem to solve. This is now the
strongest argument for the per-pulse `k` question P23 asks `aim` — the two legs do not merely
prefer different `k`, they are pushed apart by a term that changes sign between them.

#### A mechanism correction we owe on leg 1, in the same family as P17

**The "cup closed at the ship end" is not a loss-cone mirror. It is a magnetic wall.** A mirror
holds particles by `mu` conservation, and P11 shows there is no `mu`; `(J x B) . B = 0`, so the
field exerts no force along a field line on an isotropic fluid. What actually stops the plume is
**magnetic pressure**: at low `beta` the field is a wall the flow cannot push through.

The good news is that the conclusion survives the mechanism change, and with room:

| `k` | `T_stag` | stagnation pressure | `B` needed to stand it off | vs the 12 T chamber | `beta` |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2.00 | 24 515 K | 16.21 MPa | 6.38 T | **1.9×** | 0.283 |
| 4.00 | 20 726 K | 11.83 MPa | 5.45 T | 2.2× | 0.206 |
| **8.52** | 15 147 K | 7.19 MPa | 4.25 T | **2.8×** | 0.125 |
| 12.00 | 10 503 K | 4.71 MPa | 3.44 T | 3.5× | 0.082 |

**The wall holds on every row**, and — usefully — the retired loss-cone picture pointed the same
way (its required mirror ratio rises from 1.45 to 2.11 across the same span), so this is a
correction to the *reason*, not to the answer. But note the third cost of lowering `k` that this
exposes: **the margin falls from 2.8× to 1.9× in field**, because a lower `k` is a hotter, higher
pressure plume pushing harder on the wall that has to turn it.

#### Is the bounce elastic?

Seth's own caution — "assume the reflection is not that efficient". At `beta` = 0.08–0.28 the field
is a **spring**: it is compressed by the arriving plume and pushes back, storing and returning the
energy, so the reflection itself should be close to elastic. **The loss channel is not the bounce,
it is the dwell** — a plume held at the closed end radiates for longer, and that is `phi` (P21),
which has not been computed for the reflection specifically. That is the thing to check before
calling the bounce free, and it is cheap: it is the same quadrature with a longer residence.

#### One loose end we cannot close from here

Our `eta_geom` is computed in the **expansion frame** — `expansion.py` starts "from the plume at
rest" — so it describes the thermal jet. `aim` multiplies it onto `sqrt(1+k)`, which is built from
the **total** energy including the drift. Whether the same `eta_geom` should apply to the drift's
share of that energy is not answerable here: the drift is turned by the magnetic wall, and the
thermal jet is turned by the de Laval flare, and those are different devices with different
efficiencies. **Confirm which factor multiplies which term**, or the chain is charging a nozzle
efficiency against momentum that never entered the nozzle.

---

### P25. 0.775 is a calibration to the literature, not a threshold — and the plate never beats the nozzle

**Raised by Seth, 2026-09-05: where does 0.775 come from, is falling short a *problem*, and would
a pusher plate be better below it?**

#### Where 0.775 comes from

The paper states it plainly and it is not a physical threshold. `sec:jet_efficiency`:

> Published solenoid results for drift-free plumes report collimation ratios of 0.65 to 0.85,
> which is 130 to 170% of the 0.5 available to them by reflection. **Our own target of 0.775 is
> 147% of the 0.529 available at `k` = 8.5.** We are asking a solenoid for what solenoids are
> already reported to deliver.

So 0.775 is **the literature's own collimation performance, rescaled to our reflection baseline**.
It is a statement about what a solenoid ought to manage, not about what the architecture needs.
The paper elsewhere concedes it sits above the best directly comparable simulation: *"the 0.775 we
require sits above the 0.6 to 0.7 the best case returns"* (Inatomi's solenoid).

#### The values that *are* thresholds, and where we sit against them

| threshold | value | what fails below it | where we are |
| --- | ---: | --- | --- |
| **forward thrust**, `1/sqrt(1+k)` | **0.324** | the nozzle pushes the craft **backward**; the chain returns nothing | 0.68–0.87 — **2.1 to 2.7× clear** |
| growth break-even | between 0.324 and ~0.45 (`tab:mass_interest_growth`'s zero rows sit at `eta_geom` = 0.25–0.30) | the craft runs forward but **loses mass** — "underwater on a mass interest" | `eta_geom` = 0.58–0.98, well clear |
| scale-case economics | ~0.70 | the cost rebuttal at fleet scale | cleared above ~52 km/s |
| early-generation economics | ~0.89 | the \$3200/kg worst-case rebuttal | not cleared anywhere |
| **the 0.775 target** | 0.775 | **nothing** — it is a calibration | cleared above ~56 km/s |

**So falling short of 0.775 is not a failure, it is a worse result.** It moves the design down
`tab:mass_interest_growth` — slower doubling, smaller return — and it does not approach anything
that breaks. The number that would break something is **0.324**, and the chain clears it by more
than a factor of two everywhere, on both legs, on both branches, with or without an extension.

**What the paper should change:** say which of these is which. At present 0.775 is quoted as a
requirement, and a reader cannot tell it from the 0.324 that actually is one. Given P24, the
honest framing is: *the architecture works above 0.324, grows above roughly 0.45, and matches the
solenoid literature above 0.775.*

#### Would a pusher plate be better below 0.775? No — and mostly it is not even available

Worth doing because this repository owns `f`: `CONCLUSION.md` puts it at **0.8**, defined as the
axial momentum delivered as a fraction of full-capture-perfect-bounce, so a plate's impulse is
`2f m w` = **1.600 m w**.

| `eta_jet` | overtake nozzle | vs the plate | head-on nozzle | vs the plate |
| ---: | ---: | ---: | ---: | ---: |
| 0.324 (thrust floor) | 2.000 | 1.25× | 0.000 | — |
| **0.679 (our 45.58 value)** | **3.095** | **1.93×** | 1.095 | 0.68× |
| 0.775 (the target) | 3.391 | 2.12× | 1.391 | 0.87× |
| 0.846 (our 75 value) | 3.610 | 2.26× | 1.610 | 1.01× |

**On the overtake the nozzle beats an `f` = 0.8 plate for any `eta_jet` above 0.194**, and we are
3.5× above that. Even at our shortfall the nozzle delivers **93% more impulse** than a perfect-as-
built plate. The plate never wins.

**On the head-on departure a plate cannot be used at all**, and the reason is structural rather
than numerical. The projectile arrives with momentum `-m w`; a plate returns it at `+f m w`, so the
impulse is `-(1+f) m w` — **it decelerates the craft.** Forward thrust on that leg requires
ejecting *more* momentum than arrived, and **only added mass can do that.** That is the reason the
departure burn needs a slug at all, and it is worth one sentence in the paper because it makes the
whole nozzle-versus-plate question disappear on that leg.

**And at these speeds the plate is not an option regardless.** This repository's `f` study runs
3.2–16 km/s; the nozzle legs run 45.58–75 km/s, where `sec:minimum_nozzle`'s "no wall survives the
pass" applies — the front's own radiative flux is 4.5 TW/m² against graphite's 13.1 MW/m²
(`snowplow.py`). The two architectures do not overlap in velocity, so the comparison above is a
sanity check rather than a live design choice.

---

### P26. Why our regime beats the literature's 0.34 — and it is a *regime* argument, not a geometry one

**Raised by Seth, 2026-09-05, and it is the right thing to want in the paper.** The paper already
gives four reasons our nozzle should beat Schilling's 0.34, and all four indict **his geometry**:
his field is a 32-strut cage with gaps (Ampère's law on an interior loop encloses no net current),
its maximum sits at the *closed* apex, plasma trapped there enters a resistive heating loop, and
his run stops with performance still climbing. The paper's summary is sharp and correct — *"the
distinction that matters is residence rather than resistivity"* — and it even concedes that per
interaction **his field grips more cleanly than ours**, `Rm` of 300–3e4 against our 39–650.

**What that argument cannot do is explain why we should beat the *good* nozzles in the same
literature** — Ahedo's 0.63–0.83 or Inatomi's 0.6–0.7, neither of which has a cage. The paper's
0.775 sits above both, and at present it is justified only by a scaling of collimation ratios.
**This batch supplies the missing reason, and it is fundamental rather than geometric.**

#### The advantage: our plume is a fluid, so the field is a wall rather than a maze

Every nozzle in that comparison is modelled **collisionless**. Ahedo's is explicitly "collisionless,
electron-magnetized, current-free and low `beta`". In that regime a plasma is a population of
*particles*, each conserving its own `mu`, and the field's *topology* governs everything: particles
mirror, trap in local `|B|` minima, and leak through gaps. Trapping is the dominant failure mode,
and Schilling's apex bottle is one instance of it.

**Our plume is not in that regime, and P11 measures the distance: `Kn` = 2.5e-7**, on the longest of
three collision channels. A parcel collides two to ten million times crossing the bore. It is a
**fluid**, and three things follow that have no counterpart in the collisionless literature:

1. **It cannot be magnetically trapped. At all.** `(J x B) . B = 0` identically, so for a
   scalar-pressure fluid the field exerts **no force along a field line**. The mirror force of
   kinetic theory is proportional to `p_perp - p_par`, and collisions erase that anisotropy in
   ~1e-10 s against a 2 ms transit. **Every trapping failure mode in that literature — loss cones,
   magnetic bottles, apex traps, ripple wells — is structurally unavailable to us** (P17 reaches
   the same conclusion from the ripple side, and withdraws P8's coil-count requirement because of
   it). The paper's "every gram has a downhill path out" is not a property of the graded profile.
   It is a property of the *regime*, and it would hold even if the profile were not monotone.
2. **The field is a wall, not a fence.** At `beta` = 0.013–0.15 the magnetic pressure exceeds the
   plasma's by 7 to 75×, so a continuous winding confines the plume the way a pipe confines gas.
   P24 shows the same thing doing the harder job at the closed end, standing off the full
   stagnation pressure with 2.8× margin in field.
3. **The bore is a pipe.** Traced through the real field, `<cos theta>` inside the winding is
   **0.999** — there is no divergence loss at all until the plume passes the last coil. The entire
   remaining loss is a detachment-surface question, which is a much smaller and better-posed
   problem than "how much plasma does the topology lose".

#### What the paper should say, and why it matters more than the cage argument

> The published numbers this section brackets against were all computed for collisionless plumes,
> where the field is a topology that particles navigate and trapping is the dominant loss. Our
> plume is not in that regime. At the flown bag density the Knudsen number is `2.5e-7`: a parcel
> collides millions of times crossing the bore, so the plasma is a continuum fluid with a scalar
> pressure. For such a fluid the Lorentz force has no component along a field line, so the field
> cannot hold the plume back at all — it can only wall it in, and at `beta` of 0.01 to 0.15 it walls
> it in with an order of magnitude to spare. The collimation figures in the literature are
> therefore not a ceiling on this device. They are measurements of a different problem, in which
> the nozzle has to win an argument with the field's topology that ours never has.

**Three reasons to prefer this to the cage argument**, which should be kept but demoted:

1. It explains the gap against the **good** nozzles, not only against the cage.
2. It is a statement about our own measured state rather than about someone else's design, so it
   cannot be answered by "he should have used a solenoid".
3. It is **falsifiable and already tested**: if `Kn` were of order 1 the argument collapses, and
   P11 is the run that checks it.

#### The honest cost of the argument, and it must be stated with it

**If their regime does not bound us from above, it does not support us from below either.** The
paper currently draws comfort from "we are asking a solenoid for what solenoids are already
reported to deliver". That sentence has to go if this argument is adopted: a collisionless
solenoid's 0.65–0.85 is neither a ceiling nor a floor for a continuum one. **Our `eta_geom` then
stands on our own solve and nothing else** — 0.58–0.86 as flown, 0.86–0.98 with R11's extension.
That is a stronger position but a lonelier one, and the paper should say so rather than keep both
arguments.

#### And it makes the argon case in one line

The same density that makes us collisional is what makes the water dissociation toll bite (P15).
**Our regime buys the nozzle and pays for it in chemistry.** Argon keeps the first and drops the
second, which is why P15 recommends it as the next build: it is not an incremental efficiency, it
is the only change that keeps the regime advantage without the regime's cost.

---

## What the replies got wrong

Three items, and two of them are retired by the replies' own amendment box.

### 1. R4's premise does not survive R1's amendment

R4 reprices our deferred N7 from "worth ~3% on `f_d`" to "worth a factor of 2 to 3 on `eta_geom`",
on the strength of a mass-versus-station profile changing the station-weighted `<alpha>`.

**Station weighting is a `mu`-framework argument.** It only bites if a parcel is stuck with the
field ratio of its *birth station* — which is exactly what `mu` conservation asserts and what R1's
own amendment retires. In a collisional nozzle a parcel's conversion is set by the total pressure
drop it falls through, not by its birth station's `B`. R1's amendment box concedes this for the
*levels* ("read it as a ranking rather than as levels") but keeps the station weighting as "the
reason the correction was found", and R4 is built on the part that does not survive.

**So R4's factor of 2–3 is withdrawn, and its own table's spread from 0.70–0.88 down to 0.24–0.30
does not exist.** The rod-resolved plow is still worth having — it is worth about 3% on `f_d`, as
our deferred table said, and it would settle R15's question 2 — but it is not the most valuable run
on the list.

**Where station weighting *does* survive is geometry**, and P16 uses it there: which flux tube a
parcel rides is fixed by where it starts, even though how much energy it converts is not. That
distinction is the whole of it, and it is worth stating in the paper because both repositories have
been sliding between the two senses.

### 2. R10's threshold is built from the wrong `alpha`, twice over

R10 sets its trapping threshold at `1/(1 - alpha)` = 1.096 from P1's `alpha` = 0.088. Two problems,
beyond the fact that P17 retires the loss cone entirely:

- `alpha` = 0.088 is the anisotropy of a **free expansion with no nozzle**, measured to answer
  N1. It is not the pitch-angle distribution of plasma in the bore, which is collisional and
  isotropic in the comoving frame.
- Even taking the loss cone at face value, a *collisional* population is continuously scattered
  into the loss cone, so residence is set by collisional diffusion rather than by a pitch-angle
  boundary.

The direction of R10's correction — depth matters, not existence — is right, and P17 keeps it.

### 3. R1's table keeps the station weighting it concedes for the levels

R1's `<alpha>` = 0.345 row and the `eta_geom` = 0.45–0.57 it produces are in the retired framework.
The amendment box handles the levels and not the weighting. Nothing in P12 depends on either.

---

## Smaller corrections, no argument attached

- **R2's paraxial concern was justified but the sign was wrong.** Evaluating `B_r = -(r/2)
  dB_z/dz` out to 6 m against a 3.5 m coil is indeed past where paraxial is safe, but the real
  solve fans *more*, not less — so the paraxial numbers were generous.
- **`fireball.DIVERGENCE_HALF_ANGLE_DEG` = 45°** is this repository's free-jet divergence, and it
  is the same 45° the fan reaches at 2 exit radii. The agreement is not by construction.
- **P3's detachment window is doing more work than it looks.** `<cos theta>` runs 0.95 at 1.0 exit
  radii and 0.69 at 2.0, so the same bracket that P3 reported as a modest uncertainty is a 27%
  swing in `eta_geom`. Narrowing `M_A` at the exit is worth more than it appears.

---

## Deferred, with cost

| item | why | what it needs |
| --- | --- | --- |
| **R13, argon** | **a new species, not a new calculation** | an `eos_argon` beside `eos_water`: Ar I–III Saha ladder with the ionisation potentials and statistical weights, a frozen-composition variant for the ADR-0026 bracket, and a sound speed. Then the whole `expansion` → `detachment` → `fireball` → `radiance` chain re-run on it, plus `conductivity` for sub-ask 4's field leak. The Saha machinery exists and is tested; the species does not. **Recommended as the next build** — see P15. |
| **R8, the full MHD** | **not deferred — declined, with the physical reason** (P14) | Would still need a prescribed-inflow BC, an inward-normal immersed wall for `r = r_w(z)` (the IBM only does upward-facing `z = z_s(r)`), table EOS wired into `state.rs`/`riemann.rs`, and radiation the 2D track lacks. **We are not asking for it**, because the answer does not turn on it: no structure is needed for detachment, and a wall cannot act in the sub-Alfvénic stretch where the pointing loss occurs. Reopen only if the magnetic extension is refused on cost. |
| **R4 / R5, the wall-interaction half** | same machinery | as above. Reduced in value by [What the replies got wrong](#what-the-replies-got-wrong) item 1, but still owed. |
| **The `1/r`-graded profile of R15's third question** | a new field grading | not modelled here. It is a re-derivation of the 20/12/9/5 T profile against a flared wall, then this document's whole chain re-run on it. Cheap in machinery, not free in runs. See P19. |

---

## Provenance

Every figure above was produced by running the code on 2026-09-05.

| run | covers |
| --- | --- |
| `make analysis-continuum` | P11; the Knudsen number, the three channels, collisions per transit |
| `make analysis-nozzle-fluxtube` | P12's `<cos theta>`, P16; flux tubes on the Biot–Savart solve, the downstream fan |
| `make analysis-nozzle-extension` | P12, P13, P14; the continuum `eta_geom`, the cap, the extension sweep |
| `make analysis-nozzle-residence` | P17; the choking criterion and the coil counts |
| `make analysis-nozzle-front` | P18, P19's question 2; the shock closure and the contact stations |
| `make analysis-nozzle-phi` | P21; the radiated-share quadrature, and P23's radiative cost of low `k` |
| `puffsat.toll` (closed form) + `puffsat.radiance` | P15's crossover speeds and P23's `k` sweep |
| `make analysis-replies` | all of the above, in dependency order |

Artifacts: `data/results/continuum_check.csv`, `nozzle_fluxtube.csv`, `nozzle_extension.csv`,
`nozzle_residence.csv`, `nozzle_front.csv`, `nozzle_phi.csv`.

**Validation.** `extension.continued_history` reproduces the shipped `cooling_history.csv` exactly
at the flown geometry — `M` = 2.685–3.441, `beta` = 0.013–0.073, `M_A` = 0.353–0.578, which are
P3's published values — before it is run anywhere new. `front.shock_state` reproduces the paper's
own 94 600 K and 21.1 km/s to 0.03% and 0.1%. The flux-tube tracer is pinned by the analytic
invariant `B r^2 = const`, and the continuum `eta_geom` is pinned against R1's own
`gamma M^2` form to machine precision.

**Errors made and caught while doing this work**, recorded because each would have shipped a wrong
number:

1. **The front's spreading speed held fixed while the projectile decelerated.** That made
   `dr/dz = c_exp/v` steepen and reported the front touching the wall *sooner* than ADR-0012's
   cone — the opposite of the truth. `c_exp` is the shocked layer's sound speed and the shock is
   driven by the same `v`; the two nearly cancel. Caught by asking where 21.1 km/s comes from,
   which is also how P18's confirmation was found.
2. **Flux tubes launched at the chamber plane**, i.e. modelling a steady nozzle with one inlet
   rather than a column that already fills the bore. Worth a factor of three on P16, and it is the
   difference between agreeing with R5 and contradicting it.
3. **A predicted sign that the run contradicted.** The module docstring asserted that ADR-0012's
   cap must cost `eta_geom` by cutting the Mach number. It does cut the Mach number and it improves
   `eta_geom` anyway (P13). The docstring was corrected to the run rather than the reverse.
4. **The cone divergence factor applied to extensions whose flow never reaches `M_A` = 1.** R11's
   claim is that the extension moves the crossing *inside* the flare; applying its `<cos theta>`
   where the crossing does not land there would have assumed the conclusion. Now checked per case.
