# The nozzle asks N1–N7, answered

Answers to `katzseth22202/Balloon-Pulse-Propulsion` @ `36080e1`,
`docs/nozzle_asks_for_impact_sim.md` (raised 2026-09-03). Worked in
`katzseth22202/puffsat_impact_simulation` the same day.

**This document is written to be copied into the paper repository and worked there.** Section
["What the paper should change"](#what-the-paper-should-change) is self-contained: every number
needed to make every edit is stated in it, because an agent in the paper repository cannot run the
code that produced them. Every figure here came from a run, not from recall; the
[provenance table](#provenance) says which.

---

## Status

| # | verdict | does the paper change? |
| --- | --- | --- |
| **N1** | **Answered, and it overturns an assumption.** `alpha` = 0.088, not the 1/3 assumed. The cause is the column geometry, proved by a spherical control. | **yes** (P1) |
| **N2** | **Answered, and it clears.** The flare recovers a directed jet; `eta_geom` = 0.70–0.88 at the flown flare, 371% of the measured baseline. | **yes** (P2) — the good news |
| **N3** | **Partly.** `eps_b` = 1.5 is a tautology no run can move; `beta(z)` = 0.013–0.073 delivered; a new realizability finding. Liner mass fraction deferred. | **yes** (P6, P7) |
| **N4** | **Flown case answered; near-Sun declined for cause.** The radiative bracket is 10x wide, not 1.0–3.6%. Two of the paper's own numbers disagree by 10x. | **yes** (P4, P5) |
| **N5** | **Answered as a constraint on windings**, not a yes or no: ≥ 36 coils. Residence times deferred. | **yes** (P8) |
| **N6** | **Answered, and it overturns the single-station check.** `M_A` never crosses 1 inside the column. A field window exists that fixes it. | **yes** (P3, P9) |
| **N7** | **The sign correction delivered**; the realised snowplow value still owed. | **yes** (P10) |

**Ten paper edits fall out, P1–P10.** P1, P3 and P4 are the ones that change a claim rather than a
number. P2 is the only one that makes the paper *stronger*.

**Three of the asks were themselves mistaken**, and those corrections are in
[What the ask got wrong](#what-the-ask-got-wrong). They matter because two of them would have sent
this work at the wrong target.

---

## What the paper should change

### P1. `eq:reflection_baseline` assumes an isotropy the flown geometry does not have

**Locate:** `grep -n "eq:reflection_baseline" templateArxiv.tex`, the derivation paragraph — "For an
isotropic expansion $\cos\theta$ is spread evenly over $[-1,1]$".

**Now:** the thermal remainder is taken as isotropic, `alpha = <v_z^2>/<v^2>` = 1/3, giving a
drift-free baseline of 0.461 and 0.529 at `k` = 8.5.

**Should be:** state that the *flown column* is not isotropic. Simulating the merged fireball's
free expansion into vacuum gives

| body | `alpha` | drift-free baseline `sqrt(2 alpha/pi)` |
| --- | ---: | ---: |
| **23 m x 3.02 m column (flown)** | **0.088** | **0.237** |
| 5.4 m sphere (the paper's original bag) | 0.360 | 0.479 |
| isotropic, as assumed | 0.333 | 0.461 |

**The cause is the bag's shape, and the control proves it.** Same 238 kg, same density, same
energy — only the geometry differs, and `alpha` moves by a factor of **4.1**. A 23 m x 3.02 m
column is nearly an infinite cylinder: free expansion follows the steepest pressure gradient, which
is radial everywhere except near the two ends, and the ends hold little of the mass. **58% of the
expelled mass ends up within 12 degrees of the transverse plane.**

**Why this matters beyond one number.** `sec:needle_through_fog` moves the bag from a 5.4 m sphere
to a 23 m column because "length is the cheap dimension for a rocket to carry and diameter is the
expensive one" — a launch-fairing decision. Isotropy is *exact* for the sphere and false for the
column, so **a packing choice silently invalidated the assumption `eq:reflection_baseline` rests
on.** The paper should say so, because a reader who assumes the sphere's isotropy for the column
will get the baseline wrong by a factor of two.

**This does not damage the design** — see P2. It damages the *yardstick*.

**Reproduce:** `cargo test --release -p euler2d --test merge_expansion -- --nocapture report`.

---

### P2. `eta_geom` is now bounded, and the flare is what bounds it

**Priority: this is the item that makes the paper stronger, and it answers the sentence the paper
flags as its own largest gap.**

**Locate:** `grep -n "eta_geom" templateArxiv.tex`, and `sec:jet_efficiency`'s "Nothing in this
paper or in either companion repository bounds $\eta_{\mathrm{geom}}$, and it is the largest
remaining uncertainty in the growth chain."

**Should be:** that sentence can go. In a diverging field `mu = v_perp^2/2B` is an adiabatic
invariant, so as `B` falls through the flare `v_perp^2` falls in proportion and energy conservation
puts the difference into `v_par`. With flux conservation tying `B_exit/B_throat` to `1/(A/A*)`:

    alpha_exit = 1 - (1 - alpha_0)/(A/A*)

| `A/A*` | `alpha` at the exit | `eta_geom` (Gaussian spread) | `eta_geom` (directed exhaust) |
| ---: | ---: | ---: | ---: |
| 2.0 | 0.544 | 0.589 | 0.738 |
| **4.0 (flown)** | **0.772** | **0.701** | **0.879** |
| **11.3 (see P9)** | **0.919** | **0.765** | **0.959** |
| 20.0 | 0.954 | 0.780 | 0.977 |

**`alpha_0` enters only through `(1 - alpha_0)`, which is at most 1**, so a deep enough flare
recovers a directed jet from *any* starting anisotropy. **The flare is not a mitigation for P1; it
is the mechanism that undoes it.** Against the ask's own calibration — ≥130% of baseline supports
the paper, <100% is a serious problem — the flown flare returns **371% of the measured baseline**
and **191% of the isotropic one**.

**The consequence that should be in the text.** `eta_jet = eta_chem x eta_geom` clears the paper's
own 0.775 target at 75 km/s (0.638–0.800 at the flown flare) and misses at 45.58 km/s
(0.512–0.642) — but `eta_geom` is *identical* at both speeds, because the flare does not care how
hot the pulse was. **Where the chain misses, it is `eta_chem` (0.910 → 0.731) doing it.** The
largest remaining uncertainty in the growth chain moves off the nozzle and onto cold-end chemistry.

**The assumption behind it, checked:** `mu` is invariant only while the field varies slowly over a
gyroradius. `r_gyro/L` runs 1.4e-6 to 3.0e-6 across the column, so the guiding-centre picture holds
by six orders of magnitude.

**Reproduce:** `make analysis-nozzle-jet`.

---

### P3. The detachment check does not survive being done per station

**Priority: the only place the paper prints a number that a consistent calculation contradicts.**

**Locate:** `grep -n "1.63" templateArxiv.tex`, in `sec:jet_efficiency`'s Alfvén paragraph.

**Now:** "The plume leaves at $v_g\eta_{\mathrm{chem}} = \SI{10.8}{\kilo\meter\per\second}$, which is
Mach \num{1.63} in Alfvén terms. The hottest pulse gives ... \num{2.06}." Read as clearing the
detachment condition.

**Should be:** computed consistently at each station of the solved expansion, `M_A` rises
monotonically from **0.12–0.22 at the throat to 0.35–0.58 at the exit and never crosses 1.**

**Why the paper's number differs, and it reproduces exactly.** The paper divides the *exit* speed
(10.8 km/s) by an Alfvén speed built from the *bag* density (0.323 kg/m³) and the bag's standoff
field (4.1 T). The plume thins by 13x between those two states and `v_A ~ rho^{-1/2}`, so using the
pre-expansion density with the post-expansion speed inflates `M_A` by
`sqrt(0.323/0.0251)` = **3.59** — and `0.58 x 3.59 = 2.06`, the paper's hot-pulse figure to three
figures.

Underneath it, `eq:alfven_from_standoff`'s `v_A = sqrt(2 R_g T/Mbar)` follows from `beta = 1`. The
paper's two routes to `v_A` agree to 3% and it reads that as corroboration, but they are not
independent: `tab:bag_sizing`'s 4.1 T **is** the standoff field at that state, so both evaluate the
same assumption. The solved expansion runs at `beta` = 0.013–0.073, because the field is graded
against the **collision's** snowplow pressure (159 MPa at 1 m) and then asked to steer the
**expansion** that follows (10 MPa at the throat, 0.7 MPa at the exit) — roughly 15–75x
over-strength for the flow it steers, which is exactly what holds the plume sub-Alfvénic.

**What survives.** The downstream argument — past the last coil a plume's pressure falls as `R^-5`
against a vacuum field's `R^-6` — is untouched, and on its own scalings the crossing lands
**1.44–2.00 exit radii out**. What dies is the claim that detachment happens *inside* the nozzle,
and the number offered for it. **Say the crossing is downstream, and quote 1.4–2.0 exit radii.**

**Reproduce:** `make analysis-nozzle-detachment`.

---

### P4. Two of the paper's own routes to the same Jupiter heat load differ by 10x

**Locate:** `grep -n "42.6" templateArxiv.tex` (the intercepted flash) and `tab:bag_sizing`'s
radiated-share column.

**The arithmetic.** The Jupiter-only chain runs at **2 pulses per second** (stated twice, in the
shell-melt and liner-ablation appendices) and 62.9 GJ per pulse, so **125.9 GW average burn power**.

| route | value | implied product (share x sky fraction) |
| --- | ---: | ---: |
| booked intercepted flash (`tab:bag_state`) | **42.6 MW** | 3.4e-4 |
| implied by `tab:bag_sizing`'s 3.6% at "a tenth of the sky" | **453 MW** | 3.6e-3 |

**They are the same physical quantity and they disagree by 10.6x.** Backing the sky fraction out of
the booked figure gives **0.9% (hot) to 3.4% (cold)**, against the "of order a tenth" the gate is
stated with. Either the sky fraction is far below a tenth for this nozzle, or the booked flash is
low. **Say which.**

**The load-bearing conclusion is unaffected:** the booked figure sits **23x inside** the Jupiter
gate, so the flown case clears comfortably either way. The discrepancy matters because a reader can
otherwise reach the opposite conclusion — see P5.

**Reproduce:** `make analysis-nozzle-ledger`.

---

### P5. The `1e-4` passive-structure gate is a near-Sun number and should say so

**Locate:** `grep -n "passive structure" templateArxiv.tex`, `sec:minimum_nozzle`.

**Now:** "The product therefore has to sit near \num{e-4} for the nozzle to stay a passive
structure. With coils filling of order a tenth of the sky, radiative escape must remain below about
a tenth of a percent." Every marker in that paragraph is near-Sun — 4 solar radii, 618 km/s
retrograde arrival, the 1000 t ship, 155 km/s effective exhaust.

**Should be:** attach the gate explicitly to the near-Sun burn. **The gate is a power balance and
scales inversely with burn power**, and the Jupiter-only chain runs 79x below the near-Sun case:

| case | burn power | gate product | allowed radiated share |
| --- | ---: | ---: | ---: |
| near-Sun periapsis burn | 10 000 GW | 1.0e-4 | 0.10% |
| **Jupiter-only, head-on departure** | **125.9 GW** | **7.9e-3** | **7.95%** |

Without this, a reader who finds the 0.1% gate and `tab:bag_sizing`'s 1.0–3.6% in the same section
concludes the design fails by 10–36x. It does not. **The paper never makes that transfer; the ask
did, which is why it is worth closing off explicitly.**

**And the flown radiative share is a wider bracket than the paper states.** Integrating `sigma T^4`
over the *solved* cooling history rather than holding `T` fixed:

| closing speed | equilibrium branch | frozen branch | gate at that speed |
| ---: | ---: | ---: | ---: |
| 75.00 | 2.59% | 1.30% | 7.95% |
| 65.00 | 4.94% | 1.41% | 10.58% |
| **56.53** | **13.25%** | 1.52% | **13.99%** |
| 45.58 | 5.18% | 1.29% | 21.51% |

`tab:bag_sizing` quotes **1.0–3.6%**, which spans the frozen branch only; the equilibrium branch
reaches 3.7x the paper's top figure. The mechanism is recombination reheating — equilibrium exits at
16 224 K from a 26 200 K start, frozen at 5296 K, and radiation goes as `T^4`. **The 56.53 km/s
equilibrium cell clears its gate by 6%**, and is also where the adiabatic assumption is thinnest
(`min t_rad/t` = 2.7 against 7–29 elsewhere). Quote the bracket, not the frozen half of it.

*Denominator note:* the solved figures are shares of *internal* energy; the gate is on *pulse*
energy. The conversion is `(1 - f_d)` = 0.895, moving 13.25% to 11.9% and the margin from 1.06x to
1.18x.

**Reproduce:** `make analysis-expansion` (the bracket), `make analysis-nozzle-ledger` (the gates).

---

### P6. The chamber field the graded profile asks for is not realizable by a simple solenoid

**Locate:** `grep -n "20}{\\\\tesla}" templateArxiv.tex`, the graded profile at the end of
`sec:needle_through_fog`.

**First, a finding that makes the profile easier to work with:** the four stated stations are not
four independent requirements but **one power law, `B = 19.80 z^-0.4405` to within 1.7%**. Through
standoff (`B^2 = 2 mu0 p`) that says the snowplow pressure falls as `p ~ z^-0.881`, very nearly
`1/z` — which is what a front spreading one-dimensionally down a fixed bore should do, and a
self-consistency check the paper never states. It is worth stating.

**The problem:** a solenoid smooths field structure over about its own radius, and the profile
demands a 40% drop between 1 m and 3 m — a 2 m gradient — from coils about 3.5 m in radius. An
all-positive winding cannot deliver it:

| coil radius | delivered at 1 m | shortfall | `beta` at the chamber |
| ---: | ---: | ---: | ---: |
| 1.5 m | 17.23 T | 13.0% | 1.32 |
| 2.5 m | 16.26 T | 17.9% | 1.48 |
| **3.5 m** | **15.42 T** | **22.1%** | **1.65** |
| 5.0 m | 14.33 T | 27.6% | 1.91 |

The shortfall scaling with coil radius is the signature of the mechanism, not an artifact.
**Undershooting `B` is overshooting `beta`**, at the highest-pressure station in the column
(159 MPa). Matching the profile exactly requires **counter-wound coils** — which is what an
unconstrained fit reaches for, and what nobody builds.

**Scope, stated precisely because it narrows the finding:** this `beta` = 1.65 is against the
*collision's* snowplow pressure, where the field is at standoff by construction. Against the much
thinner *expansion* (P3) the same shortfall gives `beta` = 0.10, a 10x margin. **The chamber
shortfall binds during the collision, not during the expansion.**

**What helps.** Not more current (a scale cannot fix a shape: 1.28x current reaches the chamber but
costs 1.65x field energy and over-provisions 29–43% everywhere else, moving the hot-pulse virial
structure from 10–30 t to 16.5–49.4 t — and better superconductors reduce the *conductor* term, not
the structure term that grows). Not a thinner plume (it lowers demanded and delivered field
together, leaving the ratio untouched, while the wider bore it needs widens the brush). What helps
is **geometry**: a longer column moves `beta` 1.65 → 1.43 (32 m) → 1.26 (50 m) because
`eq:bore_from_length` ties length to a narrower bore, and a **converging bore** lets the coils be
small exactly where the gradient is steep. The latter is blocked on P7.

**Reproduce:** `make analysis-nozzle-field`.

---

### P7. `eq:bore_from_length`'s cylinder and the area-ratio-4 nozzle are different hardware

**This one needs a decision, not a calculation, and several other items hang off it.**

**Locate:** `eq:bore_from_length` (`r = sqrt(V/(pi \ell))`) against the expansion's area ratio.

`eq:bore_from_length` gives a **constant 3.02 m cylinder** for a 23 m column of 660 m³. The
expansion model runs an **area ratio of 1 → 4**, and P9 wants 11.3–14.8. A flux tube that expands
fourfold has twice the radius; a cylinder does not have room for it.

**Which is real decides three things:** whether the cheap fix to P6 exists (small coils need a
converging bore), what the flare in P2 can actually be, and whether P9's field window is reachable.
**State the bore's shape.**

---

### P8. The residence claim is about the on-axis profile, and the answer is a constraint on windings

**Locate:** `grep -n "no local minimum" templateArxiv.tex`, `sec:jet_efficiency`.

**Now:** "The graded profile of \SI{20}{\tesla} at the chamber falling to \SI{5}{\tesla} at the exit
has no local minimum anywhere along it, so nothing can sit in it ... Every gram has a downhill path
out."

**Should be:** the claim holds, **conditional on the winding**, and the paper specifies none. The
entire magnet specification is four field values; there is no coil count, no radii, no currents, no
turn count. Tracing `|B|` off-axis over a swept winding family:

| coils over the column | on-axis | r = 1.0 m | r = 2.0 m | r = 2.5 m |
| ---: | --- | --- | --- | --- |
| 12 | clean | **2 minima** | **7** | **7**, mirror ratio 1.36 |
| 18 | clean | clean | **5** | **11** |
| 24 | clean | clean | clean | **11** |
| **≥ 36** | clean | clean | clean | **clean** |

**At ≥ 36 coils — roughly one per metre at a 3.5 m coil radius — there is no interior `|B|` minimum
anywhere out to 2.5 m, and the claim stands.** Below that, traps appear and reach inward as the
winding thins.

**The part that should change the sentence:** the on-axis column is clean at **every** coil count,
including the ones that trap badly. The argument as written could not have detected the failure mode
even where it occurs. Say the profile is monotonic *and* that the winding is dense enough that the
off-axis field is too.

There is also a tension to resolve in the same passage: residence is defended with "a solenoid's
winding is **continuous**, so the bore is walled by field rather than fenced by it" — but a
*uniform* solenoid cannot produce a graded profile at all. Grading means varying turn density, which
is what makes a winding discrete-able in the first place.

*Caveat to carry:* this is vacuum magnetostatics. At `beta` ~ 1 the plume pushes field out, which
deepens ripple. A trap found here is real; the absence of one is not proof.

**Reproduce:** `make analysis-nozzle-field`.

---

### P9. A field window exists that fixes P3, and the flown flare is not in it

**Do this in the same pass as P3** — it is the other half of the same finding.

`M_A = M_sonic sqrt(gamma/2) sqrt(beta)`: field and density both cancel out of the ratio. At
standoff (`beta` = 1) this is `M_A = M_sonic/1.095`, **exactly the paper's own "the Alfvén surface
sits a tenth of the way past the sonic throat"**. So the paper's physics was never wrong — only its
`beta`. And `M_A` carries `sqrt(beta)`, so a 15x over-strength field costs a factor of 4.

At every exit there is a window between the least field that still stands the plume off the wall
and the most that still allows `M_A` = 1:

| leg | contain (min) | release (max) | design | `M_A` now | `M_A` at standoff |
| --- | ---: | ---: | ---: | ---: | ---: |
| 45.58 equilibrium | 0.66 T | 1.76 T | 5.00 T | 0.35 | 2.66 |
| 75 equilibrium | 1.35 T | 2.88 T | 5.00 T | 0.58 | 2.13 |
| frozen legs | 0.57–0.92 T | 1.77–2.89 T | 5.00 T | 0.35–0.58 | 3.14 |

**The window exists on every leg and the design sits above it, not inside it.** Since flux
conservation ties the exit field to the flare, the fix is geometric: **`A/A*` between 11.3 and 14.8
gives every leg containment *and* detachment**, against the flown 4. The lower bound is set by
detachment on the coldest equilibrium leg, the upper by containment on the hottest.

**The cost, which is not optional.** One static field cannot be strong for the collision and weak
for the expansion at the same station, and the paper requires it static in time. A field weakened to
the release ceiling sits **3–8x under the collision's snowplow pressure** there, so the snowplow gets
past. Whether that is affordable in ablation is a collision-phase calculation nobody has done.

**Reproduce:** `make analysis-nozzle-detachment`.

---

### P10. `eq:reflection_baseline` is applied to both legs, and the drift enters with opposite signs

**Locate:** `grep -n "f_d = 1/(1+k)\|sec:two_leg_nozzle" templateArxiv.tex`.

**Now:** `eq:reflection_baseline` is derived once in a head-on framing, `f_d = 1/(1+k)` is
substituted, and `sec:two_leg_nozzle` never cites it — but the cycle "collides twice, and the two
collisions face opposite ways."

**Should be:** the *magnitude* is the same on both legs (`|-V + u mu|` has the same distribution as
`|V + u mu|`), but the floor it is measured against is not:

| leg | `f_d` | baseline | pass-through floor | work the nozzle must do |
| --- | ---: | ---: | ---: | ---: |
| head-on (departure) | 0.105 | 0.529 | **+0.324** | 0.204 |
| overtake (growth push) | 0.105 | 0.529 | **−0.324** | **0.853** |

**A nozzle doing nothing returns `+sqrt(f_d)` head-on and `-sqrt(f_d)` on the overtake**, so the
overtake nozzle does **4.18x the work** to reach the same baseline. The paper's companion
`CONTEXT.md` has this; the paper does not.

*This does not change P2's conclusion:* `mu` conversion is indifferent to the drift's sign, and
leg 1's mirror reflects the pancake essentially completely. Including the drift in the pitch angle
drops `sin^2(theta)` from 0.912 to 0.690 — a quarter of the margin — but the minimum mirror ratio
that still turns the flown plume is **1.45**, below anything the graded column approaches. **The
overtake mirror is not a binding constraint.**

**Still owed:** the *realised* `f_d` from a snowplow merge, against the clean-merge 0.105. The run
behind P1 takes `f_d` as an input (a uniform merged body), so it structurally cannot produce it.
That needs a rod-resolved merge — see [Deferred](#deferred-with-cost).

**Reproduce:** `make analysis-nozzle-ledger`.

---

## Smaller corrections, no argument attached

- **`k` = 8.52, not 8.5.** The flown masses (213 kg / 25 kg) give 8.52. A 0.2% inconsistency in the
  paper's own numbers; harmless, but the masses are what define the physics.
- **Column length.** This repo carries `FIELD_LENGTH` = 23.8 m against the paper's 23.0 m, and a
  volume of 672.9 m³ against 660 m³. Unreconciled; pick one.
- **`THROAT_RADIUS` = 3.0 m is correctly a *radius*.** The paper's own 28 m² cross-section and
  `r = sqrt(660/(pi x 23))` = 3.02 both confirm it. Recorded so it is not re-litigated.
- **`tab:bag_sizing`'s 0.323 kg/m³ is the slug alone** (213 kg over 660 m³), before the projectile.
  The *merged* 238 kg in the same volume is 0.361.

---

## What the ask got wrong

Three items in `nozzle_asks_for_impact_sim.md` itself, two of which would have aimed this work at
the wrong target.

1. **N1's premise.** "A diagnostic on runs that probably already exist" — it is not. This repo's
   nozzle strand is quasi-1D with `v_r` identically zero, so N1, N5's off-axis check and N3's liner
   fraction are all *radial* questions posed to a model with no radial dimension. All three needed
   new machinery.
2. **N4 carries the near-Sun gate across to Jupiter numbers.** It opens with the `1e-4` product,
   then asks for the radiated share "at both ends of the burn" citing `tab:bag_sizing`'s flown bag.
   The paper never makes that transfer (P5). Had it not been caught, the flown case would have been
   reported as failing by 10–36x.
3. **"The burn's speed range of 45.58–75 km/s" is two different legs.** 45.58 is the cold end of the
   *growth push* (leg 1, an overtake, where closing speed *falls*); ~75 km/s is the *departure burn*
   (leg 2, head-on). Also, 45.58 is a **3S** number — the cold end of the slowest three-synodic
   cycle — while the ask reads as if one cadence.

*Not* a correction, and recorded because it cost time: `k` = 8.5 is right for both legs and both
cadences. The `8.21–8.69 / 7.43–7.56` column in the companion `CONTEXT.md` is a **pusher-plate
payload boost ratio** on leg 1, not a slug ratio.

---

## Deferred, with cost

| item | why | what it needs |
| --- | --- | --- |
| **N3, mass fraction reaching the liner** | needs the plume against a wall, which no current kernel does | full-bore 2D: a prescribed-inflow BC (`euler2d` has none), an inward-normal immersed wall for `r = r_w(z)` (the IBM only does upward-facing `z = z_s(r)`), table EOS wired into `state.rs`/`riemann.rs` (a rewrite, not a patch), and radiation the 2D track deliberately lacks. Mesh cost is *not* the constraint — a bore run is tens of seconds. |
| **N5, residence times and mass left in the bore** | same | same machinery |
| **N7, the realised `f_d`** | the P1 run starts from an already-merged uniform body | a rod-resolved snowplow merge: 25 kg rod into the graded column, resolving the front |
| **N4, the 100 g near-Sun case** | **declined for cause** | out of scope by three prior decisions (`puffsat_impact_sim_design.md` twice, ADR-0010). Physically ~500 km/s is ~1e7 K, an order of magnitude above the Jupiter table's 1.2e6 K ceiling and 170x above baseline; needs new tables, a re-validated Saha ladder (the current one is tuned for `Z_bar` ~ 4–4.5), and real opacity — the Kramers shape already underflows past 1e5 K. That is a study, not an item. |

**Two questions this work generated that were not in the seven:**

1. **Is the escaping snowplow affordable?** P9's field window buys detachment by letting the
   collision-phase snowplow past the weakened exit field, 3–8x over. That is a collision-phase
   ablation calculation.
2. **What happens when sub-Alfvénic plasma is forced against a wall?** If a physical diverging
   nozzle takes over downstream, the plasma does work against field lines that want to curve back —
   currents into the wall, reaction forces on the last coil. Real MHD, not modelled here.

---

## Provenance

Every figure above was produced by running the code on 2026-09-03, not by recalling it.

| run | covers |
| --- | --- |
| `make analysis-nozzle-ledger` | P4, P5's gates, P10; the reflection baseline, signed `f_d`, `eps_b`, the power balance |
| `make analysis-nozzle-field` | P6, P8; the power-law fit, the winding sweep, realizability, the column-length lever |
| `make analysis-nozzle-detachment` | P3, P9; `M_A(z)`, `beta(z)`, the field window, liner ablation depth |
| `make analysis-nozzle-jet` | P2; the adiabatic conversion, `eta_jet` per leg, the adiabaticity check |
| `make analysis-expansion` | P5's radiative bracket (`data/results/cooling_history.csv`) |
| `cargo test --release -p euler2d --test merge_expansion -- --nocapture report` | P1; `alpha`, the spherical control, the `cos(theta)` histogram |

Artifacts: `data/results/nozzle_baselines.csv`, `nozzle_field_ripple.csv`,
`nozzle_detachment.csv`, `nozzle_jet.csv`.

**Errors made and caught while doing this work**, recorded because each would have shipped a wrong
number and because the tests that caught them are the reason to trust the rest: a transverse mean
folded into the drift energy (axisymmetric `u_r >= 0` does not cancel, but true 3D transverse
momentum is zero by symmetry); `Grid2D::run_to` integrating a *duration* rather than to an absolute
time, which mislabelled a diagnostic trace; a mass-conservation check that passed because a
staircased initial condition and a boundary loss cancelled; and a least-squares winding that reached
for counter-wound end coils and manufactured its own ripple.
