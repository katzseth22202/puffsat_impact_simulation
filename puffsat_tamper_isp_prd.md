# PuffSat Tamped Nozzle — Effective Isp Study

**Product requirements document**

**Status:** approved plan, no code written.
**Date:** 2026-08-12.
**Audience:** this document is written to be self-contained. A hydrocode specialist with
no prior exposure to this project should be able to read it end to end and understand what
is being modelled, why, what the code must do, and how we will know it is right. Prior
project documents are cited for provenance, not as prerequisites.
**Relationship to the rest of the repository:** the existing study
([`puffsat_impact_sim_design.md`](puffsat_impact_sim_design.md),
[`CONCLUSION.md`](CONCLUSION.md)) computes a per-collision momentum-transfer efficiency
`f` for a gas pulse striking a pusher plate. **This is a different study** with a
different deliverable (effective specific impulse), a different device (a tamped nozzle),
and a different regime. It reuses that study's kernels, tables, and validation discipline.
Where it *departs* from a decision recorded there, §12 says so explicitly.

---

## 1. Executive summary

A spacecraft is accelerated by shooting mass at it. Projectiles ("PuffSats") are launched
from elsewhere and meet the vehicle head-on at ~75 km/s. Each projectile flies through a
hole in the centre of the vehicle's pusher plate and buries itself in a slug of ice carried
just beyond it. The collision vaporises everything; the resulting fireball expands, and the
fraction of it that comes back and strikes the plate provides thrust.

Momentum conservation carries the merged blob *away* from the plate, so only the part of
the fireball that outruns its own recoil ever returns. Ballistically that is ~31%. A
**tamper** — a dense mass placed on the far side of the slug — is proposed to turn some of
the away-going half around.

**The question this study answers:** does the tamper pay for its own mass, and if so in
what configuration?

**Why it is not settled already.** Two defensible hand models bracket the answer at
**591–965 s of effective Isp against 984 s for no tamper at all** — i.e. the bracket
straddles "clearly worth it" and "actively harmful." The disagreement is entirely about
pressure coupling in a partially-ionised, optically-thick, self-similar expansion, which is
not an algebra problem.

**What makes this tractable:** §3 establishes a single thermodynamic ceiling that every
configuration can be scored against, converting a contested multiplier into a measured
**realization fraction**. The tamper must achieve **>62.9%** of that ceiling to beat doing
nothing. A *perfect* mirror achieves 61.7%. So the entire question is whether pressure
coupling — omitted from every existing model — closes a 1.2-percentage-point gap.

**Vehicle context:** 1000 t vehicle; a 200 kg reference pulse at 1–4 Hz gives 1.69 m/s of Δv
per pulse (0.17–0.69 g) over ~5400 pulses. Pulse mass and cadence are a **free trade at
fixed thrust**, not independently pinned (§4.1).

**The plate side turned out to be simpler than first scoped.** An ablating surface is
temperature-pinned, so the substrate cannot exceed the ablator's vaporisation temperature —
a steel plate equilibrates near **750–905 K** at 4 Hz and is comfortably fine with no
material escalation (§6.5.2). The only plate-side quantity that matters is therefore
**ablator mass**, which is charged in the Isp denominator and is currently uncertain by 27×.

**Two passes.** Pass 1 (Rungs 1–5) **ignores ablator mass entirely**, isolating the tamper
question — the point of the study — from the plate question and yielding an explicit *upper
bound* on Isp. Pass 2 folds in the measured ablator cost (Rung 6) last, revising Isp down
and shifting the optimal mass ratio up. **Uncertainty reduction precedes cost refinement**:
the RT bound is inside Pass 1's bracket and is its widest term, so it runs at Rung 2.

---

## 2. Architecture background

### 2.1 The propulsion concept

The vehicle does not carry its propellant's energy. Projectiles are accelerated by
infrastructure elsewhere in the solar system and aimed at the vehicle. On the outbound leg
the vehicle is on a prograde escape trajectory and the projectiles are retrograde, so they
meet **head-on** at a closing speed `w` of 74–81 km/s. This study uses `w = 75 km/s`.

Because the encounter is head-on, the projectile's momentum is *opposite* the desired
thrust. It is therefore a **debit**, and net thrust exists only if the merged material is
ejected backwards faster than the projectile arrived forwards. This is the origin of the
`−1` in every impulse law below.

### 2.2 The device

```
        thrust  ←── +z

   vehicle │  pusher plate        vacuum gap or         tamper
           │  (concave, hole      snow interlayer       (dense ice
           │   at the vertex)     (standoff s)           shell)
           │        ▓                                       ▓
           │        ▓        ○ slug + fireball              ▓
   ════════╪════════▓══════ ( ● ) ══════════════════════════▓══
           │        ▓        ○                              ▓
           │        ▓                                       ▓
                    ↑                                       ↑
              projectile enters                    fireball's backward
              through vertex hole,                 half is turned here
              travelling −z at 75 km/s
                                    standoff to plate d ──────►
```

1. The projectile passes through a hole at the plate's vertex without touching it.
2. It buries itself in the **slug** — carried mass held beyond the plate. The collision is
   completely inelastic and everything vaporises.
3. The merged blob's centre of mass recedes from the plate at `V = w/(1+k)` while the
   fireball expands at `u = w√k/(1+k)`, where `k` = slug mass / projectile mass. Since
   `u/V = √k`, only material with `cos θ > 1/√k` outruns the recoil and reaches the plate —
   a **capture fraction of `(1 − 1/√k)/2`**, which is 31.2% at `k = 7.06` and **exactly
   zero at `k ≤ 1`**.
4. The **tamper** sits beyond the slug and turns part of the away-going half around. It is
   never struck by plate-bound material.
5. Gas reaching the plate stagnates against it and re-expands, delivering axial impulse.

### 2.3 What "tamper" means here

A tamper in the inertial-confinement sense: it does not need to *survive*, only to remain
in the way for the microseconds of the event. **It is fully vaporised** — at these
energies it is heated ~26× past sublimation. Vaporisation is not a failure mode, because
the tamper works by **areal density `σ = m/A`**, which is conserved regardless of phase
(§6.3). It is expected to end up several times cooler than the fireball.

---

## 3. The objective function

### 3.1 Effective Isp and the mass that is charged

The deliverable is **effective specific impulse in the thrust direction**:

```
Isp_eff = J / (g₀ · m_charged)
```

`m_charged` is **all expended carried mass**. It is evaluated in two passes:

| | `m_charged` | status |
|---|---|---|
| **Pass 1** (Rungs 1–5) | slug + tamper + interlayer | **an upper bound on Isp**, always labelled as such |
| **Pass 2** (Rung 6) | the above **+ ablator** | the honest number |

Pass 1 deliberately excludes the ablator so that the tamper question — the point of the
study — is isolated from the plate question and can be answered without waiting on Rung 6.
Pass 1 also reproduces the denominator prior work used, so it is directly comparable to
existing figures. **No Pass-1 Isp may be quoted without the upper-bound label**, since the
ablator's share is currently unbounded between ~2% and ~60%.

**The projectile is not charged here, and the study reports a second metric instead of
pricing it.** The vehicle spends two distinct resources per unit of Δv: mass it *carries*
(charged by the rocket equation) and *projectiles*, which it does not carry but which
someone had to manufacture, launch, and aim. These are not convertible without a
program-level economic model that lies outside this study, so no combined figure of merit is
constructed. Both are reported as curves:

| metric | definition | meaning |
|---|---|---|
| **Effective Isp** | `J / (g₀ · m_charged)` | carried-mass economy — the rocket-equation currency |
| **Projectile economy** | `β · w`, i.e. `J / m_projectile` | infrastructure throughput per unit impulse |

`β` already appears throughout §3.4; naming it as the second metric is the only change.

**What stays out of scope is only the exchange rate** between projectile consumption and
delivered payload — that needs program economics this study does not have. The mass ratio
itself is a design variable this study *does* optimise; see §3.5.

### 3.2 The ceiling

Vehicle frame, `+z` = thrust direction. Impulse on the vehicle is the momentum brought in
minus the momentum carried out:

```
J = −m_i·w − P_ejecta
```

The ejecta's total energy is `½ m_i w²` (slug and tamper start at rest in the vehicle
frame), so by Cauchy–Schwarz `|P_ejecta| ≤ √(2 M_ej E)` with `M_ej = m_i(1+K)`. Define

```
K = k(1+τ) = total CARRIED mass per projectile kg
             (τ = tamper mass per slug kg)
```

Then, per projectile kg:

```
J_max   = w·(√(1+K) − 1)
v_e,max = J_max/K = w·(√(1+K) − 1)/K
```

Equality requires **every ejecta element to end up moving −z at one common speed**, which
in turn requires the plate to catch and reverse all plate-bound material. This reproduces
the project's prior analytic results exactly (`β_ideal = √(1+k) − 1`; bare-plate optimum at
`k* = 7.057`; 1014 s at `w̄ = 77.28 km/s`).

### 3.3 Two consequences that define the study

**(a) At the ceiling, the tamper is neutral.** `v_e,max` depends only on `K`, not on how
`K` splits between slug and tamper. **The tamper can therefore never be justified as a
momentum multiplier** — only as a *realizability* device. The deliverable is a
**realization fraction**: what share of the ceiling a configuration achieves.

**(b) The tamper's backward recoil is credited, not lost.** Its `−z` momentum is the exact
accounting mirror of the `+z` impulse delivered to the plate. This reverses the prior
framing, in which tamper recoil was treated as the mechanism that kills the concept.

**Design consequence:** what is actually lost is **entropy** — energy thermalised in the
tamper, which must then be re-expanded a second time, isotropically, at a second efficiency
cost. Entropy production rises steeply with shock strength. **The ideal tamper is a
maximally isentropic piston, not a cold specular mirror.** Instrument the tamper's
entropy/thermal budget and its centre-of-mass momentum, not its reflected-speed ratio.

### 3.4 The number to beat

At `w = 75 km/s`, `k = 7.06`. `β` is axial impulse per projectile kg in units of `w`;
`β_bare` and `β_tamp` are the ballistic closed forms from prior work, re-verified here
(`β_bare(7.06) = 0.90871`, `β_bare(14.12) = 1.68353`):

| configuration | `K` | `β` | % of ceiling | Isp |
|---|---|---|---|---|
| **bare plate, k = 7.06** | 7.06 | 0.9087 | **49.4%** | **984 s** |
| bare plate, k = 14.12 — same mass spent as *slug* | 14.12 | 1.6835 | 58.3% | 912 s |
| **perfect-mirror tamper, τ = 1** | 14.12 | 1.7825 | **61.7%** | **965 s** |
| *break-even against the bare plate at k = 7.06* | 14.12 | 1.817 | **62.9%** | 984 s |
| ceiling at `K` = 7.06 / 14.12 | — | 1.839 / 2.889 | 100% | 1992 / 1565 s |

Read this row by row. A tamper at `τ = 1` **beats** spending the same mass as extra slug
(61.7% vs 58.3%), but **loses** to not spending the mass at all (61.7% vs the 62.9% needed)
— *and that is with a perfect mirror*. Any finite-mass mirror does worse. Adding the
free-plate elastic lower bound gives the honest bracket of **591–965 s against 984 s**.

**Therefore the tamper pays only if pressure coupling lifts it above 62.9%.** Every prior
model of this device is ballistic (straight-line gas elements, `ρv²` only, no `P` term), and
at `k ≈ 7` the fireball's sound speed and its recoil velocity are equal to within 5% — the
worst possible regime in which to omit pressure. That is the gap this study exists to close.

---

### 3.5 The mass ratio is a design variable, and the Isp optimum is not where it settles

`K` — carried mass per projectile kg, equivalently the split of the fixed 200 kg pulse
between projectile and carried mass (`m_i = 200/(1+K)`) — is a primary design output of this
study, not an inherited input.

**There is a genuine physics optimum.** `K* ≈ 7.06` for a bare plate is not an economic
result: it is the competition between capture fraction, which rises with `K`, and exhaust
velocity, which falls with it. Pure geometry, entirely in scope.

**But it is very flat, and the quantity it trades against is steep.** Three exact relations
govern the trade — note that the first two contain **no pulse mass and no cadence**, so they
hold under any pulse-size/rate schedule:

```
energy per unit impulse         E/J = w/(2β)     [J per N·s]
projectiles per unit impulse         = 1/(β·w)
effective Isp                        ∝ β/K
```

**Plate heat load and projectile consumption are the same function** — both scale as `1/β`
— so they are not two currencies needing an exchange rate. They are one quantity, improving
steeply with `K` while Isp degrades gently:

| `K` | `β_bare` | Isp | vs peak | **energy per unit impulse** | vs peak |
|---|---|---|---|---|---|
| 6 | 0.768 | 979 s | −0.6% | 48.8 kJ/N·s | +18% |
| **7.06** | 0.909 | **984 s** | **peak** | 41.3 kJ/N·s | — |
| 8 | 1.027 | 981 s | −0.3% | 36.5 kJ/N·s | −11% |
| 10 | 1.260 | 963 s | −2.1% | 29.8 kJ/N·s | −28% |
| 14.12 | 1.684 | 912 s | −7.4% | 22.3 kJ/N·s | **−46%** |
| 16 | 1.858 | 888 s | −9.8% | 20.2 kJ/N·s | −51% |
| 32 | 3.084 | 737 s | −25% | 12.2 kJ/N·s | −71% |

**What sets the design point.** Plate *temperature* does not — it is self-limiting at the
ablator's vaporisation point regardless of `K` or cadence (§6.5.2). What does is the
**ablator's share of the Isp denominator**, since ablation is sub-linear in fluence and so
falls with `β`. The denominator is therefore

```
m_charged / J  ∝  K/β  +  (ablator term)/β^α        with α < 1
```

The first term is minimised at `K* ≈ 7.06`; the second decreases monotonically with `K`. So
their sum has **a genuine interior minimum above `K* = 7.06`**, located by how large the
ablator's share turns out to be. This is entirely internal to Isp — no external constraint,
no economic weighting.

**Consequence for the two passes.** On **Pass 1** the ablator term is zero by construction,
so the optimum sits exactly at the ceiling optimum `K* ≈ 7.06`. On **Pass 2** it shifts
upward, and by how much is one of Rung 6's outputs. The `K` grid must therefore span well
past 7 — at least 6–32 — rather than clustering around the Pass-1 answer.

**Consequence — the tamper's case gets harder.** The tamper and extra slug are two routes to
the same end: both raise `β`. §3.4 shows extra slug reaches 58.3% of ceiling against the
tamper's 61.7%, and extra slug additionally carries no RT risk, no interlayer, no porosity
model, and no assembly complexity. **The tamper must beat a rival that is strictly simpler.**

---

## 4. Reference case

`w = 75 km/s`, `k = 7.06`, total per-pulse mass 200 kg.

| quantity | τ = 0 | τ = 1 |
|---|---|---|
| projectile `m_i` | 24.81 kg | 13.23 kg |
| slug | 175.2 kg | 93.4 kg |
| tamper | — | 93.4 kg |
| **pulse energy `½ m_i w²`** | **69.8 GJ** | **37.2 GJ** |
| blob CM recoil `V = w/(1+k)` | 9.31 km/s | 9.31 km/s |
| fireball expansion `u = w√k/(1+k)` | 24.72 km/s | 24.72 km/s |

**Blob specific internal energy** `e = ½w²·k/(1+k)² = 305.7 MJ/kg`, i.e. **57.1 eV per
H₂O molecule**. After atomisation (9.5 eV) and first ionisation of 2H + O (~40.8 eV), only
~7 eV of genuinely thermal energy remains, spread over ~6–7 particles: **T ≈ 1.2 eV
(≈14 kK)**, not 57 eV. Sound speed `c_s = √(γ(γ−1)e) ≈ 9.8 km/s` at `γ_eff = 1.25`.

**This is why the EOS is load-bearing rather than a detail.** The plume temperature is set
almost entirely by how much energy the EOS charges for dissociation and ionisation, and it
lands in the 0.5–2 eV window where the ionisation fraction changes fastest. An ideal-gamma
EOS would over-predict `T` by more than an order of magnitude and, through it, the pressure
that the whole tamper mechanism depends on.

### 4.1 Vehicle and mission context

**Vehicle mass is pinned at 1000 t** (the plate is 5% of it). **Pulse mass and cadence are
not independently pinned** — they are a free trade at fixed thrust. A 200 kg pulse is the
reference point, not a constraint:

| quantity (200 kg reference pulse) | value | derivation |
|---|---|---|
| impulse per pulse (bare, `K` = 7.06, τ = 0) | 1.69×10⁶ N·s | `β_bare · m_i · w` |
| **Δv per pulse** | **1.69 m/s** | modest, as intended |
| **vehicle acceleration** | **0.17 g at 1 Hz / 0.69 g at 4 Hz** | |
| carried mass flow | 200 kg/s / **800 kg/s** | |
| departure burn (Δv ≈ 7.06 km/s at Isp 984 s) | ~1350 s, ~5400 pulses at 4 Hz | mass ratio 2.08 |
| plate recoil per pulse (50 t plate) | 34 m/s | 0.17% of gas speed |

Two cross-checks fall out. The 800 kg/s carried-mass flow at 4 Hz matches this project's
independently-derived requirement (8.85 kg/s for a 10 t vehicle, i.e. 885 kg/s at 1000 t) to
within 10%. And the 34 m/s plate recoil against ~20 km/s gas confirms the **rigid-wall
assumption with a factor of ~600 in hand**.

**The pulse-mass / cadence trade is genuinely free, and average heat load is invariant under
it.** Since `E/J = w/(2β)` contains neither pulse mass nor cadence (§3.5), average thermal
power is `Thrust · w/(2β)` — **doubling pulse mass and halving cadence changes it not at
all.** Thrust, and therefore acceleration and gravity loss, are likewise held. What the
trade *does* buy:

- **Less total ablator**, because ablation is sub-linear in fluence — concentrating the same
  energy into fewer, larger pulses ablates less overall, and ablator mass is charged (§3.1).
- **Better optical depth**: `Σ` at the plate scales with pulse mass, moving further from the
  `τ ~ 1` regime where flux-limited diffusion is weakest and where this project has
  previously suffered a 2000× opacity error (§5.3).
- **Scale-invariance of the tamper physics**: both `Θ` (§6.1) and the RT mixed fraction
  (§6.7) are ratios and do not move, so the trade cannot rescue or break the tamper.
- Cost: higher per-pulse fluence means deeper ablation per pulse, so the layer must be
  thicker — though the longer interval gives more time to lay it down.

**Cadence is set by gravity losses alone.** Prior finite-burn integration on this project
gives departure gravity losses of 2.2% at 1 g, 5.1% at 0.5 g, 9.3% at 0.25 g, so higher
acceleration is preferred — and at fixed thrust that is a statement about pulse mass × rate,
not about cadence alone. Plate thermal rejection does **not** compete for this knob (§6.5.2).

One consequence stays out of scope but is no longer deferrable: the plate takes 34 m/s per
pulse and must return before the next one, implying metres of shock-absorber stroke.

### 4.2 The configuration under test

Two configurations were scoped as corners of a single `(slug density, standoff)` design
space. **Arm B is now provisionally foreclosed** (§6.7), so Arm D is the design and Arm B
survives only as a control case.

- **Arm D — filled (the design).** Solid ice slug → **snow or slush interlayer** occupying
  the standoff → ice tamper, as one assembled body. The interlayer transmits pressure
  continuously, so the tamper is loaded gently and early — the isentropic-piston ideal; it
  needs no support structure; and it is itself reaction mass.
- **Arm B — vacuum (control only).** Solid ice slug → **vacuum gap** → ice tamper on a
  consumable spacer. Two independent arguments close it: the plume goes ballistic before it
  is turned, so turning it is a re-thermalising ram (the largest entropy source available,
  against the §3.3 design objective); and Rayleigh–Taylor **disrupts its tamper completely**
  before the tamper finishes its job (§6.7). It is retained as a control — a code that does
  not reproduce that disruption is not to be trusted on Arm D either.

---

## 5. Physical regime

### 5.1 State space the code must cover

| axis | range | note |
|---|---|---|
| density | **10⁻⁴ – 10³ kg/m³** | solid ice at 917 down to ~4×10⁻³ at the plate |
| temperature | **300 K – ~10⁵ K** | 14 kK in the fireball; **50–80 kK at plate stagnation** |
| specific internal energy | up to ~300 MJ/kg | 57 eV/molecule |
| composition | H₂O through dissociation and first ionisation | multi-stage O ionisation at stagnation |

**Both ends exceed the project's baseline table** (`data/tables/water.json`: ρ 0.01–1198
kg/m³ over 61 points, T 300–60,000 K over 60 points). The extended-grid table built for a
prior 69 km/s scenario (`water_jupiter.json`: T to 1.2×10⁶ K, ρ to 10⁻⁴, carrying real
TOPS/OPLIB gray opacities) **does** cover it and is the intended baseline here.

### 5.2 Timescales

| event | scale |
|---|---|
| projectile transit of the slug | ~10 µs |
| slug disassembly `r_slug/u` | **11.7 µs** (ice slug) / **27.6 µs** (snow slug) |
| shock transit of the tamper | 1.4–2.9 µs (see §6.1) |
| tamper confinement (recoil-limited) | ~45 µs |
| lateral communication across the fireball `r/c_s` | ~69 µs |
| **arrival window at the plate** (velocity dispersion, d = 10 m) | **~750 µs** |
| stagnation-layer radial relief `R/c_s` | ~1 ms |
| total simulated time | **~2–3 ms** |

### 5.3 Optical depth

Areal density of the gas slab at the plate is `Σ ≈ 0.063 kg/m²`. With κ for warm dense
water uncertain across ~two decades (10–1000 m²/kg), `τ = κΣ` spans **0.63 to 63** — it
**straddles `τ ~ 1`**, which is exactly where flux-limited diffusion is weakest. Prior
experience on this project is directly relevant: at 69 km/s an interim Kramers opacity ran
~2000× low at stagnation and falsely predicted `τ ~ 1`, moving `e_eff` from 0.42 to 0.65
once real opacities were used. **Real tabulated opacity is a requirement here, not a
refinement.**

---

## 6. Findings that set the requirements

All derived from scratch; estimates are flagged as such.

### 6.1 The transit ratio — the tamper's binding timing constraint

A tamper functions only if it can communicate its inertia across its own thickness before
the event ends. Let `Θ = t_disassembly / t_shock-transit`. With
`t_dis = r_slug/u`, tamper thickness `h = m_t/(2πs²ρ_t)`, and `m_s = (4/3)πr³ρ_s`:

```
Θ = (3/2) · (ρ_tamper/ρ_slug) · (U_shock/u) · (s/r_slug)² / τ
```

`U_shock ≈ 12 km/s` (ice Hugoniot at ~70 GPa — **estimated**, worth confirming). At τ = 1:

| configuration | slug radius | tamper thickness | transit | disassembly | **Θ** |
|---|---|---|---|---|---|
| snow slug (70 kg/m³), contact ice tamper | 0.683 m | 3.5 cm | 2.9 µs | 27.6 µs | **9.5** ✓ |
| ice slug (917), contact ice tamper | 0.290 m | 19.3 cm | 16.1 µs | 11.7 µs | **0.7** ✗ |
| ice slug, tamper at 1 m standoff | 0.290 m | 1.6 cm | 1.4 µs | 11.7 µs | **8.7** ✓ |

**The tamper/slug density ratio is the dominant design parameter — for a geometric timing
reason, not an acoustic one.** It sets tamper thickness at fixed mass. There is no hard cap
on `τ` near 1: with a snow slug, `Θ > 3` holds out to `τ ≈ 3`.

### 6.2 Projectile penetration — a hard constraint on projectile geometry

Momentum-conserving snowplow along the axis; retained speed = `σ_proj/(σ_proj + σ_target)`.
A 13.23 kg ice projectile has areal density **421 kg/m²** as a 0.2 m-diameter rod, or
**95.5 kg/m²** as a 0.42 m-diameter disk. Target column densities: snow slug 95.6 kg/m²,
solid ice slug 531.7 kg/m².

| projectile | target | retained speed | **energy deposited** |
|---|---|---|---|
| 0.2 m rod | snow slug | 0.815 | **34%** — punches clean through |
| 0.42 m disk | snow slug | 0.500 | 75% |
| 0.2 m rod | solid ice slug | 0.442 | **80%** |

**Requirement:** the projectile's areal density must be ≲ the slug's column density, or the
slug becomes a spectator and the tamper becomes the fireball. The vertex hole is *not* the
binding limit — a 0.5 m hole in a 15 m-radius plate is 0.03% of its area.

**Projectile geometry is therefore a design variable inside this study, not an external
input.** The swept parameters are areal density `σ_proj` (equivalently diameter at fixed
mass), aspect ratio, and bulk density, subject to three constraints: it must pass the vertex
hole, it must survive launch and interplanetary transit as a coherent body, and it must be
deliverable to the aim tolerance. The table above already shows a genuine interior optimum
— too slender and it punches through, too broad and it fails the hole and the delivery
constraints.

**Consequence for the model:** energy deposition becomes an *output* of the near-field
simulation rather than a prescribed initial condition. That is a real scope increase —
resolving a 75 km/s hypervelocity penetration with fragmentation, jetting, and phase change
on a ~10 µs timescale is the domain of dedicated impact codes, and it carries its own
validation burden (§8). It does **not** drive resolution: 1–2 cm cells give 20–40 across a
0.42 m projectile, far coarser than the tamper's 1.75 mm requirement.

*Model caveat: the snowplow screen above ignores lateral spreading of projectile debris and
treats a disintegrating ice rod as a coherent penetrator. It is an order-of-magnitude
estimate, and replacing it is exactly what the near-field simulation is for.*

### 6.3 The tamper's operating window

Three clocks govern the tamper, and **temperature is in none of them**. At τ = 1, s = 1 m,
`σ_t = 14.9 kg/m²`:

| clock | mechanism | value |
|---|---|---|
| slug disassembly | the event it must act on ends | 11.7 µs |
| **recoil** `σ_t·Δv/P` | it comoves with the plume; no further momentum transfers | **~45 µs** |
| lateral spread `R/c_s` | area grows, `σ` falls | ~200 µs |

It completes its work ~4× before it disperses. The tamper–plume collision deposits
~72 MJ/kg — **26× past sublimation, and ~4× cooler than the fireball** (≈0.5 eV vs 1.2 eV).

Cross-check: `σ_tamper/σ_plume = 1.75` here, giving a 1-D elastic free-plate
reflected/incident ratio of **0.274**, which reproduces the prior analytic lower bound at
τ = 1 to two figures.

### 6.4 The plate sees a quasi-steady plenum, not a bounce

Plate-bound gas spans ~8–20 km/s, so at `d = 10 m` the **arrival window is ~750 µs**, while
radial relief of the stagnation layer takes `R/c_s ≈ 1 ms`. **These are comparable**, so the
layer is continuously fed and cannot relieve during the pulse. The prior study's regime is
the exact opposite — a ~µs bounce with fast relief. Consequences:

- **Survivability inverts.** Peak plate pressure is **≈ 2 MPa** against a 400 MPa
  facesheet baseline — roughly **200× margin**. (Relief sets layer thickness, not stagnation
  pressure, so this holds in the plenum regime.) **Survivability in this architecture is
  entirely thermal.**
- Whole-plate rigidity passes trivially: a 30 m plate's first flexural mode is 10–100 ms
  against a 750 µs pulse, so the rigid-wall assumption underpinning the impulse calculation
  is safe.
- **The existing kernels' plate-side conventions are both wrong here**: an initial condition
  that places a finite slug in a box, and an integration window that stops when the wall
  force decays to `10⁻³` of peak. Both need rework for a sustained feed.

### 6.5 Ablation is the plate's only thermal sink, and its mass is a forced, unbounded cost

**Scaling.** Fluence at the plate `∝ 1/d²` (independent of `R`); captured mass `∝ R²/d²`.
So **capture-per-unit-fluence `∝ R²` and nothing else — plate radius is the only lever that
breaks the standoff conflict.** In the linear regime the ablator is a fixed fractional tax,
invariant in both `R` and `d`. With vapour shielding, ablation is sub-linear in fluence
(`σ_abl ∝ Φ^α`, `α < 1`), so the tax scales as `d^{2(1−α)}` — **it grows with standoff.**
Standoff therefore trades a peak pressure that is not binding for an ablator mass that is
charged in the denominator.

**The plate has no thermal sink except ablation.** Over a 750 µs residence, against an
incident fluence of **≈12.6 MJ/m²**:

| sink | capacity | share |
|---|---|---|
| soak into steel (`√(4αt)` ≈ 173 µm → 1.36 kg/m²) | ~0.95 MJ/m² | 8% |
| re-radiation, steel at 1700 K | 3.6×10⁻⁴ MJ/m² | 0.003% |
| re-radiation, SiC at 2800 K | 2.6×10⁻³ MJ/m² | 0.021% |

Re-radiation is four orders too slow to matter and soak covers at most ~8%. **Ablation is
not a lever we choose; it is the only sink.** The only genuine design choices are ablator
material (`Q*`) and where the mass is placed (taper, §6.6).

*Note the soak row is a capacity, not a prediction — it assumes the full penetration depth
reaches melting, which happens only if the vapour curtain fails. §6.5.2 shows the physical
value sits far below it, which is why the plate temperature is not a constraint.*

The entire ablator mass is therefore set by one thing: **how much of the incident fluence
the vapour curtain blocks before it reaches the wall** — a self-consistent balance that
resists analytic estimation:

| basis | depth/pulse | mass at R = 15 m | Isp tax on 200 kg |
|---|---|---|---|
| prior measured 16 km/s ablation fraction (3.7–8.9% of pulse mass), rescaled to this Σ | 7–16 µm | 4.5–10 kg | **2–5%** |
| same, with naive `v⁸` scaling to ~20 km/s arrival | 40–190 µm | 25–121 kg | **13–60%** |

**A 27× spread**, and the study's largest single unknown. It is excluded from Pass 1 by
construction (§3.1) precisely so that it cannot contaminate the tamper verdict; Rung 6
collapses it for Pass 2.

**Why excluding the ablator is conservative for the tamper verdict.** An earlier draft of
this document claimed the opposite — that the tamper "induces roughly proportional extra
ablation," partly consuming its own gain. That is true *per pulse* and wrong *per unit
impulse*, which is what Isp measures. Since `E/J = w/(2β)` (§3.5), raising `β` lowers the
energy delivered per unit impulse and therefore the ablation per unit impulse. **The tamper
reduces ablator mass per unit impulse**, so Pass 1's exclusion of the ablator understates the
tamper rather than flattering it.

**And the ablator degrades the answer without deciding it.** Even the pessimistic branch
(121 kg against a 200 kg pulse) gives a denominator ratio of 200/321 = 0.62, taking a 984 s
Pass-1 figure to ~610 s — still ~1.6× methalox. The ablator therefore moves the absolute
number materially but does not change whether the architecture is worth pursuing, and it
nearly cancels in the *comparative* question of whether the tamper pays. That is why it is
last (Rung 6) rather than first.

#### 6.5.1 At 1–4 Hz the ablator becomes a mass-flow problem, not a coatings problem

The 27× uncertainty above translates directly into two qualitatively different vehicles:

| branch | per pulse | **at 4 Hz** | vs the 800 kg/s propellant flow | replenishment per pulse |
|---|---|---|---|---|
| optimistic | 4.5–10 kg | **18–40 kg/s** | 2–5% | 7–16 µm in <250 ms |
| pessimistic | 25–121 kg | **100–484 kg/s** | **13–60%** | 40–190 µm in <250 ms |

On the pessimistic branch the ablator is a **larger consumable stream than anything else on
the vehicle except the propellant itself**, and the surface-renewal system must lay down a
~0.2 mm film over 707 m² between pulses. Since the ablator is charged in the denominator
(§3.1), this is what Rung 6 resolves.

**Sub-linearity makes pulse size a lever on this.** Because ablation goes as `Φ^α` with
`α < 1`, concentrating the same total energy into fewer, larger pulses ablates *less* in
total. Pulse mass and cadence trade freely at fixed thrust (§4.1), so **larger pulses at
proportionally lower cadence is a direct reduction in ablator mass** — one of the few levers
that improves the denominator without touching the physics of the tamper.

#### 6.5.2 Inter-pulse balance — the plate is self-limiting, and steel is fine

Ablation is mass-transfer cooling: its enthalpy leaves with the vapour, so the plate must
reject only what **soaks in**. The decisive point is that **an ablating surface is a
temperature-pinned boundary.** While it is ablating it sits at the ablator's vaporisation
temperature `T_abl`; the substrate beneath can never exceed that, because as the bulk
approaches `T_abl` the driving ΔT — and with it the conducted flux — goes to zero. The plate
temperature is therefore bounded by the ablator, not by the plume.

Solving the steady balance for a steel plate — conducted in per pulse
`≈ k·ΔT·√(t/πα)` (k = 45 W/m/K, α = 1.2×10⁻⁵ m²/s, t = 750 µs) against `2σT⁴` radiated from
both faces, at 4 Hz:

| ablator surface temperature `T_abl` | **plate equilibrium** |
|---|---|
| 800 K | **753 K** |
| 1000 K | **905 K** |

Steel is usable to ~1000–1200 K. **It clears this comfortably, at 4 Hz, with no material
escalation and no active cooling.** This is why Orion's steel plate worked, and it is
cadence-independent in the same way §4.1's heat load is: raising cadence raises both the
conducted flux and the equilibrium radiating temperature, and the balance simply moves along
the same curve.

**An earlier draft of this document got this wrong** and should not be relied on: it took the
*capacity* of the thermal penetration depth (0.95 MJ/m², reached only if the vapour curtain
fails completely) as the actual soak, concluded the plate was ~2× short at 4 Hz, and built a
refractory-material ladder on it. The plausible soak range is 0.15%–7.5% of incident fluence
— a 50× spread — and the self-limiting argument above shows the physical answer sits near the
bottom of it.

**What actually binds is ablator burn-through, not plate temperature.** The entire argument
above assumes the ablator layer is present throughout the pulse. If it burns through partway,
the wall sees the plume directly for the remainder — at 50–80 kK, with no pinning — and the
failure is abrupt rather than graceful. So the plate-side deliverable reduces to two
questions, both answered by the same Rung 6 computation: **how much ablator is consumed per
pulse** (it is charged in the denominator, §3.1) and **what burn-through margin remains**.

#### 6.5.3 Regenerative cooling — available, but not needed

**§6.5.2 removes the motivation for this**, since the plate self-limits near 750–905 K. It is
retained as a documented option in case burn-through margin turns out to demand a colder
substrate, or a later scenario pushes `T_abl` higher.

Cooling cannot help *within* a pulse in any case: the thermal wave reaches only ~173 µm in
750 µs, so no coolant channel is in contact with the load. Ablation is the intra-pulse sink
for any material. Between pulses, water is already the propellant, so **200 kg/pulse of it
flows regardless — the cooling would be mass-free.**

| water state reached | Δh | per pulse (200 kg) | share of the 672 MJ soak |
|---|---|---|---|
| ice warmed to just under melt | 0.33 MJ/kg | 66 MJ | 10% |
| melted to liquid at 273 K | 0.66 MJ/kg | 132 MJ | 20% |
| liquid at 373 K | 1.08 MJ/kg | **216 MJ** | **32%** |
| saturated steam | 3.34 MJ/kg | 668 MJ | 99% |

**Phase change is one-way without a radiator**, so anything boiled must reject its latent
heat again before it can become a slug — and a steam slug is ~300× too dilute to stop the
projectile (§6.2). Boiled coolant therefore cannot double as propellant. Melting is
admissible: liquid water at 1000 kg/m³ against ice's 917 means a **liquid or slush slug has
ample column density.** So the free regenerative budget caps at **~32% of the upper-bound
soak.**

**If the true soak is ≲1/3 of the upper bound, propellant regenerative cooling closes the
inter-pulse gap entirely and at zero mass cost.** Above that, the options are vented water
(charged in the denominator exactly like the ablator), more radiator area, or a hotter
plate.

**The structural objection does not bind in this regime.** This project forbids voids
directly behind the hot face, because a void is a free surface that reflects the compressive
pulse as full tension and causes prompt spall. That rule was derived at 400 MPa–2 GPa
facesheet loads; **peak pressure here is ~2 MPa**, three orders of margin. Coolant channels
are structurally admissible, and they need to sit within ~3 mm of the surface (the 250 ms
diffusion length in steel) — a placement the original rule would have forbidden.

#### 6.5.4 Plate material: steel plus a thin oil ablator, and the question largely dissolves

Every candidate criterion turns out to be non-binding:

| criterion | status |
|---|---|
| shock (2 MPa gas, ~3 MPa inertial at 4600 g) | non-binding — ~200× margin for any structural material |
| spall | non-binding at 2 MPa; the rule that forbade voids was derived at 400 MPa–2 GPa |
| bending / whole-plate rigidity | passes trivially (first mode 10–100 ms vs a 750 µs pulse) |
| **steady-state temperature** | **bounded by `T_abl` ≈ 750–905 K — steel clears it (§6.5.2)** |
| atomic-oxygen attack | the renewed oil layer, not the substrate, meets the plume |

**Baseline: steel structure + a thin renewed oil ablator**, which is the Orion configuration
and is vindicated by the self-limiting argument rather than merely inherited from it.

A **thin SiC or ceramic front layer** is retained as a cheap hedge. It costs little, and it
buys margin in the one failure mode that is abrupt: if the ablator burns through mid-pulse,
the substrate briefly meets 50–80 kK plume with no temperature pinning, and a refractory face
degrades where steel would melt. It is insurance against burn-through, not a thermal
requirement.

**Escalation path, if burn-through margin proves thin:** Nb alloy (~1600 K), SiC/C-SiC or
Mo/TZM (~1900 K), or carbon-carbon (~2200–2500 K). Carbon-carbon is more available here than
in the prior study — that study rejected bare C-C because it burns in atomic oxygen, but with
oil renewed every pulse it is never bare during exposure and between pulses there is no gas
at all. At 1800 kg/m³ a 50 t plate is **14 mm thick at R = 25 m against steel's 3.2 mm**.
None of this is invoked unless Rung 6 shows burn-through margin is inadequate.

**On the Orion precedent.** The claim that a thin ablative oil protects a steel plate is
mechanistically correct and is the vapour-shielding physics this project already models: the
ablated vapour is opaque and absorbs the incoming radiation, so the substrate never sees the
full flux, largely independent of what the substrate is. Reproducing Orion's
impulse-per-pulse and ablation-per-pulse is already this project's designated keystone
validation. **What Orion validates is intra-pulse survival — the half where steel is fine.**
It does not address steady-state rejection across ~5400 pulses at 1–4 Hz, which is the half
that is marginal. *Citation caution: this project records that its Orion references are
secondary and the primary source was never read (firewalled), so any number leaned on here
must be checked against the originals.*

### 6.6 Plate shape

**The prior study's foreclosure of a deep dish is conditional on plane-wave incidence and
does not transfer.** A paraboloid focuses parallel→point and collimates point→parallel; it
cannot do both. A plane-wave cloud striking a dish gets its rebound *focused* into a hot
spot in strongly-radiating, optically-thick gas — the reason the deep dish was rejected. A
**point source at the focus** produces the opposite: a collimated beam and no focus at all.

*The fireball genuinely is a point source.* For an instantaneous expansion every element's
trajectory is `(V_cm + u·r̂)·t`, so all rays trace back to a **fixed** origin — CM recession
skews the angular distribution but does not move the apparent source. The blur is the finite
disassembly time, `u·t_dis ≈ 0.6 m`, against a focal length of 6–10 m: about 0.04 rad.

The focus-matched optimum is `d/D = R/(8d)` → **0.19** at (R = 15 m, d = 10 m), **0.31** at
(15 m, 6 m) — inside the previously foreclosed band.

**But the prize is bounded and the parabola carries a mass penalty:**

| | value |
|---|---|
| specular upper bound, parabola/flat, mass-weighted `(1+cosθ)/(2cosθ)` | 1.09 / 1.19 / **1.23** at `R/d` = 1 / 2 / 2.5 |
| prior *measured* concave lift at plane-wave incidence | `eta_capture` 0.915 → 0.977 → 0.994, **+9%** |
| paraboloid surface area vs flat disk (`F = d`, `R = 2.5d`) | **+32%** |

A ≤23% impulse gain against a ~32% area penalty in ablator and structure — with sub-linear
ablation making extra area *worse*, and more normal rim incidence collecting more flux.
**Under a charged-ablator denominator, flat may win.** The prior measurement is also direct
evidence that stagnation blunts shape by roughly an order of magnitude relative to ray
optics: the gas does not reflect, it stagnates into a subsonic plenum and re-expands, and
shape acts through *confinement* (raising the rim so the layer must relieve axially), not
through reflection.

**Taper is two separate levers, both closed-form cold-path calculations** requiring no
additional kernel runs:

- **Ablator taper.** Flat-plate flux `∝ d/(d²+r²)^{3/2}` falls 3× at `r = d` and 11× at
  `r = 2d`. Matching thickness to local fluence costs `2d²(1−cos θ_max)/R²` of the uniform
  mass: **0.59 / 0.28 / 0.20** at `R/d` = 1 / 2 / 2.5 — a **41–80% saving** on a term that
  may be 60% of the mass budget. Probably the largest single Isp lever on the plate side.
- **Structural taper.** The same profile applies to the impulse distribution; at a fixed
  plate-mass ceiling it buys up to **~2.2× the radius**. *Upper bound only* — real plates
  have a minimum gauge and bending loads do not track local pressure.

### 6.7 Rayleigh–Taylor limits what can be shaped, and threatens the mechanism

The tamper is a light fireball pushing a dense shell: Rayleigh–Taylor unstable. With
`a = P/σ ≈ 1.9×10⁸ m/s²` and `γ = √(A·k·a)`, over the 27.6 µs confinement window:

| feature wavelength | growth |
|---|---|
| 10 cm | **×20 — destroyed** |
| 1 m | ×2.6 — marginal |

The tamper is only ~0.68 m in radius, so **only gross shape survives; parabolic figuring of
the tamper is not maintainable.** Hence the tamper's design variable is angular coverage,
not curvature.

More seriously: **mixing decides whether the tamper remains a coherent inertial piston or
is entrained into the plume as ordinary reaction mass** — and that is precisely the
difference between 61.7% and 58.3% of ceiling, i.e. the whole result. Mixing is desirable at
the fireball/interlayer interface (it makes the interlayer a continuous pressure-bearing
medium) and fatal at the plume/tamper interface. **The two interfaces want opposite
outcomes.**

**An analytic bound, and it forecloses Arm B.** Using the self-similar mix growth
`h ≈ α·A·a·t²` with `α = 0.02–0.05` and Atwood number ≈ 1:

| | acceleration `a = P/σ` | window | bubble depth `h` | tamper thickness | **mixed fraction** |
|---|---|---|---|---|---|
| **Arm D** (contact, snow slug) | 1.9×10⁸ m/s² | 27.6 µs | 2.9–7.3 mm | 34.8 mm | **8–21%** (16–63% total mix width) |
| **Arm B** (vacuum standoff, 1 m) | 3.9×10⁸ m/s² | 45 µs | 16–40 mm | 16.2 mm | **100–246% — fully disrupted** |

**Arm B's tamper is shredded before it finishes its job**, across the whole plausible `α`
range — a thin sheet rammed by a fast plume is torn apart, while a thicker shell pressed
gently holds together. This is an *independent* mechanism reaching the same verdict as the
entropy argument in §3.3, and together they foreclose Arm B (§4.2).

*Bound, not prediction: the self-similar law assumes sustained acceleration and a broad
initial perturbation spectrum, whereas the acceleration here decays and the initial spectrum
is a manufacturing property. But a 2.5× overshoot does not survive a factor-of-two
correction.*

**For Arm D the bound is not decisive, and that is the problem.** A 16–63% mixed fraction
moves the realization fraction roughly between "piston" (~62%) and "just extra slug"
(~58%) — **wider than the 62.9% margin being decided.** An axisymmetric code cannot narrow
it. See the uncertainty budget (§13.1), where this is the top-ranked contributor to the
Isp bracket after the ablator.

### 6.8 Checked and dismissed

- **Pulse-to-pulse interference.** Exhaust leaves at ~20 km/s; the next projectile closes at
  75 km/s. After a ~20 s cadence the exhaust column density is ~0.3 kg/m² against a
  projectile areal density of 95–421 kg/m². Three orders down — not a constraint, even with
  a collimated rebound.
- **Vertex hole leakage.** 0.03% of plate area. Aiming tolerance is a separate problem and
  is out of scope here.

---

## 7. Modelling requirements

### 7.1 Staging is mandatory

| | scale | implied resolution |
|---|---|---|
| tamper thickness | 3.5 cm | 1.75 mm (20 cells) |
| standoff + plate radius | up to ~25 m + ~40 m | ~14,300 cells per direction |

A uniform 2-D grid spanning both is **~3×10⁸ cells** — infeasible once, let alone across a
sweep. The pipeline is therefore staged, with a control-surface hand-off:

1. **Near field** — resolve the tamper; run until the flow is ballistic. Output: mass /
   velocity / angle distribution on a control sphere, plus the entropy budget.
2. **Transport** — control sphere → plate. **Not free-streaming:** terminal Mach is ≈2.5,
   so pressure still steers the flow and a straight-line ballistic map is too crude. A
   coarse 2-D Euler continuation with spherical-inflow boundary conditions.
3. **Plate** — shaped immersed boundary, sustained-feed integration window, flat and
   parabolic families. Produces impulse, pressure field, and the flux map.
4. **Cold path** — flux map → 1-D ablating-wall model at representative radii → ablated
   mass → Isp denominator; plus both tapers and the frontier.

### 7.2 Physics tiers

- **Tier 1 — multi-material hydrodynamics with a real EOS, no radiation.** Answers the
  first-order question: does a finite-mass tamper turn the plume, or is it blown out of the
  way? If the realization fraction at τ = 1 does not approach 62.9%, the tamper is dead and
  later tiers are unnecessary.
- **Tier 2 — add radiation transport.** Required at the plate, where `τ` straddles 1 (§5.3)
  and where the thermal load is the binding constraint. Flux-limited diffusion with
  Rosseland means in the diffusion coefficient and Planck means in the emission source.
- **Explicitly not required:** material strength (everything is far past melt on µs
  timescales), chemistry beyond dissociation/ionisation, gravity, MHD.

**Geometry:** 2-D axisymmetric `(r, z)`. The problem has a genuine symmetry axis. 1-D
spherical is valid for the near field over the tamper's operating window (§10, Rung 1);
3-D buys nothing until RT is addressed, and then it buys a great deal (§6.7).

### 7.3 EOS and opacity

- **Water/ice:** equilibrium EOS with dissociation and first ionisation over ρ 10⁻⁴–10³
  kg/m³ and T 300–10⁵ K. The project's extended-grid table meets this.
- **Snow at 70 kg/m³ is not a standard EOS entry.** It must be modelled as porous ice — a
  P-α or ε-α compaction model layered on the ice EOS. Porous compaction is dissipative and
  changes how the projectile's energy is deposited; **it must not be approximated by simply
  rescaling the ice EOS density.** This is the single largest new physics item.
- **Opacity:** real tabulated per-regime opacities, not Kramers (§5.3).
- **Interlayer density is a swept axis, not a fixed choice.** Impedance matching argues for
  roughly the geometric mean of the fireball and tamper densities — **~400 kg/m³, packed
  slush rather than snow at 70.** Snow may be too light: shocked to very high velocity, it
  would arrive at the tamper as a thin fast sheet and ram rather than press, defeating the
  purpose of Arm D.

### 7.4 Initial and boundary conditions

- **Near field:** layered spherical/axisymmetric initial condition — projectile, slug,
  interlayer, tamper. **The projectile's energy deposition is resolved, not prescribed**
  (§6.2), since projectile geometry is a swept design variable and §6.2 shows the answer is
  strongly sensitive to it. A prescribed-deposition mode is retained as a cheap screening
  path and as the control that isolates deposition error from tamper physics.
- **Plate:** rigid, non-ablating immersed boundary for the impulse calculation; the ablating
  wall enters through the cold path (§7.1 stage 4).
- **Far boundary:** absorbing/outflow.

### 7.5 Resolution requirements

- ≥20 cells across the tamper thickness (**≤1.75 mm** for a 3.5 cm tamper).
- Convergence demonstrated by grid refinement on every quoted result, not assumed.
- **A cautionary precedent from this project:** a defect in an immersed-boundary
  implementation (a dish rim's side face omitted) silently biased every curved-plate result
  by ~1% and produced a physically plausible but spurious signal that survived review for
  months. Geometry handling gets its own targeted tests.

### 7.6 Diagnostics to instrument

1. **Axial impulse on the plate vs time**, `∫∫(ρv_z² + P)·(n̂·ẑ) dA dt`, with the `ρv²` and
   `P` terms reported **separately** — every prior model kept only `ρv²`, so the size of the
   pressure term is itself a headline result.
2. **Realization fraction** — delivered impulse / `w(√(1+K) − 1)`. The deliverable.
3. **Entropy budget of the tamper** — this is what the piston framing says is lost (§3.3).
4. **Tamper CM momentum vs time** — the credited quantity, and the measure of recoil.
5. **Fraction of tamper mass ending up plate-directed** vs entrained (the §6.7 mixing
   question, insofar as an axisymmetric code can bound it).
6. **Angular distribution of mass flux across the control sphere** — directly comparable to
   the analytic isotropic assumption; shows whether the tamper *turns* the plume or merely
   *stops* it.
7. **Flux and pressure maps over the plate surface** — feed the taper calculations and the
   ablation model.
8. **Mass, momentum and energy audits at every dump**, closing to <1%. The entire argument
   is an energy-to-momentum conversion, so this is the minimum bar.

---

## 8. Validation and acceptance criteria

No result is trusted before these pass. Written **before** the code, per project convention.

**Analytic anchors specific to this device:**

- [ ] **Ceiling.** In the perfectly-collimated limit, `J → w(√(1+K) − 1)` per projectile kg.
- [ ] **Bare ballistic limit.** At artificially low density (collisionless expansion),
      `β_bare(7.06) → 0.9087` and the capture fraction → 0.3118.
- [ ] **`k ≤ 1` zero-thrust floor.** Ballistically nothing reaches the plate below `k = 1`.
      A hydrocode should show *small but non-zero* thrust from pressure — **how much is
      itself a useful measurement** of how wrong the ballistic model is, and it is cheap.
- [ ] **Free-plate elastic bound.** Reflected/incident → 0.274 at τ = 1 in the ballistic
      limit (§6.3).
- [ ] **`k → 0` degeneracy.** Zero net impulse as slug mass vanishes. A code that produces
      thrust at `k = 0` is wrong.

**Standard verification:**

- [ ] Sod shock tube; Sedov blast; Noh implosion (the axisymmetric source term).
- [ ] Marshak wave for the flux-limited diffusion.
- [ ] Smooth-flow order-of-accuracy test at the scheme's formal rate.
- [ ] Two-material shock-tube with an exact interface solution (**new** — the existing
      kernels are single-material).
- [ ] Porous-compaction Hugoniot against published P-α data for porous ice or a comparable
      material (**new**).
- [ ] **Hypervelocity penetration benchmark** (**new**, required once projectile geometry is
      a resolved variable, §6.2): reproduce a published ice-on-ice or ice-on-porous-ice
      impact — crater/penetration depth and energy partition — against experimental or
      impact-code reference data. This is a distinct validation domain from the shock-tube
      family and should not be assumed to come free with them.

**Audits:**

- [ ] Mass, momentum, energy closing to <1% at every dump.
- [ ] Grid convergence on every quoted number.

---

## 9. What already exists in this repository

| component | what it is | reuse |
|---|---|---|
| `crates/hydro1d` | 1-D **planar** Lagrangian rad-hydro: staggered mesh + artificial viscosity, pluggable EOS, flux-limited diffusion (Levermore–Pomraning), coupled gas+solid conduction, two-phase condensation, **ablating wall** (`Q*` surface balance, blowing factor, vapour shield). Single material. | **Extend** to spherical geometry + multi-material (Rung 1). Its ablating wall is used unchanged in Rung 6. |
| `crates/euler2d` | 2-D axisymmetric Euler: HLLC Godunov, MUSCL-Hancock 2nd order, Strang splitting, conservative cylindrical source, **ghost-cell immersed boundary** for a shaped plate. Ideal gas, single material. | **Extend** with the table EOS and multi-material (Rung 3); its immersed boundary is reused for the plate (Rung 4). |
| `crates/tables` | Shared JSON table loader. | Unchanged. |
| `crates/sweep` | Rayon-parallel sweep driver, JSONL output. | Extend with new sweep modes. |
| `python/puffsat` | EOS/opacity table generation (CoolProp + analytic Saha/CEA-style, TOPS/OPLIB overlay), frontier extraction, analysis, plotting. | **Extend** with the porosity model; reuse the analysis and plotting pipeline. |
| `data/tables/water_jupiter.json` | Extended-grid equilibrium water table, T to 1.2×10⁶ K, ρ to 10⁻⁴, real gray opacities. | **Baseline table for this study** (§5.1). |
| Build | Top-level `Makefile` over `cargo` + `uv`; JSON tables in, JSONL results out. | Unchanged. |

**Verified working and directly transferable:** the analytic equilibrium water EOS, the
FLD implementation, the ablating-wall model, the immersed-boundary plate, the sweep/analysis
plumbing, and the validation discipline.

**Genuinely missing:** multi-material tracking, spherical 1-D geometry, a table EOS inside
the 2-D kernel, porous-ice compaction, spherical-inflow boundary conditions, and a
sustained-feed integration window.

---

## 10. Work plan

**Definition of done.** The study ends when it returns **effective Isp for Arm D within a
stated bracket**, together with a pay/don't-pay verdict on the tamper against the 62.9%
threshold (§3.4), the projectile-economy curve (§3.1), and the optimal mass ratio `K` (§3.5).
The result is explicitly **single-code and preliminary**; the independent hydrocode
cross-check is deferred and named as the outstanding validation. A bracket that straddles
62.9% is a legitimate outcome and is reported as such — the deliverable is the bracket, not a
verdict forced past the evidence.

**Two passes (§3.1), executed in strict order.** Rungs 1–5 constitute **Pass 1**, which
excludes ablator mass and returns an *upper bound* on Isp with the optimum at `K* ≈ 7.06`.
Rung 6 adds the measured ablator (**Pass 2**), revising Isp down and shifting `K*` up.

**Uncertainty reduction comes before cost refinement.** RT is *inside* Pass 1's bracket and
is its widest contributor; the ablator is *outside* Pass 1 entirely and moves the absolute
number without changing the tamper verdict (§6.5). So RT is Rung 2 — immediately after the
run that produces the acceleration history it needs — and the ablator is last.

**Rungs 1 + 2 already give a preliminary verdict**, combining the 1-D full-shell result with
the analytic capture geometry. That is the cheapest point at which to decide whether to
continue, and it precedes any 2-D work.

### Rung 1 — 1-D spherical multi-material (the new physics, isolated)

*Why 1-D is valid here, not merely cheap: lateral communication across the fireball takes
`r/c_s ≈ 69 µs`, against a 45 µs confinement time. The hemispheres are dynamically
independent for the entire window in which the tamper does its work.*

- [ ] Add spherical geometry to the 1-D Lagrangian kernel.
- [ ] Add a per-cell material index (projectile / slug / interlayer / tamper).
- [ ] Build the porous-ice (P-α or ε-α) compaction model and its table pathway.
- [ ] Add resolved projectile penetration and energy deposition, with a prescribed-deposition
      mode retained as the control (§7.4).
- [ ] Write the §8 acceptance tests **first**; make them pass.
- [ ] Sweep **Arm D** over `(ρ_interlayer, standoff, τ)`; run **Arm B as a control** and
      confirm the code reproduces its RT disruption (§6.7) rather than assuming it.
- [ ] Sweep `K` over at least 6–32 and **locate the Pass-1 design point** (expected at the
      ceiling optimum `K* ≈ 7.06`, since Pass 1 carries no ablator term), reporting Isp,
      projectile economy, and energy per unit impulse together (§3.5). Report how the tamper
      moves the optimum relative to bare.
- [ ] Sweep **projectile geometry** — areal density, aspect ratio, bulk density — against the
      vertex-hole and deliverability constraints (§6.2). Report the deposited-energy fraction
      as a first-class output, since 34% vs 80% changes the device.
- [ ] Report: realization fraction, tamper entropy budget, tamper CM momentum, `Θ` scaling
      confirmed, and `Λ` for a full-shell tamper as a bound.
- [ ] **Emit the acceleration history `a(t)` at the plume/tamper interface** — the input
      Rung 2 needs, and free from this run.
- [ ] **Gate:** if the realization fraction cannot approach 62.9%, the tamper is dead and
      the remaining rungs are descoped to a bare-plate Isp confirmation.

### Rung 2 — RT treatment for Arm D (the widest Pass-1 uncertainty)

Arm D's mixed fraction is currently bounded only at **16–63%**, spanning ~58–62% of ceiling
against a 62.9% threshold (§6.7). **Until this narrows, no amount of downstream work can
produce a verdict**, so it runs before any 2-D effort. Escalate only as far as needed:

- [ ] **Tightened analytic bound.** Integrate RT growth against Rung 1's *actual* decaying
      `a(t)` instead of the constant-acceleration self-similar law, and against a stated
      initial perturbation spectrum rather than an assumed broad one. Cheap; may alone be
      enough, since the current bound is deliberately crude in both respects.
- [ ] **Gate:** if the tightened bound no longer straddles 62.9%, stop here.
- [ ] **Mix model**, if it still straddles — a buoyancy-drag or K-L-style mix width at the
      interface in the 1-D kernel, calibrated against published RT mixing data.
- [ ] **Resolved spot-check**, only if the mix model is inconclusive: a 2-D or 3-D
      interface-resolved run at the single design point. Expensive and last.
- [ ] Report the RT contribution as an explicit ± on the realization fraction, whatever
      level it is settled at.

### Rung 3 — 2-D axisymmetric multi-material

- [ ] Port the table EOS into the 2-D kernel.
- [ ] Add multi-material tracking.
- [ ] Reproduce Rung 1 in the full-shell limit (**gate**).
- [ ] Sweep tamper **angular coverage** at fixed mass (§6.7 — coverage, not curvature).
- [ ] Emit the control-sphere mass/velocity/angle distribution.

### Rung 4 — transport and plate

- [ ] Spherical-inflow boundary condition from the control-sphere distribution.
- [ ] Coarse far-field continuation to the plate (not free-streaming, §7.1).
- [ ] Replace the decay-based integration window with a sustained-feed criterion (§6.4).
- [ ] Sweep plate radius `R` (mass-ceiling constrained), standoff `d`, and shape: flat plus
      the paraboloid family including `d/D` up to ~0.35.
- [ ] Emit plate flux and pressure maps.

### Rung 5 — Pass 1 deliverable

**This is the study's answer to the question it was created to ask.**

- [ ] Frozen-recombination bracket at 57 eV/molecule (§13.1) — both ends.
- [ ] Fold interlayer-density and projectile-deposition sensitivity into the bracket.
- [ ] Both taper calculations (closed form, §6.6).
- [ ] Assemble **effective Isp for Arm D within a stated bracket**, on carried mass
      *excluding* the ablator — explicitly labelled an upper bound (§3.1).
- [ ] Report the **realization fraction against 62.9%** and the resulting pay/don't-pay
      verdict on the tamper, or an honest "cannot distinguish" if the bracket straddles it.
- [ ] Report the projectile-economy curve and the Pass-1 optimum `K` (§3.5).

### Rung 6 — ablator and Pass 2 (last)

The ablator moves the absolute Isp but not the tamper verdict, and *excluding* it is
conservative for the tamper (§6.5). So it is refined only once the verdict exists.

- [ ] Run the existing 1-D ablating-wall model at this study's areal density, arrival
      velocity, and 750 µs residence time to collapse the **27× uncertainty** (§6.5).
- [ ] Flux map → ablating-wall model at representative radii → ablated mass per
      configuration, with its `Q*` and opacity bracket, and as a mass flow (§6.5.1).
- [ ] Report **burn-through margin** — the one abrupt failure mode, and what actually binds
      on the plate side (§6.5.2).
- [ ] Sweep **pulse mass** to quantify the sub-linearity gain: larger pulses at
      proportionally lower cadence should cut total ablator mass at fixed thrust (§4.1).
- [ ] Confirm the **self-limiting plate equilibrium** (§6.5.2) against the coupled model, so
      the steel baseline rests on a computed interface condition rather than a closed form.
- [ ] **Fold the ablator into the denominator** and re-derive Isp — the Pass 1 → Pass 2 step
      (§3.1). Report both, with Pass 1 still labelled an upper bound.
- [ ] **Relocate the optimum `K`** with the ablator term present (§3.5).
- [ ] **Gate:** if burn-through margin is inadequate, invoke the §6.5.4 escalation path
      (ceramic face → Nb/Mo → carbon-carbon) rather than changing the architecture.
- [ ] Frontier extraction, plots, and the committed deliverable artifacts.

**Deferred:** the independent hydrocode cross-check. The result is therefore reported as
**single-code and preliminary**, following this project's existing precedent that a
preliminary number is quotable provided it is labelled as such and the cross-check is named
as the outstanding validation.

---

## 11. Decisions taken, and why

| # | Decision | Rationale |
|---|---|---|
| D1 | **The tamper is an isentropic piston, not a cold mirror.** Instrument entropy and CM momentum. | "Reflects all energy, absorbs none, still recoils" is self-contradictory at finite mass. Under §3.3(b) the recoil energy is *credited*; entropy is what is actually lost. |
| D2 | **Arm D (filled interlayer) is the design; Arm B (vacuum gap) is provisionally foreclosed and retained as a control.** | Two independent arguments close Arm B: it is a re-thermalising ram rather than a piston (§3.3), and RT disrupts its tamper completely before it acts (§6.7). Effort concentrates on narrowing Arm D's bracket. |
| D2a | **The mass ratio `K` is optimised here, against Isp subject to the plate thermal budget — not against an economic weighting.** | The Isp optimum is real but flat (±0.6% over `K` = 6–8), while heat and projectile consumption per unit impulse both scale as `1/β` and improve steeply. The binding constraint is thermal, and it is in scope (§3.5). |
| D2b | **Projectiles are not priced into a combined figure of merit.** | Converting projectiles to payload-equivalent needs a program economic model outside this study; importing it would make a physics result move with someone else's assumptions. Projectile economy is reported as `β` (§3.1). |
| D3 | **Plate radius swept; plate mass a ceiling.** | §6.5: `R` is the only lever that breaks the capture-vs-fluence conflict. |
| D4 | **Isp denominator = all expended carried mass, ablator included — but evaluated in two passes**, Pass 1 excluding it as an explicit upper bound. | Isp is a rocket-equation quantity, so anything carried and expended is charged. But the ablator is uncertain by 27× and answering it is a *plate* question; excluding it in Pass 1 decouples the *tamper* question and lets Rungs 1–5 proceed without waiting, with RT narrowed first (§3.1). |
| D4a | **Pulse mass and cadence trade freely at fixed thrust; neither is independently pinned.** | `E/J = w/(2β)` contains neither, so average heat load and gravity loss are invariant under the trade — while sub-linear ablation means larger, rarer pulses cut ablator mass (§4.1, §6.5.1). |
| D4b | **Plate material: steel + thin renewed oil, with a thin ceramic hedge.** | An ablating surface is temperature-pinned, so the plate self-limits near 750–905 K; strength, spall, and bending are all non-binding at 2 MPa. Escalation is triggered only by inadequate burn-through margin (§6.5.2, §6.5.4). |
| D5 | **Plate shape sweep = flat + parabola family + both tapers.** Shape itself stays open. | §6.6: the prior foreclosure is conditional on plane-wave incidence; but the parabola's area penalty may exceed its impulse gain. Genuinely two-sided. |
| D6 | **Keep the 1-D-thermophysics × 2-D-geometry factorization**, adding a cold-path flux map. | The gas is optically thick, so radiation stays local-diffusive and the coupling is one-way. A monolithic 2-D rad-hydro across a sweep is not affordable. |
| D7 | **Ablator mass is an emergent cost per configuration, not a specified thickness.** | §6.5: ablation is the plate's only thermal sink, so it is forced by physics. Specifying a thickness the balance does not respect would silently burn through. |
| D8 | **Tamper design variable is angular coverage, not curvature.** | §6.7: RT destroys sub-metre features within the confinement window. |
| D9 | **Build order: 1-D spherical first, then 2-D.** | Isolates new physics (multi-material, porosity) from new geometry, with known-answer anchors at each step. |

---

## 12. Departures from prior decisions in this repository

Each of these reverses or narrows a recorded decision and requires its own record before the
work lands.

| Prior decision | Transfers? | Why |
|---|---|---|
| 1-D/2-D factorization | **Yes** (D6) | With a cold-path flux map added |
| Facesheet survivability ladder, peak `≈1.2ρv²` | Yes, but **not binding** | ~200× margin (§6.4) |
| Ablating-wall model | **Yes, and promoted to the critical path** | It is now the largest unknown in the denominator |
| "The ablator is not a pressure device" | Yes, and confirmed | But moot — pressure is not the constraint here |
| Frozen-recombination bracket | Yes, and **worse here** | 57 eV/molecule vs the prior study's regime |
| Rigid-during-pulse gate | Yes, passes trivially | §6.4 |
| `10⁻³`-of-peak integration window | **No** | Wrong for a 750 µs sustained feed (§6.4) |
| RT/RM deferral | **No** | Load-bearing in three places (§6.7) |
| Deep-dish foreclosure | **No** | Conditional on plane-wave incidence (§6.6) |
| Ablator, vehicle scale and cadence held out of scope | **No** | An Isp deliverable pulls all three inside the boundary |
| Inter-pulse plate thermal accumulation excluded as cadence-dependent | **Partly** | Now computed rather than excluded — but the answer is that it does not bind: the plate self-limits at `T_abl` (§6.5.2) |
| Projectile geometry treated as an unrecorded external input | **No** | §6.2 makes it a swept design variable with a real interior optimum |
| No voids behind the hot face (spall risk from free-surface reflection) | **No** | Derived at 400 MPa–2 GPa; peak here is ~2 MPa. Not currently exercised, since cooling proved unnecessary (§6.5.3) |
| Carbon-carbon rejected as a hot face (burns in atomic O) | **Weakened** | Exposure is only during the pulse, when renewed oil covers it; between pulses there is no gas. Held in reserve on the escalation path, not adopted (§6.5.4) |
| SiC + Ti as the settled hot-face stack | **Departed** | Selected there for oxidation and per-pulse thermal shock at GPa loads. Here loads are ~2 MPa and the plate self-limits, so **steel + oil is the baseline** and ceramic is a burn-through hedge (§6.5.4) |
| Existing `eta_capture` sweep data | Cross-check only | Incidence geometry differs |

**Records written.** Each clears the project's bar for an architecture decision record —
hard to reverse, surprising without context, and the result of a genuine trade-off:

- [x] **[ADR-0030](docs/adr/0030-tamper-isentropic-piston-not-mirror.md)** — the tamper is an
      isentropic piston, not a mirror; recoil is credited; the deliverable is realization
      fraction against a single ceiling. Replaces the prior `Λ > 1 + τ` criterion (§3.3).
- [x] **[ADR-0031](docs/adr/0031-isp-deliverable-pulls-scope-inward.md)** — an Isp deliverable
      pulls ablator mass, vehicle scale, and cadence inside the scope boundary that prior
      scope places outside, and inter-pulse thermal accumulation with them (§3.1).
- [x] **[ADR-0032](docs/adr/0032-deep-dish-foreclosure-is-incidence-conditional.md)** — the
      deep-dish foreclosure is conditional on plane-wave incidence and inverts under a point
      source (§6.6). ADR-0021 carries a pointer and is otherwise unchanged.
- [x] **[ADR-0033](docs/adr/0033-rt-deferral-does-not-transfer-to-the-tamper.md)** — the
      RT/RM deferral does not transfer; RT is load-bearing here in three places (§6.7).
      ADR-0020 carries a pointer and is otherwise unchanged.

Canonical terms for this study are in [`CONTEXT.md`](CONTEXT.md) under *Language — tamped-nozzle
study*; §14 below is the reading-order glossary for this document.

---

## 13. Open questions and assumptions register

State these in any write-up. Each could move the answer.

1. **RT/RM is un-modelled and, for Arm D, wider than the margin.** Mixing decides whether
   the tamper is a piston or entrained payload. The analytic bound (§6.7) closes Arm B but
   leaves Arm D at **16–63% mixed**, spanning ~58–62% of ceiling against a 62.9% threshold.
   An axisymmetric code cannot narrow it. **Top-ranked uncertainty-reduction target**
   (§13.1, Rung 2 — it now runs *before* the ablator, not after).
2. **Interlayer density may want to be ~400 kg/m³, not 70.** Impedance matching argues for
   packed slush over snow. Cheap to sweep; a real risk to Arm D as currently specified.
3. **~~Vehicle mass~~ — RESOLVED: 1000 t.** The plate is 5% of the vehicle and the prior
   10 t/50 t inconsistency is gone (§4.1).
4. **~~Cadence~~ — RESOLVED: not thermally limited.** The plate self-limits at the ablator's
   vaporisation temperature (§6.5.2), and average heat load is invariant under the
   pulse-mass/cadence trade (§4.1). Cadence is set by gravity losses alone, and larger pulses
   at lower rate is a *free* variation that also cuts ablator mass. What remains open is the
   shock-absorber stroke this implies, which is out of scope.
5. **~~Projectile geometry~~ — RESOLVED: it is in scope as a swept design variable** (§6.2).
   The residual risk is that resolving hypervelocity penetration is a distinct validation
   domain (§8) and could dominate Rung 1's effort.
6. **~~Plate material~~ — LARGELY RESOLVED: steel + thin renewed oil, with a ceramic hedge.**
   Every candidate criterion is non-binding except steady-state temperature, and that
   self-limits at `T_abl` (§6.5.2, §6.5.4). **What replaces it as the open question is
   burn-through margin** — the one abrupt failure mode, since a mid-pulse breach exposes the
   substrate to 50–80 kK plume with no temperature pinning. Rung 6 reports it. Escalation
   materials are named but not invoked.
7. **`K` and `τ` grids not pinned.** §3.5 shows the design point is expected *above* the
   `K* ≈ 7.06` Isp peak — plausibly near `K ≈ 14` if the thermal budget binds at ~2× — so
   the grid must span at least `K` = 6–32 rather than clustering around 7. `τ` spans 0–2.
   The Pass-1 grid is set now; the Pass-2 shift depends on
   Rung 6's ablator number.
8. **~~Projectile pricing~~ — RESOLVED: not priced into a combined metric.** An earlier
   draft proposed converting projectiles into payload-equivalent at an externally-supplied
   rate. Dropped: it would make a physics result depend on a program economic model outside
   this study. Projectile economy is reported as `β`, which §3.5 shows is the *same*
   function as plate heat load per unit impulse — so the trade is internal to the physics
   and needs no exchange rate.
9. **`U_shock ≈ 12 km/s` in ice at ~70 GPa is an estimate**, and `Θ` depends on it linearly.
10. **The snowplow penetration model (§6.2) is crude** — no lateral spreading, coherent
    penetrator. It is a screen, not a result.
11. **Arm B is *provisionally* foreclosed, on two analytic arguments** (§4.2, §6.7), not on
    a simulation. Reopening it requires new evidence, but the foreclosure should be
    confirmed by the control run rather than assumed.

### 13.1 Uncertainty budget for the Arm D Isp bracket

The deliverable is a bracket (§10), so its width is a first-class output. Contributors,
ranked by current width, with what collapses each:

| # | contributor | current width | affects | in Pass 1? | collapsed by |
|---|---|---|---|---|---|
| 1 | **RT mix fraction at the plume/tamper interface** | 16–63% mixed → ~58–62% of ceiling | realization fraction | **yes** | **Rung 2** — bound, mix model, or resolved spot-check |
| 2 | **Frozen recombination** at 57 eV/molecule | this project's largest quantified physics uncertainty, in a *cooler* regime than here | realization fraction | **yes** | Rung 5 bracket, both ends |
| 3 | **Projectile deposition** (34% vs 75% of energy into the slug) | changes which body is the fireball | everything | **yes** | Rung 1 — resolved penetration |
| 4 | **Interlayer density** (70 vs ~400 kg/m³) | may invert Arm D's loading mechanism | realization fraction | **yes** | Rung 1 — cheap sweep axis |
| 5 | **Ablator mass** (vapour-shielding balance) | **27×** → 2–60% of the mass budget | denominator | **no — excluded by construction** | Rung 6 (last) |
| 6 | **Opacity** near `τ ~ 1` at the plate | prior experience: 2000× error from the wrong table | flux map, denominator | partly | use the real-opacity extended table |
| 7 | **`U_shock` in ice** | linear in `Θ` | tamper timing | yes | literature lookup |

**The two-pass split (§3.1) is what makes this tractable.** The single widest contributor —
the ablator, at 27× — sits in the denominator and is *excluded* from Pass 1, so it cannot
contaminate the tamper verdict. Pass 1's bracket is set by contributors 1–4, of which
**RT alone still exceeds the 62.9% margin being decided**. Until that is narrowed, Pass 1
cannot distinguish "the tamper pays" from "the tamper is dead," and every interim number must
carry its bracket and its upper-bound label.

---

## 14. Glossary

| term | meaning |
|---|---|
| **PuffSat / projectile** | Mass accelerated by external infrastructure and aimed at the vehicle. Not carried, so not charged in Isp. |
| **Slug** | Carried mass the projectile buries itself in. Vaporised; becomes the fireball. |
| **Tamper** | Dense carried mass beyond the slug that turns part of the away-going fireball around. Fully vaporised; works by inertia (areal density), not by survival. |
| **Interlayer** | Low-density carried mass filling the standoff between slug and tamper in Arm D. Both spacer and pressure-transmitting medium. |
| **Pusher plate** | The permanent structure that receives the impulse. Has a hole at its vertex for the projectile. |
| **`k`** | Slug mass / projectile mass. ≈7 at the bare optimum. A continuous design variable. |
| **`τ`** | Tamper mass per slug kg. |
| **`K = k(1+τ)`** | Total *carried* mass per projectile kg — the only quantity the ceiling depends on. |
| **`w`** | Closing speed, 75 km/s. |
| **Capture fraction** | `(1 − 1/√k)/2` — share of an isotropic fireball that outruns its own recoil. Zero at `k ≤ 1`. |
| **Realization fraction** | Delivered impulse / `w(√(1+K) − 1)`. **The deliverable.** Must exceed 62.9% at τ = 1 for the tamper to pay. |
| **Projectile economy** | Impulse per projectile kg, `β·w`. The second reported metric, alongside effective Isp. Not combined with it — see §3.1. |
| **Arm D / Arm B** | Filled-interlayer configuration (the design) / vacuum-gap configuration (provisionally foreclosed, retained as a control). |
| **Transit ratio `Θ`** | Slug disassembly time / tamper shock-transit time. Must exceed ~3 for the tamper to function. |
| **Arrival window** | Duration over which fireball material reaches the plate, set by its velocity dispersion. ~750 µs here. |
| **Areal density `σ`** | Mass per unit area. The tamper's figure of merit, conserved through vaporisation. |
