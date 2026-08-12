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

## 0. Symbols and notation

Every symbol used later in this document is defined here, so it can be read front to back
without forward references. Values in the *reference* column are the analytic, no-interlayer
reference case of §4 (`w = 75 km/s`, `k = 7.06`, `μ = 0`, 200 kg total encounter mass) and
are quoted for orientation only — the section cited is authoritative. Filled Arm D has
`μ > 0`, so those reference masses do not describe Arm D. §0.8 lists the few deliberately
similar symbols; §14 is the glossary of named concepts, and [`CONTEXT.md`](CONTEXT.md) holds
the canonical prose terms.

**Frame and sign conventions.** Cylindrical `(r, z)`, axisymmetric about the device axis.
**`+z` is the thrust direction.** The projectile arrives travelling `−z`, so its momentum is a
debit, and useful ejecta also leaves travelling `−z`; impulse `J` is counted positive in `+z`.
This is the origin of the `−1` in every impulse law (§2.1).

**Normalisation.** `β`, `k`, `τ_t`, `μ`, `K`, `K_ej`, and `C` are dimensionless. Lower-case
`j = J/m_i` denotes impulse per projectile mass, with units `N·s/kg = m/s`; upper-case `J`
denotes per-pulse impulse. `m_charged`, `E`, and the figures in §4 are per pulse unless stated
otherwise. The two normalisations are interchangeable for Isp, since its numerator and
denominator both scale with `m_i`.

### 0.1 Masses, mass ratios, and the deliverables

| symbol | meaning | reference |
|---|---|---|
| `m_i` | projectile (PuffSat) mass per pulse — **never charged in Isp**, because it is supplied externally rather than carried by the vehicle. Written `m_projectile` where the contrast with carried mass is the point (§3.1) | 24.81 kg at `τ_t = 0`; 13.23 kg at `τ_t = 1`, both with `μ = 0` (§4) |
| `m_s`, `m_t`, `m_int` | slug, tamper, and interlayer mass per pulse | 175.2 kg / — / — at `τ_t = μ = 0`; 93.4 kg / 93.4 kg / — at `τ_t = 1, μ = 0` |
| `m_hydro` | all nonprojectile carried mass represented in the near-field hydrodynamics: at minimum `m_s + m_t + m_int`, plus any consumable spacer/support or added material region | 175.2 kg (bare reference) |
| `m_abl` | ablator mass expended per pulse; zero by construction in Pass 1 and measured per configuration in Pass 2 | — |
| `m_other` | other expended carried mass not represented in the near-field hydrodynamics. Must be named per configuration; permanent dry structure is excluded | zero in current references |
| `m_enc = m_i+m_hydro` | **encounter mass per pulse** — projectile plus near-field carried mass. This is what “200 kg reference encounter” means; it is not the Isp denominator | 200 kg in §4 |
| `m_charged = m_hydro+m_abl+m_other` | all expended **carried** mass per pulse. The Isp denominator (§3.1); a consumable Arm B spacer is charged here and included in `m_hydro` if simulated | 175.2 kg (Pass 1, bare reference) |
| `k` | **slug ratio** — slug mass per projectile kg. Continuous, not a count | 7.06; bare-plate optimum `k*` = 7.060 |
| `τ_t = m_t/m_s` | **tamper ratio** — tamper mass per slug kg. The subscript distinguishes it from optical depth `τ_opt` | swept 0–2; `τ_t = 1` is the analytic reference tamper |
| `μ = m_int/m_s` | **interlayer ratio** — interlayer mass per slug kg. Zero for the bare and vacuum-gap references; nonzero for filled Arm D | swept, not yet pinned |
| `K = m_hydro/m_i` | **hydrodynamic carried-mass ratio**. For the named slug/tamper/interlayer regions only, `K = k(1+τ_t+μ)`; additional simulated material adds its own mass ratio | 7.06 (`τ_t = μ = 0`) / 14.12 (`τ_t = 1, μ = 0`); swept 6–32 |
| `a_abl = m_abl/m_i` | ablator mass expended per projectile kg | zero in Pass 1; measured in Pass 2 |
| `a_other = m_other/m_i` | other expended carried mass per projectile kg | zero in current references |
| `C = m_charged/m_i` | **charged-mass ratio** — the denominator mass per projectile kg, `C = K + a_abl + a_other` | equals `K` in §4 |
| `K_ej = (M_ej−m_i)/m_i` | **ejecta-mass ratio** — nonprojectile mass leaving with the event's ejecta per projectile kg. The ceiling depends on `K_ej`, not on how that mass was allocated | `K_ej = K` in Pass 1; see §3.2 for Pass 2 |
| `K*`, `k*` | the value of `K` (or `k`) that maximises Isp **for a named configuration family and pass**; there is no universal `K*` | bare Pass-1 control `k*_bare = 7.060`; other optima are measured (§3.5) |
| `Isp_eff = J/(g₀·m_charged) = βw/(g₀C)` | **effective specific impulse** — the final carried-mass-economy deliverable (§3.1) | 984 s bare plate (Pass 1, upper bound) |
| `g₀` | standard gravity, 9.80665 m/s² | — |
| `r_real = J/J_max = β/(√(1+K_ej)−1)` | **realization fraction** — net vehicle impulse as a share of the ceiling for the same ejecta mass (§3.3). A hydrodynamic comparison metric feeding Isp, not a replacement for Isp. The subscript distinguishes it from radial coordinate `r` | 49.4% bare; 62.9% is only the §3.4 reference-case break-even |
| **ballistic capture fraction** | `max[0, (1 − 1/√k)/2]` — share of the isotropic, pressure-free ballistic fireball that outruns its own recoil and reaches the plate. Its zero below `k = 1` is a ballistic-model limit, not a general zero-thrust claim (§2.2, §8) | 31.2% at `k = 7.06` |

### 0.2 Speeds, impulse, and energy

| symbol | meaning | reference |
|---|---|---|
| `w` | **closing speed** of projectile and vehicle, head-on, vehicle frame | 75 km/s (envelope 74–81) |
| `w̄` | mass-weighted mean closing speed, used only when reproducing prior work's Isp figure | 77.28 km/s (§3.2) |
| `V = w/(1+k)` | merged blob's centre-of-mass recoil speed, directed *away* from the plate | 9.31 km/s |
| `u = w√k/(1+k)` | fireball expansion speed in the blob frame. `u/V = √k` sets the ballistic capture fraction | 24.72 km/s |
| `U_shock` | shock speed through the tamper (ice Hugoniot at ~70 GPa) — **an estimate** (§13.9) | ≈12 km/s |
| `c_s` | sound speed | ≈9.8 km/s in the fireball |
| `M` | Mach number | ≈2.5 terminal, at the transport hand-off (§7.1) |
| `Δv_vehicle` | vehicle velocity increment per pulse, `J/M_vehicle` (§4.1) | 1.69 m/s per pulse |
| `Δv_t` | tamper velocity change used in the recoil-time estimate `σ_t·Δv_t/P` (§6.3); not the vehicle's pulse increment | — |
| `J` | **net axial impulse delivered to the vehicle**, positive in `+z` and including the incoming projectile's `−m_i·w` momentum debit; `J = β·m_i·w` | 1.69×10⁶ N·s per 200 kg reference encounter |
| `j = J/m_i` | net axial vehicle impulse per projectile mass | `βw` |
| `P_ejecta` | signed total axial ejecta momentum; useful, plate-opposed ejecta has `P_ejecta < 0` | — |
| `M_ej = m_i(1+K_ej)` | total ejecta mass, including the projectile | — |
| `j_max = w(√(1+K_ej) − 1)`; `J_max = m_i·j_max` | **the thermodynamic ceiling**, respectively per projectile mass and per pulse (§3.2) | — |
| `v_e,max = j_max/K_ej` | net effective exhaust velocity at the ceiling, normalised by nonprojectile ejecta mass | — |
| `β = J/(m_i·w)` | **dimensionless net-vehicle-impulse coefficient**: net vehicle impulse divided by the magnitude of the incoming projectile momentum. It includes the projectile momentum debit. `β·w = J/m_i` has velocity units but is not the speed or velocity increment of any material or body. `β_bare`, `β_tamp` are the ballistic closed forms; `β_ideal = √(1+k) − 1` | `β_bare(7.06)` = 0.9087, `β_bare(14.12)` = 1.6835 |
| **projectile economy** | `β·w = J/m_i` — net vehicle impulse per projectile mass, the second reported metric (§3.1); a velocity-equivalent impulse normalisation, not a physical velocity | — |
| `E` | incoming projectile kinetic energy per pulse, `½m_i w²`; `E/J = w/(2β)` is energy per unit impulse (§3.5) | 69.8 GJ (`τ_t = 0`) / 37.2 GJ (`τ_t = 1`), both at `μ = 0` |

### 0.3 Geometry

| symbol | meaning | reference |
|---|---|---|
| `θ` | polar angle from the plate-directed `+z` axis; `θ_max` is the angle the plate rim subtends at the source (§6.6) | ballistic material with `cos θ > 1/√k` reaches the plate |
| `r`, `z` | radial and axial coordinate; `r̂`, `ẑ` are unit vectors. `n̂_wall` points from fluid into the plate (the direction of pressure force on it); `n̂_cs` points outward from a fluid control volume | — |
| `v_z` | axial fluid-velocity component, used in the control-surface momentum flux (§7.6) | — |
| `r_slug` | slug radius | 0.290 m (ice) / 0.683 m (snow) |
| `h_t` | **tamper thickness**, `h_t = m_t/(2πs²ρ_t)` (§6.1) | 3.5 cm (snow slug) / 19.3 cm (ice slug, contact) |
| `s` | **tamper stand-off radius** — radius of the tamper shell from the blob centre, so `s = r_slug` at contact. In Arm D the interlayer occupies the radial gap `r_slug < r < s` | 0.29 m (ice) / 0.68 m (snow) at contact; 1 m in the standoff cases |
| `d` | **plate standoff** — axial distance from the fireball's apparent source to the plate vertex. Distinct from `s` | 10 m reference; §7.1 spans to ~25 m |
| `R` | pusher-plate **radius** (its diameter is `2R`) | 15 m reference; up to ~25 m under a mass ceiling |
| `F` | paraboloid focal length; the focus-matched design sets `F = d` (§6.6) | 6–10 m |
| `δ/D` | paraboloid **dish depth over diameter** — the shape parameter, `δ/D = R/(8F)`. `δ` is not the plate standoff `d` | 0.19 at (`R = 15 m`, `F = 10 m`); swept to ~0.35 |
| `A` | named geometric area, as in `σ = m/A`; the reference full-shell tamper uses hemispherical area `2πs²`, while partial coverage uses its actual covered area. Distinct from Atwood number `A_RT` | — |

### 0.4 Material state and thermophysics

| symbol | meaning | reference |
|---|---|---|
| `ρ` | mass density; `ρ_slug`/`ρ_s`, `ρ_tamper`/`ρ_t`, `ρ_interlayer` | ice 917, snow 70, slush ~400 kg/m³ |
| `σ = m/A` | **areal density** — the tamper's figure of merit, conserved through vaporisation (§2.3). `σ_t` (tamper), `σ_proj`, `σ_target`, `σ_plume`, `σ_abl` (ablated per pulse) | `σ_t` = 14.9 kg/m² at `τ_t = 1`, `s = 1 m` |
| `Σ` | gas mass column per unit plate area at arrival — deliberately distinct from generic `σ`. The quoted reference is plate-area-averaged; local solver outputs are functions of radius and time | ≈0.063 kg/m² |
| `e` | specific internal energy. *Distinct from `e_eff`, §0.7* | 305.7 MJ/kg = 57.1 eV per H₂O molecule |
| `T` | temperature (K, or eV where the plasma state matters) | 14 kK fireball; 50–80 kK at stagnation |
| `P` | pressure. *Momentum in `P_ejecta` only* | ~2 MPa peak at the plate; ~70 GPa in the tamper shock |
| `γ` | ratio of specific heats, `γ_eff`; distinct from RT growth rate `γ_RT` | `γ_eff` = 1.25 |
| `κ` | gray opacity | 10–1000 m²/kg — uncertain across two decades |
| `τ_opt = κΣ` | **optical depth** at the plate (§5.3). Distinct from tamper ratio `τ_t` | 0.63–63 — straddles `τ_opt ~ 1` |

### 0.5 Timescales and the transit ratio

| symbol | meaning | reference |
|---|---|---|
| `t` | time | total simulated ~2–3 ms |
| `t_dis = r_slug/u` | **slug disassembly time** — the event the tamper must act within. Spelled `t_disassembly` in §6.1 | 11.7 µs (ice) / 27.6 µs (snow) |
| `t_shock-transit` | shock transit of the tamper, `h_t/U_shock` | 1.4–2.9 µs |
| `Θ = t_dis / t_shock-transit` | **transit ratio** — how many tamper shock-transit times fit within the slug-disassembly time. `Θ ≳ 3` is the current screening gate, not a universal instability threshold; the sweep tests sensitivity (§6.1) | 9.5 (snow slug) / 0.7 (ice slug, contact) |
| **arrival window** | time between the first and last arrivals in a stated velocity or cumulative-mass interval. The ~750 µs estimate uses the 8–20 km/s ballistic range; simulations must state their mass-flux cutoff so numerical tails do not define the interval (§6.4) | ~750 µs at `d = 10 m` |

### 0.6 Plate, ablation, and Rayleigh–Taylor

| symbol | meaning | reference |
|---|---|---|
| `Φ` | time-integrated incident energy flux per unit plate area over one pulse; state whether a quoted value is local, peak, or area-averaged | ≈12.6 MJ/m², plate-area-averaged at `R = 15 m`, `d = 10 m`. An order-of-magnitude reference whose capture-energy basis is not derived in this document; Rung 4's flux map supersedes it (§6.5) |
| `α_abl` | ablation sub-linearity exponent, `σ_abl ∝ Φ^α_abl`, `α_abl < 1` — vapour shielding (§6.5) | — |
| `Q*` | effective energy removed from the incident/wall balance per kg of ablator expelled, used in `ṁ = q_in/Q*`; it is a model parameter that may include phase change and other unresolved losses, not merely latent heat | — |
| `T_abl` | effective ablator vaporisation temperature imposed by the wall model while ablation is active. It bounds the substrate only while an unbroken ablating layer maintains that boundary (§6.5.2) | 800–1000 K → plate 753–905 K |
| `k_cond` | thermal conductivity | 45 W/m/K (steel) |
| `α_th` | thermal diffusivity | 1.2×10⁻⁵ m²/s (steel) |
| `σ_SB` | Stefan–Boltzmann constant, in the `2σ_SB T⁴` two-face radiation term (§6.5.2) | 5.67×10⁻⁸ W/m²/K⁴ |
| `Δh` | specific enthalpy rise of regenerative coolant (§6.5.3) | 0.33–3.34 MJ/kg |
| `a_RT = P_int/σ_t` | screening estimate of acceleration driving RT at the plume/tamper interface. `P_int` is the interface pressure; simulations emit the area-weighted history `a_RT(t)` rather than assuming a constant (§6.7) | 1.9×10⁸ m/s² (Arm D) / 3.9×10⁸ (Arm B) |
| `A_RT` | Atwood number — light fluid pushing dense, so ≈1 here | ≈1 |
| `k_RT = 2π/λ` | RT wavenumber for a feature of wavelength `λ` | 10 cm and 1 m features tabulated |
| `γ_RT = √(A_RT·k_RT·a_RT)` | RT linear growth rate | ×20 at 10 cm over the window |
| `h_b = α_mix·A_RT·a_RT·t²` | one-sided RT bubble depth in the self-similar screening model (§6.7) | — |
| `α_mix` | self-similar RT mix coefficient | 0.02–0.05 |
| **bubble-depth fraction** | `h_b/h_t`, the one-sided RT bubble depth divided by tamper thickness. The quoted total mix-width estimate is a separate, model-dependent quantity and must be labelled as such (§6.7) | 8–21% on Arm D's snow-slug stand-in geometry (§13.12); ≥100% Arm B |
| `α`, `ε` (compaction) | distension ratio and strain in the **P-α / ε-α** porous-ice model (§7.3). *Fourth meaning of `α`* | snow at 70 kg/m³ |

### 0.7 Symbols inherited from the per-collision `f(v)` study

These belong to the *other* study in this repository ([`puffsat_impact_sim_design.md`](puffsat_impact_sim_design.md),
[`CONTEXT.md`](CONTEXT.md)) and appear here only when prior results are cited. **None is a
deliverable of this study.**

| symbol | meaning | where it appears here |
|---|---|---|
| `f` | that study's fudge factor — delivered axial momentum as a fraction of full capture plus perfect elastic bounce | front matter only |
| `e_eff` | effective restitution, `p_rebound/p_in`. **Not** the specific internal energy `e` | §5.3, quoting the 2000× opacity error |
| `eta_capture` | capture *efficiency* — a 2-D/1-D wall-impulse ratio. Distinct from this study's capture *fraction*, which is a free-flight geometric share | §6.6, §12 |
| `Λ` | the prior **tamper multiplier**, tamped impulse over bare impulse, with criterion `Λ > 1 + τ_t`. Superseded by the realization fraction, which does not miscount recoil as loss (ADR-0030); still reported once at Rung 1 as a bound comparable to prior work | §10 Rung 1, §12 |

### 0.8 Similar symbols retained for domain conventions

Load-bearing collisions have explicit subscripts. The few remaining similar forms are retained
because they are conventional, but this table makes their scopes explicit:

| symbol | meanings | how to tell them apart |
|---|---|---|
| `τ_t`, `τ_opt` | tamper ratio / optical depth | always keep the subscript |
| `k`, `k_cond`, `k_RT` | slug ratio / thermal conductivity / RT wavenumber | always keep the subscript outside the slug ratio |
| `α_abl`, `α_th`, `α_mix`, `α` in P-α | ablation exponent / thermal diffusivity / RT mix coefficient / porous distension | the compaction-model name retains its literature notation |
| `σ`, `σ_SB` | areal density / Stefan–Boltzmann constant | always keep `SB` on the constant |
| `A`, `A_RT` | geometric area / Atwood number | always keep `RT` on the latter |
| `γ`, `γ_RT` | ratio of specific heats / RT growth rate | always keep `RT` on the latter |
| `h_t`, `h_b` | tamper thickness / one-sided RT bubble depth | their ratio `h_b/h_t` is the bubble-depth fraction |
| `d`, `δ` | plate standoff / dish depth | dish shape is `δ/D`, never `d/D` |
| `e` | specific internal energy (§4) / `e_eff`, restitution from the prior study (§5.3) | `e_eff` always carries its subscript |
| `P` | pressure (throughout) / momentum, in `P_ejecta` (§3.2) | momentum only ever as `P_ejecta` |

For example, §6.5.2 writes the heat-balance expression unambiguously as
`k_cond·ΔT·√(t/(π α_th))` against `2σ_SB T⁴`.

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

**What makes this tractable:** §3 establishes a single thermodynamic-ceiling definition that every
configuration can be scored against, converting a contested multiplier into a measured
**realization fraction**. In the no-interlayer reference comparison of §3.4, the tamper must
achieve **>62.9%** of that ceiling to beat the bare control, while a *perfect* mirror achieves
61.7%. Swept designs use the configuration-specific break-even rule in §3.4. So the reference
question is whether pressure coupling — omitted from every existing model — closes a
1.2-percentage-point gap.

**Vehicle context:** 1000 t vehicle (initial mass); a 200 kg reference encounter at 1–4 Hz gives
1.69 m/s of `Δv_vehicle`
per pulse (0.17–0.69 g) over ~2960 pulses. Encounter mass and cadence are a **free trade at
fixed thrust**, not independently pinned (§4.1).

**The plate-side working hypothesis is simpler than first scoped.** While an intact ablating
layer is active, its effective vaporisation temperature pins the boundary; the closed-form
screen predicts a steel plate near **750–905 K** at 4 Hz. Rung 6 must verify that interface
condition and burn-through margin. The remaining load-bearing plate quantity is
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
   a **ballistic capture fraction of `max[0, (1 − 1/√k)/2]`**, which is 31.2% at
   `k = 7.06` and zero at `k ≤ 1`. This is the pressure-free ballistic model's result;
   hydrodynamic pressure may still produce small nonzero thrust below that threshold (§8).
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
        = βw / (g₀ · C)
```

`m_charged` is **all expended carried mass**. It is evaluated in two passes:

| | `m_charged` | status |
|---|---|---|
| **Pass 1** (Rungs 1–5) | all expended near-field carried mass `m_hydro` + named `m_other`; ablator excluded | **an upper bound on Isp**, always labelled as such |
| **Pass 2** (Rung 6) | Pass 1 **+ ablator** | the final charged-mass number |

Pass 1 deliberately excludes the ablator so that the tamper question — the point of the
study — is isolated from the plate question and can be answered without waiting on Rung 6.
Pass 1 preserves the no-ablator convention used by prior work. Its bare and no-interlayer
reference cases are therefore directly comparable; filled Arm D additionally charges its
interlayer. **No Pass-1 Isp may be quoted without the upper-bound label**, since the ablator's
share is currently unbounded between ~2% and ~60%.

The distinction between `K`, `K_ej`, and `C` matters only when the interlayer or ablator is
present: `K` inventories the mass in the near-field hydrodynamics, `K_ej` inventories the
mass in the ceiling calculation, and `C` inventories the carried mass charged by Isp. They
are equal in the no-ablator reference cases, which is why earlier algebra could use one letter.

**The projectile is not charged here, and the study reports a second metric instead of
pricing it.** The vehicle spends two distinct resources per unit of Δv: mass it *carries*
(charged by the rocket equation) and *projectiles*, which it does not carry but which
someone had to manufacture, launch, and aim. These are not convertible without a
program-level economic model that lies outside this study, so no combined figure of merit is
constructed. Both are reported as curves:

| metric | definition | meaning |
|---|---|---|
| **Effective Isp** | `J / (g₀ · m_charged)` | carried-mass economy — the rocket-equation currency |
| **Projectile economy** | `β · w = J / m_projectile` | net vehicle impulse per projectile mass; a velocity-equivalent metric, not a physical velocity |

Here `β = J/(m_projectile·w)` is dimensionless and `J` is the net vehicle impulse after
subtracting the incoming projectile momentum debit. Naming `β·w` as the second metric is the
only change.

**What stays out of scope is only the exchange rate** between projectile consumption and
delivered payload — that needs program economics this study does not have. The mass ratio
itself is a design variable this study *does* optimise; see §3.5.

### 3.2 The ceiling

Vehicle frame, `+z` = thrust direction. Impulse on the vehicle is the momentum brought in
minus the momentum carried out:

```
J = −m_i·w − P_ejecta
```

The ejecta's total energy is at most `½ m_i w²` (all carried material starts at rest in the
vehicle frame), so by Cauchy–Schwarz `|P_ejecta| ≤ √(2 M_ej E)`. Define

```
K    = m_hydro/m_i         hydrodynamic carried mass / projectile mass
     = k(1+τ_t+μ)          when slug, tamper, and interlayer are the only regions
K_ej = (M_ej−m_i)/m_i      nonprojectile ejecta mass / projectile mass
```

Thus the per-projectile-mass and per-pulse ceilings are, respectively:

```
j_max   = J_max/m_i = w·(√(1+K_ej) − 1)
J_max   = m_i·j_max
v_e,max = j_max/K_ej
```

In a closed no-ablator calculation, `K_ej = K`. In Pass 2, `K_ej` also includes ablator mass
actually ejected across the system boundary; `C` separately charges all expended ablator mass. The wall model must
report both so Pass 2 does not silently assume they are identical.

Equality requires **every ejecta element to end up moving −z at one common speed**, which
in turn requires the plate to catch and reverse all plate-bound material. At `K_ej = k` this
recovers prior work's ideal-collimation coefficient exactly, `β_ideal = √(1+k) − 1`, so the
ceiling is a reformulation of the existing model rather than a new one. **The two figures
usually quoted beside it are *ballistic*, not ceiling, results** — the bare-plate optimum at
`k* = 7.060` and 1014 s at `w̄ = 77.28 km/s` both follow from `β_bare`, and the ceiling itself
has no interior optimum in `K_ej` (§3.5). (Rung 0's calculator reproduces the 1014 s figure to
four digits from the same `β_bare`, and locates that model's optimum at `k* = 7.0600`; prior
work's quoted 7.057 is the same optimum to optimiser tolerance on a curve that is flat to
±0.6% over `k` = 6–8.)

### 3.3 Two consequences that define the study

**(a) At the ceiling, the tamper is neutral.** `v_e,max` depends only on `K_ej`, not on how
that ejecta mass splits among slug, interlayer, tamper, and ablator. **The tamper can therefore never be justified as a
momentum multiplier** — only as a *realizability* device. The hydrodynamic comparison metric is
the **realization fraction** `r_real = J/J_max`: what share of the same-ejecta-mass ceiling a
configuration achieves. It feeds the final effective-Isp deliverable.

**(b) The tamper's backward recoil is credited, not lost.** Its `−z` momentum contributes
usefully to the full-system momentum balance. It is not, by itself, equal to plate impulse:
the projectile, slug, interlayer, all escaping material, and the vehicle share that balance.
This reverses the prior framing, in which tamper recoil was treated as the mechanism that
kills the concept.

**Design consequence:** entropy production is one important source of shortfall — energy
thermalised in the tamper must then be re-expanded a second time, isotropically, at a second
efficiency cost. Angular dispersion, ejecta velocity variance, residual internal energy,
radiation, incomplete capture, and wrong-way material can also lower `r_real`. Entropy production
rises steeply with shock strength. **The study's design hypothesis is therefore that the ideal
tamper is a maximally isentropic piston rather than a cold specular mirror** — the ceiling proof
disqualifies the mirror as a *justification* for the tamper, but it does not by itself establish
the piston as the realizable optimum; that is what the sweep tests. Instrument the tamper's
entropy/thermal budget and its centre-of-mass momentum alongside the complete shortfall and
momentum ledgers, not its reflected-speed ratio alone.

### 3.4 The reference number to beat

This subsection is an analytic **reference comparison**, not the acceptance rule for every
swept design. It holds `w = 75 km/s`, `k = 7.06`, `μ = a_abl = 0`, and compares a
`τ_t = 1` candidate with the bare `k = 7.06` control. `β = J/(m_i·w)` is the dimensionless net-vehicle-impulse
coefficient, including the incoming projectile momentum debit. Thus `β·w = J/m_i` is a
velocity-equivalent impulse per projectile mass, **not a physical axial velocity**;
`β_bare` and `β_tamp` are the ballistic closed forms from prior work, re-verified here
(`β_bare(7.06) = 0.90871`, `β_bare(14.12) = 1.68353`):

| configuration | `K` | `β` | % of ceiling | Isp |
|---|---|---|---|---|
| **bare plate, k = 7.06** | 7.06 | 0.9087 | **49.4%** | **984 s** |
| bare plate, k = 14.12 — same mass spent as *slug* | 14.12 | 1.6835 | 58.3% | 912 s |
| **perfect-mirror tamper, τ_t = 1** | 14.12 | 1.7825 | **61.7%** | **965 s** |
| *break-even against the bare plate at k = 7.06* | 14.12 | 1.817 | **62.9%** | 984 s |
| ceiling at `K` = 7.06 / 14.12 | — | 1.839 / 2.889 | 100% | 1992 / 1565 s |

Read this row by row. In this reference comparison, a tamper at `τ_t = 1` **beats** spending the same mass as extra slug
(61.7% vs 58.3%), but **loses** to not spending the mass at all (61.7% vs the 62.9% needed)
— *and that is with a perfect mirror*. Any finite-mass mirror does worse. Adding the
free-plate elastic lower bound gives the honest bracket of **591–965 s against 984 s**.

**Therefore this reference tamper pays only if pressure coupling lifts it above 62.9%.** Every prior
model of this device is ballistic (straight-line gas elements, `ρv²` only, no `P` term), and
at `k ≈ 7` the fireball's sound speed and its recoil velocity are equal to within 5% — the
worst possible regime in which to omit pressure. That is the gap this study exists to close.

For any other candidate or for Pass 2, break-even is configuration-specific:

```
β_candidate/C_candidate > β_reference/C_reference
r_real_candidate > [β_reference·C_candidate]
                    / [C_reference·(√(1+K_ej_candidate)−1)]
```

Its required realization fraction follows from that candidate's `K_ej` and `C`; 62.9% must
not be applied as a universal gate. The same-mass extra-slug competitor uses the same `C` as
the candidate, including interlayer and configuration-dependent ablator mass.

---

### 3.5 The mass ratio is a design variable, and the Isp optimum is not where it settles

`K` — near-field carried mass per projectile kg — is a primary design output of this study,
not an inherited input. Where total encounter mass is held at 200 kg,
`m_i = 200/(1+K)`; that equality includes the projectile plus near-field carried mass and
does not include the Pass-2 ablator unless the sweep explicitly holds total charged mass fixed.

**There is a genuine bare-ballistic optimum.** `k*_bare ≈ 7.06` for a bare plate is not an economic
result: it is the competition between ballistic capture fraction, which rises with `k`, and exhaust
velocity, which falls with it. Pure geometry, entirely in scope.

**But it is very flat, and the quantity it trades against is steep.** Three exact relations
govern the trade — note that the first two contain **no encounter mass and no cadence**, so they
hold under any encounter-size/rate schedule at fixed dimensionless design:

```
energy per unit impulse         E/J = w/(2β)     [J per N·s]
projectiles per unit impulse         = 1/(β·w)
effective Isp                        = βw/(g₀C)
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
| 32 | 3.068 | 733 s | −25.5% | 12.2 kJ/N·s | −70% |

**What sets the design point.** Plate *temperature* does not — on the §6.5.2 screen it is
self-limiting at the ablator's vaporisation point regardless of `K` or cadence. What does is the
**ablator's share of the Isp denominator**, since ablation is sub-linear in fluence and so
falls with `β`. The denominator is therefore

```
m_charged / J  ∝  K/β  +  (ablator term)/β^α_abl    with α_abl < 1
```

For the **bare ballistic family**, the first term is minimised at `k*_bare ≈ 7.06`; the
ablator term may shift that family's Pass-2 optimum upward. A tamped or interlayered family
has a different `β(K, τ_t, μ, geometry, …)` and must locate its own Pass-1 and Pass-2 optima
from the computed curves. There is no ceiling optimum at 7.06: at the ideal ceiling,
`v_e,max = w/(√(1+K_ej)+1)` decreases monotonically with `K_ej`.

**Consequence for the two passes.** On **Pass 1**, the bare control should reproduce
`k*_bare ≈ 7.06`; tamped and interlayered designs locate their own optima. On **Pass 2** every
family is re-optimised with its measured ablator term. The `K` grid must therefore span at
least 6–32 rather than clustering around the bare-control answer.

**Consequence — the tamper's case gets harder.** The tamper and extra slug are two routes to
the same end: both raise `β`. §3.4 shows extra slug reaches 58.3% of ceiling against the
tamper's 61.7%, and extra slug additionally carries no RT risk, no interlayer, no porosity
model, and no assembly complexity. **The tamper must beat a rival that is strictly simpler.**

---

## 4. Analytic no-interlayer reference case

`w = 75 km/s`, `k = 7.06`, `μ = 0`, total encounter mass per pulse
`m_i + m_s + m_t = 200 kg`. This table supports the analytic comparisons in §3.4; it is not
the filled Arm D mass ledger.

| quantity | `τ_t = 0` | `τ_t = 1` |
|---|---|---|
| projectile `m_i` | 24.81 kg | 13.23 kg |
| slug | 175.2 kg | 93.4 kg |
| tamper | — | 93.4 kg |
| interlayer | — | — |
| **incoming projectile energy `½ m_i w²`** | **69.8 GJ** | **37.2 GJ** |
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

**Vehicle mass is pinned at 1000 t** — the *initial* (wet) mass, from which the departure-burn
row below is derived, and the plate is 5% of it. **Encounter mass and cadence are
not independently pinned** — they are a free trade at fixed thrust. A 200 kg encounter is the
reference point, not a constraint; configuration-specific carried mass is `m_charged`, not 200 kg:

| quantity (200 kg bare reference encounter) | value | derivation |
|---|---|---|
| impulse per pulse (bare, `K = 7.06`, `τ_t = μ = 0`) | 1.69×10⁶ N·s | `β_bare · m_i · w` |
| **vehicle `Δv_vehicle` per pulse** | **1.69 m/s** | modest, as intended |
| **vehicle acceleration** | **0.17 g at 1 Hz / 0.69 g at 4 Hz** | |
| total encounter mass flow | 200 kg/s / **800 kg/s** | includes externally supplied projectiles |
| carried-mass flow, bare control | 175 kg/s / **701 kg/s** | `m_s` only |
| departure burn (Δv ≈ 7.06 km/s at Isp 984 s) | **~740 s, ~2960 pulses at 4 Hz** | mass ratio 2.08 on 1000 t *initial* mass → 519 t charged, at 175.2 kg/pulse |
| plate recoil per pulse (50 t plate) | 34 m/s | 0.17% of gas speed |

Two cross-checks fall out. Both flows land just below this project's independently-derived
carried-mass requirement (8.85 kg/s for a 10 t vehicle, i.e. 885 kg/s at 1000 t): the
bare-control carried-mass flow is ~21% under it and the total encounter flow ~10% under. They
are kept separate because the projectile is not carried, so the comparison is an
order-of-magnitude agreement rather than a match. The 34 m/s plate recoil against ~20 km/s gas
confirms the **rigid-wall assumption with a factor of ~600 in hand**.

**The encounter-mass / cadence trade is genuinely free, and average heat load is invariant under
it.** Here encounter mass means `m_enc`, scaled at fixed dimensionless mass ratios and geometry.
Since `E/J = w/(2β)` contains neither encounter mass nor cadence (§3.5), average thermal
power is `Thrust · w/(2β)` — **doubling encounter mass and halving cadence changes it not at
all.** Thrust, and therefore acceleration and gravity loss, are likewise held. What the
trade *does* buy:

- **Less total ablator**, because ablation is sub-linear in fluence — concentrating the same
  energy into fewer, larger encounters ablates less overall, and ablator mass is charged (§3.1).
- **Better optical depth**: `Σ` at the plate scales with encounter mass, moving further from the
  `τ_opt ~ 1` regime where flux-limited diffusion is weakest and where this project has
  previously suffered a 2000× opacity error (§5.3).
- **Scale-invariance of the tamper physics**: both `Θ` (§6.1) and the RT bubble-depth fraction
  (§6.7) are ratios and do not move, so the trade cannot rescue or break the tamper.
- Cost: higher per-pulse fluence means deeper ablation per pulse, so the layer must be
  thicker — though the longer interval gives more time to lay it down.

**Cadence is set by gravity losses alone.** Prior finite-burn integration on this project
gives departure gravity losses of 2.2% at 1 g, 5.1% at 0.5 g, 9.3% at 0.25 g, so higher
acceleration is preferred — and at fixed thrust that is a statement about encounter mass × rate,
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

**No existing table covers the required material domain.** The baseline
`data/tables/water.json` spans ρ 0.01–1198 kg/m³ and T 300–60,000 K, but its high-density
states do not constitute a validated solid-ice impact EOS. The prior 69 km/s table
`water_jupiter.json` extends T to 1.2×10⁶ K and low density to 10⁻⁴ kg/m³, but stops at
30 kg/m³ and is a chemical-equilibrium vapor/plasma EOS. It cannot represent 917 kg/m³ ice,
ice Hugoniots, phase boundaries under shock compression, projectile penetration, or porous
compaction. It remains useful only after material has entered its validated fluid/plasma regime.

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
water uncertain across ~two decades (10–1000 m²/kg), `τ_opt = κΣ` spans **0.63 to 63** — it
**straddles `τ_opt ~ 1`**, which is exactly where flux-limited diffusion is weakest. Prior
experience on this project is directly relevant: at 69 km/s an interim Kramers opacity ran
~2000× low at stagnation and falsely predicted `τ_opt ~ 1`, moving `e_eff` from 0.42 to 0.65
once real opacities were used. **Real tabulated opacity is a requirement here, not a
refinement.**

---

## 6. Findings that set the requirements

All derived from scratch; estimates are flagged as such.

### 6.1 The transit ratio — the tamper's binding timing constraint

A tamper functions only if it can communicate its inertia across its own thickness before
the event ends. Let `Θ = t_disassembly / t_shock-transit`. With
`t_dis = r_slug/u`, tamper thickness `h_t = m_t/(2πs²ρ_t)`, and `m_s = (4/3)πr³ρ_s`:

```
Θ = (3/2) · (ρ_tamper/ρ_slug) · (U_shock/u) · (s/r_slug)² / τ_t
```

`U_shock ≈ 12 km/s` (ice Hugoniot at ~70 GPa — **estimated**, worth confirming). The
working gate `Θ ≳ 3` requires about three traversals within `t_dis`; it is a design screen,
not a derived sharp threshold. At `τ_t = 1`:

| configuration | slug radius | tamper thickness | transit | disassembly | **Θ** |
|---|---|---|---|---|---|
| snow slug (70 kg/m³), contact ice tamper | 0.683 m | 3.5 cm | 2.9 µs | 27.6 µs | **9.5** ✓ |
| ice slug (917), contact ice tamper | 0.290 m | 19.3 cm | 16.1 µs | 11.7 µs | **0.7** ✗ |
| ice slug, tamper at 1 m standoff | 0.290 m | 1.6 cm | 1.4 µs | 11.7 µs | **8.7** ✓ |

**The tamper/slug density ratio is the dominant design parameter — for a geometric timing
reason, not an acoustic one.** It sets tamper thickness at fixed mass. There is no hard cap
on `τ_t` near 1: with a snow slug, `Θ > 3` holds out to `τ_t ≈ 3`.

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

Three clocks govern the tamper, and **temperature is in none of them**. At `τ_t = 1`, `s = 1 m`,
`σ_t = 14.9 kg/m²`:

| clock | mechanism | value |
|---|---|---|
| slug disassembly | the event it must act on ends | 11.7 µs |
| **recoil** `σ_t·Δv_t/P` | it comoves with the plume; no further momentum transfers | **~45 µs** |
| lateral spread `R/c_s` | area grows, `σ` falls | ~200 µs |

It completes its work ~4× before it disperses. The tamper–plume collision deposits
~72 MJ/kg — **26× past sublimation, and ~4× cooler than the fireball** (≈0.5 eV vs 1.2 eV).

Cross-check: `σ_tamper/σ_plume = 1.75` here, giving a 1-D elastic free-plate
reflected/incident ratio of **0.274**, which reproduces the prior analytic lower bound at
`τ_t = 1` to two figures.

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
(`σ_abl ∝ Φ^α_abl`, `α_abl < 1`), so the tax scales as `d^{2(1−α_abl)}` — **it grows with standoff.**
Standoff therefore trades a peak pressure that is not binding for an ablator mass that is
charged in the denominator.

**The plate has no thermal sink except ablation.** Over a 750 µs residence, against an
incident fluence of **≈12.6 MJ/m²** (plate-area-averaged at `R = 15 m`, `d = 10 m`; an
order-of-magnitude reference, not a derived capture-energy result — Rung 4's flux map replaces
it, and the local peak near the vertex will exceed it):

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

| basis | depth/pulse | mass at `R = 15 m` | illustrative share of a 200 kg charged-mass basis |
|---|---|---|---|
| prior study's measured 16 km/s ablation fraction (3.7–8.9% of its incident gas-pulse mass), rescaled to this Σ | 7–16 µm | 4.5–10 kg | **2–5%** |
| same, with naive `v⁸` scaling to ~20 km/s arrival | 40–190 µm | 25–121 kg | **13–60%** |

**A 27× spread**, and the study's largest single unknown. It is excluded from Pass 1 by
construction (§3.1) precisely so that it cannot contaminate the tamper verdict; Rung 6
collapses it for Pass 2.

**Why excluding the ablator is *expected* to be conservative for the tamper verdict — a
hypothesis, not a result.** An earlier draft of this document claimed the opposite — that the
tamper "induces roughly proportional extra ablation," partly consuming its own gain. That is
true *per pulse* and wrong *per unit impulse* at fixed energy partition, which is what Isp
measures: since `E/J = w/(2β)` (§3.5), raising `β` lowers the total energy delivered per unit
impulse. But `E/J` fixes only the *total* energy per unit impulse, not the ablation it causes.
Ablation also responds to how that energy arrives — the fraction that reaches the plate at all,
its radial and angular distribution, arrival velocity and spectrum, residence time, plate area
and incidence angle, and vapour-curtain formation — and a tamper changes several of those. **So
the expectation is that the tamper lowers ablator mass per unit impulse, but Pass 1's exclusion
is not proven conservative for the *relative* verdict.** Rung 6 tests it by measuring ablator
mass per configuration rather than assuming proportionality.

**And the ablator is expected to degrade the answer without deciding it.** Even the pessimistic branch
(121 kg against a 200 kg pulse) gives a denominator ratio of 200/321 = 0.62, taking a 984 s
Pass-1 figure to ~610 s — still ~1.6× methalox. The ablator therefore moves the absolute
number materially but does not change whether the architecture is worth pursuing, and it is
expected to largely cancel in the *comparative* question of whether the tamper pays — an
expectation Rung 6 checks rather than inherits. That is why it is last (Rung 6) rather than
first.

#### 6.5.1 At 1–4 Hz the ablator becomes a mass-flow problem, not a coatings problem

The 27× uncertainty above translates directly into two qualitatively different vehicles:

| branch | per pulse | **at 4 Hz** | vs the 800 kg/s total encounter flow | replenishment per pulse |
|---|---|---|---|---|
| optimistic | 4.5–10 kg | **18–40 kg/s** | 2–5% | 7–16 µm in <250 ms |
| pessimistic | 25–121 kg | **100–484 kg/s** | **13–60%** | 40–190 µm in <250 ms |

On the pessimistic branch the ablator is a **larger consumable stream than anything else on
the vehicle except the propellant itself**, and the surface-renewal system must lay down a
~0.2 mm film over 707 m² between pulses. Since the ablator is charged in the denominator
(§3.1), this is what Rung 6 resolves.

**Sub-linearity makes pulse size a lever on this.** Because ablation goes as `Φ^α_abl` with
`α_abl < 1`, concentrating the same total energy into fewer, larger encounters ablates *less* in
total. Encounter mass and cadence trade freely at fixed thrust (§4.1), so **larger encounters at
proportionally lower cadence is a direct reduction in ablator mass** — one of the few levers
that improves the denominator without touching the physics of the tamper.

#### 6.5.2 Inter-pulse balance — on this screen the plate self-limits, and steel is fine

**Everything in this subsection is a closed-form screen, not a computed result**, and is
carried as a hypothesis until Rung 6 tests it. It assumes an unbroken ablating layer holds
`T_abl` for the whole pulse, a single lumped heat partition, and a cyclic equilibrium reached
without oil pyrolysis or vapour-layer behaviour that the closed form does not represent. Rung 6
computes the sustained-inflow ablating-wall problem and verifies (or falsifies) the interface
condition, the heat partition, the burn-through margin, and the cyclic equilibrium. Until it
does, "steel plus oil is settled" and the 750–905 K equilibrium are **testable hypotheses**,
not results — and the §6.5.4 baseline rests on them.

Ablation is mass-transfer cooling: its enthalpy leaves with the vapour, so the plate must
reject only what **soaks in**. The decisive point is that **an ablating surface is a
temperature-pinned boundary.** While it is ablating it sits at the ablator's vaporisation
temperature `T_abl`; the substrate beneath can never exceed that, because as the bulk
approaches `T_abl` the driving ΔT — and with it the conducted flux — goes to zero. The plate
temperature is therefore bounded by the ablator, not by the plume.

Solving the steady balance for a steel plate — conducted in per pulse
`≈ k_cond·ΔT·√(t/(π α_th))` (`k_cond = 45 W/m/K`, `α_th = 1.2×10⁻⁵ m²/s`, `t = 750 µs`)
against `2σ_SB T⁴` radiated from
both faces, at 4 Hz:

| ablator surface temperature `T_abl` | **plate equilibrium** |
|---|---|
| 800 K | **753 K** |
| 1000 K | **905 K** |

Steel is usable to ~1000–1200 K. **On these assumptions it clears this comfortably, at 4 Hz,
with no material escalation and no active cooling.** This is why Orion's steel plate worked,
and it is cadence-independent in the same way §4.1's heat load is: raising cadence raises both the
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

#### 6.5.4 Plate material: steel plus a thin oil ablator, and — if §6.5.2 holds — the question largely dissolves

Every candidate criterion turns out to be non-binding:

| criterion | status |
|---|---|
| shock (2 MPa gas, ~3 MPa inertial at 4600 g) | non-binding — ~200× margin for any structural material |
| spall | non-binding at 2 MPa; the rule that forbade voids was derived at 400 MPa–2 GPa |
| bending / whole-plate rigidity | passes trivially (first mode 10–100 ms vs a 750 µs pulse) |
| **steady-state temperature** | **screened as bounded by `T_abl` ≈ 750–905 K — steel clears it; a hypothesis pending Rung 6 (§6.5.2)** |
| atomic-oxygen attack | the renewed oil layer, not the substrate, meets the plume |

**Baseline: steel structure + a thin renewed oil ablator**, which is the Orion configuration
and rests on the self-limiting argument rather than being merely inherited — but that argument
is a closed-form screen (§6.5.2), so the baseline is **provisional on Rung 6** rather than
settled. It is the right thing to build against; it is not yet a verified plate design.

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
It does not address steady-state rejection across ~2960 pulses at 1–4 Hz, which is the half
that is marginal. *Citation caution: this project records that its Orion references are
secondary and the primary source was never read (firewalled), so any number leaned on here
must be checked against the originals.*

### 6.6 Plate shape

**The prior study's foreclosure of a deep dish is conditional on plane-wave incidence and
does not transfer.** A paraboloid focuses parallel→point and collimates point→parallel; it
cannot do both. A plane-wave cloud striking a dish gets its rebound *focused* into a hot
spot in strongly-radiating, optically-thick gas — the reason the deep dish was rejected. A
source **at the focus** is the opposite case, so that mechanism does not foreclose the shape
here. It does not follow that a dish collimates usefully: whether it does is a simulation
result, not a corollary of ray optics (ADR-0032).

*A point-source screen, not a description of the flow.* For an instantaneous ballistic
expansion every element's trajectory is `(V_cm + u·r̂)·t`, so all rays trace back to a **fixed**
origin — CM recession skews the angular distribution but does not move the apparent source. The
blur is the finite disassembly time, `u·t_dis ≈ 0.6 m`, against a focal length of 6–10 m: about
0.04 rad. The real plume is finite-duration, spatially extended, and pressure-steered, so this
bounds the ray-optics geometry only.

The focus-matched shape is `δ/D = R/(8d)` because `F = d`: **0.19** at
(`R = 15 m`, `d = 10 m`) and **0.31** at (`R = 15 m`, `d = 6 m`) — inside the previously
foreclosed band. Here `δ` is dish depth; `d` is source-to-plate standoff.

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
`a_RT = P/σ ≈ 1.9×10⁸ m/s²` and `γ_RT = √(A_RT·k_RT·a_RT)`, over the 27.6 µs
confinement window:

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

**An analytic bound, and it provisionally forecloses Arm B.** Using the self-similar mix growth
`h_b ≈ α_mix·A_RT·a_RT·t²` with `α_mix = 0.02–0.05` and Atwood number ≈ 1:

| geometry | acceleration `a_RT = P/σ` | window | bubble depth `h_b` | tamper thickness `h_t` | **bubble-depth fraction `h_b/h_t`** |
|---|---|---|---|---|---|
| snow slug, contact tamper (`s = 0.683 m`); slug-disassembly window | 1.9×10⁸ m/s² | 27.6 µs | 2.9–7.3 mm | 34.8 mm | **8–21%** |
| **Arm B** — ice slug, vacuum standoff `s = 1 m`; tamper-recoil window | 3.9×10⁸ m/s² | 45 µs | 16–40 mm | 16.2 mm | **100–246% — fully disrupted** |

**Neither row is Arm D as §4.2 defines it** (ice slug → filled interlayer → tamper at
standoff). The first row is the snow-*slug* corner of the original `(slug density, standoff)`
scoping, carried here as Arm D's stand-in; the two rows also use different windows. Since
`h_b ∝ a_RT·t²` the substitution is not neutral and it cuts both ways: an ice slug's 11.7 µs
disassembly window alone drops the fraction to ~3–8% at equal acceleration, while a tamper at
larger `s` is thinner (16.2 mm at `s = 1 m`) and sees the higher `a_RT` that goes with a smaller
`σ_t`. **Arm D's own numbers are an output of Rung 1's measured `a_RT(t)`, not a value this
table supplies** (§13.12) — the conclusions below are stated on the stand-in geometry.

For the snow-slug row, a separate heuristic that adds spike and bubble penetration gives a
**16–63% total mix-width fraction**. That is not the same quantity as `h_b/h_t`, and the map
from either measure into realization fraction is a screening assumption to be calibrated.

**The screen predicts that Arm B's tamper is shredded before it finishes its job**, across the whole plausible `α_mix`
range — a thin sheet rammed by a fast plume is torn apart, while a thicker shell pressed
gently holds together. This is an *independent* mechanism reaching the same verdict as the
entropy argument in §3.3, and together they provisionally foreclose Arm B (§4.2).

*Bound, not prediction: the self-similar law assumes sustained acceleration and a broad
initial perturbation spectrum, whereas the acceleration here decays and the initial spectrum
is a manufacturing property. But a 2.5× overshoot does not survive a factor-of-two
correction.*

**On the stand-in geometry the bound is not decisive, and that is the problem.** A 16–63% total
mix-width estimate moves the realization fraction roughly between "piston" (~62%) and "just
extra slug" (~58%) under the current screening map. That is close enough to the reference 62.9%
comparison to be decision-limiting, and Arm D proper has not been screened at all. Axisymmetric calculations can represent some modal growth
and provide bounds, but not a general 3-D perturbation spectrum or turbulent cascade. See
the uncertainty budget (§13.1), where this is the top-ranked Pass-1 contributor.

### 6.8 Checked and dismissed

- **Pulse-to-pulse interference.** Exhaust leaves at ~20 km/s; the next projectile closes at
  75 km/s. At the 1–4 Hz reference cadence, the next encounter follows after 0.25–1 s. The
  exhaust-column screen must be recomputed at that spacing, not at the obsolete 20 s value,
  before interference is dismissed; projectile areal density 95–421 kg/m² is the comparison
  scale.
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

1. **Near field** — resolve the tamper; run until the flow is ballistic. Output: the
   **time-resolved control-sphere state** — mass, momentum, and energy flux, pressure, and
   material fraction as functions of polar angle and time — plus the entropy budget.
   **A marginal mass/velocity/angle histogram is not sufficient**: stage 2 is not
   free-streaming, so it needs the joint distribution and the pressure that steers the flow,
   and marginals discard exactly those correlations.
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
  way? For the analytic reference candidate (`τ_t = 1`, `μ = a_abl = 0`), 62.9% is the
  reference break-even. Every swept configuration is tested against its own §3.4 condition.
- **Tier 2 — add radiation transport.** Required at the plate, where `τ_opt` straddles 1 (§5.3)
  and where the thermal load is the binding constraint. Flux-limited diffusion with
  Rosseland means in the diffusion coefficient and Planck means in the emission source.
- **Explicitly not required:** material strength (everything is far past melt on µs
  timescales), chemistry beyond dissociation/ionisation, gravity, MHD.

**Geometry:** 2-D axisymmetric `(r, z)` is required for resolved projectile penetration,
jetting, angular tamper coverage, and plate transport. Spherical 1-D is a prescribed-deposition,
full-shell screen of material interfaces, pressure coupling, and acceleration history; it cannot
represent an axial projectile, projectile aspect ratio, angular deposition, or general RT.
Three-dimensional work is gated on the RT uncertainty (§6.7).

### 7.3 EOS and opacity

- **Water/ice:** equilibrium EOS with dissociation and first ionisation over ρ 10⁻⁴–10³
  kg/m³ and T 300–10⁵ K, including a defensible solid/liquid cold curve, shock-compressed
  states, phase-energy reference, and release paths. No current project table meets this;
  `water_jupiter.json` supplies only the lower-density vapor/plasma part of the eventual model.
- **Snow at 70 kg/m³ is not a standard EOS entry.** It must be modelled as porous ice — a
  P-α or ε-α compaction model layered on the ice EOS. Porous compaction is dissipative and
  changes how the projectile's energy is deposited; **it must not be approximated by simply
  rescaling the ice EOS density.** This is the single largest new physics item.
- **Handoff:** define a thermodynamically consistent transition from the dense ice/water model
  into the existing dissociating/ionising vapor-plasma EOS, with one energy zero and no pressure
  or energy discontinuity over the overlap region.
- **Opacity:** real tabulated per-regime opacities, not Kramers (§5.3).
- **Interlayer density is a swept axis, not a fixed choice.** Impedance matching argues for
  roughly the geometric mean of the fireball and tamper densities — **~400 kg/m³, packed
  slush rather than snow at 70.** Snow may be too light: shocked to very high velocity, it
  would arrive at the tamper as a thin fast sheet and ram rather than press, defeating the
  purpose of Arm D.

### 7.4 Initial and boundary conditions

- **Near field:** the 1-D spherical screen initializes slug, interlayer, and tamper with a
  bracketed prescribed-deposition profile; it cannot contain a directional projectile. The
  2-D axisymmetric near field initializes projectile, slug, interlayer, and tamper and resolves
  projectile energy deposition (§6.2). The prescribed-deposition mode remains the control that
  isolates deposition error from tamper physics.
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

1. **Axial wall impulse vs time.** At the stationary inviscid plate, impermeability removes
   advective mass flux, so wall force is pressure traction:
   `J_wall,z = ∫∫ P(n̂_wall·ẑ) dA dt`, with `n̂_wall` defined in §0.3. Separately report axial
   momentum crossing a named upstream control surface as
   `∫∫[ρv_z(v·n̂_cs)+P(n̂_cs·ẑ)] dA dt`, where `n̂_cs` points out of the fluid control volume.
   The familiar `ρv_z²+P` is only the planar special case; it is not the curved-wall integrand. The difference
   between upstream ram flux and actual wall pressure impulse is a headline diagnostic.
2. **Realization fraction** — `r_real = J/J_max = β/(√(1+K_ej)−1)`, using net vehicle impulse
   and the same configuration's ejecta-mass ceiling. A comparison metric feeding Isp.
3. **Complete ceiling-shortfall ledger** — entropy production, angular dispersion, ejecta
   velocity variance, residual internal energy, radiation, incomplete capture, and wrong-way
   material, to the extent each can be separated.
4. **Tamper CM momentum vs time** — a credited contribution to the full-system momentum
   audit, not a standalone proxy for plate impulse.
5. **Fraction of tamper mass ending up plate-directed** vs entrained (the §6.7 mixing
   question, insofar as an axisymmetric code can bound it).
6. **Time-resolved control-sphere state** — mass, momentum, and energy flux, pressure, and
   material fraction against polar angle and time (§7.1). Its angular mass-flux marginal is
   directly comparable to the analytic isotropic assumption and shows whether the tamper
   *turns* the plume or merely *stops* it, but the **joint, time-resolved** form is what the
   transport stage consumes.
7. **Flux and pressure maps over the plate surface** — feed the taper calculations and the
   ablation model.
8. **Mass, momentum and energy audits at every dump**, closing to <1%. The entire argument
   is an energy-to-momentum conversion, so this is the minimum bar.

---

## 8. Validation and acceptance criteria

No result is trusted before these pass. Written **before** the code, per project convention.

**Analytic anchors specific to this device:**

- [ ] **Ceiling.** In the perfectly-collimated limit,
      `j = J/m_i → w(√(1+K_ej) − 1)`.
- [ ] **Bare ballistic limit.** At artificially low density (collisionless expansion),
      `β_bare(7.06) → 0.9087` and the ballistic capture fraction → 0.3118.
- [ ] **`k ≤ 1` ballistic zero-capture floor.** Ballistically nothing reaches the plate below `k = 1`.
      A hydrocode should show *small but non-zero* thrust from pressure — **how much is
      itself a useful measurement** of how wrong the ballistic model is, and it is cheap.
- [ ] **Free-plate elastic bound.** Reflected/incident → 0.274 at `τ_t = 1` in the ballistic
      limit (§6.3).
- [ ] **`k → 0` degeneracy.** Zero net impulse as slug mass vanishes. A code that produces
      thrust at `k = 0` is wrong.

**Standard verification:**

- [ ] Sod shock tube; Sedov blast; Noh implosion (the axisymmetric source term).
- [ ] Marshak wave for the flux-limited diffusion.
- [ ] Smooth-flow order-of-accuracy test at the scheme's formal rate.
- [ ] Two-material shock-tube with an exact interface solution (**new** — the existing
      kernels are single-material).
- [ ] Dense-ice and porous-ice Hugoniots plus release paths against published data;
      compaction energy and the handoff to the vapor/plasma EOS must close (**new**).
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
| `crates/hydro1d` | 1-D **planar** Lagrangian rad-hydro: staggered mesh + artificial viscosity, pluggable EOS, flux-limited diffusion (Levermore–Pomraning), coupled gas+solid conduction, two-phase condensation, **ablating wall** (`Q*` surface balance, blowing factor, vapour shield). Single material. | **Extend** to spherical geometry + multi-material prescribed-deposition screening (Rung 1). Its ablating wall is used in Rung 6. |
| `crates/euler2d` | 2-D axisymmetric Euler: HLLC Godunov, MUSCL-Hancock 2nd order, Strang splitting, conservative cylindrical source, **ghost-cell immersed boundary** for a shaped plate. Ideal gas, single material. | **Extend** with the table EOS and multi-material (Rung 3); its immersed boundary is reused for the plate (Rung 4). |
| `crates/tables` | Shared JSON table loader. | Unchanged. |
| `crates/sweep` | Rayon-parallel sweep driver, JSONL output. | Extend with new sweep modes. |
| `python/puffsat` | EOS/opacity table generation (CoolProp + analytic Saha/CEA-style, TOPS/OPLIB overlay), frontier extraction, analysis, plotting. | **Extend** with the porosity model; reuse the analysis and plotting pipeline. |
| `data/tables/water_jupiter.json` | Equilibrium vapor/plasma table, T to 1.2×10⁶ K, ρ 10⁻⁴–30 kg/m³, real gray opacities. | Reuse in its validated lower-density regime after the new dense-material handoff; not a solid/porous impact EOS (§5.1, §7.3). |
| Build | Top-level `Makefile` over `cargo` + `uv`; JSON tables in, JSONL results out. | Unchanged. |

**Verified working and directly transferable within their existing regimes:** the equilibrium water vapor/plasma EOS, the
FLD implementation, the ablating-wall model, the immersed-boundary plate, the sweep/analysis
plumbing, and the validation discipline.

**Genuinely missing:** a dense solid/liquid ice EOS and its handoff to the vapor/plasma table
(§7.3 — the largest single item), porous-ice compaction, multi-material tracking, spherical 1-D
geometry, a table EOS inside the 2-D kernel, spherical-inflow boundary conditions, and a
sustained-feed integration window.

---

## 10. Work plan

**Definition of done.** The study ends when it returns **effective Isp for Arm D within a
stated bracket**, together with a pay/don't-pay verdict from configuration-specific §3.4
break-even comparisons, the projectile-economy curve (§3.1), and the optimal mass ratio `K`
for each compared configuration family (§3.5).
The result is explicitly **single-code and preliminary**; the independent hydrocode
cross-check is deferred and named as the outstanding validation. A bracket that straddles
its applicable break-even frontier is a legitimate outcome and is reported as such — the
deliverable is the bracket, not a verdict forced past the evidence.

**And a narrow win is not a verdict.** If the best tamped and best untamped designs end up
separated by only a few realization-percentage points — the scale of the §13.1 residuals, and
of the 1.2-point gap the whole study is chasing — then the *sign* of that difference is a
single-code artifact until an independent hydrocode or methodologically distinct spot-check
reproduces it. In that case the study reports the magnitude, the bracket, and the named
outstanding cross-check, and **does not report a design conclusion**. Only a separation
comfortably outside the combined uncertainty is quotable as pay/don't-pay on one code.

**Two passes (§3.1), executed in strict order.** Rung 0 is a cold-path prerequisite to both:
it fixes the ledger and the closed-form references everything else is scored against. Rungs 1–5
then constitute **Pass 1**, which excludes ablator mass and returns an *upper bound* on Isp.
Its bare control should recover
`k*_bare ≈ 7.06`; each tamped/interlayered family locates its own optimum. Rung 6 adds the
measured ablator (**Pass 2**) and re-optimises every family.

**Uncertainty reduction comes before cost refinement.** RT is *inside* Pass 1's bracket and
is its widest contributor; the ablator is *outside* Pass 1 entirely and moves the absolute
number without changing the tamper verdict (§6.5). So RT is Rung 2 — immediately after the
run that produces the acceleration history it needs — and the ablator is last.

**Rungs 1 + 2 give a preliminary prescribed-deposition screen**, combining the 1-D full-shell
result with analytic capture geometry and an RT bound. It may foreclose a clearly poor design,
but cannot establish a positive final verdict before 2-D resolves deposition and angular flow.

### Rung 0 — the analytic reference ledger (no kernel work)

A single cold-path calculator that owns **every closed-form number this document quotes**, so
the PRD, the ADRs, and the eventual analysis cannot drift apart. It is a day's work and it is
the only thing downstream rungs are permitted to quote closed-form figures from.

- [ ] Implement the mass ledger (`k`, `τ_t`, `μ`, `K`, `K_ej`, `a_abl`, `a_other`, `C`) and the
      per-pulse masses at a stated `m_enc`, with the interlayer and ablator conventions of §0.1
      applied uniformly.
- [ ] Implement `β_bare(k)`, `β_ideal`, the ballistic capture fraction, the ceiling
      `j_max`/`v_e,max`, realization fraction, effective Isp, projectile economy, and energy
      per unit impulse.
- [ ] Implement the **configuration-specific break-even** comparison of §3.4 — candidate
      against both the bare reference and the same-charged-mass extra-slug competitor — as a
      function, not a constant. 62.9% is one evaluation of it, not a gate.
- [ ] **Anchor:** reproduce prior work's 1014 s at `w̄ = 77.28 km/s` and `β_bare(7.06) = 0.9087`,
      confirming the inherited ballistic model before anything is built on it.
- [ ] Regenerate the §3.4, §3.5, and §4 tables from it and reconcile any residual disagreement
      with a quoted figure, rather than carrying both.

### Rung 1 — material qualification and 1-D spherical prescribed-deposition screen

*Why 1-D is useful here: lateral communication across the fireball takes `r/c_s ≈ 69 µs`,
against a 45 µs confinement time, so a full-shell model can cheaply screen radial loading.
This does not make an axial projectile spherical or permit 1-D to resolve penetration or RT.*

- [ ] Add spherical geometry to the 1-D Lagrangian kernel.
- [ ] Add a per-cell material index (projectile / slug / interlayer / tamper).
- [ ] Select or construct the dense ice/water EOS, define its handoff to the vapor/plasma EOS,
      and build the porous-ice (P-α or ε-α) compaction model (§7.3).
- [ ] **Material-qualification gate, before any screening number is quoted.** Dense-ice and
      porous-ice Hugoniots and release paths reproduced against published data; compaction
      energy accounted; the dense→vapor/plasma handoff continuous in pressure and energy on one
      energy zero (§8). This gates the *screen itself*, not merely Rung 3's penetration work:
      the slug and tamper start as solid ice at 917 kg/m³ and the interlayer as porous ice, so a
      screen run on an unqualified material model measures the model, not the device.
- [ ] Use bracketed prescribed-deposition profiles in 1-D; resolved axial projectile
      penetration moves to Rung 3.
- [ ] Write the §8 acceptance tests **first**; make them pass.
- [ ] **Pin Arm D's geometry and mass ledger** — `r_slug`, `s`, `μ`, `ρ_interlayer` — and emit
      its own `Θ`, `h_t`, and `a_RT(t)`. Every Arm D screening number quoted so far belongs to
      the snow-slug stand-in geometry, not to Arm D as §4.2 defines it (§6.7, §13.12).
- [ ] Sweep **Arm D** over `(ρ_interlayer, μ, standoff, τ_t)`; run **Arm B as a control** and
      produce its loading and acceleration history for the separate RT treatment (§6.7).
- [ ] Sweep `K` over at least 6–32 and **locate each Pass-1 design point**, reporting Isp,
      projectile economy, and energy per unit impulse together (§3.5). Require the bare
      control to recover `k*_bare ≈ 7.06`; measure how each tamped/interlayered optimum moves.
- [ ] Sweep deposited-energy fraction and radial deposition profile over the §6.2 analytic
      bracket; do not label these as resolved projectile-geometry results.
- [ ] Report: realization fraction, tamper entropy budget, tamper CM momentum, `Θ` scaling
      confirmed, and `Λ` for a full-shell tamper as a bound.
- [ ] **Emit the acceleration history `a_RT(t)` at the plume/tamper interface** — the input
      Rung 2 needs, and free from this run.
- [ ] **Gate:** compare each candidate's `β/C` with the appropriate bare and same-charged-mass
      extra-slug references. If its bracket cannot approach either break-even frontier, the
      tamper is dead and the remaining rungs are descoped to a bare-plate Isp confirmation.

### Rung 2 — RT treatment for Arm D (the widest Pass-1 uncertainty)

Arm D's total mix-width fraction is currently estimated only at **16–63%**, mapping under the
screening assumption to ~58–62% of ceiling near the 62.9% reference threshold (§6.7).
**Until this narrows, no amount of downstream work can produce a verdict**, so it runs before
any 2-D effort. Escalate only as far as needed:

- [ ] **Tightened analytic bound.** Integrate RT growth against Rung 1's *actual* decaying
      `a_RT(t)` instead of the constant-acceleration self-similar law, and against a stated
      initial perturbation spectrum rather than an assumed broad one. Cheap; may alone be
      enough, since the current bound is deliberately crude in both respects.
- [ ] **Gate:** if the tightened bound no longer straddles the applicable configuration-specific
      break-even frontier, stop here.
- [ ] **Mix model**, if it still straddles — a buoyancy-drag or K-L-style mix width at the
      interface in the 1-D kernel, calibrated against published RT mixing data.
- [ ] **Resolved spot-check**, only if the mix model is inconclusive: a 2-D or 3-D
      interface-resolved run at the single design point. Expensive and last.
- [ ] Report the RT contribution as an explicit ± on the realization fraction, whatever
      level it is settled at.

### Rung 3 — 2-D axisymmetric multi-material

- [ ] Port the table EOS into the 2-D kernel.
- [ ] Add multi-material tracking.
- [ ] Reproduce Rung 1 in the full-shell prescribed-deposition limit (**gate**).
- [ ] Add resolved projectile penetration and energy deposition; sweep projectile areal
      density, aspect ratio, and bulk density against the vertex-hole and deliverability
      constraints (§6.2). Report deposited-energy fraction as a first-class output.
- [ ] Sweep tamper **angular coverage** at fixed mass (§6.7 — coverage, not curvature).
- [ ] Emit the **time-resolved control-sphere state** (mass/momentum/energy flux, pressure,
      material fraction vs angle and time), not a marginal histogram (§7.1, §7.6.6).

### Rung 4 — transport and plate

- [ ] Spherical-inflow boundary condition from the control-sphere state.
- [ ] **Handoff verification (gate):** the reconstructed inflow must reproduce the near-field
      run's mass, momentum, and energy fluxes and its pressure–angle correlation across the
      control surface. An overlap-region comparison against a single-domain run at reduced
      resolution is the check; a handoff that only matches marginals fails it.
- [ ] Coarse far-field continuation to the plate (not free-streaming, §7.1).
- [ ] Replace the decay-based integration window with a sustained-feed criterion (§6.4).
- [ ] Sweep plate radius `R` (mass-ceiling constrained), standoff `d`, and shape: flat plus
      the paraboloid family including `δ/D` up to ~0.35.
- [ ] Emit plate flux and pressure maps.
- [ ] **Gate (D6):** one coupled rad-hydro spot-check at the `τ_opt ~ 1` corner, showing that
      radiation and thermal loss cannot move axial impulse or the velocity/angle distribution by
      an amount comparable with the tamper margin being decided. If it can, the one-way
      1-D-thermophysics × 2-D-geometry factorization does not license the sweep (§5.3, §12).

### Rung 5 — Pass 1 deliverable

**This is the study's answer to the question it was created to ask.**

- [ ] Frozen-recombination bracket at 57 eV/molecule (§13.1) — both ends.
- [ ] Fold interlayer-density and projectile-deposition sensitivity into the bracket.
- [ ] Both taper calculations (closed form, §6.6).
- [ ] Assemble **effective Isp for Arm D within a stated bracket**, on carried mass
      *excluding* the ablator — explicitly labelled an upper bound (§3.1).
- [ ] Report each **realization fraction and configuration-specific break-even frontier** and
      the resulting pay/don't-pay verdict on the tamper, or an honest "cannot distinguish"
      if a bracket straddles its frontier. Retain 62.9% as the labelled §3.4 reference annotation.
- [ ] Report the projectile-economy curve and the Pass-1 optimum `K` (§3.5).

### Rung 6 — ablator and Pass 2 (last)

The ablator is expected to move the absolute Isp without changing the tamper verdict, with its
exclusion expected-conservative but unproven (§6.5). So it is refined only once the verdict
exists — and this rung tests that expectation rather than assuming it.

- [ ] Run the existing 1-D ablating-wall model at this study's areal density, arrival
      velocity, and 750 µs residence time to collapse the **27× uncertainty** (§6.5).
- [ ] Flux map → ablating-wall model at representative radii → ablated mass per
      configuration, with its `Q*` and opacity bracket, and as a mass flow (§6.5.1).
- [ ] Report **burn-through margin** — the one abrupt failure mode, and what actually binds
      on the plate side (§6.5.2).
- [ ] Sweep **encounter mass `m_enc`** at fixed ratios to quantify the sub-linearity gain: larger encounters at
      proportionally lower cadence should cut total ablator mass at fixed thrust (§4.1).
- [ ] Confirm the **self-limiting plate equilibrium** (§6.5.2) against the coupled model, so
      the steel baseline rests on a computed interface condition rather than a closed form.
- [ ] **Fold the ablator into the denominator** and re-derive Isp — the Pass 1 → Pass 2 step
      (§3.1). Report both, with Pass 1 still labelled an upper bound.
- [ ] **Relocate every configuration family's optimum `K`** with the ablator term present (§3.5).
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
| D1 | **Frame the tamper as an isentropic piston, not a cold mirror**, and instrument entropy and CM momentum accordingly. A framing and instrumentation decision, plus a design hypothesis — not a result. | "Reflects all energy, absorbs none, still recoils" is self-contradictory at finite mass, so the mirror is disqualified as a justification. Under §3.3(b) recoil momentum is credited; entropy is one important part of the complete ceiling-shortfall ledger. Whether gentle, early, pressure-mediated loading actually maximises `r_real` is measured, not assumed. |
| D2 | **Arm D (filled interlayer) is the design; Arm B (vacuum gap) is provisionally foreclosed and retained as a control.** | Two independent arguments close Arm B: it is a re-thermalising ram rather than a piston (§3.3), and RT disrupts its tamper completely before it acts (§6.7). Effort concentrates on narrowing Arm D's bracket. |
| D2a | **The mass ratio `K` is optimised separately for each configuration family and pass, against Isp — not against an economic weighting.** | The bare ballistic optimum is real but flat (±0.6% over `k` = 6–8), while tampers and interlayers change `β(K)` and therefore may move it. Ablator cost is configuration-specific and in scope (§3.5). |
| D2b | **Projectiles are not priced into a combined figure of merit.** | Converting projectiles to payload-equivalent needs a program economic model outside this study; importing it would make a physics result move with someone else's assumptions. Projectile economy is reported as `β` (§3.1). |
| D3 | **Plate radius swept; plate mass a ceiling.** | §6.5: `R` is the only lever that breaks the capture-vs-fluence conflict. |
| D4 | **Isp denominator = all expended carried mass, ablator included — but evaluated in two passes**, Pass 1 excluding it as an explicit upper bound. | Isp is a rocket-equation quantity, so anything carried and expended is charged. But the ablator is uncertain by 27× and answering it is a *plate* question; excluding it in Pass 1 decouples the *tamper* question and lets Rungs 1–5 proceed without waiting, with RT narrowed first (§3.1). |
| D4a | **Encounter mass `m_enc` and cadence trade freely at fixed thrust; neither is independently pinned.** | At fixed dimensionless design, `E/J = w/(2β)` contains neither, so average heat load and gravity loss are invariant under the trade — while sub-linear ablation means larger, rarer encounters cut ablator mass (§4.1, §6.5.1). |
| D4b | **Plate material: steel + thin renewed oil, with a thin ceramic hedge — provisional on Rung 6.** | An ablating surface is temperature-pinned, so the plate self-limits near 750–905 K; strength, spall, and bending are all non-binding at 2 MPa. The temperature half of that is a *closed-form screen*, not a computed result, so the decision is a baseline to build against and Rung 6 confirms or falsifies it. Escalation is triggered only by inadequate burn-through margin (§6.5.2, §6.5.4). |
| D5 | **Plate shape sweep = flat + parabola family + both tapers.** Shape itself stays open. | §6.6: the prior foreclosure is conditional on plane-wave incidence; but the parabola's area penalty may exceed its impulse gain. Genuinely two-sided. |
| D6 | **Keep the 1-D-thermophysics × 2-D-geometry factorization**, adding a cold-path flux map — but gate it. | Radiation is local-diffusive and one-way wherever the flow is optically thick, which holds in the near field; at the plate `τ_opt` spans 0.63–63 and straddles 1 (§5.3), so the factorization is *adopted subject to* a coupled sensitivity spot-check at that corner (Rung 4) rather than assumed. A monolithic 2-D rad-hydro across the sweep is not affordable. |
| D7 | **Ablator mass is an emergent cost per configuration, not a specified thickness.** | §6.5: ablation is the plate's only thermal sink, so it is forced by physics. Specifying a thickness the balance does not respect would silently burn through. |
| D8 | **Tamper design variable is angular coverage, not curvature.** | §6.7: RT destroys sub-metre features within the confinement window. |
| D9 | **Build order: dense-material qualification and 1-D spherical screen first, then resolved 2-D.** | Isolates EOS, multi-material, porosity, and prescribed radial loading before adding directional projectile penetration, with known-answer anchors at each step. |

---

## 12. Departures from prior decisions in this repository

Each of these reverses or narrows a recorded decision and requires its own record before the
work lands.

| Prior decision | Transfers? | Why |
|---|---|---|
| 1-D/2-D factorization | **Yes** (D6) | With a cold-path flux map added, and a Rung 4 coupled spot-check gate, since `τ_opt` straddles 1 at the plate |
| Facesheet survivability ladder, peak `≈1.2ρv²` | Yes, but **not binding** | ~200× margin (§6.4) |
| Ablating-wall model | **Yes, and promoted to the critical path** | It is now the largest unknown in the denominator |
| "The ablator is not a pressure device" | Yes, and confirmed | But moot — pressure is not the constraint here |
| Frozen-recombination bracket | Yes, and **worse here** | 57 eV/molecule vs the prior study's regime |
| Rigid-during-pulse gate | Yes, passes trivially | §6.4 |
| `10⁻³`-of-peak integration window | **No** | Wrong for a 750 µs sustained feed (§6.4) |
| RT/RM deferral | **No** | Load-bearing in three places (§6.7) |
| Deep-dish foreclosure | **No** | Conditional on plane-wave incidence (§6.6) |
| Ablator, vehicle scale and cadence held out of scope | **No** | An Isp deliverable pulls all three inside the boundary |
| Inter-pulse plate thermal accumulation excluded as cadence-dependent | **Partly** | Now screened rather than excluded — and on that closed-form screen it does not bind, because the plate self-limits at `T_abl`. Rung 6 computes it (§6.5.2) |
| Projectile geometry treated as an unrecorded external input | **No** | §6.2 makes it a swept design variable with a real interior optimum |
| No voids behind the hot face (spall risk from free-surface reflection) | **No** | Derived at 400 MPa–2 GPa; peak here is ~2 MPa. Not currently exercised, since cooling proved unnecessary (§6.5.3) |
| Carbon-carbon rejected as a hot face (burns in atomic O) | **Weakened** | Exposure is only during the pulse, when renewed oil covers it; between pulses there is no gas. Held in reserve on the escalation path, not adopted (§6.5.4) |
| SiC + Ti as the settled hot-face stack | **Departed** | Selected there for oxidation and per-pulse thermal shock at GPa loads. Here loads are ~2 MPa and the plate self-limits, so **steel + oil is the baseline** and ceramic is a burn-through hedge (§6.5.4) |
| Existing `eta_capture` sweep data | Cross-check only | Incidence geometry differs |

**Records written.** Each clears the project's bar for an architecture decision record —
hard to reverse, surprising without context, and the result of a genuine trade-off:

- [x] **[ADR-0030](docs/adr/0030-tamper-isentropic-piston-not-mirror.md)** — the tamper is an
      isentropic piston, not a mirror; recoil is credited; realization fraction is the
      hydrodynamic comparison metric under one ceiling definition, while effective Isp is the
      deliverable. Replaces the prior `Λ > 1 + τ_t` criterion (§3.3).
- [x] **[ADR-0031](docs/adr/0031-isp-deliverable-pulls-scope-inward.md)** — an Isp deliverable
      pulls ablator mass, vehicle scale, and cadence inside the scope boundary that prior
      scope places outside, and inter-pulse thermal accumulation with them (§3.1).
- [x] **[ADR-0032](docs/adr/0032-deep-dish-foreclosure-is-incidence-conditional.md)** — the
      deep-dish foreclosure is conditional on plane-wave incidence and must be reopened for a
      finite-standoff source (§6.6). ADR-0021 carries a pointer and is otherwise unchanged.
- [x] **[ADR-0033](docs/adr/0033-rt-deferral-does-not-transfer-to-the-tamper.md)** — the
      RT/RM deferral does not transfer; RT is load-bearing here in three places (§6.7).
      ADR-0020 carries a pointer and is otherwise unchanged.

Canonical terms for this study are in [`CONTEXT.md`](CONTEXT.md) under *Language — tamped-nozzle
study*; §0 is the symbol table and §14 the reading-order glossary for this document.

---

## 13. Open questions and assumptions register

State these in any write-up. Each could move the answer.

1. **RT/RM is un-modelled and, for Arm D, wider than the margin.** Mixing decides whether
   the tamper is a piston or entrained payload. The analytic screen (§6.7) provisionally closes
   Arm B but leaves Arm D's stand-in geometry at a **16–63% total mix-width estimate**, mapped
   heuristically to ~58–62% of ceiling near the 62.9% reference threshold — and Arm D proper is
   unscreened (item 12). Axisymmetric calculations can bound
   some modes but not the general 3-D cascade. **Top-ranked uncertainty-reduction target**
   (§13.1, Rung 2 — it now runs *before* the ablator, not after).
2. **Interlayer density may want to be ~400 kg/m³, not 70.** Impedance matching argues for
   packed slush over snow. Cheap to sweep; a real risk to Arm D as currently specified.
3. **~~Vehicle mass~~ — RESOLVED: 1000 t.** The plate is 5% of the vehicle and the prior
   10 t/50 t inconsistency is gone (§4.1).
4. **~~Cadence~~ — RESOLVED on a screen: not thermally limited.** The plate self-limits at the
   ablator's vaporisation temperature — a closed-form screen Rung 6 must confirm
   (§6.5.2) — and average heat load is invariant under the
   encounter-mass/cadence trade (§4.1). Cadence is set by gravity losses alone, and larger encounters
   at lower rate is a *free* variation that also cuts ablator mass. What remains open is the
   shock-absorber stroke this implies, which is out of scope.
5. **~~Projectile geometry~~ — RESOLVED: it is in scope as a swept design variable** (§6.2).
   The residual risk is that resolving hypervelocity penetration is a distinct validation
   domain (§8) and could dominate Rung 3's effort.
6. **~~Plate material~~ — PROVISIONALLY RESOLVED, on a closed-form screen: steel + thin
   renewed oil, with a ceramic hedge.**
   Every candidate criterion is non-binding except steady-state temperature, and that
   self-limits at `T_abl` — which is screened, not computed, so the baseline stands or falls
   with Rung 6 (§6.5.2, §6.5.4). **What replaces it as the open question is
   burn-through margin** — the one abrupt failure mode, since a mid-pulse breach exposes the
   substrate to 50–80 kK plume with no temperature pinning. Rung 6 reports it. Escalation
   materials are named but not invoked.
7. **`K`, `μ`, and `τ_t` grids not pinned.** The bare Pass-1 control should peak near
   `k = 7.06`; tamped/interlayered optima are unknown, so `K` must span at least 6–32 rather
   than clustering around 7. `τ_t` spans 0–2, and `μ` must be set by the interlayer-mass sweep.
   Pass 2 re-optimises after Rung 6 measures the ablator term.
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
12. **Arm D has no screened geometry or mass ledger of its own yet.** §4.2 defines it as ice
    slug → filled interlayer → tamper at standoff, but every quoted Arm D number — the §6.7
    RT screen, the §6.1 `Θ` rows, and the §4 reference masses — belongs either to the
    snow-*slug* corner of the original `(slug density, standoff)` scoping or to a `μ = 0` case.
    Rung 1 must fix Arm D's `(r_slug, s, μ, ρ_interlayer)` and emit its `Θ`, `a_RT(t)`, and mass
    ledger before any Arm D screening number is quoted as such.

### 13.1 Uncertainty budget for the Arm D Isp bracket

The deliverable is a bracket (§10), so its width is a first-class output. Contributors,
ranked by current width, with what collapses each:

| # | contributor | current width | affects | in Pass 1? | collapsed by |
|---|---|---|---|---|---|
| 1 | **RT total mix-width fraction at the plume/tamper interface** | 16–63% heuristic width → ~58–62% of ceiling under the screening map, on the snow-slug stand-in geometry (§13.12) | realization fraction | **yes** | **Rung 2** — bound, mix model, or resolved spot-check |
| 2 | **Frozen recombination** at 57 eV/molecule | this project's largest quantified physics uncertainty, in a *cooler* regime than here | realization fraction | **yes** | Rung 5 bracket, both ends |
| 3 | **Projectile deposition** (34% vs 75% of energy into the slug) | changes which body is the fireball | everything | **yes** | Rung 1 bracket; Rung 3 resolved penetration |
| 4 | **Interlayer density** (70 vs ~400 kg/m³) | may invert Arm D's loading mechanism | realization fraction | **yes** | Rung 1 — cheap sweep axis |
| 5 | **Ablator mass** (vapour-shielding balance) | **27×** → 2–60% of the mass budget | denominator | **no — excluded by construction** | Rung 6 (last) |
| 6 | **Opacity** near `τ_opt ~ 1` at the plate | prior experience: 2000× error from the wrong table | flux map, denominator | partly | the real-opacity table within its validated low-density regime (§9), plus the Rung 4 coupled spot-check (D6) |
| 7 | **`U_shock` in ice** | linear in `Θ` | tamper timing | yes | literature lookup |

**The two-pass split (§3.1) is what makes this tractable.** The single widest contributor —
the ablator, at 27× — sits in the denominator and is *excluded* from Pass 1, so it cannot
contaminate the tamper verdict. Pass 1's bracket is set by contributors 1–4, of which
**RT alone still spans the reference 62.9% comparison closely enough to be decision-limiting**.
Until that is narrowed, Pass 1
cannot distinguish "the tamper pays" from "the tamper is dead," and every interim number must
carry its bracket and its upper-bound label.

---

## 14. Glossary

Named concepts, in reading order. **Symbols are defined in full in §0**; the handful repeated
here are the load-bearing ones, quoted so this table stands alone.

| term | meaning |
|---|---|
| **PuffSat / projectile** | Mass accelerated by external infrastructure and aimed at the vehicle. Not carried, so not charged in Isp. |
| **Pulse / encounter** | One projectile interaction and its associated carried reaction mass. `m_enc = m_i + m_hydro`; “200 kg” denotes this reference encounter mass, not charged carried mass. |
| **Slug** | Carried mass the projectile buries itself in. Vaporised; becomes the fireball. |
| **Tamper** | Dense carried mass beyond the slug that turns part of the away-going fireball around. Fully vaporised; works by inertia (areal density), not by survival. |
| **Interlayer** | Low-density carried mass filling the standoff between slug and tamper in Arm D. Both spacer and pressure-transmitting medium. |
| **Pusher plate** | The permanent structure that receives the impulse. Has a hole at its vertex for the projectile. |
| **`k`** | Slug mass / projectile mass. ≈7 at the bare optimum. A continuous design variable. |
| **`τ_t`** | Tamper mass per slug kg; distinct from optical depth `τ_opt`. |
| **`μ`** | Interlayer mass per slug kg. |
| **`K = m_hydro/m_i`** | All near-field hydrodynamic carried mass per projectile kg; equals `k(1+τ_t+μ)` when slug, tamper, and interlayer are the only material regions. |
| **`K_ej`** | Nonprojectile ejecta mass per projectile kg; the mass used by the thermodynamic ceiling. |
| **`C`** | All expended carried mass per projectile kg; the mass ratio used by effective Isp. |
| **`w`** | Closing speed, 75 km/s. |
| **`β = J/(m_i·w)`** | Dimensionless net-vehicle-impulse coefficient, including the incoming projectile momentum debit. `β·w` has velocity units but is not a physical velocity. |
| **Ballistic capture fraction** | `max[0, (1 − 1/√k)/2]` — share of the pressure-free isotropic ballistic fireball that outruns its own recoil. Its zero at `k ≤ 1` is not a general zero-thrust claim. |
| **Realization fraction** | `r_real = β/(√(1+K_ej)−1)`, the net vehicle impulse divided by the same-ejecta-mass ceiling. A hydrodynamic comparison metric feeding Isp; 62.9% is only the labelled §3.4 reference break-even. |
| **Projectile economy** | Net vehicle impulse per projectile mass, `β·w = J/m_i`. A velocity-equivalent metric, not a physical velocity. Reported alongside effective Isp and not combined with it — see §3.1. |
| **Arm D / Arm B** | Filled-interlayer configuration (the design) / vacuum-gap configuration (provisionally foreclosed, retained as a control). |
| **Transit ratio `Θ`** | Slug disassembly time / tamper shock-transit time. `Θ ≳ 3` is the current screening gate, not a universal sharp threshold. |
| **Arrival window** | Duration over which fireball material reaches the plate, set by its velocity dispersion. ~750 µs here. |
| **Areal density `σ`** | Mass per unit area. The tamper's figure of merit, conserved through vaporisation. |
