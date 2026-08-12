# The tamper is an isentropic piston, not a mirror: its recoil is credited and configurations share one ceiling definition

The tamped head-on nozzle (`puffsat_tamper_isp_prd.md`) inherited its objective function from a
problem statement whose comparison metric was a **tamper multiplier** `Λ(τ_t, k)` — tamped
impulse over bare impulse — with the pay/don't-pay criterion `Λ > 1 + τ_t`, since a tamper
expelled without contributing thrust divides the effective exhaust velocity by `(1+τ_t)`.
That framing is replaced.
It is not merely awkward: it makes the tamper's own recoil look like the mechanism that kills the
concept, when the accounting says the opposite.

**The ceiling.** Work the whole event in the vehicle frame, `+z` = thrust direction. Impulse on the
vehicle is the momentum brought in minus the momentum carried out, `J = −m_i·w − P_ejecta`. The
ejecta's total energy is at most `½ m_i w²` (carried material starts at rest in the vehicle
frame), so by Cauchy–Schwarz `|P_ejecta| ≤ √(2 M_ej E)`. Define the mass ratios explicitly:

```
k     = m_slug/m_i
τ_t   = m_tamper/m_slug
μ     = m_interlayer/m_slug
K     = m_hydro/m_i                # near-field carried mass / projectile mass
      = k(1+τ_t+μ)                 # if these are the only material regions
K_ej  = (M_ej−m_i)/m_i             # nonprojectile ejecta mass / projectile mass
```

Then:

```
j_max   = J_max/m_i = w·(√(1+K_ej) − 1)
v_e,max = j_max/K_ej
```

In a closed no-ablator calculation, `K_ej = K`. Pass 2 must add ablator mass actually ejected
across the system boundary to `K_ej`; the charged-mass ratio `C = m_charged/m_i` separately
includes all expended ablator mass and is the effective-Isp denominator.

Equality requires every ejecta element to end up moving `−z` at one common speed, which requires the
plate to catch and reverse all plate-bound material. At `K_ej = k` this recovers prior work's
ideal-collimation coefficient exactly, `β_ideal = √(1+k) − 1`, so it is a reformulation of the
existing model rather than a new one. The other two figures usually quoted beside it — the bare-plate
optimum at `k* = 7.060` and 1014 s at `w̄ = 77.28 km/s` — are results of the **ballistic** model
(`β_bare`); the ceiling has no interior optimum in `K_ej`, since `v_e,max = w/(√(1+K_ej)+1)` falls
monotonically.

**Two consequences, both counterintuitive.**

1. **At the ceiling the tamper is neutral.** `v_e,max` depends only on `K_ej`, not on how the
   ejecta mass splits among slug, interlayer, tamper, and ablator. **A tamper can therefore never
   be justified as a momentum multiplier** — only as a *realizability* device. The hydrodynamic
   comparison metric becomes a **realization fraction** `r_real = J/J_max`: what share of the
   same-ejecta-mass ceiling a configuration achieves. Effective Isp remains the final deliverable.
2. **The tamper's backward recoil is credited, not lost.** Its `−z` momentum contributes usefully
   to the full-system momentum audit. It is not generally equal to plate impulse on its own,
   because the projectile, slug, interlayer, all escaping material, and vehicle share the balance.

**Entropy is one important source of ceiling shortfall.** Energy thermalised in the tamper must
be re-expanded a second time, isotropically, at a second efficiency cost. Angular dispersion,
ejecta velocity variance, residual internal energy, radiation, incomplete capture, and wrong-way
material can also lower `r_real`. Entropy production rises steeply with shock strength, so gentle,
early, pressure-mediated loading is the design hypothesis to test against a single strong ram.
Instrument the **complete shortfall ledger**, including the tamper's entropy/thermal budget and
centre-of-mass momentum, rather than only a reflected-speed ratio.

**The reformulation also sharpens the reference comparison into one number.** At `w = 75 km/s`,
`k = 7.06`, with `β = J/(m_i w)`:

| configuration | `K` | `β` | % of ceiling |
|---|---|---|---|
| bare plate, `k` = 7.06 | 7.06 | 0.9087 | **49.4%** |
| bare plate, `k` = 14.12 — same mass spent as *slug* | 14.12 | 1.6835 | 58.3% |
| perfect-mirror tamper, `τ_t = 1` | 14.12 | 1.7825 | **61.7%** |
| **break-even against the bare plate at `k` = 7.06** | 14.12 | 1.817 | **62.9%** |

A *perfect* mirror misses this **reference-case** break-even by 1.2 percentage points. The 62.9%
value applies only to `k = 7.06`, `τ_t = 1`, `μ = a_abl = 0` against the bare `k = 7.06`
control. Every swept candidate instead uses `β_candidate/C_candidate > β_reference/C_reference`;
its required `r_real` depends on its own `K_ej` and `C`. So the reference tamper's case rests on
pressure coupling, which no existing model of this device contains — every one is ballistic
(straight-line elements, `ρv²` only, no `P` term), and at `k ≈ 7` the fireball's sound speed and its
recoil velocity are equal to within 5%, the worst possible regime in which to omit pressure.

## Consequences

- **The study reports a realization fraction, not a multiplier, as its hydrodynamic comparison
  metric.** Every configuration is scored with the same ceiling definition at its own `K_ej`, so
  "spend the mass as tamper" and "spend it as slug" are directly comparable — and the second
  is only 3.4 points behind while being strictly simpler
  (no interlayer, no porosity model, no assembly, no RT exposure). The tamper must beat a simpler
  rival, which the `Λ` framing hid.
- **A plausible outcome of the study is "use more slug."** That is a legitimate result and is not to
  be argued past.
- **Kernel instrumentation follows the piston framing**: the complete shortfall ledger, including
  entropy budget and tamper CM momentum, is first-class; reflected-speed ratio is not a deliverable.

## Considered Options

- **Keep the `Λ > 1 + τ_t` multiplier criterion.** Rejected: it treats recoil as loss when the
  accounting credits it, it cannot compare a tamper against the same mass spent as slug, and it
  frames a realizability question as a momentum question.
- **Cold specular mirror as the design objective** (maximise reflected speed, minimise energy
  uptake). Rejected: self-contradictory at finite mass — taking momentum `p` costs exactly `p²/2m_t`
  — and the table above shows it cannot pay in the `τ_t = 1` reference case even when perfect.
- **Purely empirical `Λ(τ_t, k)` with no mechanism claim.** Rejected: it gives no design guidance on
  standoff, thickness, or interlayer before an expensive sweep, and provides no basis for choosing
  what to instrument.
