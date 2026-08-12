# The tamper is an isentropic piston, not a mirror: its recoil is credited, and the deliverable is a realization fraction against a single ceiling

The tamped head-on nozzle (`puffsat_tamper_isp_prd.md`) inherited its objective function from a
problem statement whose deliverable was a **tamper multiplier** `Λ(τ, k)` — tamped impulse over
bare impulse — with the pay/don't-pay criterion `Λ > 1 + τ`, since a tamper expelled without
contributing thrust divides the effective exhaust velocity by `(1+τ)`. That framing is replaced.
It is not merely awkward: it makes the tamper's own recoil look like the mechanism that kills the
concept, when the accounting says the opposite.

**The ceiling.** Work the whole event in the vehicle frame, `+z` = thrust direction. Impulse on the
vehicle is the momentum brought in minus the momentum carried out, `J = −m_i·w − P_ejecta`. The
ejecta's total energy is `½ m_i w²` (slug and tamper start at rest in the vehicle frame), so by
Cauchy–Schwarz `|P_ejecta| ≤ √(2 M_ej E)` with `M_ej = m_i(1+K)`. Writing `K = k(1+τ)` for the
**total carried mass per projectile kg**, per projectile kg:

```
J_max   = w·(√(1+K) − 1)
v_e,max = J_max/K = w·(√(1+K) − 1)/K
```

Equality requires every ejecta element to end up moving `−z` at one common speed, which requires the
plate to catch and reverse all plate-bound material. This reproduces the prior analytic results
exactly (`β_ideal = √(1+k) − 1`; bare-plate optimum at `k* = 7.057`; 1014 s at `w̄ = 77.28 km/s`), so
it is a reformulation of the existing model rather than a new one.

**Two consequences, both counterintuitive.**

1. **At the ceiling the tamper is neutral.** `v_e,max` depends only on `K`, not on how `K` splits
   between slug and tamper. **A tamper can therefore never be justified as a momentum multiplier** —
   only as a *realizability* device. The deliverable becomes a **realization fraction**: what share
   of the ceiling a configuration achieves.
2. **The tamper's backward recoil is credited, not lost.** Its `−z` momentum is the exact accounting
   mirror of the `+z` impulse the plate received. Recoil is not the failure mode the `Λ` framing made
   it look like.

**What is actually lost is entropy.** Energy thermalised in the tamper must be re-expanded a second
time, isotropically, at a second efficiency cost. Entropy production rises steeply with shock
strength, so gentle, early, pressure-mediated loading beats a single strong ram. The tamper should
therefore be instrumented on its **entropy/thermal budget and its centre-of-mass momentum**, not on a
reflected-speed ratio. This is what selects a filled interlayer over a vacuum standoff, and a
thicker gently-loaded shell over a thin dense mirror.

**The reformulation also sharpens the verdict into one number.** At `w = 75 km/s`, `k = 7.06`, with
`β` in units of `w` per projectile kg:

| configuration | `K` | `β` | % of ceiling |
|---|---|---|---|
| bare plate, `k` = 7.06 | 7.06 | 0.9087 | **49.4%** |
| bare plate, `k` = 14.12 — same mass spent as *slug* | 14.12 | 1.6835 | 58.3% |
| perfect-mirror tamper, τ = 1 | 14.12 | 1.7825 | **61.7%** |
| **break-even against the bare plate at `k` = 7.06** | 14.12 | 1.817 | **62.9%** |

A *perfect* mirror misses break-even by 1.2 percentage points. So the tamper's entire case rests on
pressure coupling, which no existing model of this device contains — every one is ballistic
(straight-line elements, `ρv²` only, no `P` term), and at `k ≈ 7` the fireball's sound speed and its
recoil velocity are equal to within 5%, the worst possible regime in which to omit pressure.

## Consequences

- **The study reports a realization fraction, not a multiplier.** A tamped configuration is scored
  against the same ceiling as a bare one, so "spend the mass as tamper" and "spend it as slug" are
  directly comparable — and the second is only 3.4 points behind while being strictly simpler
  (no interlayer, no porosity model, no assembly, no RT exposure). The tamper must beat a simpler
  rival, which the `Λ` framing hid.
- **A plausible outcome of the study is "use more slug."** That is a legitimate result and is not to
  be argued past.
- **Kernel instrumentation follows the piston framing**: entropy budget and tamper CM momentum are
  first-class outputs; reflected-speed ratio is not measured as a deliverable.

## Considered Options

- **Keep the `Λ > 1 + τ` multiplier criterion.** Rejected: it treats recoil as loss when the
  accounting credits it, it cannot compare a tamper against the same mass spent as slug, and it
  frames a realizability question as a momentum question.
- **Cold specular mirror as the design objective** (maximise reflected speed, minimise energy
  uptake). Rejected: self-contradictory at finite mass — taking momentum `p` costs exactly `p²/2m_t`
  — and the table above shows it cannot pay at `τ = 1` even when perfect.
- **Purely empirical `Λ(τ, k)` with no mechanism claim.** Rejected: it gives no design guidance on
  standoff, thickness, or interlayer before an expensive sweep, and provides no basis for choosing
  what to instrument.
