# An effective-Isp deliverable pulls ablator mass, vehicle scale, and cadence inside the scope boundary that `f(v)` kept them outside

Design §11 places total-performance analysis, vehicle mass, pulse cadence, and multi-pulse plate
thermal accumulation out of scope, and design §7 explicitly reports ablation-per-pulse "as the
requirement back-propagated to the MEMS replenishment system rather than imposing it as a budget."
Those exclusions are correct **because `f` is a per-collision momentum ratio** — a property of one
collision, reusable regardless of how pulses are smoothed into vehicle motion, and blind to what the
vehicle carries.

The tamped-nozzle study (`puffsat_tamper_isp_prd.md`) does not have that property. Its deliverable is
**effective specific impulse**, which is a rocket-equation quantity: `Isp = J/(g₀·m_charged)`. Every
one of the parked items either sits in that denominator or constrains it. Reporting an Isp while
excluding the ablator would be reporting a fraction with a term missing from the bottom.

**What comes in, and why each is load-bearing rather than incidental:**

- **Ablator mass.** It is carried and expended per pulse, so it is charged. Its magnitude is
  currently uncertain by **27×** — between ~2% and ~60% of the propellant budget — which is the
  difference between a footnote and the largest term after the propellant itself.
- **Vehicle mass (1000 t) and cadence (1–4 Hz).** Needed to size the plate, to close the
  mass-flow cross-check, and to ask whether the plate can reject its absorbed heat between pulses.
- **Inter-pulse plate thermal accumulation.** Design §11 excludes it *on the grounds that it depends
  on a cadence the `f(v)` study never fixed*. With cadence pinned, the exclusion's own precondition
  is gone, so it is computed rather than deferred. (The answer turns out to be benign — an ablating
  surface is temperature-pinned, so the plate self-limits near 750–905 K — but that is a result, not
  an assumption, and it could not have been asserted without doing the work.)

**What stays out, deliberately:** the conversion rate between projectile consumption and delivered
payload (that needs program economics this repository has no basis to model — projectile economy is
reported as `β` instead, and no combined figure of merit is constructed); shock-absorber stroke;
whole-plate structural design beyond the closed-form checks ADR-0027 already permits.

**Precedent.** ADR-0027 established that a special scenario may host a bounded exception to the §11
boundary when its scale makes a parked question first-order, kept narrow so it does not erode the
boundary for the envelope study. This is the same move, with a different trigger: not scale, but a
change in the *kind* of quantity being delivered.

**Mitigation — a two-pass denominator.** The ablator's 27× uncertainty is a *plate* question, and
letting it gate the *tamper* question would be a scoping failure of its own. So Isp is evaluated in
two passes: **Pass 1 excludes ablator mass entirely** and is reported as an explicit upper bound
(it also reproduces the denominator prior work used, so it is directly comparable); **Pass 2** folds
in the measured ablator. Excluding it is *expected* to be conservative for the tamper verdict,
because `E/J = w/(2β)` means raising `β` lowers total energy delivered per unit impulse. That is an
expectation, not a proof: `E/J` fixes the energy, not the ablation it causes, and a tamper can also
change the fraction reaching the plate, its angular and radial distribution, arrival velocity,
residence time, and vapour-shield formation. PRD Rung 6 measures ablator mass per configuration and
tests the expectation rather than inheriting it.

## Consequences

- The repository now hosts a **plate thermal-balance calculation** and an **ablator mass budget** it
  previously, and correctly, excluded. Both are scenario-scoped to this study; the `f(v)` envelope
  study keeps the strict §7/§11 boundary.
- **Pass-1 Isp figures must never be quoted without their upper-bound label.** The excluded term is
  not small and not bounded.
- Ablator work is sequenced **last** (PRD Rung 6), after the tamper verdict exists, because it is
  expected to move the absolute number without changing the verdict — even the pessimistic branch
  takes a 984 s Pass-1 figure to ~610 s, still ~1.6× methalox. If Rung 6 finds the ablator term is
  configuration-dependent enough to move the comparison, the verdict is re-taken at Pass 2.

## Considered Options

- **Hold the §7/§11 boundary and report `f` only.** Rejected: `f` does not answer the question this
  study exists to ask, and a tamper's whole case is a mass trade that `f` cannot express.
- **Report Isp but treat the ablator as a plate requirement, not a budget** (the design §7
  precedent). Rejected: it hides a term potentially larger than everything except the propellant,
  inside the one quantity being delivered.
- **Import a projectile-to-payload exchange rate and report a single combined figure of merit.**
  Rejected: it would make a physics result move whenever an external economic model changed, and
  bury a program-level judgement inside a physics deliverable.
