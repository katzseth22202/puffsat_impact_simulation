# The RT/RM deferral does not transfer to the tamped nozzle: there, mixing decides the result

ADR-0020 declines to model Rayleigh–Taylor and Richtmyer–Meshkov instability, on an argument that is
correct and worth restating: `e_eff` is a pure momentum ratio, so RT "corrugates the interface (more
area, faster mixing) but neither creates nor destroys axial momentum," and can therefore reach `f`
only through a loss channel, gated on a conductive-channel watchdog.

That argument is sound **for the interface it is about** — the near-wall boundary layer, where the
unstable surface separates two parcels of gas that are both going to deliver their momentum to the
plate regardless of how corrugated the boundary between them is. It does not transfer to the tamped
nozzle (`puffsat_tamper_isp_prd.md`), and the reason is not that the argument is weaker there but
that the *quantity being disturbed is different*.

**At the plume/tamper interface, mixing changes the mass distribution, not the interface shape.** The
tamper's entire function is inertial: it is a coherent, high-areal-density body that takes backward
momentum from the plume and is credited for it (ADR-0030). If RT entrains it into the plume, it stops
being a piston and becomes ordinary reaction mass — which the ceiling framework shows is a *worse*
use of the same kilograms. This is not a corrugation with conserved momentum on both sides; it is a
change in which body is carrying what.

**And the effect is comparable with the reference margin being decided.** Using the self-similar
mix growth
`h_b ≈ α_mix·A_RT·a_RT·t²` with `α_mix = 0.02–0.05` and Atwood number ≈ 1:

| configuration | `a_RT = P/σ` | window | bubble depth `h_b` | tamper thickness `h_t` | **bubble-depth fraction `h_b/h_t`** |
|---|---|---|---|---|---|
| snow slug, contact tamper (`s = 0.683 m`); slug-disassembly window | 1.9×10⁸ m/s² | 27.6 µs | 2.9–7.3 mm | 34.8 mm | **8–21%** |
| vacuum standoff, ice slug, tamper at `s = 1 m`; tamper-recoil window | 3.9×10⁸ m/s² | 45 µs | 16–40 mm | 16.2 mm | **100–246% — fully disrupted** |

**Neither row is the filled configuration as the PRD defines it** (ice slug → filled interlayer
→ tamper at standoff), and the two rows use different windows. The first is the snow-*slug*
corner of the original `(slug density, standoff)` scoping, carried as the filled design's
stand-in; since `h_b ∝ a_RT·t²` the substitution is not neutral and cuts both ways. The filled
configuration's own numbers are an output of the PRD's Rung 1 measured `a_RT(t)`, not values
this table supplies.

The screen predicts that the standoff configuration's tamper is **shredded before it finishes
its job**, across the whole plausible `α_mix` range — a thin sheet rammed by a fast plume is torn
apart, while a thicker shell pressed gently holds together. That is an independent mechanism
reaching the same verdict as ADR-0030's entropy argument, and together they provisionally foreclose
it. It remains a *constant-acceleration self-similar estimate with an assumed broad initial
perturbation spectrum*: a strong provisional foreclosure, not a verified result, and the control
run is what confirms it. On the stand-in geometry the bound is *not* decisive: a separate 16–63% total mix-width
heuristic maps the realization fraction
roughly between "piston" (~62% of ceiling) and "just extra slug" (~58%), near the 62.9%
reference-case threshold. Bubble-depth fraction and total mix-width fraction are not interchangeable,
and the map from either to realization fraction is a screening assumption. The uncertainty is
wider than the decision.

**It is load-bearing in three distinct places**, which is why it cannot be carried as a footnote:

1. **Tamper coherence** — piston or entrained payload, i.e. the result itself.
2. **The interlayer interface**, which wants the *opposite* outcome: mixing there is desirable,
   because it makes the interlayer a continuous pressure-bearing medium. Two interfaces, opposite
   requirements, same instability.
3. **Tamper shape survival** — at these accelerations 10 cm features grow ×20 over the confinement
   window while 1 m features grow ×2.6, so on a 0.68 m-radius tamper only gross shape survives. This
   is why the tamper's design variable is angular coverage rather than curvature.

**Axisymmetry is useful but incomplete.** An axisymmetric code can represent some modal growth and
provide bounds, but cannot reproduce a general three-dimensional perturbation spectrum or turbulent
mixing cascade. Escalation beyond an axisymmetric screen is a dimensionality question, not merely a
resolution question.

## Consequences

- The tamped-nozzle study carries an **explicit RT treatment as a numbered rung**, sequenced
  *second* — immediately after the run that produces the acceleration history it needs, and before
  any 2-D work — because no downstream effort can produce a verdict while this term is wider than the
  margin. It is escalated only as far as needed: a tightened analytic bound against the actual
  decaying `a_RT(t)` and a stated initial perturbation spectrum; then a mix model; then a resolved
  spot-check.
- The RT contribution is reported as an **explicit ± on the realization fraction**, at whatever level
  it is settled.
- **ADR-0020 stands for the `f(v)` study.** Its watchdog and its scope are unchanged; what changes is
  that its reasoning is now recorded as interface-specific rather than general.

## Considered Options

- **Inherit ADR-0020 and carry RT as a stated caveat.** Rejected: the caveat would be wider than the
  quantity being decided, which makes the deliverable unfalsifiable rather than merely uncertain.
- **Go straight to a 3-D resolved simulation.** Rejected as the *first* step: expensive, and the
  present bound is deliberately crude in two fixable respects (it assumes constant acceleration, and
  a broad initial perturbation spectrum), so the cheap refinement may settle it.
- **Design the tamper to be RT-immune.** Not available: the instability is driven by a light fluid
  accelerating a dense one, which is the tamper's operating principle, not an incidental feature of
  it.
