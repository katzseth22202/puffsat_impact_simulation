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

**And the effect is larger than the margin being decided.** Using the self-similar mix growth
`h ≈ α·A·a·t²` with `α = 0.02–0.05` and Atwood number ≈ 1:

| configuration | `a = P/σ` | window | bubble depth `h` | tamper thickness | **mixed fraction** |
|---|---|---|---|---|---|
| filled interlayer, contact tamper | 1.9×10⁸ m/s² | 27.6 µs | 2.9–7.3 mm | 34.8 mm | **8–21%** (16–63% total mix width) |
| vacuum standoff, tamper at 1 m | 3.9×10⁸ m/s² | 45 µs | 16–40 mm | 16.2 mm | **100–246% — fully disrupted** |

The standoff configuration's tamper is **shredded before it finishes its job**, across the whole
plausible `α` range — a thin sheet rammed by a fast plume is torn apart, while a thicker shell
pressed gently holds together. That is an independent mechanism reaching the same verdict as
ADR-0030's entropy argument, and together they foreclose it. For the filled configuration the bound
is *not* decisive: a 16–63% mixed fraction moves the realization fraction roughly between "piston"
(~62% of ceiling) and "just extra slug" (~58%), **against a 62.9% threshold.** The uncertainty is
wider than the decision.

**It is load-bearing in three distinct places**, which is why it cannot be carried as a footnote:

1. **Tamper coherence** — piston or entrained payload, i.e. the result itself.
2. **The interlayer interface**, which wants the *opposite* outcome: mixing there is desirable,
   because it makes the interlayer a continuous pressure-bearing medium. Two interfaces, opposite
   requirements, same instability.
3. **Tamper shape survival** — at these accelerations 10 cm features grow ×20 over the confinement
   window while 1 m features grow ×2.6, so on a 0.68 m-radius tamper only gross shape survives. This
   is why the tamper's design variable is angular coverage rather than curvature.

**An axisymmetric code cannot represent it.** Both kernels in this repository are axisymmetric by
design, so this is not a resolution question but a dimensionality one.

## Consequences

- The tamped-nozzle study carries an **explicit RT treatment as a numbered rung**, sequenced
  *second* — immediately after the run that produces the acceleration history it needs, and before
  any 2-D work — because no downstream effort can produce a verdict while this term is wider than the
  margin. It is escalated only as far as needed: a tightened analytic bound against the actual
  decaying `a(t)` and a stated initial perturbation spectrum; then a mix model; then a resolved
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
