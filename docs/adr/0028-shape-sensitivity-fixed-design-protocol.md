# Pulse-shape sensitivity is measured as raw `f` at a fixed design, fixed mass, on one fixed grid

The shape-sensitivity study (design §13) asks whether slight pulse-shape changes produce only
slight changes in `f` and `J_wall` — the smoothness property that makes delivery dispersion a
linear guidance correction for the pushed vehicle. This ADR records the three protocol decisions
the quoted sensitivities `S_x = (Δf/f)/(Δx/x)` depend on. They are worth recording because each
inverts an existing repo convention or an obvious default, and because redoing the study under a
different protocol would invalidate the committed `shape_sensitivity.*` artifacts.

## Decision 1: raw `f(shape)` at a fixed design point — not the survivability-constrained frontier

The deliverable everywhere else in the repo is the *best survivable* `f`: an argmax over geometry
subject to the peak-pressure gate. That object is intrinsically discontinuous — a small shape
change can push `peak_wall_pressure` across the gate and knock a corner out of the feasible set,
jumping the constrained optimum. Measuring its sensitivity would manufacture a cliff that has
nothing to do with guidance, where the plate is already built and only the arriving pulse varies.
So the study freezes the design at the best-survivable baseline nominal (`d/D = 0.1`, `L/D = 0.3`,
`r_foot/R = 0.5`, M = 20; dip and 16 km/s anchors) and perturbs only the pulse. Survivability is
kept as a **separate margin check** (`peak_wall_pressure` under the gate over the whole shape box),
so the one genuine step a dispersion could trigger is reported as a margin number, not smuggled
into the smoothness claim.

## Decision 2: fixed pulse mass and speed — shape moves `f` through *both* ADR-0003 factors

Delivery dispersion distorts the cloud; it does not change how much gas PuffSat brought. Holding
`m` and `v` fixed pins `p_in` and makes "same pulse, different shape" well-defined — but it means
reshaping changes the areal density `Σ`, so `e_eff` moves too. The protocol therefore re-runs the
1D pipeline fresh at each sample's `Σ` (equilibrium chemistry, the headline convention) rather
than interpolating frontier data that was never gridded in `Σ`. For the taper axis, where `Σ`
becomes radius-dependent, `e_eff` is evaluated at the mass-weighted mean `Σ` with a `Σ`-min/max
bound on the profile effect — deliberately *not* un-deferring the Rung-D `Σ`-resolved `e_eff(ρ)`
work; if the bound comes back non-small, that deferred work has become load-bearing and the study
halts there as a finding. A three-point frozen-chemistry spot-check at the dip anchor keeps the
smoothness claim from silently inheriting an equilibrium-only slope (ADR-0026).

## Decision 3: one fixed grid for the whole box — inverting the geometry sweep's scaling convention

The existing D-cc geometry sweep sizes the domain to the cloud (`r_foot = 1` WLOG, plate scaled to
match). Under that convention the immersed plate's grid representation changes between samples,
and IBM stair-step jitter lands at exactly the `Δf ~ 0.005–0.02` differences the study resolves.
The sensitivity runs invert it: the plate and grid are frozen once, and the cloud varies on top.
Companion noise rules: the `f` assembly normalizes by the *measured* initialized `p_in` (never the
analytic value), the noise floor `σ_noise` comes from ≥ 3 refined-grid repeats, structure below
`2σ_noise` is reported as below-noise (a pass), and a flagged cliff must survive grid refinement
before it is called physical.

## Scope boundary: symmetry-breaking modes are bounded, not computed

Lateral offset, tilt, and sideways drift are outside the axisymmetric kernel and stay outside
(§11 keeps spatially-resolved off-center `f` out of scope; ADR-0002/0023 economics rule out a 3D
kernel for this question). They are covered by the §13 analytic linear bound — slope constants
plus the rim-clip margin (`δ = R − r_foot`, itself C¹ not a step) — which is consistent with §11,
not an exception to it. The **shape box (±20%, one-sided taper/divergence) is an assumption, not
derived dispersion**; the deferred cloud-schedule study owns the real delivery numbers, and every
quoted `S` carries that caveat.

## Considered Options

- **Sensitivity of the constrained frontier.** Rejected: measures optimizer discontinuities, not
  the guidance-relevant physics (Decision 1).
- **Interpolate existing sweep/frontier data instead of fresh runs.** Rejected: the existing grids
  were domain-rescaled per case (grid noise between points) and never gridded in `Σ`; the 1D/2D
  runs are cheap at ~13 samples × 2 anchors × 2 plate shapes.
- **Fixed pulse *density* instead of fixed mass.** Rejected: would decouple `e_eff` from shape
  entirely, understating sensitivity — dispersion conserves mass, not density.
- **Full frozen-chemistry box.** Rejected as the default: `S` is a relative derivative, largely
  level-insensitive; the dip spot-check is the cheap guard, and a wild frozen slope would upgrade
  the ADR-0026 caveat rather than re-run the box.
- **3D kernel for offset/tilt.** Rejected: a discipline step-change for modes an analytic linear
  bound already answers at the "slight → slight" level the mission argument needs.

## Consequences

- The committed `shape_sensitivity.csv` / `.png` are protocol-bound: comparing their `f` values
  against the domain-rescaled geometry-sweep grids point-for-point is invalid; only the *nominal*
  sample is expected to reconcile (within the noise floor) with the frontier value.
- The `Σ`-min/max taper bound is a tripwire that can re-prioritize the deferred `Σ`-resolved
  `e_eff(ρ)` work; the frozen spot-check can upgrade the ADR-0026 caveat. Both are named halt
  conditions, not silent approximations.
- New `SlugConfig` fields (edge taper, divergence `α`) must default to the current behavior with
  identity-regression tests, so every pre-existing result stays valid unchanged.
