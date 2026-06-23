# The 2D `eta_capture` kernel is an axisymmetric finite-volume HLLC Godunov solver (MUSCL-Hancock, Strang-split, cylindrical source)

The geometry track (ADR-0003) needs a second hot-path kernel: a **radiation-free 2D axisymmetric
Euler solver** for `eta_capture`. This ADR pins its numerics — the 2D sibling of the 1D kernel's
ADR-0022 (which chose a Lagrangian staggered-AV scheme). The 2D kernel is a **from-scratch Eulerian
finite-volume Godunov** solver (`crates/euler2d`), deliberately a *different* scheme from the 1D
kernel so that the two agreeing on the bounce is an independent cross-check, not a tautology.

**Equations.** 2D axisymmetric (cylindrical `r`–`z`, `∂/∂θ = 0`) compressible Euler, conserved
vector `U = [ρ, ρu_z, ρu_r, E]`. Axisymmetry is the **conservative cylindrical radial update**
`(1/r)∂(r F_r)/∂r` (radius-weighted faces) plus the pressure source `p/r` in the radial-momentum
equation — so mass, axial momentum, and energy conserve exactly and only radial momentum carries a
geometric source. The axis `r = 0` needs no special case: the inner face area `r_{1/2} = 0` kills
the flux through it.

**Scheme.** Second-order **MUSCL-Hancock** (van-Leer-limited primitive-conserved slopes, a Hancock
half-step predictor, a positivity fallback to first order where a face would lose positive density
or pressure) with the **HLLC** interface flux (Davis wave speeds; HLLC, not HLL, to resolve the
cloud/vacuum contact). Time advance is **Strang dimensional splitting** `Z(dt/2)·R(dt)·Z(dt/2)`. A
density/pressure **vacuum floor** resets evacuated cells (the rarefaction tail) to a quiescent state
— the 2D analogue of the 1D tail guard (ADR-0001).

**EOS.** A calibrated effective-γ ideal gas (ADR-0008), baked into the kernel — no equilibrium-table
lookup, since `eta_capture` is a ratio in which EOS error is largely common-mode.

**Verification (the rung's exit criteria — analytic + convergence, CLAUDE.md).** Sod shock tube
(embedded in 2D, against `hydro1d`'s exact Riemann oracle); the cylindrical **Sedov** blast (self-
similar `R_s ∝ t^{1/2}`) and **Noh** implosion (exact post-shock density 16 for γ = 5/3) — both
exercising the axisymmetric source; a smooth-advection **order-of-accuracy** test (≈ 2nd order in
L1); and the confined plane-wave bounce reproducing the independent 1D kernel's `1+e_eff`.

## Considered Options

- **Reuse the 1D Lagrangian staggered-AV scheme (ADR-0022) in 2D.** Rejected: Lagrangian remap in 2D
  is far more complex than Eulerian finite volume for the capture geometry, and using the *same*
  scheme would forfeit the independent-kernel cross-check the Eulerian solver provides.
- **First-order Godunov only.** Rejected: too diffusive for the rebound geometry; the order-of-
  accuracy gate requires a second-order scheme. (First-order was the D0 stepping stone.)
- **Unsplit (corner-transport) integration.** Deferred: Strang splitting with MUSCL-Hancock is
  second-order, simpler, and adequate for the verification suite; CTU can be revisited if the
  splitting error ever shows up against the cross-code (rung F).
- **Body-fitted / staircase concave boundary.** Out of scope here (flat plate only); when the
  shallow-concave plate lands, a **cut-cell immersed** reflecting boundary is preferred — a
  staircase would bias the axial-momentum ratio `eta_capture` measures.
