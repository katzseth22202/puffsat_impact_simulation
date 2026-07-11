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

## Amendment (2026-06): the shallow-concave reflecting boundary is a ghost-cell immersed boundary (true-normal mirror), not a literal cut-cell

When the shallow-concave plate landed (Rung D follow-on), the dished surface `z_s(r) =
depth·(r/r_plate)²` cuts diagonally across the square `(z,r)` mesh, so the flat plate's "reflect at
the `z=0` grid line" trick no longer applies. The boundary method chosen is a **ghost-cell immersed
boundary (IBM) with a true-normal mirror**, not the cut-cell named below:

- Cells whose centers lie *under* the surface are masked **solid**; before each of the three Strang
  sub-sweeps, `apply_immersed_bc` refills each solid cell by mirroring the fluid across the **true
  local surface normal** `n̂` (image point `x_i = x_c + 2·dist·n̂`, sampled from the nearest fluid
  cell): copy ρ and p, reflect velocity `u → u − 2(u·n̂)n̂`. The square mesh and the verified D0–D3
  sweeps are untouched (the pass is gated on `plate_profile.is_some()`), and the timestep stays the
  regular CFL — there is **no cut-cell small-cell instability** to manage.
- **Why this honors the anti-staircase rationale below.** The staircase bias this ADR warned against
  comes from snapping the wall normal to the grid axes; the true-normal mirror uses the *exact* `n̂`,
  so a free-slip wall exerts a purely normal impulse. This is gated, not assumed: a **D4b normal-only
  wall-impulse test** confirms `|Δp·t̂|/|Δp·n̂| < 0.15` off an inclined wall (a staircase gives ≈ 0.5).
- **Why this is safe for the `eta_capture` ratio.** `eta_capture`'s flat denominator (the confined
  plane wave) is grid-aligned, so the IBM's boundary error does **not** cancel against it. Mitigation:
  the **flat plate is also run through the IBM** (raised to `z₀ > 0`), so flat and concave differ
  *only* by the surface profile; and a **consistency gate (D4c)** ties the IBM flat wall to the
  verified grid-aligned flat `eta_capture` to rel < 0.10 — load-bearing, not cosmetic, since it is
  what makes the curvature gain trustworthy.
- **Axial-force simplification:** `J_wall^2D = Σ_r P(r)·(n̂·ẑ)·dA = Σ_r P(r)·r·dr` exactly — the
  `√(1+s²)` arc-length factor cancels the `1/√(1+s²)` of `n̂·ẑ`, so the axial impulse is the surface
  pressure on the *projected* annulus regardless of slope (design §8).

The literal **cut-cell / body-fitted boundary stays the path only if a future rung needs higher
boundary fidelity** than the IBM delivers; nothing measured here motivated it.

## Correction (2026-07-10): rim side face added to the dish SDF; face states guarded by a signal-speed envelope

Chasing an unphysical `eta_capture = 7.6` in the Jupiter sweep's M = 40 anchor
(d/D = 0.10, L/D = 0.3, r_foot/R = 0.5) exposed **two kernel defects**, both fixed:

1. **The dish rim's vertical side face was missing from the IBM geometry.** `PlateProfile::Dish`
   modeled only the top surface `z_s(r)`, so solid cells just inboard of the rim were mirrored
   across the faraway *top* instead of the nearby side face — spurious radial fluxes at the rim
   corner acted as an energy source, bounded at M ≲ 20 but self-exciting at M = 40. Fix:
   `normal(z, r)` returns the side-face normal `(0, 1)` where the side is nearest, and
   `signed_distance` is the intersection SDF `max(top_distance, r − r_plate)`.
2. **Unguarded MUSCL-Hancock face states at vacuum-floor-adjacent faces.** The reconstructed face
   density can approach zero while the momentum slope stays finite, so `u = m/ρ` diverges even
   though the ρ/p positivity fallback passes; the garbage fluxes plus the vacuum-floor reset then
   create mass and energy non-conservatively. Fix: each Hancock face pair is validated against a
   **signal-speed envelope** (1.5× the stencil's `max(|u_n|, |u_t|) + 2c/(γ−1)` vacuum-front bound,
   Toro §4.6), falling back to first order at that face.

A third symptom — a phantom plateau force that never decayed, running the impulse integration to
the step cap — was downstream of these. Regression gates: a rim-side-face unit test (plate.rs) and
an M = 40 dished-bounce termination/restitution test (bounce.rs). **Consequence:** the pre-fix
M = 10 concave corner etas were ~1–1.6 % high and M = 20 systematically low; post-fix M = 10/20/40
agree to < 1 % (Mach insensitivity restored). Best survivable `f` moves dip 0.777 → **0.768**,
16 km/s 0.805 → **0.798** (pure `eta_capture` shift; the 1D `e_eff` anchors are untouched).

## Considered Options

- **Ghost-cell IBM with a true-normal mirror (chosen for the concave boundary, amendment above).**
  Keeps the verified square-grid scheme and CFL timestep; captures the true normal (no staircase
  bias); flat-vs-concave run apples-to-apples through the same boundary operator.
- **Reuse the 1D Lagrangian staggered-AV scheme (ADR-0022) in 2D.** Rejected: Lagrangian remap in 2D
  is far more complex than Eulerian finite volume for the capture geometry, and using the *same*
  scheme would forfeit the independent-kernel cross-check the Eulerian solver provides.
- **First-order Godunov only.** Rejected: too diffusive for the rebound geometry; the order-of-
  accuracy gate requires a second-order scheme. (First-order was the D0 stepping stone.)
- **Unsplit (corner-transport) integration.** Deferred: Strang splitting with MUSCL-Hancock is
  second-order, simpler, and adequate for the verification suite; CTU can be revisited if the
  splitting error ever shows up against the cross-code (rung F).
- **Staircase concave boundary.** Rejected: a grid-axis-snapped wall biases the axial-momentum ratio
  `eta_capture` measures (the D4b test catches it). Superseded by the ghost-cell IBM (amendment above).
- **Cut-cell / body-fitted concave boundary.** Deferred: highest boundary fidelity, but reshaping
  cells brings the small-cell timestep hazard, and the IBM's true-normal mirror already clears the
  anti-staircase gate. Revisit only if a later rung needs sub-IBM boundary fidelity.
