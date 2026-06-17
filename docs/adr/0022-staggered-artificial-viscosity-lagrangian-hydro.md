# The 1D kernel uses staggered artificial-viscosity Lagrangian hydro, not Godunov

The 1D rad-hydro kernel (`crates/hydro1d`, rung A onward) discretizes the Lagrangian
equations on a **staggered mesh with artificial viscosity** (von Neumann–Richtmyer /
Wilkins): velocity at nodes, density/energy/pressure at cell centers, shocks captured by a
scalar artificial-viscosity term `q` rather than by a Riemann solver at each interface. Rung
A is ideal-gas (effective-gamma EOS), but the scheme is chosen for where it has to go —
rungs B/C add a tabulated equilibrium EOS (ADR-0007), gray flux-limited radiation diffusion
(ADR-0006), and a conducting/ablating wall (ADR-0005, ADR-0014).

**Extensibility is the deciding criterion**, because rung A is scaffolding for the real
kernel, not the deliverable. An artificial-viscosity scheme swallows an arbitrary tabulated
EOS with a bare `p(ρ,e)` call, takes the radiation step as a clean operator-split implicit
diffusion solve on the same mesh, and reads the **wall impulse** (ADR-0001's restitution
measure) directly off a boundary node. Two secondary pulls agree: the validation keystone is
reproducing **Project Orion** (ADR-0013's V&V), whose rad-hydro lineage *is* staggered-AV
Lagrangian, so this is the most faithful path; and the order-of-accuracy test runs on a
*smooth* problem regardless of scheme (any shock-capturing method is first-order at the
shock), so Godunov's sharper shocks buy rung A little.

**Consequence:** the exact Riemann solver written as the Sod oracle stays **test-only** — it
does not double as the kernel's flux function. That is accepted, and is in fact preferred: it
keeps the oracle a small, standalone, independently-testable closed-form helper. The cost of
artificial viscosity — one or two `q` coefficients to tune and shocks smeared over a few
cells — is a well-trodden, one-dimensional affair for ideal-gas Sod.

## Considered Options

- **Godunov-type Lagrangian** (Riemann problem solved at each interface) — sharper shocks,
  no `q`-tuning, and the exact Riemann solver could double as the interface flux. Rejected:
  the reuse synergy evaporates at rung B, where a *general* tabulated EOS forces the exact
  solver off onto an approximate/iterative Riemann solver with EOS calls inside the star-state
  iteration — exactly when you'd want to cash it in — while radiation coupling and the
  conducting-wall BC are both more involved in a Godunov framework.
- **Eulerian (fixed-grid) hydro** — not considered seriously: the problem is a gas column
  bouncing off a moving interface with a contact-tracked wall, and restitution is a wall
  impulse (ADR-0001); a Lagrangian mesh tracks the material and the wall natively. ADR-0002's
  factorization already commits the restitution leg to a 1D Lagrangian treatment.
