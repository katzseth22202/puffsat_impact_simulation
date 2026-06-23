# Wall BC for the 1D rad-hydro: cold black absorber + semi-infinite conducting solid

The wall boundary condition fixes two of the loss channels (radiative-to-plate, channel 1a, and
conductive-to-plate, channel 2; see ADR-0016 for the five-channel decomposition), so it is pinned
explicitly rather than left to the solver default.

**Radiation — cold black absorber.** During the ~µs pulse the surface stays cold (~300–600 K,
per ADR-0004) while the gas radiates at kK temperatures, so the wall's re-emission `σ·T_wall⁴` is
~10³–10⁶× smaller than the incident flux. The wall is modeled with absorptivity ≈ 1 and negligible
re-emission: the radiative loss channel is simply the flux-limited-diffusion flux that reaches the
wall. No wall albedo or radiative back-reaction during the pulse.

**Conduction — semi-infinite conducting solid, not isothermal Dirichlet.** The gas couples to a 1D
heat-conduction solution in the solid, so the conductive-to-plate loss is computed self-consistently
from the solid's effusivity `√(kρc)`, and the interface temperature *emerges* — the cold value of
ADR-0004 becomes a verifiable output, not an assumption. Consequence (accepted, arguably a feature):
the conductive loss becomes facesheet-material-dependent, which later enables "does a lower-effusivity
ablator cut conductive loss?" studies.

**Validity guard — single layer vs layered stack.** Baseline = a single semi-infinite solid with the
*ablator's* properties. But the per-pulse thermal penetration depth `δ = √(ατ)` is ~1–11 µm
(ablator `α ~ 1.3e-7 m²/s`, `τ ~ 1 µs–1 ms`), comparable to a renewed ablator layer thickness. If `δ`
reaches the SiC beneath, the model resolves the two-layer ablator-on-SiC stack — which gives a colder
interface (SiC effusivity ~30× the ablator's) and hence *more* condensation/conduction loss. So the
single-layer ablator baseline is the *less* conservative choice and must be checked against the actual
ablator thickness, not assumed valid.

**Turbulent (RT) enhancement — out of scope, but watched here.** This conductive solve is *laminar*.
Rayleigh–Taylor mixing at the near-wall boundary layer would enhance wall heat transfer above this
value, one-sidedly, and could push the true conductive loss above the floor's laminar number. RT is
not modeled (ADR-0020), but the conductive-to-plate channel (1D output) is its watchdog: if it grows
to a material slice of `(1−f)` at any anchor, a bounding turbulent-conduction correction is applied
there before `f` is quoted.

**Amendment (Rung B, 2026-06): gas-side thermal resistance — the conductive coupling's missing
half.** Coupling this semi-infinite solid to the *inviscid* gas (Euler + artificial viscosity,
ADR-0022) at the 16 km/s anchor over-drains the near-wall gas. The kernel has no physical gas-side
conduction, so no thermal boundary layer forms in the gas and the interface flux is limited *only* by
the solid's effusivity. The (physically-correct) first-step semi-infinite flux integral
`effusivity·ΔT·2√(Δt/π)` then exceeds the thin wall *gas* cell's heat content, zeroing it and
collapsing the bounce — a real physics/numerics gap, not a tuning issue. A faithful conductive
channel therefore needs a **gas-side conduction operator** (a parabolic diffusion term on the gas
energy + flux-continuous interface coupling), with its own coupled-conduction analytic acceptance
test. **Decision (2026-06): defer the conductive channel for the high-v `e_eff` pass** — the sweep
runs `wall = None` (`loss_conductive = 0`). Justified: `e_eff` is loss-insensitive (0.63 with no
conduction vs 0.64 lossless; ≤1.6 % over a 100× opacity swing — ADR-0007 amendment), so the
deliverable is unaffected. The fix is **bundled with the wall-flux/survivability rung** (design §10
"B-flux") alongside the real per-regime opacity table, since channels 1a + 2 are the same plate heat
load (ADR-0010/0011).

**Amendment (B-flux, 2026-06): gas-side conduction operator landed.** The gas now carries its own
thermal conductivity `k_gas` (a table field, CoolProp transport at low-v; ADR-0007), and conduction
is solved as a **combined gas+solid backward-Euler system** — one tridiagonal over the union
`[gas_{n−1}…gas_0 | solid_0…solid_{m−1}]` on the same Thomas solver (`Solid::step_coupled`). Every
face uses the series resistance of its two half-cells `G = 1/(½w_L/k_L + ½w_R/k_R)`, which carries the
non-uniform Lagrangian gas mesh, variable `k_gas`, and the gas|solid material jump in one formula
(reducing to `k/dx` for the uniform solid, matching the original step-temperature solve). The
**interface temperature now emerges from flux continuity** instead of being pinned to the bulk gas
temperature, so the over-drain cannot recur: the gas-side resistance bounds the interface flux.
Verified test-first against the **two-semi-infinite-media contact** analytic (interface jumps to the
effusivity-weighted `T_i = (e_g·T_g + e_s·T_s)/(e_g+e_s)`, erf profile each side), plus an
order-of-accuracy refinement test and an energy-closure check. Wired into both `CoupledBounce` and
`CondensingBounce`; it engages only where the table provides `k_gas`. **Still deferred:** the high-v
**plasma transport** `k_gas` (Spitzer-like) — the high-v plasma table carries none, so the 16 km/s
`e_eff` pass keeps `wall = None`; that, the real opacity table, and the survivability report remain
the B-flux high-v sibling. (Finding: at low-v, conduction with the real `k_gas` is negligible over the
µs bounce — see ADR-0004 amendment and design §10 Rung C.)

## Considered Options

- **Isothermal Dirichlet wall.** Rejected: imposes the interface temperature and makes the conductive
  flux a near-wall mesh artifact rather than a physical number.
- **Radiative wall with albedo / re-emission.** Rejected as unnecessary during the pulse: the cold
  surface makes re-emission ≤10⁻³ of the incident flux.
