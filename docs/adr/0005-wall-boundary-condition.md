# Wall BC for the 1D rad-hydro: cold black absorber + semi-infinite conducting solid

The wall boundary condition fixes two of the four loss channels (radiative-to-plate and
conductive-to-plate), so it is pinned explicitly rather than left to the solver default.

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

## Considered Options

- **Isothermal Dirichlet wall.** Rejected: imposes the interface temperature and makes the conductive
  flux a near-wall mesh artifact rather than a physical number.
- **Radiative wall with albedo / re-emission.** Rejected as unnecessary during the pulse: the cold
  surface makes re-emission ≤10⁻³ of the incident flux.
