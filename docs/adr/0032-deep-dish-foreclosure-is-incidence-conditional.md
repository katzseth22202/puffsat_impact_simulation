# The deep-dish foreclosure is conditional on plane-wave incidence and must be reopened for a finite-standoff source

ADR-0021 forecloses the deep dish and states it is "never reopened." That foreclosure carries an
unstated precondition — **the incidence geometry of the arriving gas** — and under the tamped-nozzle
study's geometry (`puffsat_tamper_isp_prd.md`) the precondition fails, so the conclusion must be
retested rather than inherited.
This ADR records the precondition so that ADR-0021 can stand unchanged for the study it was written
for, and so a reader meeting a deep dish in the tamper work does not read it as ADR-0021 being
ignored.

**ADR-0021's mechanism, restated precisely.** A deep dish is rejected because "its focal hot spot
concentrates the rebound in the strongest-radiating (`flux ~ v⁸`), `τ_opt≫1` gas and couples radiation
into the otherwise radiation-free 2D geometry track." That is a statement about what a dish does to a
**plane-wave cloud**: a paraboloid maps parallel rays onto its focus, so a plane-wave rebound is
focused, and the hot spot is real.

**Ray optics motivate reopening, not the answer.** A paraboloid focuses *parallel → point* and
collimates *point → parallel*. In the tamper study the gas arrives from a **compact fireball at
finite standoff**, not as a plane wave, so the old plane-wave focal-hot-spot argument no longer
forecloses the shape. But a finite-duration, spatially extended, pressure-steered plume is not an
optical point source; whether it collimates usefully is a simulation result.

**Point-source screen.** For an instantaneous ballistic expansion every element's
trajectory is `(V_cm + u·r̂)·t`, so all rays trace back to a **fixed** origin; the receding centre of
mass skews the angular distribution but does not move the apparent source. The blur is the finite
disassembly time, `u·t_dis ≈ 25 km/s × 25 µs ≈ 0.6 m`, against a focal length of 6–10 m — about
0.04 rad.

**And the focus-matched ray-optics shape lands inside the foreclosed band.** For a paraboloid
with its focus at the fireball, dish depth-to-diameter is `δ/D = R/(8d)`: **0.19** at
`R = 15 m, d = 10 m` and **0.31** at
`R = 15 m, d = 6 m`, against ADR-0021's shallow band of `δ/D ≤ 0.15`.

**This does not mean the deep dish wins there.** Two effects cut against it, and the study treats
flat and parabolic as genuinely two-sided:

- **The prize is bounded at ≤ 23%.** The specular upper bound, mass-weighted over solid angle, is
  `(1+cosθ)/(2cosθ)` for a focus-matched parabola against `2cosθ` for a flat plate: 1.09 / 1.19 /
  1.23 at `R/d` = 1 / 2 / 2.5.
- **Stagnation blunts shape by roughly an order of magnitude.** A hypersonic gas does not reflect; it
  shocks, stagnates into a subsonic plenum, and re-expands, so curvature acts through *confinement*
  (raising the rim so the layer must relieve axially), not reflection. ADR-0021's own measured
  concave lift — `eta_capture` 0.915 → 0.977 → 0.994 over `δ/D` 0 / 0.10 / 0.15, about **+9%** — is
  the direct evidence for this, and the readings above 1.0 that once suggested otherwise turned out
  to be a rim-corner boundary artifact.
- **A focus-matched paraboloid has ~32% more surface area** than the flat disk it replaces
  (`F = d`, `R = 2.5d`). Since ablator mass is charged in the tamper study's denominator (ADR-0031)
  and ablation is *sub-linear* in fluence, extra area costs more total ablator, not less — and the
  parabola also collects flux at more normal incidence at the rim, where a flat plate takes a cosine
  discount.

So the parabola buys ≤ 23% of impulse for ~32% more area, and **flat may still win on mass.**

## Consequences

- **ADR-0021 stands unamended for the `f(v)` envelope study**, whose cloud is delivered as a shaped
  pulse onto the plate — plane-wave incidence, where its reasoning is correct and its foreclosure
  holds.
- **Its foreclosure now carries a stated precondition.** Any future study should check its incidence
  geometry before inheriting it. The general form: *a dish focuses whatever is parallel and
  collimates whatever is at its focus; the hot-spot objection applies only to the former.*
- The tamped-nozzle study sweeps flat plus the paraboloid family out to `δ/D ≈ 0.35`, and reports the
  shape result as an outcome rather than assuming the parabola.

## Considered Options

- **Amend ADR-0021 in place.** Rejected: ADR-0021's amendments are all internal to the `f(v)` study,
  and importing a second study's geometry into it would muddy a record that is correct as written. A
  pointer is added there instead.
- **Inherit the foreclosure and cap depth at `δ/D ≤ 0.15`.** Rejected: off-optimum by construction
  under point-source incidence, and it would inherit a conclusion whose stated mechanism does not
  apply — the worst kind of borrowed decision.
- **Adopt the focus-matched parabola as the baseline.** Rejected as premature: the ≤ 23% specular
  prize, the ~9% measured stagnation reality, and the ~32% area penalty do not obviously net
  positive, so shape is swept rather than chosen.
