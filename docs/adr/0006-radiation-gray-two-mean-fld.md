# Radiation: gray two-mean flux-limited diffusion (Rosseland + Planck), multigroup gated

The 1D rad-hydro uses gray flux-limited diffusion with **two** opacity means, spanning the full
velocity envelope, with multigroup deferred behind a gate.

**Two means, each in its correct regime.** Frequency-dependent opacity is reduced to two gray means
that weight it differently: **Rosseland** (harmonic, `∂B/∂T`-weighted) in the **diffusion
coefficient**, correct optically thick (high-v, `τ≫1`); **Planck** (arithmetic, `B`-weighted) in the
**emission/absorption source**, correct optically thin (low-v, `τ≪1`). Both are tabulated. Using
Rosseland alone — the common shortcut — gets optically-thin emission wrong and therefore mis-states
the **low-v optical fraction**, which is exactly the quantity that gates the dark-oil seed study (§7).
So even though the low-v optical loss is small in magnitude, the second mean is required to keep that
go/no-go decision honest.

**Levermore–Pomraning flux limiter** ties the regimes together: the single FLD solver reduces to
`c·E` free-streaming at `τ≪1` and Fickian diffusion at `τ≫1`, so one model spans 16 km/s (thick) to
3.2 km/s (thin) without switching. Without it, gray diffusion over-transports radiation in the thin
low-v gas (unphysical superluminal flux).

**Gray baseline, multigroup gated.** Multigroup is built only if the emergent-flux *spectrum* proves
to matter — band-selective ablator absorption, or a gray survivability flux that is sensitive to the
mean. This is also where the Castro MGFLD cross-check (rung F) would apply.

## Numerics (settled at B3c)

The discretization choices that fall out of the Su–Olson Marshak-wave acceptance (the Rung-B exit
criterion). All three were validated against the published `ε=0.1, τ=1` benchmark table.

- **Energy-based matter coupling.** The linearized-implicit substep returns the matter
  internal-energy change `Δe` (energy/volume), not a temperature change `δT`; the caller advances
  `e_mat += Δe` and inverts the EOS for `T`. A `T += Δe/c_v` update divides finite absorbed energy
  by the vanishing `c_v = αT³` at the cold wave front and overflows. This is a real robustness fix
  for steep heat capacities, not a benchmark hack — it matters for the production water table near
  phase change / ionization, where `c_v` swings hard.
- **Marshak (incident-current Robin) boundary.** A radiation source surface enters as
  `F(0) = (c/2)(e_inc − E₀)` — a surface conductance `(c/2)/dx` independent of the interior
  diffusion `D`. Pure Dirichlet wrongly pins `E(0) → e_inc` (the benchmark surface value is `~0.55`,
  not `1`). The same incident-flux BC feeds the cold-black-absorber wall in ADR-0005 (B4).
- **Limiter is verification-switchable.** The production default is Levermore–Pomraning, but the
  solver also exposes a **Fick (`λ≡1/3`) mode**. Su–Olson is a *pure-diffusion* benchmark, so the
  tight table match runs in Fick mode; LP deliberately departs from pure diffusion near the steep
  front and is verified separately by the free-streaming-cap test (`|F| ≤ cE`, → `cE` as the
  gradient steepens). This keeps the LP default honest without contaminating the diffusion oracle.

## Considered Options

- **Single-mean (Rosseland-only) gray diffusion.** Rejected: wrong in the optically-thin limit;
  corrupts the low-v optical fraction and the seed-study gate.
- **Multigroup from the start.** Deferred: gray two-mean spans the envelope; pay for groups only if
  the spectrum matters.
- **Temperature-based coupling (`δT` return).** Rejected at B3c: overflows as `c_v → 0` at a cold
  front; the energy-based form is conservative against the radiation update and stays finite.
