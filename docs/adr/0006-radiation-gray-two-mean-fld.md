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

## Considered Options

- **Single-mean (Rosseland-only) gray diffusion.** Rejected: wrong in the optically-thin limit;
  corrupts the low-v optical fraction and the seed-study gate.
- **Multigroup from the start.** Deferred: gray two-mean spans the envelope; pay for groups only if
  the spectrum matters.
