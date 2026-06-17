# O2 and ANFO-mix spot-check tables: full pipeline re-run scoped to the optimum, compared at the water optimum

ADR-0007's table pipeline is explicitly water-only. The §4 confirmatory spot-checks — oxygen (LOX)
and an ANFO-product mix (N2/H2O/CO2/CO) — need their own EOS/opacity. This ADR fixes how they are
built and how they are compared, without water-scaling the answer.

**Full pipeline re-run per species, scoped to the optimum's `(ρ,T)` neighborhood only.** The
spot-checks are run *at the optimum, not swept* (§4), so the table is needed over a narrow `(ρ,T)`
band, not the full envelope — a handful of points, cheap even for line-by-line opacity.

- **O2: a trivial re-run.** Single species: CoolProp (real fluid + two-phase), CEA/Saha (O2→2O
  dissociation + ionization), HITEMP/ExoMol (O2 bands) + TOPS/OPLIB (atomic O/O+ plasma). Same
  pipeline, different input.
- **ANFO mix: lean on CEA's native multi-species equilibrium** for the EOS and the `(ρ,T)`-dependent
  equilibrium composition (CEA is built for mixtures — arguably easier than pure water). Opacity is
  the real work: the **equilibrium-composition-weighted sum of per-component spectral opacities**
  (HITEMP covers all of N2/H2O/CO2/CO; TOPS/OPLIB per element for the plasma), **summed spectrally
  first, then Rosseland/Planck-mean-ed** — the means do not add linearly across species. Composition
  shifts with `(ρ,T)`, so opacity couples to the CEA equilibrium at each point.

**Compared at the water optimum, not each species' own optimum.** The spot-checks run at the water
optimum's cloud conditions (same ρ, shape, footprint), an apples-to-apples test of "is water's
`f(v)` representative?" Re-optimize a species only if it fails there.

**Why water is expected to remain the conservative baseline** (the claim the spot-check defends):
heavier species (O2=32, CO2=44) run cooler for given KE → less ~T⁴ radiation → higher `e_eff`; H2O
condenses most readily (2.26 MJ/kg) → most condensation loss, the dominant low-v sink; O2 barely
condenses in the bounce (90 K). The most likely thing that could *break* the claim is the ANFO mix's
CO2 being a strong IR radiator, so the spot-check reports `e_eff` + the full loss decomposition
(ADR-0016) and watches that channel.

## Considered Options

- **Water-equivalent scaling** (treat O2/ANFO as water via mean-MW / effective-γ). Rejected: the
  spot-check exists to confirm water is representative and conservative; approximating them as water
  assumes the answer.
- **Full swept tables per species.** Rejected: the spot-checks are confirmatory at the optimum, not
  swept (§4); a narrow `(ρ,T)` band suffices.
- **Spot-check at each species' own re-optimized optimum.** Rejected as the baseline: it confounds
  the species effect with re-optimization; the water optimum isolates the comparison.
