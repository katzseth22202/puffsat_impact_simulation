# The transitional anchor is the least-certain point: compute τ, deploy transport there, expect a possible e_eff minimum

The velocity sweep is anchored at the two well-modeled ends (3.2 and 16 km/s) plus one or more transitional points whose velocity is **located by measurement, not assumed**: a dense sweep finds where `τ` crosses ~1 and where `e_eff` dips. The transitional anchor
(partial ionization spanning ~6–10 km/s) is treated specially because it sits at the
confluence of the two model weaknesses and may be the physical worst case for radiative loss.

**Modeled with the high-v package.** It needs the equilibrium ionization EOS and plasma opacity that
the low-v package lacks, so the transitional anchor is the high-v package run at low `v`, not a third
code.

**Flagged as the least-certain anchor, for two stacked reasons:**
- `τ ~ 1` is where flux-limited diffusion is least accurate — Levermore–Pomraning is exact only at
  `τ≪1` and `τ≫1`; the transition is the interpolation gap.
- Partial ionization (10–20 kK) is where opacity is hardest (the §4 "not Kramers" regime) and where
  the molecular↔plasma table seam (ADR-0007) falls.

**`τ` is computed from the real tables at the anchor, not assumed.** Partial-ionization line opacity
may put it at `τ≫1` (FLD fine) or `τ~1` (FLD weakest); that measurement decides whether there is a
problem at all.

**If `τ~1`, the transport-level check is deployed here, not at the high-v end.** FLD is already exact
at `τ≫1`, so an independent radiation *model* (Quokka M1, or an Sn/short-characteristics solve) earns
its keep only near the transition. At minimum, FLD is bracketed between an optically-thin-emission
bound and an optically-thick-diffusion bound. The transitional `f` is reported with radiation-model
error bars.

**`e_eff` may have a local minimum at the transition.** Radiative loss can peak here: too cold to
radiate at 3.2 km/s, radiating-but-*trapped* at 16 km/s (`τ≫1`), but at the transition hot enough to
radiate strongly *and* thin enough (`τ~1`) for that radiation to escape to the wall and sideways. So
velocity is swept densely (~5–9 km/s) around the transition rather than interpolated, to catch the dip. The number of transitional anchors is itself an outcome of this sweep: one by default, but if the data reveal two distinct features at different velocities — a `τ~1` radiative-leak dip and a separate dissociation/ionization (EOS specific-heat) feature — both are adopted as anchors.

## Considered Options

- **Trust FLD at the transition like the other anchors.** Rejected: `τ~1` is exactly where FLD is
  weakest — the one point that most needs a transport check.
- **Interpolate `f(v)` between the well-modeled ends.** Rejected: `e_eff` may dip at the transition,
  which interpolation would miss.
