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

## Amendment (2026-06): the EOS dip is found; the radiative dip stays pending

The dense sweep landed (`crates/sweep --transitional`, `analysis.py --axis v`; design §10
"Transitional anchor"). It separates the two anticipated features by the **two pieces of physics that
are computable now vs not**, and reports both as one decomposition:

- **EOS-only curve (`run_bounce`, radiation+conduction off) — computable now, and it dips.** This is
  the dissociation/ionization specific-heat feature of the ADR's "two distinct features" clause,
  isolated from radiation. The ρ-mean `e_eff(v)` has a clear **interior minimum ≈ 0.57 near 11 km/s**,
  below both the 0.74 (3.2 km/s) and 0.64 (16 km/s) endpoints — so the worst case really is the
  transition, not an endpoint (vindicating the "located by measurement" stance and the C/B-flux
  finding that neither endpoint is the floor). The dip velocity (~11 km/s) sits at/above the a-priori
  ~5–9 km/s window: it bottoms where the equilibrium chemistry most actively absorbs the stagnation
  enthalpy.
- **Radiation-on curve (`CoupledBounce`, interim opacity, `wall = None`) — the comparison, not the
  answer at `τ~1`.** It sits only ≈ 0.004 below the EOS curve (the interim-opacity band). But the
  interim Kramers opacity (`κ_R ∝ T⁻³·⁵`, B5c-2) is **structurally wrong at the transition** — it
  yields `τ ≫ 1` where reality is `τ ~ 1` — so it does **not** resolve the separate `τ~1` radiative-
  leak dip this ADR predicted. That dip is gated on the **real per-regime opacity table** (the deferred
  B-flux sibling), where `τ` is computed from the real tables (as this ADR requires) and, if `τ~1`
  there, the transport-level check (Quokka M1 / Sₙ) is deployed.
- **Decision on sequencing.** The EOS-only floor (≈ 0.57) is the trustworthy worst-case-so-far and a
  lower bound on `f`'s EOS side; the radiative leak can only push it deeper. So the real-opacity rung
  is a **refinement of a known floor, not a gate** — if ≈ 0.57 already clears a useful `f`, the radiative
  dip is a quantification exercise, not a blocker.
- **Package seam.** The sweep starts at ~5 km/s: below it the high-v `eos_water` dissociation chemistry
  degrades and the low-v CoolProp/two-phase package (Rung C) takes over. The transitional sweep is the
  high-v package run down to that ~5 km/s validity seam.

**Amendment (Rung E, 2026-06): the dip is *not* radiatively fillable by the ablating shield.** Rung E
(ADR-0014) tested the ablating wall as the candidate lever to fill this dip. It does **not**: the
0.57 worst case is an **EOS specific-heat sink** (enthalpy into dissociation/ionization), not a wall
radiative loss, so the vapor shield has almost nothing to recover there — ablating lifts the 11 km/s
`e_eff` only to **[0.570, 0.580]** (recovery ≤ +0.014). This *sharpens* the "refinement, not a gate"
decision above: the real-opacity radiative-leak quantification matters for the rigid floor's **loss
accounting** (channel 1a), but it cannot lift the dip's `e_eff` via ablation, because the dip's loss
is in the EOS, not in the radiation field. The 0.57 EOS floor stands as the worst case; the geometry
(concave `eta_capture`) lever, not the wall, is what clears a useful `f` there (Rung D, `f ≈ 0.83`).
