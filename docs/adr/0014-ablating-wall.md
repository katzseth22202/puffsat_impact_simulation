# The ablating wall is a Phase-2 e_eff lever, not just protection; the rigid wall is the conservative floor

The survivability frontier runs on a passive ablating/transpiring wall (Phase 2). Counter to the
intuition that ablation only "costs" energy, it tends to *raise* `e_eff` — so the rigid wall
(ADR-0005) is the conservative floor and the ablating wall is the best-estimate refinement.

**Why ablation raises `e_eff`.** The vapor acts on two loss channels: **blowing** (injected vapor
thickens the boundary layer and cuts the conductive flux into the solid) and **vapor shielding**
(the vapor absorbs incoming radiation in the near-wall layer before it reaches the cold solid).
Both keep energy *in the gas* instead of losing it to the wall. The only true energy sink is the
ablation enthalpy `Q*·ṁ`, and `ṁ` is small — a few µm of ablator ≈ 0.4 kg vs the 25 kg pulse,
~1.5%. The recovery outweighs the cost, so `e_eff` rises relative to the rigid wall.

**Ties to the dual-curve and the transition.** Rigid wall = conservative floor for *both* `e_eff`
and survivability (full wall losses, full flux to the SiC); ablating wall = best-estimate. Per
ADR-0013, ablation is a gated refinement (rigid → lower-bound curve, ablating → best-estimate curve),
run where the floor dips below `f = 0.8` (ADR-0009). Its `e_eff` leverage should *peak at the
transitional anchor*: at 3.2 km/s the gas barely radiates, at 16 km/s the radiation is already
trapped (`τ≫1`), but at `τ~1` radiation reaches the wall and vapor shielding has the most to recover
— so the ablating wall is the prime candidate to fill the transitional `e_eff` dip (ADR-0012).

**Model.** Quasi-steady surface energy balance: incoming flux → ablation rate via an effective heat
of ablation `Q*`, with vapor injected as a **mass/momentum/energy source at the 1D wall boundary**.
Three coupled effects, all self-consistent in the existing 1D rad-hydro: blowing correction on the
conductive flux; vapor opacity in the FLD for shielding; injected mass in the bounce. **Passive
ablation** of the MEMS-renewed sacrificial layer — not pumped coolant; renewal is between pulses,
passive during the pulse.

**Vapor approximated, not speciated (baseline).** The vapor is treated as a `Q*` enthalpy sink + a
gray near-wall absorber + a mass source, without full silicone/ionic-liquid pyrolysis chemistry —
even though the ablating wall injects ablator species the water-only table pipeline (ADR-0007) does
not cover. Full ablator-vapor speciation is gated on whether spot-checks show it moves `e_eff`. `Q*`
is parameterized (silicone ~2–10 MJ/kg literature), with sensitivity reported.

## Considered Options

- **Transient / full-pyrolysis ablation** (char layer, decomposition kinetics, gas-phase products)
  from the start. Deferred: quasi-steady `Q*` captures the first-order energy interception; full
  pyrolysis only if spot-checks show vapor chemistry matters.
- **Treat ablation mass-loading as a first-order perturbation needing full vapor chemistry.**
  Rejected as baseline: `ṁ ~ 1.5%` of the pulse; the dominant effect is on the loss channels
  (blowing, shielding), captured with an effective `Q*` + gray absorber.
