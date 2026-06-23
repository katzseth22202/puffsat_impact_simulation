# 2D geometry sweep uses a calibrated effective-gamma EOS, spot-checked against equilibrium EOS

`eta_capture` is the lossless 2D/1D wall-impulse ratio (ADR-0003), and the 2D geometry sweep runs
hundreds of (curvature × cloud-shape) cases. That sweep uses a **calibrated effective-gamma** EOS,
not the full tabulated equilibrium EOS, with an equilibrium-EOS spot-check guarding the approximation.

**Rationale.** Because `eta_capture` is a ratio with the same EOS in numerator and denominator, EOS
error is largely **common-mode and cancels** — a calibrated effective-gamma can give the right
*ratio* even where it gets neither `J_wall` individually right. So the full equilibrium EOS (a table
lookup per cell per step, × hundreds of runs) is an unnecessary cost for the sweep.

**The imperfection it guards against.** Cancellation is not exact: in 2D the gas expands *sideways*
into growing volume and samples a different thermodynamic path (more cooling/recombination) than the
1D axial re-expansion, so the EOS does not perfectly divide out. The guard is a **full-equilibrium-EOS
spot-check at the optimum, plus one flat and one shallow-concave case**. If the ratio shifts beyond
the reporting tolerance, the sideways-path effect is non-negligible and the sweep falls back to the
equilibrium EOS.

**Non-negotiable.** The 2D run and its 1D-lossless denominator use the *identical* EOS (whichever),
radiation-free and adiabatic, with `J_wall` extracted by the same `10⁻³` force-decay cutoff as
ADR-0001. The cost optimization is licensed by the spot-check, not assumed.

## Amendment (2026-06): the effective-gamma kernel landed; spot-check still deferred

The 2D kernel (`crates/euler2d`, ADR-0023) is built with the **calibrated effective-gamma ideal gas
baked in** (`p = (γ−1)ρe`, `c = √(γ(γ−1)e)`), exactly as this ADR licenses — no per-cell equilibrium
table lookup. The flat-plate `eta_capture` (0.81 → 0.92 over `r_foot/L = 0.5 → 2`, ADR-0003) is
reported under this EOS. Both the free and the confined-2D denominator runs use the *identical* γ-law,
radiation-free and adiabatic, with `J_wall` extracted by the same `10⁻³` cutoff — so the EOS error is
common-mode by construction at this stage.

**The equilibrium-EOS spot-check remains deferred to the optimum**, per this ADR: the sideways-path
residual (2D expansion samples a different thermodynamic path than the 1D denominator) is only worth
measuring at the operating point, and the concave plate + parametric sweep that the spot-check guards
are themselves the follow-on rung. Nothing landed here changes that deferral.

## Considered Options

- **Full equilibrium EOS across the whole 2D sweep.** Rejected as baseline: pays a large per-run cost
  to remove doubt that the ratio's common-mode cancellation already removes for all but the
  sideways-path residual — which the spot-check measures directly.
- **Effective-gamma with no spot-check.** Rejected: leaves the sideways-path cancellation error
  unquantified.
