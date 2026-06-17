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

## Considered Options

- **Full equilibrium EOS across the whole 2D sweep.** Rejected as baseline: pays a large per-run cost
  to remove doubt that the ratio's common-mode cancellation already removes for all but the
  sideways-path residual — which the spot-check measures directly.
- **Effective-gamma with no spot-check.** Rejected: leaves the sideways-path cancellation error
  unquantified.
