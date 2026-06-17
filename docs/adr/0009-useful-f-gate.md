# The "useful f" gate is the paper's f = 0.8, evaluated against the conservative baseline

Several downstream decisions are gated on whether `f` is "useful": whether to build the kinetic
condensation model (ADR-0004), whether to pull the dark-oil seed lever (§7), and whether the low-v
worst case is settled. This ADR fixes that gate at a number.

**Threshold = 0.8, the paper's own assumed value.** The sim exists to replace the paper's constant
`f = 0.8` with a real `f(v)`, so the scope-respecting gate is: does the *conservative* (equilibrium)
baseline already clear 0.8?
- Conservative `f ≥ 0.8` at an anchor → that case is settled; recovery levers (kinetic condensation,
  seed) are not built.
- Conservative `f < 0.8` → the paper's assumption is at risk; pull the recovery levers to claw back
  toward 0.8 and report the residual.

**Why not a performance threshold.** The natural definition — "high enough that the propulsion
concept closes" — is a mass-ratio/mission question, and all vehicle-performance analysis is out of
scope (§11). Anchoring to 0.8 keeps the gate inside scope: it is the existing claim the analysis
defends, not a new requirement. §6 already ties 0.8 to `eta_capture ~ 0.95`, `e_eff ~ 0.68`.

**Always report the sensitivity, not just the verdict.** The mass-ratio law `m_r/m_p = 2f/ln(...)`
is linear in `f`, so `f(v)` is reported with the rule "an X% shortfall in `f` is an X% increase in
required mass ratio." That hands any future performance study the lever without this sim defining
the mission.

## Considered Options

- **A performance/mission-derived threshold.** Rejected: drags in the out-of-scope vehicle analysis
  (§11).
- **A margin above the 0.5 inelastic floor, or a 0.7–0.8 band.** Rejected as the primary anchor: 0.8
  is the specific existing claim the analysis is meant to defend, which makes the verdict directly
  meaningful; the floor and band are conveyed via the reported sensitivity instead.
