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

## Amendment (2026-06, Rung S): the gate evaluated against the *survivable* frontier

The survivability frontier (design §7, ADR-0010/0011) now lets the gate be evaluated against the
*survivable* `f`, not the unconstrained best — the honest reading, since an operating point that
exceeds the facesheet pressure limit is not available:

- **At the transitional dip (the conservative worst case, `e_eff ≈ 0.57`, 11 km/s):** the best
  survivable `f ≈ 0.80` — it **clears the gate**, on a moderate-footprint shallow-concave disk that
  peaks at ~370 MPa (under the 400 MPa baseline) even after a 1.78× concave focusing penalty. The
  unconstrained `f`-max corner (short disk + tight footprint) is foreclosed: it peaks at ~2.3 GPa.
- **At 16 km/s (`e_eff ≈ 0.63`):** the best survivable `f ≈ 0.78` under the 400 MPa **baseline** —
  *just below* the gate — rising to `≈ 0.84` at the **relaxed 900 MPa** high-v limit (the §7 sweep
  to 700/900). So at 16 km/s the verdict is limit-conditional: it clears 0.8 only if the SiC+Ti stack
  carries the relaxed pressure, which is exactly the kind of one-sided refinement the gate triggers.

**Consequence for the lever-gating logic above:** the conservative baseline clears 0.8 at the dip
(settled there) but sits marginally below at 16 km/s under the baseline pressure limit — so the
16 km/s anchor is where the recovery levers (ablating wall, ADR-0014; the relaxed `P_limit`) earn
their keep. This is still a conservative-floor reading (rigid wall, worst-case `e_eff`, the concave
focusing penalty applied); the best-estimate curve (ADR-0013) only moves it up. No `f` is quoted
externally until validation (design §9). *(Superseded — see the 2026-06 paper-draft amendment
below: a preliminary, conservative-floor `f` is quotable in the white paper.)*

## Amendment (2026-06, paper draft): a preliminary `f` is quotable in the white paper

The original stance — "no `f` is quoted externally until the hydrocode cross-check lands"
(design §9) — is **relaxed for the white paper**. The paper is an explicitly speculative concept
paper, not a validated engineering deliverable, so a preliminary `f` is quotable *provided it is
framed as preliminary*. The bar for "don't mislead readers" is met by labeling, not by withholding.

Conditions on any externally quoted `f`:
- Labeled **preliminary** and attributed to a **single in-house code**, with the FLASH cross-check
  named as the outstanding independent validation (ADR-0009 above / design §9).
- Quoted as the **conservative-floor** reading *with its sensitivity* (the linear mass-ratio rule
  above; the two-curve deliverable, ADR-0013) — never as a guaranteed, optimized, or
  margin-carrying value.
- Carries the velocity-dependence caveat: the conservative worst case is the **interior ~11 km/s
  dip**, not an endpoint (ADR-0012), and 16 km/s is limit-conditional (Rung S amendment above).

The FLASH cross-check **remains the open validation gate**: it upgrades the quoted number from
*preliminary* to *validated*. It no longer blocks external quotation — it gates only the
"validated" label.

## Considered Options

- **A performance/mission-derived threshold.** Rejected: drags in the out-of-scope vehicle analysis
  (§11).
- **A margin above the 0.5 inelastic floor, or a 0.7–0.8 band.** Rejected as the primary anchor: 0.8
  is the specific existing claim the analysis is meant to defend, which makes the verdict directly
  meaningful; the floor and band are conveyed via the reported sensitivity instead.
