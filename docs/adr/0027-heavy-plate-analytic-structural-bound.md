# First-cut analytic plate-structural bound: the heavy-plate scenario carries a closed-form go/no-go, decoupled from `f(v)`

The **heavy-plate** special scenario (100 kg water pulses, `R = 15 m` / 30 m-diameter plate,
plate mass **≤ 40 t** as a ceiling, swept 16–28 km/s at 0.5 km/s) makes a question the envelope
study deliberately parked: **can the plate structure itself survive the per-pulse impulse?** At
these parameters one pulse delivers `≈ (1 + e_eff)·m·v ≈ 5×10⁶ N·s` into a 30 m plate — a
structural ask an order of magnitude beyond the design's 5 m / 3–4 t reference body, so "it's a
Phase-2 detail" no longer defers cleanly. This ADR records the decision to answer it with a
**closed-form first-cut bound carried as a companion to `f(v)`**, and — as important — the explicit
limits of what that bound does and does not claim.

## The scope boundary this deviates from

Design §5, §11, and ADR-0011 scope this repository to **gas-dynamics `f(v)` plus *local facesheet*
survivability**. The whole-plate structural body (Ti truss/corrugated core, tensioned
high-strength-fiber back-face, first-mode-vs-pulse rigidity, SiC-Ti spall accumulation) is named
there as a **parallel Phase-2 track that "does not change `f(v)`," fed by the solver's
peak-load/impulse/footprint outputs** — explicitly *not* built here. The repo has two hydro
kernels and no structural solver, by design.

The heavy-plate scenario does not overturn that boundary; it carves a **narrow, scenario-scoped
exception**: a handful of closed-form checks, decoupled from the `f(v)` pipeline, that answer
"is the ≤ 40 t plate plausibly buildable and rigid-during-pulse at this scale?" without pulling a
structural-dynamics discipline (FEA / plate-shell dynamics) into the codebase.

## Method: three closed-form checks, off the solver's load outputs

Fed by the 1D kernel's `peak_wall_pressure`, per-pulse `wall_impulse`, bounce duration, and the
Σ-contract footprint — no new kernel, no mesh:

1. **Rigid-during-pulse check (first-mode period vs pulse duration).** The candidate plate's
   fundamental flexural period must stay `≫` the ~µs–100s-of-µs bounce time, or the face is not
   locally rigid during the collision and the rigid-wall assumption behind `e_eff`/`f` breaks.
   This check therefore does double duty: it is a structural gate *and* the validity test for
   `f` at this much larger span (design §5, constraint 2). Fails here invalidate the `f(v)`
   result, not just the structure.
2. **Areal-impulse → membrane/bending stress.** The per-unit-area impulse drives a candidate
   **Ti-truss core + tensioned fiber (Vectran) back-face** at the ≤ 40 t / `R = 15 m` mass
   budget (areal mass ≈ 40 t / π·15² ≈ 57 kg/m²); compare peak membrane + bending stress to the
   material allowables. The 40 t is treated as a **ceiling** — the check asks whether *some*
   admissible construction survives at ≤ 40 t, reporting the implied minimum, not sizing a final
   design.
3. **SiC-Ti spall reflection (ADR-0011).** The peak facesheet compressive load reflects at the
   lower-impedance Ti backing as tension (`R ≈ −0.15`); confirm it stays sub-spall at the
   heavy-plate peak pressures, reusing the ADR-0011 model.

## What it claims — and what it does not

- **It is a go/no-go feasibility bound with documented candidate-construction assumptions**, not a
  validated structural design. A "go" means a plausible ≤ 40 t construction clears the checks with
  margin; a "no-go" is a real negative finding for the scenario *even if `f(v)` clears the gate*,
  and would be reported as such.
- **It is closed-form, not FEA / plate-shell dynamics.** Buckling, dynamic amplification beyond
  the first mode, joint/fatigue detail, and thermal-structural coupling are out of it (and out of
  this repo).
- **It lives beside `f(v)`, not inside it.** It consumes solver outputs and emits a separate
  verdict; it never feeds back into `e_eff`, `eta_capture`, or `f` — except through check (1),
  whose *failure* would flag the `f` result as invalid at this scale (a one-way validity gate,
  not a value coupling).
- **It is scenario-scoped to heavy-plate**, not a general repo capability. Other anchors keep the
  strict §11 boundary (facesheet + back-propagated loads only).

## Considered Options

- **Strict scope: facesheet survival + back-propagate the impulse/peak-load/footprint only.**
  Rejected *for this scenario*: at 30 m / ≤ 40 t / 5×10⁶ N·s the structural question is
  first-order, not a downstream detail — deferring it would leave "does the plate survive the
  impulse?" (the user's actual question) unanswered while reporting a healthy `f`.
- **Full structural-dynamics track (FEA or reduced plate-shell model).** Rejected: a discipline
  step-change and almost certainly a separate tool/repo, disproportionate to a first-cut
  feasibility question and coupling a heavyweight new dependency into a gas-dynamics codebase.
  Left as the named refinement if the closed-form bound lands marginal.
- **New ADR vs. amend ADR-0011.** A new ADR was chosen because this is a *scope/architecture*
  decision (what disciplines the repo will and won't host, and under what scenario-scoped
  exception), broader than ADR-0011's SiC-Ti spall physics, which it merely reuses.

## Consequences

- Sets the precedent that a **special scenario may host a bounded closed-form structural check**
  when its scale makes the Phase-2 structural question first-order — kept narrow (closed-form,
  scenario-scoped, decoupled) so it does not erode the §11 boundary for the envelope study.
- The rigid-during-pulse check becomes a **named validity gate on `f`** at large spans: if it
  fails, the heavy-plate `f(v)` is reported as provisional-pending-structure, not clean.
- Implementation (deferred): a Python cold-path module (e.g. `python/puffsat/structure.py`) fed by
  the heavy-plate sweep outputs, reported alongside the `f(v)` / facesheet-survivability frontier.
