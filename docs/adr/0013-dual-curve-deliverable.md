# The deliverable is two curves: a conservative lower-bound f(v) and a best-estimate f(v) with error bars

`f(v)` is reported as **two curves**, not one, because almost every modeling approximation in the
study was deliberately chosen one-sided (conservative). That structure lets the strongest claim land
early and keeps one-sided conservatism cleanly separated from genuine two-sided uncertainty.

**Conservative lower-bound `f(v)`** — all one-sided approximations stacked pessimistically:
- chemical energy zeroed (§4): `f ≤ 1` strictly;
- equilibrium condensation (ADR-0004): over-condenses, lowers `e_eff`;
- `P_limit = 400 MPa` floor (frontier / ADR-0009): worst-case facesheet strength;
- effective-gamma 2D sweep (ADR-0008): bounded by the equilibrium-EOS spot-check;
- factorization separability (ADR-0003): bounded by the cross-code.

This is the headline, defensible number — "even at worst, `f ≥` this" — and it is **quotable as a
floor as soon as the smoke tests (ADR-0001 limits) and the Orion validation pass**, without first
resolving every two-sided uncertainty.

**Best-estimate `f(v)` with error bars** — from a one-at-a-time sensitivity study over the genuine
two-sided uncertainties: opacity at partial ionization (ADR-0006/0007), FLD at `τ~1` (ADR-0012,
transport check at the transition), EOS-seam continuity (ADR-0007), and the `J_wall` integration
cutoff (ADR-0001).

**The gates are the error-budget logic.** The gated refinements — kinetic condensation (ADR-0004),
multigroup (ADR-0006), equilibrium-EOS spot-check (ADR-0008), transport check (ADR-0012) — are
pulled **only where the conservative floor dips below the `f = 0.8` useful-f line** (ADR-0009).
Modeling effort is spent exactly where the pessimistic floor fails to clear the bar, and nowhere else.

## Considered Options

- **A single best-estimate `f(v)` with symmetric error bars.** Rejected: it discards the defensible
  floor — the strongest and earliest-available claim, needing only the one-sided approximations plus
  Orion — and blurs one-sided conservatism into two-sided uncertainty that it is not.
