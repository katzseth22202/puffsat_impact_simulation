"""The Sn transport check's verdict on flux-limited diffusion (Q9/Q21, ADR-0012).

`hydro1d::transport` audits FLD's escape-to-space channel against a gray discrete-ordinates solve
on the same states, and `sweep --transport-check` writes one row per `(v, rho)`. What comes back is
a **relative bias on a loss channel**, which is not yet an answer to anything: a 40% bias on a
channel carrying 0.01% of the energy budget cannot move `f`, while a 5% bias on a channel carrying
half of it certainly can.

This module does two things the raw rows cannot:

1. **Translates the bias into `delta_f`**, the currency the study reports and the same one the Q18
   opacity bracket lands in, so the two radiation uncertainties can be compared directly rather
   than described.
2. **Applies the escalation gate** from design SS12.1 step 5 -- a bias above 10% escalates to a
   coupled M1 solver -- verbatim on the bias, so the decision cannot drift into "whatever we judge
   small once weighted".

Both numbers are reported. They answer different questions and neither substitutes for the other.

**Limitation, carried into every output:** the audit is one-way. It reports what transport would
have said about the escaping flux on states FLD produced; it never feeds back, so `delta_f` here is
a *bias estimate*, not a correction. Removing that caveat requires replacing the radiation operator,
which is precisely what the gate decides for or against.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

# Design SS12.1 step 5: above this relative disagreement the FLD closure is no longer adequate and
# the study escalates to a coupled M1 solver. Held as a named constant so the gate is one decision
# in one place rather than a literal scattered through the reporting.
ESCALATION_THRESHOLD = 0.10

# eta_capture along the survivable cloud schedule (design SS7; the geometry sweep gives 0.967-0.992
# there). Same nominal the opacity bracket uses, so the two delta_f numbers are commensurable.
ETA_NOMINAL = 0.98

DEFAULT_SWEEP_PATH = Path("data/results/sweep_transport_check.jsonl")
DEFAULT_SUMMARY_PATH = Path("data/results/transport_check.csv")

CSV_HEADER = (
    "v,rho_impact,length,e_eff,tau_flux_weighted,escape_share_of_ke,"
    "bias_rosseland,bias_planck,mean_selection_spread,worst_bias,worst_tau,"
    "delta_e_eff,delta_f,verdict"
)


@dataclass(frozen=True)
class AuditRow:
    """One audited `(v, rho)` case, as the Rust sweep wrote it."""

    v: float
    rho_impact: float
    length: float
    e_eff: float
    loss_escape_space: float
    transport_escape_rosseland: float
    transport_escape_planck: float
    relative_bias: float
    relative_bias_planck: float
    worst_relative_difference: float
    worst_optical_depth: float
    flux_weighted_optical_depth: float
    escape_share_of_ke: float
    converged: bool

    @property
    def mean_selection_spread(self) -> float:
        """How far the Planck-mean tally sits from the Rosseland one.

        A gray transport solve carries a single extinction coefficient while FLD is a two-mean
        model (ADR-0006: Planck for emission/absorption, Rosseland for diffusion), so **no
        single-mean Sn run is "FLD minus the closure"**. Running both brackets that ambiguity: a
        small spread means the mean choice is immaterial here and the Rosseland number is a clean
        closure verdict; a large one means part of the disagreement is about opacity means, and
        blaming the closure for it would be wrong.
        """
        return abs(self.relative_bias_planck - self.relative_bias)


@dataclass(frozen=True)
class FImpact:
    """A bias on the escape channel, translated into the study's deliverable."""

    delta_e_eff: float
    delta_f: float


def delta_f_from_bias(
    bias: float, escape_share: float, e_eff: float, eta_capture: float = ETA_NOMINAL
) -> FImpact:
    """Translate a relative bias on the escape channel into a band on `f`.

    A bias `beta` on a channel carrying fraction `phi` of the incident kinetic energy misplaces
    `phi*beta` of the energy budget. The rebound carries `e_eff**2` of the incident KE, so

    ```
    d(e_eff)/e_eff = (1/2) * d(KE)/KE  =>  d(e_eff) = phi*|beta| / (2*e_eff),
    ```

    and `f = eta*(1 + e_eff)/2` carries that to `d(f) = eta*d(e_eff)/2`.

    **First-order, and an estimate.** It books the misplaced energy but not its effect on the
    *pressure history* -- radiation lost earlier or later changes when the wall is loaded, not only
    how much energy remains. It is the right order and the right sign; it is not a correction.

    A zero or negative `e_eff` (no rebound at all) has no restitution to perturb, so the translation
    is undefined and returns zero rather than dividing by it.
    """
    if e_eff <= 0.0:
        return FImpact(delta_e_eff=0.0, delta_f=0.0)
    delta_e = escape_share * abs(bias) / (2.0 * e_eff)
    return FImpact(delta_e_eff=delta_e, delta_f=eta_capture * delta_e / 2.0)


def verdict(bias: float) -> str:
    """`ESCALATE` when the disagreement exceeds the gate, `PASS` otherwise.

    Read on the raw bias, as design SS12.1 step 5 states it, and on its magnitude: FLD
    over-reporting the escape is exactly as much a closure failure as under-reporting it.
    """
    return "ESCALATE" if abs(bias) > ESCALATION_THRESHOLD else "PASS"


@dataclass(frozen=True)
class RankedRow:
    """An audited row with its translated impact on `f`."""

    row: AuditRow
    impact: FImpact
    verdict: str


def load(path: Path) -> tuple[list[AuditRow], list[AuditRow]]:
    """Read the audit JSONL, returning `(usable, rejected)`.

    A row whose bounce never converged is rejected, not judged: since Q6 the contract is that
    `converged = false` means *no result*, not a low one. Its escape integral is truncated wherever
    the run stopped, so any bias computed from it is a solver artifact wearing a physics answer's
    clothes.
    """
    usable: list[AuditRow] = []
    rejected: list[AuditRow] = []
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        d = json.loads(line)
        row = AuditRow(
            v=float(d["v"]),
            rho_impact=float(d["rho_impact"]),
            length=float(d["length"]),
            e_eff=float(d["e_eff"]),
            loss_escape_space=float(d["loss_escape_space"]),
            transport_escape_rosseland=float(d["transport_escape_rosseland"]),
            transport_escape_planck=float(d["transport_escape_planck"]),
            relative_bias=float(d["relative_bias"]),
            relative_bias_planck=float(d["relative_bias_planck"]),
            worst_relative_difference=float(d["worst_relative_difference"]),
            worst_optical_depth=float(d["worst_optical_depth"]),
            flux_weighted_optical_depth=float(d["flux_weighted_optical_depth"]),
            escape_share_of_ke=float(d["escape_share_of_ke"]),
            converged=bool(d["converged"]),
        )
        (usable if row.converged else rejected).append(row)
    return usable, rejected


def rank(rows: list[AuditRow], eta_capture: float = ETA_NOMINAL) -> list[RankedRow]:
    """Rank by translated `delta_f`, worst first.

    Deliberately **not** by raw bias. Ranking on the bias would head the table with whichever
    energetically irrelevant channel happens to disagree most and bury the row that actually
    threatens `f`.
    """
    ranked = [
        RankedRow(
            row=r,
            impact=delta_f_from_bias(r.relative_bias, r.escape_share_of_ke, r.e_eff, eta_capture),
            verdict=verdict(r.relative_bias),
        )
        for r in rows
    ]
    ranked.sort(key=lambda x: x.impact.delta_f, reverse=True)
    return ranked


def write_summary(
    sweep_path: Path = DEFAULT_SWEEP_PATH, out_path: Path = DEFAULT_SUMMARY_PATH
) -> list[RankedRow]:
    """Reduce the audit to the verdict table; return the ranked rows."""
    rows, rejected = load(sweep_path)
    ranked = rank(rows)
    lines = [CSV_HEADER]
    for item in ranked:
        r = item.row
        lines.append(
            f"{r.v},{r.rho_impact},{r.length},{r.e_eff:.6f},"
            f"{r.flux_weighted_optical_depth:.6e},{r.escape_share_of_ke:.6e},"
            f"{r.relative_bias:.6e},{r.relative_bias_planck:.6e},"
            f"{r.mean_selection_spread:.6e},{r.worst_relative_difference:.6e},"
            f"{r.worst_optical_depth:.6e},{item.impact.delta_e_eff:.6e},"
            f"{item.impact.delta_f:.6e},{item.verdict}"
        )
    out_path.write_text("\n".join(lines) + "\n")
    if rejected:
        for r in rejected:
            print(
                f"python: WARNING rejected non-converged row v={r.v / 1000:.0f} km/s "
                f"rho={r.rho_impact:.3f} (Q6 contract: no result, not a low one)"
            )
    return ranked


def main() -> None:
    """CLI: `python -m puffsat.transport_check` -> data/results/transport_check.csv."""
    ranked = write_summary()
    if not ranked:
        print("python: no usable audit rows")
        return

    escalate = [x for x in ranked if x.verdict == "ESCALATE"]
    worst_tau = min(x.row.flux_weighted_optical_depth for x in ranked)
    print("python: Sn transport check of the FLD escape channel (Q9/Q21)")
    print(f"    rows: {len(ranked)}   escape-weighted tau range down to {worst_tau:.2f}")
    print(
        f"    max |bias| (Rosseland): {max(abs(x.row.relative_bias) for x in ranked):.2%}"
        f"   gate: {ESCALATION_THRESHOLD:.0%}"
    )
    print(f"    max mean-selection spread: {max(x.row.mean_selection_spread for x in ranked):.2%}")
    top = ranked[0]
    print(
        f"    worst by delta_f: v={top.row.v / 1000:.0f} km/s rho={top.row.rho_impact:.3f} "
        f"-> delta_f={top.impact.delta_f:.3e} "
        f"(bias {top.row.relative_bias:+.2%} on {top.row.escape_share_of_ke:.3%} of KE)"
    )
    if escalate:
        print(f"    VERDICT: ESCALATE — {len(escalate)} row(s) above the gate:")
        for x in escalate:
            print(
                f"      v={x.row.v / 1000:.0f} km/s rho={x.row.rho_impact:.3f} "
                f"bias={x.row.relative_bias:+.2%} tau={x.row.flux_weighted_optical_depth:.2f}"
            )
    else:
        print("    VERDICT: PASS — FLD stands; the bias is the radiation-model error bar.")
    print(f"python: wrote {DEFAULT_SUMMARY_PATH}")


if __name__ == "__main__":
    main()


__all__ = [
    "CSV_HEADER",
    "ESCALATION_THRESHOLD",
    "ETA_NOMINAL",
    "AuditRow",
    "FImpact",
    "RankedRow",
    "delta_f_from_bias",
    "load",
    "rank",
    "verdict",
    "write_summary",
]
