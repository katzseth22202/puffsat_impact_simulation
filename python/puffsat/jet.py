"""N2 -- `eta_jet`, and whether the flare actually turns the pancake (both legs).

Rung 4 found the merged plume is a **pancake**: `alpha = <v_z^2>/<v^2>` = 0.088 against the 1/3
`eq:reflection_baseline` assumes, because the bag became a 23 m column for launch-fairing reasons
and free expansion follows the short axis. That halves the *mirror* baseline, to 0.237.

**A mirror is not the design.** `sec:jet_efficiency` says so itself -- "a real nozzle turns flow
instead of merely reflecting it ... that is why `eq:reflection_baseline` is a reference point and
not a ceiling". This module computes what the turning actually buys, which is N2, and which the
paper says nothing in it or either companion bounds.

# The mechanism, and why a pancake is not the disaster the baseline implies

A mirror is helpless against transverse motion: it can only reverse what already points at it. A
**diverging field is not**, and the reason is the adiabatic invariant

    mu = v_perp^2 / (2 B) = const,

so as the plume runs down the flare and `B` falls, `v_perp^2` falls *in proportion*, and energy
conservation puts the difference into `v_par`. With `alpha_0` the axial share at the throat and
`B_exit/B_throat = 1/(A/A*)` by flux conservation,

    alpha_exit = 1 - (1 - alpha_0) * B_exit/B_throat.

**The flare is therefore the exact remedy for the pancake**, and the deeper the flare the more
complete the conversion. This is the standard magnetic-nozzle result, not a new claim; what is new
is applying it to a plume that is measured rather than assumed.

# The two legs are different devices

`sec:two_leg_nozzle`: "the cycle collides twice, and the two collisions face opposite ways."

- **Leg 2, head-on departure** -- a *nozzle*. The drift already points out the back, and the flare
  converts `v_perp` to `v_par` on the way out. This is the case above.
- **Leg 1, overtake** -- a *mirror*, "a cup closed at the ship end". The drift points prograde, the
  wrong way, and the field must reverse it before it is worth anything. A mirror of ratio `R`
  reflects everything except its loss cone, `sin^2(theta) < 1/R`.

And here the pancake **helps**: mirror reflection acts on the perpendicular component, so a plume
with `alpha = 0.088` has `sin^2(theta) ~ 0.91` and sits far outside any plausible loss cone. The
same anisotropy that costs leg 2 its baseline buys leg 1 its reflection.

# What this is not

A guiding-centre argument on a fluid. It assumes the invariant holds (checked by
`adiabaticity_parameter`), that the flow detaches rather than following field lines home (Rung 2's
question, and the reason `A/A*` has to reach the field window), and that the expansion is
collisionless enough for `mu` to mean anything. It replaces an unbounded parameter with a bounded
one; it is not a solved MHD nozzle.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from puffsat import field

MU0 = field.MU0

DEFAULT_OUTPUT = Path("data/results/nozzle_jet.csv")

ALPHA_MEASURED = 0.088
"""Axial share of the freely expanded merged fireball (Rung 4, `euler2d::merge`).

The 23 m column's value. The paper's original 5.4 m spherical bag gives 0.360, near the isotropic
1/3 it assumes -- the difference is shape alone.
"""

ALPHA_ISOTROPIC = 1.0 / 3.0
"""What `eq:reflection_baseline` assumes, and what a spherical bag would deliver."""

FLOWN_AREA_RATIO = 4.0
"""`A/A*` the flown design carries (`expansion.AREA_RATIO_EXIT`), i.e. 20 T -> 5 T."""

WINDOW_AREA_RATIO = 11.3
"""Lower edge of the field window from Rung 2 -- the flare that also buys detachment."""

TARGET_ETA_JET = 0.775
"""The paper's own target (`sec:jet_efficiency`), 147% of its 0.529 mirror baseline."""

ETA_CHEM = {45.58: 0.731, 56.53: 0.835, 65.0: 0.864, 75.0: 0.910}
"""Frozen-chemistry ceiling per closing speed (`toll.eta_chem`, tab:eta_chem_speed)."""

Leg = Literal["head-on", "overtake"]


# ---- The adiabatic conversion ---------------------------------------------------------------


def alpha_after_expansion(alpha_0: float, area_ratio: float) -> float:
    """Axial share after a flare of `A/A*`, by conservation of `mu = v_perp^2/2B`.

    `B_exit/B_throat = 1/(A/A*)` (flux conservation), so `v_perp^2` is divided by the area ratio
    and everything it loses becomes axial:

        alpha_exit = 1 - (1 - alpha_0)/(A/A*).

    Note what this says about the pancake: `alpha_0` enters only through `(1 - alpha_0)`, which is
    at most 1, so **a deep enough flare recovers a directed jet from any starting anisotropy**.
    The flare is not a mitigation for the pancake, it is the mechanism that undoes it.
    """
    if not 0.0 <= alpha_0 <= 1.0:
        raise ValueError("alpha_0 must lie in [0, 1]")
    if area_ratio < 1.0:
        raise ValueError("area ratio must be at least 1")
    return 1.0 - (1.0 - alpha_0) / area_ratio


def area_ratio_for_alpha(alpha_0: float, alpha_target: float) -> float:
    """The flare needed to reach a stated axial share -- the inverse of the above."""
    if alpha_target >= 1.0:
        return math.inf
    if alpha_target <= alpha_0:
        return 1.0
    return (1.0 - alpha_0) / (1.0 - alpha_target)


def eta_geom_directed(alpha: float) -> float:
    """`eta_geom` for a **directed** exhaust of axial share `alpha`: `sqrt(alpha)`.

    A nozzle's exhaust all leaves the same way, so the axial speeds do not cancel. With
    `v_g = sqrt(<v^2>)` and every element at axial share `alpha`, `<v_z>/v_g = sqrt(alpha)`.

    The optimistic bound: it assumes a single speed rather than a spread.
    """
    return math.sqrt(max(alpha, 0.0))


def eta_geom_spread(alpha: float) -> float:
    """`eta_geom` for a Gaussian spread of the same axial share: `sqrt(2 alpha/pi)`.

    The same family `eq:reflection_baseline` uses, so at `alpha = 1/3` it returns the paper's
    0.461. The pessimistic bound: a real nozzle's exhaust is directed, so the truth sits between
    this and `eta_geom_directed`, and the gap is the `sqrt(2/pi)` = 0.798 of a Maxwellian spread.
    """
    return math.sqrt(2.0 * max(alpha, 0.0) / math.pi)


def adiabaticity_parameter(b_field: float, length_scale: float, temp_k: float = 15_000.0) -> float:
    """`r_gyro/L` for a proton at the stated field -- `mu` is invariant only while this is small.

    The invariant survives when the field varies slowly over a gyroradius. A thermal proton at
    `T` has `r_g = m v_th/(q B)`; dividing by the length over which `B` changes gives the
    adiabaticity parameter. Values well below 1 mean the guiding-centre picture holds; approaching
    1 means it does not and this whole module stops applying.
    """
    m_p, q_e, k_b = 1.67262192e-27, 1.602176634e-19, 1.380649e-23
    v_th = math.sqrt(2.0 * k_b * temp_k / m_p)
    r_gyro = m_p * v_th / (q_e * b_field)
    return r_gyro / length_scale


# ---- The overtake leg: a mirror, not a nozzle -------------------------------------------------


def loss_cone_fraction(mirror_ratio: float) -> float:
    """Share of an isotropic population that escapes a mirror of ratio `R`: `sqrt(1 - 1/R)`.

    Particles escape when their pitch angle lies inside the loss cone, `sin^2(theta) < 1/R`.
    The complement is what reflects, which on the overtake leg is what becomes thrust.
    """
    if mirror_ratio <= 1.0:
        return 1.0
    return math.sqrt(1.0 - 1.0 / mirror_ratio)


def mean_sin2_pitch(alpha: float, drift_fraction: float) -> float:
    """Mean `sin^2(theta)` of the plume against an axial field, **including the bulk drift**.

    The drift matters and cuts against the plume: it is motion *along* the field, so it pushes the
    pitch angle toward the loss cone, which is the one direction a mirror cannot reverse. With
    `V/u = sqrt(f_d/(1-f_d))` and the thermal part split by `alpha`,

        v_par ~ V + sqrt(alpha) u,   v_perp ~ sqrt(1-alpha) u.

    At the flown `alpha` = 0.088 and `f_d` = 0.105 this gives `sin^2(theta)` ~ 0.69, still far
    outside any plausible loss cone -- but 0.91 would have been the answer had the drift been
    ignored, so the correction is worth a quarter of the margin.
    """
    a = max(min(alpha, 1.0), 0.0)
    f = max(min(drift_fraction, 0.999), 0.0)
    v_drift = math.sqrt(f / (1.0 - f))
    v_par = v_drift + math.sqrt(a)
    v_perp2 = 1.0 - a
    return v_perp2 / (v_perp2 + v_par * v_par)


def reflected_fraction(alpha: float, mirror_ratio: float, drift_fraction: float = 0.105) -> float:
    """Share of the plume a mirror of ratio `R` turns around, on the overtake leg.

    Mirror reflection acts on the *perpendicular* component, so a pancake reflects readily: it
    escapes only inside the loss cone, `sin^2(theta) < 1/R`.

    **The anisotropy that costs leg 2 its baseline buys leg 1 its reflection.** Modelled on the
    mean pitch angle rather than the full distribution, so it saturates at 1 -- it is a bound
    showing the leg-1 mirror is not the binding constraint, not a mirror design.
    """
    if mirror_ratio <= 1.0:
        return 0.0
    return 1.0 if mean_sin2_pitch(alpha, drift_fraction) > 1.0 / mirror_ratio else 0.0


# ---- Assembled per-leg efficiency -------------------------------------------------------------


@dataclass(frozen=True)
class JetEfficiency:
    """`eta_jet` on one leg at one flare, with the bracket its model choice leaves."""

    leg: Leg
    closing_speed_km_s: float
    area_ratio: float
    alpha_in: float
    alpha_out: float
    eta_geom_lo: float
    eta_geom_hi: float
    eta_chem: float
    eta_jet_lo: float
    eta_jet_hi: float
    reflected: float
    """Leg 1 only: share the mirror turns around. 1.0 on leg 2, which needs no reversal."""

    @property
    def clears_target(self) -> bool:
        """Whether the optimistic bound reaches the paper's own 0.775."""
        return self.eta_jet_hi >= TARGET_ETA_JET


def jet_efficiency(
    *,
    leg: Leg = "head-on",
    closing_speed: float = 75.0,
    area_ratio: float = FLOWN_AREA_RATIO,
    alpha_0: float = ALPHA_MEASURED,
    mirror_ratio: float = 4.0,
) -> JetEfficiency:
    """`eta_jet = eta_chem * eta_geom` on one leg, with `eta_geom` from the adiabatic conversion.

    The two `eta_geom` bounds are the same axial share read through a directed exhaust
    (`sqrt(alpha)`) and through a Gaussian spread (`sqrt(2 alpha/pi)`). A real nozzle is directed,
    so the truth sits nearer the upper bound; both are reported because nothing here resolves the
    exhaust's speed distribution.
    """
    alpha_out = alpha_after_expansion(alpha_0, area_ratio)
    lo, hi = eta_geom_spread(alpha_out), eta_geom_directed(alpha_out)
    chem = ETA_CHEM.get(round(closing_speed, 2), 0.910)
    reflected = 1.0 if leg == "head-on" else reflected_fraction(alpha_0, mirror_ratio)
    return JetEfficiency(
        leg=leg,
        closing_speed_km_s=closing_speed,
        area_ratio=area_ratio,
        alpha_in=alpha_0,
        alpha_out=alpha_out,
        eta_geom_lo=lo * reflected,
        eta_geom_hi=hi * reflected,
        eta_chem=chem,
        eta_jet_lo=chem * lo * reflected,
        eta_jet_hi=chem * hi * reflected,
        reflected=reflected,
    )


def flare_sweep(
    area_ratios: tuple[float, ...] = (2.0, 4.0, 6.0, 8.0, 11.3, 14.8, 20.0),
    *,
    alpha_0: float = ALPHA_MEASURED,
) -> list[tuple[float, float, float, float]]:
    """`(A/A*, alpha_out, eta_geom_lo, eta_geom_hi)` -- the N2 answer in one table."""
    return [
        (
            ar,
            alpha_after_expansion(alpha_0, ar),
            eta_geom_spread(alpha_after_expansion(alpha_0, ar)),
            eta_geom_directed(alpha_after_expansion(alpha_0, ar)),
        )
        for ar in area_ratios
    ]


def write_jet(rows: list[JetEfficiency], path: Path = DEFAULT_OUTPUT) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "leg",
                "closing_speed_km_s",
                "area_ratio",
                "alpha_in",
                "alpha_out",
                "eta_geom_lo",
                "eta_geom_hi",
                "eta_chem",
                "eta_jet_lo",
                "eta_jet_hi",
                "reflected",
            ]
        )
        for r in rows:
            writer.writerow(
                [
                    r.leg,
                    f"{r.closing_speed_km_s:g}",
                    f"{r.area_ratio:g}",
                    f"{r.alpha_in:.4f}",
                    f"{r.alpha_out:.4f}",
                    f"{r.eta_geom_lo:.4f}",
                    f"{r.eta_geom_hi:.4f}",
                    f"{r.eta_chem:.4f}",
                    f"{r.eta_jet_lo:.4f}",
                    f"{r.eta_jet_hi:.4f}",
                    f"{r.reflected:.4f}",
                ]
            )


def main() -> None:
    """Does the flare turn the pancake, on both legs?"""
    parser = argparse.ArgumentParser(description="N2: eta_jet and the flare")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    print("== N2: does the flare turn the pancake? ==")
    print(f"measured alpha at the throat: {ALPHA_MEASURED} (Rung 4; a sphere would give 0.360)")
    print("mu = v_perp^2/2B is invariant, so alpha_exit = 1 - (1 - alpha_0)/(A/A*)\n")
    print(f"{'A/A*':>7} {'alpha_out':>10} {'eta_geom (spread)':>18} {'eta_geom (directed)':>20}")
    for ar, a_out, lo, hi in flare_sweep():
        mark = (
            "  <- flown"
            if ar == FLOWN_AREA_RATIO
            else ("  <- field window" if abs(ar - WINDOW_AREA_RATIO) < 0.1 else "")
        )
        print(f"{ar:7.1f} {a_out:10.4f} {lo:18.4f} {hi:20.4f}{mark}")

    print("\n== eta_jet per leg, at the flown flare and at the field window ==")
    print("The two legs share eta_geom: mu conversion does not care which way the drift points,")
    print("and the leg-1 mirror reflects the pancake essentially completely. The leg asymmetry")
    print("lives in the impulse ledger's +1 (overtake) against -1 (head-on), which is")
    print("`aim_is_all_you_need`'s bookkeeping and deliberately not re-derived here.")
    rows: list[JetEfficiency] = []
    print(
        f"{'leg':>10} {'w':>7} {'A/A*':>6} {'alpha_out':>10} {'eta_chem':>9} "
        f"{'eta_jet lo-hi':>16} {'vs 0.775':>10}"
    )
    for leg in ("head-on", "overtake"):
        for ar in (FLOWN_AREA_RATIO, WINDOW_AREA_RATIO):
            for w in (45.58, 75.0):
                r = jet_efficiency(leg=leg, closing_speed=w, area_ratio=ar)
                rows.append(r)
                verdict = "clears" if r.clears_target else "short"
                print(
                    f"{r.leg:>10} {w:7.2f} {ar:6.1f} {r.alpha_out:10.4f} {r.eta_chem:9.3f} "
                    f"{r.eta_jet_lo:7.3f}-{r.eta_jet_hi:<8.3f} {verdict:>10}"
                )

    print("\n== the ask's calibration: % of the reflection baseline ==")
    base_pancake = eta_geom_spread(ALPHA_MEASURED)
    base_isotropic = eta_geom_spread(ALPHA_ISOTROPIC)
    for ar in (FLOWN_AREA_RATIO, WINDOW_AREA_RATIO):
        hi = eta_geom_directed(alpha_after_expansion(ALPHA_MEASURED, ar))
        print(
            f"  A/A* = {ar:4.1f}: eta_geom {hi:.3f} is {100 * hi / base_pancake:5.0f}% of the "
            f"measured baseline ({base_pancake:.3f})"
        )
        print(
            f"{'':16} and {100 * hi / base_isotropic:5.0f}% of the isotropic one "
            f"({base_isotropic:.3f})"
        )
    print("  The ask: >=130% supports the paper, <100% is a serious problem.")

    print("\n== is mu actually invariant here? ==")
    print(f"{'station':>12} {'B [T]':>7} {'L [m]':>7} {'r_gyro/L':>10}")
    for label, b, length in (("throat", 20.0, 3.0), ("mid-column", 9.0, 6.0), ("exit", 5.0, 23.0)):
        print(f"{label:>12} {b:7.1f} {length:7.1f} {adiabaticity_parameter(b, length):10.2e}")
    print("  Far below 1 everywhere, so the guiding-centre picture holds comfortably.")

    write_jet(rows, args.output)
    print(f"\nwrote {len(rows)} rows -> {args.output}")


if __name__ == "__main__":
    main()
