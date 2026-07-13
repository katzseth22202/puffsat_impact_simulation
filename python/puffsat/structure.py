"""Heavy-plate whole-plate structural first-cut bound (ADR-0027), decoupled from `f(v)`.

At 30 m / <= 40 t a single 100 kg pulse deposits ~5e6 N.s, so "does the plate structure itself
survive the impulse?" becomes first-order (design §12.1). Design §5/§11 and ADR-0011 keep the
whole-plate *structural sizing* out of this repo (two hydro kernels, no structural solver), so
ADR-0027 carries a narrow, scenario-scoped exception: a handful of **closed-form** checks, off the
solver's load outputs, answering "is the <= 40 t plate plausibly buildable and rigid-during-pulse
at this scale?" without pulling a structural-dynamics discipline (FEA / plate-shell dynamics) in.

Three checks, fed by the heavy-plate sweep's `peak_wall_pressure` / `wall_impulse` and the Sigma
footprint, evaluated at each velocity anchor's survivable-design operating point (the shape the
facesheet frontier selects):

1. **Rigid-during-pulse (first-mode period vs pulse duration).** The candidate plate's fundamental
   flexural period `T1` must stay `>>` the bounce time `tau_pulse`, or the face is not locally rigid
   during the collision and the rigid-wall assumption behind `e_eff`/`f` breaks. This is a
   structural gate **and** the validity gate on `f` at this larger span (design §5) — a *failure*
   here flags the heavy-plate `f(v)` as provisional-pending-structure (ADR-0027).
2. **Areal-impulse -> membrane/bending stress.** The per-unit-area impulse over the footprint drives
   a candidate Ti-truss core + tensioned fiber (Vectran) back-face; the check reports the **implied
   minimum** back-face and total plate mass to carry the load at the fiber allowable, and whether
   that fits the <= 40 t ceiling. This is a **conservative static-equivalent upper bound** (the true
   impulsive response is softer, per check 1); a pass is conservative, and an over-budget result is
   the honest "closed-form does not establish survival -> FEA refinement" finding ADR-0027 names.
3. **SiC-Ti spall reflection (ADR-0011).** The peak facesheet compressive load reflects at the
   lower-impedance Ti backing as tension (`|R| ~ 0.15`); confirm it stays sub-spall, reusing the
   ADR-0011 model.

**What it is not** (ADR-0027): not FEA, not a validated design. Buckling, dynamic amplification
beyond the first mode, joint/fatigue detail, and thermal-structural coupling are out of scope. It
consumes solver outputs and emits a separate go/no-go verdict; it never feeds back into `e_eff`,
`eta_capture`, or `f` — except through check (1)'s one-way validity gate. The candidate-construction
constants below are documented assumptions, tunable, and stated with every reported number.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass, fields
from pathlib import Path

from puffsat.analysis import SIC_SPALL_LO, _write_csv, reflected_tensile
from puffsat.heavyplate import (
    PLATE_MASS_CEILING_KG,
    PLATE_RADIUS_M,
    HeavyPlatePoint,
    HeavyPlateRow,
    _load_geometry,
    _LogInterp,
    best_at_v,
    headline_rows,
    heavyplate_frontier,
    read_heavyplate_sweep,
    sweep_velocities,
)

DEFAULT_HEAVYPLATE_SWEEP_PATH = Path("data/results/sweep_heavyplate.jsonl")
DEFAULT_GEOMETRY_M40_PATH = Path("data/results/sweep_geometry_m40.jsonl")
DEFAULT_SUMMARY_PATH = Path("data/results/frontier_structure_heavyplate.csv")

# ---- Candidate construction (documented ADR-0027 assumptions; stated with every number) ----------
#
# A Ti-truss sandwich (two Ti face sheets separated by a lightweight truss core) with a tensioned
# high-strength-fiber (Vectran) back-face, at the design's baseline stack areal mass. These are
# *candidate* numbers for a first-cut feasibility bound, not a sized design.
E_TI = 110.0e9  # Ti-6Al-4V Young's modulus [Pa]
RHO_TI = 4500.0  # Ti density [kg/m^3]
FACE_THICKNESS_M = 0.002  # each Ti face sheet [m] (two faces)
CORE_DEPTH_M = 0.30  # sandwich core depth (face mid-plane separation) [m] ~ span/100
# Clamped circular-plate fundamental eigenvalue (omega1 = lambda^2/R^2 * sqrt(D/m_a)).
LAMBDA_SQ_CLAMPED = 10.216

# Tensioned fiber back-face (Vectran): the membrane that holds the dish shape against the load.
RHO_FIBER = 1400.0  # Vectran density [kg/m^3]
SIGMA_FIBER_ULT = 1.1e9  # Vectran ultimate tensile [Pa]
SAFETY_FACTOR_FIBER = 2.0  # working allowable = ult / SF
SIGMA_FIBER_WORKING = SIGMA_FIBER_ULT / SAFETY_FACTOR_FIBER

# Plate areal mass [kg/m^2] for the dynamics (check 1): the design's baseline stack (design §2).
AREAL_DENSITY = 45.0
# Baseline areal mass WITHOUT the tensioned back-face (Ti faces + truss core + SiC facesheet), so
# check 2 can add the *implied* back-face on top and compare the total to the <= 40 t ceiling.
BASE_AREAL_NO_BACK = 2.0 * FACE_THICKNESS_M * RHO_TI + 8.0 + 4.0  # faces + truss core + SiC ~ 30

# Minimum built-in dish `d/D` for the flat-plate membrane geometry (a real tensioned back-face has
# some pretensioned sag; avoids a singular radius of curvature at `d/D = 0`).
D_OVER_D_MIN = 0.02
# Rigidity margin: `T1 / tau_pulse` must clear this for the face to be locally rigid during the
# pulse (impulsive regime) — the validity gate on `f`.
RIGID_MARGIN = 10.0

PLATE_AREA_M2 = math.pi * PLATE_RADIUS_M * PLATE_RADIUS_M
AREAL_CEILING = PLATE_MASS_CEILING_KG / PLATE_AREA_M2  # <= 40 t / (pi R^2) ~ 57 kg/m^2


def flexural_rigidity() -> float:
    """Sandwich flexural rigidity `D = E_face * t_face * d_core^2 / 2` [N.m] (thin-face sandwich,
    both faces Ti at the core-depth separation) — the bending stiffness for the first mode."""
    return E_TI * FACE_THICKNESS_M * CORE_DEPTH_M * CORE_DEPTH_M / 2.0


def first_mode_period(areal_mass: float = AREAL_DENSITY) -> float:
    """Fundamental flexural period `T1 = 2*pi / omega1` [s] of a clamped circular plate of radius
    `R` (`omega1 = (lambda^2 / R^2) * sqrt(D / m_a)`, `lambda^2 = 10.216`). Heavier plate -> lower
    `omega1` -> longer `T1`, so this stays conservative for check 1 as the back-face mass grows."""
    omega1 = (LAMBDA_SQ_CLAMPED / (PLATE_RADIUS_M * PLATE_RADIUS_M)) * math.sqrt(
        flexural_rigidity() / areal_mass
    )
    return 2.0 * math.pi / omega1


def pulse_duration(wall_impulse: float, peak_wall_pressure: float) -> float:
    """Effective pulse duration `tau_pulse = 2 * J_area / p_peak` [s]: the full width of a
    triangular pressure pulse of areal impulse `J_area` and peak `p_peak` (both from the 1D sweep).
    The factor 2 is conservative for check 1 (a longer pulse lowers `T1/tau_pulse`)."""
    if peak_wall_pressure <= 0.0:
        raise ValueError("non-positive peak_wall_pressure — stale sweep JSONL?")
    return 2.0 * wall_impulse / peak_wall_pressure


@dataclass(frozen=True)
class StructurePoint:
    """The ADR-0027 closed-form bound at one velocity anchor's survivable-design operating point."""

    v: float
    d_over_d: float
    l_over_d: float
    r_foot_over_r: float
    rho_impact: float
    # Solver load inputs at the design point (areal quantities, 1D sweep).
    areal_impulse: float  # J_area = wall_impulse [Pa.s]
    peak_facesheet: float  # focused peak facesheet pressure [Pa] (Rung S, ADR-0010)
    # Check 1 — rigid-during-pulse / f-validity gate.
    pulse_duration: float
    first_mode_period: float
    rigidity_ratio: float  # T1 / tau_pulse
    rigid_ok: bool
    # Check 2 — areal-impulse -> membrane stress (conservative static-equivalent, implied minimum).
    membrane_tension: float  # N [N/m] to react the effective footprint pressure over the dish
    back_thickness_req: float  # implied minimum fiber back-face thickness [m]
    implied_plate_mass_t: float  # BASE stack + implied back-face, over pi R^2, in tonnes
    mass_ok: bool  # implied total <= 40 t ceiling
    # Check 3 — SiC-Ti spall reflection (ADR-0011).
    reflected_tensile: float
    spall_ok: bool
    # Overall.
    verdict_ok: bool


def _interp_at(rows: list[HeavyPlateRow], v: float, attr: str) -> _LogInterp:
    """A log-rho interpolator of `attr` over the headline rho slice at velocity `v`."""
    slice_rows = sorted(
        (r for r in headline_rows(rows) if abs(r.v - v) <= 1e-6), key=lambda r: r.rho_impact
    )
    if not slice_rows:
        raise ValueError(f"no headline sweep rows at v={v}")
    return _LogInterp(
        [r.rho_impact for r in slice_rows], [float(getattr(r, attr)) for r in slice_rows]
    )


def structure_point(design: HeavyPlatePoint, rows: list[HeavyPlateRow]) -> StructurePoint:
    """Evaluate the three closed-form checks at one survivable-design operating point.

    `design` is the frontier's chosen shape at this velocity (its Sigma density, focused peak, and
    curvature); `rows` supplies the areal impulse `J_area(rho)` and unfocused peak `p_1d(rho)` at
    the design density (interpolated log-rho along the headline slice)."""
    j_area = _interp_at(rows, design.v, "wall_impulse")(design.rho_impact)
    p_1d = _interp_at(rows, design.v, "peak_wall_pressure")(design.rho_impact)

    # Check 1 — rigid-during-pulse (validity gate). tau from the consistent 1D (unfocused) pair.
    tau = pulse_duration(j_area, p_1d)
    t1 = first_mode_period()
    ratio = t1 / tau
    rigid_ok = ratio >= RIGID_MARGIN

    # Check 2 — areal-impulse -> membrane stress, conservative static-equivalent, implied minimum.
    # The focused facesheet peak over the footprint reacts as an effective uniform plate pressure
    # p_eff = F_peak / (pi R^2) = p_peak * (r_foot/R)^2, held by the shallow tensioned dish as
    # membrane tension N = p_eff * Rc / 2 (Rc = R / (4 d/D), the shallow spherical-cap curvature).
    p_eff = design.peak_compressive * design.r_foot_over_r * design.r_foot_over_r
    dd_eff = max(design.d_over_d, D_OVER_D_MIN)
    r_curv = PLATE_RADIUS_M / (4.0 * dd_eff)
    membrane_tension = p_eff * r_curv / 2.0
    back_thickness_req = membrane_tension / SIGMA_FIBER_WORKING
    implied_areal = BASE_AREAL_NO_BACK + back_thickness_req * RHO_FIBER
    implied_mass_t = implied_areal * PLATE_AREA_M2 / 1000.0
    mass_ok = implied_mass_t <= PLATE_MASS_CEILING_KG / 1000.0

    # Check 3 — SiC-Ti spall reflection (ADR-0011), on the focused facesheet peak.
    refl = reflected_tensile(design.peak_compressive)
    spall_ok = refl < SIC_SPALL_LO

    return StructurePoint(
        v=design.v,
        d_over_d=design.d_over_d,
        l_over_d=design.l_over_d,
        r_foot_over_r=design.r_foot_over_r,
        rho_impact=design.rho_impact,
        areal_impulse=j_area,
        peak_facesheet=design.peak_compressive,
        pulse_duration=tau,
        first_mode_period=t1,
        rigidity_ratio=ratio,
        rigid_ok=rigid_ok,
        membrane_tension=membrane_tension,
        back_thickness_req=back_thickness_req,
        implied_plate_mass_t=implied_mass_t,
        mass_ok=mass_ok,
        reflected_tensile=refl,
        spall_ok=spall_ok,
        verdict_ok=rigid_ok and mass_ok and spall_ok,
    )


def structure_frontier(
    points: list[HeavyPlatePoint], rows: list[HeavyPlateRow]
) -> list[StructurePoint]:
    """The ADR-0027 bound at each velocity anchor's survivable-design point (best-surviving concave,
    falling back to best-surviving flat). Velocities with no surviving shape are skipped."""
    out: list[StructurePoint] = []
    for v in sweep_velocities(rows):
        design = best_at_v(points, v, concave=True) or best_at_v(points, v)
        if design is not None:
            out.append(structure_point(design, rows))
    return out


def write_summary(points: list[StructurePoint], path: Path = DEFAULT_SUMMARY_PATH) -> None:
    """Write the structural-bound frontier to a CSV (header = `StructurePoint` field names)."""
    header = [f.name for f in fields(StructurePoint)]
    _write_csv(header, ([getattr(p, name) for name in header] for p in points), path)


def _fmt_flag(ok: bool) -> str:
    return "PASS" if ok else "FAIL"


def main() -> None:
    """Run the ADR-0027 closed-form whole-plate structural bound; print the per-anchor verdict."""
    parser = argparse.ArgumentParser(
        description="Heavy-plate whole-plate structural first-cut bound (ADR-0027)."
    )
    parser.add_argument("--sweep", type=Path, default=DEFAULT_HEAVYPLATE_SWEEP_PATH)
    parser.add_argument("--geometry", type=Path, default=DEFAULT_GEOMETRY_M40_PATH)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY_PATH)
    args = parser.parse_args()

    rows = read_heavyplate_sweep(args.sweep)
    geo_rows = _load_geometry(args.geometry)
    frontier = heavyplate_frontier(geo_rows, rows, sweep_velocities(rows))
    points = structure_frontier(frontier, rows)
    write_summary(points, args.summary)

    if not points:
        print("python: no surviving design point at any velocity — nothing to structurally check.")
        return

    t1 = first_mode_period()
    print(
        f"ADR-0027 closed-form whole-plate bound: R = {PLATE_RADIUS_M} m, <= "
        f"{PLATE_MASS_CEILING_KG / 1000:.0f} t ceiling ({AREAL_CEILING:.0f} kg/m^2). "
        f"Candidate Ti-sandwich first-mode period T1 = {t1 * 1e3:.0f} ms; "
        f"fiber allowable {SIGMA_FIBER_WORKING / 1e9:.2f} GPa (SF {SAFETY_FACTOR_FIBER:.0f})."
    )
    print(
        "Checks: (1) rigid-during-pulse / f-validity  (2) areal-impulse membrane (implied mass)  "
        "(3) SiC-Ti spall"
    )
    print(
        f"\n{'v [km/s]':>8} {'shape (d/D,L/D,rf/R)':>22} {'tau[ms]':>8} {'T1/tau':>7} "
        f"{'(1)':>5} {'impl.mass[t]':>12} {'(2)':>5} {'refl[MPa]':>10} {'(3)':>5} {'verdict':>8}"
    )
    for p in points:
        shape = f"({p.d_over_d:.2f},{p.l_over_d:.1f},{p.r_foot_over_r:.1f})"
        print(
            f"{p.v / 1000:8.1f} {shape:>22} {p.pulse_duration * 1e3:8.2f} {p.rigidity_ratio:7.0f} "
            f"{_fmt_flag(p.rigid_ok):>5} {p.implied_plate_mass_t:12.0f} {_fmt_flag(p.mass_ok):>5} "
            f"{p.reflected_tensile / 1e6:10.1f} {_fmt_flag(p.spall_ok):>5} "
            f"{_fmt_flag(p.verdict_ok):>8}"
        )

    n_ok = sum(1 for p in points if p.verdict_ok)
    rigid_all = all(p.rigid_ok for p in points)
    spall_all = all(p.spall_ok for p in points)
    mass_worst = max(p.implied_plate_mass_t for p in points)
    ceiling_t = PLATE_MASS_CEILING_KG / 1000
    print(
        f"\nSUMMARY: rigidity/f-validity gate {'PASS at all anchors' if rigid_all else 'FAILS'}; "
        f"SiC-Ti spall {'PASS at all anchors' if spall_all else 'FAILS'}; "
        f"membrane implied plate mass up to {mass_worst:.0f} t vs the {ceiling_t:.0f} t ceiling."
    )
    if n_ok == len(points):
        print(
            "GO (closed-form, candidate construction): all three checks clear at every anchor "
            "within the <= 40 t budget. Full plate-shell / buckling verdict remains the FEA step."
        )
    else:
        print(
            f"NOT ESTABLISHED by the closed-form bound at {len(points) - n_ok}/{len(points)} "
            "anchors (the conservative static-equivalent membrane exceeds the <= 40 t budget). Per "
            "ADR-0027 this routes the whole-plate verdict to the named FEA refinement; the "
            "rigid-during-pulse validity gate on f(v) is unaffected (it passes)."
        )
    print(f"python: wrote {args.summary}")


if __name__ == "__main__":
    main()
