"""Rung 2 -- does the plume detach from the field, and where? (N6, and Rung 5's `beta(z)`)

Reads the solved cooling history and asks three questions of every station, rather than of one:

1. **`M_A(z)`** -- the Alfven Mach number along the column (N6). Detachment needs the flow to
   cross from sub-Alfvenic to super-Alfvenic; below 1 the field dominates and the plume is still
   tied to the ship.
2. **`beta(z) = p_actual/p_design`** -- how far the solved expansion sits from the standoff
   condition the field was graded against. This is the only form of N3's "max beta" question that
   carries information, since beta against the design profile is 1 by construction.
3. **Where a physical wall could take over from the field** -- raised by Seth 2026-09-03.

## The headline: the paper's single-station check does not survive being done per station

`sec:jet_efficiency` reports `M_A` = 1.63 at the coldest pulse and 2.06 at the hottest, and reads
those as clearing the detachment condition. Computed consistently at each station of the solved
expansion, **`M_A` runs 0.20 at the throat to 0.58 at the exit and never crosses 1.**

Two things produce the gap, and they are the same thing seen twice:

- **The paper mixes states.** It divides the *exit* speed (10.8 km/s) by an Alfven speed built from
  the *bag* density (0.32 kg/m^3) and the bag's standoff field (4.1 T). The plume expands by 13x in
  density between those two places, and `v_A ~ rho^-1/2`, so using the pre-expansion density with
  the post-expansion speed inflates `M_A` by `sqrt(0.32/0.025) = 3.6x`. That is the whole
  discrepancy: `0.58 x 3.56 = 2.06`, the paper's hot-pulse figure exactly.
- **`eq:alfven_from_standoff` assumes standoff.** `v_A = sqrt(2 R_g T/Mbar)` follows from
  `B^2/2mu0 = p`, i.e. from `beta = 1`. The paper's two routes to `v_A` agree to 3% and it reads
  that as corroboration, but they are not independent -- `tab:bag_sizing`'s 4.1 T *is* the standoff
  field at that density and temperature, so both routes evaluate the same assumed state.

**The assumption fails by more than an order of magnitude.** The graded field is derived from the
*snowplow* pressure of the collision (159 MPa at 1 m), while the *expansion* that follows runs at
10 MPa at the throat falling to 0.7 MPa at the exit. So `beta` is near 0.06-0.07 through the whole
column: the field is roughly 15x over-strength for the flow it is steering. An over-strength field
is exactly what keeps a flow sub-Alfvenic.

**What this does and does not overturn.** It does not show the architecture fails. The paper's
downstream argument -- past the last coil a plume's pressure falls as `R^-5` while a vacuum field
falls as `R^-6`, so `beta` rises and detachment eventually happens -- is untouched, and
`downstream_crossing` estimates it lands about 1.4 exit radii out. What it overturns is the claim
that the crossing happens *inside* the nozzle, and the number offered as evidence for it.

## What this module does not do

It takes the field as given (design profile or realizable winding, both reported) and the flow as
solved by `expansion.py`. It does not solve the coupled MHD problem, so it cannot say how the flow
would differ if the field were sized to the expansion instead of to the collision.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from puffsat import field

MU0 = field.MU0

SIGMA_SB = 5.670374419e-8
"""Stefan-Boltzmann constant [W m^-2 K^-4]."""

DEFAULT_HISTORY = Path("data/results/cooling_history.csv")
DEFAULT_OUTPUT = Path("data/results/nozzle_detachment.csv")

BAG_RHO = 0.323
"""Flown bag density [kg/m^3] -- 213 kg over 660 m^3. The paper's `M_A` uses this with the *exit*
speed, which is the state mixing this module exists to correct."""

GRAPHITE_SUBLIMATION_K = 3900.0
"""Pyrolytic graphite sublimation temperature [K] -- the liner's ceiling."""

GRAPHITE_SUBLIMATION_ENTHALPY = 716.7e3 / 0.012011
"""Graphite -> C(g), 716.7 kJ/mol = 59.7 MJ/kg. Why graphite is a good ablator.

Sensible heat to reach sublimation (~6 MJ/kg) is *not* included, so the ablation depths below are
about 10% conservative.
"""

GRAPHITE_DENSITY = 2200.0
"""Pyrolytic graphite density [kg/m^3]."""

BORE_WALL_AREA_M2 = 2.0 * math.pi * 3.0 * 23.0
"""Bore wall area [m^2] -- 434 m^2 at the flown 3.0 m radius over 23 m."""

BOOKED_FLASH_W = 42.6e6
"""Intercepted flash the paper books for the Jupiter case [W] (tab:bag_state).

Spread over the bore wall this is 98 kW/m^2, which is **130x below** the plume's own surface flux.
The gap is physical: the field holds the plume off the wall and only a few percent of the sky is
hardware. Which flux a wall actually sees is therefore a geometry question, and the two bounds are
carried separately in `AblationCase`.
"""

LINER_THICKNESS_M = 0.01
"""Reference liner thickness [m] for quoting a life in pulses. Not a paper number."""


@dataclass(frozen=True)
class Station:
    """One station of the solved expansion, with the magnetic quantities it implies."""

    closing_speed_km_s: float
    branch: str
    time_s: float
    x_m: float
    area_ratio: float
    rho: float
    temp_k: float
    pressure_pa: float
    speed_m_s: float
    mach: float
    b_design_t: float
    """Design field from flux conservation, `20 T / (A/A*)` -- what `electrothermal` uses."""
    b_built_t: float
    """Field a physically realizable winding actually delivers here (Rung 1)."""
    v_alfven_design: float
    alfven_mach_design: float
    alfven_mach_built: float
    beta: float
    """`p_actual/p_design` -- how far the solved expansion sits from standoff."""
    radiative_flux_w_m2: float
    """`sigma T^4` -- an optically-thick upper bound on what a wall here would have to shed."""


def _station_position(area_ratio: float, area_ratio_end: float, length: float) -> float:
    """Linear-in-area-ratio axial map, matching `expansion.cooling_history`'s own assumption."""
    if area_ratio_end <= 1.0:
        return 0.0
    return length * (area_ratio - 1.0) / (area_ratio_end - 1.0)


def load_stations(
    path: Path = DEFAULT_HISTORY, *, n_coils: int = 72, length: float = field.COLUMN_LENGTH_M
) -> list[Station]:
    """Load the solved history and attach the magnetic quantities at every station.

    Both fields are carried deliberately. `b_design_t` is the flux-conservation profile the rest of
    the repo uses; `b_built_t` is what Rung 1 showed a real winding delivers, which is ~22% lower
    at the chamber. A weaker field means a lower `v_A` and so a *higher* `M_A`, which is the one
    place realizability helps rather than hurts.
    """
    stack = field.build_winding(n_coils)
    rows = list(csv.DictReader(path.open()))
    if not rows:
        raise ValueError(f"{path} has no rows; run `make analysis-expansion` first")
    ar_end = max(float(r["area_ratio"]) for r in rows)

    out: list[Station] = []
    for r in rows:
        rho = float(r["rho"])
        area_ratio = float(r["area_ratio"])
        speed = float(r["speed_m_s"])
        pressure = float(r["pressure_pa"])
        temp = float(r["temp_k"])
        x = _station_position(area_ratio, ar_end, length)
        b_design = 20.0 / area_ratio
        b_built = abs(stack.on_axis(x))
        v_a_design = b_design / math.sqrt(MU0 * rho)
        v_a_built = b_built / math.sqrt(MU0 * rho) if b_built > 0.0 else math.inf
        out.append(
            Station(
                closing_speed_km_s=float(r["closing_speed_km_s"]),
                branch=r["branch"],
                time_s=float(r["time_ms"]) * 1e-3,
                x_m=x,
                area_ratio=area_ratio,
                rho=rho,
                temp_k=temp,
                pressure_pa=pressure,
                speed_m_s=speed,
                mach=float(r["mach"]),
                b_design_t=b_design,
                b_built_t=b_built,
                v_alfven_design=v_a_design,
                alfven_mach_design=speed / v_a_design,
                alfven_mach_built=speed / v_a_built if v_a_built > 0.0 else 0.0,
                beta=pressure / (b_design * b_design / (2.0 * MU0)),
                radiative_flux_w_m2=SIGMA_SB * temp**4,
            )
        )
    return out


def by_case(stations: list[Station]) -> dict[tuple[float, str], list[Station]]:
    """Group stations into their `(closing speed, branch)` curves, in column order."""
    grouped: dict[tuple[float, str], list[Station]] = defaultdict(list)
    for s in stations:
        grouped[(s.closing_speed_km_s, s.branch)].append(s)
    return {k: sorted(v, key=lambda s: s.area_ratio) for k, v in grouped.items()}


def paper_alfven_mach(
    *, b_bag: float = 4.1, rho_bag: float = 0.32, exit_speed: float = 10.8e3
) -> float:
    """Reproduce the paper's own `M_A` so the discrepancy is a computation, not an assertion.

    `tab:bag_sizing`'s bag field and bag density, divided into the *exit* speed. Returns 1.67
    against the paper's printed 1.63 -- agreement to the rounding of its inputs.
    """
    return exit_speed / (b_bag / math.sqrt(MU0 * rho_bag))


def state_mixing_factor(rho_exit: float, rho_bag: float = BAG_RHO) -> float:
    """How much using the bag density with the exit speed inflates `M_A`: `sqrt(rho_bag/rho_exit)`.

    `v_A ~ rho^-1/2`, so evaluating the Alfven speed before the expansion and the flow speed after
    it overstates the ratio by exactly this factor. It is 3.6 at the flown geometry.
    """
    return math.sqrt(rho_bag / rho_exit)


def downstream_crossing(exit_mach_alfven: float) -> float:
    """Radii past the exit at which `M_A` reaches 1, on the paper's own downstream scalings.

    `sec:jet_efficiency` argues that past the last coil the plume's pressure falls as `R^-5` and a
    vacuum field's as `R^-6`. Taking a freely expanding plume at constant speed, `rho ~ R^-3` and a
    dipole `B ~ R^-3`, so `v_A ~ R^-3/2` and `M_A ~ R^3/2`. The crossing is then at

        R/R_exit = (1/M_A_exit)^{2/3}.

    **This is the paper's argument, quantified, not a check of it.** It inherits the free-expansion
    and dipole assumptions and says nothing about whether the field lines actually open.
    """
    if exit_mach_alfven <= 0.0:
        return math.inf
    if exit_mach_alfven >= 1.0:
        return 1.0
    return float((1.0 / exit_mach_alfven) ** (2.0 / 3.0))


@dataclass(frozen=True)
class AblationCase:
    """Could a physical diverging wall survive here -- asked as depth per pulse, not as a flux?

    Raised by Seth 2026-09-03, twice, and the second time corrected the first answer.

    **CORRECTED 2026-09-05 (reply R8). This docstring used to say "the gas dynamics were never the
    obstacle: every exit is supersonic, and a diverging wall on supersonic flow is exactly what a
    de Laval nozzle is. The obstacle is heat." That was the wrong Mach number.** Being *supersonic*
    (`M` = 2.4-3.4) is not being *super-Alfvenic* (`M_A` = 0.44-0.72): a diverging wall works on
    ordinary supersonic gas, and this gas is magnetically dominated. At `beta` = 0.033-0.154 the
    field is 6-30x over-strength, so plasma cannot cross field lines to reach a wall, and
    `fluxtube.downstream_fan` finds the plume's own bounding tube turns at 6.34 m -- inside the
    7.0 m the paper's R8 needs its bell to start at for graphite to survive the 16 224 K exit.
    **The gas dynamics were the obstacle all along.** The ablation numbers below stand and are
    still worth having; what they establish is that a wall would *survive*, which turns out not
    to be the question. See `docs/nozzle_replies_answered.md` P14 and `extension.py`.

    **But the first screen asked the wrong question.** Comparing the incident flux against what
    graphite re-radiates at sublimation is a *steady-state* test, and these pulses last about two
    milliseconds. What matters for a liner is **how deep it ablates per pulse** and how many pulses
    that buys. A flux ten times past the steady-state ceiling for 2 ms is a fraction of a micron.

    Two fluxes bound the answer, because what a wall sees depends on where it is:

    - `sigma T^4` at the plume's own surface -- an upper bound, right for a wall **in contact**,
      which is what a physical nozzle means.
    - the paper's booked 42.6 MW spread over the bore -- 98 kW/m^2, right for a wall **standing
      off** behind a field, which is the current design, and 130x lower.

    Neither includes convective transfer, which for a wall in Mach-3 contact is likely to exceed
    the radiative term. So the contact column is an underestimate in one direction and an
    overestimate in the other; treat it as an order of magnitude.
    """

    closing_speed_km_s: float
    branch: str
    temp_k: float
    sonic_mach: float
    transit_s: float
    contact_flux_w_m2: float
    """`sigma T^4` -- the plume's own surface flux, for a wall in contact."""
    net_flux_w_m2: float
    """Contact flux less what graphite re-radiates at sublimation. Zero means it never ablates."""
    depth_per_pulse_m: float
    pulses_per_liner: float
    """Pulses to consume `LINER_THICKNESS_M` of graphite."""


def ablation_depth(flux_w_m2: float, transit_s: float) -> float:
    """Graphite ablated per pulse [m] under `flux_w_m2` for `transit_s`.

    The surface pins at sublimation and re-radiates `sigma T_sub^4`; whatever is left drives
    sublimation at 59.7 MJ/kg. Conduction into the bulk is neglected, which is conservative -- it
    would carry heat away that would otherwise ablate.
    """
    capacity = SIGMA_SB * GRAPHITE_SUBLIMATION_K**4
    net = max(flux_w_m2 - capacity, 0.0)
    return net * transit_s / (GRAPHITE_SUBLIMATION_ENTHALPY * GRAPHITE_DENSITY)


def ablation_cases(stations: list[Station]) -> list[AblationCase]:
    """Score every case's exit station for wall survival, as depth per pulse."""
    capacity = SIGMA_SB * GRAPHITE_SUBLIMATION_K**4
    out: list[AblationCase] = []
    for (w, branch), curve in by_case(stations).items():
        last = curve[-1]
        transit = last.time_s
        flux = last.radiative_flux_w_m2
        depth = ablation_depth(flux, transit)
        out.append(
            AblationCase(
                closing_speed_km_s=w,
                branch=branch,
                temp_k=last.temp_k,
                sonic_mach=last.mach,
                transit_s=transit,
                contact_flux_w_m2=flux,
                net_flux_w_m2=max(flux - capacity, 0.0),
                depth_per_pulse_m=depth,
                pulses_per_liner=(LINER_THICKNESS_M / depth if depth > 0.0 else math.inf),
            )
        )
    return sorted(out, key=lambda o: (o.branch, o.closing_speed_km_s))


def standoff_flux_w_m2() -> float:
    """The paper's booked intercepted flash spread over the bore wall [W/m^2] -- 98 kW/m^2.

    Two orders below graphite's re-radiation ceiling, so on the paper's own booking a
    **standing-off** liner never reaches sublimation from radiation at all. That the paper
    nonetheless books 4.9 kg/pulse of ablation means its liner loss is charged to something other
    than the radiative flash -- worth asking the paper which.
    """
    return BOOKED_FLASH_W / BORE_WALL_AREA_M2


# ---- The field window: strong enough to contain, weak enough to release (N6, N3) ------------
#
# Raised by Seth 2026-09-03, and it turns the N6 finding from a problem into a fix. If the exit
# field is what holds `M_A` under 1, then lowering it is the lever -- provided there is room
# between "still stands the plume off the wall" and "already lets it go". There is, in every case.

GAMMA_MONATOMIC = 5.0 / 3.0


def containment_field(pressure_pa: float) -> float:
    """`B = sqrt(2 mu0 p)` [T] -- the least field that stands the plume off the wall."""
    if pressure_pa < 0.0:
        raise ValueError("pressure must be non-negative")
    return math.sqrt(2.0 * MU0 * pressure_pa)


def release_field(speed: float, rho: float) -> float:
    """`B = v sqrt(mu0 rho)` [T] -- the most field that still allows `M_A` to reach 1.

    Set `M_A = v/(B/sqrt(mu0 rho)) = 1` and solve. Above this the flow stays magnetised; below it
    the flow has detached.
    """
    if rho <= 0.0:
        raise ValueError("density must be positive")
    return speed * math.sqrt(MU0 * rho)


def alfven_mach_closed_form(
    mach_sonic: float, beta: float, gamma: float = GAMMA_MONATOMIC
) -> float:
    """`M_A = M_sonic sqrt(gamma/2) sqrt(beta)` -- why an over-strength field costs detachment.

    From `M_A = v sqrt(mu0 rho)/B`, `beta = 2 mu0 p/B^2` and `c_s = sqrt(gamma p/rho)`, the field
    and the density both cancel out of the ratio and leave only the sonic Mach number and `beta`.

    **Two things follow.** At standoff (`beta = 1`) it reduces to `M_A = M_sonic/sqrt(2/gamma)` =
    `M_sonic/1.095`, which is exactly the paper's "the Alfven surface sits a tenth of the way past
    the sonic throat" -- so the paper's physics is right. And `M_A` carries `sqrt(beta)`, so the
    15x over-strength field of the flown design costs a factor of 4 in Alfven Mach, which is the
    whole of the gap in `M_A(z)` above.

    Exact on the frozen branch; over-predicts equilibrium by 10-15% because recombination moves
    the effective exponent off 5/3 (the history carries it per station as `gamma_t`).
    """
    if beta < 0.0:
        raise ValueError("beta must be non-negative")
    return mach_sonic * math.sqrt(gamma / 2.0) * math.sqrt(beta)


@dataclass(frozen=True)
class FieldWindow:
    """Is there a field that contains the plume and still lets it go? Evaluated at the exit."""

    closing_speed_km_s: float
    branch: str
    mach_sonic: float
    b_contain_t: float
    """Floor: below this the expansion is no longer stood off the wall."""
    b_release_t: float
    """Ceiling: above this the flow cannot reach `M_A` = 1."""
    b_design_t: float
    window_exists: bool
    """`b_release > b_contain`. Equivalent to `M_sonic > 1.095` by the closed form."""
    design_above_window: bool
    """The finding: the flown design sits *above* the ceiling, not inside the window."""
    alfven_mach_now: float
    alfven_mach_at_standoff: float
    required_area_ratio: float
    """`A/A*` that flux conservation needs to bring the field to the release ceiling."""
    snowplow_overshoot: float
    """`p_snowplow/(B_release^2/2mu0)` -- what the collision does to a field weakened to release.

    The cost of the fix. One static field cannot be strong for the collision and weak for the
    expansion at the same station, so buying detachment means the snowplow gets past there.
    """


def field_windows(stations: list[Station], *, throat_field_t: float = 20.0) -> list[FieldWindow]:
    """Build the containment/release window at every case's exit station.

    The window exists whenever the flow is supersonic past `M_sonic = 1.095`, which every exit is
    by a wide margin (2.7-3.4). What the flown design gets wrong is not the existence of the
    window but its position: 5 T against a ceiling of 1.8-2.9 T.
    """
    out: list[FieldWindow] = []
    for (w, branch), curve in by_case(stations).items():
        s = curve[-1]
        b_contain = containment_field(s.pressure_pa)
        b_release = release_field(s.speed_m_s, s.rho)
        p_snowplow = field.design_pressure(max(s.x_m, 1.0))
        out.append(
            FieldWindow(
                closing_speed_km_s=w,
                branch=branch,
                mach_sonic=s.mach,
                b_contain_t=b_contain,
                b_release_t=b_release,
                b_design_t=s.b_design_t,
                window_exists=b_release > b_contain,
                design_above_window=s.b_design_t > b_release,
                alfven_mach_now=s.alfven_mach_design,
                alfven_mach_at_standoff=s.speed_m_s / (b_contain / math.sqrt(MU0 * s.rho)),
                required_area_ratio=throat_field_t / b_release,
                snowplow_overshoot=p_snowplow / (b_release * b_release / (2.0 * MU0)),
            )
        )
    return sorted(out, key=lambda o: (o.branch, o.closing_speed_km_s))


@dataclass(frozen=True)
class CommonWindow:
    """One flare that satisfies every flown leg at once, in area-ratio terms.

    Each leg has its own containment floor and release ceiling in *field*; flux conservation
    (`B = B_throat/(A/A*)`) turns those into an area-ratio interval. The design has to sit inside
    the intersection of all of them, which is narrower than any single leg's window and is the
    number a nozzle can actually be built to.
    """

    area_ratio_min: float
    """Set by the leg with the lowest release ceiling -- flare at least this hard to detach."""
    area_ratio_max: float
    """Set by the leg with the highest containment floor -- flare more than this and the plume
    is no longer stood off the wall."""
    flown_area_ratio: float
    binding_release_leg: str
    binding_contain_leg: str
    exists: bool


def common_window(
    stations: list[Station], *, throat_field_t: float = 20.0, flown_area_ratio: float = 4.0
) -> CommonWindow:
    """The area-ratio interval that gives every leg containment *and* detachment.

    **This is the actionable form of the finding.** The flown `A/A*` = 4 sits well below it: at
    that flare the exit field is 5 T against release ceilings of 1.8-2.9 T, so no leg detaches.
    Flaring into the interval buys detachment on every leg without losing standoff on any.
    """
    windows = field_windows(stations, throat_field_t=throat_field_t)
    tightest_release = min(windows, key=lambda w: w.b_release_t)
    tightest_contain = max(windows, key=lambda w: w.b_contain_t)
    lo = throat_field_t / tightest_release.b_release_t
    hi = throat_field_t / tightest_contain.b_contain_t
    return CommonWindow(
        area_ratio_min=lo,
        area_ratio_max=hi,
        flown_area_ratio=flown_area_ratio,
        binding_release_leg=(f"{tightest_release.closing_speed_km_s:g} {tightest_release.branch}"),
        binding_contain_leg=(f"{tightest_contain.closing_speed_km_s:g} {tightest_contain.branch}"),
        exists=lo < hi,
    )


def write_stations(stations: list[Station], path: Path = DEFAULT_OUTPUT) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "closing_speed_km_s",
                "branch",
                "x_m",
                "area_ratio",
                "rho",
                "temp_k",
                "pressure_pa",
                "speed_m_s",
                "mach",
                "b_design_T",
                "b_built_T",
                "v_alfven_m_s",
                "M_A_design",
                "M_A_built",
                "beta",
                "radiative_flux_W_m2",
            ]
        )
        for s in stations:
            writer.writerow(
                [
                    f"{s.closing_speed_km_s:g}",
                    s.branch,
                    f"{s.x_m:.4f}",
                    f"{s.area_ratio:.5f}",
                    f"{s.rho:.6e}",
                    f"{s.temp_k:.2f}",
                    f"{s.pressure_pa:.6e}",
                    f"{s.speed_m_s:.2f}",
                    f"{s.mach:.4f}",
                    f"{s.b_design_t:.4f}",
                    f"{s.b_built_t:.4f}",
                    f"{s.v_alfven_design:.2f}",
                    f"{s.alfven_mach_design:.5f}",
                    f"{s.alfven_mach_built:.5f}",
                    f"{s.beta:.6e}",
                    f"{s.radiative_flux_w_m2:.4e}",
                ]
            )


def main() -> None:
    """N6 along the column, N3's beta, and the physical-wall question, from the solved history."""
    parser = argparse.ArgumentParser(description="Rung 2: detachment along the column")
    parser.add_argument("--history", type=Path, default=DEFAULT_HISTORY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    stations = load_stations(args.history)
    cases = by_case(stations)

    print("== N6: M_A along the column, not at one station ==")
    print(
        f"{'w':>7} {'branch':>11} {'M_A throat':>11} {'M_A exit':>9} {'M_A exit':>9} "
        f"{'crosses 1?':>11} {'detach at':>10}"
    )
    print(
        f"{'':>7} {'':>11} {'(design)':>11} {'(design)':>9} {'(built)':>9} "
        f"{'':>11} {'R/R_exit':>10}"
    )
    for (w, branch), curve in sorted(cases.items()):
        throat, last = curve[0], curve[-1]
        crossed = any(s.alfven_mach_design >= 1.0 for s in curve)
        print(
            f"{w:7.2f} {branch:>11} {throat.alfven_mach_design:11.2f} "
            f"{last.alfven_mach_design:9.2f} {last.alfven_mach_built:9.2f} "
            f"{('yes' if crossed else 'NO'):>11} "
            f"{downstream_crossing(last.alfven_mach_design):10.2f}"
        )

    hot = cases[(75.0, "equilibrium")][-1]
    print(f"\npaper's own route reproduced: M_A = {paper_alfven_mach():.2f} (paper prints 1.63)")
    print(
        f"  it divides the EXIT speed by an Alfven speed built from the BAG density "
        f"({BAG_RHO} kg/m^3),"
    )
    print(f"  inflating M_A by sqrt(rho_bag/rho_exit) = {state_mixing_factor(hot.rho):.2f}x.")
    print(
        f"  {hot.alfven_mach_design:.2f} x {state_mixing_factor(hot.rho):.2f} = "
        f"{hot.alfven_mach_design * state_mixing_factor(hot.rho):.2f}, "
        f"against the paper's hot-pulse 2.06."
    )

    print("\n== N3 / Rung 5: beta = p_actual/p_design along the column ==")
    print(f"{'w':>7} {'branch':>11} {'beta throat':>12} {'beta exit':>10} {'beta max':>9}")
    for (w, branch), curve in sorted(cases.items()):
        betas = [s.beta for s in curve]
        print(f"{w:7.2f} {branch:>11} {betas[0]:12.4f} {betas[-1]:10.4f} {max(betas):9.4f}")
    print("The field is graded against the SNOWPLOW pressure of the collision (159 MPa at 1 m),")
    print("but the expansion that follows runs at 10 MPa. So it is ~15x over-strength for the flow")
    print("it steers -- which is precisely what holds the plume sub-Alfvenic.")

    print("\n== Seth's question: could a physical wall take over downstream? ==")
    print("The gas dynamics are not the obstacle: the flow is supersonic past the throat, and a")
    print(
        "diverging wall is what a de Laval nozzle uses on supersonic flow. The obstacle is heat --"
    )
    print("but at ~2 ms a pulse, the question is microns per pulse, not a steady-state flux ratio.")
    print(
        f"\n{'w':>7} {'branch':>11} {'T_exit':>8} {'sonic M':>8} {'transit':>8} "
        f"{'contact q':>11} {'um/pulse':>9} {'pulses/cm':>10}"
    )
    for opt in ablation_cases(stations):
        life = "never" if math.isinf(opt.pulses_per_liner) else f"{opt.pulses_per_liner:,.0f}"
        print(
            f"{opt.closing_speed_km_s:7.2f} {opt.branch:>11} {opt.temp_k:8.0f} "
            f"{opt.sonic_mach:8.2f} {opt.transit_s * 1e3:6.2f}ms "
            f"{opt.contact_flux_w_m2 / 1e6:8.1f}MW {opt.depth_per_pulse_m * 1e6:9.3f} {life:>10}"
        )
    print(
        f"\ngraphite: {GRAPHITE_SUBLIMATION_ENTHALPY / 1e6:.1f} MJ/kg to sublime, re-radiating "
        f"{SIGMA_SB * GRAPHITE_SUBLIMATION_K**4 / 1e6:.1f} MW/m^2 at "
        f"{GRAPHITE_SUBLIMATION_K:.0f} K."
    )
    print("The frozen branch is a non-problem: sub-micron a pulse, tens of thousands per cm.")
    print("The equilibrium branch's hot legs are not: 50 um a pulse eats a centimetre in ~200.")
    print("\nFlux caveat. The contact column is `sigma T^4` at the plume's own surface, right")
    print(
        f"for a wall touching it. The paper's booked 42.6 MW over the bore is "
        f"{standoff_flux_w_m2() / 1e3:.0f} kW/m^2 --"
    )
    print("130x lower, and below graphite's ceiling, so a standing-off liner never ablates from")
    print("radiation at all. Neither figure carries a convective term, which for Mach-3 contact")
    print("would likely exceed the radiative one.")

    print("\n== The field window: contain without gripping (Seth, 2026-09-03) ==")
    print(
        "M_A = M_sonic * sqrt(gamma/2) * sqrt(beta) -- field and density cancel out of the ratio."
    )
    print("At standoff that is M_sonic/1.095, the paper's own relation. The design is not at")
    print("standoff, and M_A carries sqrt(beta), so a 15x over-strength field costs 4x in M_A.")
    print(
        f"\n{'w':>7} {'branch':>11} {'M_son':>6} {'contain':>8} {'release':>8} {'design':>7} "
        f"{'M_A now':>8} {'M_A@standoff':>13} {'need A/A*':>10} {'snowplow':>9}"
    )
    for win in field_windows(stations):
        print(
            f"{win.closing_speed_km_s:7.2f} {win.branch:>11} {win.mach_sonic:6.2f} "
            f"{win.b_contain_t:7.2f}T {win.b_release_t:7.2f}T {win.b_design_t:6.2f}T "
            f"{win.alfven_mach_now:8.2f} {win.alfven_mach_at_standoff:13.2f} "
            f"{win.required_area_ratio:10.1f} {win.snowplow_overshoot:8.1f}x"
        )
    print("\nThe window exists in every case, and the design sits ABOVE it rather than inside.")
    print("Flux conservation ties the exit field to the flare, so reaching the release ceiling")
    print("means an area ratio near 7 rather than the flown 4 -- flare harder.")
    print("The cost is the last column: one static field cannot be strong for the collision and")
    print("weak for the expansion, so a field weakened to release lets the snowplow past. That is")
    print("what a physical wall in the diverging section would be there to take.")

    cw = common_window(stations)
    print(
        f"\nOne flare for the whole fleet: A/A* between {cw.area_ratio_min:.1f} and "
        f"{cw.area_ratio_max:.1f}, against the flown {cw.flown_area_ratio:.0f}."
    )
    print(
        f"  lower bound set by detachment on the {cw.binding_release_leg} leg; "
        f"upper by containment on {cw.binding_contain_leg}."
    )

    write_stations(stations, args.output)
    print(f"\nwrote {len(stations)} stations -> {args.output}")


if __name__ == "__main__":
    main()
