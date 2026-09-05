"""Does a magnetic extension move detachment inside the nozzle, and what does it buy? (R11, R1)

Three questions land on one run, because they are all properties of the same continued expansion:

- **R1** -- `eta_geom` in the regime that actually applies. `continuum.py` retires `jet.py`'s
  `mu`-conservation model; what replaces it is the ordinary nozzle decomposition, and it needs the
  solved Mach number, which this module already has at every station.
- **R11** -- `M_A(z)` through a prescribed extension past the bag, on our solver rather than the
  paper's 1-D isentropic area-Mach relation, and where the crossing lands.
- **R8's cheap partial** -- the same `M_A(z)` under a *prescribed wall contour* instead of free
  expansion, which is a boundary condition rather than a solver.

# `eta_geom` in a collisional nozzle

`eta_geom` is the share of the plume's speed that points down the thrust axis, `<v_z>/v_g` with
`v_g = sqrt(<v^2>)`. In a continuum exhaust a parcel has a **directed** speed `u` along its
streamline plus **isotropic thermal** motion about it. Thermal motion averages to zero in the
mean and contributes fully to the mean square, so with streamlines inclined at `theta`:

    <v_z> = u <cos theta>,    <v^2> = u^2 + <v_th^2> = u^2 + 3 R_sp T

    eta_geom = <cos theta> / sqrt(1 + 3 R_sp T / u^2)

which is the paper's `<cos theta>/sqrt(1 + 3/(gamma M^2))` with `R_sp T = c_s^2/gamma` substituted
back out. **It is evaluated here without `gamma`**, from `R_sp = k_B/m_bar` with `m_bar` the mean
*heavy-particle* mass out of `eos_water.composition`, because the mixture's effective `gamma`
varies through the expansion and the mass carriers are what set `v_g`. The dissociation and
ionisation stores are not in `<v^2>`; they are `eta_chem`'s business, and double-counting them
here is the error the factorisation exists to prevent.

**Two terms, and they pull opposite ways under every lever.** A deeper flare raises `M` (more of
the thermal store becomes directed) and lowers `<cos theta>` (the streamlines fan). That is why a
magnetic nozzle has an optimum rather than being monotone in `A/A*`, which is R1's correction to
`jet.py`'s docstring and is reproduced here as an output rather than asserted.

# Why the cap costs area ratio, which nobody has priced

`A/A*` is the field ratio by flux conservation. The paper's ADR-0012 caps the peak at 12 T while
P9 holds the exit at 5 T, so the **magnet's own area ratio falls from 4.0 to 2.40**. That is a 40%
cut in the expansion, and the obvious expectation is that it costs `eta_geom`: a smaller area ratio
is a lower exit Mach number, and less of the thermal store gets turned into directed motion.

**It does cost that, and it buys back more.** Both profiles are run rather than argued, and the
cap comes out *ahead*:

- `eta_thermal` falls, as expected -- the exit Mach number drops from 3.23 to 2.42 on the cold
  equilibrium leg, and `eta_thermal` with it, 0.909 to 0.869.
- **`<cos theta>` rises by much more.** A smaller area ratio leaves the plume denser at the exit
  at the same 5 T, so `v_A = B/sqrt(mu0 rho)` is lower and `M_A` is *higher* -- 0.44 against 0.35.
  A higher `M_A` detaches sooner (1.73 exit radii instead of 2.00), and sooner means less fanning.

So the cap is not the compromise ADR-0012 prices it as. Its own table has `eta_geom` moving by
under 0.005 across the cap; in the regime that applies it moves by **+0.06 to +0.15, in the
paper's favour**. The two frameworks disagree because in the `mu` picture `alpha` depends on the
field *ratio*, which the cap preserves at the exit, while in the continuum picture the cap moves
the detachment surface.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

import numpy as np

from puffsat import eos_water, expansion, field, fluxtube

MU0 = field.MU0
K_B = eos_water.K_B

DEFAULT_OUTPUT = Path("data/results/nozzle_extension.csv")

EXIT_FIELD_T = 5.0
"""Field at the magnet exit [T]. P9 calls this a collision requirement, not a standoff one."""

FORWARD_THRUST_FLOOR = 0.324
"""`1/sqrt(1+k)` at `k` = 8.52 -- `sec:mass_interest`'s floor. Below it the growth claim fails."""

SWEPT_TARGET = 0.775
"""`eta_jet` the growth tables sweep to (`sec:jet_efficiency`)."""

ETA_CHEM = {45.58: 0.731, 56.53: 0.835, 65.0: 0.864, 75.0: 0.910}
"""Frozen-chemistry ceiling per closing speed. **A water figure** -- see R13."""

EXTENSION_RATIOS: tuple[tuple[str, float], ...] = (
    ("none", 1.0),
    ("hot-favourable", 2.63),
    ("hot-pessimistic", 2.90),
    ("cold-favourable", 5.60),
    ("cold-pessimistic", 7.45),
)
"""R11's own area ratios past the magnet exit, so the answer is comparable to its table."""


# ---- The continuum `eta_geom` ------------------------------------------------------------------


def thermal_speed_squared(rho: float, temp: float) -> float:
    """`<v_th^2> = 3 k T/m_bar` [m^2/s^2] over the **heavy** particles.

    Electrons are excluded: they carry the thermal energy but essentially none of the mass, and
    `v_g` is defined by energy per unit mass. Including them raises `<v^2>` by the electron share
    of the thermal energy, which is under a percent of the heavy-particle term at these
    ionisation fractions; `main` prints the sensitivity rather than hiding the choice.
    """
    comp = eos_water.composition(rho, temp)
    n_heavy = comp.n_neutral_heavy + comp.n_hp + sum(comp.n_o_ions)
    if n_heavy <= 0.0:
        raise ValueError("no heavy particles; the composition solve failed")
    m_bar = rho / n_heavy
    return 3.0 * K_B * temp / m_bar


def eta_thermal(rho: float, temp: float, speed: float) -> float:
    """`1/sqrt(1 + <v_th^2>/u^2)` -- the directed share of a perfectly aligned exhaust.

    The `<cos theta>` = 1 limit of `eta_geom`, and the ceiling any nozzle geometry can reach at a
    given Mach number.
    """
    if speed <= 0.0:
        return 0.0
    return 1.0 / math.sqrt(1.0 + thermal_speed_squared(rho, temp) / (speed * speed))


def eta_geom(rho: float, temp: float, speed: float, mean_cos_theta: float) -> float:
    """The full continuum `eta_geom`: divergence times thermal alignment."""
    return mean_cos_theta * eta_thermal(rho, temp, speed)


def cone_mean_cos(half_angle_rad: float) -> float:
    """`(1 + cos theta_max)/2` -- the conical-nozzle divergence factor, for a uniform cone.

    The standard result for a source flow filling a cone. Carried because it is what R11's own
    probes used, so the two sides can be compared before the traced value replaces it.
    """
    return 0.5 * (1.0 + math.cos(half_angle_rad))


# ---- One continued expansion --------------------------------------------------------------------


@dataclass(frozen=True)
class Station:
    """One station of a continued expansion, with its magnetic and geometric quantities."""

    closing_speed_km_s: float
    branch: str
    area_ratio: float
    x_m: float
    radius_m: float
    rho: float
    temp_k: float
    pressure_pa: float
    speed_m_s: float
    mach: float
    b_t: float
    v_alfven: float
    alfven_mach: float
    beta: float
    eta_thermal: float


def continued_history(
    closing_speed: float,
    temp_0: float,
    *,
    frozen: bool,
    chamber_field_t: float,
    area_ratio_end: float,
    length_m: float,
    chamber_radius_m: float = expansion.CHAMBER_RADIUS,
    steps: int = 320,
    stride: int = 8,
) -> list[Station]:
    """Solve the isentropic expansion out to `area_ratio_end` and attach the magnetic quantities.

    The field at a station is `B = B_chamber/(A/A*)`, which is flux conservation and the same
    relation the rest of the repo uses -- but with the chamber field as a parameter, because
    ADR-0012's cap moves it and `detachment.load_stations` hardcodes 20 T.
    """
    if frozen:
        y = eos_water.frozen_composition(expansion.BAG_RHO, temp_0)
        eos: expansion.Eos = lambda r, t: eos_water.pressure_energy_frozen(r, t, y)  # noqa: E731
        c_s: expansion.SoundSpeed = lambda r, t: eos_water.sound_speed_frozen(r, t, y)  # noqa: E731
    else:
        eos, c_s = eos_water.pressure_energy, eos_water.sound_speed

    rows = expansion.cooling_history(
        expansion.BAG_RHO,
        temp_0,
        eos,
        c_s,
        area_ratio_end,
        length_m,
        steps=steps,
        expansion_ratio=4096.0,
    )
    sampled = [*rows[::stride], rows[-1]]

    out: list[Station] = []
    for row in sampled:
        b = chamber_field_t / row.area_ratio
        v_a = b / math.sqrt(MU0 * row.rho)
        out.append(
            Station(
                closing_speed_km_s=closing_speed,
                branch="frozen" if frozen else "equilibrium",
                area_ratio=row.area_ratio,
                x_m=row.x,
                radius_m=chamber_radius_m * math.sqrt(row.area_ratio),
                rho=row.rho,
                temp_k=row.temp,
                pressure_pa=row.pressure,
                speed_m_s=row.speed,
                mach=row.mach,
                b_t=b,
                v_alfven=v_a,
                alfven_mach=row.speed / v_a if v_a > 0.0 else math.inf,
                beta=row.pressure / (b * b / (2.0 * MU0)),
                eta_thermal=eta_thermal(row.rho, row.temp, row.speed),
            )
        )
    return out


def alfven_crossing(stations: list[Station]) -> Station | None:
    """The first station at which `M_A` reaches 1, or `None` if the expansion never gets there."""
    for lo, hi in pairwise(stations):
        if lo.alfven_mach <= 1.0 <= hi.alfven_mach:
            return hi
    return None


# ---- The sweep R11 asks for ---------------------------------------------------------------------


@dataclass(frozen=True)
class Case:
    """One (leg, branch, extension) combination, reduced to what the paper needs."""

    closing_speed_km_s: float
    branch: str
    profile: str
    extension: str
    area_ratio_exit: float
    length_m: float
    exit_field_t: float
    exit_radius_m: float
    exit_mach: float
    exit_beta: float
    exit_alfven_mach: float
    detachment_radii: float
    live_fraction: float
    crossing_area_ratio: float | None
    crossing_x_m: float | None
    eta_thermal_exit: float
    mean_cos_theta: float
    eta_geom: float
    eta_jet_water: float


def half_angle(length_m: float, r_start: float, r_end: float) -> float:
    """Flare half-angle [rad] of a straight-walled extension between two radii."""
    if length_m <= 0.0:
        return 0.0
    return math.atan((r_end - r_start) / length_m)


def run_case(
    closing_speed: float,
    temp_0: float,
    *,
    frozen: bool,
    chamber_field_t: float,
    profile_name: str,
    extension_name: str,
    extension_ratio: float,
    fan: Fan,
    magnet_length_m: float = fluxtube.COLUMN_LENGTH_M,
    flare_half_angle_deg: float = 15.0,
) -> tuple[Case, list[Station]]:
    """One leg, one branch, one extension. Returns the summary and the stations behind it.

    **Where `<cos theta>` comes from depends on whether there is an extension**, and that is the
    whole of R11's argument rather than a modelling convenience:

    - *No extension.* The plume leaves the last coil sub-Alfvenic, so it is still tied to the
      field and follows the opening vacuum lines until `M_A` reaches 1. The divergence is then
      whatever those lines have done by the detachment surface, which `fluxtube.downstream_fan`
      traces and `Fan.at` reads off.
    - *With an extension.* The crossing is reached **inside** a controlled flare, so the exhaust
      direction is set by the flare rather than by free fanning, and the conical divergence factor
      applies.
    """
    magnet_ratio = chamber_field_t / EXIT_FIELD_T
    area_ratio_end = magnet_ratio * extension_ratio

    r_magnet_exit = expansion.CHAMBER_RADIUS * math.sqrt(magnet_ratio)
    r_end = expansion.CHAMBER_RADIUS * math.sqrt(area_ratio_end)
    ext_length = (
        0.0
        if extension_ratio <= 1.0
        else (r_end - r_magnet_exit) / math.tan(math.radians(flare_half_angle_deg))
    )
    length = magnet_length_m + ext_length

    stations = continued_history(
        closing_speed,
        temp_0,
        frozen=frozen,
        chamber_field_t=chamber_field_t,
        area_ratio_end=area_ratio_end,
        length_m=length,
    )
    exit_station = stations[-1]
    crossing = alfven_crossing(stations)

    detach_radii = detachment_radii(exit_station.alfven_mach)
    if extension_ratio > 1.0 and exit_station.alfven_mach >= 1.0:
        # The crossing is reached **inside** the flare, so the exhaust direction is set by the
        # flare wall. This is R11's whole claim, and it is checked rather than assumed.
        cos_theta = cone_mean_cos(math.radians(flare_half_angle_deg))
        live = 1.0
    else:
        # Still sub-Alfvenic where the winding stops, extension or not: the plume remains tied to
        # the field and fans freely past the last coil until it can leave. The fan is read in
        # exit radii past the winding's end, which is where a solenoid end's geometry is
        # self-similar, so the magnet's traced fan carries over to an extended one.
        cos_theta = fan.at(detach_radii)
        live = fan.live_at(detach_radii)

    eta_g = eta_geom(exit_station.rho, exit_station.temp_k, exit_station.speed_m_s, cos_theta)
    chem = ETA_CHEM.get(closing_speed, 1.0)

    return (
        Case(
            closing_speed_km_s=closing_speed,
            branch="frozen" if frozen else "equilibrium",
            profile=profile_name,
            extension=extension_name,
            area_ratio_exit=area_ratio_end,
            length_m=length,
            exit_field_t=chamber_field_t / area_ratio_end,
            exit_radius_m=r_end,
            exit_mach=exit_station.mach,
            exit_beta=exit_station.beta,
            exit_alfven_mach=exit_station.alfven_mach,
            detachment_radii=detach_radii,
            live_fraction=live,
            crossing_area_ratio=None if crossing is None else crossing.area_ratio,
            crossing_x_m=None if crossing is None else crossing.x_m,
            eta_thermal_exit=exit_station.eta_thermal,
            mean_cos_theta=cos_theta,
            eta_geom=eta_g,
            eta_jet_water=eta_g * chem,
        ),
        stations,
    )


# ---- Tying the detachment surface to the traced fan --------------------------------------------


@dataclass(frozen=True)
class Fan:
    """`<cos theta>` against distance past the exit, from the real field trace."""

    radii_out: tuple[float, ...]
    mean_cos_theta: tuple[float, ...]
    live_fraction: tuple[float, ...]

    def at(self, radii_out: float) -> float:
        """`<cos theta>` at a detachment surface, in exit radii past the exit plane."""
        return float(np.interp(radii_out, np.array(self.radii_out), np.array(self.mean_cos_theta)))

    def live_at(self, radii_out: float) -> float:
        """Share of tubes still running forward there -- the rest have turned."""
        return float(np.interp(radii_out, np.array(self.radii_out), np.array(self.live_fraction)))


def traced_fan(chamber_field_t: float, *, n_radial: int = 8) -> Fan:
    """Trace the plume's fan past the last coil for a given chamber field.

    The exit radius follows from the magnet's own area ratio, `r_exit = r_chamber sqrt(B*/B_exit)`,
    so the cap moves it: 6.00 m uncapped, 4.70 m at 12 T.
    """
    stacks = fluxtube.build_stacks(cap_t=chamber_field_t)
    stack, _profile = stacks["flared-capped"]
    r_exit = expansion.CHAMBER_RADIUS * math.sqrt(chamber_field_t / EXIT_FIELD_T)
    rows = fluxtube.downstream_fan(stack, r_exit=r_exit, n_radial=n_radial)
    usable = [r for r in rows if not math.isnan(r.mean_cos_theta)]
    return Fan(
        radii_out=tuple(r.radii_out for r in usable),
        mean_cos_theta=tuple(r.mean_cos_theta for r in usable),
        live_fraction=tuple(r.live_fraction for r in usable),
    )


def detachment_radii(exit_alfven_mach: float) -> float:
    """Exit radii past the exit plane at which `M_A` reaches 1, on the paper's own scalings.

    `detachment.downstream_crossing`, restated here so `extension.py` does not import a module
    that hardcodes the 20 T chamber. Free expansion at constant speed gives `rho ~ R^-3`, a
    dipole gives `B ~ R^-3`, so `M_A ~ R^{3/2}` and the crossing is `(1/M_A_exit)^{2/3}`.
    """
    if exit_alfven_mach <= 0.0:
        return math.inf
    if exit_alfven_mach >= 1.0:
        return 1.0
    return float((1.0 / exit_alfven_mach) ** (2.0 / 3.0))


# ---- The sweep -----------------------------------------------------------------------------------

PROFILES: tuple[tuple[str, float], ...] = (
    ("flown-20T", 20.0),
    ("capped-12T", field.ADR0012_CAP_T),
)
"""The chamber field before and after ADR-0012's cap. The cap moves `A/A*`, so it moves `M`."""

LEGS: tuple[tuple[float, float], ...] = (
    (45.58, 14700.0),
    (75.0, 26200.0),
)
"""The two legs the paper flies: leg 1's cold end (overtake) and leg 2 (head-on departure)."""


CSV_HEADER = (
    "closing_speed_km_s",
    "branch",
    "profile",
    "extension",
    "area_ratio_exit",
    "length_m",
    "exit_field_T",
    "exit_radius_m",
    "exit_mach",
    "exit_beta",
    "exit_alfven_mach",
    "detachment_radii",
    "live_fraction",
    "eta_thermal_exit",
    "mean_cos_theta",
    "eta_geom",
    "eta_jet_water",
)


def write_cases(rows: list[Case], path: Path = DEFAULT_OUTPUT) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(CSV_HEADER)
        for c in rows:
            writer.writerow(
                [
                    f"{c.closing_speed_km_s:g}",
                    c.branch,
                    c.profile,
                    c.extension,
                    f"{c.area_ratio_exit:.4f}",
                    f"{c.length_m:.2f}",
                    f"{c.exit_field_t:.4f}",
                    f"{c.exit_radius_m:.3f}",
                    f"{c.exit_mach:.4f}",
                    f"{c.exit_beta:.6f}",
                    f"{c.exit_alfven_mach:.4f}",
                    f"{c.detachment_radii:.4f}",
                    f"{c.live_fraction:.4f}",
                    f"{c.eta_thermal_exit:.4f}",
                    f"{c.mean_cos_theta:.4f}",
                    f"{c.eta_geom:.4f}",
                    f"{c.eta_jet_water:.4f}",
                ]
            )


def main() -> None:
    """R1's replacement `eta_geom`, and R11's extension, on the solved expansion."""
    parser = argparse.ArgumentParser(description="The magnetic extension and eta_geom (R11, R1)")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    fans = {name: traced_fan(b) for name, b in PROFILES}

    print("== The traced fan past the last coil (R2 downstream; R1's divergence term) ==")
    print("Inside the winding the bore is axial and <cos theta> is within 0.2% of 1. All of")
    print("the divergence is downstream, so it is a detachment-surface question.\n")
    print(f"{'profile':>12} {'radii out':>10} {'<cos theta>':>12} {'tubes still forward':>20}")
    for name, _b in PROFILES:
        for radii in (0.0, 0.5, 1.0, 1.44, 2.0):
            print(
                f"{name:>12} {radii:10.2f} {fans[name].at(radii):12.4f} "
                f"{100 * fans[name].live_at(radii):19.0f}%"
            )

    rows: list[Case] = []
    for profile_name, b_chamber in PROFILES:
        for speed, temp_0 in LEGS:
            for frozen in (False, True):
                for ext_name, ext_ratio in EXTENSION_RATIOS:
                    if profile_name == "flown-20T" and ext_name != "none":
                        continue  # the cap is flown; extensions are only priced on it
                    case, _stations = run_case(
                        speed,
                        temp_0,
                        frozen=frozen,
                        chamber_field_t=b_chamber,
                        profile_name=profile_name,
                        extension_name=ext_name,
                        extension_ratio=ext_ratio,
                        fan=fans[profile_name],
                    )
                    rows.append(case)

    print("\n== R1: eta_geom in the regime that applies, and what the cap costs ==")
    print(
        f"{'leg':>7} {'branch':>7} {'profile':>11} {'A/A*':>6} {'M':>6} {'M_A':>6} "
        f"{'detach':>7} {'<cos>':>6} {'eta_th':>7} {'eta_geom':>9} {'eta_jet':>8}"
    )
    for c in rows:
        if c.extension != "none":
            continue
        print(
            f"{c.closing_speed_km_s:7.2f} {c.branch[:6]:>7} {c.profile:>11} "
            f"{c.area_ratio_exit:6.2f} {c.exit_mach:6.3f} {c.exit_alfven_mach:6.3f} "
            f"{c.detachment_radii:7.2f} {c.mean_cos_theta:6.3f} {c.eta_thermal_exit:7.4f} "
            f"{c.eta_geom:9.4f} {c.eta_jet_water:8.4f}"
        )

    print("\n== R11: what the extension buys, on the capped profile ==")
    print(
        f"{'leg':>7} {'branch':>7} {'extension':>17} {'total L':>8} {'A/A*':>6} {'B_exit':>7} "
        f"{'M':>6} {'M_A':>6} {'eta_geom':>9} {'eta_jet':>8} {'vs 0.775':>9}"
    )
    for c in rows:
        if c.profile != "capped-12T":
            continue
        verdict = "clears" if c.eta_jet_water >= SWEPT_TARGET else "misses"
        print(
            f"{c.closing_speed_km_s:7.2f} {c.branch[:6]:>7} {c.extension:>17} "
            f"{c.length_m:7.1f}m {c.area_ratio_exit:6.2f} {c.exit_field_t:6.2f}T "
            f"{c.exit_mach:6.3f} {c.exit_alfven_mach:6.3f} {c.eta_geom:9.4f} "
            f"{c.eta_jet_water:8.4f} {verdict:>9}"
        )

    floor_fails = [c for c in rows if c.eta_jet_water < FORWARD_THRUST_FLOOR]
    print(f"\nforward-thrust floor {FORWARD_THRUST_FLOOR}: {len(floor_fails)} of {len(rows)} fail")
    best = max(rows, key=lambda c: c.eta_geom)
    print(
        f"best eta_geom  : {best.eta_geom:.4f} ({best.extension}, {best.branch}, "
        f"{best.closing_speed_km_s:g} km/s)"
    )

    write_cases(rows, args.output)
    print(f"\nwrote {len(rows)} rows -> {args.output}")


if __name__ == "__main__":
    main()
