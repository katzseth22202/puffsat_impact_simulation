"""Result 17: the lossless gas spring behind the growth-push plate -- its loss, heat and ring-down.

    PYTHONPATH=python uv run python -m puffsat.xray_collision.gas_spring

The parent's S46 (`docs/fast_growth_wave_asks_for_impact_sim.md`, commit 5d1379c) drops ADR-0064's
absorber-set plate (`I^2 / (2 s M a)`, a plate stopped within each pulse, which would have to shed
2.2 GJ a pulse). The paper's `sec:pusher_plate_damping` spring replaces it. The plate swings on a
nearly lossless gas spring, each pulse reverses its swing, and the law fixes the plate's mass times
its stroke: about 507 t m for 7.35 MN s at 2 Hz. On an 8 m stroke that is a 63 t plate.

Four things are computed here.

1. **The cycle, solved nonlinearly.** The paper's law is a linear spring. A gas column compressed
   by most of its stroke is not linear, so the periodic orbit is shot for an adiabatic column
   `F = F0 (H / (H - x))^gamma` on the plate-craft reduced mass, kicked by `J / m` each pulse.
   Out come the column height `H`, the stroke, the force swing and the pressure swing. The
   column's height is fixed by the dynamics. Pressure sets only the footprint, so the bag's fibre
   and gas masses (both `~ P V`) are the same at any pressure.
2. **The loss per cycle**, by three channels:
   - *Thermal hysteresis.* The gas heats on compression and loses some of that heat to the bag
     wall before it re-expands. The lumped gas is integrated cycle by cycle against a wall held at
     `T_w`, with a heat-transfer coefficient bracketed three ways: the laminar oscillating
     boundary layer (`k / delta_kappa`, Swift), turbulent forced convection at the gas's own speed
     (Colburn analogy), and Kornhauser & Smith's measured gas-spring correlation
     `Nu = 0.56 Pe_omega^0.69` (extrapolated far past its data). The largest is the design case.
   - *Fibre hysteresis.* The hoop fibres cycle with the pressure. Loss `pi eta sigma_a^2 / E` per
     unit volume, with the loss factor `eta` an assumed 0.01-0.05 bracket (no measured value for
     Vectran fabric at 2 Hz was found).
   - *Wall friction* of the gas's breathing flow (turbulent skin friction, the same velocity).
3. **The bags' heat over the push.** The losses go into the bag wall, fibre and gas. The rise is
   given with no cooling (an upper bound) and with radiation off the bundle's outside.
4. **The plate's ring-down** (S46 item 2). A free axisymmetric Kirchhoff plate (Rayleigh-Ritz on
   Legendre polynomials in `(r/R)^2`) with the skirt as a rim mass. It is loaded by the sheet's
   2-D impulse profile (Result 15's k = 7.07 run), and the residual after rigid-body motion rings
   in the plate's own modes. The static bending from where the bags push is also given, and the
   sandwich depth a ring of bags would need.

Gas properties are CoolProp's at the bag's mean state. Bag-fabric, bladder and loss-factor values
are datasheet-class or assumed; each is marked where it is set.
"""

from __future__ import annotations

import csv
import math
from dataclasses import asdict, dataclass
from functools import cache
from pathlib import Path

import CoolProp.CoolProp as CP
import numpy as np
from numpy.polynomial import legendre
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.linalg import eigh
from scipy.optimize import brentq, fsolve

from puffsat.xray_collision import charge, coated_sheet, fast_wave

G0 = 9.80665
SIGMA_SB = 5.670374419e-8

# --- The push (ADR-0064 amendment; parent S45-S46).
PULSE = fast_wave.PULSE  # [N s], 7.35 MN s
CADENCE = fast_wave.CADENCE  # [Hz]
CRAFT = coated_sheet.LAUNCH_UNIT / 1.330  # [kg], the unit at the push's end (mass ratio 1.330)
PLATE = 63.0e3  # [kg], the parent's estimate at an 8 m stroke
STROKE = 8.0  # [m], Seth's
THETA_PARENT = 2.0 * math.pi / 3.0
#: The push's pulse count: the parent's 4.4 TJ / 2.2 GJ.
CYCLES = 2000
T_WALL = 300.0  # [K], the bag at the start

# --- Bag materials. Datasheet-class, not sourced here.
#: Vectran HT yarn: ~3.0 GPa and ~65 GPa (the repo's dry_wrap.py value), 1400 kg/m^3, c_p taken
#: as 1200 J/kg/K (liquid-crystal polyesters run 1.0-1.3 kJ/kg/K).
FIBRE_RHO = 1400.0
FIBRE_E = 65.0e9
FIBRE_ULT = 3.0e9
FIBRE_SF = 4.0  # yarn-to-fabric translation, seams and 2000-cycle fatigue together (assumed)
FIBRE_CP = 1200.0
#: Fibre loss factor bracket (assumed; aramid and LCP fibres are usually quoted at 1-5%).
FIBRE_ETA = (0.01, 0.05)
#: A gas-tight film bladder inside the fabric: ~140 um of polyimide-class film (assumed).
BLADDER_AREAL = 0.2  # [kg/m^2]
BLADDER_CP = 1090.0
#: The bundle's outside radiates to a cold sink (assumed emissivity and sink).
EMISSIVITY = 0.85
T_SINK = 250.0
#: Service limits for the bag's materials (datasheet-class): Vectran's long-term rating is
#: usually put near 150 C, aramid's near 180-200 C, a polyimide bladder's near 400 C.
T_LIMITS = {"Vectran": 423.0, "aramid": 463.0, "polyimide bladder": 673.0}

GASES = {"nitrogen": "Nitrogen", "helium": "Helium", "hydrogen": "Hydrogen", "argon": "Argon"}


@dataclass(frozen=True)
class Gas:
    name: str
    r: float  # [J/kg/K]
    gamma: float
    cp: float
    mu: float  # [Pa s]
    k: float  # [W/m/K]
    rho: float  # [kg/m^3] at the mean state


@cache
def gas(name: str, pressure: float, temp: float = T_WALL) -> Gas:
    """CoolProp's properties at `pressure` [Pa] and `temp` [K]."""
    fluid = GASES[name]

    def prop(key: str) -> float:
        return float(CP.PropsSI(key, "P", pressure, "T", temp, fluid))

    cp, cv = prop("Cpmass"), prop("Cvmass")
    return Gas(name, cp - cv, cp / cv, cp, prop("V"), prop("L"), prop("Dmass"))


# --- 1. The cycle ---------------------------------------------------------------------------


def reduced(m: float, ship: float = CRAFT) -> float:
    return m * ship / (m + ship)


def linear_stroke(m: float, theta: float) -> float:
    """The paper's linear-spring stroke [m]: `(Ma / m nu^2) tan(theta/4) / (2 theta)` with
    `Ma = J nu`, i.e. `(J / m) tan(theta/4) / (2 omega)`; `theta` counts on the reduced mass."""
    return PULSE * CADENCE / (m * CADENCE**2) * math.tan(theta / 4.0) / (2.0 * theta)


@dataclass(frozen=True)
class Orbit:
    """The steady adiabatic cycle. `x` is the column's compression from where its force is the
    craft's mean force `f0`; the kick adds `J / m` to the compression rate."""

    m: float
    column: float  # [m], H at x = 0
    gamma: float
    f0: float
    x0: float  # just after the kick
    v0: float
    x_min: float
    x_max: float
    f_min: float
    f_max: float
    f_mean: float
    theta_lin: float  # the linearized spring's theta at x = 0

    @property
    def stroke(self) -> float:
        return self.x_max - self.x_min

    @property
    def peak_over_mean(self) -> float:
        return self.f_max / self.f_mean

    @property
    def pressure_ratio(self) -> float:
        return self.f_max / self.f_min


def orbit(m: float, column: float, gamma: float, ship: float = CRAFT) -> Orbit:
    """Shoot the periodic orbit: `(x0, v0)` after a kick returns to itself one period later."""
    mu = reduced(m, ship)
    f0 = PULSE * CADENCE * ship / (m + ship)
    period = 1.0 / CADENCE
    kick = PULSE / m

    def rhs(t: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        f = f0 * (column / (column - y[0])) ** gamma
        return np.array([y[1], -f / mu])

    def residual(z: NDArray[np.float64]) -> NDArray[np.float64]:
        end = solve_ivp(rhs, (0.0, period), z, rtol=1e-10, atol=1e-10).y[:, -1]
        return np.array([end[0] - z[0], end[1] + kick - z[1]])

    k_lin = gamma * f0 / column
    theta = math.sqrt(k_lin / mu) / CADENCE
    # Linear guess: velocity +kick/2 after the kick, at the cycle's centre.
    z = fsolve(residual, np.array([-0.0, 0.5 * kick]), xtol=1e-12)
    ts = np.linspace(0.0, period, 4001)
    xs = solve_ivp(rhs, (0.0, period), z, rtol=1e-10, atol=1e-10, t_eval=ts).y[0]
    fs = f0 * (column / (column - xs)) ** gamma
    return Orbit(
        m,
        column,
        gamma,
        f0,
        float(z[0]),
        float(z[1]),
        float(xs.min()),
        float(xs.max()),
        float(fs.min()),
        float(fs.max()),
        float(np.trapezoid(fs, ts) / period),
        theta,
    )


def column_for_stroke(m: float, stroke: float, gamma: float) -> Orbit:
    """The column height whose steady orbit sweeps `stroke` at plate mass `m`."""
    lo = 0.6 * stroke  # the column must be longer than the compression it takes
    hi = 400.0
    h = brentq(lambda c: orbit(m, c, gamma).stroke - stroke, lo, hi, xtol=1e-4)
    return orbit(m, h, gamma)


def column_for_theta(m: float, theta: float, gamma: float) -> Orbit:
    """The column whose linearized spring has `theta` at the mean force."""
    mu = reduced(m)
    f0 = PULSE * CADENCE * CRAFT / (m + CRAFT)
    k = mu * (theta * CADENCE) ** 2
    return orbit(m, gamma * f0 / k, gamma)


# --- 2. The bag -----------------------------------------------------------------------------


@dataclass(frozen=True)
class Bag:
    """`cells` cylinders side by side under the plate, mean pressure `pressure`."""

    label: str
    pressure: float  # [Pa] at x = 0
    cells: int
    gas: str


BAGS = (
    Bag("whole floor, nitrogen", 52.0e3, 19, "nitrogen"),
    Bag("1 MPa, 4 cells, nitrogen", 1.0e6, 4, "nitrogen"),
    Bag("whole floor, helium", 52.0e3, 19, "helium"),
    Bag("1 MPa, 4 cells, helium", 1.0e6, 4, "helium"),
    Bag("whole floor, hydrogen", 52.0e3, 19, "hydrogen"),
    Bag("whole floor, argon", 52.0e3, 19, "argon"),
)


@dataclass(frozen=True)
class Geometry:
    area: float  # [m^2], total footprint
    cell_radius: float
    volume: float  # [m^3] at x = 0
    wall: float  # [m^2], all cells' side walls at x = 0 (ends are the plate and the deck)
    ends: float  # [m^2]
    outside: float  # [m^2], the bundle's outer side, which sees space
    fibre: float  # [kg]
    bladder: float  # [kg]
    gas: float  # [kg]


def geometry(bag: Bag, o: Orbit) -> Geometry:
    area = o.f0 / bag.pressure
    r = math.sqrt(area / (bag.cells * math.pi))
    h_max = o.column - o.x_min
    volume = area * o.column
    wall = bag.cells * 2.0 * math.pi * r * o.column
    ends = 2.0 * area
    outside = 2.0 * math.pi * math.sqrt(area / math.pi) * o.column
    # Hoop fibres only (the ends are the plate and the deck): 2 P V rho / sigma, at the peak
    # pressure over the longest column.
    sigma = FIBRE_ULT / FIBRE_SF
    fibre = 2.0 * (bag.pressure * o.f_max / o.f0) * area * h_max * FIBRE_RHO / sigma
    bladder = BLADDER_AREAL * (bag.cells * 2.0 * math.pi * r * h_max + ends)
    g = gas(bag.gas, bag.pressure)
    return Geometry(area, r, volume, wall, ends, outside, fibre, bladder, g.rho * volume)


# --- 3. The loss per cycle ------------------------------------------------------------------


def h_laminar(g: Gas) -> float:
    """`k / delta_kappa` at the pulse frequency: the oscillating thermal boundary layer."""
    alpha = g.k / (g.rho * g.cp)
    delta = math.sqrt(2.0 * alpha / (2.0 * math.pi * CADENCE))
    return g.k / delta


def gas_speed(o: Orbit) -> float:
    """The breathing flow's mean speed [m/s]: the gas moves at `xdot z / H`, so its column-mean
    speed is `|xdot| / 2`; `|xdot|` is taken at the orbit's RMS, `~ kick / (2 sqrt 3)`."""
    return 0.5 * PULSE / o.m / (2.0 * math.sqrt(3.0))


def skin_friction(g: Gas, u: float, length: float) -> float:
    """Turbulent flat-plate mean `C_f = 0.074 Re^-0.2`."""
    return float(0.074 * (g.rho * u * length / g.mu) ** -0.2)


def h_turbulent(g: Gas, o: Orbit) -> float:
    """Colburn analogy, `St Pr^(2/3) = C_f / 2`, at the breathing speed along the column."""
    u = gas_speed(o)
    pr = g.cp * g.mu / g.k
    return float(g.rho * g.cp * u * 0.5 * skin_friction(g, u, o.column) * pr ** (-2.0 / 3.0))


def h_kornhauser(g: Gas, geo: Geometry) -> float:
    """Kornhauser & Smith (1994): `Nu = 0.56 Pe_omega^0.69` on `D_h = 4 V / S`, with
    `Pe_omega = omega D_h^2 / (4 alpha)`. Their data reach Pe ~ 1e4; this is far past it."""
    alpha = g.k / (g.rho * g.cp)
    d_h = 4.0 * geo.volume / (geo.wall + geo.ends)
    pe = 2.0 * math.pi * CADENCE * d_h**2 / (4.0 * alpha)
    return float(0.56 * pe**0.69 * g.k / d_h)


def thermal_loss(bag: Bag, o: Orbit, h: float, cycles: int = 400) -> float:
    """Mechanical energy lost per cycle [J] to heat exchange with a wall at `T_WALL`.

    The lumped gas (mass `m_g`, volume `A (H - x)`) and the plate are integrated together:
    `m_g c_v dT/dt = -P dV/dt - h S (T - T_w)`, `F = P A`. In the periodic state the heat the
    wall takes per cycle is the work the kicks put in beyond the adiabatic orbit's."""
    g = gas(bag.gas, bag.pressure)
    geo = geometry(bag, o)
    area = geo.area
    cv = g.cp / g.gamma
    mu = reduced(o.m)
    m_g = bag.pressure * geo.volume / (g.r * T_WALL)
    s = geo.wall + geo.ends

    def rhs(t: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        x, v, temp, _ = y
        vol = area * (o.column - x)
        p = m_g * g.r * temp / vol
        dvol = -area * v
        return np.array(
            [
                v,
                -p * area / mu,
                (-p * dvol - h * s * (temp - T_WALL)) / (m_g * cv),
                h * s * (temp - T_WALL),
            ]
        )

    # Start on the adiabatic orbit; the gas there is at T_w where x = 0.
    temp0 = T_WALL * (o.column / (o.column - o.x0)) ** (g.gamma - 1.0)
    y = np.array([o.x0, o.v0, temp0, 0.0])
    period = 1.0 / CADENCE
    last = math.nan
    for i in range(cycles):
        sol = solve_ivp(rhs, (0.0, period), y, rtol=1e-9, atol=1e-9)
        end = sol.y[:, -1]
        q_cycle = float(end[3] - y[3])
        y = np.array([end[0], end[1] + PULSE / o.m, end[2], end[3]])
        if i > 20 and abs(q_cycle - last) <= 1e-3 * abs(q_cycle):
            return q_cycle
        last = q_cycle
    return last


def fibre_loss(bag: Bag, o: Orbit, eta: float) -> float:
    """[J] per cycle: `pi eta sigma_a^2 / E` over the fibres' volume, the stress swinging with
    the pressure from its least to its peak (the fibres are rated at the peak)."""
    geo = geometry(bag, o)
    sigma = FIBRE_ULT / FIBRE_SF
    amp = 0.5 * sigma * (1.0 - o.f_min / o.f_max)
    return math.pi * eta * amp**2 / FIBRE_E * geo.fibre / FIBRE_RHO


def friction_loss(bag: Bag, o: Orbit) -> float:
    """[J] per cycle: `1/2 rho C_f u^3` over the side walls, with the column's `u = xdot z / H`
    (its cube averages to `xdot^3 / 4`), `xdot` cubed and averaged on a triangle wave."""
    g = gas(bag.gas, bag.pressure)
    geo = geometry(bag, o)
    kick = PULSE / o.m
    mean_cube = (0.5 * kick) ** 3 / 4.0  # |v|^3 on a velocity ramp from -kick/2 to +kick/2
    cf = skin_friction(g, gas_speed(o), o.column)
    return 0.5 * g.rho * cf * mean_cube / 4.0 * geo.wall / CADENCE


# --- 4. The bags' heat over the push --------------------------------------------------------


def heat_capacity(bag: Bag, o: Orbit) -> float:
    g = gas(bag.gas, bag.pressure)
    geo = geometry(bag, o)
    return geo.fibre * FIBRE_CP + geo.bladder * BLADDER_CP + geo.gas * g.cp / g.gamma


def bag_temperature(bag: Bag, o: Orbit, per_cycle: float, radiate: bool) -> float:
    """[K] at the push's end: the losses into the bag's lumped heat capacity, with or without
    radiation off the bundle's outside to `T_SINK` (the start is held at T_WALL in both)."""
    c = heat_capacity(bag, o)
    power = per_cycle * CADENCE
    out = geometry(bag, o).outside if radiate else 0.0
    q0 = EMISSIVITY * SIGMA_SB * out * (T_WALL**4 - T_SINK**4)  # balanced by sunlight at start

    def rhs(t: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        loss = EMISSIVITY * SIGMA_SB * out * (y[0] ** 4 - T_SINK**4) - q0
        return np.array([(power - loss) / c])

    sol = solve_ivp(rhs, (0.0, CYCLES / CADENCE), [T_WALL], rtol=1e-8)
    return float(sol.y[0, -1])


@dataclass(frozen=True)
class SpringRow:
    bag: str
    gas: str
    pressure_kpa: float
    column_m: float
    stroke_m: float
    theta_lin: float
    peak_over_mean: float
    pressure_ratio: float
    footprint_m2: float
    cell_radius_m: float
    gas_t: float
    fibre_t: float
    bladder_t: float
    h_laminar: float  # [W/m^2/K]
    h_turbulent: float
    h_kornhauser: float
    swing_mj: float  # J^2 / 8m, the parent's "0.1 GJ"
    thermal_mj: float  # at the largest h
    thermal_laminar_mj: float
    fibre_mj_lo: float
    fibre_mj_hi: float
    friction_mj: float
    total_mj: float  # largest h, high fibre eta
    loss_share: float  # of the swing, per cycle
    push_gj: float  # over CYCLES
    rise_k_no_cooling: float
    rise_k_radiating: float
    gas_pressure_rise: float  # the gas's P at fixed volume, from the end temperature


def spring_rows(o: Orbit) -> list[SpringRow]:
    out = []
    for bag in BAGS:
        oo = column_for_stroke(o.m, STROKE, gas(bag.gas, bag.pressure).gamma)
        g = gas(bag.gas, bag.pressure)
        geo = geometry(bag, oo)
        hs = (h_laminar(g), h_turbulent(g, oo), h_kornhauser(g, geo))
        th = thermal_loss(bag, oo, max(hs))
        th_lam = thermal_loss(bag, oo, hs[0])
        fl = [fibre_loss(bag, oo, e) for e in FIBRE_ETA]
        fr = friction_loss(bag, oo)
        total = th + fl[1] + fr
        swing = PULSE**2 / (8.0 * oo.m)
        t_hot = bag_temperature(bag, oo, total, radiate=False)
        t_rad = bag_temperature(bag, oo, total, radiate=True)
        out.append(
            SpringRow(
                bag.label,
                bag.gas,
                bag.pressure / 1e3,
                oo.column,
                oo.stroke,
                oo.theta_lin,
                oo.peak_over_mean,
                oo.pressure_ratio,
                geo.area,
                geo.cell_radius,
                geo.gas / 1e3,
                geo.fibre / 1e3,
                geo.bladder / 1e3,
                hs[0],
                hs[1],
                hs[2],
                swing / 1e6,
                th / 1e6,
                th_lam / 1e6,
                fl[0] / 1e6,
                fl[1] / 1e6,
                fr / 1e6,
                total / 1e6,
                total / swing,
                total * CYCLES / 1e9,
                t_hot - T_WALL,
                t_rad - T_WALL,
                t_hot / T_WALL,
            )
        )
    return out


@dataclass(frozen=True)
class DesignRow:
    """The plate-mass / column trade on nitrogen: the paper's linear law against the solve."""

    plate_t: float
    theta_lin: float
    stroke_linear_m: float
    stroke_m: float
    column_m: float
    peak_over_mean: float
    pressure_ratio: float
    gas_and_fibre_pv_mj: float  # P V at x = 0 = f0 H


def design_rows() -> list[DesignRow]:
    gamma = gas("nitrogen", 52.0e3).gamma
    out = []
    for m in (40.0e3, 63.0e3, 100.0e3):
        for theta in (math.pi / 2.0, THETA_PARENT, 5.0 * math.pi / 6.0, math.pi):
            o = column_for_theta(m, theta, gamma)
            out.append(
                DesignRow(
                    m / 1e3,
                    theta,
                    linear_stroke(m, theta),
                    o.stroke,
                    o.column,
                    o.peak_over_mean,
                    o.pressure_ratio,
                    o.f0 * o.column / 1e6,
                )
            )
    return out


# --- 5. The plate's ring-down ---------------------------------------------------------------

STEEL_RHO = 7850.0
STEEL_E = 205.0e9
STEEL_NU = 0.29
PLATE_R = 9.5  # [m], the flown push plate (Result 13)
SKIRT = 2.0  # [m]
#: The flown sheet at 100 um: 7.35 MN s / (k v_e) = 28.3 kg of Kapton over ~8 m radius.
SHEET_FLOWN = 8.0
SHEET_K = 7.07  # the 9.86 R_sun headline's push
#: How long the pulse loads the face: the merged gas crosses the 2 m standoff at its `u_rms`,
#: `w sqrt(k) / (1 + k)` = 32.7 km/s at 99.2 km/s and k = 7.07, so ~61 us. Modes faster than
#: this are loaded quasi-statically, so the rectangular pulse's `sinc` cuts them off.
PULSE_TIME = 2.0 / (99.2e3 * math.sqrt(SHEET_K) / (1.0 + SHEET_K))
#: Damping on the bending modes between pulses (assumed bracket: welded steel against the bags).
ZETAS = (0.005, 0.02)
N_BASIS = 24
_GAUSS_X, _GAUSS_W = legendre.leggauss(600)
_R = 0.5 * PLATE_R * (_GAUSS_X + 1.0)
_W = 0.5 * PLATE_R * _GAUSS_W * 2.0 * math.pi * _R

PROFILE_JSONL = charge.OUT_DIR / "plate_profile_2d.jsonl"


@cache
def _profiles() -> dict[float, dict[str, object]]:
    return {float(str(r["sheet_radius"])): r for r in charge.load_jsonl(PROFILE_JSONL)}


def capture(sheet: float) -> float:
    """`J / (m v)` of a sheet of radius `sheet` on the flown plate (k = 7.07, 2 m standoff)."""
    return float(str(_profiles()[sheet]["impulse_per_mv"]))


def face_coeff(sheet: float) -> float:
    """`P_peak / (sigma v^2)` of that sheet on the flown plate."""
    return float(str(_profiles()[sheet]["peak_per_mv2"])) * math.pi * sheet**2


@cache
def impulse_profile(sheet: float) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """`(r_edges, j)` on the floor: the 2-D cumulative profile differenced on its 0.25 m rings,
    normalized to PULSE inside PLATE_R (the skirt's axial share, < 2%, is dropped)."""
    prof = _profiles()[sheet]["profile"]
    assert isinstance(prof, list)
    rs = np.array([0.0] + [float(p[0]) for p in prof])
    cum = np.array([0.0] + [float(p[1]) for p in prof])
    keep = rs <= PLATE_R + 1e-9
    rs, cum = rs[keep], cum[keep]
    j = np.diff(cum) / (math.pi * np.diff(rs**2))
    return rs, j * PULSE / cum[-1]


def j_at(r: NDArray[np.float64], sheet: float) -> NDArray[np.float64]:
    """Impulse per area [N s/m^2] at `r`, piecewise constant on the 0.25 m rings."""
    e, j = impulse_profile(sheet)
    idx = np.clip(np.searchsorted(e, r, side="right") - 1, 0, len(j) - 1)
    out: NDArray[np.float64] = j[idx]
    return out


@dataclass(frozen=True)
class Section:
    """The plate's floor, its mass spread evenly: solid steel (`depth` 0) or a sandwich of two
    steel faces `depth` apart (all the mass in the faces, so the stiffness is a lower bound).
    With `skirt_on_floor` the skirt shares the gauge and rides the floor's edge as a rim mass;
    otherwise it is held on the craft's frame (the parent's argon plate, its P10) and all of
    `plate` is floor."""

    label: str
    plate: float  # [kg], the sprung mass
    depth: float
    skirt_on_floor: bool

    @property
    def rim(self) -> float:
        if not self.skirt_on_floor:
            return 0.0
        return self.plate * 2.0 * PLATE_R * SKIRT / (PLATE_R**2 + 2.0 * PLATE_R * SKIRT)

    @property
    def areal(self) -> float:
        return (self.plate - self.rim) / (math.pi * PLATE_R**2)

    @property
    def rigidity(self) -> float:
        if self.depth == 0.0:
            h = self.areal / STEEL_RHO
            return STEEL_E * h**3 / (12.0 * (1.0 - STEEL_NU**2))
        tf = 0.5 * self.areal / STEEL_RHO
        return STEEL_E * tf * self.depth**2 / (2.0 * (1.0 - STEEL_NU**2))

    @property
    def stress_per_moment(self) -> float:
        """Surface (face) stress per unit bending moment per width [1/m^2]."""
        if self.depth == 0.0:
            return 6.0 / (self.areal / STEEL_RHO) ** 2
        return 1.0 / (0.5 * self.areal / STEEL_RHO * self.depth)

    @property
    def face_mm(self) -> float:
        """The solid floor, or each sandwich face [mm]."""
        t = self.areal / STEEL_RHO
        return 1e3 * (t if self.depth == 0.0 else 0.5 * t)


def _basis(n: int, r: NDArray[np.float64]) -> tuple[NDArray[np.float64], ...]:
    """`phi, phi'', phi'/r` of `P_n(2 (r/R)^2 - 1)` in r."""
    c = np.zeros(n + 1)
    c[n] = 1.0
    rho = r / PLATE_R
    s = 2.0 * rho**2 - 1.0
    p = legendre.legval(s, c)
    dp = legendre.legval(s, legendre.legder(c))
    ddp = legendre.legval(s, legendre.legder(c, 2))
    ds = 4.0 * rho / PLATE_R
    return tuple(
        np.asarray(a, dtype=np.float64)
        for a in (p, ddp * ds**2 + dp * 4.0 / PLATE_R**2, dp * 4.0 / PLATE_R**2)
    )


def _edge(n: int) -> NDArray[np.float64]:
    return np.array([_basis(i, np.array([PLATE_R]))[0][0] for i in range(n)])


@dataclass(frozen=True)
class PlateModel:
    section: Section
    omega: NDArray[np.float64]  # [rad/s], elastic modes
    modes: NDArray[np.float64]  # coefficients, M-normalized, one column per elastic mode


@cache
def plate_model(sec: Section, n_basis: int = N_BASIS) -> PlateModel:
    """Rayleigh-Ritz on a free axisymmetric Kirchhoff plate with the rim as an edge mass. The
    bags are left out: their stiffness (`k / A` ~ 4 kN/m^3) is far below the plate's."""
    b = [_basis(i, _R) for i in range(n_basis)]
    phi = np.array([bi[0] for bi in b])
    dd = np.array([bi[1] for bi in b])
    d1 = np.array([bi[2] for bi in b])
    edge = _edge(n_basis)
    mass = (phi * (_W * sec.areal)) @ phi.T + sec.rim * np.outer(edge, edge)
    wd = _W * sec.rigidity
    stiff = (dd * wd) @ dd.T + (d1 * wd) @ d1.T + STEEL_NU * ((dd * wd) @ d1.T + (d1 * wd) @ dd.T)
    vals, vecs = eigh(stiff, mass)
    # Drop the rigid translation (the least eigenvalue, ~0).
    return PlateModel(sec, np.sqrt(np.maximum(vals[1:], 0.0)), vecs[:, 1:])


def stress_shapes(model: PlateModel, r: NDArray[np.float64]) -> NDArray[np.float64]:
    """Surface stress [Pa] per unit modal coordinate, `(len(r), modes)`: the larger of the
    radial and hoop bending stress."""
    n = model.modes.shape[0]
    b = [_basis(i, r) for i in range(n)]
    dd = np.array([bi[1] for bi in b]).T @ model.modes
    d1 = np.array([bi[2] for bi in b]).T @ model.modes
    sec = model.section
    m_r = sec.rigidity * (dd + STEEL_NU * d1)
    m_t = sec.rigidity * (d1 + STEEL_NU * dd)
    out: NDArray[np.float64] = np.maximum(np.abs(m_r), np.abs(m_t)) * sec.stress_per_moment
    return out


@dataclass(frozen=True)
class Kick:
    ringing_share: float  # of the kick's kinetic energy, in the bending modes
    ringing: float  # [J]
    peak: float  # [Pa], undamped, over one pulse interval
    bound: float  # [Pa], every mode at its peak together


def kick(model: PlateModel, sheet: float) -> Kick:
    """The plate's bending response to one pulse of the sheet `sheet`, applied over PULSE_TIME."""
    sec = model.section
    n = model.modes.shape[0]
    v0 = j_at(_R, sheet) / sec.areal
    phi = np.array([_basis(i, _R)[0] for i in range(n)])
    qdot = model.modes.T @ (phi @ (_W * sec.areal * v0))  # the rim starts at rest
    ke_total = 0.5 * float(np.sum(_W * sec.areal * v0**2))
    ke_ring = 0.5 * float(np.sum(qdot**2))
    x = 0.5 * model.omega * PULSE_TIME
    amp = qdot / model.omega * np.abs(np.sinc(x / math.pi))
    r_out = np.linspace(0.0, PLATE_R, 191)
    shapes = stress_shapes(model, r_out)
    keep = model.omega * PULSE_TIME < 20.0  # the sinc has cut the rest below 10%
    w_top = float(model.omega[keep].max())
    step = 0.1 * math.pi / w_top
    ts = np.arange(0.0, min(1.0 / CADENCE, 20000 * step), step)
    sig = np.abs(shapes[:, keep] @ (amp[keep, None] * np.sin(np.outer(model.omega[keep], ts))))
    return Kick(
        ke_ring / ke_total,
        ke_ring,
        float(sig.max()),
        float((shapes[:, keep] @ np.abs(amp[keep])).max()),
    )


def static_stress(model: PlateModel, force: float, ring: tuple[float, float] | None) -> float:
    """[Pa] Peak bending stress when `force` pushes on the plate through bags, against the
    plate's own inertia. `ring = (r1, r2)` puts the bags under that annulus; `None` spreads them
    under the whole floor, which bends it by the rim's share only."""
    sec = model.section
    n = model.modes.shape[0]
    floor = sec.areal * math.pi * PLATE_R**2
    acc = force / (floor + sec.rim)
    if ring is None:
        bag = np.full_like(_R, force / (math.pi * PLATE_R**2))
    else:
        r1, r2 = ring
        inside = (r1 <= _R) & (r2 >= _R)
        bag = np.where(inside, force / (math.pi * (r2**2 - r1**2)), 0.0)
    phi = np.array([_basis(i, _R)[0] for i in range(n)])
    load = phi @ (_W * (bag - sec.areal * acc)) - sec.rim * acc * _edge(n)
    q = (model.modes.T @ load) / model.omega**2
    r_out = np.linspace(0.0, PLATE_R, 191)
    return float(np.abs(stress_shapes(model, r_out) @ q).max())


#: `(section, sheets)`: ADR-0064's plate and an even 63 t plate under the flown 8 m sheet, then
#: the remedy -- a sheet as wide as the floor, with the skirt on the frame -- and its jitter.
CASES = (
    (Section("12.2 t solid, skirt on floor (ADR-0064)", 12.2e3, 0.0, True), (SHEET_FLOWN,)),
    (Section("63 t solid, skirt on floor", PLATE, 0.0, True), (SHEET_FLOWN,)),
    (Section("63 t sandwich 0.3 m, skirt on floor", PLATE, 0.3, True), (SHEET_FLOWN, 9.5)),
    (
        Section("63 t sandwich 0.3 m, skirt on frame", PLATE, 0.3, False),
        (SHEET_FLOWN, 9.0, 9.5, 10.0),
    ),
    (Section("63 t sandwich 1.0 m, skirt on frame", PLATE, 1.0, False), (9.0, 9.5, 10.0)),
)


@dataclass(frozen=True)
class RingRow:
    section: str
    sheet_m: float
    kapton_um: float  # the pulse's 28 kg of Kapton spread over the sheet
    kapton_holds_u8: bool  # Result 13's rule at the 9.86 R_sun headline, u = 8 km/s
    capture: float  # J/mv against the flown 8 m sheet's
    face_coeff: float  # P_peak / (sigma v^2) on that floor
    face_mm: float
    f1_hz: float
    rigid_dv: float  # [m/s], J / m
    ringing_share: float
    ringing_mj: float  # per pulse
    ringing_mw: float  # at 2 Hz, if each pulse's ringing is damped out before the next
    peak_mpa: float
    peak_coarse_mpa: float  # on two-thirds of the basis: the convergence check
    bound_mpa: float
    buildup_lo: float  # steady state over many pulses at ZETAS (coherent worst case)
    buildup_hi: float
    static_floor_mpa: float  # bags under the whole floor, at the cycle's peak force
    static_ring_mpa: float  # bags under a 5-8 m ring


def ring_rows() -> list[RingRow]:
    peak_force = column_for_stroke(PLATE, STROKE, 1.4).f_max
    headline = fast_wave.CLOSURES[0]
    sheet_mass = PULSE / (headline.k * headline.v_e)
    out = []
    for sec, sheets in CASES:
        model = plate_model(sec)
        coarse = plate_model(sec, 2 * N_BASIS // 3)
        build = [1.0 / (1.0 - math.exp(-z * float(model.omega[0]) / CADENCE)) for z in ZETAS]
        for s in sheets:
            k = kick(model, s)
            kapton = sheet_mass / (math.pi * s**2) / coated_sheet.KAPTON_RHO
            out.append(
                RingRow(
                    sec.label,
                    s,
                    1e6 * kapton,
                    fast_wave.kapton_holds(headline, 8.0e3, kapton),
                    capture(s) / capture(SHEET_FLOWN),
                    face_coeff(s),
                    sec.face_mm,
                    float(model.omega[0]) / (2.0 * math.pi),
                    PULSE / sec.plate,
                    k.ringing_share,
                    k.ringing / 1e6,
                    k.ringing * CADENCE / 1e6,
                    k.peak / 1e6,
                    kick(coarse, s).peak / 1e6,
                    k.bound / 1e6,
                    build[0],
                    build[1],
                    static_stress(model, peak_force, None) / 1e6,
                    static_stress(model, peak_force, (5.0, 8.0)) / 1e6,
                )
            )
    return out


# --- Output ---------------------------------------------------------------------------------

OUT_DESIGN = charge.OUT_DIR / "gas_spring_design.csv"
OUT_SPRING = charge.OUT_DIR / "gas_spring_loss.csv"
OUT_RING = charge.OUT_DIR / "gas_spring_ring.csv"


def _write(path: Path, rows: list[object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    dicts = [asdict(r) for r in rows]  # type: ignore[call-overload]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(dicts[0]))
        w.writeheader()
        w.writerows(dicts)


def main() -> None:
    print(
        f"  the paper's law: m s = {PLATE * linear_stroke(PLATE, THETA_PARENT) / 1e3:.0f} t m"
        f" at theta 2 pi / 3 -> {linear_stroke(PLATE, THETA_PARENT):.2f} m at {PLATE / 1e3:.0f} t"
    )
    des = design_rows()
    for d in des:
        print(
            f"  design {d.plate_t:5.1f} t theta {d.theta_lin:.3f}: stroke linear"
            f" {d.stroke_linear_m:5.2f} m, solved {d.stroke_m:5.2f} m | column {d.column_m:5.1f} m"
            f" | peak/mean {d.peak_over_mean:.2f}, P ratio {d.pressure_ratio:.2f}"
            f" | P V {d.gas_and_fibre_pv_mj:4.0f} MJ"
        )
    spr = spring_rows(column_for_stroke(PLATE, STROKE, 1.4))
    for s in spr:
        print(
            f"  {s.bag:<26} column {s.column_m:4.1f} m, theta {s.theta_lin:.2f},"
            f" peak/mean {s.peak_over_mean:.2f}, P x{s.pressure_ratio:.2f} | {s.footprint_m2:5.1f}"
            f" m^2 ({s.cell_radius_m:.2f} m cells) | gas {s.gas_t:.2f} t, fibre {s.fibre_t:.2f} t,"
            f" bladder {s.bladder_t:.2f} t"
        )
        print(
            f"      h lam/turb/K&S {s.h_laminar:5.1f}/{s.h_turbulent:5.1f}/{s.h_kornhauser:6.1f}"
            f" W/m^2/K | per cycle: thermal {s.thermal_mj:.3f}"
            f" (laminar {s.thermal_laminar_mj:.4f}),"
            f" fibre {s.fibre_mj_lo:.3f}-{s.fibre_mj_hi:.3f}, friction {s.friction_mj:.4f} MJ"
            f" -> {s.total_mj:.3f} MJ = {s.loss_share:.2%} of the {s.swing_mj:.0f} MJ swing"
            f" | push {s.push_gj:.2f} GJ: +{s.rise_k_no_cooling:.0f} K bare,"
            f" +{s.rise_k_radiating:.0f} K radiating | gas P x{s.gas_pressure_rise:.2f}"
        )
    ring = ring_rows()
    print(f"  the pulse loads the face for {1e6 * PULSE_TIME:.0f} us")
    for rr in ring:
        print(
            f"  ring {rr.section:<40} sheet {rr.sheet_m:4.1f} m, {rr.kapton_um:3.0f} um"
            f" ({'holds' if rr.kapton_holds_u8 else 'fails'} at u 8) (capture x{rr.capture:.3f},"
            f" face {rr.face_coeff:.2f} sigma v^2) | {rr.face_mm:4.1f} mm, f1 {rr.f1_hz:5.1f} Hz,"
            f" kick {rr.rigid_dv:4.0f} m/s | ringing {rr.ringing_share:5.1%} = {rr.ringing_mj:5.1f}"
            f" MJ ({rr.ringing_mw:4.0f} MW) | peak {rr.peak_mpa:6.0f} MPa (coarser"
            f" {rr.peak_coarse_mpa:6.0f}, bound {rr.bound_mpa:6.0f}) x{rr.buildup_lo:.2f}-"
            f"{rr.buildup_hi:.2f} | static: floor bags"
            f" {rr.static_floor_mpa:4.0f}, 5-8 m ring {rr.static_ring_mpa:5.0f} MPa"
        )
    _write(OUT_DESIGN, list(des))
    _write(OUT_SPRING, list(spr))
    _write(OUT_RING, list(ring))
    print(f"python: wrote {OUT_DESIGN.parent}/gas_spring_*.csv")


if __name__ == "__main__":
    main()
