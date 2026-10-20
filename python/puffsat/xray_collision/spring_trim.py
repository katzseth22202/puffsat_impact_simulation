"""Result 19: keeping the lossless spring on its cycle -- the pulse scatter, the drift and the trim.

    PYTHONPATH=python uv run python -m puffsat.xray_collision.spring_trim

The parent's S48 (`docs/spring_trim_asks_for_impact_sim.md`, commit 94f8808): Result 17's gas
spring has almost no loss, so it never forgets an impulse error, and the 63 t plate's swing
random-walks off its stroke. Two levers are taken (Seth, 2026-10-20). The first is the **next
charge's k**: the craft deploys each polyethylene charge, so it can size it from the swing it has
just measured. The second is the craft's **conventional thrusters**, which take up whatever is left
of the gap between the push it expected and the push it got. That budget is small and not yet
charged. Four things are computed here.

1. **The slopes** (`sweep --spring-trim-2d`). These are 2-D runs of the floor-wide 9.5 m sheet on
   the flown plate (ADR-0068). They give the impulse at k = 5-9 for a fixed sheet mass, and at a
   standoff of 1.5-3 m. Result 17's sheet-radius runs give the loss from a sheet that overhangs.
2. **The scatter budget.** Each source's assumed 1-sigma error, times its slope. There are two
   outputs: the heave, which is the fraction of the pulse's impulse, and the rocking. Rocking is
   the swing an off-centre impulse puts into the plate's edge, as a fraction of the heave kick.
3. **The drift, flown.** Result 17's nonlinear adiabatic column, `F = F0 (H / (H - x))^gamma` on
   the reduced mass, kicked every 0.5 s for 2000 pulses with a random impulse error. Its loss is a
   viscous drag sized to Result 17's 0.14-0.37% of the swing a cycle (`gas_spring_loss.csv`). The
   equivalence is an assumption, since the real loss is thermal hysteresis. The integration is RK4,
   vectorized over the random seeds.
4. **The trim.** A discrete LQR on the orbit's one-period map (the monodromy, from finite
   differences of the same RK4 flow). Its state is the swing's position and velocity just after a
   kick. It can also carry the summed impulse error, so the k trim holds the craft's delta-v too
   and the thrusters take only the remainder. Its input is the next charge's impulse, clipped to
   the k range the 2-D runs cover. It works with a measurement error and with a lag of 0-2 pulses
   between committing a charge and its pulse. Rocking cannot be trimmed by k. It is given the same
   controller acting on a charge with a slight wedge, which shifts the impulse's centroid.

The spring is flown at theta = 2 pi / 3 (the flown 18.6 m column) and at theta = pi (8.3 m), as
S48 item 4 asks. The pulse's slow taper as propellant burns (a programmed `J(t)`) is not modelled:
it moves the steady orbit slowly and predictably, which a feed-forward on k follows.
"""

from __future__ import annotations

import csv
import math
from dataclasses import asdict, dataclass
from functools import cache
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.linalg import solve_discrete_are
from scipy.optimize import brentq, fsolve

from puffsat.xray_collision import charge, gas_spring

G0 = 9.80665
PULSE = gas_spring.PULSE  # [N s]
CADENCE = gas_spring.CADENCE  # [Hz]
PERIOD = 1.0 / CADENCE
PLATE = gas_spring.PLATE  # [kg]
CRAFT = gas_spring.CRAFT  # [kg], the unit at the push's end
PLATE_R = gas_spring.PLATE_R  # [m]
CYCLES = gas_spring.CYCLES
K0 = gas_spring.SHEET_K  # 7.07, the 9.86 R_sun headline's push
W = 99.2e3  # [m/s], that push's closing speed (fast_wave.CLOSURES[0])
STANDOFF = 2.0  # [m]
SHEET_R = 9.5  # [m], ADR-0068's floor-wide sheet
#: The flown 18.6 m column, the stiffest of Result 17's columns that stays stable (11.9 m) and
#: the parent's gas-column fix (8.3 m).
THETAS = {"2pi/3": 2.0 * math.pi / 3.0, "5pi/6": 5.0 * math.pi / 6.0, "pi": math.pi}

JSONL = charge.OUT_DIR / "spring_trim_2d.jsonl"
JSONL_COARSE = charge.OUT_DIR / "spring_trim_2d_coarse.jsonl"
LOSS_CSV = charge.OUT_DIR / "gas_spring_loss.csv"


# --- 1. The slopes from the 2-D runs --------------------------------------------------------


@cache
def _runs(grid: str = "x1") -> tuple[dict[str, object], ...]:
    return tuple(charge.load_jsonl(JSONL if grid == "x1" else JSONL_COARSE))


def floor_capture(row: dict[str, object]) -> float:
    """`J / (m v)` landing inside the floor's 9.5 m. The skirt hangs on the craft's frame
    (ADR-0068), so its share does not kick the plate."""
    prof = row["profile"]
    assert isinstance(prof, list)
    return float(np.interp(PLATE_R, [float(p[0]) for p in prof], [float(p[1]) for p in prof]))


def ladder(grid: str = "x1") -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """`(k, J/(m v))` at the flown 2 m standoff, sorted by k."""
    rows = sorted(
        (float(str(r["k"])), floor_capture(r))
        for r in _runs(grid)
        if float(str(r["height"])) == STANDOFF
    )
    return np.array([r[0] for r in rows]), np.array([r[1] for r in rows])


def capture0(grid: str = "x1") -> float:
    ks, cs = ladder(grid)
    return float(np.interp(K0, ks, cs))


def slope_k(grid: str = "x1") -> float:
    """`d ln J / d k` at K0 and fixed sheet mass: a straight-line fit to the ladder's runs at
    k = 6, 7.07 and 8."""
    ks, cs = ladder(grid)
    near = np.abs(ks - K0) <= 1.0
    a = np.polyfit(ks[near], cs[near], 1)[0]
    return float(a / capture0(grid))


def beta(grid: str = "x1") -> float:
    """`d ln J / d ln k` at fixed sheet mass: the charge's share of the impulse's response."""
    return K0 * slope_k(grid)


def slope_standoff(grid: str = "x1") -> float:
    """`d ln J / d h` [1/m] at K0: a straight-line fit over standoffs of 1.5-3 m."""
    rows = sorted(
        (float(str(r["height"])), floor_capture(r)) for r in _runs(grid) if float(str(r["k"])) == K0
    )
    hs = np.array([r[0] for r in rows])
    cs = np.array([r[1] for r in rows])
    return float(np.polyfit(hs, cs, 1)[0] / capture0(grid))


def overhang_loss() -> float:
    """The impulse a sheet loses per unit of its mass that lands outside the floor. It comes from
    Result 17's 10 m sheet against the 9.5 m sheet on the same floor (fixed sheet mass)."""
    lost = 1.0 - gas_spring.capture(10.0) / gas_spring.capture(9.5)
    outside = 1.0 - (9.5 / 10.0) ** 2
    return lost / outside


def sheet_mass() -> float:
    """[kg] the floor-wide sheet's mass for one 7.35 MN s pulse at K0 and W."""
    return PULSE / (capture0() * W)


def trim_range() -> tuple[float, float]:
    """The fractional impulse trim the 2-D ladder covers (k = 5 to 9)."""
    _, cs = ladder()
    c0 = capture0()
    return float(cs[0] / c0 - 1.0), float(cs[-1] / c0 - 1.0)


def k_for(u: NDArray[np.float64]) -> NDArray[np.float64]:
    """The charge's k that gives a fractional impulse trim `u`, interpolated on the ladder."""
    ks, cs = ladder()
    out: NDArray[np.float64] = np.interp(1.0 + u, cs / capture0(), ks)
    return out


@dataclass(frozen=True)
class SlopeRow:
    quantity: str
    value: float
    coarse: float
    unit: str


def slope_rows() -> list[SlopeRow]:
    lo, hi = trim_range()
    return [
        SlopeRow("J/(m v) at k 7.07, 2 m (floor only)", capture0(), capture0("x2"), "-"),
        SlopeRow("d ln J / d k (fixed sheet)", slope_k(), slope_k("x2"), "1/k"),
        SlopeRow("d ln J / d ln k (fixed sheet), beta", beta(), beta("x2"), "-"),
        SlopeRow("d ln J / d standoff", slope_standoff(), slope_standoff("x2"), "1/m"),
        SlopeRow("impulse lost per overhanging mass share", overhang_loss(), math.nan, "-"),
        SlopeRow("sheet mass", sheet_mass(), math.nan, "kg"),
        SlopeRow("charge mass at k 7.07", K0 * sheet_mass(), math.nan, "kg"),
        SlopeRow("trim floor (k = 5)", lo, math.nan, "of J"),
        SlopeRow("trim ceiling (k = 9)", hi, math.nan, "of J"),
    ]


# --- 2. The scatter budget ------------------------------------------------------------------

#: Assumed 1-sigma errors. None is measured. Each is a design tolerance the trim must cover.
SHEET_MASS_ERR = 0.01  # the sheet's mass on arrival: weighed at launch, eroded en route
CHARGE_MASS_ERR = 0.005  # the craft meters each charge
STANDOFF_ERR = 0.05  # [m], the charge's placement against where the plate will be
AIM_ERR = 0.10  # [m] per axis, the sheet's centre against the plate's
SPEED_ERR = 1.0e-4  # the sheet's closing speed (10 m/s at 99 km/s)
SHEET_AREAL_ERR = 0.05  # the film's areal density over each m^2 (gauge tolerance)
CHARGE_AREAL_ERR = 0.02  # the charge's areal density over each m^2
AREAL_CELL = 1.0  # [m^2], the scale those two vary on
TIMING_ERR = 1.0e-4  # [s], the sheet's arrival time (10 m along its track)
#: [m/s] the slowing the last pulse's plume gives the sheet (ADR-0069 item 6: 3-5 m/s). It is
#: predictable, so only its variation scatters; its whole size is taken as the bound.
PLUME_SLOWING = 5.0


def centroid_per_areal(share: float) -> float:
    """[m per unit areal error] the impulse centroid's 1-sigma shift on one axis. The cause is
    independent areal errors on `AREAL_CELL` patches, each of which moves the local impulse by
    `share` of its own error: `share sqrt(a / 4 pi)` for a uniform disk."""
    return share * math.sqrt(AREAL_CELL / (4.0 * math.pi))


def rock_share(offset: float) -> float:
    """The edge's kick as a fraction of the heave kick, for an impulse centroid `offset` off the
    axis. On a uniform floor, `I = m R^2 / 4`, so the edge moves at `4 (J / m) (offset / R)`."""
    return 4.0 * offset / PLATE_R


def aim_loss_stats(sigma: float) -> tuple[float, float]:
    """`(mean, std)` of the impulse a centred-in-error sheet loses off the floor's edge. The
    offset is Rayleigh with `sigma` per axis, and the lune that overhangs is `2 d / (pi R)` of
    the sheet. This is an upper bound, since the far side's gap gains a little (Result 17's 9 m
    sheet captures more than the 9.5 m)."""
    scale = overhang_loss() * 2.0 / (math.pi * SHEET_R)
    return scale * sigma * math.sqrt(math.pi / 2.0), scale * sigma * math.sqrt(2.0 - math.pi / 2.0)


@dataclass(frozen=True)
class ScatterRow:
    source: str
    error: float
    unit: str
    heave: float  # 1-sigma fraction of J
    rock: float  # 1-sigma edge kick per axis, fraction of the heave kick
    note: str


def scatter_rows() -> list[ScatterRow]:
    b = beta()
    mean_loss, std_loss = aim_loss_stats(AIM_ERR)
    rows = [
        ScatterRow(
            "sheet mass on arrival",
            SHEET_MASS_ERR,
            "of m",
            (1.0 - b) * SHEET_MASS_ERR,
            0.0,
            "fixed charge: k falls as m rises",
        ),
        ScatterRow(
            "charge mass", CHARGE_MASS_ERR, "of k m", b * CHARGE_MASS_ERR, 0.0, "metered by craft"
        ),
        ScatterRow(
            "standoff",
            STANDOFF_ERR,
            "m",
            abs(slope_standoff()) * STANDOFF_ERR,
            0.0,
            "2-D slope at 1.5-3 m",
        ),
        ScatterRow(
            "aim offset",
            AIM_ERR,
            "m per axis",
            std_loss,
            rock_share(AIM_ERR),
            f"also a mean loss {mean_loss:.2%}, taken out by the mean k",
        ),
        ScatterRow("closing speed", SPEED_ERR, "of w", SPEED_ERR, 0.0, "J ~ w at fixed k"),
        ScatterRow(
            "plume drag on the sheet",
            PLUME_SLOWING / W,
            "of w",
            PLUME_SLOWING / W,
            0.0,
            "ADR-0069: 3-5 m/s, predictable; bounded by its whole size",
        ),
        ScatterRow(
            "sheet areal density",
            SHEET_AREAL_ERR,
            "per m^2",
            0.0,
            rock_share(centroid_per_areal(1.0 - b) * SHEET_AREAL_ERR),
            "its total is the sheet-mass line",
        ),
        ScatterRow(
            "charge areal density",
            CHARGE_AREAL_ERR,
            "per m^2",
            0.0,
            rock_share(centroid_per_areal(b) * CHARGE_AREAL_ERR),
            "its total is the charge-mass line",
        ),
        ScatterRow(
            "arrival time",
            TIMING_ERR,
            "s",
            0.0,
            0.0,
            f"a phase error: {PULSE / PLATE * TIMING_ERR * 100:.1f} cm of swing",
        ),
    ]
    heave = math.sqrt(sum(r.heave**2 for r in rows))
    rock = math.sqrt(sum(r.rock**2 for r in rows))
    rows.append(ScatterRow("total (rss)", math.nan, "-", heave, rock, ""))
    return rows


def scatter_total() -> tuple[float, float]:
    t = scatter_rows()[-1]
    return t.heave, t.rock


# --- 3. The spring and its one-period map ---------------------------------------------------

STEPS = 64  # RK4 steps a period (hw ~ 0.05-0.1; tests hold the orbit to 1e-6 of gas_spring's)


def loss_shares() -> tuple[float, float]:
    """Result 17's least and largest loss a cycle, as a share of the swing's energy."""
    with LOSS_CSV.open() as fh:
        shares = [float(r["loss_share"]) for r in csv.DictReader(fh)]
    return min(shares), max(shares)


def _gamma() -> float:
    return gas_spring.gas("nitrogen", 52.0e3).gamma


@dataclass(frozen=True)
class Spring:
    label: str
    theta: float
    loss: float  # share of the swing's energy a cycle
    column: float  # [m]
    gamma: float
    f0: float  # [N]
    mu: float  # [kg]
    kick: float  # [m/s], J / m
    damp: float  # [1/s], the drag's c / mu
    x0: float  # [m], just after the kick, on the (damped) orbit
    v0: float
    x_lo: float
    x_hi: float
    a: tuple[tuple[float, float], tuple[float, float]]  # the flow's one-period Jacobian

    @property
    def mono(self) -> NDArray[np.float64]:
        return np.array(self.a)

    @property
    def stroke(self) -> float:
        return self.x_hi - self.x_lo


def _accel(
    x: NDArray[np.float64],
    v: NDArray[np.float64],
    column: NDArray[np.float64],
    damp: NDArray[np.float64],
    gamma: float,
    f0: float,
    mu: float,
) -> NDArray[np.float64]:
    out: NDArray[np.float64] = -f0 * (column / (column - x)) ** gamma / mu - damp * v
    return out


def flow(
    x: NDArray[np.float64],
    v: NDArray[np.float64],
    column: NDArray[np.float64],
    damp: NDArray[np.float64],
    gamma: float,
    f0: float,
    mu: float,
    record: bool = False,
) -> tuple[NDArray[np.float64], ...]:
    """One period of RK4. Returns `(x, v, x_min, x_max)`, plus `(xs, vs)` at each step when
    `record` is set."""
    h = PERIOD / STEPS
    lo, hi = x.copy(), x.copy()
    xs, vs = [], []
    for _ in range(STEPS):
        a1 = _accel(x, v, column, damp, gamma, f0, mu)
        x2, v2 = x + 0.5 * h * v, v + 0.5 * h * a1
        a2 = _accel(x2, v2, column, damp, gamma, f0, mu)
        x3, v3 = x + 0.5 * h * v2, v + 0.5 * h * a2
        a3 = _accel(x3, v3, column, damp, gamma, f0, mu)
        x4, v4 = x + h * v3, v + h * a3
        a4 = _accel(x4, v4, column, damp, gamma, f0, mu)
        x = x + h / 6.0 * (v + 2.0 * v2 + 2.0 * v3 + v4)
        v = v + h / 6.0 * (a1 + 2.0 * a2 + 2.0 * a3 + a4)
        lo, hi = np.minimum(lo, x), np.maximum(hi, x)
        if record:
            xs.append(x.copy())
            vs.append(v.copy())
    if record:
        return x, v, lo, hi, np.array(xs), np.array(vs)
    return x, v, lo, hi


def _one(
    z: NDArray[np.float64], column: float, damp: float, consts: tuple[float, float, float]
) -> tuple[NDArray[np.float64], ...]:
    g, f0, mu = consts
    return flow(np.array([z[0]]), np.array([z[1]]), np.array([column]), np.array([damp]), g, f0, mu)


def swing_energy() -> float:
    """Result 17's swing, `J^2 / 8 m`."""
    return PULSE**2 / (8.0 * PLATE)


@cache
def spring(label: str, loss: float) -> Spring:
    """The orbit at `THETAS[label]` with a drag that takes `loss` of the swing a cycle."""
    theta = THETAS[label]
    g = _gamma()
    o = gas_spring.column_for_theta(PLATE, theta, g)
    mu = gas_spring.reduced(PLATE)
    kick = PULSE / PLATE
    consts = (g, o.f0, mu)
    # The drag: c int v^2 dt over the undamped orbit = loss x the swing.
    rec = flow(
        np.array([o.x0]), np.array([o.v0]), np.array([o.column]), np.array([0.0]), *consts, True
    )
    vs = rec[5][:, 0]
    damp = loss * swing_energy() / (mu * float(np.sum(vs**2)) * PERIOD / STEPS)

    def residual(z: NDArray[np.float64]) -> NDArray[np.float64]:
        x, v, _, _ = _one(z, o.column, damp, consts)
        return np.array([x[0] - z[0], v[0] + kick - z[1]])

    z = fsolve(residual, np.array([o.x0, o.v0]), xtol=1e-10)
    _, _, lo, hi = _one(z, o.column, damp, consts)
    # The flow's Jacobian over one period, by central differences.
    jac = np.zeros((2, 2))
    for j, eps in enumerate((1e-5, 1e-5)):
        dz = np.zeros(2)
        dz[j] = eps
        p = _one(z + dz, o.column, damp, consts)
        m = _one(z - dz, o.column, damp, consts)
        jac[0, j] = (p[0][0] - m[0][0]) / (2 * eps)
        jac[1, j] = (p[1][0] - m[1][0]) / (2 * eps)
    a = ((float(jac[0, 0]), float(jac[0, 1])), (float(jac[1, 0]), float(jac[1, 1])))
    return Spring(
        label,
        theta,
        loss,
        o.column,
        g,
        o.f0,
        mu,
        kick,
        damp,
        float(z[0]),
        float(z[1]),
        float(lo[0]),
        float(hi[0]),
        a,
    )


@cache
def substep_maps(label: str, loss: float) -> NDArray[np.float64]:
    """`(STEPS, 2, 2)`: the variational map from just after a kick to each RK4 step. It is used
    for the rocking modes' excursions within a period."""
    sp = spring(label, loss)
    consts = (sp.gamma, sp.f0, sp.mu)
    out = np.zeros((STEPS, 2, 2))
    for j in range(2):
        dz = np.zeros(2)
        dz[j] = 1e-5
        args = (np.array([sp.column]), np.array([sp.damp]), *consts, True)
        p = flow(np.array([sp.x0 + dz[0]]), np.array([sp.v0 + dz[1]]), *args)
        m = flow(np.array([sp.x0 - dz[0]]), np.array([sp.v0 - dz[1]]), *args)
        out[:, 0, j] = (p[4][:, 0] - m[4][:, 0]) / 2e-5
        out[:, 1, j] = (p[5][:, 0] - m[5][:, 0]) / 2e-5
    return out


def rotation(sp: Spring) -> tuple[float, float]:
    """`(angle, modulus)` of the one-period map's eigenvalues: the swing's phase advance per
    pulse and its decay."""
    lam = np.linalg.eigvals(sp.mono)
    return float(abs(np.angle(lam[0]))), float(np.max(np.abs(lam)))


def variational_trace(theta: float) -> float:
    """The trace of the undamped orbit's one-period map at `theta`, from the exact variational
    equation `d'' = -(gamma F / (H - x)) d / mu` along the orbit (solve_ivp, independent of the
    RK4 above). A trace below -2 means a real eigenvalue below -1: the steady cycle is unstable
    and an error grows every pulse with alternating sign (Hill's first tongue, at a phase of pi)."""
    g = _gamma()
    o = gas_spring.column_for_theta(PLATE, theta, g)
    mu = gas_spring.reduced(PLATE)

    def rhs(t: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        x, v, a, b, c, d = y
        f = o.f0 * (o.column / (o.column - x)) ** g
        k = g * f / (o.column - x)
        return np.array([v, -f / mu, b, -k / mu * a, d, -k / mu * c])

    y = solve_ivp(rhs, (0.0, PERIOD), [o.x0, o.v0, 1.0, 0.0, 0.0, 1.0], rtol=1e-11, atol=1e-12)
    end = y.y[:, -1]
    return float(end[2] + end[5])


def stability_edge() -> float:
    """The linearized `theta` at which the nonlinear cycle goes unstable (trace = -2)."""
    return float(brentq(lambda t: variational_trace(t) + 2.0, 0.75 * math.pi, math.pi, xtol=1e-6))


def two_pulse_phase_trim(sp: Spring, error: float = 0.1) -> float:
    """The largest of the two impulse trims (fraction of J) that cancel a pure position error of
    `error` m, deadbeat, in two pulses. It measures how well the next charge reaches the swing's
    phase."""
    a = sp.mono
    b = np.array([0.0, sp.kick])
    # z2 = A^2 e + A b u1 + b u2 = 0, with e the post-kick error.
    reach = np.column_stack([a @ b, b])
    u = np.linalg.solve(reach, -(a @ a) @ np.array([error, 0.0]))
    return float(np.max(np.abs(u)))


# --- 4. Flying the push ---------------------------------------------------------------------

#: The swing's measurement: laser ranging of the plate over each 0.5 s interval (assumed).
MEAS_X = 0.01  # [m]
MEAS_V = 0.01  # [m/s]
#: LQR weights (chosen): position in units of 0.5 m, velocity of omega x 0.5 m, the summed
#: impulse in units of 0.2 pulse, the trim in units of 5% of a pulse.
Q_X = 0.5
Q_S = 0.2
R_U = 0.05
SEEDS = 100
#: [m] an excursion past either end of the steady stroke that ends the push: the bags would
#: have to be built half again as long (chosen).
LOST = 4.0
#: Thruster specific impulse (assumed, hydrazine class) for the uncharged trim budget.
THRUSTER_ISP = 300.0  # [s]
#: The rocking trim: the charge's areal density may tilt up to this fraction across its radius.
WEDGE_MAX = 0.30


def gains(sp: Spring, integrate: bool) -> NDArray[np.float64]:
    """LQR on `z' = A z + B u`, with the summed impulse error as a third state if `integrate`."""
    omega = sp.theta * CADENCE
    a = np.zeros((3, 3))
    a[:2, :2] = sp.mono
    a[2, 2] = 1.0
    b = np.array([[0.0], [sp.kick], [1.0]])
    q = np.diag([1.0 / Q_X**2, 1.0 / (omega * Q_X) ** 2, 1.0 / Q_S**2])
    r = np.array([[1.0 / R_U**2]])
    if not integrate:
        a, b, q = a[:2, :2], b[:2], q[:2, :2]
    p = solve_discrete_are(a, b, q, r)
    k = np.linalg.solve(r + b.T @ p @ b, b.T @ p @ a)[0]
    return np.concatenate([k, [0.0]]) if not integrate else k


@dataclass(frozen=True)
class Case:
    spring: str
    loss: float
    scatter: float
    timing: float = 0.0  # [s] 1-sigma arrival error
    trim: bool = False
    integrate: bool = True
    latency: int = 0  # extra pulses between committing a charge and its pulse


@dataclass(frozen=True)
class Flight:
    excursion_p50: float  # [m] beyond the steady orbit's range, over the push
    excursion_p95: float
    first_metre_p50: float  # pulses to a 1 m excursion (inf if never)
    lost_p50: float  # pulses to an excursion of LOST (inf if never)
    lost_share: float  # of seeds, by the push's end
    final_dv_rms: float  # [m/s] the swing's velocity error at the end
    trim_rms: float  # fraction of J
    k_rms: float
    k_lo_p01: float
    k_hi_p99: float
    saturated: float  # share of pulses at the ladder's edge
    craft_dv_p95: float  # [m/s] left for the thrusters at the push's end


def fly(cases: list[Case], pulses: int = CYCLES, seed: int = 48) -> list[Flight]:
    """All `cases` together, `SEEDS` seeds each, on the nonlinear spring. Every case in a call
    must share its latency."""
    lat = {c.latency for c in cases}
    assert len(lat) == 1
    latency = lat.pop()
    rng = np.random.default_rng(seed)
    springs = [spring(c.spring, c.loss) for c in cases]
    rep = np.repeat
    n = len(cases) * SEEDS
    col = rep([s.column for s in springs], SEEDS)
    damp = rep([s.damp for s in springs], SEEDS)
    x_ref, v_ref = rep([s.x0 for s in springs], SEEDS), rep([s.v0 for s in springs], SEEDS)
    x_lo, x_hi = rep([s.x_lo for s in springs], SEEDS), rep([s.x_hi for s in springs], SEEDS)
    sig = rep([c.scatter for c in cases], SEEDS)
    sig_t = rep([c.timing for c in cases], SEEDS)
    on = rep([c.trim for c in cases], SEEDS)
    kk = np.repeat(
        np.array([gains(s, c.integrate) for s, c in zip(springs, cases, strict=True)]), SEEDS, 0
    )
    kk[~on] = 0.0
    aa = np.zeros((n, 3, 3))
    aa[:, :2, :2] = np.repeat(np.array([s.mono for s in springs]), SEEDS, 0)
    aa[:, 2, 2] = 1.0
    sp0 = springs[0]
    gamma, f0, mu, kick = sp0.gamma, sp0.f0, sp0.mu, sp0.kick
    bb = np.array([0.0, kick, 1.0])
    u_lo, u_hi = trim_range()

    x, v = x_ref.copy(), v_ref.copy()
    s_true = np.zeros(n)
    s_meas = np.zeros(n)
    queue = [np.zeros(n) for _ in range(latency)]
    lo_all, hi_all = x.copy(), x.copy()
    first = np.full(n, np.inf)
    lost_at = np.full(n, np.inf)
    trims = np.zeros((pulses, n))
    u_next = np.zeros(n)
    for i in range(pulses):
        x, v, lo, hi = flow(x, v, col, damp, gamma, f0, mu)
        lo_all, hi_all = np.minimum(lo_all, lo), np.maximum(hi_all, hi)
        exc = np.maximum(np.maximum(hi_all - x_hi, x_lo - lo_all), 0.0)
        first = np.where((exc >= 1.0) & np.isinf(first), i + 1, first)
        # A seed that leaves the stroke by LOST is out of the push. It is put back on the orbit
        # so the batch stays finite, and its excursion is held at infinity.
        gone = ~np.isfinite(exc) | (exc >= LOST)
        lost_at = np.where(gone & np.isinf(lost_at), i + 1, lost_at)
        x, v = np.where(gone, x_ref, x), np.where(gone, v_ref, v)
        lo_all, hi_all = np.where(gone, x_ref, lo_all), np.where(gone, x_ref, hi_all)
        u = queue.pop(0) if latency else u_next
        trims[i] = u
        w = rng.normal(0.0, 1.0, n) * sig
        v = v + kick * (1.0 + w + u)
        x = x - kick * rng.normal(0.0, 1.0, n) * sig_t
        s_true = s_true + w + u
        s_meas = s_meas + w + u + rng.normal(0.0, 1.0, n) * (math.sqrt(2.0) * MEAS_V / kick)
        z = np.stack(
            [
                x - x_ref + rng.normal(0.0, 1.0, n) * MEAS_X,
                v - v_ref + rng.normal(0.0, 1.0, n) * MEAS_V,
                s_meas,
            ],
            axis=1,
        )
        for q in queue:
            z = np.einsum("nij,nj->ni", aa, z) + bb * q[:, None]
        u_new = np.clip(-np.einsum("ni,ni->n", kk, z), u_lo, u_hi)
        if latency:
            queue.append(u_new)
        else:
            u_next = u_new
    exc = np.maximum(np.maximum(hi_all - x_hi, x_lo - lo_all), 0.0)
    exc = np.where(np.isfinite(lost_at), np.inf, exc)
    ks = k_for(trims)
    out = []
    for j, c in enumerate(cases):
        sl = slice(j * SEEDS, (j + 1) * SEEDS)
        t = trims[:, sl]
        kj = ks[:, sl]
        sat = float(np.mean((t <= u_lo + 1e-12) | (t >= u_hi - 1e-12))) if c.trim else 0.0
        alive = np.isinf(lost_at[sl])
        dv = (v[sl] - v_ref[sl])[alive]
        out.append(
            Flight(
                float(np.quantile(exc[sl], 0.50, method="higher")),
                float(np.quantile(exc[sl], 0.95, method="higher")),
                float(np.quantile(first[sl], 0.5, method="higher")),
                float(np.quantile(lost_at[sl], 0.5, method="higher")),
                float(np.mean(~alive)),
                float(np.sqrt(np.mean(dv**2))) if dv.size else math.nan,
                float(np.sqrt(np.mean(t**2))),
                float(np.sqrt(np.mean((kj - K0) ** 2))) if c.trim else 0.0,
                float(np.percentile(kj, 1)) if c.trim else K0,
                float(np.percentile(kj, 99)) if c.trim else K0,
                sat,
                float(np.percentile(np.abs(s_true[sl]), 95) * PULSE / CRAFT),
            )
        )
    return out


def lost_linear(label: str, loss: float, scatter: float, seed: int = 4849) -> float:
    """Median pulses to an excursion of LOST on the linearized map. It is the same walk as `fly`
    without the column's nonlinearity, and shows how much of the drift's pace that adds."""
    sp = spring(label, loss)
    sub = substep_maps(label, loss)
    rng = np.random.default_rng(seed)
    n = 4 * SEEDS
    z = np.zeros((n, 2))
    worst = np.zeros(n)
    first = np.full(n, np.inf)
    for i in range(CYCLES):
        worst = np.maximum(worst, np.abs(sub[:, 0, :] @ z.T).max(axis=0))
        first = np.where((worst >= LOST) & np.isinf(first), i + 1, first)
        z = z @ sp.mono.T
        z[:, 1] += sp.kick * rng.normal(0.0, scatter, n)
    return float(np.quantile(first, 0.5, method="higher"))


def thruster_propellant(impulse: float) -> float:
    """[kg] at `THRUSTER_ISP`."""
    return impulse / (G0 * THRUSTER_ISP)


@dataclass(frozen=True)
class DriftRow:
    spring: str
    column_m: float
    stroke_m: float
    loss_share: float
    scatter: float
    excursion_p50_m: float
    excursion_p95_m: float
    pulses_to_1m_p50: float
    pulses_to_lost_p50: float  # an excursion of LOST m
    pulses_to_lost_linear: float  # the same on the linearized map
    lost_share: float
    swing_error_rms: float  # [m/s] at the end, seeds still flying
    swing_error_estimate: float  # the parent's s J / m / sqrt(d) (inf when lossless)
    pulses_parent: float  # the parent's: a velocity error of half the kick, (1 / 2 s)^2


SCATTERS = (0.01, 0.02, 0.05, 0.10)


def drift_rows() -> list[DriftRow]:
    lo, hi = loss_shares()
    cases = [Case(lab, d, s) for lab in THETAS for d in (0.0, lo, hi) for s in SCATTERS]
    flights = fly(cases)
    out = []
    for c, f in zip(cases, flights, strict=True):
        sp = spring(c.spring, c.loss)
        est = c.scatter * sp.kick / math.sqrt(c.loss) if c.loss > 0 else math.inf
        out.append(
            DriftRow(
                c.spring,
                sp.column,
                sp.stroke,
                c.loss,
                c.scatter,
                f.excursion_p50,
                f.excursion_p95,
                f.first_metre_p50,
                f.lost_p50,
                lost_linear(c.spring, c.loss, c.scatter),
                f.lost_share,
                f.final_dv_rms,
                est,
                (0.5 / c.scatter) ** 2,
            )
        )
    return out


@dataclass(frozen=True)
class TrimRow:
    spring: str
    scatter: float
    timing_ms: float
    integrate: bool
    latency: int
    excursion_p50_m: float
    excursion_p95_m: float
    trim_rms: float
    k_rms: float
    k_p01: float
    k_p99: float
    charge_rms_kg: float
    saturated: float
    lost_share: float
    craft_dv_p95: float  # [m/s] for the thrusters
    thruster_kg: float  # at THRUSTER_ISP on the end-of-push craft


def trim_cases() -> list[Case]:
    """The budget's scatter and the parent's 2-10%, on the least loss (the hardest case)."""
    lo, _ = loss_shares()
    two, five, pi = list(THETAS)
    budget = round(scatter_total()[0], 4)
    cases = [Case(two, lo, s, trim=True) for s in (budget, 0.02, 0.05, 0.10)]
    cases += [Case(two, lo, budget, trim=True, integrate=False)]
    cases += [Case(two, lo, budget, TIMING_ERR, trim=True)]
    cases += [Case(two, lo, s, trim=True, latency=lag) for s in (budget, 0.02) for lag in (1, 2)]
    cases += [Case(five, lo, s, trim=True) for s in (budget, 0.02)]
    cases += [Case(five, lo, budget, trim=True, latency=2)]
    cases += [Case(pi, lo, s, trim=True) for s in (budget, 0.02, 0.05)]
    cases += [Case(pi, lo, budget, TIMING_ERR, trim=True)]
    cases += [Case(pi, lo, budget, trim=True, latency=lag) for lag in (1, 2)]
    return cases


def trim_rows() -> list[TrimRow]:
    cases = trim_cases()
    flights: dict[Case, Flight] = {}
    for latency in sorted({c.latency for c in cases}):
        group = [c for c in cases if c.latency == latency]
        flights.update(zip(group, fly(group), strict=True))
    m = sheet_mass()
    out = []
    for c in cases:
        f = flights[c]
        out.append(
            TrimRow(
                c.spring,
                c.scatter,
                c.timing * 1e3,
                c.integrate,
                c.latency,
                f.excursion_p50,
                f.excursion_p95,
                f.trim_rms,
                f.k_rms,
                f.k_lo_p01,
                f.k_hi_p99,
                f.k_rms * m,
                f.saturated,
                f.lost_share,
                f.craft_dv_p95,
                thruster_propellant(f.craft_dv_p95 * CRAFT),
            )
        )
    return out


@dataclass(frozen=True)
class LeverRow:
    """What each lever would spend to hold the spring alone, over the push."""

    lever: str
    scatter: float
    impulse_mns: float  # [MN s] summed over the push
    propellant_kg: float
    note: str


def lever_rows() -> list[LeverRow]:
    out = []
    mean_abs = math.sqrt(2.0 / math.pi)
    for s in (0.01, 0.05):
        imp = CYCLES * s * PULSE * mean_abs
        out.append(
            LeverRow(
                "thrusters on the plate cancel each pulse's error",
                s,
                imp / 1e6,
                thruster_propellant(imp),
                "the whole of each error, every pulse",
            )
        )
        net = math.sqrt(CYCLES) * s * PULSE
        out.append(
            LeverRow(
                "thrusters take only the summed push error",
                s,
                net / 1e6,
                thruster_propellant(net),
                "1-sigma random walk; the spring is not held",
            )
        )
    _, rock = scatter_total()
    # Per axis the angular impulse is J x offset = J x (share R / 4); a couple at the rim pays it
    # as that over R in thruster impulse.
    imp = 2.0 * CYCLES * rock * PULSE / 4.0 * mean_abs
    out.append(
        LeverRow(
            "thruster couples at the rim cancel each pulse's rocking",
            rock,
            imp / 1e6,
            thruster_propellant(imp),
            "both axes; angular impulse J x offset / R",
        )
    )
    # Alignment guides (the parent's absorber budget carries them) that hold the plate's tilt
    # hand each pulse's angular impulse J x offset to the craft. The craft has no spring, so only
    # the sum matters: a random walk, cleared by couples at an arm of PLATE_R (assumed).
    net = 2.0 * math.sqrt(CYCLES) * PULSE * AIM_ERR / PLATE_R
    out.append(
        LeverRow(
            "tilt guides; thrusters take only the summed angular impulse",
            AIM_ERR,
            net / 1e6,
            thruster_propellant(net),
            "aim error in m per axis; 1-sigma walk, both axes, couples at R",
        )
    )
    return out


# --- 5. Rocking -----------------------------------------------------------------------------


@dataclass(frozen=True)
class RockRow:
    spring: str
    aim_m: float
    rock_share: float  # 1-sigma edge kick per axis, of the heave kick
    trimmed: bool
    latency: int
    edge_p50_m: float  # the edge's excursion over the push, both axes
    edge_p95_m: float
    wedge_rms: float  # the charge's areal tilt across its radius
    wedge_p99: float
    lost_share: float


def rock_share_at(aim: float) -> float:
    """The rss edge share with the aim error set to `aim` and the areal terms as budgeted."""
    rows = [r for r in scatter_rows()[:-1] if r.source != "aim offset"]
    return math.sqrt(sum(r.rock**2 for r in rows) + rock_share(aim) ** 2)


def fly_rocking(
    label: str, loss: float, share: float, trim: bool, latency: int, seed: int = 4848
) -> tuple[float, float, float, float, float]:
    """Both rocking axes, linearized about the heave orbit. On a uniform floor over uniform bags,
    rocking obeys the heave's variational equation. Edge kicks are `share x kick`, and the trim
    is a wedge on the charge: `beta x wedge` of the kick at the edge."""
    sp = spring(label, loss)
    a = sp.mono
    sub = substep_maps(label, loss)
    b = np.array([0.0, sp.kick])
    k = gains(sp, integrate=False)[:2]
    rng = np.random.default_rng(seed)
    n = 2 * SEEDS * 4  # two axes per seed
    z = np.zeros((n, 2))
    worst = np.zeros(n // 2)
    queue = [np.zeros(n) for _ in range(latency)]
    u_max = beta() * WEDGE_MAX
    u_next = np.zeros(n)
    wedges = []
    lost = np.zeros(n // 2, dtype=bool)
    for _ in range(CYCLES):
        xs = sub[:, 0, :] @ z.T  # (STEPS, n): each axis's edge displacement this period
        edge = np.sqrt(xs[:, : n // 2] ** 2 + xs[:, n // 2 :] ** 2).max(axis=0)
        worst = np.maximum(worst, edge)
        # As in `fly`: a seed whose edge leaves the stroke by LOST is out, and is reset.
        lost |= worst >= LOST
        both = np.concatenate([lost, lost])
        z[both] = 0.0
        z = z @ a.T
        u = queue.pop(0) if latency else u_next
        wedges.append(u / beta())
        z = z + np.outer(rng.normal(0.0, share, n) + u, b)
        zm = z + rng.normal(0.0, 1.0, (n, 2)) * np.array([MEAS_X, MEAS_V])
        for q in queue:
            zm = zm @ a.T + np.outer(q, b)
        u_new = np.clip(-(zm @ k), -u_max, u_max) if trim else np.zeros(n)
        if latency:
            queue.append(u_new)
        else:
            u_next = u_new
    wd = np.abs(np.array(wedges))
    worst = np.where(lost, np.inf, worst)
    return (
        float(np.quantile(worst, 0.50, method="higher")),
        float(np.quantile(worst, 0.95, method="higher")),
        float(np.sqrt(np.mean(wd**2))),
        float(np.percentile(wd, 99)),
        float(np.mean(lost)),
    )


AIMS = (0.01, 0.03, 0.10, 0.30)


def rock_rows() -> list[RockRow]:
    lo, _ = loss_shares()
    out = []
    for lab in THETAS:
        for aim in AIMS:
            share = rock_share_at(aim)
            for trim, lat in ((False, 0), (True, 0), (True, 2)):
                p50, p95, wr, w99, gone = fly_rocking(lab, lo, share, trim, lat)
                out.append(RockRow(lab, aim, share, trim, lat, p50, p95, wr, w99, gone))
    return out


# --- 6. The spring's map --------------------------------------------------------------------


@dataclass(frozen=True)
class MapRow:
    spring: str
    loss_share: float
    column_m: float
    stroke_m: float
    trace: float  # of the one-period map; below -2 the cycle is unstable
    phase_per_pulse: float  # [rad] the one-period map's rotation
    decay_per_pulse: float  # its largest eigenvalue's modulus (> 1: unstable)
    phase_trim_10cm: float  # fraction of J, two-pulse deadbeat on a 10 cm phase error
    k_gain_x: float  # LQR: trim per m of position error
    k_gain_v: float  # per m/s


def map_rows() -> list[MapRow]:
    lo, hi = loss_shares()
    out = []
    for lab in THETAS:
        for d in (0.0, lo, hi):
            sp = spring(lab, d)
            ang, mod = rotation(sp)
            tr = float(np.trace(sp.mono))
            g = gains(sp, integrate=False)
            out.append(
                MapRow(
                    lab,
                    d,
                    sp.column,
                    sp.stroke,
                    tr,
                    ang,
                    mod,
                    two_pulse_phase_trim(sp),
                    float(g[0]),
                    float(g[1]),
                )
            )
    return out


# --- Output ---------------------------------------------------------------------------------

OUT_SLOPES = charge.OUT_DIR / "spring_trim_slopes.csv"
OUT_SCATTER = charge.OUT_DIR / "spring_trim_scatter.csv"
OUT_MAP = charge.OUT_DIR / "spring_trim_map.csv"
OUT_EDGE = charge.OUT_DIR / "spring_trim_edge.csv"
OUT_DRIFT = charge.OUT_DIR / "spring_trim_drift.csv"
OUT_TRIM = charge.OUT_DIR / "spring_trim_trim.csv"
OUT_LEVERS = charge.OUT_DIR / "spring_trim_levers.csv"
OUT_ROCK = charge.OUT_DIR / "spring_trim_rocking.csv"


def _write(path: Path, rows: list[object]) -> None:
    dicts = [asdict(r) for r in rows]  # type: ignore[call-overload]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(dicts[0]))
        w.writeheader()
        w.writerows(dicts)


def main() -> None:
    for sl in slope_rows():
        print(f"python: {sl.quantity:<44} {sl.value:9.4f}  (coarse {sl.coarse:.4f}) {sl.unit}")
    _write(OUT_SLOPES, list(slope_rows()))
    sc = scatter_rows()
    for sr in sc:
        print(
            f"python: scatter {sr.source:<22} heave {sr.heave:7.3%}  rock {sr.rock:7.3%}  {sr.note}"
        )
    _write(OUT_SCATTER, list(sc))
    mp = map_rows()
    for mr in mp:
        print(
            f"python: map {mr.spring:<5} d {mr.loss_share:.4f} | trace {mr.trace:.3f} phase"
            f" {mr.phase_per_pulse:.3f} rad"
            f" |lam| {mr.decay_per_pulse:.5f} | 10 cm phase trim {mr.phase_trim_10cm:.2%}"
        )
    _write(OUT_MAP, list(mp))
    edge = stability_edge()
    eo = gas_spring.column_for_theta(PLATE, edge, _gamma())
    print(
        f"python: the cycle goes unstable above theta = {edge / math.pi:.3f} pi:"
        f" column {eo.column:.2f} m, stroke {eo.stroke:.2f} m, peak/mean {eo.peak_over_mean:.3f}"
    )
    with OUT_EDGE.open("w") as fh:
        fh.write("theta_edge,theta_edge_over_pi,column_m,stroke_m,peak_over_mean\n")
        fh.write(f"{edge},{edge / math.pi},{eo.column},{eo.stroke},{eo.peak_over_mean}\n")
    dr = drift_rows()
    for dr_ in dr:
        print(
            f"python: drift {dr_.spring:<5} d {dr_.loss_share:.4f} s {dr_.scatter:4.0%} | exc p50"
            f" {dr_.excursion_p50_m:6.2f} p95 {dr_.excursion_p95_m:6.2f} m | 1 m after"
            f" {dr_.pulses_to_1m_p50:6.0f}, {LOST:.0f} m after {dr_.pulses_to_lost_p50:6.0f}"
            f" (lost {dr_.lost_share:.0%}) | dv {dr_.swing_error_rms:6.1f}"
            f" (est {dr_.swing_error_estimate:6.1f})"
        )
    _write(OUT_DRIFT, list(dr))
    tr = trim_rows()
    for tr_ in tr:
        print(
            f"python: trim {tr_.spring:<5} s {tr_.scatter:4.0%} t {tr_.timing_ms:.1f} ms int"
            f" {tr_.integrate!s:<5} lag {tr_.latency} | exc p50 {tr_.excursion_p50_m:5.2f} p95"
            f" {tr_.excursion_p95_m:5.2f} m | trim {tr_.trim_rms:6.2%} k rms {tr_.k_rms:.2f}"
            f" [{tr_.k_p01:.2f}, {tr_.k_p99:.2f}] sat {tr_.saturated:.1%} | craft dv"
            f" {tr_.craft_dv_p95:.3f} m/s {tr_.thruster_kg:.0f} kg | lost {tr_.lost_share:.0%}"
        )
    _write(OUT_TRIM, list(tr))
    lv = lever_rows()
    for lr in lv:
        print(
            f"python: lever {lr.lever:<52} s {lr.scatter:.3%} {lr.impulse_mns:8.2f} MN s"
            f" {lr.propellant_kg:8.0f} kg"
        )
    _write(OUT_LEVERS, list(lv))
    rk = rock_rows()
    for rr in rk:
        print(
            f"python: rock {rr.spring:<5} aim {rr.aim_m:.2f} m share {rr.rock_share:6.2%} trim"
            f" {rr.trimmed!s:<5} lag {rr.latency} | edge p50 {rr.edge_p50_m:6.2f} p95"
            f" {rr.edge_p95_m:6.2f} m | wedge rms {rr.wedge_rms:.1%} p99 {rr.wedge_p99:.1%}"
            f" | lost {rr.lost_share:.0%}"
        )
    _write(OUT_ROCK, list(rk))
    print(f"python: wrote {OUT_SLOPES.parent}/spring_trim_*.csv")


if __name__ == "__main__":
    main()
