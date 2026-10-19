"""Result 13b: can a coated, thicker Kapton sheet fly through the old plume at 2 Hz and 1 g?

    PYTHONPATH=python uv run python -m puffsat.xray_collision.coated_sheet

Seth asked whether a JWST-style sheet -- an emissive doped-silicon face over a reflecting
aluminium (or silver) film, on ~100 um of Kapton -- survives the push at 2 Hz and 1 g mean.

**The plume** is the parent's expanding-cloud model (`templateArxiv.tex`,
`eq:jovian_plume_density`, `sec:charge_pusher_plume`), re-implemented here and checked against the
parent's printed doses. The k-th earlier pulse's cloud, of mass `m_g`, was launched `k tau` before
the sheet's impact. At cloud age `t` the sheet is `V (k tau - t)` from the plate, the cloud's
centre is at `U t`, and the gas there moves ballistically from the plate:

    rho = m_g / ((2 pi)^1.5 (u t)^3) exp(-[V k tau - (V + U) t]^2 / (2 u^2 t^2)),  w = V k tau / t,

and the sheet's face takes the incident kinetic flux `q = rho w^3 / 2` (full absorption, the
parent's bound). `U = (I - m_s V) / m_g` is what the pulse's impulse requires. The clouds add
independently.

**The sheet** is 1-D conduction through the Kapton (backward Euler), with the coating films as
lumped heat capacity at the faces. Each face radiates `eps sigma T^4`. The coatings change only
the faces' emissivity: a 50-150 nm film conducts across itself in nanoseconds, so a metal layer
under the emissive one is at the face's temperature and cannot reflect deposited heat. It reflects
*photons*; the plume's heat arrives as ~50 eV/amu particles that stop in the top ~0.1 um.

**The JWST lesson is the gap**, so a two-film sheet is also run: a thin sacrificial shield film a
short gap ahead of the main sheet. Until the shield reaches `T_FAIL` it stops all the plume gas,
and the main sheet sees only its infrared through an aluminised face. After that the main sheet
takes the plume directly (the shield's vapour is assumed to stop nothing).
"""

from __future__ import annotations

import csv
import math
from collections.abc import Callable
from dataclasses import asdict, dataclass

import numpy as np
from numpy.typing import NDArray

from puffsat.xray_collision import charge

SIGMA = 5.670374419e-8
G0 = 9.80665
V_CLOSE = 100.0e3  # [m/s]
LAUNCH_UNIT = 1500.0e3  # [kg]
K_PUSH = 6.74
VE_PUSH = 38.4e3  # [m/s] per kg of charge (calcs law.push_ve(6.74, 100))
SPREADS = (5.0e3, 8.0e3, 12.0e3)  # [m/s], the parent's one-axis u
T_START = 300.0  # [K], the parent's start
#: (label, cadence [Hz], mean thrust [g on LAUNCH_UNIT]).
PUSHES = (
    ("4 Hz, 12 MN s (priced)", 4.0, 12.0e6 * 4.0 / (LAUNCH_UNIT * G0)),
    ("1 Hz, 1 g", 1.0, 1.0),
    ("2 Hz, 1 g", 2.0, 1.0),
    ("2 Hz, 2 g", 2.0, 2.0),
)
#: The parent's printed doses [J/m^2] (u = 5 / 8 / 12) that the model must reproduce.
PARENT_DOSES = {
    "4 Hz, 12 MN s (priced)": (635.0e3, 251.0e3, 114.0e3),
    "1 Hz, 1 g": (47.0e3, 19.0e3, 8.0e3),
}

# --- Kapton (DuPont HN datasheet; c_p held at its room value, which overstates the rise).
KAPTON_RHO = 1420.0
KAPTON_CP = 1090.0
KAPTON_K = 0.12  # [W/m/K]
T_RATING = 673.15  # [K], the 400 C continuous rating
T_FAIL = 800.0  # [K], taken as rapid decomposition (TGA onset ~500-550 C; assumed, not verified)
# --- Coating films, lumped [J/m^2/K]: 50 nm Si, 100 nm Al (JWST layer 1, per NASA / Wikipedia).
SI_FILM = 2330.0 * 700.0 * 50e-9
AL_FILM = 2700.0 * 900.0 * 100e-9
#: Faces' infrared emissivity: (front toward the plume, back). Bracketing values, not measured
#: data: bare aluminium ~0.04, Kapton and doped silicon ~0.8, and the ideal 1.0 as the bound on
#: what any coating can do.
STACKS = {
    "bare Kapton": (0.8, 0.8, 0.0),
    "Si / Al / Kapton / Al (JWST layer 1)": (0.8, 0.04, SI_FILM + 2 * AL_FILM),
    "Si / Al / Kapton (Seth's)": (0.8, 0.8, SI_FILM + AL_FILM),
    "ideal emitter both faces": (1.0, 1.0, 0.0),
    "no radiation (parent's bound)": (0.0, 0.0, 0.0),
}
THICKNESSES = (100e-6, 150e-6, 200e-6, 250e-6)
SHIELD = 25e-6  # [m], the sacrificial front film of the two-film sheet
WINDOW = 3.0  # [s] of approach integrated (older clouds add < 1% of the dose)
DT = 2.0e-4  # [s]
CELLS = 40

OUT_DOSE = charge.OUT_DIR / "coated_sheet_dose.csv"
OUT_SHEET = charge.OUT_DIR / "coated_sheet.csv"


@dataclass(frozen=True)
class Push:
    label: str
    cadence: float  # [Hz]
    pulse: float  # [N s]
    sheet: float  # [kg] of impactor per pulse
    cloud: float  # [kg], m_g: the sheet plus its charge
    recession: float  # [m/s], U
    v_close: float = V_CLOSE  # [m/s], the sheets' closing speed


def push(
    label: str,
    cadence: float,
    g_mean: float,
    v_close: float = V_CLOSE,
    k: float = K_PUSH,
    v_e: float = VE_PUSH,
) -> Push:
    pulse = LAUNCH_UNIT * g_mean * G0 / cadence
    sheet = pulse / (k * v_e)
    cloud = sheet * (1.0 + k)
    return Push(label, cadence, pulse, sheet, cloud, (pulse - sheet * v_close) / cloud, v_close)


def flux(p: Push, u: float, before: NDArray[np.float64]) -> NDArray[np.float64]:
    """Incident kinetic flux [W/m^2] on the sheet `before` seconds ahead of its impact."""
    tau = 1.0 / p.cadence
    q = np.zeros_like(before)
    for k in range(1, int(WINDOW * p.cadence) + 40):
        t = k * tau - before  # the k-th cloud's age
        ok = t > 0.0
        tt = t[ok]
        rho = (
            p.cloud
            / ((2.0 * math.pi) ** 1.5 * (u * tt) ** 3)
            * np.exp(
                -((p.v_close * k * tau - (p.v_close + p.recession) * tt) ** 2) / (2 * (u * tt) ** 2)
            )
        )
        w = p.v_close * k * tau / tt
        q[ok] += 0.5 * rho * w**3
    return q


def _times() -> NDArray[np.float64]:
    """`before` from WINDOW down to 0, at the solver's step."""
    return np.linspace(WINDOW, 0.0, round(WINDOW / DT) + 1)


def dose(p: Push, u: float) -> float:
    s = _times()
    q = flux(p, u, s)
    return float(np.sum(0.5 * (q[1:] + q[:-1]) * DT))


@dataclass(frozen=True)
class Heating:
    face_max: float  # [K]
    mean_max: float  # [K]
    rating_before: float  # [s] before impact the face passes T_RATING (0 if never)
    fail_before: float  # [s] before impact it passes T_FAIL (0 if never)
    fail_depth: float  # [m], the deepest the sheet ever passes T_FAIL: the decomposed skin


@dataclass(frozen=True)
class Solid:
    """A sheet material: density, through-thickness conductivity, specific heat (a function of
    temperature) and the two face temperatures `heat` reports crossings of."""

    name: str
    rho: float  # [kg/m^3]
    k: float  # [W/m/K], through the thickness
    cp: Callable[[NDArray[np.float64]], NDArray[np.float64]]  # [J/kg/K] at T [K]
    t_rating: float  # [K]
    t_fail: float  # [K]


KAPTON = Solid(
    "Kapton",
    KAPTON_RHO,
    KAPTON_K,
    lambda t: np.full_like(t, KAPTON_CP),
    T_RATING,
    T_FAIL,
)


def heat(
    q: NDArray[np.float64],
    thickness: float,
    eps_front: float,
    eps_back: float,
    films: float,
    t0: float = T_START,
    solid: Solid = KAPTON,
) -> tuple[Heating, NDArray[np.float64]]:
    """Backward-Euler conduction through the sheet under face flux `q` (on `_times()`), with
    radiation linearised about the last step and the heat capacity taken at the last step's
    temperature. Returns the summary and the face history."""
    n = CELLS
    dx = thickness / n
    films_at = np.zeros(n)
    films_at[0] += 0.5 * films if eps_back > 0 else films
    films_at[-1] += 0.5 * films if eps_back > 0 else 0.0
    g = solid.k / dx
    temp = np.full(n, t0)
    face = np.empty(len(q))
    face[0] = t0
    mean_max = t0
    s = _times()
    cross = {solid.t_rating: 0.0, solid.t_fail: 0.0}
    deepest = 0
    conduct = np.zeros((n, n))
    for j in range(n):
        if j > 0:
            conduct[j, j] += g
            conduct[j, j - 1] -= g
        if j < n - 1:
            conduct[j, j] += g
            conduct[j, j + 1] -= g
    for i in range(1, len(q)):
        cap = solid.rho * solid.cp(temp) * dx + films_at  # [J/m^2/K] per cell
        a = conduct + np.diag(cap / DT)
        b = cap / DT * temp
        # Radiation, linearised: eps sigma T^4 ~ eps sigma (4 T0^3 T - 3 T0^4).
        for j, eps in ((0, eps_front), (n - 1, eps_back)):
            t_old = temp[j]
            a[j, j] += 4.0 * eps * SIGMA * t_old**3
            b[j] += 3.0 * eps * SIGMA * t_old**4
        b[0] += q[i]
        temp = np.linalg.solve(a, b)
        face[i] = temp[0]
        mean_max = max(mean_max, float(np.mean(temp)))
        deepest = max(deepest, int(np.sum(temp >= solid.t_fail)))
        for lim in cross:
            if cross[lim] == 0.0 and temp[0] >= lim:
                cross[lim] = float(s[i])
    return (
        Heating(
            float(face.max()),
            mean_max,
            cross[solid.t_rating],
            cross[solid.t_fail],
            deepest * dx,
        ),
        face,
    )


@dataclass(frozen=True)
class SheetRow:
    push: str
    u: float  # [km/s]
    dose: float  # [kJ/m^2]
    stack: str
    kapton_um: float
    face_max: float
    mean_max: float
    rating_before_ms: float
    fail_before_ms: float
    skin_um: float  # the decomposed skin (past T_FAIL)
    holds: bool  # face stays under T_RATING until the last 10 ms


LAST = 0.010  # [s]: what happens inside the last 10 ms is taken as harmless (S37.2 is open)


def two_film(
    q: NDArray[np.float64], main: float, eps_shield: float = 0.8, eps_main: float = 0.04
) -> tuple[Heating, float]:
    """A SHIELD-thick film ahead of the main sheet. Until the shield passes T_FAIL it takes all
    the flux (one lump: it is thin) and the main sheet sees its infrared through an aluminised
    face of `eps_main` (parallel plates). After that the main sheet takes the plume, conducted
    through its thickness by :func:`heat`. Returns the main sheet's heating and how long before
    impact the shield failed (0 if it never did)."""
    s = _times()
    c_shield = KAPTON_RHO * KAPTON_CP * SHIELD + SI_FILM + AL_FILM
    eff = 1.0 / (1.0 / eps_shield + 1.0 / eps_main - 1.0)
    ts = T_START
    q_main = np.zeros_like(q)
    failed_at = 0.0
    for i in range(1, len(q)):
        if failed_at == 0.0:
            ir = eff * SIGMA * (ts**4 - T_START**4)
            ts += DT * (q[i] - eps_shield * SIGMA * ts**4 - ir) / c_shield
            q_main[i] = ir
            if ts >= T_FAIL:
                failed_at = float(s[i])
        else:
            q_main[i] = q[i]
    h, _ = heat(q_main, main, eps_main, eps_main, 2 * AL_FILM)
    return h, failed_at


def main() -> None:
    s = _times()
    pushes = [push(label, f, gm) for label, f, gm in PUSHES]
    dose_rows = []
    for p in pushes:
        doses = [dose(p, u) for u in SPREADS]
        parent = PARENT_DOSES.get(p.label)
        q8 = flux(p, 8.0e3, s)
        peak_at = float(s[int(np.argmax(q8))])
        print(
            f"  {p.label:<22} {p.pulse / 1e6:5.2f} MN s, sheet {p.sheet:5.1f} kg, cloud"
            f" {p.cloud:5.0f} kg, U {p.recession / 1e3:4.1f} km/s | dose u = 5/8/12:"
            f" {'/'.join(f'{d / 1e3:.0f}' for d in doses)} kJ/m^2"
            + (f" (parent {'/'.join(f'{d / 1e3:.0f}' for d in parent)})" if parent else "")
            + f" | peak flux {q8.max() / 1e3:.0f} kW/m^2 at {1e3 * peak_at:.0f} ms out (u = 8)"
        )
        dose_rows.append(
            {
                "push": p.label,
                "pulse_MNs": p.pulse / 1e6,
                "sheet_kg": p.sheet,
                "cloud_kg": p.cloud,
                "U_kms": p.recession / 1e3,
                **{f"dose_u{u / 1e3:.0f}_kJ": d / 1e3 for u, d in zip(SPREADS, doses, strict=True)},
                "peak_flux_u8_kW": q8.max() / 1e3,
                "peak_before_u8_ms": 1e3 * peak_at,
            }
        )
    rows: list[SheetRow] = []
    for p in pushes[1:]:
        for u in SPREADS:
            q = flux(p, u, s)
            d = float(np.sum(0.5 * (q[1:] + q[:-1]) * DT))
            for name, (ef, eb, films) in STACKS.items():
                for th in THICKNESSES:
                    h, _ = heat(q, th, ef, eb, films)
                    rows.append(
                        SheetRow(
                            p.label,
                            u / 1e3,
                            d / 1e3,
                            name,
                            1e6 * th,
                            h.face_max,
                            h.mean_max,
                            1e3 * h.rating_before,
                            1e3 * h.fail_before,
                            1e6 * h.fail_depth,
                            h.rating_before <= LAST,
                        )
                    )
            for th in THICKNESSES:
                h, gone = two_film(q, th)
                rows.append(
                    SheetRow(
                        p.label,
                        u / 1e3,
                        d / 1e3,
                        f"two-film: {SHIELD * 1e6:.0f} um shield + gap (shield gone"
                        f" {1e3 * gone:.0f} ms out)",
                        1e6 * th,
                        h.face_max,
                        h.mean_max,
                        1e3 * h.rating_before,
                        1e3 * h.fail_before,
                        1e6 * h.fail_depth,
                        h.rating_before <= LAST,
                    )
                )
    for r in rows:
        if r.push == "2 Hz, 1 g" or r.kapton_um == 100.0:
            print(
                f"  {r.push:<10} u {r.u:2.0f} {r.dose:4.0f} kJ/m^2 | {r.stack:<38}"
                f" {r.kapton_um:3.0f} um | face {r.face_max:5.0f} K, mean {r.mean_max:5.0f} K"
                f" | past 400 C {r.rating_before_ms:4.0f} ms out, past {T_FAIL:.0f} K"
                f" {r.fail_before_ms:4.0f} ms, skin {r.skin_um:3.0f} um"
                f" | {'holds' if r.holds else 'fails'}"
            )
    OUT_DOSE.parent.mkdir(parents=True, exist_ok=True)
    with OUT_DOSE.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(dose_rows[0]))
        w.writeheader()
        w.writerows(dose_rows)
    with OUT_SHEET.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(asdict(rows[0])))
        w.writeheader()
        w.writerows(asdict(r) for r in rows)
    print(f"python: wrote {OUT_DOSE}, {OUT_SHEET}")


if __name__ == "__main__":
    main()
