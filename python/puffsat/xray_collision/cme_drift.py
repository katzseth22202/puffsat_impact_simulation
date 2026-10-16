"""CME density at depth, non-radial CME flow, and the drift of an edge-on Sun-wave disk.

    PYTHONPATH=python uv run python -m puffsat.xray_collision.cme_drift

Answers item 4 (asks a and b) of the parent's `docs/solar_environment_answers_for_impact_sim.md`
(`Balloon-Pulse-Propulsion` 53ca8e8), and checks its item 5. Every input is a cited literature
value. The scaling to 4 and 9.86 R_sun, the drift and the sunlight integral are computed here.

**(a) Density at depth** (`cme_density.csv`). Three estimates, high to low:

- *Coronagraph, per event.* Temmer et al. 2021 (JGR 126, e2020JA028380, Table A1) give, for 29
  well-observed CMEs of 2008-2014, the deprojected mass at ~15 R_sun (STEREO COR2) over the GCS
  flux-rope volume at 21.5 R_sun. They assume self-similar expansion (volume ~ r^3) and a uniform
  plasma. Here that density is scaled inward as r^-3, times the share of the mass already lifted.
  Bein et al. 2013 (ApJ 768, 31) find the mass still grows by 2-6% of its final value per R_sun,
  so the mass at r is bracketed as `1 - f (15 - r)` of the 15 R_sun mass, with f = 0.02 and 0.06.
- *In situ, extrapolated.* Temmer & Bothmer 2022 (A&A 665, A70) fit Helios 0.3-1 au ICMEs:
  magnetic ejecta `7.1 R^-2.4` and sheath `22.3 R^-1.7` cm^-3, with R in au. The sheath fit meets
  the ambient wind at ~13 R_sun, which they read as where sheaths start to form. They note the
  coronagraph densities at 15-30 R_sun lie above the inward extrapolation, but within one standard
  deviation of it.
- *Parker, measured.* Romeo et al. 2023 (ApJ 954, 168) on the 2022-09-05 CME at 13.3-15 R_sun:
  electron density 1962 -> 3989 cm^-3 across the shock (their Table 2), ~1e4 cm^-3 in the
  closed-field interval, and below 1000 cm^-3 in most of the rest. Parker crossed the flank and
  the wake, not the core.

Peak to mean: the measured shock compression is 2.0. The sheath-to-ejecta ratio of Temmer et al.'s
29 events at 1 au has a median of 2.3, computed here from their table. A strong shock bounds the
sheath at 4 times the gas ahead of it.

**(b) Non-radial flow** (`cme_tilt.csv`).

- *Near the Sun:* Parker's shock downstream (Romeo et al. Table 2) is v_RTN = (455, -358, 251)
  km/s. The normal component is 0.40 of the speed, a 23 degree tilt out of the RT plane, and the
  shock normal is 23 degrees out of it too. This is the flank sheath, where plasma is pushed
  around the ejecta.
- *At 1 au:* Al-Haddad et al. 2022 (ApJ 927, 68; STEREO-A PLASTIC, 2007-2019). Inside
  ejecta, the normal component has a median of 7 km/s. That is ~0.9 degrees on the 1 au ejecta
  speed of Temmer et al.'s events (median computed here). The largest CME deflection in 13 years
  is 18.4 degrees, in a sheath; only ~2% of all points lie more than 5 degrees north or south.

**The drift** (`cme_drift.csv`). This takes the parent's disk (1 kg, 0.5 m, 2.8 mm graphite) and
its passage times (1 h at 4 R_sun, 2.5 h at 9.86 R_sun). It flies face normal out of the orbital
plane at periapsis speed, tangentially. The plasma flows at speed v, tilted phi out of the plane.
The relative flow tilts `asin(v sin phi / sqrt(v^2 + v_disk^2))` onto the face. The disk is
absorbing, so its force is `rho v_rel^2 (A_face sin + A_edge cos)` along the flow. The CME's own
speed carries most of the ram: at 1000 km/s it is ~91% of `v_rel^2` at 4 R_sun (disk 307 km/s)
and ~96% at 9.86 R_sun (`cme_share_of_ram`). The face term goes as `rho v_rel v sin(phi)`. The
disk's own motion dilutes the tilt but raises the ram, so it barely changes the face push.
The flow is held for the whole passage, so `dv = F t/m` and `drift = F t^2/(2m)`, as in the
parent's table. Cold gas is priced at the parent's 60 s.

**Item 5, edge-on sunlight** (`edge_sunlight`). The share of the face-on outward push that an
edge-on disk still feels, integrated over the finite solar disk. Uniform, it is `~4 alpha/(3 pi)`
for angular radius alpha. Limb darkening (linear, u = 0.6) concentrates the light toward the centre
and lowers it. The edge's own `D t` is added.
"""

from __future__ import annotations

import csv
import math
import statistics
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from puffsat.xray_collision import charge, mission

M_P = 1.67262192e-27  # [kg]; Temmer et al. convert mass density with the proton mass
DEPTHS_RSUN = (4.0, 9.86)
#: The parent's CME passage time at each depth (answers doc, E2).
PASSAGE = {4.0: 3600.0, 9.86: 2.5 * 3600.0}  # [s]
#: The parent's disk: 1 kg, 0.5 m graphite, 2.8 mm thick (answers doc, E2).
DISK_MASS = 1.0  # [kg]
DISK_D = 0.5  # [m]
DISK_T = 2.8e-3  # [m]
COLD_GAS_ISP = 60.0  # [s], the parent's cold gas
G0 = 9.80665
#: The parent's 2022-09-05 CME speed (Patel et al. 2023), the extreme row.
EXTREME_SPEED = 2500.0e3  # [m/s]
#: CME speeds swept at the coronagraph median density: the push goes as the speed squared.
SPEED_SWEEP = (500.0e3, 1000.0e3, 2500.0e3)  # [m/s]
#: The parent's assumed CME densities at 4 R_sun (answers doc, E2), [cm^-3] -> pressure check.
PARENT_ASSUMED = {"typical": (1.0e6, 1.0e3), "extreme": (1.0e7, 2.45e3)}  # (n, v km/s)
#: Bein et al. 2013: the CME mass still grows by 2-6% of its final mass per R_sun.
MASS_GROWTH = (0.02, 0.06)
MASS_REFERENCE_RSUN = 15.0  # Temmer et al.'s mass is measured at ~15 R_sun
DENSITY_REFERENCE_RSUN = 21.5  # ...and its density is quoted at 21.5 R_sun

#: Temmer et al. 2021, Table A1: (event, deprojected mass [g], CME speed at 21.5 R_sun [km/s],
#: mass density at 21.5 R_sun [kg/m^3], 1 au ME speed [km/s], 1 au sheath and ME density
#: [cm^-3]).
TEMMER_2021 = (
    ("2008/12/12", 4.50e15, 470.0, 1.36e-17, 340.0, 15.7, 3.9),
    ("2009/12/16", 2.20e15, 370.0, 4.81e-18, 380.0, 3.0, 3.7),
    ("2010/04/03", 6.04e15, 900.0, 1.42e-17, 650.0, 9.4, 4.2),
    ("2010/04/08", 7.52e15, 500.0, 3.08e-17, 410.0, 9.7, 10.1),
    ("2010/05/24", 3.20e15, 390.0, 4.36e-18, 360.0, 19.2, 7.2),
    ("2010/06/16", 2.35e15, 430.0, 7.97e-18, 360.0, 5.9, 5.9),
    ("2010/10/26", 7.14e15, 450.0, 1.86e-17, 340.0, 13.4, 7.8),
    ("2011/01/30", 5.60e15, 310.0, 1.84e-17, 410.0, 10.6, 18.0),
    ("2011/02/15", 6.83e15, 700.0, 8.95e-18, 470.0, 11.6, 1.4),
    ("2011/03/03", 3.13e15, 510.0, 5.83e-18, 440.0, 4.2, 5.1),
    ("2011/06/02", 3.80e15, 830.0, 5.78e-18, 510.0, 34.3, 14.9),
    ("2011/06/14", 1.03e16, 770.0, 2.14e-17, 490.0, 5.6, 7.0),
    ("2011/08/04", 6.84e15, 1300.0, 6.93e-18, 540.0, 7.5, 1.2),
    ("2011/09/13", 4.20e15, 600.0, 5.98e-18, 450.0, 13.3, 4.3),
    ("2011/10/22", 1.23e16, 610.0, 1.01e-17, 470.0, 25.6, 11.0),
    ("2011/10/27", 3.01e15, 700.0, 5.76e-18, 370.0, 7.5, 3.7),
    ("2011/11/26", 1.00e16, 1350.0, 7.67e-18, 450.0, 9.2, 13.4),
    ("2012/01/19", 9.17e15, 1000.0, 9.60e-18, 450.0, 26.6, 6.7),
    ("2012/03/13", 1.00e16, 850.0, 6.02e-18, 710.0, 9.1, 4.7),
    ("2012/06/14", 8.46e15, 1000.0, 8.39e-18, 450.0, 37.7, 22.6),
    ("2012/07/12", 1.84e16, 1150.0, 1.75e-17, 490.0, 15.8, 2.3),
    ("2012/09/28", 9.61e15, 1000.0, 1.34e-17, 370.0, 20.3, 10.8),
    ("2012/10/05", 6.40e15, 600.0, 9.05e-18, 400.0, 15.7, 5.2),
    ("2012/11/09", 5.19e15, 600.0, 5.77e-18, 380.0, 22.3, 6.5),
    ("2012/11/23", 3.45e15, 680.0, 4.06e-18, 450.0, 7.3, 2.6),
    ("2013/04/11", 1.51e16, 700.0, 2.18e-17, 410.0, 13.4, 4.1),
    ("2013/07/09", 2.93e15, 550.0, 4.57e-18, 410.0, 5.7, 3.7),
    ("2013/09/29", 1.37e16, 1000.0, 1.65e-17, 470.0, 12.2, 1.8),
    ("2014/08/15", 1.60e15, 450.0, 3.68e-18, 360.0, 18.2, 8.1),
)
#: Temmer & Bothmer 2022, Fig. 4: Helios fits n = c R[au]^p, [cm^-3].
HELIOS_FITS = {"magnetic ejecta": (7.1, -2.4), "sheath": (22.3, -1.7), "ambient wind": (7.0, -2.1)}
#: Romeo et al. 2023, Table 2: Parker's shock at 2022-09-05, upstream and downstream.
PSP_UPSTREAM_N = 1962.0  # [cm^-3]
PSP_DOWNSTREAM_N = 3989.0
PSP_DOWNSTREAM_V = (455.0, -358.0, 251.0)  # [km/s], RTN
PSP_CLOSED_FIELD_N = 1.0e4  # [cm^-3], the bidirectional-electron interval
PSP_RSUN = 13.3
#: Al-Haddad et al. 2022: the normal flow inside ejecta at 1 au (median), and the largest CME
#: deflection of 2007-2019 (a sheath).
ALHADDAD_EJECTA_VN = 7.0  # [km/s]
ALHADDAD_SHEATH_MAX_DEG = 18.4
STRONG_SHOCK_COMPRESSION = 4.0
LIMB_DARKENING_U = 0.6

OUT_DENSITY = charge.OUT_DIR / "cme_density.csv"
OUT_TILT = charge.OUT_DIR / "cme_tilt.csv"
OUT_DRIFT = charge.OUT_DIR / "cme_drift.csv"
OUT_SUNLIGHT = charge.OUT_DIR / "cme_edge_sunlight.csv"


def coronagraph_density(rho_21: float, r_sun: float, growth: float) -> float:
    """Number density [cm^-3] at r_sun from Temmer et al.'s 21.5 R_sun density: self-similar
    r^-3, times the share of the 15 R_sun mass already lifted, `1 - growth (15 - r)`."""
    lifted = 1.0 - growth * max(MASS_REFERENCE_RSUN - r_sun, 0.0)
    n_21 = rho_21 / M_P * 1e-6
    return n_21 * (DENSITY_REFERENCE_RSUN / r_sun) ** 3 * lifted


def helios_density(structure: str, r_sun: float) -> float:
    c, p = HELIOS_FITS[structure]
    return float(c * (r_sun * charge.R_SUN / charge.AU) ** p)


def peak_to_mean() -> dict[str, float]:
    ratios = [sh / me for *_, sh, me in TEMMER_2021]
    return {
        "psp_shock_compression": PSP_DOWNSTREAM_N / PSP_UPSTREAM_N,
        "sheath_over_ejecta_1au_median": statistics.median(ratios),
        "strong_shock_bound": STRONG_SHOCK_COMPRESSION,
    }


@dataclass(frozen=True)
class DensityRow:
    r_sun: float
    estimate: str
    n_median: float  # [cm^-3]
    n_high: float  # [cm^-3]: the densest of 29 events, or the same value for a fit
    speed_median: float  # [km/s]
    pressure_median: float  # [Pa], rho v^2
    pressure_high: float  # [Pa], densest event x shock compression at EXTREME_SPEED


def density_rows() -> list[DensityRow]:
    comp = peak_to_mean()["psp_shock_compression"]
    rows: list[DensityRow] = []
    rs = [*DEPTHS_RSUN, PSP_RSUN]
    for r in rs:
        for f in MASS_GROWTH:
            ns = [coronagraph_density(rho, r, f) for _, _, _, rho, *_ in TEMMER_2021]
            ps = [
                coronagraph_density(rho, r, f) * 1e6 * M_P * (v * 1e3) ** 2
                for _, _, v, rho, *_ in TEMMER_2021
            ]
            rows.append(
                DensityRow(
                    r_sun=r,
                    estimate=f"coronagraph, mass grows {f:.0%}/R_sun",
                    n_median=statistics.median(ns),
                    n_high=max(ns),
                    speed_median=statistics.median(v for _, _, v, *_ in TEMMER_2021),
                    pressure_median=statistics.median(ps),
                    pressure_high=max(ns) * comp * 1e6 * M_P * EXTREME_SPEED**2,
                )
            )
        for s in ("magnetic ejecta", "sheath", "ambient wind"):
            n = helios_density(s, r)
            v = statistics.median(v for _, _, v, *_ in TEMMER_2021)
            rows.append(
                DensityRow(
                    r_sun=r,
                    estimate=f"Helios fit, {s}",
                    n_median=n,
                    n_high=n,
                    speed_median=v,
                    pressure_median=n * 1e6 * M_P * (v * 1e3) ** 2,
                    pressure_high=n * comp * 1e6 * M_P * EXTREME_SPEED**2,
                )
            )
    return rows


@dataclass(frozen=True)
class TiltRow:
    source: str
    tilt_deg: float  # out of the orbital (RT) plane, Sun frame
    out_of_plane_share: float  # v_N / |v|


def tilt_rows() -> list[TiltRow]:
    vr, vt, vn = PSP_DOWNSTREAM_V
    speed = math.sqrt(vr * vr + vt * vt + vn * vn)
    v_me = statistics.median(vme for *_, vme, _, _ in TEMMER_2021)
    sheath = math.sin(math.radians(ALHADDAD_SHEATH_MAX_DEG))
    return [
        TiltRow("radial (the parent's item 1)", 0.0, 0.0),
        TiltRow(
            "1 au ejecta, median normal flow (Al-Haddad)",
            math.degrees(math.asin(ALHADDAD_EJECTA_VN / v_me)),
            ALHADDAD_EJECTA_VN / v_me,
        ),
        TiltRow("1 au sheath, largest of 2007-2019 (Al-Haddad)", ALHADDAD_SHEATH_MAX_DEG, sheath),
        TiltRow(
            "13-15 R_sun flank sheath, Parker (Romeo)",
            math.degrees(math.asin(vn / speed)),
            vn / speed,
        ),
    ]


@dataclass(frozen=True)
class Drift:
    r_sun: float
    case: str
    tilt_source: str
    plasma_tilt_deg: float
    n: float  # [cm^-3]
    speed: float  # [km/s]
    relative_tilt_deg: float
    disk_speed: float  # [km/s], tangential at periapsis
    cme_share_of_ram: float  # v^2 / (v^2 + v_disk^2): the CME's own flow's share of rho v_rel^2
    pressure: float  # [Pa]
    force: float  # [N]
    dv: float  # [m/s]
    drift: float  # [m]
    cold_gas: float  # [kg]


def disk_areas() -> tuple[float, float]:
    return math.pi * DISK_D**2 / 4.0, DISK_D * DISK_T


def drift(r_sun: float, n: float, speed: float, tilt_deg: float) -> tuple[float, ...]:
    """(relative tilt [deg], rho v_rel^2 [Pa], force [N], dv [m/s], drift [m]) for an edge-on
    disk at periapsis speed, held in plasma of density n [cm^-3], speed [m/s], tilted out of the
    orbital plane by tilt_deg, for the parent's passage time."""
    v_disk = mission.periapsis(r_sun).periapsis_speed
    phi = math.radians(tilt_deg)
    v_rel = math.hypot(speed, v_disk)
    rel = math.asin(speed * math.sin(phi) / v_rel)
    face, edge = disk_areas()
    pressure = n * 1e6 * M_P * v_rel**2
    force = pressure * (face * abs(math.sin(rel)) + edge * math.cos(rel))
    t = PASSAGE[r_sun]
    return (
        math.degrees(rel),
        pressure,
        force,
        force * t / DISK_MASS,
        0.5 * force * t * t / DISK_MASS,
    )


def drift_rows() -> list[Drift]:
    comp = peak_to_mean()["psp_shock_compression"]
    cases: list[tuple[float, str, float, float]] = []
    v_med = 1e3 * statistics.median(v for _, _, v, *_ in TEMMER_2021)
    for r in DEPTHS_RSUN:
        hi = max(coronagraph_density(rho, r, MASS_GROWTH[0]) for _, _, _, rho, *_ in TEMMER_2021)
        med = statistics.median(
            coronagraph_density(rho, r, MASS_GROWTH[0]) for _, _, _, rho, *_ in TEMMER_2021
        )
        cases += [
            (
                r,
                "in situ (Helios ejecta fit), median speed",
                helios_density("magnetic ejecta", r),
                v_med,
            ),
            (r, "coronagraph median, median speed", med, v_med),
            (r, "densest of 29 x shock 2.0, at 2500 km/s", hi * comp, EXTREME_SPEED),
        ]
        for v_sweep in SPEED_SWEEP:
            cases.append((r, f"coronagraph median, {v_sweep / 1e3:.0f} km/s", med, v_sweep))
    out: list[Drift] = []
    for r, case, n, v in cases:
        v_disk = mission.periapsis(r).periapsis_speed
        for tr in tilt_rows():
            rel, p, f, dv, x = drift(r, n, v, tr.tilt_deg)
            out.append(
                Drift(
                    r_sun=r,
                    case=case,
                    tilt_source=tr.source,
                    plasma_tilt_deg=tr.tilt_deg,
                    n=n,
                    speed=v / 1e3,
                    relative_tilt_deg=rel,
                    disk_speed=v_disk / 1e3,
                    cme_share_of_ram=v * v / (v * v + v_disk * v_disk),
                    pressure=p,
                    force=f,
                    dv=dv,
                    drift=x,
                    cold_gas=DISK_MASS * dv / (COLD_GAS_ISP * G0),
                )
            )
    return out


def edge_sunlight(r_sun: float, limb_u: float, n: int = 400) -> float:
    """Outward push on an edge-on disk (face normal out of the orbital plane) as a share of the
    face-on push, integrated over a solar disk of angular radius asin(1/r_sun) with linear limb
    darkening I = 1 - u (1 - mu). Absorbing faces; the edge's projected `D t` is included."""
    alpha = math.asin(1.0 / r_sun)
    rho = (np.arange(n) + 0.5) / n * alpha
    psi = (np.arange(n) + 0.5) / n * 2.0 * math.pi
    rr, pp = np.meshgrid(rho, psi, indexing="ij")
    mu = np.sqrt(np.clip(1.0 - (np.sin(rr) / math.sin(alpha)) ** 2, 0.0, 1.0))
    weight = (1.0 - limb_u * (1.0 - mu)) * np.sin(rr)  # intensity x dOmega / (drho dpsi)
    radial = np.cos(rr)
    normal = np.abs(np.sin(rr) * np.sin(pp))
    face_on = float(np.sum(weight * radial * radial))
    faces = float(np.sum(weight * normal * radial))
    face, edge = disk_areas()
    edge_share = float(np.sum(weight * radial * radial)) * edge / face
    return (faces + edge_share) / face_on


def sunlight_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for r in DEPTHS_RSUN:
        alpha = math.asin(1.0 / r)
        rows.append(
            {
                "r_sun": r,
                "sun_angular_radius_deg": math.degrees(alpha),
                "uniform_small_angle": 4.0 * alpha / (3.0 * math.pi),
                "uniform": edge_sunlight(r, 0.0),
                "limb_darkened": edge_sunlight(r, LIMB_DARKENING_U),
            }
        )
    return rows


def _write(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    dens = density_rows()
    for d in dens:
        print(
            f"  {d.r_sun:5.2f} R_sun {d.estimate:<38} n {d.n_median:8.2e} (high {d.n_high:8.2e})"
            f" cm^-3, rho v^2 {d.pressure_median:8.2e} Pa (high {d.pressure_high:8.2e})"
        )
    print("  peak to mean:", {k: round(v, 2) for k, v in peak_to_mean().items()})
    for n_par, v_par in PARENT_ASSUMED.values():
        print(
            f"  parent assumed {n_par:.0e} cm^-3 at {v_par:.0f} km/s:"
            f" {n_par * 1e6 * M_P * (v_par * 1e3) ** 2:.2e} Pa"
        )
    tilts = tilt_rows()
    for t in tilts:
        print(f"  tilt {t.tilt_deg:5.1f} deg ({t.out_of_plane_share:.2f}): {t.source}")
    drifts = drift_rows()
    for x in drifts:
        print(
            f"  {x.r_sun:5.2f} R_sun {x.case[:42]:<42} {x.plasma_tilt_deg:5.1f} deg:"
            f" dv {x.dv:8.3f} m/s, drift {x.drift / 1e3:8.3f} km, gas {1e3 * x.cold_gas:6.1f} g"
        )
    sun = sunlight_rows()
    for s in sun:
        print(
            f"  edge-on sunlight at {s['r_sun']} R_sun: uniform {float(str(s['uniform'])):.3f},"
            f" limb-darkened {float(str(s['limb_darkened'])):.3f} of face-on"
        )
    _write(OUT_DENSITY, [asdict(d) for d in dens])
    _write(OUT_TILT, [asdict(t) for t in tilts])
    _write(OUT_DRIFT, [asdict(x) for x in drifts])
    _write(OUT_SUNLIGHT, sun)
    print(f"python: wrote {OUT_DENSITY}, {OUT_TILT}, {OUT_DRIFT}, {OUT_SUNLIGHT}")


if __name__ == "__main__":
    main()
