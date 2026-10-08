"""Multigroup escape for the one-zone fireball: is the light of a dilute argon collision trapped?

The gray fireball (`fireball.collide`) uses one Planck and one Rosseland mean. That cannot say
whether the light sits in a band the gas absorbs. Cold argon is opaque above its 15.76 eV
ionization edge (3p photoionization, ~6e4 m^2/kg) and transparent below its 11.6 eV first
resonance line, so a gray mean hides exactly the question asked of a ~4 eV argon fireball.

Here the same homologous sphere radiates group by group from a TOPS (LANL OPLIB) multigroup pull,
`data/tables/tops/tops_argon_groups.html` (0.5-10 eV, 1e-6 to 0.1 kg/m^3, 120 log groups over
1-300 eV). In group g the luminosity uses the gray model's thin/thick bridge:

    L_g = 4 pi R^2 sigma T^4 b_g / (1 + 3/4 kR_g rho 2R + 1/(4/3 k_g rho R))

with b_g the Planck fraction in the group. Light below 1 eV is lumped into the first group and
light above 300 eV into the last. `k_g` is the group's Planck mean (`emission="planck"`) or its
Rosseland mean (`"rosseland"`). The Planck mean overstates emission where thick lines fill the
group, and the Rosseland mean understates it, so the pair brackets line saturation.

The light must also cross argon that has not been shocked yet. `f_through_cold_target`
attenuates the emitted spectrum by the whole cold target column `rho0 R_target` at the TOPS floor
(0.5 eV). That is a lower bound: it counts absorbed light as lost, when in fact it heats and
ionizes gas that is about to be shocked anyway.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Literal

import numpy as np

from puffsat.eos_water import Vec
from puffsat.tops import CM2G_TO_SI, GCC_TO_SI
from puffsat.xray_collision.eos import EV, K_B, build_table
from puffsat.xray_collision.fireball import (
    SIGMA_SB,
    STEP_FRACTION,
    T_STOP_EV,
    Kinematics,
    planck_above,
    shocked_state,
)
from puffsat.xray_collision.materials import ARGON

ARGON_GROUPS_PULL = Path("data/tables/tops/tops_argon_groups.html")
V_REL_ARGON_BAGS = 66.0e3  # [m/s] the rocket case as asked for the argon bags
AR_RESONANCE_EV = 11.6  # Ar I 3p^5 4s resonance lines (11.62, 11.83 eV): cold argon clear below
AR_EDGE_EV = 15.76  # Ar I ionization: cold argon opaque above
EUV_EV = 30.0
COLD_KT_EV = 0.5  # the TOPS floor, standing in for unshocked argon

Emission = Literal["planck", "rosseland"]

_HEADER = re.compile(r"Energy\s+Ross mg\s+Planck mg\s+for T, density =\s+(\S+)\s+(\S+)")
_NUMBER = re.compile(r"^\d\.\d+E[-+]\d+$")


@dataclass(frozen=True)
class TopsGroups:
    """TOPS multigroup means on a (kT, rho) grid, SI, logs held for interpolation."""

    edges_ev: Vec  # n_g + 1 group boundaries [eV]
    kt_ev: Vec  # [eV]
    rho: Vec  # [kg/m^3]
    ln_kr: Vec  # (n_T, n_rho, n_g), ln of [m^2/kg]
    ln_kp: Vec

    def kappa(self, kt_ev: float, rho: float) -> tuple[Vec, Vec]:
        """Group Rosseland and Planck means [m^2/kg], log-log bilinear, held at the grid edges."""
        x, y = np.log(self.kt_ev), np.log(self.rho)
        xq = min(max(math.log(kt_ev), float(x[0])), float(x[-1]))
        yq = min(max(math.log(rho), float(y[0])), float(y[-1]))
        i = int(np.clip(np.searchsorted(x, xq) - 1, 0, len(x) - 2))
        j = int(np.clip(np.searchsorted(y, yq) - 1, 0, len(y) - 2))
        a = (xq - x[i]) / (x[i + 1] - x[i])
        b = (yq - y[j]) / (y[j + 1] - y[j])

        def bil(f: Vec) -> Vec:
            return np.asarray(
                np.exp(
                    (1 - a) * (1 - b) * f[i, j]
                    + a * (1 - b) * f[i + 1, j]
                    + (1 - a) * b * f[i, j + 1]
                    + a * b * f[i + 1, j + 1]
                ),
                dtype=np.float64,
            )

        return bil(self.ln_kr), bil(self.ln_kp)


def parse_tops_groups(html: str, eg_high_kev: float = 0.3) -> TopsGroups:
    """Parse a TOPS multigroup results page. Each (T, rho) block lists every group's lower edge
    with its Rosseland and Planck means. The top edge is the request's `eg_high_kev`."""
    text = html.replace("<br/>", "\n").replace("&nbsp;", " ")
    text = re.sub(r"<[^>]+>", "", text)
    heads = list(_HEADER.finditer(text))
    blocks: dict[tuple[float, float], Vec] = {}
    for h, nxt in zip(heads, [*heads[1:], None], strict=True):
        body = text[h.end() : nxt.start() if nxt else len(text)]
        rows = []
        for line in body.splitlines():
            cols = line.split()
            if not cols:
                continue
            if len(cols) == 3 and all(_NUMBER.match(c) for c in cols):
                rows.append([float(c) for c in cols])
            elif rows:
                break  # the group list ends at the first non-row
        blocks[(float(h.group(1)), float(h.group(2)))] = np.array(rows)
    ts = sorted({k[0] for k in blocks})
    rhos = sorted({k[1] for k in blocks})
    lower = blocks[(ts[0], rhos[0])][:, 0]
    kr = np.empty((len(ts), len(rhos), len(lower)))
    kp = np.empty_like(kr)
    for i, t in enumerate(ts):
        for j, r in enumerate(rhos):
            kr[i, j] = blocks[(t, r)][:, 1]
            kp[i, j] = blocks[(t, r)][:, 2]
    floor = 1e-30
    return TopsGroups(
        edges_ev=np.append(lower, eg_high_kev) * 1e3,
        kt_ev=np.array(ts) * 1e3,
        rho=np.array(rhos) * GCC_TO_SI,
        ln_kr=np.log(np.maximum(kr * CM2G_TO_SI, floor)),
        ln_kp=np.log(np.maximum(kp * CM2G_TO_SI, floor)),
    )


@lru_cache(maxsize=2)
def load_tops_groups(path: Path = ARGON_GROUPS_PULL) -> TopsGroups:
    return parse_tops_groups(path.read_text())


def group_fractions(edges_ev: Vec, kt_ev: float) -> Vec:
    """Share of sigma T^4 in each group; the tails below and above the grid join the end groups."""
    above = np.array([planck_above(e / kt_ev) for e in edges_ev])
    frac = above[:-1] - above[1:]
    frac[0] += 1.0 - above[0]
    frac[-1] += above[-1]
    return np.asarray(frac, dtype=np.float64)


@dataclass(frozen=True)
class MultigroupResult:
    rho0: float  # both bags' density at the moment of impact [kg/m^3]
    m_impactor: float  # [kg]
    mass_ratio: float  # k: target mass at rest per kg of impactor
    v_rel: float  # [m/s]
    emission: str
    lowering: bool
    r_impactor: float  # [m] radius of the impactor's cloud at rho0
    r_target: float  # [m]
    heat: float  # Q [J]
    heat_fraction: float  # Q / impactor KE = k/(1+k)
    rho_shock: float
    kt_start_ev: float
    f_rad: float  # radiated share of Q
    f_below_resonance: float  # photons under 11.6 eV: cold argon is clear
    f_resonance_to_edge: float  # 11.6-15.76 eV: clear except in the lines
    f_edge_to_euv: float  # 15.76-30 eV: photoionizes cold argon
    f_above_euv: float  # above 30 eV
    f_through_cold_target: float  # f_rad after the whole cold target column (lower bound)
    t_half: float  # [s] time by which half the light is out
    t_end: float  # [s] time to cool to T_STOP_EV
    closure: float  # (E_rad + E_kin + E_int) / Q

    @property
    def light_of_ship_ke(self) -> float:
        return self.f_rad * self.heat_fraction


def collide_groups(
    rho0: float,
    m_impactor: float = 10.0,
    mass_ratio: float = 2.0,
    v_rel: float = V_REL_ARGON_BAGS,
    *,
    emission: Emission = "planck",
    lowering: bool = False,
    groups: TopsGroups | None = None,
) -> MultigroupResult:
    """Argon cloud of `m_impactor` at `rho0` into `mass_ratio` times its mass at rest, at the
    same density. One-zone fireball as in `fireball.collide`, with per-group escape (RK4)."""
    grp = load_tops_groups() if groups is None else groups
    table = build_table(ARGON, lowering, False)
    kin = Kinematics(m_impactor, mass_ratio * m_impactor, v_rel)
    mass, q = kin.m_total, kin.heat
    shock = shocked_state(table, 0.5 * v_rel, rho0)
    r0 = (3.0 * mass / (4.0 * math.pi * shock.rho)) ** (1.0 / 3.0)

    def derivs(y: Vec) -> tuple[Vec, float, float, float]:
        r, rdot, e_int = float(y[0]), float(y[1]), float(y[2])
        rho = mass / (4.0 / 3.0 * math.pi * r**3)
        temp = table.temp(rho, max(e_int, 0.0) / mass)
        p = table.state(rho, temp)[0]
        kt = K_B * temp / EV
        kr, kp = grp.kappa(kt, rho)
        k_emit = kp if emission == "planck" else kr
        g = 1.0 / (1.0 + 0.75 * kr * rho * 2.0 * r + 1.0 / (4.0 / 3.0 * k_emit * rho * r))
        lum_g = 4.0 * math.pi * r**2 * SIGMA_SB * temp**4 * group_fractions(grp.edges_ev, kt) * g
        lum = float(lum_g.sum())
        head = [
            rdot,
            20.0 * math.pi / 3.0 * r**2 * p / mass,
            -p * 4.0 * math.pi * r**2 * rdot - lum,
        ]
        t_hydro = r / max(abs(rdot), math.sqrt(p / rho))
        return np.concatenate((head, lum_g)), temp, lum, t_hydro

    n_g = len(grp.edges_ev) - 1
    y = np.concatenate(([r0, 0.0, q], np.zeros(n_g)))
    t = 0.0
    history: list[tuple[float, float]] = [(0.0, 0.0)]
    while True:
        k1, temp, lum, t_hydro = derivs(y)
        if temp < T_STOP_EV * EV / K_B or y[2] < 1e-6 * q:
            break
        dt = STEP_FRACTION * min(t_hydro, float(y[2]) / max(lum, 1e-300))
        k2 = derivs(y + 0.5 * dt * k1)[0]
        k3 = derivs(y + 0.5 * dt * k2)[0]
        k4 = derivs(y + dt * k3)[0]
        y = y + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        t += dt
        history.append((t, float(y[3:].sum())))

    spec = y[3:] / q
    lower = grp.edges_ev[:-1]
    e_rad = float(y[3:].sum())
    times, cum = zip(*history, strict=True)
    t_half = float(np.interp(0.5 * e_rad, cum, times))

    r_imp = (3.0 * m_impactor / (4.0 * math.pi * rho0)) ** (1.0 / 3.0)
    r_tgt = r_imp * mass_ratio ** (1.0 / 3.0)
    kr_cold, kp_cold = grp.kappa(COLD_KT_EV, rho0)
    k_cold = kp_cold if emission == "planck" else kr_cold

    def band(lo: float, hi: float) -> float:
        return float(spec[(lower >= lo) & (lower < hi)].sum())

    return MultigroupResult(
        rho0=rho0,
        m_impactor=m_impactor,
        mass_ratio=mass_ratio,
        v_rel=v_rel,
        emission=emission,
        lowering=lowering,
        r_impactor=r_imp,
        r_target=r_tgt,
        heat=q,
        heat_fraction=kin.heat_fraction,
        rho_shock=shock.rho,
        kt_start_ev=K_B * table.temp(shock.rho, q / mass) / EV,
        f_rad=e_rad / q,
        f_below_resonance=band(0.0, AR_RESONANCE_EV),
        f_resonance_to_edge=band(AR_RESONANCE_EV, AR_EDGE_EV),
        f_edge_to_euv=band(AR_EDGE_EV, EUV_EV),
        f_above_euv=band(EUV_EV, math.inf),
        f_through_cold_target=float((spec * np.exp(-k_cold * rho0 * r_tgt)).sum()),
        t_half=t_half,
        t_end=t,
        closure=(e_rad + 0.3 * mass * float(y[1]) ** 2 + float(y[2])) / q,
    )
