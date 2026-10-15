"""Freeze-out in the charge pusher's fireball: do the energies locked in ions and bonds return?

    PYTHONPATH=python uv run python -m puffsat.xray_collision.freeze

Result 10's speeds assume every joule of heat ends as expansion. Two stores could keep some:

- **Ionization.** About 18% of the fireball's energy starts as ionization. It returns only if ions
  recapture electrons while the gas is still dense enough. The store is followed as a relaxation,
  `d e_ion / dt = -(e_ion - e_ion,eq) / tau_rec`, against Saha equilibrium on the expansion. The
  rate is three-body plus radiative recombination, using `puffsat.recombination`'s hydrogenic
  coefficients (Zel'dovich & Raizer; Osterbrock & Ferland). They are conservative for multiply
  charged carbon, whose capture is faster.
- **Bonds.** The CH2 table counts atoms only, so polyethylene's 84 MJ/kg of bonds is never paid
  back in it. Hydrogen's share re-forms as H2 if `H + H + H -> H2 + H` beats the expansion where
  H2 becomes favoured. That rate is Jacobs, Giedt & Cohen 1967 (`walled_nozzle.rates.H_ATOMIC_H`),
  reported with its margin in decades. Carbon's share has no rate to take
  (`walled_nozzle.rates.CARBON_CHANNEL_NOTE`: it is a nucleation problem), so it is bracketed
  between frozen and condensed.

The expansion is the one-zone homologous sphere of `fireball.py`, run on the CH2 table from its
1284 kg/m^3 ceiling until the gas reaches the plate. The race is a Damkohler number,
`Da = tau_exp / tau_rec` (`recombination.verdict`).
"""

from __future__ import annotations

import csv
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from puffsat import recombination as rc
from puffsat.walled_nozzle import rates
from puffsat.water_plate import spray_eos
from puffsat.xray_collision import charge
from puffsat.xray_collision.eos import EV, K_B, EosTable

CH2 = spray_eos.AtomMixture((("C", 1.0), ("H", 2.0)))
OUT = Path("data/results/xray_collision/charge_freeze.csv")
EV_J = 1.602176634e-19
AMU = 1.66053906660e-27
H2_BOND = 4.478 * EV_J  # D0 of H2
CH2_ATOMIZATION = (2.0 * 4.3 + 3.6) * EV_J  # = charge.PE_BONDS per CH2 unit
H2_WINDOW = (2500.0, 4500.0)  # [K] where H2 becomes favoured at these densities


def ion_energy(rho: float, temp: float) -> tuple[float, float, float]:
    """`(e_ion [J/kg], n_e, n_heavy)` in Saha equilibrium for CH2 at `(rho, temp)`."""
    s = spray_eos.state(CH2, rho, temp)
    e = 0.0
    for name, n in s.n_element.items():
        e += n * float(np.dot(s.fractions[name], spray_eos._cumulative_ip(name)))
    return e / rho, s.n_e, sum(s.n_element.values())


@dataclass(frozen=True)
class FreezeResult:
    k: float
    t_plate: float  # [s] when the fireball's edge reaches the plate's focal distance
    rho_plate: float
    kt_plate_ev: float
    ion_share_start: float  # of the heat
    ion_share_frozen: float  # still locked at the plate, of the heat
    da_ion_min: float  # smallest Damkohler number met on the way
    t_h2_window: float  # [K] where the H2 race is read (the window's middle, if reached)
    da_h2: float
    h2_margin_decades: float  # log10(Da), against the rate's stated uncertainty
    h2_extrapolated: bool
    bond_share: float  # bond toll as a share of the heat
    velocity_loss_ions: float  # 1 - sqrt(1 - frozen ion share)
    velocity_loss_bonds_carbon_frozen: float  # H2 re-formed, carbon left atomic
    velocity_loss_bonds_all_frozen: float  # 1 once the bonds outweigh the heat (no exhaust)
    velocity_loss_worst: float  # ions frozen as found, plus every bond


def freeze(
    k: float, table: EosTable, plate_distance: float = 4.17, v_rel: float = charge.V_REL
) -> FreezeResult:
    """Expand the k fireball (1 kg projectile closing at `v_rel`) and race both stores
    against it."""
    m = 1.0 + k
    heat = 0.5 * k / (1.0 + k) * v_rel**2  # J per kg of projectile
    e = heat / m
    rho = float(np.exp(table.ln_rho[-1]))
    r = (3.0 * m / (4.0 * math.pi * rho)) ** (1.0 / 3.0)
    rdot, t = 0.0, 0.0
    e_ion0 = ion_energy(rho, table.temp(rho, e))[0]
    x = e_ion0  # the ionization store actually held
    da_min = math.inf
    h2 = (math.nan, math.nan)
    t_mid = 0.5 * (H2_WINDOW[0] + H2_WINDOW[1])
    while True:
        rho = m / (4.0 / 3.0 * math.pi * r**3)
        temp = table.temp(rho, e)
        p = table.state(rho, temp)[0]
        x_eq, n_e, _ = ion_energy(rho, temp)
        tau_exp = r / (3.0 * rdot) if rdot > 0.0 else math.inf
        alpha = rc.three_body_recombination(temp, n_e) + rc.radiative_recombination(temp)
        tau_rec = 1.0 / (alpha * n_e) if n_e > 0.0 else math.inf
        if math.isfinite(tau_exp) and x > 1e-6 * e:
            da_min = min(da_min, tau_exp / tau_rec)
        if math.isnan(h2[0]) and temp < t_mid and math.isfinite(tau_exp):
            n_h = rho * CH2.per_kg["H"]
            tau_h2 = 1.0 / (rates.H_ATOMIC_H(temp) * n_h * n_h)
            h2 = (temp, tau_exp / tau_h2)
        edge = r * math.sqrt(5.0 / 3.0)  # a uniform sphere's rms radius is sqrt(3/5) of its edge
        if edge >= plate_distance or temp < H2_WINDOW[0]:
            break
        # One step: homologous sphere (fireball.py's Rddot), adiabatic energy, relaxing store.
        # The relaxation below is the exact exponential over the step, so a recombination time far
        # shorter than the step (solid density) just pins the store to equilibrium.
        dt = 0.02 * r / max(rdot, math.sqrt(p / rho))
        rddot = 20.0 * math.pi / 3.0 * r**2 * p / m
        dv = 4.0 * math.pi * r**2 * rdot * dt
        e -= p * dv / m
        x += -(x - x_eq) * (1.0 - math.exp(-dt / tau_rec)) if math.isfinite(tau_rec) else 0.0
        rdot += rddot * dt
        r += rdot * dt
        t += dt
    frozen = max(x - ion_energy(rho, temp)[0], 0.0) / (heat / m)
    bond = charge.PE_BONDS * k / (k + 1.0) / (heat / m)
    # Carbon's share is bracketed as wholly frozen (atoms) or wholly condensed (returned).
    carbon_part = (CH2_ATOMIZATION - H2_BOND) / CH2_ATOMIZATION
    da_h2 = h2[1]
    return FreezeResult(
        k=k,
        t_plate=t,
        rho_plate=rho,
        kt_plate_ev=temp * K_B / EV,
        ion_share_start=e_ion0 / (heat / m),
        ion_share_frozen=frozen,
        da_ion_min=da_min,
        t_h2_window=h2[0],
        da_h2=da_h2,
        h2_margin_decades=math.log10(da_h2) if da_h2 > 0 else math.nan,
        h2_extrapolated=rates.H_ATOMIC_H.extrapolated(h2[0]) if not math.isnan(h2[0]) else True,
        bond_share=bond,
        velocity_loss_ions=1.0 - math.sqrt(1.0 - frozen),
        velocity_loss_bonds_carbon_frozen=1.0 - math.sqrt(max(1.0 - bond * carbon_part, 0.0)),
        velocity_loss_bonds_all_frozen=1.0 - math.sqrt(max(1.0 - bond, 0.0)),
        velocity_loss_worst=1.0 - math.sqrt(max(1.0 - frozen - bond, 0.0)),
    )


def main() -> None:
    tab = charge.ch2_table()
    rows = []
    for k in (20.0, 30.0, 40.0, 60.0, 80.0, 100.0):
        f = freeze(k, tab)
        rows.append(asdict(f))
        print(
            f"  k {f.k:3.0f}: plate at {1e6 * f.t_plate:5.1f} us, rho {f.rho_plate:.2e} kg/m^3, "
            f"kT {f.kt_plate_ev:4.2f} eV | ions {f.ion_share_start:.1%} -> frozen "
            f"{f.ion_share_frozen:.2%} (Da_min {f.da_ion_min:.1e}) | "
            + (
                "H2 window not reached before the plate"
                if math.isnan(f.da_h2)
                else f"H2 at {f.t_h2_window:.0f} K Da {f.da_h2:.1e}"
            )
            + " | velocity lost: ions "
            f"{f.velocity_loss_ions:.2%}, bonds C-frozen {f.velocity_loss_bonds_carbon_frozen:.1%}"
            f" / all {f.velocity_loss_bonds_all_frozen:.1%}"
        )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"python: wrote {len(rows)} rows -> {OUT}")


if __name__ == "__main__":
    main()
