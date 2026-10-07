"""The near-term chamber's volume trade: a bigger chamber softens the wall's spike, at what Isp?

The blast solves (`chamber_blast.py`, `vessel_blast.py`) put the reflected spike at the wall at a
roughly fixed multiple of the quasi-static pressure, independent of size. At a fixed pulse
energy the spike therefore falls as `1/V`. This module prices that on the nozzle side, with
`near_term.py`'s solved physics: the real equilibrium gas, the nozzle isentrope, the choked
blowdown and the `H + H + M` freeze clock. Rod, pairings and chamber temperature are fixed.

At each volume the throat grows with the chamber, holding `V/A*` at the near-term 134 m, so
every chamber empties on the same clock. The chamber temperature is set by energy per kilogram,
so it does not move. What moves is the pressure, and with it the recombination race in the
nozzle.

Spike scaling, from the converged runs:
- centred blast, 20 m^3: ~1.8 GPa methane, ~2.95 GPa hydrogen (1-D sphere extrapolated);
- the long plug (10 kg over 2.5 m, back-weighted): 0.61 GPa worst at 1 cm cells in the methane
  chamber, taken as ~1-1.5 GPa converged.

Both are scaled as `20/V`. The wall-area factor `(V/20)^(2/3)` is what fixed-thickness layers
(liner, minimum steel gauge, skins) grow by; the pressure-carrying mass is constant at fixed `pV`.
"""

from __future__ import annotations

import csv
import dataclasses
import math
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

from puffsat.walled_nozzle import near_term as nt

OUTPUT = Path("data/results/walled_nozzle/near_term/volume_trade.csv")
VOLUMES = (20.0, 40.0, 60.0, 80.0, 120.0, 160.0, 240.0)
V_PER_THROAT = nt.VOLUME / nt.THROAT_AREAS[0]
AREA_RATIOS = (100.0, 300.0)
SPIKE_20 = {"hydrogen_5500": (2.95e9, math.nan), "methane_7000": (1.8e9, 1.25e9)}
"""(centred, long plug) converged spike at 20 m^3 [Pa]; the long plug is methane only."""


@dataclass(frozen=True)
class Row:
    pairing: str
    volume_m3: float
    charge_kg: float
    pressure_bar: float
    throat_m2: float
    blowdown_efold_ms: float
    area_ratio: float
    exit_diameter_m: float
    eta_blowdown: float
    eta_net: float
    store_at_risk: float
    isp_effective_s: float
    spike_centred_gpa: float
    spike_long_plug_gpa: float
    wall_area_factor: float


def efold_time(sol: nt.ChamberSolution, br: nt.Branch, tab: nt.Isentrope, throat: float) -> float:
    """Time for the chamber to fall to 1/e of its mass [s]."""
    rows = nt.blowdown(sol, br, tab, throat)
    for prev, row in pairwise(rows):
        if row.mass_fraction <= math.exp(-1.0):
            f = (prev.mass_fraction - math.exp(-1.0)) / (prev.mass_fraction - row.mass_fraction)
            return prev.time + f * (row.time - prev.time)
    return math.nan


def run(report: nt.Report = print) -> list[Row]:
    out: list[Row] = []
    for pairing in (p for p in nt.PAIRINGS if p.name in SPIKE_20):
        for volume in VOLUMES:
            sol = nt.solve_charge(pairing, "janaf", volume)
            br = nt.branch("equilibrium", sol)
            tab = nt.isentrope(br, sol.rho, pairing.temp)
            throat = volume / V_PER_THROAT
            effs = {e.area_ratio: e for e in nt.efficiencies(sol, br, tab, AREA_RATIOS)}
            efold = efold_time(sol, br, tab, throat)
            centred, plug = SPIKE_20[pairing.name]
            for ar in AREA_RATIOS:
                eff = effs[ar]
                fc = dataclasses.replace(
                    nt.freeze_check(sol, tab, throat, ar), eta_blowdown=eff.eta_blowdown
                )
                row = Row(
                    pairing=pairing.name,
                    volume_m3=volume,
                    charge_kg=sol.charge_mass,
                    pressure_bar=sol.pressure / 1e5,
                    throat_m2=throat,
                    blowdown_efold_ms=efold * 1e3,
                    area_ratio=ar,
                    exit_diameter_m=math.sqrt(4.0 * throat * ar / math.pi),
                    eta_blowdown=eff.eta_blowdown,
                    eta_net=fc.eta_net,
                    store_at_risk=fc.store_at_risk,
                    isp_effective_s=eff.isp.isp_effective,
                    spike_centred_gpa=centred * 20.0 / volume / 1e9,
                    spike_long_plug_gpa=plug * 20.0 / volume / 1e9,
                    wall_area_factor=(volume / 20.0) ** (2.0 / 3.0),
                )
                out.append(row)
                report(
                    f"{pairing.name:14s} V={volume:5.0f} p={row.pressure_bar:5.0f} bar "
                    f"A*={throat:5.3f} efold={row.blowdown_efold_ms:5.1f} ms AR={ar:4.0f} "
                    f"D_exit={row.exit_diameter_m:4.1f} m eta_bd={row.eta_blowdown:.3f} "
                    f"eta_net={row.eta_net:.3f} Isp_eq={row.isp_effective_s:5.0f} s "
                    f"spike={row.spike_centred_gpa:.2f}/{row.spike_long_plug_gpa:.2f} GPa "
                    f"area x{row.wall_area_factor:.2f}"
                )
    return out


def main() -> None:
    rows = run()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fields = list(Row.__dataclass_fields__)
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fields)
        for r in rows:
            writer.writerow([getattr(r, f) for f in fields])
    print(f"\nwrote {len(rows)} rows -> {OUTPUT}")


if __name__ == "__main__":
    main()
