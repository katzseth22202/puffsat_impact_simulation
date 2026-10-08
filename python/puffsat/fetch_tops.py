"""Fetch the water gray-opacity table from TOPS (aphysics2.lanl.gov) — provenance regenerator.

The TOPS web form is a **two-stage** submission: posting the request form returns an intermediate
confirmation page whose own form must be submitted to get the results table (a direct POST to
`/submit` 500s). The flow below mirrors pyTOPSScrape (Boudreaux, pypi.org/project/pyTOPSScrape),
which carries T-1's permission for automated queries; the server also 500s transiently under
load, hence the retry loop.

Requires the `fetch` extra (`uv run --extra fetch python -m puffsat.fetch_tops`). The saved HTML
is parsed by `puffsat.tops` when `tables --jupiter --tops <file>` builds the scenario table;
keep the pull verbatim (it is the citable provenance artifact, ADR-0007/ADR-0025).
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import mechanize  # untyped (pyproject per-module override); wrapped here — only bytes leave

TOPS_URL = "https://aphysics2.lanl.gov/apps/"
DEFAULT_OUT = Path("data/tables/tops/tops_water_gray.html")

# Request extents: the OPLIB floor up past the table's 1.2e6 K top (keV), and the Jupiter table's
# full rho range (g/cc), log-spaced at the table's own resolution.
T_LOW_KEV = "0.0005"
T_UP_KEV = "0.125"
RHO_LOW_GCC = "1.0e-7"
RHO_UP_GCC = "3.0e-2"
N_RHO = "48"


def fetch_water_gray(
    timeout_s: float = 120.0,
    mixture: str = "2. h 1. o",
    mixname: str = "water",
    t_up_kev: str = T_UP_KEV,
    rho_low_gcc: str = RHO_LOW_GCC,
    rho_up_gcc: str = RHO_UP_GCC,
    n_rho: str = N_RHO,
    n_groups: str | None = None,
    eg_low_kev: str = "0.001",
    eg_high_kev: str = "0.3",
) -> bytes:
    """One two-stage TOPS submission, gray means. Returns HTML.

    With `n_groups` it asks for multigroup means instead: Rosseland and Planck means over
    `n_groups` log-spaced photon-energy boundaries from `eg_low_kev` to `eg_high_kev` (TOPS
    returns one fewer group than boundaries), after the gray table for each temperature.

    `mixture` is TOPS's number-fraction string; the default is water (2 H : 1 O atomic).
    `t_up_kev` must be one of the form's listed temperatures. Pull one mixture at a time: TOPS
    keeps request state per server session, and two concurrent submissions can both come back
    as whichever mixture reached it last."""
    br = mechanize.Browser()
    br.set_handle_robots(False)  # per-pyTOPSScrape: T-1 permits automated queries
    br.open(TOPS_URL, timeout=timeout_s)
    br.select_form(nr=0)
    br.form["mixture"] = mixture
    br.form["mixname"] = mixname
    br.form.find_control(name="tlow", type="select").get(T_LOW_KEV).selected = True
    br.form.find_control(name="tup", type="select").get(t_up_kev).selected = True
    br.form["rlow"] = rho_low_gcc
    br.form["rup"] = rho_up_gcc
    br.form["nr"] = n_rho
    if n_groups is None:
        br.form.find_control(name="datype").value = ["gray"]
    else:
        br.form.find_control(name="egrid").value = ["range"]
        br.form["egplow"] = eg_low_kev
        br.form["egphigh"] = eg_high_kev
        br.form["ngpengs"] = n_groups
        br.form.find_control(name="espace").value = ["log"]
        br.form.find_control(name="datype").value = ["groups"]
    br.submit()  # stage 1: request form -> confirmation page
    br.select_form(nr=0)
    response = br.submit()  # stage 2: confirmation -> results table
    html = bytes(response.read())
    br.close()
    return html


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch the TOPS water gray-opacity table.")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--attempts", type=int, default=5)
    parser.add_argument("--mixture", default="2. h 1. o", help="TOPS number-fraction string")
    parser.add_argument(
        "--mixname",
        default="water",
        help="letters and digits only: TOPS returns HTTP 500 for a name with an underscore",
    )
    parser.add_argument("--t-up-kev", default=T_UP_KEV, help="a temperature the form lists")
    parser.add_argument("--rho-low-gcc", default=RHO_LOW_GCC)
    parser.add_argument("--rho-up-gcc", default=RHO_UP_GCC)
    parser.add_argument("--n-rho", default=N_RHO)
    parser.add_argument("--n-groups", default=None, help="multigroup boundaries; omit for gray")
    parser.add_argument("--eg-low-kev", default="0.001")
    parser.add_argument("--eg-high-kev", default="0.3")
    args = parser.parse_args()

    for attempt in range(1, args.attempts + 1):
        try:
            html = fetch_water_gray(
                mixture=args.mixture,
                mixname=args.mixname,
                t_up_kev=args.t_up_kev,
                rho_low_gcc=args.rho_low_gcc,
                rho_up_gcc=args.rho_up_gcc,
                n_rho=args.n_rho,
                n_groups=args.n_groups,
                eg_low_kev=args.eg_low_kev,
                eg_high_kev=args.eg_high_kev,
            )
            break
        except Exception as exc:  # broad on purpose — the server 500s transiently; retry them all
            print(f"python: TOPS attempt {attempt}/{args.attempts} failed: {exc}", file=sys.stderr)
            time.sleep(5.0)
    else:
        sys.exit("python: TOPS unreachable — kept the existing pull (if any)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_bytes(html)
    kind = "gray" if args.n_groups is None else f"{args.n_groups}-boundary multigroup"
    print(f"python: wrote TOPS {args.mixname} {kind} pull -> {args.out} ({len(html)} bytes)")


if __name__ == "__main__":
    main()
