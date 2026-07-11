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


def fetch_water_gray(timeout_s: float = 120.0) -> bytes:
    """One two-stage TOPS submission for water (2 H : 1 O atomic), gray means. Returns HTML."""
    br = mechanize.Browser()
    br.set_handle_robots(False)  # per-pyTOPSScrape: T-1 permits automated queries
    br.open(TOPS_URL, timeout=timeout_s)
    br.select_form(nr=0)
    br.form["mixture"] = "2. h 1. o"
    br.form["mixname"] = "water"
    br.form.find_control(name="tlow", type="select").get(T_LOW_KEV).selected = True
    br.form.find_control(name="tup", type="select").get(T_UP_KEV).selected = True
    br.form["rlow"] = RHO_LOW_GCC
    br.form["rup"] = RHO_UP_GCC
    br.form["nr"] = N_RHO
    br.form.find_control(name="datype").value = ["gray"]
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
    args = parser.parse_args()

    for attempt in range(1, args.attempts + 1):
        try:
            html = fetch_water_gray()
            break
        except Exception as exc:  # broad on purpose — the server 500s transiently; retry them all
            print(f"python: TOPS attempt {attempt}/{args.attempts} failed: {exc}", file=sys.stderr)
            time.sleep(5.0)
    else:
        sys.exit("python: TOPS unreachable — kept the existing pull (if any)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_bytes(html)
    print(f"python: wrote TOPS water gray pull -> {args.out} ({len(html)} bytes)")


if __name__ == "__main__":
    main()
