"""Boundary plumbing smoke test (the Python ends of the round-trip).

`write` emits a 2x2 toy `(rho, T) -> p` JSON table; the Rust `smoke` crate interpolates it and
appends a JSONL row; `check` reads that row back and asserts the value. Driven by `make smoke`.
Stdlib-only on purpose — this validates the JSON/JSONL contract (ADR-0019), not the science.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

TABLE_PATH = Path("data/tables/smoke.json")
RESULT_PATH = Path("data/results/smoke.jsonl")
EXPECTED_P = 15.0  # bilinear midpoint of corners [0, 10, 20, 30]
TOL = 1e-9


def write_table() -> None:
    """Write the toy table; `p` is flattened row-major over `(rho, T)`."""
    table = {
        "shape": [2, 2],
        "rho_grid": [0.0, 1.0],
        "T_grid": [0.0, 1.0],
        "p": [0.0, 10.0, 20.0, 30.0],
        "provenance": {"source": "smoke-test toy table", "generated_by": "puffsat.smoke"},
    }
    TABLE_PATH.parent.mkdir(parents=True, exist_ok=True)
    TABLE_PATH.write_text(json.dumps(table, indent=2))
    print(f"python: wrote toy table -> {TABLE_PATH}")


def check_result() -> None:
    """Read the last JSONL record written by Rust and assert the interpolated value."""
    if not RESULT_PATH.exists():
        sys.exit(f"FAIL: {RESULT_PATH} missing (did the Rust step run?)")
    lines = RESULT_PATH.read_text().splitlines()
    if not lines:
        sys.exit(f"FAIL: {RESULT_PATH} is empty")
    p = float(json.loads(lines[-1])["p_interp"])
    if abs(p - EXPECTED_P) > TOL:
        sys.exit(f"FAIL: p_interp={p}, expected {EXPECTED_P}")
    print(f"python: read p_interp={p} (expected {EXPECTED_P}) -> SMOKE OK")


def main(argv: list[str]) -> None:
    if len(argv) != 2 or argv[1] not in {"write", "check"}:
        sys.exit("usage: smoke.py [write|check]")
    (write_table if argv[1] == "write" else check_result)()


if __name__ == "__main__":
    main(sys.argv)
