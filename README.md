# PuffSat Impact Simulation

A physics simulation to compute the per-collision coefficient of restitution of PuffSat gas against a pusher plate — the fudge factor `f` in the paper [*Aim Is All You Need: A Speculative White Paper on PuffSat Pulsed Propulsion*](https://github.com/katzseth22202/Balloon-Pulse-Propulsion).

## What This Computes

The core output is `f(v)`, the fudge factor as a function of impact speed across the mission envelope (3.2–16 km/s):

```
f = eta_capture * (1 + e_eff) / 2
```

- `eta_capture` — fraction of axial momentum that lands and rebounds usefully (geometry)
- `e_eff` — effective restitution after radiative, conductive, and condensation losses (thermophysics)

This backs paper §3.2, which currently assumes a constant `f = 0.8`. The sim replaces that with a defensible `f(v)` curve and a full loss budget.

## Architecture

Two physics tracks feed a combined result:

| Track | Tool | Output |
|---|---|---|
| Restitution + wall load | 1D Lagrangian rad-hydro (Rust) | `e_eff`, peak flux, peak pressure |
| Geometric capture | 2D axisymmetric Euler, radiation-free (Rust) | `eta_capture` vs plate shape |
| Confirmation | FLASH (deferred) | independent cross-check |

Python handles EOS/opacity table generation (Cantera/CoolProp/CEA) and all analysis/plotting. The boundary between Rust and Python is a plain-text file format — JSON for the EOS/opacity tables, JSONL for the sweep results (ADR-0019) — exchanged through the gitignored `data/` directory. No FFI.

## Getting Started

The repo is a cargo workspace (Rust hot path) plus a [uv](https://docs.astral.sh/uv/)-managed Python project (cold path), glued by a `Makefile` (ADR-0018).

### Prerequisites

- **Rust** (stable) via [rustup](https://rustup.rs) — `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh`
- **A C linker** — Rust links through the system `cc`. On Debian/Ubuntu: `sudo apt-get install -y gcc` (or `build-essential`).
- **uv** (manages the Python 3.12+ interpreter and the virtualenv) — `curl -LsSf https://astral.sh/uv/install.sh | sh`

### Set up and verify

```bash
git clone <repo-url> && cd puffsat_impact_simulation
uv sync          # create .venv and install dev tools (ruff, mypy)
make smoke       # end-to-end plumbing test; prints "SMOKE OK"
```

`make smoke` runs the full cross-language round-trip — Python writes a JSON table, the Rust `smoke` crate reads it and appends a JSONL result, Python reads it back and asserts the value — so a green run confirms both toolchains, the workspace build, and the file boundary together.

### Common tasks

| Command | What it does |
|---|---|
| `make smoke` | Plumbing round-trip test (Python ↔ Rust via `data/`) |
| `make build` | Compile the Rust workspace |
| `make test` | Run all tests (`cargo test`; pytest later) |
| `make lint` | ruff + mypy + clippy + fmt checks (the CI gate) |
| `make fmt` | Auto-format Python and Rust |
| `make tables` / `make sweep` / `make analysis` | Physics pipeline (stubs until build rungs B+) |

All generated tables and results land in the gitignored `data/` directory.

## Build Order

- **A.** 1D ideal-gas smoke test — Sod shock tube, elastic/inelastic momentum limits
- **B.** 1D high-v package — equilibrium EOS, real opacity tables, flux-limited diffusion at 16 km/s
- **C.** 1D low-v package — cool-gas, optically thin, condensation loss; worst case is water at 3.2 km/s
- **D.** 2D Euler geometry — `eta_capture` for flat and shallow-concave plates
- **E.** Ablating wall and levers — transpiring wall, dark-oil opacity seed study
- **F.** Validation — FLASH cross-check and Project Orion impulse/ablation reproduction

## Reference

- The paper this simulation backs: [*Aim Is All You Need: A Speculative White Paper on PuffSat Pulsed Propulsion*](https://github.com/katzseth22202/Balloon-Pulse-Propulsion)
- Full physics rationale, gas model, plate construction, design variables, and validation plan: [`puffsat_impact_sim_design.md`](puffsat_impact_sim_design.md)
