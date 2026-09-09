# PuffSat Impact Simulation — Project Guidance

Per-collision study computing the paper's fudge factor `f(v)` — the momentum-transfer
efficiency of a PuffSat gas pulse bouncing off the pusher plate.

- **Glossary (canonical terms):** [`CONTEXT.md`](CONTEXT.md) — read before using domain words.
- **Full design + rationale:** [`puffsat_impact_sim_design.md`](puffsat_impact_sim_design.md).
- **Decisions:** [`docs/adr/`](docs/adr/).

The **water-injected overtake plate study** is separate from both `f(v)` and the
head-on tamper study. Scope and current results:
[`docs/water_injected_overtake_plate_study.md`](docs/water_injected_overtake_plate_study.md).
Its Python is `python/puffsat/water_plate/`, outputs are under
`data/results/water_plate/`, and targets are `make water-plate-ledger`,
`make water-plate-profiles`, `make water-plate-thermo`, `make water-plate-chemistry`,
`make water-plate-flow-table`, `make water-plate-flow`, `make water-plate-flow-analysis`,
`make water-plate-flow-report`, `make water-plate-freeze-parent`,
`make water-plate-freeze-prepare`, `make water-plate-freeze`,
`make water-plate-freeze-report`, `make water-plate-molecular-prepare`,
`make water-plate-molecular`, `make water-plate-molecular-report`,
`make water-plate-selective-report`, `make water-plate-kinetics-audit`,
`make water-plate-pressure-kinetics-audit`,
`make water-plate-thermal-kinetics-audit`,
`make water-plate-parcels`, `make water-plate-parcel-report`,
`make water-plate-parcel-check`,
and `make water-plate-test`.
The same-state reaction-switch comparison is scoped in
[ADR-0044](docs/adr/0044-water-plate-parcel-frozen-restarts.md).
Its ordered molecular/ionization separation is scoped in
[ADR-0045](docs/adr/0045-water-plate-selective-chemistry.md).
The local neutral-network rate/applicability screen is scoped in
[ADR-0046](docs/adr/0046-water-plate-kinetics-applicability.md).
The direct water-association pressure-source coverage and channel-removal
sensitivity are scoped in
[ADR-0047](docs/adr/0047-water-plate-pressure-kinetics.md).
Rate-control sensitivities and fixed-energy local temperature feedback are
scoped in [ADR-0048](docs/adr/0048-water-plate-local-thermal-kinetics.md).
Nonlinear neutral chemistry on prescribed parcel density histories is scoped
in [ADR-0049](docs/adr/0049-water-plate-nonlinear-parcels.md). It solves temperature
from energy conservation, but does not recompute hydrodynamics or wall impulse.
The flow reference is planar, stratified and
prepositioned, not a distributed injector; its provisional cold/dense EOS is
scoped in [ADR-0043](docs/adr/0043-water-plate-cold-source-and-provisional-dense-eos.md).
Its incoming momentum is a credit, and its injected
water ratio and geometry must not inherit head-on-study conventions.

## Architecture (settled)

Two languages, one repo. The split is ADR-0002.

- **Rust** — the hot path: 1D Lagrangian rad-hydro kernel, 2D axisymmetric Euler kernel,
  shared table loader, and the rayon sweep driver. Organized as a **cargo workspace**
  (`crates/`).
- **Python** — the cold path only: EOS/opacity table generation (CoolProp / CEA / Saha)
  and all analysis, frontier extraction, and plotting. Managed with **uv**
  (`pyproject.toml`, `python/puffsat/`).
- **Boundary = plain-text files**, no FFI/PyO3 and no binary-format deps (ADR-0002; the
  scale does not justify binary — ADR-0019):
  - **JSON** for the gridded EOS/opacity tables (ADR-0007, amended): one object with grids,
    flattened field arrays, and nested provenance. Python writes (`json`), Rust reads the
    whole table into RAM once (`serde_json`) and interpolates in-memory. The loader lives
    behind a swappable module; revisit only if a production table exceeds ~25 MB / ~1s parse.
  - **JSONL** for sweep results (one JSON object per run, appendable + crash-resilient):
    Rust writes (`serde_json`), Python reads (`pandas.read_json(lines=True)`).
  - Both live under the gitignored `data/` directory — except the small deliverable
    artifacts `CONCLUSION.md` cites as evidence (`data/results/frontier*.csv` and the
    result `*.png` figures), which are committed in place; the bulk tables/JSONL are not
    (ADR-0025).

## Build

A top-level **Makefile** is the single entry point (ADR-0018); it delegates to `cargo`
and `uv` and tracks the `tables → results → analysis` file-dependency graph.

```
make smoke      # boundary round-trip plumbing test (Python -> JSON -> Rust -> JSONL -> Python)
make tables     # generate EOS/opacity tables (uv)
make build      # cargo build --release
make sweep      # run the sweep (cargo, rayon); depends on tables
make analysis   # frontier extraction + plots (uv); depends on sweep
make test       # cargo test + pytest
make lint       # ruff + mypy + cargo clippy + cargo fmt --check
make fmt        # ruff format + cargo fmt
```

**The `tamper-*` targets belong to the other study.** `puffsat_tamper_isp_prd.md` (effective Isp
for a tamped nozzle) is a *separate study* sharing this repo's kernels, tables, and validation
discipline — not its deliverable, regime, or plate-side conventions. Its targets are prefixed
`tamper-`, its artifacts live under `data/results/tamper/`, and its Python is its own package
(`python/puffsat/tamper/`). Departures from decisions recorded for the `f(v)` study are in PRD §12
and ADR-0030..0034.

```
make tamper-ledger   # Rung 0: the analytic reference ledger (owns every closed-form number)
make tamper-test     # that study's tests alone
```

**The `walled-nozzle-*` targets belong to a third study.** The walled thermal nozzle answers asks
N9-N11 of the companion repository against its ADR-0016: a 3 m bore, **no field**, a 200-673 m^3
chamber and 25 kg at 75 km/s into methane. It shares the EOS machinery and the Bray criterion and
nothing else -- not the `f(v)` bag, not the magnetic nozzle's geometry, not the plate-side
conventions. Scope is [ADR-0050](docs/adr/0050-walled-thermal-nozzle-methane-chamber.md); results
are in [`docs/walled_nozzle_asks_answered.md`](docs/walled_nozzle_asks_answered.md). Its Python is
`python/puffsat/walled_nozzle/` plus the top-level `python/puffsat/eos_methane.py`, and its
outputs are under `data/results/walled_nozzle/`.

```
make walled-nozzle-chamber    # N10.4b + N9.0: the solved chamber charge and the sealed-vessel clock
make walled-nozzle-hydrogen   # N10: the Project 242 validation and the pure-hydrogen rung
make walled-nozzle-test       # that study's tests alone
```

## Rust coding standards

**Safe Rust only.** The workspace denies `unsafe_code`. Use safe Rust exclusively.

- If a task seems to genuinely require `unsafe` (FFI, raw-pointer perf, etc.), **stop and
  ask Seth first.** Do not introduce `unsafe` unilaterally.
- If approved, the `unsafe` block must carry an explicit `#[allow(unsafe_code)]` on the
  narrowest possible scope plus a `// SAFETY:` comment justifying every invariant relied
  on. This is `deny`, not `forbid`, precisely so this reviewed escape hatch exists — but it
  is never taken without sign-off.

**Lints.** Workspace-wide via `[workspace.lints]`; each crate opts in with
`[lints] workspace = true`.

```toml
[workspace.lints.rust]
unsafe_code = "deny"          # safe Rust only; opt-in requires sign-off (see above)
missing_debug_implementations = "warn"
rust_2018_idioms = "warn"

[workspace.lints.clippy]
all      = { level = "deny",  priority = -1 }   # correctness/style: must be clean
pedantic = { level = "warn",  priority = -1 }   # guidance, not a gate (noisy on numerics)
cast_precision_loss      = "allow"              # index/count → f64 is exact below 2^53; false positive here
cast_possible_truncation = "deny"               # f64 → int can corrupt an index; gate it
cast_sign_loss           = "deny"               # ditto (provably-safe sites carry a narrow #[allow])
```

- `cargo clippy --all-targets --all-features` must be clean for the `all` group, which
  `[workspace.lints]` sets to `deny` — so a plain `cargo clippy` (what `make lint` / CI runs)
  already fails the build on any `all`-group lint. Do **not** add `-D warnings`: it would
  escalate `pedantic` too. `pedantic` is `warn` so it advises without blocking; bump or relax
  individual pedantic lints (priority `0`, so they override the group) as the code stabilizes.
- **Cast lints (settled).** `cast_precision_loss` is **allowed**: every site is a grid index or
  cell count widened to `f64` (`i as f64 * dx`), which is *exact* below 2⁵³ (~9e15) for our ~10³
  grids and has no float-preserving alternative — it was a false positive burying real warnings.
  It guards integer representability, **not** the `f64` arithmetic that governs solver accuracy
  (that is the job of the order-of-accuracy/convergence tests). Conversely the `f64 → int` casts
  (`cast_possible_truncation`, `cast_sign_loss`) are **denied**, since they can silently corrupt
  an index; the rare provably-safe site takes a narrow `#[allow(...)]` with a `// SAFE:` comment.
- `cargo fmt` is mandatory (checked in `make lint`).
- This is a correctness-critical from-scratch solver: prefer clarity and explicit numeric
  types over cleverness. Document the physics/equation a function implements.

## Python coding standards

- **uv** for environments and running (`uv run`, `uv sync`). Never call a bare `python`.
- **ruff** for lint + format (`make lint` / `make fmt`).
- **mypy strict.** `strict = true` (plus `warn_unreachable`, `warn_redundant_casts`).
  Untyped scientific libraries (CoolProp, Cantera, CEA/rocketcea) are isolated behind
  **per-module** `ignore_missing_imports` overrides — **never a global one**, so any new
  untyped import is still flagged. Wrap every untyped-library call in a small typed function
  at the boundary so `Any` cannot leak into the codebase (cast the return explicitly,
  e.g. `return float(PropsSI(...))`).
- `numpy` and `pandas` ship types — keep them clean under strict.

## Testing (TDD via analytic + convergence)

Kernels are built **test-first** (red-green-refactor), but the "unit" is *a kernel reproducing
a known solution to tolerance*, not micro-functions:

- Write the **analytic-solution acceptance test first** — it is the rung's exit criterion and
  the answer is known before the solver exists: Sod / Noh / Sedov for the Euler track, a
  Marshak wave for flux-limited diffusion, and the `f -> 1` (elastic) / `f -> 0.5` (stick)
  momentum limits (design §9–10).
- Add an **order-of-accuracy test**: the error norm must shrink at the scheme's formal rate
  under grid refinement — the strongest single correctness signal for a solver.
- Pure helpers with a closed form (Riemann flux, table interpolation, Thomas solve) get unit
  tests too, but don't force ceremony on exploratory numerics.

Run with `make test` (`cargo test`; `pytest` once Python tests exist).

## Workflow

- Read relevant files (and the relevant ADR) before changing anything.
- Run `make test` after changes; run `make lint` before committing.
- Follow existing structure and the canonical terms in `CONTEXT.md`.
- Keep `CONTEXT.md` a glossary only — implementation decisions go in `docs/adr/`.

## Scratch space (`todos/`)

The gitignored `todos/` directory is disposable scratch space for things we're working on
together — task checklists, working notes, intermediate plans. Claude may freely create files
there to track or stage a single task. Nothing in it is checked in or load-bearing; clean it up
when the task is done.
