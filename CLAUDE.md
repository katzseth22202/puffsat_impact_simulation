# PuffSat Impact Simulation — Project Guidance

Per-collision study computing the paper's fudge factor `f(v)` — the momentum-transfer
efficiency of a PuffSat gas pulse bouncing off the pusher plate.

- **Glossary (canonical terms):** [`CONTEXT.md`](CONTEXT.md) — read before using domain words.
- **Full design + rationale:** [`puffsat_impact_sim_design.md`](puffsat_impact_sim_design.md).
- **Decisions:** [`docs/adr/`](docs/adr/).

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
  - Both live under the gitignored `data/` directory.

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
```

- `cargo clippy --all-targets --all-features` must be **warning-free** for the `all` group
  (it runs `-D warnings` in `make lint` / CI). `pedantic` is `warn` so it advises without
  blocking — bump individual pedantic lints to `deny` as the code stabilizes if useful.
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
