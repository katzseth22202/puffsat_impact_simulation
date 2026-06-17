# Monorepo: a cargo workspace and a uv project, orchestrated by one Makefile

The Rust hot path (ADR-0002) is a **cargo workspace** under `crates/` — the 1D rad-hydro
kernel, the 2D Euler kernel, a shared table loader, and the sweep driver as separate crates
so they grow on independent lifecycles. The Python cold path is a sibling **uv project**
(`pyproject.toml`, `python/puffsat/`). The two languages meet only through the gitignored
`data/` directory (ADR-0019). A top-level **Makefile** is the single build entry point,
delegating to `cargo` and `uv` (`make tables | build | sweep | analysis | smoke | test |
lint | fmt`).

Make — not `just` or `cargo-xtask` — because the pipeline **is** a file-dependency DAG:
`data/tables/*.json` (built by uv) → `data/results/*.jsonl` (built by cargo) → analysis and
plots (uv). Make's timestamp model expresses that natively, so `make analysis` will not
rebuild a multi-hour sweep whose inputs are unchanged. Make is also already installed, adding
no new dependency.

## Considered Options

- **just** — cleaner command-runner ergonomics, but a *task runner* with no file-timestamp
  dependency tracking (it always re-runs targets) and an extra install. Rejected: loses the
  incremental-rebuild property that is the whole reason to use a build tool here.
- **cargo-xtask** — keeps orchestration in one Rust toolchain, but driving `uv`/Python and
  expressing a cross-language file DAG from a Rust binary is clumsier and heavier to write.
  Rejected.
- **Single Rust crate, split later** — rejected: the 1D kernel, 2D kernel, and table loader
  have distinct lifecycles and the sweep driver binary links the kernels; a workspace keeps
  those boundaries clean from rung A, avoiding a later refactor.
