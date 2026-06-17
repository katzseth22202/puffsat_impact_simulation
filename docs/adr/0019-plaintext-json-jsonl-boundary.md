# Plain-text JSON/JSONL across the Python↔Rust boundary — no binary formats

§8 and ADR-0002 originally specified **HDF5** for the EOS/opacity tables and **Parquet** for
the sweep results. The actual data scale does not justify either, so the boundary is plain
text instead:

- **Tables → JSON.** One object per table: axis grids, flattened field arrays
  (`p, e, c_s, κ_Rosseland, κ_Planck`), and a nested `provenance` object. The table is loaded
  **whole into RAM once** (§8), which negates HDF5's chunked-random-access advantage; at
  single-digit-to-tens of MB the parse is sub-second and paid once. `serde_json` round-trips
  f64 **losslessly**, so there is no precision cost, and nested JSON carries ADR-0007's
  per-regime provenance more naturally than HDF5/safetensors string attributes.
- **Results → JSONL.** One JSON object per run. Sweep output is "hundreds, not thousands" of
  runs (~400 rows × ~25 fields ≈ 60 KB) — far below where Parquet's columnar typing earns its
  dependency. JSONL is **appendable** (each parallel run appends a line, so a crash mid-sweep
  preserves completed runs), greppable, and typed (numbers stay numbers, unlike CSV).

This removes every binary-format dependency — `hdf5`/`safetensors`/`arrow`/`parquet` on the
Rust side, `h5py`/`pyarrow` on the Python side — leaving the boundary as `serde_json` ↔ stdlib
`json`. The table loader stays behind a swappable module.

## Consequences

- Revisit **only the tables** if a production table ever exceeds ~25 MB as JSON or ~1 s to
  parse; the drop-in upgrade is `safetensors` (compact, load-whole, minimal deps) behind the
  unchanged loader interface. Results remain JSONL at this study's scale regardless of count.
- Supersedes the HDF5/Parquet specifics in §8 and ADR-0002's stack note. ADR-0002's core
  decision — Rust core, Python only at the table/analysis boundary, **no FFI/PyO3** — is
  unchanged; the boundary is still a file format, just a plain-text one.

## Considered Options

- **Keep HDF5 + Parquet** — the heaviest dependency set for the smallest data, plus real
  HDF5-from-Rust build friction (the original `hdf5` crate is unmaintained; this environment
  has no system `libhdf5` and no `pkg-config`). Rejected.
- **safetensors tables + JSONL results** — a reasonable hedge *if* tables prove large, but it
  adds a binary dependency now to solve a problem not yet measured. Deferred to the revisit
  trigger above.
