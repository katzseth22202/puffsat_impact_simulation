# Rust for the solver core, Python only at the table/analysis boundary

The 1D rad-hydro kernel, the 2D Euler kernel, and the sweep driver are written in **Rust**.
**Python** is used only off the hot path: generating EOS/opacity tables (Cantera / CoolProp /
CEA) and doing all analysis, frontier extraction, and plotting. The boundary is a file format
(HDF5 / Parquet), no FFI or PyO3.

Rust is chosen for **compile-time memory- and data-race safety in a from-scratch,
correctness-critical solver**, plus `cargo`/`rayon` ergonomics. It is *not* chosen for the
parameter-sweep parallelism: the sweep is embarrassingly parallel (independent runs), which is
process-level and language-agnostic — C++, Rust, even Python `multiprocessing` all do it equally
well. So the multithreading argument is not the deciding factor; defensible correctness is.

## Considered Options

- **C/C++.** Its decisive advantage is the mature implicit-solver ecosystem
  (PETSc / HYPRE / SUNDIALS) for stiff multi-dimensional linear systems. **The factored
  architecture neutralizes that advantage:** the 1D radiation flux-limited-diffusion step is a
  *tridiagonal* solve (Thomas algorithm, no library), and the 2D track is *explicit* Godunov
  Euler (no implicit solve at all). So C++ would bring build/dependency overhead (CMake / Conan /
  vcpkg) for almost none of its ecosystem benefit. **Switch to C++ only if** the plan changes to
  heavy implicit multi-D solves, extending an existing C++ hydro code, or tight MPI cluster
  scaling — none of which are current.

## Stack

`ndarray` (arrays), `faer` (dense/sparse linear algebra if ever needed), `rayon` (sweep
data-parallelism), `hdf5` (the Python↔Rust table boundary). A NumPy prototype of the 1D kernel
for rung A remains acceptable — debug the physics before fighting the borrow checker, then port.
