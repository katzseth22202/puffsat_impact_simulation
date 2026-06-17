//! Rung-A 1D Lagrangian ideal-gas hydrodynamics.
//!
//! The kernel itself (staggered artificial-viscosity Lagrangian, ADR-0022) is built
//! test-first against analytic solutions. This crate currently exposes the [`riemann`]
//! module — the *exact* Riemann solver for the 1D Euler equations, which serves as the
//! analytic oracle for the Sod shock-tube acceptance test (it is test-only and does **not**
//! double as a flux function; ADR-0022).

pub mod riemann;
