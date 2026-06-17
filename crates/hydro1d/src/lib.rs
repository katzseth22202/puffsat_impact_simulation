//! Rung-A 1D Lagrangian ideal-gas hydrodynamics.
//!
//! - [`kernel`] — the solver: a staggered-grid Lagrangian scheme with von Neumann–Richtmyer
//!   artificial viscosity (ADR-0022), built test-first against analytic solutions. Generic over
//!   the [`eos`] it carries.
//! - [`eos`] — the equation of state the kernel calls (`p(ρ, e)`; ADR-0022): rung A's analytic
//!   [`eos::IdealGas`] and rung B's tabulated [`eos::TableEos`] (ADR-0007).
//! - [`riemann`] — the *exact* Riemann solver for the 1D Euler equations, used as the analytic
//!   oracle for the Sod shock-tube acceptance test. It is test-only and does **not** double as
//!   a flux function (the AV kernel has no Riemann solver; ADR-0022).

pub mod eos;
pub mod kernel;
pub mod radiation;
pub mod riemann;

/// A primitive fluid state `(ρ, u, p)` — the shared currency between the kernel and the oracle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Primitive {
    /// Mass density `ρ`.
    pub rho: f64,
    /// Velocity `u`.
    pub u: f64,
    /// Pressure `p`.
    pub p: f64,
}

impl Primitive {
    /// Construct a primitive state.
    #[must_use]
    pub fn new(rho: f64, u: f64, p: f64) -> Self {
        Self { rho, u, p }
    }
}
