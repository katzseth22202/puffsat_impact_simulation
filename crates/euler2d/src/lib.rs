//! 2D axisymmetric Euler kernel for the radiation-free `eta_capture` geometry track (ADR-0003).
//!
//! `f = eta_capture · (1 + e_eff)/2`. The 1D rad-hydro kernel ([`hydro1d`], a separate crate)
//! delivers `e_eff`; this kernel delivers `eta_capture`, the **lossless 2D/1D wall-impulse ratio**
//! (ADR-0003) — pure capture geometry, with the common gas-dynamic re-expansion divided out. The
//! sweep is radiation-free and uses a calibrated effective-γ EOS (ADR-0008).
//!
//! - [`state`] — conserved/primitive state and the γ-law gas.
//! - [`riemann`] — the HLLC interface flux for each dimensional sweep (Toro §10.4).
//! - [`kernel`] — the finite-volume Godunov solver ([`kernel::Grid2D`]), built test-first against
//!   Sod (D0), then Sedov/Noh + order-of-accuracy (D1), then the flat-plate bounce (D2).
//!
//! Numerics are ADR-0023 (the 2D sibling of the 1D kernel's ADR-0022).

pub mod kernel;
pub mod riemann;
pub mod state;
