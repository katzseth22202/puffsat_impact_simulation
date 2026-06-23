//! Conserved/primitive fluid state and the γ-law gas for the 2D axisymmetric Euler kernel.
//!
//! The conserved vector per cell is `U = [ρ, ρu_z, ρu_r, E]`, with `E` the **total** energy
//! density (internal + kinetic). The gas is a calibrated effective-γ ideal gas (ADR-0008):
//!
//! ```text
//! p = (γ − 1) ρ e_int,    e_int = E/ρ − ½(u_z² + u_r²),    c = √(γ p / ρ).
//! ```
//!
//! `eta_capture` is a lossless 2D/1D ratio (ADR-0003), so the same γ divides out of numerator and
//! denominator — the kernel never needs the full equilibrium EOS for the geometry sweep.

/// Conserved state: density, axial & radial momentum density, total energy density.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cons {
    /// Mass density `ρ`.
    pub rho: f64,
    /// Axial momentum density `ρu_z`.
    pub mz: f64,
    /// Radial momentum density `ρu_r`.
    pub mr: f64,
    /// Total energy density `E` (internal + kinetic).
    pub e_tot: f64,
}

/// Primitive state: density, axial & radial velocity, pressure.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Prim {
    /// Mass density `ρ`.
    pub rho: f64,
    /// Axial velocity `u_z`.
    pub uz: f64,
    /// Radial velocity `u_r`.
    pub ur: f64,
    /// Pressure `p`.
    pub p: f64,
}

impl Cons {
    /// Construct a conserved state from its components.
    #[must_use]
    pub fn new(rho: f64, mz: f64, mr: f64, e_tot: f64) -> Self {
        Self { rho, mz, mr, e_tot }
    }

    /// Conserved state for a primitive `(ρ, u_z, u_r, p)` under the γ-law (`E = p/(γ−1) + ½ρ|u|²`).
    #[must_use]
    pub fn from_prim(w: Prim, gamma: f64) -> Self {
        let ke = 0.5 * w.rho * (w.uz * w.uz + w.ur * w.ur);
        Self {
            rho: w.rho,
            mz: w.rho * w.uz,
            mr: w.rho * w.ur,
            e_tot: w.p / (gamma - 1.0) + ke,
        }
    }
}

impl Prim {
    /// Construct a primitive state from its components.
    #[must_use]
    pub fn new(rho: f64, uz: f64, ur: f64, p: f64) -> Self {
        Self { rho, uz, ur, p }
    }

    /// Primitive state from a conserved `U` under the γ-law.
    #[must_use]
    pub fn from_cons(u: Cons, gamma: f64) -> Self {
        let uz = u.mz / u.rho;
        let ur = u.mr / u.rho;
        let ke = 0.5 * u.rho * (uz * uz + ur * ur);
        let p = (gamma - 1.0) * (u.e_tot - ke);
        Self {
            rho: u.rho,
            uz,
            ur,
            p,
        }
    }

    /// Adiabatic sound speed `c = √(γ p / ρ)`.
    #[must_use]
    pub fn sound_speed(&self, gamma: f64) -> f64 {
        (gamma * self.p / self.rho).max(0.0).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::{Cons, Prim};
    use approx::assert_relative_eq;

    const GAMMA: f64 = 1.4;

    #[test]
    fn prim_cons_roundtrip() {
        let w = Prim::new(1.3, -0.7, 0.4, 2.1);
        let back = Prim::from_cons(Cons::from_prim(w, GAMMA), GAMMA);
        assert_relative_eq!(back.rho, w.rho, epsilon = 1e-14);
        assert_relative_eq!(back.uz, w.uz, epsilon = 1e-14);
        assert_relative_eq!(back.ur, w.ur, epsilon = 1e-14);
        assert_relative_eq!(back.p, w.p, epsilon = 1e-14);
    }

    #[test]
    fn sound_speed_is_gamma_p_over_rho() {
        let w = Prim::new(2.0, 0.0, 0.0, 5.0);
        assert_relative_eq!(
            w.sound_speed(GAMMA),
            (GAMMA * 5.0 / 2.0).sqrt(),
            epsilon = 1e-14
        );
    }
}
