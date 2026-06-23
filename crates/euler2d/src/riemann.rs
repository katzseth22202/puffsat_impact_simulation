//! HLLC approximate Riemann solver for the directional 1D Euler problem (Toro §10.4), used as the
//! Godunov interface flux for each dimensional sweep (ADR-0023).
//!
//! Dimensional splitting reduces each face to a 1D Riemann problem in the sweep-**normal**
//! direction; the **transverse** velocity is passively advected (it rides along with the contact).
//! So the solver works on a directional primitive `(ρ, u_n, u_t, p)` and returns the flux of
//! `[ρ, ρu_n, ρu_t, E]`. HLLC (not HLL) resolves the contact discontinuity — necessary because the
//! cloud/vacuum interface of the bounce problem *is* a contact.
//!
//! Wave-speed estimate: the **Davis** bounds `S_L = min(u_nL−c_L, u_nR−c_R)`,
//! `S_R = max(u_nL+c_L, u_nR+c_R)` (Toro §10.5.1) — simple, positively conservative, and adequate
//! for the Sod/Sedov/Noh acceptance suite. The contact speed `S*` and star states follow Toro
//! eqs. 10.37–10.39.

/// Directional primitive at an interface: `u_n` is the sweep-normal velocity, `u_t` transverse.
#[derive(Debug, Clone, Copy)]
pub struct DirState {
    /// Mass density `ρ`.
    pub rho: f64,
    /// Sweep-normal velocity `u_n`.
    pub un: f64,
    /// Transverse velocity `u_t` (passively advected).
    pub ut: f64,
    /// Pressure `p`.
    pub p: f64,
}

/// Flux of the directional conserved vector `[ρ, ρu_n, ρu_t, E]` across the interface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DirFlux {
    /// Mass flux.
    pub rho: f64,
    /// Normal-momentum flux.
    pub mn: f64,
    /// Transverse-momentum flux.
    pub mt: f64,
    /// Total-energy flux.
    pub e: f64,
}

/// Directional **conserved** vector `[ρ, ρu_n, ρu_t, E]` — the currency of the MUSCL reconstruction
/// (slope-limited per component) and the Hancock half-step predictor (kernel.rs).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DirCons {
    /// Mass density `ρ`.
    pub rho: f64,
    /// Normal momentum density `ρu_n`.
    pub mn: f64,
    /// Transverse momentum density `ρu_t`.
    pub mt: f64,
    /// Total energy density `E`.
    pub e: f64,
}

impl DirCons {
    /// Directional conserved vector for a directional primitive state.
    #[must_use]
    pub fn from_state(s: DirState, gamma: f64) -> Self {
        Self {
            rho: s.rho,
            mn: s.rho * s.un,
            mt: s.rho * s.ut,
            e: total_energy(s, gamma),
        }
    }

    /// Directional primitive state for this conserved vector.
    #[must_use]
    pub fn to_state(self, gamma: f64) -> DirState {
        let un = self.mn / self.rho;
        let ut = self.mt / self.rho;
        let p = (gamma - 1.0) * (self.e - 0.5 * self.rho * (un * un + ut * ut));
        DirState {
            rho: self.rho,
            un,
            ut,
            p,
        }
    }

    /// `self + scale · other` (component-wise) — slope/extrapolation arithmetic.
    #[must_use]
    pub fn axpy(self, scale: f64, other: DirCons) -> Self {
        Self {
            rho: self.rho + scale * other.rho,
            mn: self.mn + scale * other.mn,
            mt: self.mt + scale * other.mt,
            e: self.e + scale * other.e,
        }
    }

    /// `self + scale · flux` (component-wise) — the Hancock half-step predictor (`DirFlux` shares
    /// the `[ρ, ρu_n, ρu_t, E]` layout).
    #[must_use]
    pub fn add_flux(self, scale: f64, f: DirFlux) -> Self {
        Self {
            rho: self.rho + scale * f.rho,
            mn: self.mn + scale * f.mn,
            mt: self.mt + scale * f.mt,
            e: self.e + scale * f.e,
        }
    }
}

/// Total energy density `E = p/(γ−1) + ½ρ(u_n² + u_t²)` for a directional state.
fn total_energy(s: DirState, gamma: f64) -> f64 {
    s.p / (gamma - 1.0) + 0.5 * s.rho * (s.un * s.un + s.ut * s.ut)
}

/// The physical Euler flux of `[ρ, ρu_n, ρu_t, E]` for a single state.
pub(crate) fn phys_flux(s: DirState, gamma: f64) -> DirFlux {
    let e = total_energy(s, gamma);
    DirFlux {
        rho: s.rho * s.un,
        mn: s.rho * s.un * s.un + s.p,
        mt: s.rho * s.un * s.ut,
        e: (e + s.p) * s.un,
    }
}

/// HLLC star-region contribution `F_K + S_K (U*_K − U_K)` on side `K` (Toro eq. 10.38–10.39),
/// where `sk` is that side's outer wave speed and `s_star` the contact speed.
// star_mn/star_mt are the normal/transverse components of one star-state momentum; the near-
// identical names mirror the physics.
#[allow(clippy::similar_names)]
fn star_flux(s: DirState, sk: f64, s_star: f64, gamma: f64) -> DirFlux {
    let e = total_energy(s, gamma);
    let f = phys_flux(s, gamma);
    // U*_K = ρ_K (S_K − u_n)/(S_K − S*) · [1, S*, u_t, E/ρ + (S* − u_n)(S* + p/(ρ(S_K − u_n)))].
    let factor = s.rho * (sk - s.un) / (sk - s_star);
    let star_rho = factor;
    let star_mn = factor * s_star;
    let star_mt = factor * s.ut;
    let star_e = factor * (e / s.rho + (s_star - s.un) * (s_star + s.p / (s.rho * (sk - s.un))));
    DirFlux {
        rho: f.rho + sk * (star_rho - s.rho),
        mn: f.mn + sk * (star_mn - s.rho * s.un),
        mt: f.mt + sk * (star_mt - s.rho * s.ut),
        e: f.e + sk * (star_e - e),
    }
}

/// The HLLC numerical flux across an interface with left/right directional states.
#[must_use]
pub fn hllc_flux(l: DirState, r: DirState, gamma: f64) -> DirFlux {
    let cl = (gamma * l.p / l.rho).max(0.0).sqrt();
    let cr = (gamma * r.p / r.rho).max(0.0).sqrt();

    // Davis wave-speed bounds (Toro §10.5.1).
    let sl = (l.un - cl).min(r.un - cr);
    let sr = (l.un + cl).max(r.un + cr);

    if sl >= 0.0 {
        return phys_flux(l, gamma);
    }
    if sr <= 0.0 {
        return phys_flux(r, gamma);
    }

    // Contact speed S* (Toro eq. 10.37).
    let s_star = (r.p - l.p + l.rho * l.un * (sl - l.un) - r.rho * r.un * (sr - r.un))
        / (l.rho * (sl - l.un) - r.rho * (sr - r.un));

    if s_star >= 0.0 {
        star_flux(l, sl, s_star, gamma)
    } else {
        star_flux(r, sr, s_star, gamma)
    }
}

#[cfg(test)]
mod tests {
    use super::{DirState, hllc_flux, phys_flux};
    use approx::assert_relative_eq;

    const GAMMA: f64 = 1.4;

    fn assert_flux_eq(a: super::DirFlux, b: super::DirFlux, eps: f64) {
        assert_relative_eq!(a.rho, b.rho, epsilon = eps);
        assert_relative_eq!(a.mn, b.mn, epsilon = eps);
        assert_relative_eq!(a.mt, b.mt, epsilon = eps);
        assert_relative_eq!(a.e, b.e, epsilon = eps);
    }

    /// Consistency: equal left/right states ⇒ the HLLC flux is the exact physical flux.
    #[test]
    fn equal_states_give_physical_flux() {
        let s = DirState {
            rho: 1.2,
            un: 0.3,
            ut: -0.5,
            p: 0.9,
        };
        assert_flux_eq(hllc_flux(s, s, GAMMA), phys_flux(s, GAMMA), 1e-13);
    }

    /// Supersonic to the right (both wave bounds positive) ⇒ pure left-upwind flux.
    #[test]
    fn supersonic_right_is_left_upwind() {
        let l = DirState {
            rho: 1.0,
            un: 5.0,
            ut: 0.0,
            p: 1.0,
        };
        let r = DirState {
            rho: 0.5,
            un: 5.0,
            ut: 0.0,
            p: 0.5,
        };
        assert_flux_eq(hllc_flux(l, r, GAMMA), phys_flux(l, GAMMA), 1e-13);
    }

    /// Supersonic to the left (both wave bounds negative) ⇒ pure right-upwind flux.
    #[test]
    fn supersonic_left_is_right_upwind() {
        let l = DirState {
            rho: 1.0,
            un: -5.0,
            ut: 0.0,
            p: 1.0,
        };
        let r = DirState {
            rho: 0.5,
            un: -5.0,
            ut: 0.0,
            p: 0.5,
        };
        assert_flux_eq(hllc_flux(l, r, GAMMA), phys_flux(r, GAMMA), 1e-13);
    }

    /// The transverse velocity is passively advected: it leaves the mass and normal-momentum
    /// fluxes unchanged, and the transverse-momentum flux is the mass flux times the upwind `u_t`.
    /// These states put the contact moving right (`S* > 0`), so the left value is advected.
    #[test]
    fn transverse_velocity_is_passive() {
        let base_l = DirState {
            rho: 1.0,
            un: 0.4,
            ut: 0.0,
            p: 1.0,
        };
        let base_r = DirState {
            rho: 0.3,
            un: 0.2,
            ut: 0.0,
            p: 0.4,
        };
        let f0 = hllc_flux(base_l, base_r, GAMMA);
        let l = DirState { ut: 2.0, ..base_l };
        let r = DirState { ut: 7.0, ..base_r };
        let f = hllc_flux(l, r, GAMMA);
        // Mass and normal-momentum fluxes are independent of the transverse velocity.
        assert_relative_eq!(f.rho, f0.rho, epsilon = 1e-13);
        assert_relative_eq!(f.mn, f0.mn, epsilon = 1e-13);
        // Contact moves right ⇒ left transverse value advected: F_mt = F_rho · u_t,L.
        assert!(f.rho > 0.0);
        assert_relative_eq!(f.mt, f.rho * l.ut, epsilon = 1e-13);
    }
}
