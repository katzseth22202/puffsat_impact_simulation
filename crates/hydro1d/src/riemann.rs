//! Exact Riemann solver for the 1D Euler equations with an ideal-gas (constant-γ) EOS.
//!
//! This is the analytic **oracle** for the Sod shock-tube acceptance test, following Toro,
//! *Riemann Solvers and Numerical Methods for Fluid Dynamics*, 3rd ed., Ch. 4. It is
//! deliberately test-only: the staggered artificial-viscosity kernel (ADR-0022) captures
//! shocks with a `q` term and never solves a Riemann problem, so this solver exists purely
//! to grade the kernel against a known solution.
//!
//! The procedure: find the star-region pressure `p*` by Newton iteration on the pressure
//! function `f(p) = f_L(p) + f_R(p) + (u_R − u_L) = 0` (Toro eq. 4.5), recover the star
//! velocity `u*` (eq. 4.9), then sample the self-similar solution `W(x/t)` (the `sample`
//! routine, Toro §4.5).

/// A primitive fluid state `(ρ, u, p)`.
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

/// Ideal-gas sound speed `c = sqrt(γ p / ρ)`.
fn sound_speed(state: Primitive, gamma: f64) -> f64 {
    (gamma * state.p / state.rho).sqrt()
}

/// One side's pressure function and its derivative, `(f_K(p), f_K'(p))` (Toro eq. 4.6–4.7).
///
/// A shock (`p > p_K`) and a rarefaction (`p ≤ p_K`) take different closed forms; the
/// derivative is supplied for the Newton iteration in [`solve`].
fn pressure_function(p: f64, state: Primitive, gamma: f64) -> (f64, f64) {
    let c_s = sound_speed(state, gamma);
    if p > state.p {
        // Shock branch.
        let a_k = 2.0 / ((gamma + 1.0) * state.rho);
        let b_k = (gamma - 1.0) / (gamma + 1.0) * state.p;
        let root = (a_k / (b_k + p)).sqrt();
        let f = (p - state.p) * root;
        let df = root * (1.0 - (p - state.p) / (2.0 * (b_k + p)));
        (f, df)
    } else {
        // Rarefaction branch.
        let power = (gamma - 1.0) / (2.0 * gamma);
        let f = (2.0 * c_s / (gamma - 1.0)) * ((p / state.p).powf(power) - 1.0);
        let df = (1.0 / (state.rho * c_s)) * (p / state.p).powf(-(gamma + 1.0) / (2.0 * gamma));
        (f, df)
    }
}

/// The solved self-similar structure of a Riemann problem, ready to sample at any `ξ = x/t`.
#[derive(Debug, Clone, Copy)]
pub struct RiemannSolution {
    left: Primitive,
    right: Primitive,
    gamma: f64,
    p_star: f64,
    u_star: f64,
}

impl RiemannSolution {
    /// Star-region pressure `p*` (constant across the contact discontinuity).
    #[must_use]
    pub fn p_star(&self) -> f64 {
        self.p_star
    }

    /// Star-region velocity `u*` (the contact-discontinuity speed).
    #[must_use]
    pub fn u_star(&self) -> f64 {
        self.u_star
    }

    /// Sample the self-similar solution at `ξ = x/t`, returning the primitive state there
    /// (Toro §4.5). `ξ ≤ u*` is governed by the left wave, `ξ > u*` by the right wave; each
    /// wave is a shock or a rarefaction depending on whether `p*` exceeds the far pressure.
    #[must_use]
    pub fn sample(&self, xi: f64) -> Primitive {
        let gamma = self.gamma;
        let gm1 = gamma - 1.0;
        let gp1 = gamma + 1.0;
        let ps = self.p_star;
        let us = self.u_star;

        if xi <= us {
            // Left of the contact: left-going wave.
            let state = self.left;
            let c_s = sound_speed(state, gamma);
            if ps > state.p {
                // Left shock.
                let s_shock = state.u
                    - c_s * ((gp1 / (2.0 * gamma)) * (ps / state.p) + gm1 / (2.0 * gamma)).sqrt();
                if xi < s_shock {
                    state
                } else {
                    let rho = state.rho * ((ps / state.p) + gm1 / gp1)
                        / ((gm1 / gp1) * (ps / state.p) + 1.0);
                    Primitive::new(rho, us, ps)
                }
            } else {
                // Left rarefaction.
                let c_star = c_s * (ps / state.p).powf(gm1 / (2.0 * gamma));
                let s_head = state.u - c_s;
                let s_tail = us - c_star;
                if xi < s_head {
                    state
                } else if xi > s_tail {
                    let rho = state.rho * (ps / state.p).powf(1.0 / gamma);
                    Primitive::new(rho, us, ps)
                } else {
                    let factor = 2.0 / gp1 + (gm1 / (gp1 * c_s)) * (state.u - xi);
                    let rho = state.rho * factor.powf(2.0 / gm1);
                    let vel = 2.0 / gp1 * (c_s + 0.5 * gm1 * state.u + xi);
                    let pres = state.p * factor.powf(2.0 * gamma / gm1);
                    Primitive::new(rho, vel, pres)
                }
            }
        } else {
            // Right of the contact: right-going wave.
            let state = self.right;
            let c_s = sound_speed(state, gamma);
            if ps > state.p {
                // Right shock.
                let s_shock = state.u
                    + c_s * ((gp1 / (2.0 * gamma)) * (ps / state.p) + gm1 / (2.0 * gamma)).sqrt();
                if xi > s_shock {
                    state
                } else {
                    let rho = state.rho * ((ps / state.p) + gm1 / gp1)
                        / ((gm1 / gp1) * (ps / state.p) + 1.0);
                    Primitive::new(rho, us, ps)
                }
            } else {
                // Right rarefaction.
                let c_star = c_s * (ps / state.p).powf(gm1 / (2.0 * gamma));
                let s_head = state.u + c_s;
                let s_tail = us + c_star;
                if xi > s_head {
                    state
                } else if xi < s_tail {
                    let rho = state.rho * (ps / state.p).powf(1.0 / gamma);
                    Primitive::new(rho, us, ps)
                } else {
                    let factor = 2.0 / gp1 - (gm1 / (gp1 * c_s)) * (state.u - xi);
                    let rho = state.rho * factor.powf(2.0 / gm1);
                    let vel = 2.0 / gp1 * (-c_s + 0.5 * gm1 * state.u + xi);
                    let pres = state.p * factor.powf(2.0 * gamma / gm1);
                    Primitive::new(rho, vel, pres)
                }
            }
        }
    }
}

/// Solve the Riemann problem between a `left` and `right` ideal-gas state.
///
/// Finds `p*` by Newton iteration on `f_L(p) + f_R(p) + (u_R − u_L) = 0` (Toro eq. 4.5),
/// seeded with the two-rarefaction/PVRS guess (eq. 4.47), then recovers `u*` (eq. 4.9).
#[must_use]
pub fn solve(left: Primitive, right: Primitive, gamma: f64) -> RiemannSolution {
    let c_l = sound_speed(left, gamma);
    let c_r = sound_speed(right, gamma);
    let du = right.u - left.u;

    // PVRS initial guess (Toro eq. 4.47), floored to stay positive.
    let p_pv = 0.5 * (left.p + right.p) - 0.125 * du * (left.rho + right.rho) * (c_l + c_r);
    let mut p = p_pv.max(1e-6);

    for _ in 0..100 {
        let (f_l, df_l) = pressure_function(p, left, gamma);
        let (f_r, df_r) = pressure_function(p, right, gamma);
        let p_next = (p - (f_l + f_r + du) / (df_l + df_r)).max(1e-12);
        let converged = (p_next - p).abs() / (0.5 * (p_next + p)) < 1e-12;
        p = p_next;
        if converged {
            break;
        }
    }

    let (f_l, _) = pressure_function(p, left, gamma);
    let (f_r, _) = pressure_function(p, right, gamma);
    let u_star = 0.5 * (left.u + right.u) + 0.5 * (f_r - f_l);

    RiemannSolution {
        left,
        right,
        gamma,
        p_star: p,
        u_star,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const GAMMA: f64 = 1.4;

    fn sod() -> (Primitive, Primitive) {
        (
            Primitive::new(1.0, 0.0, 1.0),
            Primitive::new(0.125, 0.0, 0.1),
        )
    }

    fn assert_state(got: Primitive, rho: f64, u: f64, p: f64, eps: f64) {
        assert_relative_eq!(got.rho, rho, epsilon = eps);
        assert_relative_eq!(got.u, u, epsilon = eps);
        assert_relative_eq!(got.p, p, epsilon = eps);
    }

    /// The standard Sod shock tube (Toro Table 4.1, test 1): the star state is a known
    /// reference, `p* = 0.30313`, `u* = 0.92745`, for γ = 1.4.
    #[test]
    fn sod_star_state_matches_toro() {
        let (left, right) = sod();
        let sol = solve(left, right, GAMMA);
        assert_relative_eq!(sol.p_star(), 0.303_130, epsilon = 1e-5);
        assert_relative_eq!(sol.u_star(), 0.927_453, epsilon = 1e-5);
    }

    /// Far outside the wave fan, sampling must return the untouched initial states.
    #[test]
    fn sod_sample_far_field_is_initial_state() {
        let (left, right) = sod();
        let sol = solve(left, right, GAMMA);
        assert_state(sol.sample(-5.0), left.rho, left.u, left.p, 1e-12);
        assert_state(sol.sample(5.0), right.rho, right.u, right.p, 1e-12);
    }

    /// Between the left rarefaction tail and the contact, the star-left density is the known
    /// `ρ*_L = 0.42632` at `u*`, `p*`. (`ξ = 0` lies in this region for Sod.)
    #[test]
    fn sod_sample_star_left_density() {
        let (left, right) = sod();
        let sol = solve(left, right, GAMMA);
        assert_state(sol.sample(0.0), 0.426_320, 0.927_453, 0.303_130, 1e-4);
    }

    /// Between the contact and the right shock, the star-right density is the known
    /// `ρ*_R = 0.26557` at `u*`, `p*`. (`ξ = 1.0` lies in this region for Sod.)
    #[test]
    fn sod_sample_star_right_density() {
        let (left, right) = sod();
        let sol = solve(left, right, GAMMA);
        assert_state(sol.sample(1.0), 0.265_574, 0.927_453, 0.303_130, 1e-4);
    }
}
