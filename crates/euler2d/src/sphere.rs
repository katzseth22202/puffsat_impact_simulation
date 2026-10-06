//! A 1D spherically symmetric Euler solver for a blast inside a closed spherical vessel.
//!
//! The walled-nozzle chamber is a sphere with the energy released near its centre. Its wall load
//! is the blast's reflected spike, then reverberations that settle to the quasi-static pressure.
//! That is a 1D spherical problem. The 2D kernel's immersed boundary takes only single-valued
//! surfaces `z = z_s(r)`, and in 1D the wall spike resolves far more finely.
//!
//! **Scheme.** The same building blocks as the 2D kernel (ADR-0023): MUSCL-Hancock with
//! minmod-limited slopes on the conserved variables, and the HLLC flux of [`crate::riemann`]. The
//! spherical update is conservative in mass and energy:
//!
//! `V_i dU_i/dt = A_{i-1/2} F_{i-1/2} - A_{i+1/2} F_{i+1/2} + S_i`,
//!
//! with `A = r^2`, `V = (r_+^3 - r_-^3)/3` (the common `4 pi` dropped) and the momentum source
//! `S = p_i (A_{i+1/2} - A_{i-1/2})`. That source balances a uniform pressure exactly, so gas at
//! rest stays at rest. The centre face has zero area. The outer face is a rigid wall, and the
//! momentum flux through it is the wall pressure.
//!
//! The Hancock predictor is planar (it drops the geometric terms for the half step). That is
//! second order away from the centre, and the wall is where the answer is read.

use crate::riemann::{DirCons, DirFlux, DirState, hllc_flux, phys_flux};

/// Density floor [kg/m^3] and pressure floor [Pa] that keep the post-blast cavity positive.
const RHO_FLOOR: f64 = 1e-8;
const P_FLOOR: f64 = 1e-3;

/// The closed sphere and its gas.
#[derive(Debug, Clone)]
pub struct Sphere1D {
    n: usize,
    dr: f64,
    gamma: f64,
    cfl: f64,
    u: Vec<DirCons>,
    wall_pressure: f64,
}

impl Sphere1D {
    /// A sphere of radius `radius` in `n` equal shells, γ-law gas, filled from `f(r) -> (ρ, u, p)`.
    ///
    /// # Panics
    /// If `n` is zero or `radius` is not positive.
    #[must_use]
    pub fn new(n: usize, radius: f64, gamma: f64, f: impl Fn(f64) -> (f64, f64, f64)) -> Self {
        assert!(n > 0 && radius > 0.0, "sphere must be non-empty");
        let dr = radius / n as f64;
        let u = (0..n)
            .map(|i| {
                let (rho, un, p) = f((i as f64 + 0.5) * dr);
                DirCons::from_state(
                    DirState {
                        rho,
                        un,
                        ut: 0.0,
                        p,
                    },
                    gamma,
                )
            })
            .collect();
        Self {
            n,
            dr,
            gamma,
            cfl: 0.4,
            u,
            wall_pressure: 0.0,
        }
    }

    /// Shell volume (without `4 pi`) of cell `i`.
    fn volume(&self, i: usize) -> f64 {
        let (lo, hi) = (i as f64 * self.dr, (i as f64 + 1.0) * self.dr);
        (hi.powi(3) - lo.powi(3)) / 3.0
    }

    /// Face area (without `4 pi`) of face `j` (face `j` is the inner face of cell `j`).
    fn area(&self, j: usize) -> f64 {
        (j as f64 * self.dr).powi(2)
    }

    /// Total mass [kg] (with `4 pi`).
    #[must_use]
    pub fn mass(&self) -> f64 {
        4.0 * std::f64::consts::PI
            * (0..self.n)
                .map(|i| self.u[i].rho * self.volume(i))
                .sum::<f64>()
    }

    /// Total energy [J] (with `4 pi`).
    #[must_use]
    pub fn energy(&self) -> f64 {
        4.0 * std::f64::consts::PI
            * (0..self.n)
                .map(|i| self.u[i].e * self.volume(i))
                .sum::<f64>()
    }

    /// Primitive state of cell `i` as `(ρ, u, p)`.
    #[must_use]
    pub fn prim(&self, i: usize) -> (f64, f64, f64) {
        let s = self.u[i].to_state(self.gamma);
        (s.rho, s.un, s.p)
    }

    /// Pressure on the wall from the last step [Pa].
    #[must_use]
    pub fn wall_pressure(&self) -> f64 {
        self.wall_pressure
    }

    /// Add `energy` [J] of internal energy and `mass` [kg] spread evenly by volume over `r < radius`.
    /// This is a charge going off at the centre (a rod striking its plug).
    pub fn deposit(&mut self, energy: f64, mass: f64, radius: f64) {
        let cells: Vec<usize> = (0..self.n)
            .filter(|&i| (i as f64 + 0.5) * self.dr < radius)
            .collect();
        let total: f64 =
            cells.iter().map(|&i| self.volume(i)).sum::<f64>() * 4.0 * std::f64::consts::PI;
        for i in cells {
            self.u[i].e += energy / total;
            // Added mass arrives at rest in this frame: no momentum, no kinetic energy.
            self.u[i].rho += mass / total;
        }
    }

    /// CFL-limited time step [s].
    #[must_use]
    pub fn stable_dt(&self) -> f64 {
        let fastest = self
            .u
            .iter()
            .map(|c| {
                let s = c.to_state(self.gamma);
                s.un.abs() + (self.gamma * s.p.max(P_FLOOR) / s.rho.max(RHO_FLOOR)).sqrt()
            })
            .fold(0.0_f64, f64::max);
        self.cfl * self.dr / fastest
    }

    fn mirror(c: DirCons) -> DirCons {
        DirCons { mn: -c.mn, ..c }
    }

    fn limited(a: f64, b: f64) -> f64 {
        if a * b <= 0.0 {
            0.0
        } else if a.abs() < b.abs() {
            a
        } else {
            b
        }
    }

    fn slope(l: DirCons, c: DirCons, r: DirCons) -> DirCons {
        DirCons {
            rho: Self::limited(c.rho - l.rho, r.rho - c.rho),
            mn: Self::limited(c.mn - l.mn, r.mn - c.mn),
            mt: 0.0,
            e: Self::limited(c.e - l.e, r.e - c.e),
        }
    }

    fn state(&self, c: DirCons) -> DirState {
        let mut s = c.to_state(self.gamma);
        s.rho = s.rho.max(RHO_FLOOR);
        s.p = s.p.max(P_FLOOR);
        s.ut = 0.0;
        s
    }

    /// Advance one step of `dt` [s].
    pub fn step(&mut self, dt: f64) {
        let n = self.n;
        let g = self.gamma;
        let half = 0.5 * dt / self.dr;
        // Reconstruct and predict each cell's left and right face states.
        let mut left = Vec::with_capacity(n);
        let mut right = Vec::with_capacity(n);
        for i in 0..n {
            let c = self.u[i];
            let l = if i == 0 {
                Self::mirror(c)
            } else {
                self.u[i - 1]
            };
            let r = if i + 1 == n {
                Self::mirror(c)
            } else {
                self.u[i + 1]
            };
            let d = Self::slope(l, c, r);
            let ul = c.axpy(-0.5, d);
            let ur = c.axpy(0.5, d);
            let fl = phys_flux(self.state(ul), g);
            let fr = phys_flux(self.state(ur), g);
            let diff = DirFlux {
                rho: fl.rho - fr.rho,
                mn: fl.mn - fr.mn,
                mt: 0.0,
                e: fl.e - fr.e,
            };
            left.push(ul.add_flux(half, diff));
            right.push(ur.add_flux(half, diff));
        }
        // Face fluxes: face j sits between cell j-1 (its right state) and cell j (its left).
        let mut flux = Vec::with_capacity(n + 1);
        flux.push(DirFlux {
            rho: 0.0,
            mn: 0.0,
            mt: 0.0,
            e: 0.0,
        });
        for j in 1..n {
            flux.push(hllc_flux(self.state(right[j - 1]), self.state(left[j]), g));
        }
        let wall_in = self.state(right[n - 1]);
        let wall_out = DirState {
            un: -wall_in.un,
            ..wall_in
        };
        let wall = hllc_flux(wall_in, wall_out, g);
        self.wall_pressure = wall.mn;
        flux.push(DirFlux {
            rho: 0.0,
            mn: wall.mn,
            mt: 0.0,
            e: 0.0,
        });
        for i in 0..n {
            let (a_lo, a_hi) = (self.area(i), self.area(i + 1));
            let v = self.volume(i);
            let p_mid = 0.5 * (self.state(left[i]).p + self.state(right[i]).p);
            let c = &mut self.u[i];
            c.rho += dt / v * (a_lo * flux[i].rho - a_hi * flux[i + 1].rho);
            c.mn += dt / v * (a_lo * flux[i].mn - a_hi * flux[i + 1].mn + p_mid * (a_hi - a_lo));
            c.e += dt / v * (a_lo * flux[i].e - a_hi * flux[i + 1].e);
            c.rho = c.rho.max(RHO_FLOOR);
            let kinetic = 0.5 * c.mn * c.mn / c.rho;
            c.e = c.e.max(kinetic + P_FLOOR / (g - 1.0));
        }
    }

    /// Run to `t_end` [s], calling `record(t, wall_pressure)` after every step.
    pub fn run_to(&mut self, t0: f64, t_end: f64, mut record: impl FnMut(f64, f64)) -> f64 {
        let mut t = t0;
        while t < t_end {
            let dt = self.stable_dt().min(t_end - t);
            self.step(dt);
            t += dt;
            record(t, self.wall_pressure);
        }
        t
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn uniform(n: usize) -> Sphere1D {
        Sphere1D::new(n, 1.0, 1.4, |_| (1.0, 0.0, 1.0))
    }

    #[test]
    fn gas_at_rest_stays_at_rest() {
        let mut s = uniform(100);
        s.run_to(0.0, 0.5, |_, _| {});
        for i in 0..100 {
            let (rho, u, p) = s.prim(i);
            assert!((rho - 1.0).abs() < 1e-12 && u.abs() < 1e-12 && (p - 1.0).abs() < 1e-12);
        }
        assert!((s.wall_pressure() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn closed_sphere_conserves_mass_and_energy() {
        let mut s = uniform(200);
        s.deposit(50.0, 0.0, 0.1);
        let (m0, e0) = (s.mass(), s.energy());
        s.run_to(0.0, 0.3, |_, _| {});
        assert!((s.mass() - m0).abs() / m0 < 1e-10);
        assert!((s.energy() - e0).abs() / e0 < 1e-6);
    }

    /// Sedov-Taylor: `R = xi0 (E t^2 / rho)^(1/5)` with xi0 = 1.0328 at gamma 1.4, spherical
    /// (`alpha` = 0.851; Kamm and Timmes, LA-UR-07-2849).
    #[test]
    fn blast_radius_follows_sedov() {
        let mut s = Sphere1D::new(800, 1.0, 1.4, |_| (1.0, 0.0, 1e-6));
        let energy = 1.0;
        s.deposit(energy, 0.0, 0.01);
        let t_end = 0.15;
        s.run_to(0.0, t_end, |_, _| {});
        // Shock radius: the outermost cell well above ambient density.
        let shock = (0..800)
            .rev()
            .find(|&i| s.prim(i).0 > 2.0)
            .map_or(0.0, |i| (i as f64 + 0.5) / 800.0);
        let expected = 1.0328 * (energy * t_end * t_end).powf(0.2);
        assert!(
            (shock - expected).abs() / expected < 0.01,
            "shock {shock} vs {expected}"
        );
    }
}
