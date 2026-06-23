//! The 2D axisymmetric Euler kernel: a finite-volume HLLC Godunov solver on an `(z, r)` mesh
//! (ADR-0023).
//!
//! **D0 (this slice) is Cartesian** — first-order Godunov with dimensional (Godunov) splitting,
//! verified on the Sod shock tube embedded in 2D. The cylindrical geometric source
//! `(1/r)∂(rF_r)/∂r` + `p/r` and the second-order MUSCL-Hancock reconstruction land in D1; the
//! flat-plate bounce and `eta_capture` extraction in D2.
//!
//! Layout: conserved cells in row-major order `idx(iz, ir) = iz·nr + ir`, `z` increasing with `iz`
//! (the plate sits at `z = 0`, `iz = 0`), `r` increasing with `ir` (the axis at `r = 0`, `ir = 0`).
//! Each step sweeps `z` then `r`; the transverse velocity is passively advected (see [`riemann`]).

use crate::riemann::{DirFlux, DirState, hllc_flux};
use crate::state::{Cons, Prim};

/// A domain-edge boundary condition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Bc {
    /// Reflecting (rigid wall / symmetry axis): the sweep-normal velocity flips in the ghost cell.
    Reflect,
    /// Transmissive (zero-gradient outflow): the ghost cell copies the boundary cell.
    Transmissive,
}

/// Which axis a 1D sweep runs along — fixes which momentum component is "normal".
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Axis {
    Z,
    R,
}

/// The 2D Euler mesh and its state.
#[derive(Debug, Clone)]
pub struct Grid2D {
    nz: usize,
    nr: usize,
    dz: f64,
    dr: f64,
    gamma: f64,
    cfl: f64,
    /// Conserved cells, `idx(iz, ir) = iz·nr + ir`.
    u: Vec<Cons>,
    /// `z = 0` (plate) and `z = z_max` boundaries.
    pub bc_zlo: Bc,
    pub bc_zhi: Bc,
    /// `r = 0` (axis) and `r = r_max` boundaries.
    pub bc_rlo: Bc,
    pub bc_rhi: Bc,
}

impl Grid2D {
    /// A mesh of `nz × nr` cells with the given spacings and γ, all cells set to a unit placeholder
    /// state; call [`Self::init`] to set the initial condition. Boundaries default to transmissive.
    ///
    /// # Panics
    /// If `nz` or `nr` is zero (an empty mesh).
    #[must_use]
    pub fn new(nz: usize, nr: usize, dz: f64, dr: f64, gamma: f64) -> Self {
        assert!(nz > 0 && nr > 0, "grid must be non-empty");
        let placeholder = Cons::from_prim(Prim::new(1.0, 0.0, 0.0, 1.0), gamma);
        Self {
            nz,
            nr,
            dz,
            dr,
            gamma,
            cfl: 0.4,
            u: vec![placeholder; nz * nr],
            bc_zlo: Bc::Transmissive,
            bc_zhi: Bc::Transmissive,
            bc_rlo: Bc::Transmissive,
            bc_rhi: Bc::Transmissive,
        }
    }

    #[inline]
    fn idx(&self, iz: usize, ir: usize) -> usize {
        iz * self.nr + ir
    }

    /// Number of axial / radial cells.
    #[must_use]
    pub fn nz(&self) -> usize {
        self.nz
    }
    #[must_use]
    pub fn nr(&self) -> usize {
        self.nr
    }

    /// Set the initial condition cell-by-cell from a primitive-state function of `(iz, ir)`.
    pub fn init(&mut self, f: impl Fn(usize, usize) -> Prim) {
        for iz in 0..self.nz {
            for ir in 0..self.nr {
                let k = self.idx(iz, ir);
                self.u[k] = Cons::from_prim(f(iz, ir), self.gamma);
            }
        }
    }

    /// Primitive state of cell `(iz, ir)`.
    #[must_use]
    pub fn prim(&self, iz: usize, ir: usize) -> Prim {
        Prim::from_cons(self.u[self.idx(iz, ir)], self.gamma)
    }

    /// Conserved state of cell `(iz, ir)`.
    #[must_use]
    pub fn cons(&self, iz: usize, ir: usize) -> Cons {
        self.u[self.idx(iz, ir)]
    }

    /// Stable time step `dt = cfl / max_cell max((|u_z|+c)/dz, (|u_r|+c)/dr)` (the dimensional-split
    /// CFL bound: each 1D sweep must satisfy its own CFL).
    #[must_use]
    pub fn stable_dt(&self) -> f64 {
        let mut inv = 0.0_f64;
        for &c in &self.u {
            let w = Prim::from_cons(c, self.gamma);
            let cs = w.sound_speed(self.gamma);
            inv = inv
                .max((w.uz.abs() + cs) / self.dz)
                .max((w.ur.abs() + cs) / self.dr);
        }
        self.cfl / inv
    }

    /// Advance one step of size `dt` (z-sweep then r-sweep; Godunov splitting, first order).
    pub fn step(&mut self, dt: f64) {
        self.sweep(Axis::Z, dt);
        self.sweep(Axis::R, dt);
    }

    /// Advance to `t_end`, clamping the final step to land exactly on it. Returns the step count.
    ///
    /// # Panics
    /// If the run exceeds a large step cap without reaching `t_end` (a non-converging time step,
    /// e.g. a NaN sound speed) — a guard against silent infinite loops.
    pub fn run_to(&mut self, t_end: f64) -> usize {
        let mut t = 0.0;
        let mut steps = 0;
        while t < t_end {
            let dt = self.stable_dt().min(t_end - t);
            self.step(dt);
            t += dt;
            steps += 1;
            assert!(steps < 10_000_000, "run_to did not converge");
        }
        steps
    }

    /// One directional sweep over every line in `axis`, updating the conserved state in place.
    fn sweep(&mut self, axis: Axis, dt: f64) {
        let (n_lines, n_cells, dx, bc_lo, bc_hi) = match axis {
            Axis::Z => (self.nr, self.nz, self.dz, self.bc_zlo, self.bc_zhi),
            Axis::R => (self.nz, self.nr, self.dr, self.bc_rlo, self.bc_rhi),
        };
        for line in 0..n_lines {
            // Gather the line's directional states.
            let cells: Vec<DirState> = (0..n_cells)
                .map(|i| {
                    let (iz, ir) = match axis {
                        Axis::Z => (i, line),
                        Axis::R => (line, i),
                    };
                    dir_state(self.cons(iz, ir), self.gamma, axis)
                })
                .collect();

            // Pad with ghost cells, then compute the n_cells+1 face fluxes.
            let ghost_lo = ghost(cells[0], bc_lo);
            let ghost_hi = ghost(cells[n_cells - 1], bc_hi);
            let mut fluxes: Vec<DirFlux> = Vec::with_capacity(n_cells + 1);
            for k in 0..=n_cells {
                let left = if k == 0 { ghost_lo } else { cells[k - 1] };
                let right = if k == n_cells { ghost_hi } else { cells[k] };
                fluxes.push(hllc_flux(left, right, self.gamma));
            }

            // Conservative update U_i -= (dt/dx)(F_{i+1/2} − F_{i−1/2}).
            for i in 0..n_cells {
                let (iz, ir) = match axis {
                    Axis::Z => (i, line),
                    Axis::R => (line, i),
                };
                let k = self.idx(iz, ir);
                self.u[k] = advance(self.u[k], fluxes[i], fluxes[i + 1], dt / dx, axis);
            }
        }
    }
}

/// Build a directional primitive `(ρ, u_n, u_t, p)` from a conserved cell for the given sweep axis.
fn dir_state(c: Cons, gamma: f64, axis: Axis) -> DirState {
    let w = Prim::from_cons(c, gamma);
    let (un, ut) = match axis {
        Axis::Z => (w.uz, w.ur),
        Axis::R => (w.ur, w.uz),
    };
    DirState {
        rho: w.rho,
        un,
        ut,
        p: w.p,
    }
}

/// The ghost state outside a boundary cell: a transmissive copy, or a reflection (normal velocity
/// flipped) for a rigid wall / symmetry axis.
fn ghost(boundary: DirState, bc: Bc) -> DirState {
    match bc {
        Bc::Transmissive => boundary,
        Bc::Reflect => DirState {
            un: -boundary.un,
            ..boundary
        },
    }
}

/// Apply the conservative flux difference to a cell, mapping the directional flux
/// `[ρ, ρu_n, ρu_t, E]` back onto `[ρ, ρu_z, ρu_r, E]` for the sweep axis.
// The momentum-component deltas (d_mn/d_mt/d_mz/d_mr) are deliberately near-identical: they are the
// four components of one momentum vector, so the names mirror the physics rather than disambiguate.
#[allow(clippy::similar_names)]
fn advance(c: Cons, f_lo: DirFlux, f_hi: DirFlux, inv_dx_dt: f64, axis: Axis) -> Cons {
    let d_rho = inv_dx_dt * (f_hi.rho - f_lo.rho);
    let d_mn = inv_dx_dt * (f_hi.mn - f_lo.mn);
    let d_mt = inv_dx_dt * (f_hi.mt - f_lo.mt);
    let d_e = inv_dx_dt * (f_hi.e - f_lo.e);
    let (d_mz, d_mr) = match axis {
        Axis::Z => (d_mn, d_mt),
        Axis::R => (d_mt, d_mn),
    };
    Cons {
        rho: c.rho - d_rho,
        mz: c.mz - d_mz,
        mr: c.mr - d_mr,
        e_tot: c.e_tot - d_e,
    }
}

#[cfg(test)]
mod tests {
    use super::{Bc, Grid2D};
    use crate::state::Prim;

    const GAMMA: f64 = 1.4;

    /// A uniform, motionless state is a fixed point of the scheme to round-off (no spurious fluxes).
    #[test]
    fn uniform_state_is_steady() {
        let mut g = Grid2D::new(8, 6, 0.1, 0.1, GAMMA);
        g.init(|_, _| Prim::new(1.3, 0.0, 0.0, 2.0));
        g.step(0.01);
        for iz in 0..g.nz() {
            for ir in 0..g.nr() {
                let w = g.prim(iz, ir);
                assert!((w.rho - 1.3).abs() < 1e-13);
                assert!((w.p - 2.0).abs() < 1e-13);
                assert!(w.uz.abs() < 1e-13 && w.ur.abs() < 1e-13);
            }
        }
    }

    /// Total mass and total axial momentum are conserved across a step with reflecting walls all
    /// round (a closed box) to round-off — the conservative-form bookkeeping check.
    #[test]
    fn closed_box_conserves_mass_and_momentum() {
        let mut g = Grid2D::new(10, 8, 0.1, 0.1, GAMMA);
        g.bc_zlo = Bc::Reflect;
        g.bc_zhi = Bc::Reflect;
        g.bc_rlo = Bc::Reflect;
        g.bc_rhi = Bc::Reflect;
        // A pressure blob that will drive flow but cannot leave the closed box.
        g.init(|iz, ir| {
            let hot = (4..6).contains(&iz) && ir < 3;
            Prim::new(1.0, 0.0, 0.0, if hot { 5.0 } else { 1.0 })
        });
        let mass0: f64 = (0..g.nz())
            .flat_map(|iz| (0..g.nr()).map(move |ir| (iz, ir)))
            .map(|(iz, ir)| g.cons(iz, ir).rho)
            .sum();
        for _ in 0..20 {
            let dt = g.stable_dt();
            g.step(dt);
        }
        let mass1: f64 = (0..g.nz())
            .flat_map(|iz| (0..g.nr()).map(move |ir| (iz, ir)))
            .map(|(iz, ir)| g.cons(iz, ir).rho)
            .sum();
        assert!(
            (mass1 - mass0).abs() < 1e-10,
            "mass drift {}",
            mass1 - mass0
        );
    }
}
