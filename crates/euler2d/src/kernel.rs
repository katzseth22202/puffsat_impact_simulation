//! The 2D axisymmetric Euler kernel: a finite-volume HLLC Godunov solver on an `(z, r)` mesh
//! (ADR-0023).
//!
//! **Scheme (D1):** second-order **MUSCL-Hancock** (van-Leer-limited slopes, a Hancock half-step
//! predictor) with the **HLLC** interface flux, advanced by **Strang dimensional splitting**
//! `Z(dt/2) · R(dt) · Z(dt/2)`. With [`Grid2D::set_axisymmetric`] the radial sweep uses the
//! conservative **cylindrical** update — radius-weighted faces `(1/r)∂(rF_r)/∂r` plus the pressure
//! source `p/r` in the radial-momentum equation — so mass, axial momentum, and energy conserve
//! exactly and only radial momentum carries a geometric source. The axis `r=0` needs no special
//! case: the inner face area `r_{1/2}=0` kills the flux through it.
//!
//! Layout: conserved cells in row-major order `idx(iz, ir) = iz·nr + ir`, `z` increasing with `iz`
//! (the plate sits at `z = 0`, `iz = 0`), `r` increasing with `ir` (the axis at `r = 0`, `ir = 0`).
//! The transverse velocity is passively advected (see [`riemann`]).

use crate::plate::PlateProfile;
use crate::riemann::{DirCons, DirFlux, DirState, hllc_flux, phys_flux};
use crate::state::{Cons, Prim};

/// A domain-edge boundary condition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Bc {
    /// Reflecting (rigid wall / symmetry axis): the sweep-normal velocity flips in the ghost cells.
    Reflect,
    /// Transmissive (zero-gradient outflow): the ghost cells copy the boundary cell.
    Transmissive,
    /// Periodic: the ghost cells wrap to the opposite end (used by the smooth order-of-accuracy
    /// test).
    Periodic,
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
    /// When true the radial sweep uses the cylindrical (axisymmetric) update; otherwise Cartesian.
    axisymmetric: bool,
    /// Finite plate at `z = 0`: when `Some(r_plate)` the `z`-lo boundary reflects for `r ≤ r_plate`
    /// and is transmissive beyond (gas past the plate edge escapes — §7). `None` ⇒ use `bc_zlo`
    /// uniformly.
    plate_radius: Option<f64>,
    /// Shallow-concave (or inclined) plate imposed as a ghost-cell immersed boundary (true-normal
    /// mirror, ADR-0023 amendment). When `Some`, solid cells under the surface are refilled before
    /// each sweep and the wall impulse is taken at the surface cell. `None` ⇒ the grid-aligned path.
    plate_profile: Option<PlateProfile>,
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
    /// state; call [`Self::init`] to set the initial condition. Boundaries default to transmissive,
    /// geometry to Cartesian.
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
            axisymmetric: false,
            plate_radius: None,
            plate_profile: None,
            u: vec![placeholder; nz * nr],
            bc_zlo: Bc::Transmissive,
            bc_zhi: Bc::Transmissive,
            bc_rlo: Bc::Transmissive,
            bc_rhi: Bc::Transmissive,
        }
    }

    /// Enable the cylindrical (axisymmetric) radial update. `r = 0` (the axis, `ir = 0`'s inner
    /// face) should carry [`Bc::Reflect`]; the zero face area there makes it a no-op regardless.
    pub fn set_axisymmetric(&mut self, on: bool) {
        self.axisymmetric = on;
    }

    /// Set a finite plate radius at `z = 0`: the boundary reflects for `r ≤ r_plate` and is
    /// transmissive beyond (the flat-plate bounce; gas past the edge escapes). `None` restores the
    /// uniform [`Self::bc_zlo`].
    pub fn set_plate_radius(&mut self, r_plate: Option<f64>) {
        self.plate_radius = r_plate;
    }

    /// Impose a curved/inclined reflecting plate as a ghost-cell immersed boundary (ADR-0023
    /// amendment): solid cells under the surface are refilled by mirroring the adjacent fluid across
    /// the true local normal before each sweep, and the axial wall impulse is taken at the surface
    /// cell. `None` restores the grid-aligned path. The `z`-lo boundary should be transmissive (the
    /// immersed surface, not the grid edge, is the wall).
    pub fn set_plate_profile(&mut self, profile: Option<PlateProfile>) {
        self.plate_profile = profile;
    }

    /// Radial coordinate of cell `ir`'s center.
    #[inline]
    fn r_center(&self, ir: usize) -> f64 {
        (ir as f64 + 0.5) * self.dr
    }

    /// Axial coordinate of cell `iz`'s center.
    #[inline]
    fn z_center(&self, iz: usize) -> f64 {
        (iz as f64 + 0.5) * self.dz
    }

    /// Whether cell `(iz, ir)` lies inside an immersed plate (always false without a profile).
    #[must_use]
    pub fn is_solid(&self, iz: usize, ir: usize) -> bool {
        self.plate_profile
            .is_some_and(|p| p.is_solid(self.z_center(iz), self.r_center(ir)))
    }

    /// Net **axial** force on the plate `= ∫ p (n̂·ẑ) dA`. For a surface `z = z_s(r)` the tilt
    /// cancels — `(n̂·ẑ)·dA = (1/√(1+s²))·(2π r √(1+s²) dr) = 2π r dr`, the projected annulus — so the
    /// axial force is `Σ p(surface, r)·r·dr` whatever the slope, with the common `2π` dropped (it
    /// cancels in the `J_wall/p_in` restitution ratio `eta_capture` uses). With an immersed
    /// [`PlateProfile`] the pressure is taken at the lowest fluid cell in each column (the surface
    /// cell); otherwise at the grid-aligned plate (`iz = 0`, `r ≤ r_plate`).
    #[must_use]
    pub fn plate_force(&self) -> f64 {
        if let Some(profile) = self.plate_profile {
            (0..self.nr)
                .filter(|&ir| profile.covers(self.r_center(ir)))
                .filter_map(|ir| {
                    self.surface_cell(ir)
                        .map(|iz| self.prim(iz, ir).p * self.r_center(ir) * self.dr)
                })
                .sum()
        } else {
            let rp = self.plate_radius.unwrap_or(f64::INFINITY);
            (0..self.nr)
                .filter(|&ir| self.r_center(ir) <= rp)
                .map(|ir| self.prim(0, ir).p * self.r_center(ir) * self.dr)
                .sum()
        }
    }

    /// The lowest fluid cell in column `ir` — the cell whose pressure acts on the immersed surface —
    /// or `None` if the entire column is solid.
    fn surface_cell(&self, ir: usize) -> Option<usize> {
        (0..self.nz).find(|&iz| !self.is_solid(iz, ir))
    }

    /// Total axial momentum `Σ ρu_z dV` (with the common `2π` dropped, as in [`Self::plate_force`]).
    /// Its magnitude at `t = 0` is the incident axial momentum `p_in`.
    #[must_use]
    pub fn axial_momentum(&self) -> f64 {
        (0..self.nz)
            .flat_map(|iz| (0..self.nr).map(move |ir| (iz, ir)))
            .map(|(iz, ir)| self.cons(iz, ir).mz * self.r_center(ir) * self.dr * self.dz)
            .sum()
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
        for iz in 0..self.nz {
            for ir in 0..self.nr {
                // Solid (under-plate) cells hold mirror/garbage values between sweeps; they do not
                // constrain the physical time step.
                if self.is_solid(iz, ir) {
                    continue;
                }
                let w = self.prim(iz, ir);
                let cs = w.sound_speed(self.gamma);
                inv = inv
                    .max((w.uz.abs() + cs) / self.dz)
                    .max((w.ur.abs() + cs) / self.dr);
            }
        }
        self.cfl / inv
    }

    /// Advance one step of size `dt` by Strang splitting `Z(dt/2) · R(dt) · Z(dt/2)` (second order
    /// in time when each sweep is the second-order MUSCL-Hancock solve). With an immersed
    /// [`PlateProfile`] the reflecting ghost values are re-imposed before each sub-sweep; without
    /// one [`Self::apply_immersed_bc`] is a no-op, so the verification suite is untouched.
    pub fn step(&mut self, dt: f64) {
        self.apply_immersed_bc();
        self.sweep(Axis::Z, 0.5 * dt);
        self.apply_immersed_bc();
        self.sweep(Axis::R, dt);
        self.apply_immersed_bc();
        self.sweep(Axis::Z, 0.5 * dt);
    }

    /// Refill solid cells under an immersed plate surface by mirroring the adjacent fluid across the
    /// true local normal — the ghost-cell immersed boundary (ADR-0023 amendment). For each solid
    /// cell: reflect its center across the surface to an image point in the fluid, sample the fluid
    /// state there (nearest fluid cell), copy `ρ` and `p`, and reverse the wall-normal velocity
    /// component (`u → u − 2(u·n̂)n̂`). No-op when no profile is set.
    // n_z/n_r are the two components of one normal vector — intrinsically similar names.
    #[allow(clippy::similar_names)]
    fn apply_immersed_bc(&mut self) {
        let Some(profile) = self.plate_profile else {
            return;
        };
        for iz in 0..self.nz {
            for ir in 0..self.nr {
                let z = self.z_center(iz);
                let r = self.r_center(ir);
                if !profile.is_solid(z, r) {
                    continue;
                }
                let dist = profile.signed_distance(z, r); // < 0 inside the solid
                let (n_z, n_r) = profile.normal(r);
                // Image point = the cell center reflected across the surface (into the fluid).
                let z_img = z - 2.0 * dist * n_z;
                let r_img = (r - 2.0 * dist * n_r).max(0.0);
                let src = self.sample_fluid(z_img, r_img, profile);
                // Mirror: reverse the wall-normal velocity component, leaving ρ and p.
                let u_n = src.uz * n_z + src.ur * n_r;
                let mirrored = Prim::new(
                    src.rho,
                    src.uz - 2.0 * u_n * n_z,
                    src.ur - 2.0 * u_n * n_r,
                    src.p,
                );
                let k = self.idx(iz, ir);
                self.u[k] = Cons::from_prim(mirrored, self.gamma);
            }
        }
    }

    /// Sample the fluid primitive nearest the image point `(z, r)` for the immersed-boundary mirror.
    /// The image lies in the fluid by construction; should rounding land on a solid cell (a near-rim
    /// corner) march up in `z` to the first fluid cell.
    fn sample_fluid(&self, z: f64, r: f64, profile: PlateProfile) -> Prim {
        let ir = nearest_index(r, self.dr, self.nr);
        let mut iz = nearest_index(z, self.dz, self.nz);
        while iz + 1 < self.nz && profile.is_solid(self.z_center(iz), self.r_center(ir)) {
            iz += 1;
        }
        self.prim(iz, ir)
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

    /// One MUSCL-Hancock sweep over every line along `axis`, updating the conserved state in place.
    // The four momentum-flux/state components (mn/mt and the mapped mz/mr) are intrinsically
    // similarly named — they are one momentum vector — so the lint is silenced for this routine.
    #[allow(clippy::similar_names)]
    fn sweep(&mut self, axis: Axis, dt: f64) {
        let (n_lines, n, dx, bc_lo, bc_hi) = match axis {
            Axis::Z => (self.nr, self.nz, self.dz, self.bc_zlo, self.bc_zhi),
            Axis::R => (self.nz, self.nr, self.dr, self.bc_rlo, self.bc_rhi),
        };
        let axisym_r = self.axisymmetric && axis == Axis::R;
        let gamma = self.gamma;

        for line in 0..n_lines {
            let cell = |i: usize| match axis {
                Axis::Z => (i, line),
                Axis::R => (line, i),
            };

            // Directional conserved state, padded with two ghost layers each side: physical cell i
            // sits at padded index i + 2.
            let mut p = vec![DirCons::default(); n + 4];
            for i in 0..n {
                let (iz, ir) = cell(i);
                p[i + 2] = dir_cons(self.cons(iz, ir), axis);
            }
            // The finite plate makes the z-lo boundary depend on the line's radius: reflect on the
            // plate (r ≤ r_plate), transmissive past its edge. Other boundaries are uniform.
            let bc_lo_line = match (axis, self.plate_radius) {
                (Axis::Z, Some(rp)) if self.r_center(line) > rp => Bc::Transmissive,
                (Axis::Z, Some(_)) => Bc::Reflect,
                _ => bc_lo,
            };
            fill_ghosts(&mut p, n, bc_lo_line, bc_hi);

            // MUSCL slope reconstruction + Hancock half-step predictor for every padded index that
            // borders a physical face (1 ..= n+2). The predictor uses the plane flux (the geometric
            // weighting enters only the corrector below).
            let half = 0.5 * dt / dx;
            let mut face_l = vec![DirState::default(); n + 4];
            let mut face_r = vec![DirState::default(); n + 4];
            for j in 1..=n + 2 {
                let back = p[j].axpy(-1.0, p[j - 1]);
                let fwd = p[j + 1].axpy(-1.0, p[j]);
                let slope = limited_slope(back, fwd, p[j], gamma);
                let ql = p[j].axpy(-0.5, slope);
                let qr = p[j].axpy(0.5, slope);
                let fl = phys_flux(ql.to_state(gamma), gamma);
                let fr = phys_flux(qr.to_state(gamma), gamma);
                // q* = q ± ½ Δ + ½(dt/dx)(F(q_L) − F(q_R)).
                face_l[j] = ql.add_flux(half, fl).add_flux(-half, fr).to_state(gamma);
                face_r[j] = qr.add_flux(half, fl).add_flux(-half, fr).to_state(gamma);
            }

            // Face fluxes: face j (1 ..= n+1) sits between padded cells j and j+1.
            let mut flux = vec![DirFlux::default(); n + 2];
            for j in 1..=n + 1 {
                flux[j] = hllc_flux(face_r[j], face_l[j + 1], gamma);
            }

            // Conservative corrector. Cartesian cells use F_hi − F_lo; the cylindrical radial sweep
            // weights faces by their radius and adds the pressure source p/r to radial momentum.
            for i in 0..n {
                let (iz, ir) = cell(i);
                let k = self.idx(iz, ir);
                let f_lo = flux[i + 1];
                let f_hi = flux[i + 2];
                let dc = dir_cons(self.u[k], axis);
                let updated = if axisym_r {
                    let r_lo = i as f64 * dx;
                    let r_hi = (i + 1) as f64 * dx;
                    let r_c = (i as f64 + 0.5) * dx;
                    let inv_vol = dt / (r_c * dx);
                    let p_src = dt * Prim::from_cons(self.u[k], gamma).p / r_c;
                    DirCons {
                        rho: dc.rho - inv_vol * (r_hi * f_hi.rho - r_lo * f_lo.rho),
                        mn: dc.mn - inv_vol * (r_hi * f_hi.mn - r_lo * f_lo.mn) + p_src,
                        mt: dc.mt - inv_vol * (r_hi * f_hi.mt - r_lo * f_lo.mt),
                        e: dc.e - inv_vol * (r_hi * f_hi.e - r_lo * f_lo.e),
                    }
                } else {
                    let inv = dt / dx;
                    DirCons {
                        rho: dc.rho - inv * (f_hi.rho - f_lo.rho),
                        mn: dc.mn - inv * (f_hi.mn - f_lo.mn),
                        mt: dc.mt - inv * (f_hi.mt - f_lo.mt),
                        e: dc.e - inv * (f_hi.e - f_lo.e),
                    }
                };
                self.u[k] = floored(from_dir_cons(updated, axis), gamma);
            }
        }
    }
}

/// Density / pressure floor below which a cell is treated as vacuum. The bounce's rarefaction tail
/// evacuates to near-vacuum; such cells carry negligible mass and momentum, so resetting them to a
/// quiescent floor keeps the scheme robust without affecting the wall impulse (the same spirit as
/// the 1D tail guard, ADR-0001).
const RHO_FLOOR: f64 = 1.0e-9;
const P_FLOOR: f64 = 1.0e-12;

/// Pass a healthy cell through unchanged; reset a vacuum/invalid cell (non-positive or non-finite
/// density or pressure) to a quiescent floor state.
fn floored(c: Cons, gamma: f64) -> Cons {
    let w = Prim::from_cons(c, gamma);
    if w.rho.is_finite() && w.p.is_finite() && w.rho >= RHO_FLOOR && w.p >= P_FLOOR {
        c
    } else {
        Cons::from_prim(Prim::new(RHO_FLOOR, 0.0, 0.0, P_FLOOR), gamma)
    }
}

/// Index of the cell whose center is nearest coordinate `x` along an axis of spacing `dx` with `n`
/// cells, clamped to `[0, n−1]`.
fn nearest_index(x: f64, dx: f64, n: usize) -> usize {
    let f = (x / dx - 0.5).round().max(0.0);
    // SAFE: `f` is clamped to ≥ 0 above and to ≤ n−1 below; no truncation or sign loss.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let i = f as usize;
    i.min(n - 1)
}

/// Map a full conserved cell to the directional conserved vector for a sweep axis (reorder the
/// momentum so the sweep-normal component is `mn`).
fn dir_cons(c: Cons, axis: Axis) -> DirCons {
    let (mn, mt) = match axis {
        Axis::Z => (c.mz, c.mr),
        Axis::R => (c.mr, c.mz),
    };
    DirCons {
        rho: c.rho,
        mn,
        mt,
        e: c.e_tot,
    }
}

/// Inverse of [`dir_cons`]: map a directional conserved vector back to a full conserved cell.
fn from_dir_cons(dc: DirCons, axis: Axis) -> Cons {
    let (mz, mr) = match axis {
        Axis::Z => (dc.mn, dc.mt),
        Axis::R => (dc.mt, dc.mn),
    };
    Cons {
        rho: dc.rho,
        mz,
        mr,
        e_tot: dc.e,
    }
}

/// Reflection of a directional conserved cell: the sweep-normal momentum flips (rigid wall / axis).
fn reflect(dc: DirCons) -> DirCons {
    DirCons { mn: -dc.mn, ..dc }
}

/// Fill the two ghost layers each side of the padded line `p` (length `n + 4`, physical cells at
/// `2 ..= n+1`) per the low/high boundary conditions.
fn fill_ghosts(p: &mut [DirCons], n: usize, bc_lo: Bc, bc_hi: Bc) {
    match bc_lo {
        Bc::Reflect => {
            p[1] = reflect(p[2]);
            p[0] = reflect(p[3]);
        }
        Bc::Transmissive => {
            p[1] = p[2];
            p[0] = p[2];
        }
        Bc::Periodic => {
            p[1] = p[n + 1];
            p[0] = p[n];
        }
    }
    match bc_hi {
        Bc::Reflect => {
            p[n + 2] = reflect(p[n + 1]);
            p[n + 3] = reflect(p[n]);
        }
        Bc::Transmissive => {
            p[n + 2] = p[n + 1];
            p[n + 3] = p[n + 1];
        }
        Bc::Periodic => {
            p[n + 2] = p[2];
            p[n + 3] = p[3];
        }
    }
}

/// The van Leer limited slope, per conserved component (zero where the differences disagree in
/// sign — a local extremum). Falls back to a zero slope (first order) if the extrapolated face
/// states would lose positivity of density or pressure.
fn limited_slope(back: DirCons, fwd: DirCons, center: DirCons, gamma: f64) -> DirCons {
    let slope = DirCons {
        rho: vanleer(back.rho, fwd.rho),
        mn: vanleer(back.mn, fwd.mn),
        mt: vanleer(back.mt, fwd.mt),
        e: vanleer(back.e, fwd.e),
    };
    let lo = center.axpy(-0.5, slope).to_state(gamma);
    let hi = center.axpy(0.5, slope).to_state(gamma);
    if lo.rho <= 0.0 || lo.p <= 0.0 || hi.rho <= 0.0 || hi.p <= 0.0 {
        DirCons::default()
    } else {
        slope
    }
}

/// Van Leer slope limiter for one component: `2ab/(a+b)` when `a, b` share a sign, else `0`.
fn vanleer(a: f64, b: f64) -> f64 {
    if a * b <= 0.0 {
        0.0
    } else {
        2.0 * a * b / (a + b)
    }
}

impl Default for DirCons {
    fn default() -> Self {
        Self {
            rho: 0.0,
            mn: 0.0,
            mt: 0.0,
            e: 0.0,
        }
    }
}

impl Default for DirState {
    fn default() -> Self {
        Self {
            rho: 0.0,
            un: 0.0,
            ut: 0.0,
            p: 0.0,
        }
    }
}

impl Default for DirFlux {
    fn default() -> Self {
        Self {
            rho: 0.0,
            mn: 0.0,
            mt: 0.0,
            e: 0.0,
        }
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

    /// Total mass is conserved across many steps with reflecting walls all round (a closed box) to
    /// round-off — the conservative-form bookkeeping check.
    #[test]
    fn closed_box_conserves_mass() {
        let mut g = Grid2D::new(10, 8, 0.1, 0.1, GAMMA);
        g.bc_zlo = Bc::Reflect;
        g.bc_zhi = Bc::Reflect;
        g.bc_rlo = Bc::Reflect;
        g.bc_rhi = Bc::Reflect;
        g.init(|iz, ir| {
            let hot = (4..6).contains(&iz) && ir < 3;
            Prim::new(1.0, 0.0, 0.0, if hot { 5.0 } else { 1.0 })
        });
        let mass = |g: &Grid2D| -> f64 {
            (0..g.nz())
                .flat_map(|iz| (0..g.nr()).map(move |ir| (iz, ir)))
                .map(|(iz, ir)| g.cons(iz, ir).rho)
                .sum()
        };
        let mass0 = mass(&g);
        for _ in 0..20 {
            let dt = g.stable_dt();
            g.step(dt);
        }
        assert!((mass(&g) - mass0).abs() < 1e-10, "mass drift");
    }
}
