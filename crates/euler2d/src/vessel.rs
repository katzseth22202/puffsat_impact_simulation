//! An axisymmetric closed chamber, a rocket-engine-like vessel, as an immersed boundary.
//!
//! The chamber is fluid where `r < r_w(z)`, for `0 <= z <= z_hi`. The port face is the grid's
//! reflecting `z = 0` boundary. The throat at `z_hi` opens onto the grid's transmissive `z`-hi
//! boundary. The wall contour is:
//!
//! - a cylinder of radius `r_c` for `0 <= z <= z_c`;
//! - a convergence from `r_c` to the throat radius `r_t` over `z_c <= z <= z_hi`, with
//!   `r_w = r_t + (r_c - r_t) g(s)`, `s = (z - z_c)/(z_hi - z_c)`, and
//!   `g(s) = (1 - curve)(1 - s) + curve sqrt(1 - s^2)`.
//!
//! `curve = 0` is a straight cone. `curve = 1` is an elliptical (convex, bulging) nose, which turns
//! the flow gradually and keeps the wall nearly parallel to the axis until close to the throat.
//! `z_c = 0` drops the cylinder, leaving a pure cone or a curved cone.
//!
//! An optional [`Bulge`] widens the wall locally by `dr`, with half-cosine ramps, to stand the wall
//! off from a blast at one station (the plug's) without widening the whole chamber.
//!
//! The solid is filled by the same ghost-cell mirror as the plate (ADR-0023 amendment). Each solid
//! cell takes the state at its mirror image across the wall, with the wall-normal velocity
//! reversed.

/// The chamber contour.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vessel {
    /// Cylinder radius [m].
    pub r_c: f64,
    /// End of the cylinder [m] (0 for no cylinder).
    pub z_c: f64,
    /// Throat plane [m].
    pub z_hi: f64,
    /// Throat radius [m].
    pub r_t: f64,
    /// 0 = straight cone, 1 = elliptical nose.
    pub curve: f64,
    /// Depth of the port-end head [m]; 0 for a flat port face. A 2:1 elliptical head is `r_c/2`.
    pub h_head: f64,
    /// Radius of the flat port opening the head closes down to [m].
    pub r_port: f64,
    /// A local widening of the wall ([`Bulge::NONE`] for none).
    pub bulge: Bulge,
}

/// A local widening: `+dr` over `z0 <= z <= z1`, ramping up over `[z0 - ramp, z0]` and down over
/// `[z1, z1 + ramp_out]` as `dr (1 - cos(pi x))/2`, so the wall and its slope stay continuous.
/// A long `ramp_out` makes the closing shoulder glancing to a blast travelling downstream.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bulge {
    /// Extra radius at the plateau [m] (>= 0).
    pub dr: f64,
    /// Plateau start [m].
    pub z0: f64,
    /// Plateau end [m].
    pub z1: f64,
    /// Opening ramp length [m].
    pub ramp: f64,
    /// Closing ramp length [m].
    pub ramp_out: f64,
}

impl Bulge {
    /// No bulge.
    pub const NONE: Self = Self {
        dr: 0.0,
        z0: 0.0,
        z1: 0.0,
        ramp: 1.0,
        ramp_out: 1.0,
    };

    /// Extra radius at `z` [m].
    #[must_use]
    pub fn extra(&self, z: f64) -> f64 {
        if self.dr == 0.0 {
            return 0.0;
        }
        let x = if z < self.z0 {
            1.0 - (self.z0 - z) / self.ramp
        } else if z > self.z1 {
            1.0 - (z - self.z1) / self.ramp_out
        } else {
            1.0
        };
        let x = x.clamp(0.0, 1.0);
        self.dr * 0.5 * (1.0 - (std::f64::consts::PI * x).cos())
    }
}

impl Vessel {
    /// Wall radius `r_w(z)` [m], bulge included.
    #[must_use]
    pub fn r_wall(&self, z: f64) -> f64 {
        self.r_base(z) + self.bulge.extra(z)
    }

    /// This contour with `bulge` added (the volume grows by the bulge's).
    #[must_use]
    pub fn with_bulge(self, bulge: Bulge) -> Self {
        Self { bulge, ..self }
    }

    /// The largest wall radius [m]: the cylinder plus any bulge.
    #[must_use]
    pub fn r_max(&self) -> f64 {
        self.r_c + self.bulge.dr
    }

    fn r_base(&self, z: f64) -> f64 {
        if self.h_head > 0.0 && z < self.h_head {
            let u = 1.0 - z.max(0.0) / self.h_head;
            return self.r_port + (self.r_c - self.r_port) * (1.0 - u * u).sqrt();
        }
        if z <= self.z_c {
            return self.r_c;
        }
        let s = ((z - self.z_c) / (self.z_hi - self.z_c)).clamp(0.0, 1.0);
        let g = (1.0 - self.curve) * (1.0 - s) + self.curve * (1.0 - s * s).sqrt();
        self.r_t + (self.r_c - self.r_t) * g
    }

    /// `dr_w/dz` by central difference.
    #[must_use]
    pub fn slope(&self, z: f64) -> f64 {
        let h = 1e-4 * self.z_hi;
        (self.r_wall(z + h) - self.r_wall(z - h)) / (2.0 * h)
    }

    /// Whether `(z, r)` is wall material.
    #[must_use]
    pub fn is_solid(&self, z: f64, r: f64) -> bool {
        r > self.r_wall(z)
    }

    /// Signed distance to the wall [m], positive in the fluid (linearized about the local slope).
    #[must_use]
    pub fn signed_distance(&self, z: f64, r: f64) -> f64 {
        let s = self.slope(z);
        (self.r_wall(z) - r) / (1.0 + s * s).sqrt()
    }

    /// Unit wall normal `(n_z, n_r)` pointing into the fluid.
    #[must_use]
    pub fn normal(&self, z: f64) -> (f64, f64) {
        let s = self.slope(z);
        let inv = 1.0 / (1.0 + s * s).sqrt();
        (s * inv, -inv)
    }

    /// Chamber volume [m^3] by midpoint quadrature.
    #[must_use]
    pub fn volume(&self) -> f64 {
        let n = 4000;
        let dz = self.z_hi / f64::from(n);
        (0..n)
            .map(|k| {
                let r = self.r_wall((f64::from(k) + 0.5) * dz);
                std::f64::consts::PI * r * r * dz
            })
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cylinder_and_cone_volume_matches_closed_form() {
        let v = Vessel {
            r_c: 1.0,
            z_c: 2.0,
            z_hi: 3.0,
            r_t: 0.2,
            curve: 0.0,
            h_head: 0.0,
            r_port: 0.0,
            bulge: Bulge::NONE,
        };
        let cone = std::f64::consts::PI / 3.0 * 1.0 * (1.0 + 0.2 + 0.04);
        let expected = std::f64::consts::PI * 2.0 + cone;
        assert!((v.volume() - expected).abs() / expected < 1e-4);
    }

    #[test]
    fn bulge_volume_matches_closed_form() {
        // A closed cylinder with a bulge whose ramps are negligibly short: a step to r_c + dr
        // over the plateau.
        let base = Vessel {
            r_c: 1.0,
            z_c: 4.0,
            z_hi: 4.0,
            r_t: 1.0,
            curve: 0.0,
            h_head: 0.0,
            r_port: 0.0,
            bulge: Bulge::NONE,
        };
        let v = Vessel {
            bulge: Bulge {
                dr: 0.5,
                z0: 1.0,
                z1: 2.0,
                ramp: 1e-6,
                ramp_out: 1e-6,
            },
            ..base
        };
        let pi = std::f64::consts::PI;
        let expected = pi * 4.0 + pi * (1.5 * 1.5 - 1.0) * 1.0;
        assert!((v.volume() - expected).abs() / expected < 1e-3);
        assert!((v.r_wall(1.5) - 1.5).abs() < 1e-12);
        assert!((v.r_wall(0.5) - 1.0).abs() < 1e-12);
        // Ramp midpoint is half the rise, and the wall is continuous across the ramp ends.
        let smooth = Vessel {
            bulge: Bulge {
                dr: 0.5,
                z0: 1.0,
                z1: 2.0,
                ramp: 0.4,
                ramp_out: 1.0,
            },
            ..base
        };
        assert!((smooth.r_wall(0.8) - 1.25).abs() < 1e-12);
        assert!((smooth.r_wall(2.5) - 1.25).abs() < 1e-12);
        assert!((smooth.r_wall(0.6 + 1e-9) - 1.0).abs() < 1e-6);
        assert!((smooth.r_wall(3.0 - 1e-9) - 1.0).abs() < 1e-6);
    }

    // SAFE: test radii and spacings are small positive numbers, so the cell count fits a usize.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    fn chamber_grid(v: Vessel, n: usize, zhi_reflect: bool) -> crate::kernel::Grid2D {
        use crate::kernel::{Bc, Grid2D};
        let dz = v.z_hi / n as f64;
        let nr = (v.r_max() / dz).ceil() as usize + 3;
        let mut g = Grid2D::new(n, nr, dz, dz, 1.4);
        g.set_axisymmetric(true);
        g.bc_zlo = Bc::Reflect;
        g.bc_zhi = if zhi_reflect {
            Bc::Reflect
        } else {
            Bc::Transmissive
        };
        g.bc_rlo = Bc::Reflect;
        g.bc_rhi = Bc::Reflect;
        g.set_vessel(Some(v));
        g
    }

    #[test]
    fn gas_at_rest_stays_at_rest_in_a_curved_chamber() {
        use crate::state::Prim;
        let v = Vessel {
            r_c: 1.0,
            z_c: 0.5,
            z_hi: 2.0,
            r_t: 0.2,
            curve: 0.6,
            h_head: 0.5,
            r_port: 0.1,
            bulge: Bulge::NONE,
        };
        let mut g = chamber_grid(v, 80, true);
        g.init(|_, _| Prim::new(1.0, 0.0, 0.0, 1.0));
        g.run_to(0.5);
        for (_, _, p) in g.vessel_wall_pressures() {
            assert!((p - 1.0).abs() < 1e-10, "wall pressure drifted to {p}");
        }
    }

    #[test]
    fn gas_at_rest_stays_at_rest_in_a_bulged_chamber() {
        use crate::state::Prim;
        let v = Vessel {
            r_c: 1.0,
            z_c: 1.4,
            z_hi: 2.0,
            r_t: 0.2,
            curve: 0.0,
            h_head: 0.3,
            r_port: 0.1,
            bulge: Bulge {
                dr: 0.4,
                z0: 0.7,
                z1: 0.9,
                ramp: 0.3,
                ramp_out: 0.3,
            },
        };
        let mut g = chamber_grid(v, 80, true);
        g.init(|_, _| Prim::new(1.0, 0.0, 0.0, 1.0));
        g.run_to(0.5);
        for (_, _, p) in g.vessel_wall_pressures() {
            assert!((p - 1.0).abs() < 1e-10, "wall pressure drifted to {p}");
        }
    }

    #[test]
    fn closed_cylinder_blast_keeps_its_mass() {
        use crate::state::Prim;
        // A cylinder closed at both ends: z_c = z_hi makes r_w = r_c everywhere.
        let v = Vessel {
            r_c: 1.0,
            z_c: 2.0,
            z_hi: 2.0,
            r_t: 1.0,
            curve: 0.0,
            h_head: 0.0,
            r_port: 0.0,
            bulge: Bulge::NONE,
        };
        let mut g = chamber_grid(v, 80, true);
        g.init(|iz, ir| {
            let (z, r) = ((iz as f64 + 0.5) * 0.025, (ir as f64 + 0.5) * 0.025);
            let hot = (z - 1.0).hypot(r) < 0.1;
            Prim::new(1.0, 0.0, 0.0, if hot { 100.0 } else { 1.0 })
        });
        let mass = |g: &crate::kernel::Grid2D| -> f64 {
            let mut m = 0.0;
            for iz in 0..g.nz() {
                for ir in 0..g.nr() {
                    if !v.is_solid((iz as f64 + 0.5) * 0.025, (ir as f64 + 0.5) * 0.025) {
                        m += g.prim(iz, ir).rho * g.cell_volume(ir);
                    }
                }
            }
            m
        };
        let m0 = mass(&g);
        g.run_to(0.6);
        assert!((mass(&g) - m0).abs() / m0 < 0.02);
    }

    /// The multithreaded sweep must reproduce the serial one bit for bit: a blast in a domed,
    /// curved chamber with an open throat, run both ways, compared cell by cell with `==`.
    #[test]
    fn parallel_sweeps_are_bit_identical_to_serial() {
        use crate::state::Prim;
        let v = Vessel {
            r_c: 1.0,
            z_c: 1.2,
            z_hi: 2.0,
            r_t: 0.2,
            curve: 1.0,
            h_head: 0.5,
            r_port: 0.05,
            bulge: Bulge::NONE,
        };
        let run = |parallel: bool| {
            let mut g = chamber_grid(v, 60, false);
            g.set_parallel(parallel);
            let cell = v.z_hi / 60.0;
            g.init(|iz, ir| {
                let (z, r) = ((iz as f64 + 0.5) * cell, (ir as f64 + 0.5) * cell);
                let hot = (z - 0.8).hypot(r) < 0.12;
                Prim::new(1.0, 0.0, 0.0, if hot { 500.0 } else { 1.0 })
            });
            g.run_to(0.15);
            g
        };
        let (serial, parallel) = (run(false), run(true));
        for iz in 0..serial.nz() {
            for ir in 0..serial.nr() {
                assert!(
                    serial.cons(iz, ir) == parallel.cons(iz, ir),
                    "cell ({iz}, {ir}) differs between serial and parallel runs"
                );
            }
        }
    }

    #[test]
    fn normal_points_inward_on_the_cylinder() {
        let v = Vessel {
            r_c: 1.0,
            z_c: 2.0,
            z_hi: 3.0,
            r_t: 0.2,
            curve: 1.0,
            h_head: 0.0,
            r_port: 0.0,
            bulge: Bulge::NONE,
        };
        let (nz, nr) = v.normal(1.0);
        assert!(nz.abs() < 1e-9 && (nr + 1.0).abs() < 1e-9);
        assert!(v.signed_distance(1.0, 0.9) > 0.0 && v.is_solid(1.0, 1.1));
    }
}
