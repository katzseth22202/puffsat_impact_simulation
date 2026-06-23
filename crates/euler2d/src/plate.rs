//! Plate surface descriptors for the immersed-boundary reflecting wall (ADR-0023 amendment, D4).
//!
//! The shallow-concave plate's surface `z = z_s(r)` cuts diagonally across the square `(z, r)` mesh,
//! so it is imposed as a **ghost-cell immersed boundary** (a true-normal mirror) rather than the
//! grid-aligned `z = 0` reflecting BC the flat plate uses. A staircase of full square cells would
//! bias the rebound angle — and `eta_capture` *is* a rebound-angle measurement — so capturing the
//! true surface normal is the point (ADR-0023). A [`PlateProfile`] answers the three questions the
//! immersed-boundary pass and the axial wall-impulse integral need: is a point inside the solid
//! plate, what is the signed distance to the surface (negative inside), and what is the unit normal
//! `n̂` pointing into the fluid.
//!
//! Two shapes:
//! - [`PlateProfile::InclinedPlane`] — a flat wall `z = z0 + slope·r` spanning the whole domain (no
//!   edge); the immersed-boundary acceptance tests (free-slip tangency, specular normal rebound) use
//!   it because a constant tilt has a constant analytic normal.
//! - [`PlateProfile::Dish`] — the axisymmetric shallow-concave plate `z = z0 + depth·(r/r_plate)²`
//!   for `r ≤ r_plate` (parabolic — the shallow limit; a spherical cap is the alternative).
//!   `depth = (d/D)·2·r_plate` from the depth-to-diameter ratio (ADR-0021); `z0` raises the whole
//!   dish a few cells off the domain floor so a solid layer always underlies it. Gas past the rim
//!   (`r > r_plate`) is over no plate and escapes (§7).

/// A reflecting plate surface `z = z_s(r)`, imposed as an immersed boundary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PlateProfile {
    /// A flat wall `z = z0 + slope·r`, unbounded in `r` (used by the planar acceptance tests).
    InclinedPlane {
        /// Surface height on the axis `r = 0`.
        z0: f64,
        /// `dz_s/dr` (constant).
        slope: f64,
    },
    /// The axisymmetric shallow-concave dish `z = z0 + depth·(r/r_plate)²` for `r ≤ r_plate`.
    Dish {
        /// Plate radius (the rim; gas past it escapes).
        r_plate: f64,
        /// Surface height on the axis (the dish floor, raised off the domain bottom).
        z0: f64,
        /// Rim-to-floor depth `d` (so `depth = (d/D)·2·r_plate`).
        depth: f64,
    },
}

impl PlateProfile {
    /// Surface height `z_s(r)`. Past the dish rim the parabola is clamped at the rim value (the
    /// plate does not extend there — [`Self::covers`] gates the solid region).
    #[must_use]
    pub fn z_surface(&self, r: f64) -> f64 {
        match *self {
            Self::InclinedPlane { z0, slope } => z0 + slope * r,
            Self::Dish { r_plate, z0, depth } => {
                let rr = (r / r_plate).min(1.0);
                z0 + depth * rr * rr
            }
        }
    }

    /// Surface slope `dz_s/dr` at radius `r`.
    #[must_use]
    pub fn slope(&self, r: f64) -> f64 {
        match *self {
            Self::InclinedPlane { slope, .. } => slope,
            Self::Dish { r_plate, depth, .. } => {
                if r >= r_plate {
                    0.0
                } else {
                    2.0 * depth * r / (r_plate * r_plate)
                }
            }
        }
    }

    /// Whether the plate is present at radius `r` (the dish ends at its rim; the plane is unbounded).
    #[must_use]
    pub fn covers(&self, r: f64) -> bool {
        match *self {
            Self::InclinedPlane { .. } => true,
            Self::Dish { r_plate, .. } => r <= r_plate,
        }
    }

    /// Outward unit normal `(n_z, n_r)` pointing into the fluid (the `+z` side). For the surface
    /// `F = z − z_s(r) = 0`, `∇F = (1, −z_s′)`, normalized.
    #[must_use]
    pub fn normal(&self, r: f64) -> (f64, f64) {
        let s = self.slope(r);
        let inv = 1.0 / (1.0 + s * s).sqrt();
        (inv, -s * inv)
    }

    /// Whether the point `(z, r)` lies inside the solid plate.
    #[must_use]
    pub fn is_solid(&self, z: f64, r: f64) -> bool {
        self.covers(r) && z < self.z_surface(r)
    }

    /// Signed perpendicular distance to the surface (linearized about the foot of the normal):
    /// the vertical gap `z − z_s(r)` projected onto the normal, negative inside the solid. Exact for
    /// the plane and a shallow-curve approximation for the dish.
    #[must_use]
    pub fn signed_distance(&self, z: f64, r: f64) -> f64 {
        let (nz, _) = self.normal(r);
        (z - self.z_surface(r)) * nz
    }
}

#[cfg(test)]
mod tests {
    use super::PlateProfile;
    use approx::assert_relative_eq;

    #[test]
    fn dish_floor_and_rim_heights() {
        let p = PlateProfile::Dish {
            r_plate: 2.0,
            z0: 0.1,
            depth: 0.6,
        };
        assert_relative_eq!(p.z_surface(0.0), 0.1, epsilon = 1e-14); // floor on the axis
        assert_relative_eq!(p.z_surface(2.0), 0.7, epsilon = 1e-14); // rim = z0 + depth
        assert_relative_eq!(p.z_surface(1.0), 0.1 + 0.6 * 0.25, epsilon = 1e-14); // parabolic
    }

    #[test]
    fn solid_below_surface_and_only_within_rim() {
        let p = PlateProfile::Dish {
            r_plate: 2.0,
            z0: 0.1,
            depth: 0.6,
        };
        assert!(p.is_solid(0.05, 0.0)); // below the floor, on the axis
        assert!(!p.is_solid(0.5, 0.0)); // above the floor
        assert!(!p.is_solid(0.05, 3.0)); // past the rim — no plate there
    }

    #[test]
    fn normal_is_unit_and_tilts_toward_axis_on_a_rising_dish() {
        let p = PlateProfile::Dish {
            r_plate: 2.0,
            z0: 0.0,
            depth: 0.6,
        };
        let (nz, nr) = p.normal(1.5);
        assert_relative_eq!(nz * nz + nr * nr, 1.0, epsilon = 1e-14);
        assert!(nz > 0.0, "normal points into the fluid (+z)");
        assert!(
            nr < 0.0,
            "a rising dish tilts its normal toward the axis (−r)"
        );
    }

    #[test]
    fn signed_distance_sign_and_plane_exactness() {
        // For a flat inclined plane the linearized distance is exact: a point at vertical gap Δz
        // above the surface sits at perpendicular distance Δz·n_z.
        let p = PlateProfile::InclinedPlane {
            z0: 0.2,
            slope: 0.5,
        };
        let (nz, _) = p.normal(0.0);
        let zs = p.z_surface(1.0); // = 0.7
        assert_relative_eq!(p.signed_distance(zs + 0.3, 1.0), 0.3 * nz, epsilon = 1e-14);
        assert!(
            p.signed_distance(0.0, 1.0) < 0.0,
            "below the surface ⇒ negative"
        );
    }
}
