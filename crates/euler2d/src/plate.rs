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
    /// A [`Self::Dish`] floor with a wall rising from its rim (ADR-0055): the spray plate's skirt
    /// or flared bell. The wall's inner face runs `r_in(z) = r_plate + flare·(z − z_rim)` from the
    /// rim height `z_rim = z0 + depth` up to `z_rim + skirt_height`, and the wall is `thickness`
    /// thick in `r`. Below the rim it is a straight continuation, so floor and wall join. A zero
    /// `skirt_height` is exactly the dish.
    Cup {
        /// Floor (dish) radius, where the wall starts.
        r_plate: f64,
        /// Floor height on the axis.
        z0: f64,
        /// Floor rim-to-axis depth.
        depth: f64,
        /// Wall height above the rim.
        skirt_height: f64,
        /// Outward flare of the wall, `dr/dz` (0 is a straight cylinder).
        flare: f64,
        /// Wall thickness in `r` (at least two cells, so ghost cells have solid behind them).
        thickness: f64,
    },
}

impl PlateProfile {
    /// Surface height `z_s(r)`. Past the dish rim the parabola is clamped at the rim value (the
    /// plate does not extend there — [`Self::covers`] gates the solid region).
    #[must_use]
    pub fn z_surface(&self, r: f64) -> f64 {
        match *self {
            Self::InclinedPlane { z0, slope } => z0 + slope * r,
            Self::Dish { r_plate, z0, depth }
            | Self::Cup {
                r_plate, z0, depth, ..
            } => {
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
            Self::Dish { r_plate, depth, .. } | Self::Cup { r_plate, depth, .. } => {
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
            Self::Cup {
                r_plate,
                skirt_height,
                flare,
                thickness,
                ..
            } => {
                if skirt_height > 0.0 {
                    r <= r_plate + flare.max(0.0) * skirt_height + thickness
                } else {
                    r <= r_plate
                }
            }
        }
    }

    /// Outward unit normal `(n_z, n_r)` pointing into the fluid, from the plate face **nearest**
    /// the point. For the top surface `F = z − z_s(r) = 0`, `∇F = (1, −z_s′)`, normalized. The
    /// dish's solid body also ends at its rim (`r = r_plate`): a point nearer that vertical side
    /// face than the top surface takes the radial normal `(0, 1)` instead. Ignoring the side face
    /// (the pre-fix behavior) mirrored rim-adjacent solid cells across the *top* surface, feeding
    /// spurious radial fluxes into the fluid past the rim — a bounded error at M ≲ 20 that becomes
    /// a self-exciting energy source at the rim corner for very strong shocks (found at M = 40).
    #[must_use]
    pub fn normal(&self, z: f64, r: f64) -> (f64, f64) {
        match *self {
            Self::InclinedPlane { .. } => self.top_normal(r),
            Self::Dish { r_plate, .. } => self.dish_normal(z, r, r_plate),
            Self::Cup { r_plate, .. } => {
                if self.floor_solid(z, r) || !self.wall_solid(z, r) {
                    self.dish_normal(z, r, r_plate)
                } else {
                    self.wall_normal(z, r)
                }
            }
        }
    }

    /// The dish's normal: its top surface, or its rim side face where that is nearer.
    fn dish_normal(&self, z: f64, r: f64, r_plate: f64) -> (f64, f64) {
        let (d_top, d_side) = (self.top_distance(z, r), r - r_plate);
        if d_side > d_top {
            (0.0, 1.0)
        } else {
            self.top_normal(r)
        }
    }

    /// The cup wall's geometry `(r_in(z), wall top, flare, thickness)`; `None` without a wall.
    fn wall(&self, z: f64) -> Option<(f64, f64, f64, f64)> {
        match *self {
            Self::Cup {
                r_plate,
                z0,
                depth,
                skirt_height,
                flare,
                thickness,
            } if skirt_height > 0.0 => {
                let z_rim = z0 + depth;
                let r_in = r_plate + flare * (z - z_rim).max(0.0);
                Some((r_in, z_rim + skirt_height, flare, thickness))
            }
            _ => None,
        }
    }

    /// The dish floor's height span `(axis, rim)`. `None` for a flat floor or a plane.
    #[must_use]
    pub fn floor_span(&self) -> Option<(f64, f64)> {
        match *self {
            Self::Dish { z0, depth, .. } | Self::Cup { z0, depth, .. } if depth > 0.0 => {
                Some((z0, z0 + depth))
            }
            _ => None,
        }
    }

    /// The wall's height span `(rim, top)`, from the floor's rim to the wall's lip. `None` without
    /// a wall.
    #[must_use]
    pub fn wall_span(&self) -> Option<(f64, f64)> {
        match *self {
            Self::Cup {
                z0,
                depth,
                skirt_height,
                ..
            } if skirt_height > 0.0 => Some((z0 + depth, z0 + depth + skirt_height)),
            _ => None,
        }
    }

    /// The floor part of the solid (the dish).
    fn floor_solid(&self, z: f64, r: f64) -> bool {
        let r_plate = match *self {
            Self::Dish { r_plate, .. } | Self::Cup { r_plate, .. } => r_plate,
            Self::InclinedPlane { .. } => return z < self.z_surface(r),
        };
        r <= r_plate && z < self.z_surface(r)
    }

    /// The wall's signed distance (negative inside): the slab between its inner and outer faces,
    /// capped by its lip. `None` without a wall.
    fn wall_distance(&self, z: f64, r: f64) -> Option<(f64, f64, f64)> {
        let (r_in, z_top, flare, thickness) = self.wall(z)?;
        let scale = 1.0 / (1.0 + flare * flare).sqrt();
        let d_inner = (r_in - r) * scale;
        let d_outer = (r - r_in - thickness) * scale;
        Some((d_inner, d_outer, z - z_top))
    }

    /// The wall part of the solid.
    fn wall_solid(&self, z: f64, r: f64) -> bool {
        self.wall_distance(z, r)
            .is_some_and(|(a, b, c)| a.max(b).max(c) < 0.0)
    }

    /// The normal of the wall face nearest a point inside the wall: inner (toward the axis, tilted
    /// up by the flare), outer, or the lip.
    fn wall_normal(&self, z: f64, r: f64) -> (f64, f64) {
        let Some((d_inner, d_outer, d_top)) = self.wall_distance(z, r) else {
            return (1.0, 0.0);
        };
        let flare = self.wall(z).map_or(0.0, |w| w.2);
        let scale = 1.0 / (1.0 + flare * flare).sqrt();
        if d_top >= d_inner && d_top >= d_outer {
            (1.0, 0.0)
        } else if d_inner >= d_outer {
            (flare * scale, -scale)
        } else {
            (-flare * scale, scale)
        }
    }

    /// Unit normal of the top surface `z = z_s(r)` (pointing into the fluid above).
    fn top_normal(&self, r: f64) -> (f64, f64) {
        let s = self.slope(r);
        let inv = 1.0 / (1.0 + s * s).sqrt();
        (inv, -s * inv)
    }

    /// Signed perpendicular distance to the top surface (linearized about the foot of the normal):
    /// the vertical gap `z − z_s(r)` projected onto the normal, negative below the surface. Exact
    /// for the plane and a shallow-curve approximation for the dish.
    fn top_distance(&self, z: f64, r: f64) -> f64 {
        let (nz, _) = self.top_normal(r);
        (z - self.z_surface(r)) * nz
    }

    /// Whether the point `(z, r)` lies inside the solid plate.
    #[must_use]
    pub fn is_solid(&self, z: f64, r: f64) -> bool {
        match *self {
            Self::Cup { .. } => self.floor_solid(z, r) || self.wall_solid(z, r),
            _ => self.covers(r) && z < self.z_surface(r),
        }
    }

    /// Signed distance to the solid's boundary, negative inside. The dish is the intersection of
    /// two half-spaces — below the top surface *and* within the rim radius — so its signed distance
    /// is the max of the two face distances (exact away from the rim corner). The plane has only
    /// the top face.
    #[must_use]
    pub fn signed_distance(&self, z: f64, r: f64) -> f64 {
        match *self {
            Self::InclinedPlane { .. } => self.top_distance(z, r),
            Self::Dish { r_plate, .. } => self.top_distance(z, r).max(r - r_plate),
            Self::Cup { r_plate, .. } => {
                let floor = self.top_distance(z, r).max(r - r_plate);
                match self.wall_distance(z, r) {
                    Some((a, b, c)) => floor.min(a.max(b).max(c)),
                    None => floor,
                }
            }
        }
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
        // A point just under the surface at mid-radius: the top face is nearest.
        let z = p.z_surface(1.5) - 0.01;
        let (nz, nr) = p.normal(z, 1.5);
        assert_relative_eq!(nz * nz + nr * nr, 1.0, epsilon = 1e-14);
        assert!(nz > 0.0, "normal points into the fluid (+z)");
        assert!(
            nr < 0.0,
            "a rising dish tilts its normal toward the axis (−r)"
        );
    }

    #[test]
    fn rim_side_face_owns_nearby_solid_cells() {
        // The dish body ends at its rim: a solid cell just inside `r_plate` but far below the top
        // surface must mirror across the vertical side face (radial normal, side distance), not
        // across the faraway top surface — the M=40 rim-corner blow-up was exactly this.
        let p = PlateProfile::Dish {
            r_plate: 2.0,
            z0: 0.1,
            depth: 0.6,
        };
        let (z, r) = (0.3, 1.95); // deep under the rim: top face ~0.32 away, side face 0.05 away
        let (nz, nr) = p.normal(z, r);
        assert_relative_eq!(nz, 0.0, epsilon = 1e-14);
        assert_relative_eq!(nr, 1.0, epsilon = 1e-14);
        assert_relative_eq!(p.signed_distance(z, r), -0.05, epsilon = 1e-12);

        // On the axis the side face is 2.0 away and the floor 0.05 above: the top face owns it.
        let (nz, nr) = p.normal(0.05, 0.0);
        assert_relative_eq!(nz, 1.0, epsilon = 1e-14);
        assert_relative_eq!(nr, 0.0, epsilon = 1e-14);
        assert_relative_eq!(p.signed_distance(0.05, 0.0), -0.05, epsilon = 1e-12);
    }

    #[test]
    fn signed_distance_sign_and_plane_exactness() {
        // For a flat inclined plane the linearized distance is exact: a point at vertical gap Δz
        // above the surface sits at perpendicular distance Δz·n_z.
        let p = PlateProfile::InclinedPlane {
            z0: 0.2,
            slope: 0.5,
        };
        let (nz, _) = p.normal(0.0, 0.0);
        let zs = p.z_surface(1.0); // = 0.7
        assert_relative_eq!(p.signed_distance(zs + 0.3, 1.0), 0.3 * nz, epsilon = 1e-14);
        assert!(
            p.signed_distance(0.0, 1.0) < 0.0,
            "below the surface ⇒ negative"
        );
    }

    fn cup(skirt_height: f64, flare: f64) -> PlateProfile {
        PlateProfile::Cup {
            r_plate: 2.0,
            z0: 0.1,
            depth: 0.2,
            skirt_height,
            flare,
            thickness: 0.2,
        }
    }

    #[test]
    fn cup_wall_is_solid_and_its_interior_is_fluid() {
        let c = cup(1.0, 0.0);
        // Rim at z = 0.3; the wall spans r in [2.0, 2.2] up to z = 1.3.
        assert!(c.is_solid(0.8, 2.1)); // inside the wall
        assert!(!c.is_solid(0.8, 1.9)); // inside the cup
        assert!(!c.is_solid(0.8, 2.3)); // outside the wall
        assert!(!c.is_solid(1.4, 2.1)); // above the lip
        assert!(c.is_solid(0.05, 1.0)); // the floor still there
    }

    #[test]
    fn cup_inner_face_points_at_the_axis_and_up_when_flared() {
        let straight = cup(1.0, 0.0);
        let (nz, nr) = straight.normal(0.8, 2.01);
        assert_relative_eq!(nz, 0.0, epsilon = 1e-12);
        assert_relative_eq!(nr, -1.0, epsilon = 1e-12);
        let flared = cup(1.0, 0.5);
        // Inner face at z = 0.8 sits at r = 2.0 + 0.5 * 0.5 = 2.25.
        let (nz, nr) = flared.normal(0.8, 2.26);
        assert!(nz > 0.0 && nr < 0.0, "({nz}, {nr})");
        assert_relative_eq!(nz * nz + nr * nr, 1.0, epsilon = 1e-12);
    }

    #[test]
    fn a_cup_without_a_wall_is_the_dish() {
        let c = cup(0.0, 0.0);
        let d = PlateProfile::Dish {
            r_plate: 2.0,
            z0: 0.1,
            depth: 0.2,
        };
        for &(z, r) in &[
            (0.05, 0.5),
            (0.25, 1.9),
            (0.4, 2.05),
            (0.2, 2.5),
            (0.15, 1.0),
        ] {
            assert_eq!(c.is_solid(z, r), d.is_solid(z, r), "solid at ({z}, {r})");
            assert_eq!(c.covers(r), d.covers(r), "covers at {r}");
            if d.is_solid(z, r) {
                assert_eq!(c.normal(z, r), d.normal(z, r), "normal at ({z}, {r})");
                assert_eq!(c.signed_distance(z, r), d.signed_distance(z, r));
            }
        }
    }
}
