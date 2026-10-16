//! The canted departure plate (X-ray study, Result 12): a wide impactor sheet sweeping *past* a
//! flat plate, in planar `(x, z)` geometry.
//!
//! In the 3S Earth -> Jupiter -> Sun -> Earth cycle the departure burn's thrust axis meets the
//! returning stream at ~100 degrees (the calcs repo's `departure_nozzle_ledger`). The stream moves
//! almost parallel to the plate, ~10 degrees *away* from its face. A wide Kapton sheet strikes a
//! polyethylene strip held off a flat plate, and the merged fireball explodes beside the plate
//! while sliding along it. A frictionless flat plate takes no sideways force, so only the normal
//! impulse counts. Conserving the tangential momentum and spending all the energy on the normal
//! one bounds it at
//!
//! `J_n <= m w (cos theta + sqrt(1 + k - sin^2 theta))`,
//!
//! the **canted ceiling** (`theta` from the thrust axis; `theta = 0` is the overtake ceiling).
//!
//! Geometry: the plate lies along `z = 0` for `0 <= x <= plate_length` (the grid's `r` direction
//! is `x`, Cartesian), gas above it. Thrust is `-z`. The stream direction is
//! `s = (sin theta, -cos theta)` in `(x, z)`: for `theta > 90` it climbs away from the plate. Both
//! bodies are rectangles aligned with `s`: the impactor a thin slab of unit mass (per unit depth)
//! at speed 1, the charge a block of mass `k` at rest just downstream of it, so their collision is
//! a clean 1-D slab collision along `s`. Units: mass, speed and impulse per unit depth, read
//! directly in `m w`.
//!
//! The sheet spans exactly the across-stream extent of what it is aimed at, so all of its unit
//! mass strikes. (Until 2026-10-15 it was always `width` wide about the charge's centre, and ~41%
//! of it passed beside the two half-disks lying on a face.) A [`PlateLayer`] splits the charge
//! across the stream into a dense plate-side layer and the foam above it; the sheet is then aimed
//! at the foam alone, so the layer is driven only by the foam's fireball.

use crate::kernel::{Bc, Grid2D};
use crate::state::Prim;

/// The sheet's Mach number at the start: its pressure sets both bodies' cold pressure.
const COLD_MACH: f64 = 20.0;

/// The charge's cross-section (planar, so a "hemisphere" is a half-cylinder).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChargeShape {
    /// A rectangle `charge_length` along the stream and `width` across it.
    Block,
    /// A disk of diameter `width`.
    Round,
    /// A half-disk of radius `width/2` with its flat face toward the incoming sheet.
    HalfFlatUpstream,
    /// A half-disk with its flat face toward the plate, dome up.
    HalfDomeUp,
    /// A half-disk with its dome toward the plate.
    HalfDomeDown,
}

/// A dense layer on the charge's plate side (a solid methane face under polyethylene foam): the
/// charge cells nearest the plate across the stream, holding `mass_share` of `k` at
/// `density_ratio` times the foam's density.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlateLayer {
    pub mass_share: f64,
    pub density_ratio: f64,
}

/// One canted strip run.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CantedStrip {
    /// Angle between the stream and the thrust axis [deg] (~101 for the 3S departure).
    pub theta_deg: f64,
    /// Charge-to-impactor mass ratio.
    pub k: f64,
    /// Height of the charge's centre above the plate [m].
    pub standoff: f64,
    /// Distance of the charge's centre downstream of the plate's upstream edge [m].
    pub charge_x: f64,
    /// Plate length along the stream [m].
    pub plate_length: f64,
    /// Width of both bodies across the stream [m], and the charge's length along it [m].
    pub width: f64,
    pub charge_length: f64,
    /// Cell size [m] and effective gamma.
    pub dx: f64,
    pub gamma: f64,
    /// Ambient density as a fraction of the charge's.
    pub ambient: f64,
    /// `Some(R)` replaces the sheet and charge with the merged fireball they make: a hot round
    /// column of radius R and mass 1 + k, holding the collision's heat `k/(1+k)/2` and drifting
    /// along the stream at `1/(1+k)` (the planar point fireball, as Results 10-11 use in 2-D).
    pub point_radius: Option<f64>,
    /// The charge's cross-section (the slab and square runs are [`ChargeShape::Block`]).
    pub shape: ChargeShape,
    /// `Some` splits the charge into a plate-side layer and the struck foam (not with
    /// `point_radius`).
    pub layer: Option<PlateLayer>,
}

/// What the plate took.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CantedResult {
    /// Normal impulse on the plate per `m w` (lossless gamma gas).
    pub impulse: f64,
    /// The canted ceiling `cos theta + sqrt(1 + k - sin^2 theta)`.
    pub ceiling: f64,
    /// Where the impulse acts: its centroid along the plate [m], and that centroid's offset
    /// downstream of the charge's foot (the torque-free placement of the ship's centre of mass).
    pub centroid_x: f64,
    pub centroid_offset: f64,
    /// Peak plate pressure per `m w^2` (per unit depth, per m^2) and where it fell [m].
    pub peak_pressure: f64,
    pub peak_x: f64,
    /// Ambient mass as a share of the merged mass.
    pub ambient_mass_share: f64,
    /// The layer's density over the foam's as gridded (the masses are exact, so this differs a
    /// little from the requested ratio), and the layer's share of `k`; 1 and 0 with no layer.
    pub layer_density_ratio: f64,
    pub layer_mass_share: f64,
    pub steps: usize,
}

impl CantedStrip {
    fn s(&self) -> (f64, f64) {
        let t = self.theta_deg.to_radians();
        (t.sin(), -t.cos())
    }

    /// The canted ceiling for this strip.
    #[must_use]
    pub fn ceiling(&self) -> f64 {
        let t = self.theta_deg.to_radians();
        t.cos() + (1.0 + self.k - t.sin() * t.sin()).sqrt()
    }
}

/// Is `(x, z)` inside the rectangle centred at `c`, `along` long in `s` and `across` wide?
fn in_rect(x: f64, z: f64, c: (f64, f64), s: (f64, f64), along: f64, across: f64) -> bool {
    let (dx, dz) = (x - c.0, z - c.1);
    let a = dx * s.0 + dz * s.1;
    let b = -dx * s.1 + dz * s.0;
    a.abs() <= 0.5 * along && b.abs() <= 0.5 * across
}

/// Is `(x, z)` inside the charge of `c` centred at `centre`?
fn in_charge(x: f64, z: f64, centre: (f64, f64), s: (f64, f64), c: &CantedStrip) -> bool {
    let r = 0.5 * c.width;
    let (dx, dz) = (x - centre.0, z - centre.1);
    let along = dx * s.0 + dz * s.1;
    let inside = dx.hypot(dz) < r;
    match c.shape {
        ChargeShape::Block => in_rect(x, z, centre, s, c.charge_length, c.width),
        ChargeShape::Round => inside,
        ChargeShape::HalfFlatUpstream => inside && along >= 0.0,
        ChargeShape::HalfDomeUp => inside && dz >= 0.0,
        ChargeShape::HalfDomeDown => inside && dz <= 0.0,
    }
}

/// Run one strip: the sheet strikes the charge, and the merged gas expands past the plate.
#[must_use]
pub fn run_canted_strip(c: &CantedStrip) -> CantedResult {
    let s = c.s();
    let disk = 4.0 * c.dx;
    let gap = 2.0 * c.dx;
    // Domain: the plate spans [0, plate_length] from the open upstream edge, with 4 m of open
    // floor downstream. The sheet starts just upstream of the charge, inside the plate's span.
    let x0 = 0.0;
    let x_len = c.plate_length + 4.0;
    let z_len = c.standoff + 0.5 * c.width + 8.0;
    let mut nr = 0usize;
    while (nr as f64) * c.dx < x_len {
        nr += 1;
    }
    let mut nz = 0usize;
    while (nz as f64) * c.dx < z_len {
        nz += 1;
    }
    let mut g = Grid2D::new(nz, nr, c.dx, c.dx, c.gamma);
    g.set_axisymmetric(false);
    g.bc_rlo = Bc::Transmissive;
    g.bc_rhi = Bc::Transmissive;
    g.bc_zhi = Bc::Transmissive;
    g.bc_zlo = Bc::Transmissive;
    // The kernel's finite plate: the z-lo boundary reflects for x <= plate_length, open beyond.
    g.set_plate_radius(Some(c.plate_length));
    let plate = (x0, x0 + c.plate_length);
    let charge = (x0 + c.charge_x, c.standoff);
    let back = 0.5 * c.charge_length + gap + 0.5 * disk;
    let centre = |iz: usize, ir: usize| ((ir as f64 + 0.5) * c.dx, (iz as f64 + 0.5) * c.dx);
    // Across-stream coordinate about the charge's centre (positive away from the plate).
    let across = |x: f64, z: f64| -(x - charge.0) * s.1 + (z - charge.1) * s.0;
    let point = c.point_radius.is_some();
    assert!(
        !(point && c.layer.is_some()),
        "a point fireball has no layers"
    );
    let mut charge_b: Vec<f64> = Vec::new();
    if !point {
        for iz in 0..nz {
            for ir in 0..nr {
                let (x, z) = centre(iz, ir);
                if in_charge(x, z, charge, s, c) {
                    charge_b.push(across(x, z));
                }
            }
        }
        charge_b.sort_by(f64::total_cmp);
    }
    // The layer is the `n_layer` charge cells nearest the plate across the stream; `n_layer`
    // makes `r n / (r n + N - n)` the requested mass share at density ratio `r`.
    let b_cut = c.layer.map_or(f64::NEG_INFINITY, |l| {
        let n = charge_b.len() as f64;
        let want = l.mass_share * n / (l.density_ratio * (1.0 - l.mass_share) + l.mass_share);
        let mut n_layer = 1usize;
        while (n_layer as f64) + 0.5 < want && n_layer + 1 < charge_b.len() {
            n_layer += 1;
        }
        0.5 * (charge_b[n_layer - 1] + charge_b[n_layer])
    });
    // The sheet spans the struck cells' across extent, kept at the nominal edge wherever the
    // struck part reaches it (so full-width charges are unchanged).
    let (sheet_lo, sheet_hi) = if point {
        (-0.5 * c.width, 0.5 * c.width)
    } else {
        let struck: Vec<f64> = charge_b.iter().copied().filter(|&b| b > b_cut).collect();
        let (lo, hi) = (struck[0], struck[struck.len() - 1]);
        (
            if lo - c.dx <= -0.5 * c.width {
                -0.5 * c.width
            } else {
                lo - 0.5 * c.dx
            },
            if hi + c.dx >= 0.5 * c.width {
                0.5 * c.width
            } else {
                hi + 0.5 * c.dx
            },
        )
    };
    let mid = 0.5 * (sheet_lo + sheet_hi);
    let sheet = (
        charge.0 - back * s.0 - mid * s.1,
        charge.1 - back * s.1 + mid * s.0,
    );
    let sheet_width = sheet_hi - sheet_lo;
    let region = |iz: usize, ir: usize| -> u8 {
        let (x, z) = centre(iz, ir);
        if let Some(radius) = c.point_radius {
            return if (x - charge.0).hypot(z - charge.1) < radius {
                3
            } else {
                0
            };
        }
        if in_charge(x, z, charge, s, c) {
            if across(x, z) < b_cut { 4 } else { 2 }
        } else if in_rect(x, z, sheet, s, disk, sheet_width) {
            1
        } else {
            0
        }
    };
    let cell = c.dx * c.dx;
    let (mut v_sheet, mut v_foam, mut v_layer, mut v_all) = (0.0, 0.0, 0.0, 0.0);
    for iz in 0..nz {
        for ir in 0..nr {
            v_all += cell;
            match region(iz, ir) {
                1 => v_sheet += cell,
                2 | 3 => v_foam += cell,
                4 => v_layer += cell,
                _ => {}
            }
        }
    }
    let v_charge = v_foam + v_layer;
    // The point fireball has no sheet; its mass 1 + k sits in the charge's region.
    let rho_sheet = if point { 1.0 } else { 1.0 / v_sheet };
    let rho_charge = if point {
        (1.0 + c.k) / v_charge
    } else {
        c.k / v_charge
    };
    let layer_share = c.layer.map_or(0.0, |l| l.mass_share);
    let (rho_foam, rho_layer) = if c.layer.is_some() {
        (
            c.k * (1.0 - layer_share) / v_foam,
            c.k * layer_share / v_layer,
        )
    } else {
        (rho_charge, rho_charge)
    };
    let drift = 1.0 / (1.0 + c.k);
    let p_hot = (c.gamma - 1.0) * 0.5 * c.k / (1.0 + c.k) / v_charge;
    // Both bodies start cold at the sheet's pressure (Mach `COLD_MACH` sheet). Holding the charge
    // at the sheet's *temperature* instead would give it k times the sheet's specific heat, which
    // at k = 30 rivals the collision's own heat.
    let p_cold = rho_sheet / (c.gamma * COLD_MACH * COLD_MACH);
    let rho_amb = c.ambient * rho_charge;
    g.init(|iz, ir| match region(iz, ir) {
        // Prim is (rho, u_z, u_r, p); x is the grid's r direction.
        1 => Prim::new(rho_sheet, s.1, s.0, p_cold),
        2 => Prim::new(rho_foam, 0.0, 0.0, p_cold * (v_sheet / v_charge).min(1.0)),
        4 => Prim::new(rho_layer, 0.0, 0.0, p_cold * (v_sheet / v_charge).min(1.0)),
        3 => Prim::new(rho_charge, drift * s.1, drift * s.0, p_hot),
        _ => Prim::new(rho_amb, 0.0, 0.0, p_cold * 1.0e-6),
    });

    // The bodies must start cold (the collision makes the heat); the point fireball carries
    // exactly the collision's heat.
    let heat = initial_heat(&g, c.gamma);
    if point {
        assert!((heat / (0.5 * c.k / (1.0 + c.k)) - 1.0).abs() < 0.02);
    } else {
        assert!(
            heat < 0.02 * 0.5,
            "the sheet and charge must start cold against the sheet's kinetic energy of 1/2"
        );
    }
    let on_plate = |ir: usize| {
        let x = (ir as f64 + 0.5) * c.dx;
        x >= plate.0 && x <= plate.1
    };
    let force = |g: &Grid2D| -> (f64, f64, f64, f64) {
        let (mut f, mut m, mut pk, mut px) = (0.0, 0.0, 0.0_f64, 0.0);
        for ir in (0..nr).filter(|&ir| on_plate(ir)) {
            let p = g.prim(0, ir).p;
            let x = (ir as f64 + 0.5) * c.dx;
            f += p * c.dx;
            m += p * c.dx * (x - plate.0);
            if p > pk {
                pk = p;
                px = x - plate.0;
            }
        }
        (f, m, pk, px)
    };
    let (mut impulse, mut moment) = (0.0, 0.0);
    let (mut f_old, mut m_old, _, _) = force(&g);
    let (mut peak, mut peak_x, mut f_peak) = (0.0_f64, 0.0, 0.0_f64);
    let (mut steps, mut quiet) = (0usize, 0usize);
    let max_steps = 400 * nz.max(nr) + 50_000;
    while steps < max_steps {
        let dt = g.stable_dt();
        g.step(dt);
        let (f, m, pk, px) = force(&g);
        impulse += 0.5 * dt * (f_old + f);
        moment += 0.5 * dt * (m_old + m);
        (f_old, m_old) = (f, m);
        if pk > peak {
            peak = pk;
            peak_x = px;
        }
        f_peak = f_peak.max(f);
        quiet = if f_peak > 0.0 && f < 1.0e-3 * f_peak {
            quiet + 1
        } else {
            0
        };
        steps += 1;
        if quiet >= 40 {
            break;
        }
    }
    let centroid = if impulse > 0.0 { moment / impulse } else { 0.0 };
    CantedResult {
        impulse,
        ceiling: c.ceiling(),
        centroid_x: centroid,
        centroid_offset: centroid - c.charge_x,
        peak_pressure: peak,
        peak_x,
        ambient_mass_share: rho_amb * (v_all - v_sheet - v_charge) / (1.0 + c.k),
        layer_density_ratio: rho_layer / rho_foam,
        layer_mass_share: layer_share,
        steps,
    }
}

/// Internal energy on the grid per unit depth, `sum p/(gamma - 1) dV` (planar cells).
fn initial_heat(g: &Grid2D, gamma: f64) -> f64 {
    let cell = g.dr() * g.dz();
    (0..g.nz())
        .flat_map(|iz| (0..g.nr()).map(move |ir| (iz, ir)))
        .map(|(iz, ir)| g.prim(iz, ir).p / (gamma - 1.0) * cell)
        .sum()
}

#[cfg(test)]
mod tests {
    use super::{CantedStrip, run_canted_strip};

    fn small(theta_deg: f64, k: f64) -> CantedStrip {
        CantedStrip {
            theta_deg,
            k,
            standoff: 1.0,
            charge_x: 2.0,
            plate_length: 6.0,
            width: 0.8,
            charge_length: 0.3,
            dx: 0.1,
            gamma: 1.4,
            ambient: 1.0e-9,
            point_radius: None,
            shape: super::ChargeShape::Block,
            layer: None,
        }
    }

    /// Conservation bounds the canted plate: a frictionless flat plate keeps the tangential
    /// momentum out of reach and can at best spend all the energy normally, so the normal
    /// impulse stays under `cos theta + sqrt(1 + k - sin^2 theta)` -- both head-on (the overtake
    /// ceiling) and canted -- and is positive.
    #[test]
    fn canted_impulse_stays_under_its_ceiling() {
        for (theta, point) in [(0.0, None), (101.0, None), (101.0, Some(0.4))] {
            let c = CantedStrip {
                point_radius: point,
                ..small(theta, 10.0)
            };
            let r = run_canted_strip(&c);
            assert!(r.impulse > 0.0, "theta {theta}: impulse {}", r.impulse);
            assert!(
                r.impulse < r.ceiling,
                "theta {theta}: impulse {} over ceiling {}",
                r.impulse,
                r.ceiling
            );
            assert!(r.ambient_mass_share < 1.0e-3);
        }
    }

    /// A plate-side layer holds its share of `k` at close to the requested density ratio, the
    /// half-disk lying on its face is struck by the whole sheet, and both stay under the ceiling.
    #[test]
    fn layered_and_half_charges_stay_under_their_ceiling() {
        use super::{ChargeShape, PlateLayer};
        for (shape, layer) in [
            (ChargeShape::HalfDomeUp, None),
            (
                ChargeShape::Block,
                Some(PlateLayer {
                    mass_share: 0.5,
                    density_ratio: 5.0,
                }),
            ),
        ] {
            let c = CantedStrip {
                shape,
                layer,
                charge_length: 0.8,
                ..small(101.0, 10.0)
            };
            let r = run_canted_strip(&c);
            assert!(r.impulse > 0.0 && r.impulse < r.ceiling);
            if layer.is_some() {
                assert!((r.layer_density_ratio / 5.0 - 1.0).abs() < 0.35);
                assert!((r.layer_mass_share - 0.5).abs() < 1.0e-12);
            }
        }
    }
}
