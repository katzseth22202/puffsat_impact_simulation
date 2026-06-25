//! Shared loader for the unified EOS/opacity table (ADR-0007).
//!
//! The Python cold path generates **one** gridded object mapping `(ρ, T)` to the material
//! properties the 1D rad-hydro kernel needs:
//!
//! ```text
//! (ρ, T) → (p, e, c_s, κ_Rosseland, κ_Planck)
//! ```
//!
//! serialized as JSON (ADR-0019). This crate reads the whole table into RAM once and
//! interpolates in-memory (ADR-0002) — no per-cell disk I/O in the hot loop.
//!
//! # Interpolation: bilinear in `(log ρ, log T)` with log fields
//!
//! A real water table spans many decades in both axes (ρ ~ 1e-1 … and T from ~300 K to ~43 kK)
//! and in every field, so a linear interpolant on the raw values would carry large error between
//! widely-spaced nodes. Interpolating **bilinearly in `(log ρ, log T)` on the log of each field**
//! instead is near-optimal across that dynamic range, and is **exact for any power-law EOS**
//! `field = K · ρ^a · T^b` (which becomes a plane in log-log space) — including the ideal gas
//! `p = (γ−1) ρ e` with `e = c_v T`. That exactness is what lets rung B regress the table path
//! against rung A's analytic ideal-gas results to machine precision (B2).
//!
//! All five fields are positive-definite, so the log is always defined; the loader validates
//! this on read.

use std::path::Path;

use serde::Deserialize;

/// An error reading or validating a [`Table`].
#[derive(Debug)]
pub enum TableError {
    /// The file could not be read.
    Io(std::io::Error),
    /// The JSON could not be parsed into the table schema.
    Parse(serde_json::Error),
    /// The table failed a structural/physical invariant (shape, monotonicity, positivity).
    Invalid(String),
}

impl std::fmt::Display for TableError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(e) => write!(f, "reading EOS/opacity table: {e}"),
            Self::Parse(e) => write!(f, "parsing EOS/opacity table JSON: {e}"),
            Self::Invalid(m) => write!(f, "invalid EOS/opacity table: {m}"),
        }
    }
}

impl std::error::Error for TableError {}

impl From<std::io::Error> for TableError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}

impl From<serde_json::Error> for TableError {
    fn from(e: serde_json::Error) -> Self {
        Self::Parse(e)
    }
}

/// The on-disk JSON schema (ADR-0007/0019): grids, a checksum `shape`, the flattened field
/// arrays nested under `fields`, and a free-form `provenance` object (kept for humans; not read
/// into the hot path).
#[derive(Debug, Deserialize)]
struct RawTable {
    rho_grid: Vec<f64>,
    #[serde(rename = "T_grid")]
    t_grid: Vec<f64>,
    shape: [usize; 2],
    fields: RawFields,
    #[serde(default)]
    #[allow(dead_code)]
    provenance: serde_json::Value,
}

/// The flattened field arrays, each row-major over `(ρ, T)` with length `n_ρ · n_T`.
#[derive(Debug, Deserialize)]
struct RawFields {
    p: Vec<f64>,
    e: Vec<f64>,
    c_s: Vec<f64>,
    kappa_rosseland: Vec<f64>,
    kappa_planck: Vec<f64>,
    /// Optional condensed mass fraction `∈ [0, 1]` (Rung C low-v tables). Absent in the high-v table.
    #[serde(default)]
    liquid_frac: Option<Vec<f64>>,
    /// Optional gas thermal conductivity `k_gas > 0` (B-flux low-v tables). Absent in the high-v table.
    #[serde(default)]
    k_gas: Option<Vec<f64>>,
}

/// The five material properties at one `(ρ, T)` point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Material {
    /// Pressure `p`.
    pub p: f64,
    /// Specific internal energy `e`.
    pub e: f64,
    /// Adiabatic sound speed `c_s`.
    pub c_s: f64,
    /// Rosseland-mean opacity `κ_R` (diffusion coefficient; ADR-0006).
    pub kappa_rosseland: f64,
    /// Planck-mean opacity `κ_P` (emission/absorption source; ADR-0006).
    pub kappa_planck: f64,
}

/// A unified EOS/opacity table, loaded and ready for in-memory interpolation.
///
/// Holds the natural-log of both axes and of every field, so a query is a bilinear blend in
/// log-log space followed by one `exp`.
#[derive(Debug, Clone)]
pub struct Table {
    rho_grid: Vec<f64>,
    t_grid: Vec<f64>,
    log_rho: Vec<f64>,
    log_t: Vec<f64>,
    log_p: Vec<f64>,
    log_e: Vec<f64>,
    log_cs: Vec<f64>,
    log_kr: Vec<f64>,
    log_kp: Vec<f64>,
    /// Raw (not log) condensed mass fraction, interpolated **linearly** (it is `[0, 1]`-valued and
    /// legitimately `0`). `None` for tables that omit it — then [`Table::liquid_fraction`] is `0`.
    liquid_frac: Option<Vec<f64>>,
    /// Natural-log of the gas thermal conductivity `k_gas > 0` (the B-flux conduction operator's
    /// gas-side property; ADR-0005). On the positive log-interp path like the opacities. `None` for
    /// tables that omit it — then [`Table::k_gas`] is `None` and the gas gets no conduction.
    log_k_gas: Option<Vec<f64>>,
}

impl Table {
    /// Parse and validate a table from a JSON string.
    ///
    /// # Errors
    /// Returns [`TableError::Parse`] on malformed JSON and [`TableError::Invalid`] if the table
    /// violates a structural or physical invariant (shape mismatch, non-ascending grid, or a
    /// non-positive grid value or field value — all required for the log-log interpolation).
    pub fn from_json(json: &str) -> Result<Self, TableError> {
        let raw: RawTable = serde_json::from_str(json)?;
        Self::from_raw(raw)
    }

    /// Load and validate a table from a JSON file.
    ///
    /// # Errors
    /// As [`Self::from_json`], plus [`TableError::Io`] if the file cannot be read.
    pub fn load(path: impl AsRef<Path>) -> Result<Self, TableError> {
        let raw = std::fs::read_to_string(path)?;
        Self::from_json(&raw)
    }

    fn from_raw(raw: RawTable) -> Result<Self, TableError> {
        let n_rho = raw.rho_grid.len();
        let n_t = raw.t_grid.len();
        let invalid = |m: String| Err(TableError::Invalid(m));

        if n_rho < 2 || n_t < 2 {
            return invalid(format!(
                "need at least 2 points per axis, got {n_rho}×{n_t}"
            ));
        }
        if raw.shape != [n_rho, n_t] {
            return invalid(format!(
                "shape {:?} disagrees with grids ({n_rho}×{n_t})",
                raw.shape
            ));
        }
        for (name, grid) in [("rho_grid", &raw.rho_grid), ("T_grid", &raw.t_grid)] {
            if grid[0] <= 0.0 {
                return invalid(format!("{name} must be positive for log interpolation"));
            }
            if !grid.windows(2).all(|w| w[0] < w[1]) {
                return invalid(format!("{name} must be strictly ascending"));
            }
        }

        let n = n_rho * n_t;
        let log_field = |name: &str, v: &[f64]| -> Result<Vec<f64>, TableError> {
            if v.len() != n {
                return Err(TableError::Invalid(format!(
                    "field `{name}` has {} entries, expected {n} (= {n_rho}×{n_t})",
                    v.len()
                )));
            }
            // Catch ≤ 0 *and* NaN (NaN is not `Greater` than 0) without a negated float compare.
            if let Some(&bad) = v
                .iter()
                .find(|&&x| x.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater))
            {
                return Err(TableError::Invalid(format!(
                    "field `{name}` has a non-positive value {bad} (log interpolation requires > 0)"
                )));
            }
            Ok(v.iter().map(|x| x.ln()).collect())
        };

        // `liquid_frac` (optional, Rung C) is interpolated linearly, so it is *not* logged; it must
        // be fully populated and bounded in `[0, 1]` (this also rejects NaN, which fails `contains`).
        let liquid_frac = match raw.fields.liquid_frac {
            Some(v) => {
                if v.len() != n {
                    return invalid(format!(
                        "field `liquid_frac` has {} entries, expected {n}",
                        v.len()
                    ));
                }
                if let Some(&bad) = v.iter().find(|&&x| !(0.0..=1.0).contains(&x)) {
                    return invalid(format!("field `liquid_frac` value {bad} outside [0, 1]"));
                }
                Some(v)
            }
            None => None,
        };

        // `k_gas` (optional, B-flux) is strictly positive and spans decades, so it rides the same
        // positive log-interp path as the opacities (`log_field` validates length and `> 0`).
        let log_k_gas = raw
            .fields
            .k_gas
            .as_deref()
            .map(|v| log_field("k_gas", v))
            .transpose()?;

        Ok(Self {
            log_rho: raw.rho_grid.iter().map(|x| x.ln()).collect(),
            log_t: raw.t_grid.iter().map(|x| x.ln()).collect(),
            log_p: log_field("p", &raw.fields.p)?,
            log_e: log_field("e", &raw.fields.e)?,
            log_cs: log_field("c_s", &raw.fields.c_s)?,
            log_kr: log_field("kappa_rosseland", &raw.fields.kappa_rosseland)?,
            log_kp: log_field("kappa_planck", &raw.fields.kappa_planck)?,
            liquid_frac,
            log_k_gas,
            rho_grid: raw.rho_grid,
            t_grid: raw.t_grid,
        })
    }

    /// The density grid (ascending).
    #[must_use]
    pub fn rho_grid(&self) -> &[f64] {
        &self.rho_grid
    }

    /// The temperature grid (ascending).
    #[must_use]
    pub fn t_grid(&self) -> &[f64] {
        &self.t_grid
    }

    /// Pressure `p(ρ, T)`.
    #[must_use]
    pub fn pressure(&self, rho: f64, t: f64) -> f64 {
        self.interp(&self.log_p, rho, t)
    }

    /// Specific internal energy `e(ρ, T)`.
    #[must_use]
    pub fn energy(&self, rho: f64, t: f64) -> f64 {
        self.interp(&self.log_e, rho, t)
    }

    /// Adiabatic sound speed `c_s(ρ, T)`.
    #[must_use]
    pub fn sound_speed(&self, rho: f64, t: f64) -> f64 {
        self.interp(&self.log_cs, rho, t)
    }

    /// Rosseland-mean opacity `κ_R(ρ, T)`.
    #[must_use]
    pub fn kappa_rosseland(&self, rho: f64, t: f64) -> f64 {
        self.interp(&self.log_kr, rho, t)
    }

    /// Planck-mean opacity `κ_P(ρ, T)`.
    #[must_use]
    pub fn kappa_planck(&self, rho: f64, t: f64) -> f64 {
        self.interp(&self.log_kp, rho, t)
    }

    /// A copy of this table with **both** opacities (Rosseland + Planck) scaled by `factor` — the
    /// opacity-scale knob (B5d-3 sensitivity scan; Rung E τ-bracket). Multiplying `κ` by `factor`
    /// shifts its natural-log field by `ln factor`, so the log-log interpolation is preserved exactly
    /// and only the optical depth moves. `factor = 1` is the identity. The EOS fields (`p`, `e`,
    /// `c_s`) and any `liquid_frac` / `k_gas` are untouched, isolating the radiative-transport regime.
    ///
    /// # Panics
    /// Panics unless `factor > 0` (`κ` must stay positive for the log interpolation).
    #[must_use]
    pub fn with_opacity_scale(&self, factor: f64) -> Self {
        assert!(factor > 0.0, "opacity scale must be positive");
        let shift = factor.ln();
        let mut scaled = self.clone();
        for v in &mut scaled.log_kr {
            *v += shift;
        }
        for v in &mut scaled.log_kp {
            *v += shift;
        }
        scaled
    }

    /// Condensed mass fraction `liquid_frac(ρ, T) ∈ [0, 1]` — the Rung C low-v two-phase tables carry
    /// this for the wall-sticking condensation sink; tables without it (e.g. the high-v table) return
    /// `0`. Interpolated **linearly** (it is `[0, 1]`-valued and legitimately `0`, unlike the log
    /// fields), on the same `(log ρ, log T)` axes.
    #[must_use]
    pub fn liquid_fraction(&self, rho: f64, t: f64) -> f64 {
        match &self.liquid_frac {
            Some(field) => self.interp_linear(field, rho, t),
            None => 0.0,
        }
    }

    /// Gas thermal conductivity `k_gas(ρ, T) > 0` [W/m/K] — the B-flux conduction operator's gas-side
    /// property (ADR-0005). `Some` for the low-v B-flux tables that carry it (log-interpolated like
    /// the opacities), `None` for tables that omit it (e.g. the high-v table), where the caller gives
    /// the gas no conduction.
    #[must_use]
    pub fn k_gas(&self, rho: f64, t: f64) -> Option<f64> {
        self.log_k_gas
            .as_ref()
            .map(|field| self.interp(field, rho, t))
    }

    /// Specific heat capacity `c_v = ∂e/∂T` at fixed `ρ`, by centered finite difference on the
    /// tabulated `e(ρ, T)`. The radiation step (B5) multiplies this by `ρ` to get the volumetric
    /// heat capacity `ρ c_v` it needs. The step `δ` is a small fraction of `T`, matching the EOS's
    /// `dp_de` central difference.
    #[must_use]
    pub fn cv(&self, rho: f64, t: f64) -> f64 {
        let dt = 1e-4 * t.max(1.0);
        (self.energy(rho, t + dt) - self.energy(rho, t - dt)) / (2.0 * dt)
    }

    /// All five properties at `(ρ, T)` in one query.
    #[must_use]
    pub fn material(&self, rho: f64, t: f64) -> Material {
        Material {
            p: self.pressure(rho, t),
            e: self.energy(rho, t),
            c_s: self.sound_speed(rho, t),
            kappa_rosseland: self.kappa_rosseland(rho, t),
            kappa_planck: self.kappa_planck(rho, t),
        }
    }

    /// Bilinear interpolation of one log-field at `(ρ, T)`, in `(log ρ, log T)` space, returning
    /// the de-logged value. Queries outside the grid clamp to the nearest edge.
    fn interp(&self, log_field: &[f64], rho: f64, t: f64) -> f64 {
        let (i, fr) = locate(&self.log_rho, rho.ln());
        let (j, ft) = locate(&self.log_t, t.ln());
        let n_t = self.t_grid.len();
        let at = |a: usize, b: usize| log_field[a * n_t + b];
        let lo = at(i, j).mul_add(1.0 - ft, at(i, j + 1) * ft);
        let hi = at(i + 1, j).mul_add(1.0 - ft, at(i + 1, j + 1) * ft);
        lo.mul_add(1.0 - fr, hi * fr).exp()
    }

    /// Bilinear interpolation of a **raw** (non-log) field at `(ρ, T)`, on the `(log ρ, log T)` axes
    /// — for fields like `liquid_frac` that are `[0, 1]`-valued and may be `0`. No `exp`.
    fn interp_linear(&self, field: &[f64], rho: f64, t: f64) -> f64 {
        let (i, fr) = locate(&self.log_rho, rho.ln());
        let (j, ft) = locate(&self.log_t, t.ln());
        let n_t = self.t_grid.len();
        let at = |a: usize, b: usize| field[a * n_t + b];
        let lo = at(i, j).mul_add(1.0 - ft, at(i, j + 1) * ft);
        let hi = at(i + 1, j).mul_add(1.0 - ft, at(i + 1, j + 1) * ft);
        lo.mul_add(1.0 - fr, hi * fr)
    }
}

/// Locate the lower bracketing index `i` and the fractional weight `∈ [0, 1]` of `x` on a sorted
/// ascending `grid` (at least two points). Clamps to the grid ends.
fn locate(grid: &[f64], x: f64) -> (usize, f64) {
    let n = grid.len();
    if x <= grid[0] {
        return (0, 0.0);
    }
    if x >= grid[n - 1] {
        return (n - 2, 1.0);
    }
    let i = grid.partition_point(|&g| g <= x) - 1;
    let frac = (x - grid[i]) / (grid[i + 1] - grid[i]);
    (i, frac)
}

#[cfg(test)]
mod tests {
    use super::{Table, TableError};
    use approx::assert_relative_eq;

    /// Build a JSON table on log-spaced grids whose fields are exact power laws of `(ρ, T)`.
    /// `p = Kp·ρ^a·T^b`, etc. Log-log bilinear interpolation must reproduce these exactly.
    fn power_law_json(n_rho: usize, n_t: usize) -> String {
        // Geometric (log-spaced) grids.
        let rho_grid: Vec<f64> = (0..n_rho)
            .map(|i| 0.1 * 10f64.powf(i as f64 / (n_rho - 1) as f64 * 3.0)) // 0.1 … 100
            .collect();
        let t_grid: Vec<f64> = (0..n_t)
            .map(|j| 300.0 * 10f64.powf(j as f64 / (n_t - 1) as f64 * 2.0)) // 300 … 30000
            .collect();
        // field(ρ,T) = k · ρ^a · T^b, row-major over (ρ, T).
        let field = |k: f64, a: f64, b: f64| -> Vec<f64> {
            let mut v = Vec::with_capacity(n_rho * n_t);
            for &r in &rho_grid {
                for &t in &t_grid {
                    v.push(k * r.powf(a) * t.powf(b));
                }
            }
            v
        };
        let table = serde_json::json!({
            "rho_grid": rho_grid,
            "T_grid": t_grid,
            "shape": [n_rho, n_t],
            "fields": {
                "p": field(2.0, 1.0, 1.0),          // p ∝ ρ T  (ideal-gas-like)
                "e": field(3.0, 0.0, 1.0),          // e ∝ T
                "c_s": field(5.0, 0.0, 0.5),        // c_s ∝ √T
                "kappa_rosseland": field(7.0, 2.0, -3.5),
                "kappa_planck": field(11.0, 1.5, -2.0),
            },
            "provenance": {"source": "power-law unit-test table"},
        });
        table.to_string()
    }

    /// The opacity-scale knob multiplies both opacities by `factor` at every `(ρ, T)`, leaving the
    /// EOS fields (`p`, `e`, `c_s`) untouched, and `factor = 1` is the identity.
    #[test]
    fn opacity_scale_multiplies_both_opacities() {
        let table = Table::from_json(&power_law_json(6, 5)).unwrap();
        let scaled = table.with_opacity_scale(10.0);
        let identity = table.with_opacity_scale(1.0);
        for &(rho, t) in &[(0.37, 1234.0), (3.3, 555.0), (42.0, 9001.0)] {
            assert_relative_eq!(
                scaled.kappa_rosseland(rho, t),
                10.0 * table.kappa_rosseland(rho, t),
                max_relative = 1e-12
            );
            assert_relative_eq!(
                scaled.kappa_planck(rho, t),
                10.0 * table.kappa_planck(rho, t),
                max_relative = 1e-12
            );
            // EOS fields are untouched by the opacity scale.
            assert_relative_eq!(
                scaled.pressure(rho, t),
                table.pressure(rho, t),
                max_relative = 1e-12
            );
            assert_relative_eq!(
                scaled.energy(rho, t),
                table.energy(rho, t),
                max_relative = 1e-12
            );
            // factor = 1 reproduces the original opacity exactly.
            assert_relative_eq!(
                identity.kappa_rosseland(rho, t),
                table.kappa_rosseland(rho, t),
                max_relative = 1e-12
            );
        }
    }

    /// log-log bilinear interpolation is exact for power-law fields, at arbitrary off-grid points.
    #[test]
    fn power_law_interpolation_is_exact() {
        let table = Table::from_json(&power_law_json(6, 5)).unwrap();
        for &(rho, t) in &[(0.37, 1234.0), (3.3, 555.0), (42.0, 9001.0), (0.5, 4000.0)] {
            assert_relative_eq!(table.pressure(rho, t), 2.0 * rho * t, max_relative = 1e-12);
            assert_relative_eq!(table.energy(rho, t), 3.0 * t, max_relative = 1e-12);
            assert_relative_eq!(
                table.sound_speed(rho, t),
                5.0 * t.sqrt(),
                max_relative = 1e-12
            );
            assert_relative_eq!(
                table.kappa_rosseland(rho, t),
                7.0 * rho.powf(2.0) * t.powf(-3.5),
                max_relative = 1e-12
            );
            assert_relative_eq!(
                table.kappa_planck(rho, t),
                11.0 * rho.powf(1.5) * t.powf(-2.0),
                max_relative = 1e-12
            );
        }
    }

    /// A JSON table with a **curved** energy field `e = ke·T²` (other fields power-law). The log-log
    /// interpolation is still exact (e is a power law of T), so `∂e/∂T = 2·ke·T` is known — and is
    /// distinguishable from `e/T = ke·T`, so a `cv` that confused the two would be caught.
    fn quadratic_energy_json(n_rho: usize, n_t: usize, ke: f64) -> String {
        let rho_grid: Vec<f64> = (0..n_rho)
            .map(|i| 0.1 * 10f64.powf(i as f64 / (n_rho - 1) as f64 * 3.0))
            .collect();
        let t_grid: Vec<f64> = (0..n_t)
            .map(|j| 300.0 * 10f64.powf(j as f64 / (n_t - 1) as f64 * 2.0))
            .collect();
        let field = |k: f64, a: f64, b: f64| -> Vec<f64> {
            let mut v = Vec::with_capacity(n_rho * n_t);
            for &r in &rho_grid {
                for &t in &t_grid {
                    v.push(k * r.powf(a) * t.powf(b));
                }
            }
            v
        };
        serde_json::json!({
            "rho_grid": rho_grid,
            "T_grid": t_grid,
            "shape": [n_rho, n_t],
            "fields": {
                "p": field(2.0, 1.0, 1.0),
                "e": field(ke, 0.0, 2.0),       // e = ke·T²  ⇒  ∂e/∂T = 2·ke·T
                "c_s": field(5.0, 0.0, 0.5),
                "kappa_rosseland": field(7.0, 2.0, -3.5),
                "kappa_planck": field(11.0, 1.5, -2.0),
            },
            "provenance": {"source": "quadratic-energy unit-test table"},
        })
        .to_string()
    }

    /// `Table::cv` returns `∂e/∂T` (not `e/T`): the centered FD recovers `2·ke·T` for `e = ke·T²`.
    #[test]
    fn cv_matches_de_dt() {
        let ke = 0.5;
        let table = Table::from_json(&quadratic_energy_json(6, 5, ke)).unwrap();
        for &(rho, t) in &[(0.37, 1234.0), (3.3, 555.0), (42.0, 9001.0)] {
            assert_relative_eq!(table.cv(rho, t), 2.0 * ke * t, max_relative = 1e-6);
        }
    }

    /// Querying exactly on a grid node returns the stored value (the interpolation is consistent).
    #[test]
    fn grid_nodes_recover_stored_values() {
        let table = Table::from_json(&power_law_json(4, 4)).unwrap();
        let rho = table.rho_grid()[2];
        let t = table.t_grid()[1];
        assert_relative_eq!(table.pressure(rho, t), 2.0 * rho * t, max_relative = 1e-13);
    }

    /// Queries past the grid edge clamp to the boundary value (no extrapolation blow-up).
    #[test]
    fn out_of_range_clamps_to_edge() {
        let table = Table::from_json(&power_law_json(4, 4)).unwrap();
        let rmin = table.rho_grid()[0];
        let tmin = table.t_grid()[0];
        // Far below both grid mins clamps to the corner node value.
        assert_relative_eq!(
            table.pressure(rmin * 1e-6, tmin * 1e-6),
            2.0 * rmin * tmin,
            max_relative = 1e-13
        );
    }

    /// `material` agrees with the individual accessors.
    #[test]
    fn material_matches_accessors() {
        let table = Table::from_json(&power_law_json(5, 5)).unwrap();
        let (rho, t) = (1.7, 2500.0);
        let m = table.material(rho, t);
        // `material` calls the same accessors, so the results are bit-identical.
        assert_eq!(m.p.to_bits(), table.pressure(rho, t).to_bits());
        assert_eq!(
            m.kappa_planck.to_bits(),
            table.kappa_planck(rho, t).to_bits()
        );
    }

    #[test]
    fn rejects_shape_mismatch() {
        let json = r#"{
            "rho_grid": [1.0, 2.0],
            "T_grid": [1.0, 2.0],
            "shape": [2, 3],
            "fields": {"p":[1,1,1,1],"e":[1,1,1,1],"c_s":[1,1,1,1],
                       "kappa_rosseland":[1,1,1,1],"kappa_planck":[1,1,1,1]}
        }"#;
        assert!(matches!(
            Table::from_json(json),
            Err(TableError::Invalid(_))
        ));
    }

    #[test]
    fn rejects_non_positive_field() {
        let json = r#"{
            "rho_grid": [1.0, 2.0],
            "T_grid": [1.0, 2.0],
            "shape": [2, 2],
            "fields": {"p":[1,1,1,0],"e":[1,1,1,1],"c_s":[1,1,1,1],
                       "kappa_rosseland":[1,1,1,1],"kappa_planck":[1,1,1,1]}
        }"#;
        assert!(matches!(
            Table::from_json(json),
            Err(TableError::Invalid(_))
        ));
    }

    /// `liquid_frac` (optional, Rung C): interpolated linearly, `0` allowed (unlike the log fields),
    /// corners recovered exactly and the geometric midpoint is the bilinear corner average.
    #[test]
    fn liquid_fraction_linear_interp_and_zero_allowed() {
        let json = r#"{
            "rho_grid": [1.0, 2.0],
            "T_grid": [1.0, 2.0],
            "shape": [2, 2],
            "fields": {"p":[1,1,1,1],"e":[1,2,1,2],"c_s":[1,1,1,1],
                       "kappa_rosseland":[1,1,1,1],"kappa_planck":[1,1,1,1],
                       "liquid_frac":[0.0, 0.2, 0.4, 0.6]}
        }"#;
        let table = Table::from_json(json).unwrap();
        // Corner with value 0.0 loads and is returned exactly (not a positive-log field).
        assert_relative_eq!(table.liquid_fraction(1.0, 1.0), 0.0, epsilon = 1e-12);
        assert_relative_eq!(table.liquid_fraction(2.0, 2.0), 0.6, max_relative = 1e-12);
        // Geometric midpoint (√2, √2) → fr = ft = ½ → average of the four corners.
        let m = 2f64.sqrt();
        assert_relative_eq!(table.liquid_fraction(m, m), 0.3, max_relative = 1e-12);
    }

    /// A table that omits `liquid_frac` (e.g. the high-v table) loads and reports `0` everywhere.
    #[allow(clippy::float_cmp)] // exact: the absent-field branch returns the literal 0.0
    #[test]
    fn liquid_fraction_absent_is_zero() {
        let table = Table::from_json(&power_law_json(4, 4)).unwrap();
        assert_eq!(table.liquid_fraction(1.7, 2500.0), 0.0);
    }

    #[test]
    fn rejects_liquid_frac_out_of_range() {
        let json = r#"{
            "rho_grid": [1.0, 2.0],
            "T_grid": [1.0, 2.0],
            "shape": [2, 2],
            "fields": {"p":[1,1,1,1],"e":[1,1,1,1],"c_s":[1,1,1,1],
                       "kappa_rosseland":[1,1,1,1],"kappa_planck":[1,1,1,1],
                       "liquid_frac":[0.0, 0.5, 1.0, 1.5]}
        }"#;
        assert!(matches!(
            Table::from_json(json),
            Err(TableError::Invalid(_))
        ));
    }

    /// `k_gas` (optional, B-flux): strictly positive, log-interpolated like the opacities — a
    /// power-law field `k_gas = ρ·T` is reproduced exactly; a table that omits it returns `None`
    /// (no gas conduction).
    #[test]
    fn k_gas_log_interp_and_absent_is_none() {
        let json = r#"{
            "rho_grid": [1.0, 2.0],
            "T_grid": [1.0, 2.0],
            "shape": [2, 2],
            "fields": {"p":[1,1,1,1],"e":[1,1,1,1],"c_s":[1,1,1,1],
                       "kappa_rosseland":[1,1,1,1],"kappa_planck":[1,1,1,1],
                       "k_gas":[1.0, 2.0, 2.0, 4.0]}
        }"#;
        let table = Table::from_json(json).unwrap();
        assert_relative_eq!(table.k_gas(1.0, 1.0).unwrap(), 1.0, max_relative = 1e-12);
        assert_relative_eq!(table.k_gas(2.0, 2.0).unwrap(), 4.0, max_relative = 1e-12);
        // Geometric midpoint (√2, √2): log-log interp of ρ·T is exact ⇒ √2·√2 = 2.
        let m = 2f64.sqrt();
        assert_relative_eq!(table.k_gas(m, m).unwrap(), 2.0, max_relative = 1e-12);
        // A table without the field (e.g. the high-v table) has no gas conduction.
        assert!(
            Table::from_json(&power_law_json(4, 4))
                .unwrap()
                .k_gas(1.7, 2500.0)
                .is_none()
        );
    }

    #[test]
    fn rejects_non_positive_k_gas() {
        let json = r#"{
            "rho_grid": [1.0, 2.0],
            "T_grid": [1.0, 2.0],
            "shape": [2, 2],
            "fields": {"p":[1,1,1,1],"e":[1,1,1,1],"c_s":[1,1,1,1],
                       "kappa_rosseland":[1,1,1,1],"kappa_planck":[1,1,1,1],
                       "k_gas":[1.0, 0.0, 2.0, 4.0]}
        }"#;
        assert!(matches!(
            Table::from_json(json),
            Err(TableError::Invalid(_))
        ));
    }

    #[test]
    fn rejects_non_ascending_grid() {
        let json = r#"{
            "rho_grid": [2.0, 1.0],
            "T_grid": [1.0, 2.0],
            "shape": [2, 2],
            "fields": {"p":[1,1,1,1],"e":[1,1,1,1],"c_s":[1,1,1,1],
                       "kappa_rosseland":[1,1,1,1],"kappa_planck":[1,1,1,1]}
        }"#;
        assert!(matches!(
            Table::from_json(json),
            Err(TableError::Invalid(_))
        ));
    }
}
