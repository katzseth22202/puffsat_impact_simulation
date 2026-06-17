//! Rung-0 plumbing smoke test (the Rust middle of the boundary round-trip).
//!
//! Reads the toy JSON table written by `puffsat.smoke`, bilinearly interpolates one point,
//! and appends a JSONL result row. Proves: cargo toolchain, the workspace build, `serde_json`
//! both directions, and the cross-language JSON/JSONL schema contract (ADR-0019).

use std::fs;
use std::io::Write as _;
use std::path::Path;

use serde::{Deserialize, Serialize};

const TABLE_PATH: &str = "data/tables/smoke.json";
const RESULT_PATH: &str = "data/results/smoke.jsonl";

/// A gridded `(rho, T) -> field` table; `p` is flattened row-major over `(rho, T)`.
#[derive(Debug, Deserialize)]
struct Table {
    shape: [usize; 2],
    rho_grid: Vec<f64>,
    #[serde(rename = "T_grid")]
    t_grid: Vec<f64>,
    p: Vec<f64>,
}

/// One interpolated query, the JSONL result schema.
#[derive(Debug, Serialize)]
struct Record {
    rho: f64,
    #[serde(rename = "T")]
    t: f64,
    p_interp: f64,
}

/// Locate the lower bracketing index `i` and fractional weight in `[0, 1]` for `x` on `grid`.
/// Clamps to the grid ends. Assumes `grid` is sorted ascending with at least two points.
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

/// Bilinear interpolation of the `p` field at `(rho, t)`.
fn bilinear(table: &Table, rho: f64, t: f64) -> f64 {
    let n_t = table.shape[1];
    let (i, fr) = locate(&table.rho_grid, rho);
    let (j, ft) = locate(&table.t_grid, t);
    let at = |a: usize, b: usize| table.p[a * n_t + b];
    let lo = at(i, j) * (1.0 - ft) + at(i, j + 1) * ft;
    let hi = at(i + 1, j) * (1.0 - ft) + at(i + 1, j + 1) * ft;
    lo * (1.0 - fr) + hi * fr
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let raw = fs::read_to_string(TABLE_PATH)?;
    let table: Table = serde_json::from_str(&raw)?;

    let (rho, t) = (0.5, 0.5);
    let record = Record {
        rho,
        t,
        p_interp: bilinear(&table, rho, t),
    };
    let line = serde_json::to_string(&record)?;

    if let Some(parent) = Path::new(RESULT_PATH).parent() {
        fs::create_dir_all(parent)?;
    }
    let mut out = fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(RESULT_PATH)?;
    writeln!(out, "{line}")?;

    println!(
        "rust: p({rho}, {t}) = {} -> appended to {RESULT_PATH}",
        record.p_interp
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{Table, bilinear};

    fn toy() -> Table {
        Table {
            shape: [2, 2],
            rho_grid: vec![0.0, 1.0],
            t_grid: vec![0.0, 1.0],
            p: vec![0.0, 10.0, 20.0, 30.0],
        }
    }

    #[test]
    fn midpoint_is_corner_average() {
        assert!((bilinear(&toy(), 0.5, 0.5) - 15.0).abs() < 1e-12);
    }

    #[test]
    fn corners_are_exact() {
        assert!((bilinear(&toy(), 0.0, 0.0) - 0.0).abs() < 1e-12);
        assert!((bilinear(&toy(), 1.0, 1.0) - 30.0).abs() < 1e-12);
    }
}
