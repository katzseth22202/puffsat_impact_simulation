//! Equation of state for the Lagrangian kernel.
//!
//! ADR-0022's whole point is that a staggered artificial-viscosity scheme "swallows an arbitrary
//! tabulated EOS with a bare `p(ρ, e)` call." This module is that interface. The kernel evolves
//! the specific internal energy `e`, so the EOS must answer, given a cell's `(ρ, e)`:
//!
//! - the pressure `p(ρ, e)` (the force on the mesh),
//! - the adiabatic sound speed `c_s(ρ, e)` (the CFL signal speed),
//! - `∂p/∂e|_ρ` (the Newton slope of the implicit energy update), and
//! - the inverse `e(ρ, p)` (to seed `e` from an initial pressure).
//!
//! Two implementations:
//! - [`IdealGas`] — rung A's analytic `p = (γ−1) ρ e`.
//! - [`TableEos`] — rung B's tabulated equilibrium EOS (ADR-0007), wrapping a [`tables::Table`]
//!   indexed by `(ρ, T)`. Because the kernel carries `e`, every query first inverts the monotone
//!   `e(ρ, T)=e` (or `p(ρ, T)=p`) for `T` — by a `T`-grid binary search plus an analytic inverse
//!   of the log-log interpolant (see [`invert_field`]) — then reads off the tabulated property.

use tables::Table;

/// The EOS interface the kernel calls, per cell, with `e` the specific internal energy the
/// Lagrangian scheme conserves. Implementors are `Debug + Clone` so [`crate::kernel::Tube`] can
/// derive both.
pub trait Eos: core::fmt::Debug + Clone {
    /// Pressure `p(ρ, e)`.
    fn pressure(&self, rho: f64, e: f64) -> f64;
    /// Adiabatic sound speed `c_s(ρ, e)`.
    fn sound_speed(&self, rho: f64, e: f64) -> f64;
    /// Specific internal energy from a pressure, `e(ρ, p)` — seeds `e` from the initial `p`.
    fn energy_from_pressure(&self, rho: f64, p: f64) -> f64;
    /// `∂p/∂e` at fixed `ρ`: the slope of the Newton step in the implicit energy update. Only its
    /// accuracy near the root sets convergence *speed*; the root itself is fixed by `pressure`.
    fn dp_de(&self, rho: f64, e: f64) -> f64;
    /// Temperature `T(ρ, e)` — the radiation coupling (B5) needs it to set the Planck emission
    /// `aT⁴` and to index the opacity table. For the ideal gas this is `e` in its reduced units
    /// (`c_v = 1`); for the table it inverts the monotone `e(ρ, T)`.
    fn temperature(&self, rho: f64, e: f64) -> f64;
}

/// Ideal gas `p = (γ−1) ρ e` — rung A's EOS, the analytic baseline every table path regresses
/// against (B2).
#[derive(Debug, Clone, Copy)]
pub struct IdealGas {
    gamma: f64,
}

impl IdealGas {
    /// Construct an ideal gas with adiabatic index `γ`.
    #[must_use]
    pub fn new(gamma: f64) -> Self {
        Self { gamma }
    }

    /// The adiabatic index `γ`.
    #[must_use]
    pub fn gamma(&self) -> f64 {
        self.gamma
    }
}

impl Eos for IdealGas {
    fn pressure(&self, rho: f64, e: f64) -> f64 {
        (self.gamma - 1.0) * rho * e
    }

    fn sound_speed(&self, _rho: f64, e: f64) -> f64 {
        // c² = γ p / ρ = γ (γ−1) e — independent of ρ for an ideal gas.
        (self.gamma * (self.gamma - 1.0) * e).max(0.0).sqrt()
    }

    fn energy_from_pressure(&self, rho: f64, p: f64) -> f64 {
        p / ((self.gamma - 1.0) * rho)
    }

    fn dp_de(&self, rho: f64, _e: f64) -> f64 {
        (self.gamma - 1.0) * rho
    }

    fn temperature(&self, _rho: f64, e: f64) -> f64 {
        e // reduced units: e = c_v T with c_v = 1 (matches the rung-A / ideal-gas table convention)
    }
}

/// Tabulated equilibrium EOS (ADR-0007): a [`tables::Table`] of `(ρ, T) → (p, e, c_s, …)` queried
/// at the kernel's `(ρ, e)` by inverting for `T` on the monotone temperature axis.
#[derive(Debug, Clone)]
pub struct TableEos {
    table: Table,
}

impl TableEos {
    /// Wrap a loaded table as an EOS.
    #[must_use]
    pub fn new(table: Table) -> Self {
        Self { table }
    }

    /// Borrow the underlying table (opacity lookups, grids).
    #[must_use]
    pub fn table(&self) -> &Table {
        &self.table
    }

    /// Temperature such that `e(ρ, T) = e_target` (`e` is monotone increasing in `T` at fixed `ρ`,
    /// since heat capacity is positive).
    fn temperature_from_energy(&self, rho: f64, e_target: f64) -> f64 {
        invert_field(self.table.t_grid(), e_target, |t| self.table.energy(rho, t))
    }

    /// Temperature such that `p(ρ, T) = p_target` (`p` is monotone increasing in `T` at fixed `ρ`).
    fn temperature_from_pressure(&self, rho: f64, p_target: f64) -> f64 {
        invert_field(self.table.t_grid(), p_target, |t| {
            self.table.pressure(rho, t)
        })
    }

    /// Condensed mass fraction `liquid_frac(ρ, e) ∈ [0, 1]` at the kernel's `(ρ, e)`: invert the
    /// energy for `T`, then read the table's (linearly-interpolated) field — the Rung C wall-sticking
    /// sink reads this. `0` for tables without the field (e.g. the high-v table).
    #[must_use]
    pub fn liquid_fraction(&self, rho: f64, e: f64) -> f64 {
        let t = self.temperature_from_energy(rho, e);
        self.table.liquid_fraction(rho, t)
    }

    /// Gas thermal conductivity `k_gas(ρ, e) > 0` [W/m/K] at the kernel's `(ρ, e)`: invert the energy
    /// for `T`, then read the table's (log-interpolated) field — the B-flux conduction operator reads
    /// this (ADR-0005). `None` for tables without the field (e.g. the high-v table), where the gas
    /// gets no conduction.
    #[must_use]
    pub fn k_gas(&self, rho: f64, e: f64) -> Option<f64> {
        let t = self.temperature_from_energy(rho, e);
        self.table.k_gas(rho, t)
    }
}

impl Eos for TableEos {
    fn pressure(&self, rho: f64, e: f64) -> f64 {
        let t = self.temperature_from_energy(rho, e);
        self.table.pressure(rho, t)
    }

    fn sound_speed(&self, rho: f64, e: f64) -> f64 {
        let t = self.temperature_from_energy(rho, e);
        self.table.sound_speed(rho, t)
    }

    fn energy_from_pressure(&self, rho: f64, p: f64) -> f64 {
        let t = self.temperature_from_pressure(rho, p);
        self.table.energy(rho, t)
    }

    fn dp_de(&self, rho: f64, e: f64) -> f64 {
        // ∂p/∂e|_ρ = (∂p/∂T)/(∂e/∂T) at fixed ρ, by central difference about the inverted T.
        let t = self.temperature_from_energy(rho, e);
        let dt = 1e-4 * t.max(1.0);
        let dp = self.table.pressure(rho, t + dt) - self.table.pressure(rho, t - dt);
        let de = self.table.energy(rho, t + dt) - self.table.energy(rho, t - dt);
        dp / de
    }

    fn temperature(&self, rho: f64, e: f64) -> f64 {
        self.temperature_from_energy(rho, e)
    }
}

/// Invert a monotone-increasing tabulated field for the temperature where `f(T) = target`,
/// **consistently with the table's own log-log interpolation**. Binary-search the `T` grid for the
/// bracketing interval, then invert analytically inside it: at fixed `ρ` the bilinear-in-`(log ρ,
/// log T)` interpolant is log-log-linear in `T` within a cell, so this returns the exact inverse
/// of the forward interpolation (not merely a tolerance-converged root), in a few cheap
/// evaluations. Out-of-table targets clamp to the nearest grid edge, matching the table's clamp.
fn invert_field(t_grid: &[f64], target: f64, f: impl Fn(f64) -> f64) -> f64 {
    let n = t_grid.len();
    if target <= f(t_grid[0]) {
        return t_grid[0];
    }
    if target >= f(t_grid[n - 1]) {
        return t_grid[n - 1];
    }
    // Binary search for the interval [lo, hi] with f(t_grid[lo]) ≤ target < f(t_grid[hi]).
    let (mut lo, mut hi) = (0usize, n - 1);
    while hi - lo > 1 {
        let mid = usize::midpoint(lo, hi);
        if f(t_grid[mid]) <= target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let (t0, t1) = (t_grid[lo], t_grid[hi]);
    let (f0, f1) = (f(t0), f(t1));
    if f1 <= f0 {
        return t0; // degenerate (field flat over the interval); shouldn't happen for a real EOS.
    }
    // log-log-linear inverse: w = (ln target − ln f0)/(ln f1 − ln f0); ln T = ln t0 + w·(ln t1 − ln t0).
    let w = (target.ln() - f0.ln()) / (f1.ln() - f0.ln());
    t0.ln().mul_add(1.0 - w, t1.ln() * w).exp()
}

#[cfg(test)]
mod tests {
    use super::{Eos, IdealGas, TableEos};
    use approx::assert_relative_eq;
    use tables::Table;

    const GAMMA: f64 = 1.4;

    /// A table that *encodes the ideal gas* `p=(γ−1)ρe` with `e = T` (heat capacity `c_v = 1`):
    /// `e = T`, `p = (γ−1)ρT`, `c_s = √(γ(γ−1)T)`. Opacities are placeholders (unused by the EOS).
    /// Power laws in `(ρ, T)`, so the table's log-log interpolation is exact.
    fn ideal_gas_table() -> Table {
        let n = 8;
        let rho_grid: Vec<f64> = (0..n)
            .map(|i| 0.01 * 1000f64.powf(i as f64 / (n - 1) as f64)) // 0.01 … 10
            .collect();
        let t_grid: Vec<f64> = (0..n)
            .map(|j| 0.05 * 4000f64.powf(j as f64 / (n - 1) as f64)) // 0.05 … 200
            .collect();
        let mut p = Vec::new();
        let mut e = Vec::new();
        let mut cs = Vec::new();
        for &r in &rho_grid {
            for &t in &t_grid {
                p.push((GAMMA - 1.0) * r * t);
                e.push(t);
                cs.push((GAMMA * (GAMMA - 1.0) * t).sqrt());
            }
        }
        let one = vec![1.0; n * n];
        let json = serde_json::json!({
            "rho_grid": rho_grid,
            "T_grid": t_grid,
            "shape": [n, n],
            "fields": {
                "p": p, "e": e, "c_s": cs,
                "kappa_rosseland": one, "kappa_planck": one,
            },
        });
        Table::from_json(&json.to_string()).unwrap()
    }

    /// The tabulated EOS reproduces the analytic ideal gas it encodes, to the bisection tolerance.
    #[test]
    fn table_eos_matches_ideal_gas() {
        let ideal = IdealGas::new(GAMMA);
        let table = TableEos::new(ideal_gas_table());
        for &(rho, e) in &[(1.0, 2.5), (0.125, 2.0), (3.0, 0.7), (0.05, 50.0)] {
            assert_relative_eq!(
                table.pressure(rho, e),
                ideal.pressure(rho, e),
                max_relative = 1e-9
            );
            assert_relative_eq!(
                table.sound_speed(rho, e),
                ideal.sound_speed(rho, e),
                max_relative = 1e-9
            );
            assert_relative_eq!(
                table.dp_de(rho, e),
                ideal.dp_de(rho, e),
                max_relative = 1e-6
            );
        }
    }

    /// `energy_from_pressure` inverts `pressure`: round-tripping `(ρ, e) → p → e` recovers `e`.
    #[test]
    fn table_eos_pressure_energy_roundtrip() {
        let table = TableEos::new(ideal_gas_table());
        for &(rho, e) in &[(1.0, 2.5), (0.3, 12.0)] {
            let p = table.pressure(rho, e);
            assert_relative_eq!(table.energy_from_pressure(rho, p), e, max_relative = 1e-8);
        }
    }

    /// `TableEos::liquid_fraction` (Rung C): a table with `e = T` and a `liquid_frac` field returns
    /// that field at the inverted `T`; a table without the field returns `0`.
    #[allow(clippy::float_cmp)] // exact: the absent-field branch returns the literal 0.0
    #[test]
    fn liquid_fraction_reads_the_field_or_zero() {
        // e = T table (so temperature_from_energy is the identity) with a constant condensed fraction.
        let n: usize = 4;
        let rho_grid: Vec<f64> = (0..n).map(|i| 0.1 * 10f64.powf(i as f64)).collect(); // 0.1 … 100
        let t_grid: Vec<f64> = (0..n).map(|j| 300.0 * 2f64.powf(j as f64)).collect(); // 300 … 2400
        let e_field: Vec<f64> = (0..n * n).map(|idx| t_grid[idx % n]).collect(); // e = T, row-major
        let lin = |k: f64| vec![k; n * n];
        let json = serde_json::json!({
            "rho_grid": rho_grid, "T_grid": t_grid, "shape": [n, n],
            "fields": {
                "p": lin(1.0), "e": e_field,
                "c_s": lin(1.0), "kappa_rosseland": lin(1.0), "kappa_planck": lin(1.0),
                "liquid_frac": lin(0.3),
            },
        });
        let eos = TableEos::new(Table::from_json(&json.to_string()).unwrap());
        // e = T, so any e in range inverts to T = e and reads the constant 0.3.
        assert_relative_eq!(eos.liquid_fraction(1.0, 600.0), 0.3, max_relative = 1e-12);

        // The ideal-gas table omits liquid_frac → 0.
        let dry = TableEos::new(ideal_gas_table());
        assert_eq!(dry.liquid_fraction(1.0, 2.5), 0.0);
    }

    /// `TableEos::k_gas` (B-flux): a table with `e = T` and a `k_gas` field returns that field at the
    /// inverted `T`; a table without it returns `None` (the gas gets no conduction).
    #[test]
    fn k_gas_reads_the_field_or_none() {
        // e = T table (so temperature_from_energy is the identity) with a constant conductivity.
        let n: usize = 4;
        let rho_grid: Vec<f64> = (0..n).map(|i| 0.1 * 10f64.powf(i as f64)).collect(); // 0.1 … 100
        let t_grid: Vec<f64> = (0..n).map(|j| 300.0 * 2f64.powf(j as f64)).collect(); // 300 … 2400
        let e_field: Vec<f64> = (0..n * n).map(|idx| t_grid[idx % n]).collect(); // e = T, row-major
        let lin = |k: f64| vec![k; n * n];
        let json = serde_json::json!({
            "rho_grid": rho_grid, "T_grid": t_grid, "shape": [n, n],
            "fields": {
                "p": lin(1.0), "e": e_field,
                "c_s": lin(1.0), "kappa_rosseland": lin(1.0), "kappa_planck": lin(1.0),
                "k_gas": lin(0.05),
            },
        });
        let eos = TableEos::new(Table::from_json(&json.to_string()).unwrap());
        // e = T, so any e in range inverts to T = e and reads the constant 0.05.
        assert_relative_eq!(eos.k_gas(1.0, 600.0).unwrap(), 0.05, max_relative = 1e-12);

        // The ideal-gas table omits k_gas → None.
        let dry = TableEos::new(ideal_gas_table());
        assert!(dry.k_gas(1.0, 2.5).is_none());
    }

    /// `Eos::temperature` (B5a): the table inverts `e(ρ, T)` for `T` (here `e = T`, so it returns
    /// `e`), and the ideal gas returns `e` directly (reduced units `c_v = 1`). Both agree.
    #[test]
    fn temperature_inverts_energy() {
        let table = TableEos::new(ideal_gas_table());
        let ideal = IdealGas::new(GAMMA);
        for &(rho, e) in &[(1.0, 2.5), (0.3, 12.0), (5.0, 0.7)] {
            assert_relative_eq!(table.temperature(rho, e), e, max_relative = 1e-9);
            assert_relative_eq!(ideal.temperature(rho, e), e, max_relative = 1e-15);
        }
    }
}
