//! Parcel-frozen gas thermodynamics retaining the same dense Helmholtz residual.
//! No kinetics, radiation or new hydrodynamic discretization is introduced here.

use hydro1d::eos::{Eos, TableEos};
use serde::Deserialize;
use std::sync::Arc;

#[derive(Debug, Clone, Deserialize)]
pub struct ResidualKnot {
    pub log_rho: f64,
    /// Residual energy, entropy and heat capacity at 1000 K.
    pub values: [f64; 3],
    /// Derivatives with respect to log density, from pressure identities.
    pub slopes: [f64; 3],
}

#[derive(Debug, Clone, Deserialize)]
pub struct ResidualCurve {
    pub knots: Vec<ResidualKnot>,
}

impl ResidualCurve {
    /// One cubic Hermite potential and its derivatives; pressure and acoustics
    /// are derived from it, not independently interpolated residual fields.
    pub fn evaluate(&self, rho: f64) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let x = rho.ln();
        let n = self.knots.len();
        assert!(n >= 2 && x >= self.knots[0].log_rho && x <= self.knots[n - 1].log_rho);
        let hi = self
            .knots
            .partition_point(|k| k.log_rho < x)
            .clamp(1, n - 1);
        let (a, b) = (&self.knots[hi - 1], &self.knots[hi]);
        let h = b.log_rho - a.log_rho;
        let u = (x - a.log_rho) / h;
        let mut value = [0.0; 3];
        let mut first = [0.0; 3];
        let mut second = [0.0; 3];
        for j in 0..3 {
            let (y0, y1, m0, m1) = (a.values[j], b.values[j], h * a.slopes[j], h * b.slopes[j]);
            // Horner form reduces cancellation for nearby states.
            let c2 = 3.0 * (y1 - y0) - 2.0 * m0 - m1;
            let c3 = 2.0 * (y0 - y1) + m0 + m1;
            value[j] = y0 + u * (m0 + u * (c2 + u * c3));
            first[j] = (m0 + u * (2.0 * c2 + 3.0 * u * c3)) / h;
            second[j] = (2.0 * c2 + 6.0 * u * c3) / (h * h);
        }
        (value, first, second)
    }
}

#[derive(Debug, Clone, Deserialize)]
pub struct Vibration {
    pub theta_k: f64,
    pub amplitude_j_kg: f64,
}

#[derive(Debug, Clone, Deserialize)]
pub struct FrozenGas {
    pub gas_constant: f64,
    pub linear_cv: f64,
    pub chemical_energy: f64,
    pub vibrations: Vec<Vibration>,
    pub bond_store: f64,
    pub ion_store: f64,
}

impl FrozenGas {
    pub fn thermal(&self, temp: f64) -> (f64, f64) {
        let mut energy = self.linear_cv * temp + self.chemical_energy;
        let mut cv = self.linear_cv;
        for mode in &self.vibrations {
            let x = mode.theta_k / temp;
            let denominator = -(-x).exp_m1();
            let q = (-x).exp();
            energy += mode.amplitude_j_kg * q / denominator;
            cv +=
                mode.amplitude_j_kg * mode.theta_k * q / (temp * temp * denominator * denominator);
        }
        (energy, cv)
    }
}

#[derive(Debug, Clone)]
pub struct FrozenCell {
    pub gas: FrozenGas,
    pub residual: Arc<ResidualCurve>,
    pub storage_shift: f64,
    /// Constant energy and entropy-pressure corrections to match the parent
    /// table at the exact switching state; reported separately by the runner.
    pub energy_offset: f64,
    pub pressure_r_offset: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct Thermo {
    pub pressure: f64,
    pub energy: f64,
    pub cv: f64,
    pub pressure_t: f64,
    pub pressure_r: f64,
    pub energy_r: f64,
}

impl Thermo {
    pub fn sound_speed(&self, rho: f64) -> f64 {
        let cs2 = self.pressure_r
            + self.pressure_t * (self.pressure / (rho * rho) - self.energy_r) / self.cv;
        assert!(
            self.cv > 0.0 && cs2 > 0.0 && cs2.is_finite(),
            "unstable frozen closure: {self:?}"
        );
        cs2.sqrt()
    }
}

impl FrozenCell {
    pub fn at_temperature(&self, rho: f64, temp: f64) -> Thermo {
        let (v, d, dd) = self.residual.evaluate(rho);
        let t0 = 1000.0;
        let a = t0 * (1.0 - t0 / temp);
        let b = a - 0.5 * temp * (1.0 - (t0 / temp).powi(2));
        let bt = -0.5 * (1.0 - (t0 / temp).powi(2));
        let fx = d[0] - temp * d[1] + d[2] * b;
        let fxx = dd[0] - temp * dd[1] + dd[2] * b;
        let (gas_e, gas_cv) = self.gas.thermal(temp);
        let r = self.gas.gas_constant + self.pressure_r_offset;
        Thermo {
            pressure: rho * (r * temp + fx),
            energy: gas_e + v[0] + v[2] * a + self.storage_shift + self.energy_offset,
            cv: gas_cv + v[2] * (t0 / temp).powi(2),
            pressure_t: rho * (r - d[1] + d[2] * bt),
            pressure_r: r * temp + fx + fxx,
            energy_r: (d[0] + d[2] * a) / rho,
        }
    }

    /// A safeguarded solve; callers check the temperature domain after every
    /// hydro step. Clamping an intermediate Newton query is not accepted as a
    /// completed in-domain trajectory.
    pub fn invert_energy(&self, rho: f64, energy: f64) -> f64 {
        let (mut lo, mut hi) = (1000.0, 2e6);
        if energy <= self.at_temperature(rho, lo).energy {
            return lo;
        }
        if energy >= self.at_temperature(rho, hi).energy {
            return hi;
        }
        let mut t = ((energy - self.gas.chemical_energy - self.storage_shift - self.energy_offset)
            / self.gas.linear_cv)
            .clamp(lo, hi);
        for _ in 0..40 {
            let s = self.at_temperature(rho, t);
            let error = s.energy - energy;
            if error.abs() < 1e-11 * energy.abs().max(1.0) {
                return t;
            }
            if error > 0.0 {
                hi = t;
            } else {
                lo = t;
            }
            let next = t - error / s.cv;
            t = if next > lo && next < hi {
                next
            } else {
                0.5 * (lo + hi)
            };
        }
        t
    }
}

#[derive(Debug, Clone)]
pub enum ParcelEos {
    Equilibrium(Arc<TableEos>),
    Frozen(FrozenCell),
    Indexed(Vec<ParcelEos>),
}

impl Eos for ParcelEos {
    fn for_cell(&self, cell: usize) -> &Self {
        match self {
            Self::Indexed(cells) => &cells[cell],
            _ => self,
        }
    }
    fn pressure(&self, rho: f64, e: f64) -> f64 {
        match self {
            Self::Equilibrium(v) => v.pressure(rho, e),
            Self::Frozen(v) => v.at_temperature(rho, v.invert_energy(rho, e)).pressure,
            Self::Indexed(_) => panic!("indexed EOS requires a cell index"),
        }
    }
    fn sound_speed(&self, rho: f64, e: f64) -> f64 {
        match self {
            Self::Equilibrium(v) => v.sound_speed(rho, e),
            Self::Frozen(v) => v
                .at_temperature(rho, v.invert_energy(rho, e))
                .sound_speed(rho),
            Self::Indexed(_) => panic!("indexed EOS requires a cell index"),
        }
    }
    fn temperature(&self, rho: f64, e: f64) -> f64 {
        match self {
            Self::Equilibrium(v) => v.temperature(rho, e),
            Self::Frozen(v) => v.invert_energy(rho, e),
            Self::Indexed(_) => panic!("indexed EOS requires a cell index"),
        }
    }
    fn dp_de(&self, rho: f64, e: f64) -> f64 {
        match self {
            Self::Equilibrium(v) => v.dp_de(rho, e),
            Self::Frozen(v) => {
                let s = v.at_temperature(rho, v.invert_energy(rho, e));
                s.pressure_t / s.cv
            }
            Self::Indexed(_) => panic!("indexed EOS requires a cell index"),
        }
    }
    fn energy_from_pressure(&self, rho: f64, p: f64) -> f64 {
        match self {
            Self::Equilibrium(v) => v.energy_from_pressure(rho, p),
            _ => panic!("frozen restart must specify conserved energy explicitly"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use hydro1d::kernel::{Boundary, LagrangianState, Tube, Viscosity};

    fn monatomic(r: f64) -> FrozenCell {
        FrozenCell {
            gas: FrozenGas {
                gas_constant: r,
                linear_cv: 1.5 * r,
                chemical_energy: 5e7,
                vibrations: vec![],
                bond_store: 5e7,
                ion_store: 0.0,
            },
            residual: Arc::new(ResidualCurve {
                knots: vec![
                    ResidualKnot {
                        log_rho: -20.0,
                        values: [0.0; 3],
                        slopes: [0.0; 3],
                    },
                    ResidualKnot {
                        log_rho: 20.0,
                        values: [0.0; 3],
                        slopes: [0.0; 3],
                    },
                ],
            }),
            storage_shift: 4e6,
            energy_offset: 0.0,
            pressure_r_offset: 0.0,
        }
    }

    #[test]
    fn frozen_monatomic_limit_has_inert_store_and_gamma_five_thirds() {
        let cell = monatomic(1000.0);
        let eos = ParcelEos::Frozen(cell.clone());
        for (rho, t) in [(1.0, 10000.0), (100.0, 5000.0), (0.01, 100000.0)] {
            let state = cell.at_temperature(rho, t);
            assert_relative_eq!(state.pressure, rho * 1000.0 * t, max_relative = 1e-12);
            assert_relative_eq!(eos.temperature(rho, state.energy), t, max_relative = 1e-10);
            assert_relative_eq!(
                state.sound_speed(rho).powi(2),
                5.0 / 3.0 * 1000.0 * t,
                max_relative = 1e-12
            );
        }
    }

    #[test]
    fn distinct_parcel_closures_preserve_a_stationary_pressure_contact() {
        let a = monatomic(1000.0);
        let mut b = monatomic(2000.0);
        b.gas.linear_cv = 6000.0; // Different gamma, so incorrect cell dispatch cannot pass.
        let e1 = a.at_temperature(1.0, 10000.0).energy;
        let e2 = b.at_temperature(1.0, 5000.0).energy;
        let mut tube = Tube::from_lagrangian(
            LagrangianState {
                positions: vec![0.0, 0.5, 1.0],
                node_velocities: vec![0.0; 3],
                cell_masses: vec![0.5; 2],
                specific_energies: vec![e1, e2],
            },
            ParcelEos::Indexed(vec![ParcelEos::Frozen(a), ParcelEos::Frozen(b)]),
            Boundary::Wall,
            Boundary::Wall,
            Viscosity::VON_NEUMANN_RICHTMYER,
        );
        let energy = tube.total_energy();
        tube.run_to(1e-3);
        assert_relative_eq!(tube.pressure(0), 1e7, max_relative = 1e-12);
        assert_relative_eq!(tube.pressure(1), 1e7, max_relative = 1e-12);
        assert_relative_eq!(tube.total_energy(), energy, max_relative = 1e-12);
        assert!(tube.total_momentum().abs() < 1e-10);
    }

    #[test]
    fn residual_potential_obeys_maxwell_and_its_analytic_derivatives() {
        let mut cell = monatomic(1000.0);
        cell.residual = Arc::new(ResidualCurve {
            knots: [0.0, 2.0]
                .into_iter()
                .map(|x: f64| ResidualKnot {
                    log_rho: x,
                    values: [
                        2e6 * (1.0 + x + 0.5 * x * x),
                        100.0 * (x + 0.25 * x * x),
                        200.0 * (1.0 + 0.2 * x),
                    ],
                    slopes: [2e6 * (1.0 + x), 100.0 * (1.0 + 0.5 * x), 40.0],
                })
                .collect(),
        });
        let rho = 0.7_f64.exp();
        let temp = 5000.0;
        let dr = 1e-5 * rho;
        let dt = 1e-5 * temp;
        let s = cell.at_temperature(rho, temp);
        let a = cell.at_temperature(rho - dr, temp);
        let b = cell.at_temperature(rho + dr, temp);
        let c = cell.at_temperature(rho, temp - dt);
        let d = cell.at_temperature(rho, temp + dt);
        assert_relative_eq!(
            s.pressure_r,
            (b.pressure - a.pressure) / (2.0 * dr),
            max_relative = 1e-8
        );
        assert_relative_eq!(
            s.energy_r,
            (b.energy - a.energy) / (2.0 * dr),
            max_relative = 1e-8
        );
        assert_relative_eq!(
            s.pressure_t,
            (d.pressure - c.pressure) / (2.0 * dt),
            max_relative = 1e-8
        );
        assert_relative_eq!(
            s.cv,
            (d.energy - c.energy) / (2.0 * dt),
            max_relative = 1e-8
        );
        assert_relative_eq!(
            rho * rho * s.energy_r,
            s.pressure - temp * s.pressure_t,
            max_relative = 1e-12
        );
        assert!(s.sound_speed(rho) > 0.0);
    }

    #[test]
    fn vibrational_heat_capacity_is_the_derivative_of_its_energy() {
        let mut gas = monatomic(1000.0).gas;
        gas.vibrations.push(Vibration {
            theta_k: 3000.0,
            amplitude_j_kg: 1.2e6,
        });
        for t in [1000.0, 5000.0, 100000.0] {
            let dt = 1e-4 * t;
            let cv = (gas.thermal(t + dt).0 - gas.thermal(t - dt).0) / (2.0 * dt);
            assert_relative_eq!(gas.thermal(t).1, cv, max_relative = 1e-8);
        }
    }
}
