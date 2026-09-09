//! Fixed molecules with charge-neutral atomic Saha equilibrium.
//! Analytic derivatives follow implicit charge conservation, not finite differences.

use crate::{FrozenGas, Thermo};
use serde::Deserialize;

#[derive(Debug, Clone, Deserialize)]
pub struct IonLadder {
    pub abundance: f64,
    pub energies: Vec<f64>,
    pub log_prefactors: Vec<f64>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct IonizingGas {
    /// Fixed molecules and heavy atoms; excludes free electrons and ion energy.
    pub base: FrozenGas,
    pub ladders: Vec<IonLadder>,
    pub boltzmann: f64,
    pub formula_mass: f64,
}

#[derive(Debug, Default)]
struct Moments {
    charge: f64,
    variance: f64,
    charge_t: f64,
    energy: f64,
    energy_charge: f64,
    energy_t: f64,
}

impl IonizingGas {
    fn moments(&self, temp: f64, ln_ne: f64) -> Moments {
        let mut result = Moments::default();
        for a in &self.ladders {
            if a.abundance == 0.0 {
                continue;
            }
            let logs: Vec<f64> = a
                .energies
                .iter()
                .enumerate()
                .map(|(z, e)| {
                    a.log_prefactors[z] + 1.5 * z as f64 * temp.ln()
                        - e / (self.boltzmann * temp)
                        - z as f64 * ln_ne
                })
                .collect();
            let peak = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let weights: Vec<f64> = logs.iter().map(|x| (x - peak).exp()).collect();
            let norm: f64 = weights.iter().sum();
            let mut q = 0.0;
            let mut energy = 0.0;
            let mut at = 0.0;
            for (z, &w) in weights.iter().enumerate() {
                let probability = w / norm;
                q += probability * z as f64;
                energy += probability * a.energies[z];
                at += probability
                    * (1.5 * z as f64 / temp + a.energies[z] / (self.boltzmann * temp * temp));
            }
            result.charge += a.abundance * q;
            result.energy += a.abundance * energy;
            // Centered moments avoid cancellation in nearly saturated ladders.
            for (z, &w) in weights.iter().enumerate() {
                let probability = a.abundance * w / norm;
                let dz = z as f64 - q;
                let de = a.energies[z] - energy;
                let da =
                    1.5 * z as f64 / temp + a.energies[z] / (self.boltzmann * temp * temp) - at;
                result.variance += probability * dz * dz;
                result.charge_t += probability * dz * da;
                result.energy_charge += probability * de * dz;
                result.energy_t += probability * de * da;
            }
        }
        result
    }

    fn equilibrium(&self, rho: f64, temp: f64) -> Moments {
        let maximum: f64 = self
            .ladders
            .iter()
            .map(|a| a.abundance * (a.energies.len() - 1) as f64)
            .sum();
        if maximum == 0.0 {
            return Moments::default();
        }
        let ln_nf = (rho / self.formula_mass).ln();
        let mut hi = ln_nf + maximum.ln();
        let mut lo = hi - 700.0;
        let mut x = hi;
        for _ in 0..100 {
            let s = self.moments(temp, x);
            let q = s.charge.max(1e-300);
            let error = x - ln_nf - q.ln();
            if error.abs() < 2e-13 {
                return s;
            }
            if error > 0.0 {
                hi = x;
            } else {
                lo = x;
            }
            let next = x - error / (1.0 + s.variance / q);
            x = if next > lo && next < hi {
                next
            } else {
                0.5 * (lo + hi)
            };
        }
        panic!("atomic charge-neutrality solve failed at rho={rho}, T={temp}");
    }

    pub fn at_temperature(&self, rho: f64, temp: f64) -> Thermo {
        let s = self.equilibrium(rho, temp);
        let denominator = (s.charge + s.variance).max(1e-300);
        let ln_ne_r = s.charge / denominator / rho;
        let ln_ne_t = s.charge_t / denominator;
        let q_r = -s.variance * ln_ne_r;
        let q_t = s.charge_t - s.variance * ln_ne_t;
        let factor = self.boltzmann / self.formula_mass;
        let r = self.base.gas_constant + factor * s.charge;
        let (base_energy, base_cv) = self.base.thermal(temp);
        Thermo {
            pressure: rho * r * temp,
            energy: base_energy + 1.5 * factor * s.charge * temp + s.energy / self.formula_mass,
            cv: base_cv
                + 1.5 * factor * (s.charge + temp * q_t)
                + (s.energy_t - s.energy_charge * ln_ne_t) / self.formula_mass,
            pressure_t: rho * (r + temp * factor * q_t),
            pressure_r: temp * (r + rho * factor * q_r),
            energy_r: 1.5 * factor * temp * q_r - s.energy_charge * ln_ne_r / self.formula_mass,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn hydrogen() -> IonizingGas {
        // Dimensionless one-electron ladder: at rho=T=k=m=1, K=0.5
        // gives ion fraction exactly 1/2 by alpha^2/(1-alpha)=K/n.
        IonizingGas {
            base: FrozenGas {
                gas_constant: 1.0,
                linear_cv: 1.5,
                chemical_energy: 7.0,
                vibrations: vec![],
                bond_store: 7.0,
                ion_store: 0.0,
            },
            ladders: vec![IonLadder {
                abundance: 1.0,
                energies: vec![0.0, 2.0],
                log_prefactors: vec![0.0, 0.5_f64.ln() + 2.0],
            }],
            boltzmann: 1.0,
            formula_mass: 1.0,
        }
    }

    #[test]
    fn single_stage_matches_quadratic_saha_solution() {
        let gas = hydrogen();
        let s = gas.at_temperature(1.0, 1.0);
        assert_relative_eq!(s.pressure, 1.5, max_relative = 1e-12);
        assert_relative_eq!(s.energy, 10.25, max_relative = 1e-12);
    }

    #[test]
    fn constrained_derivatives_obey_maxwell_and_finite_differences() {
        let mut gas = hydrogen();
        gas.ladders.push(IonLadder {
            abundance: 0.4,
            energies: vec![0.0, 3.0, 8.0],
            log_prefactors: vec![0.0, 1.0, 2.0],
        });
        gas.base.gas_constant += 0.4;
        gas.base.linear_cv += 0.6;
        for (rho, t) in [(0.01, 0.2), (1.0, 1.0), (100.0, 5.0), (0.1, 1000.0)] {
            let dr = 1e-4 * rho;
            let dt = 1e-4 * t;
            let s = gas.at_temperature(rho, t);
            let a = gas.at_temperature(rho - dr, t);
            let b = gas.at_temperature(rho + dr, t);
            let c = gas.at_temperature(rho, t - dt);
            let d = gas.at_temperature(rho, t + dt);
            assert_relative_eq!(
                s.pressure_r,
                (b.pressure - a.pressure) / (2.0 * dr),
                max_relative = 1e-7
            );
            assert_relative_eq!(
                s.pressure_t,
                (d.pressure - c.pressure) / (2.0 * dt),
                max_relative = 1e-7
            );
            assert_relative_eq!(
                s.cv,
                (d.energy - c.energy) / (2.0 * dt),
                max_relative = 1e-7
            );
            assert_relative_eq!(
                s.energy_r,
                (b.energy - a.energy) / (2.0 * dr),
                max_relative = 1e-5,
                epsilon = 1e-10
            );
            assert_relative_eq!(
                rho * rho * s.energy_r,
                s.pressure - t * s.pressure_t,
                max_relative = 1e-8,
                epsilon = 1e-10
            );
            assert!(s.sound_speed(rho) > 0.0);
        }
    }
}
