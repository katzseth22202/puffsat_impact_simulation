//! Energy-conserving nonlinear neutral kinetics; prescribed density, NOT hydro.

use crate::ResidualCurve;
use serde::{Deserialize, Serialize};

const H: [f64; 6] = [1.0, 0.0, 2.0, 1.0, 0.0, 2.0];
const O: [f64; 6] = [0.0, 1.0, 0.0, 1.0, 2.0, 1.0];

#[derive(Debug, Clone, Deserialize)]
pub struct Rate {
    pub reactants: [i32; 6],
    pub products: [i32; 6],
    /// SI molecule-based coefficient, including the selected scenario scale.
    pub a: f64,
    pub exponent: f64,
    pub activation_k: f64,
    pub colliders: Option<[f64; 6]>,
}

#[derive(Debug, Clone, Deserialize)]
pub struct Chemistry {
    pub kb: f64,
    pub molecule_mass: f64,
    pub reference_temp: f64,
    pub reference_log_n: [f64; 6],
    pub ground_k: [f64; 6],
    pub linear_cv: [f64; 6],
    pub vibrations: [Vec<f64>; 6],
    pub reactions: Vec<Rate>,
}

impl Chemistry {
    pub fn validate(&self) -> Result<(), String> {
        if !self.kb.is_finite()
            || self.kb <= 0.0
            || !self.molecule_mass.is_finite()
            || self.molecule_mass <= 0.0
            || !self.reference_temp.is_finite()
            || self.reference_temp <= 0.0
            || self
                .ground_k
                .iter()
                .chain(&self.reference_log_n)
                .any(|x| !x.is_finite())
            || self.linear_cv.iter().any(|x| !x.is_finite() || *x <= 0.0)
            || self
                .vibrations
                .iter()
                .flatten()
                .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err("invalid species thermochemistry".into());
        }
        for rate in &self.reactions {
            let nu = std::array::from_fn(|i| f64::from(rate.products[i] - rate.reactants[i]));
            if !rate.a.is_finite()
                || rate.a < 0.0
                || !rate.exponent.is_finite()
                || !rate.activation_k.is_finite()
                || rate
                    .reactants
                    .iter()
                    .chain(&rate.products)
                    .any(|x| !(0..=3).contains(x))
                || dot(H, nu) != 0.0
                || dot(O, nu) != 0.0
                || rate
                    .colliders
                    .is_some_and(|c| c.iter().any(|v| !v.is_finite() || *v < 0.0))
            {
                return Err("invalid or nonconservative reaction".into());
            }
        }
        Ok(())
    }
    pub fn species_thermo(&self, temp: f64) -> ([f64; 6], [f64; 6]) {
        let mut u = self.ground_k;
        let mut cv = self.linear_cv;
        for i in 0..6 {
            u[i] += self.linear_cv[i] * temp;
            for &theta in &self.vibrations[i] {
                let x = theta / temp;
                let q = (-x).exp();
                let d = -(-x).exp_m1();
                u[i] += theta * q / d;
                cv[i] += x * x * q / (d * d);
            }
        }
        (u, cv)
    }

    fn log_partition(&self, i: usize, temp: f64) -> f64 {
        self.linear_cv[i] * temp.ln()
            - self.ground_k[i] / temp
            - self.vibrations[i]
                .iter()
                .map(|t| (-(-t / temp).exp_m1()).ln())
                .sum::<f64>()
    }

    pub fn log_reference(&self, temp: f64) -> [f64; 6] {
        let d: [f64; 6] = std::array::from_fn(|i| {
            self.log_partition(i, temp) - self.log_partition(i, self.reference_temp)
        });
        std::array::from_fn(|i| self.reference_log_n[i] + d[i] - H[i] * d[0] - O[i] * d[1])
    }

    /// dy/dt and its analytic derivatives at fixed density. The collider
    /// derivative multiplies NET progress; it cancels at detailed balance.
    pub fn source(&self, y: [f64; 6], rho: f64, temp: f64) -> ([f64; 6], [[f64; 6]; 6], [f64; 6]) {
        let nf = rho / self.molecule_mass;
        let ln_n = y.map(|v| (nf * v).ln());
        let reference = self.log_reference(temp);
        let (u, _) = self.species_thermo(temp);
        let mut f = [0.0; 6];
        let mut jac = [[0.0; 6]; 6];
        let mut ft = [0.0; 6];
        for rate in &self.reactions {
            if rate.a == 0.0 {
                continue;
            }
            let nu: [f64; 6] =
                std::array::from_fn(|i| f64::from(rate.products[i] - rate.reactants[i]));
            let lk = rate.a.ln() + rate.exponent * temp.ln() - rate.activation_k / temp;
            let mut lf = lk;
            let mut lr = lk;
            let mut dlogk = 0.0;
            for i in 0..6 {
                lf += f64::from(rate.reactants[i]) * ln_n[i];
                lr += f64::from(rate.products[i]) * ln_n[i] - nu[i] * reference[i];
                dlogk += nu[i] * u[i] / (temp * temp);
            }
            let collider = rate.colliders.map_or(1.0, |c| {
                c.iter().zip(y).map(|(a, b)| a * b * nf).sum::<f64>()
            });
            if collider == 0.0 {
                continue;
            }
            lf += collider.ln();
            lr += collider.ln();
            // Already per original H2O formula unit, avoiding huge volumetric f.
            lf -= nf.ln();
            lr -= nf.ln();
            let qf = lf.exp();
            let qr = lr.exp();
            let net = if lf >= lr {
                -qf * (lr - lf).exp_m1()
            } else {
                qr * (lf - lr).exp_m1()
            };
            let qt = net * (rate.exponent / temp + rate.activation_k / (temp * temp)) + qr * dlogk;
            for i in 0..6 {
                f[i] += nu[i] * net;
                ft[i] += nu[i] * qt;
                for j in 0..6 {
                    let derivative = (qf * f64::from(rate.reactants[j])
                        - qr * f64::from(rate.products[j]))
                        / y[j]
                        + rate.colliders.map_or(0.0, |c| net * c[j] * nf / collider);
                    jac[i][j] += nu[i] * derivative;
                }
            }
        }
        (f, jac, ft)
    }
}

#[derive(Debug, Clone, Copy, Deserialize, Serialize)]
pub struct State {
    pub molecules: [f64; 4],
    pub temperature: f64,
}

#[derive(Debug)]
pub struct Context<'a> {
    pub chem: &'a Chemistry,
    pub residual: &'a ResidualCurve,
    pub atoms: [f64; 2],
    pub spectators: f64,
    pub energy_offset: f64,
    pub pressure_r_offset: f64,
}

impl Context<'_> {
    pub fn populations(&self, s: State) -> [f64; 6] {
        let m = s.molecules;
        [
            self.atoms[0] - 2.0 * m[0] - m[1] - 2.0 * m[3],
            self.atoms[1] - m[1] - 2.0 * m[2] - m[3],
            m[0],
            m[1],
            m[2],
            m[3],
        ]
    }

    /// Specific internal energy, pressure and frozen-composition heat capacity.
    pub fn thermo(&self, s: State, rho: f64) -> Result<(f64, f64, f64), String> {
        let y = self.populations(s);
        let t = s.temperature;
        if !t.is_finite()
            || !(1000.0..=2e6).contains(&t)
            || y.iter().any(|v| !v.is_finite() || *v <= 0.0)
            || !rho.is_finite()
            || rho <= 0.0
            || rho.ln() < self.residual.knots[0].log_rho
            || rho.ln() > self.residual.knots.last().unwrap().log_rho
        {
            return Err("parcel outside positive-population/thermo domain".into());
        }
        let (u, cv) = self.chem.species_thermo(t);
        let (v, d, _) = self.residual.evaluate(rho);
        let r = self.chem.kb / self.chem.molecule_mass;
        let a = 1000.0 * (1.0 - 1000.0 / t);
        let b = a - 0.5 * t * (1.0 - (1000.0 / t).powi(2));
        let energy =
            r * (dot(y, u) + 1.5 * self.spectators * t) + v[0] + v[2] * a + self.energy_offset;
        let pressure = rho
            * ((r * (y.iter().sum::<f64>() + self.spectators) + self.pressure_r_offset) * t + d[0]
                - t * d[1]
                + d[2] * b);
        let capacity = r * (dot(y, cv) + 1.5 * self.spectators) + v[2] * (1000.0 / t).powi(2);
        if !energy.is_finite() || !pressure.is_finite() || pressure <= 0.0 || capacity <= 0.0 {
            return Err("nonpositive or nonfinite parcel pressure/heat capacity".into());
        }
        Ok((energy, pressure, capacity))
    }

    fn pressure_t(&self, s: State, rho: f64) -> f64 {
        let (_, d, _) = self.residual.evaluate(rho);
        rho * (self.chem.kb / self.chem.molecule_mass
            * (self.populations(s).iter().sum::<f64>() + self.spectators)
            + self.pressure_r_offset
            - d[1]
            - 0.5 * d[2] * (1.0 - (1000.0 / s.temperature).powi(2)))
    }

    pub fn bond_change(&self, start: State, end: State) -> f64 {
        let a = self.populations(start);
        let b = self.populations(end);
        self.chem.kb / self.chem.molecule_mass
            * dot(std::array::from_fn(|i| a[i] - b[i]), self.chem.ground_k)
    }
}

fn dot(a: [f64; 6], b: [f64; 6]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn solve(mut a: [[f64; 5]; 5], mut b: [f64; 5]) -> Result<[f64; 5], String> {
    // Row scaling and partial pivoting; no unsafe code or external solver FFI.
    for i in 0..5 {
        let scale = a[i].iter().fold(0.0_f64, |v, x| v.max(x.abs()));
        if !scale.is_finite() || scale == 0.0 {
            return Err("singular Newton row".into());
        }
        for v in &mut a[i] {
            *v /= scale;
        }
        b[i] /= scale;
    }
    for j in 0..5 {
        let pivot = (j..5)
            .max_by(|&a0, &b0| a[a0][j].abs().total_cmp(&a[b0][j].abs()))
            .unwrap();
        a.swap(j, pivot);
        b.swap(j, pivot);
        if a[j][j].abs() < 1e-16 {
            return Err("unresolved Newton pivot".into());
        }
        for i in j + 1..5 {
            let factor = a[i][j] / a[j][j];
            let pivot_row = a[j];
            for (value, pivot_value) in a[i][j..].iter_mut().zip(&pivot_row[j..]) {
                *value -= factor * pivot_value;
            }
            b[i] -= factor * b[j];
        }
    }
    let mut x = [0.0; 5];
    for i in (0..5).rev() {
        x[i] = (b[i] - (i + 1..5).map(|j| a[i][j] * x[j]).sum::<f64>()) / a[i][i];
    }
    if x.iter().any(|v| !v.is_finite()) {
        return Err("nonfinite Newton correction".into());
    }
    Ok(x)
}

/// Fully implicit reaction/temperature update. Work uses the new pressure,
/// matching the energy residual exactly. Failed solves are rejected by caller.
pub fn backward_step(
    c: &Context<'_>,
    old: State,
    rho0: f64,
    rho1: f64,
    dt: f64,
    heat: f64,
) -> Result<(State, f64), String> {
    if !dt.is_finite() || dt <= 0.0 || !heat.is_finite() {
        return Err("invalid step".into());
    }
    let e0 = c.thermo(old, rho0)?.0;
    let dv = 1.0 / rho1 - 1.0 / rho0;
    let mut s = old;
    let tref = c.chem.reference_temp;
    for _ in 0..32 {
        let (e, p, cv) = c.thermo(s, rho1)?;
        let y = c.populations(s);
        let (f, jac, ft) = c.chem.source(y, rho1, s.temperature);
        let (u, _) = c.chem.species_thermo(s.temperature);
        let mut matrix = [[0.0; 5]; 5];
        let mut rhs = [0.0; 5];
        for i in 0..4 {
            rhs[i] = old.molecules[i] - s.molecules[i] + dt * f[i + 2];
            for j in 0..4 {
                matrix[i][j] = if i == j { 1.0 } else { 0.0 }
                    - dt * (jac[i + 2][j + 2]
                        - H[j + 2] * jac[i + 2][0]
                        - O[j + 2] * jac[i + 2][1]);
            }
            matrix[i][4] = -dt * ft[i + 2] * tref;
        }
        rhs[4] = e0 + heat - e - p * dv;
        let gas_r = c.chem.kb / c.chem.molecule_mass;
        for j in 0..4 {
            matrix[4][j] = gas_r
                * (u[j + 2] - H[j + 2] * u[0] - O[j + 2] * u[1]
                    + rho1 * s.temperature * (1.0 - H[j + 2] - O[j + 2]) * dv);
        }
        matrix[4][4] = (cv + c.pressure_t(s, rho1) * dv) * tref;
        let delta = solve(matrix, rhs)?;
        let small = delta[..4].iter().all(|d| d.abs() < 1e-11)
            && (delta[4] * tref).abs() < 1e-10 * s.temperature
            && rhs[4].abs() < 1e-10 * e0.abs().max(1.0);
        if small {
            return Ok((s, -p * dv));
        }
        let mut fraction = 1.0;
        let mut next = s;
        let mut found = false;
        for _ in 0..40 {
            next = State {
                molecules: std::array::from_fn(|i| s.molecules[i] + fraction * delta[i]),
                temperature: s.temperature + fraction * delta[4] * tref,
            };
            if c.thermo(next, rho1).is_ok() {
                found = true;
                break;
            }
            fraction *= 0.5;
        }
        if !found {
            return Err("Newton positivity line search failed".into());
        }
        s = next;
    }
    Err("Newton iteration limit".into())
}

#[derive(Debug, Clone, Deserialize)]
pub struct Knot {
    pub time_s: f64,
    pub rho: f64,
    /// Cumulative sampled parent delta(e)+mean(p)*delta(v), NOT chemical heat.
    pub forcing_j_kg: f64,
}

#[derive(Debug, Clone, Copy, Deserialize, Serialize)]
pub struct Settings {
    pub relative_tolerance: f64,
    pub absolute_population_tolerance: f64,
    pub maximum_step_s: f64,
    pub parent_forcing: bool,
}

#[derive(Debug, Serialize)]
pub struct Sample {
    pub time_s: f64,
    pub rho: f64,
    pub state: State,
    pub net_bond_return_j_kg: f64,
    pub work_j_kg: f64,
    pub forcing_j_kg: f64,
    pub energy_error_j_kg: f64,
}

#[derive(Debug, Serialize)]
pub struct Trajectory {
    pub status: String,
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    pub samples: Vec<Sample>,
}

pub fn integrate(
    c: &Context<'_>,
    initial: State,
    knots: &[Knot],
    settings: Settings,
) -> Trajectory {
    let mut result = Trajectory {
        status: "completed".into(),
        accepted_steps: 0,
        rejected_steps: 0,
        samples: vec![],
    };
    let outcome = integrate_inner(c, initial, knots, settings, &mut result);
    if let Err(error) = outcome {
        result.status = error;
    }
    result
}

fn integrate_inner(
    c: &Context<'_>,
    initial: State,
    knots: &[Knot],
    settings: Settings,
    result: &mut Trajectory,
) -> Result<(), String> {
    if knots.len() < 2
        || !settings.relative_tolerance.is_finite()
        || settings.relative_tolerance <= 0.0
        || !settings.absolute_population_tolerance.is_finite()
        || settings.absolute_population_tolerance <= 0.0
        || !settings.maximum_step_s.is_finite()
        || settings.maximum_step_s <= 0.0
    {
        return Err("invalid parcel settings/history".into());
    }
    c.chem.validate()?;
    let e0 = c.thermo(initial, knots[0].rho)?.0;
    let mut state = initial;
    let mut work = 0.0;
    let mut heat = 0.0;
    let mut dt = settings.maximum_step_s;
    result.samples.push(Sample {
        time_s: knots[0].time_s,
        rho: knots[0].rho,
        state,
        net_bond_return_j_kg: 0.0,
        work_j_kg: 0.0,
        forcing_j_kg: 0.0,
        energy_error_j_kg: 0.0,
    });
    for pair in knots.windows(2) {
        let (a, b) = (&pair[0], &pair[1]);
        let span = b.time_s - a.time_s;
        if !span.is_finite()
            || span <= 0.0
            || b.rho <= 0.0
            || !b.rho.is_finite()
            || !a.forcing_j_kg.is_finite()
            || !b.forcing_j_kg.is_finite()
        {
            return Err("invalid parcel density/forcing history".into());
        }
        let rho_at = |t: f64| (a.rho.ln() + (b.rho / a.rho).ln() * (t - a.time_s) / span).exp();
        let qdot = if settings.parent_forcing {
            (b.forcing_j_kg - a.forcing_j_kg) / span
        } else {
            0.0
        };
        let mut t = a.time_s;
        while b.time_s - t > 1e-15 * span {
            dt = dt.min(settings.maximum_step_s).min(b.time_s - t);
            if dt < 1e-20 || result.accepted_steps + result.rejected_steps > 200_000 {
                return Err("parcel step-size/work limit".into());
            }
            let r0 = rho_at(t);
            let rm = rho_at(t + 0.5 * dt);
            let r1 = rho_at(t + dt);
            let trial = (|| -> Result<(State, f64, f64), String> {
                let coarse = backward_step(c, state, r0, r1, dt, qdot * dt)?.0;
                let (mid, w1) = backward_step(c, state, r0, rm, 0.5 * dt, 0.5 * qdot * dt)?;
                let (fine, w2) = backward_step(c, mid, rm, r1, 0.5 * dt, 0.5 * qdot * dt)?;
                let mut error = (fine.temperature - coarse.temperature).abs()
                    / (settings.relative_tolerance * fine.temperature);
                for i in 0..4 {
                    error = error.max(
                        (fine.molecules[i] - coarse.molecules[i]).abs()
                            / (settings.absolute_population_tolerance
                                + settings.relative_tolerance
                                    * fine.molecules[i].abs().max(state.molecules[i].abs())),
                    );
                }
                Ok((fine, w1 + w2, error))
            })();
            if let Ok((fine, dw, error)) = trial {
                if error <= 1.0 && error.is_finite() {
                    state = fine;
                    work += dw;
                    heat += qdot * dt;
                    t += dt;
                    result.accepted_steps += 1;
                    dt *= (0.85 / error.max(1e-8).sqrt()).clamp(0.3, 2.0);
                    continue;
                }
                dt *= (0.7 / error.sqrt()).clamp(0.1, 0.5);
            } else {
                dt *= 0.25;
            }
            result.rejected_steps += 1;
        }
        let energy = c.thermo(state, b.rho)?.0;
        result.samples.push(Sample {
            time_s: b.time_s,
            rho: b.rho,
            state,
            net_bond_return_j_kg: c.bond_change(initial, state),
            work_j_kg: work,
            forcing_j_kg: heat,
            energy_error_j_kg: energy - e0 - work - heat,
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ResidualKnot;

    fn inert() -> (Chemistry, ResidualCurve) {
        (
            Chemistry {
                kb: 1.0,
                molecule_mass: 1.0,
                reference_temp: 6000.0,
                reference_log_n: [0.0; 6],
                ground_k: [0.0; 6],
                linear_cv: [1.5; 6],
                vibrations: std::array::from_fn(|_| vec![]),
                reactions: vec![],
            },
            ResidualCurve {
                knots: [-20.0, 20.0]
                    .map(|x| ResidualKnot {
                        log_rho: x,
                        values: [0.0; 3],
                        slopes: [0.0; 3],
                    })
                    .to_vec(),
            },
        )
    }

    #[test]
    fn frozen_expansion_converges_to_analytic_adiabat() {
        let (chem, residual) = inert();
        let context = Context {
            chem: &chem,
            residual: &residual,
            atoms: [2.0, 1.0],
            spectators: 0.0,
            energy_offset: 0.0,
            pressure_r_offset: 0.0,
        };
        let start = State {
            molecules: [0.1; 4],
            temperature: 6000.0,
        };
        let mut errors = vec![];
        for steps in [20, 40, 80] {
            let mut state = start;
            let mut work = 0.0;
            for i in 0..steps {
                let r0 = (-f64::from(i) / f64::from(steps)).exp();
                let r1 = (-f64::from(i + 1) / f64::from(steps)).exp();
                let (next, dw) =
                    backward_step(&context, state, r0, r1, 1.0 / f64::from(steps), 0.0).unwrap();
                state = next;
                work += dw;
            }
            assert_eq!(state.molecules, start.molecules);
            let expected = 6000.0 * (-2.0_f64 / 3.0).exp();
            errors.push((state.temperature - expected).abs());
            let delta = context.thermo(state, (-1.0_f64).exp()).unwrap().0
                - context.thermo(start, 1.0).unwrap().0;
            assert!((delta - work).abs() < 1e-7);
        }
        assert!((errors[0] / errors[1] - 2.0).abs() < 0.08);
        assert!((errors[1] / errors[2] - 2.0).abs() < 0.08);
    }

    #[test]
    fn reversible_isomer_relaxes_analytically_and_conserves_energy() {
        let (mut chem, residual) = inert();
        // Toy isomers O2 <=> 2 O would change particle count. Instead use
        // H2 + O2 <=> 2 OH: identical toy energies, exact scalar Riccati law
        // at equal H2/O2 populations, constant total particle count.
        chem.reactions.push(Rate {
            reactants: [0, 0, 1, 0, 1, 0],
            products: [0, 0, 0, 2, 0, 0],
            a: 1.0,
            exponent: 0.0,
            activation_k: 0.0,
            colliders: None,
        });
        let context = Context {
            chem: &chem,
            residual: &residual,
            atoms: [2.0, 2.0],
            spectators: 0.0,
            energy_offset: 0.0,
            pressure_r_offset: 0.0,
        };
        // H2=O2=x, OH=1-2x. dx/dt=(1-2x)^2-x^2;
        // (x-1/3)/(x-1) decays as exp(-2t).
        let start = State {
            molecules: [0.45, 0.1, 0.45, 0.1],
            temperature: 6000.0,
        };
        let mut errors = vec![];
        for steps in [40, 80, 160] {
            let mut state = start;
            for _ in 0..steps {
                state = backward_step(&context, state, 1.0, 1.0, 1.0 / f64::from(steps), 0.0)
                    .unwrap()
                    .0;
            }
            let z = (0.45 - 1.0 / 3.0) / (0.45 - 1.0) * (-2.0_f64).exp();
            let exact = (1.0 / 3.0 - z) / (1.0 - z);
            errors.push((state.molecules[0] - exact).abs());
            assert!((state.temperature - 6000.0).abs() < 1e-7);
            assert!(context.populations(state).iter().all(|y| *y > 0.0));
        }
        assert!((errors[0] / errors[1] - 2.0).abs() < 0.08);
        assert!((errors[1] / errors[2] - 2.0).abs() < 0.08);
    }
}
