//! Rung-B exit criterion (B3c): the coupled gray-FLD substep reproduces the **Su–Olson
//! non-equilibrium Marshak wave** — the analytic two-temperature radiation-diffusion benchmark
//! Seth chose for the Marshak acceptance (memory: rung-b-radiation-decisions).
//!
//! ## The benchmark
//!
//! A cold, purely-absorbing half-space with **constant opacity** and a material heat capacity
//! `c_v = α T³` (which makes the material energy `e_m = (α/4) T⁴` linear in `a T⁴`, so the coupled
//! radiation/material system has a closed-form solution). A radiation source surface drives the
//! `x = 0` boundary through the incident-current (**Marshak**) condition. Su & Olson tabulate the
//! dimensionless radiation energy density `U` and material energy density `V` versus position `x`
//! and retarded time `τ`.
//!
//! - B. Su & G. L. Olson, *Benchmark results for the non-equilibrium Marshak diffusion problem*,
//!   JQSRT **56**, 337 (1996). Tables 1 (`U`) and 2 (`V`) for `ε = 0.1`.
//! - Reference values below are the published table, cross-checked against Frank Timmes'
//!   `suo02.f` / LANL ExactPack `suolson` (LA-CC-05-101).
//!
//! ## Mapping the benchmark onto the solver
//!
//! Dimensionless ↔ physical (Su & Olson): `x = √3 χ z`, `τ = 4 a c χ t / α`, `ε = 4 a / α`. We pick
//! consistent units `a = c = χ = 1` and `T_bc = 1`, so the incident energy density `e_inc = a T_bc⁴
//! = 1`, the dimensionless solutions are `U = E_rad` and `V = a T⁴ = T⁴`, the diffusion coefficient
//! is `c/(3χ) = 1/3`, and `τ = t/10` at `ε = 0.1` (`α = 40`).
//!
//! Two solver features this exercises (both added for this acceptance, see ADR-0006):
//! - the **Marshak** boundary `F(0) = (c/2)(e_inc − E)` — pure Dirichlet would pin `U(0) → 1`
//!   instead of the correct transient surface value;
//! - the **energy-based matter coupling** — `fld_substep` returns `Δe`, and we advance the internal
//!   energy and invert the `e_m = (α/4) T⁴` EOS for `T`. A `T += Δe/C_v` update would divide by the
//!   vanishing `C_v = α T³` at the cold wave front and blow up.
//!
//! Run in **Fick (`λ = 1/3`) mode**: Su–Olson is a *pure-diffusion* benchmark, and the production
//! Levermore–Pomraning limiter deliberately departs from it near the steep front. The limiter's
//! thin-gas behaviour is verified separately by the free-streaming-cap test.

use hydro1d::radiation::{Limiter, Medium, RadBc, RadConstants, fld_substep};

const RT3: f64 = 1.732_050_807_568_877_2;

/// Su & Olson 1996, `ε = 0.1`, `τ = 1.0`. Positions `x` and the tabulated `U` (Table 1) and `V`
/// (Table 2). (The full tables run to `x = 20`; we stop where the wave has died to `< 1e-4`.)
const XPOS: [f64; 10] = [0.0, 0.1, 0.25, 0.5, 0.75, 1.0, 2.5, 5.0, 7.5, 10.0];
const U_REF: [f64; 10] = [
    0.55182, 0.51412, 0.46198, 0.38541, 0.32046, 0.26564, 0.08147, 0.00961, 0.00097, 0.00009,
];
const V_REF: [f64; 10] = [
    0.32030, 0.29429, 0.25915, 0.20925, 0.16862, 0.13563, 0.03539, 0.00334, 0.00028, 0.00003,
];

/// Linear interpolation of a cell-centered field at dimensionless position `xq` (cell `j` sits at
/// `x = √3 (j+½) dz`). Clamps to the end cells, which is exact enough at `x = 0` and `x = x_max`.
fn interp_at_x(values: &[f64], dz: f64, xq: f64) -> f64 {
    let zq = xq / RT3;
    let fpos = zq / dz - 0.5;
    if fpos <= 0.0 {
        return values[0];
    }
    let j = fpos.floor() as usize;
    if j >= values.len() - 1 {
        return values[values.len() - 1];
    }
    let frac = fpos - j as f64;
    values[j] * (1.0 - frac) + values[j + 1] * frac
}

#[test]
fn fld_reproduces_suolson_marshak_wave() {
    // Dimensionless ε = 0.1, τ = 1.0  ->  α = 4a/ε = 40, t_final = 10τ = 10.
    let (a, c, chi) = (1.0, 1.0, 1.0);
    let eps = 0.1;
    let alpha = 4.0 * a / eps;
    let e_inc = a * 1.0_f64.powi(4); // T_bc = 1
    let consts = RadConstants { c, a };

    let n = 1200;
    let zmax = 8.0;
    let dz = zmax / n as f64;
    let t_final = 10.0;
    let nsteps = 2400;
    let dt = t_final / nsteps as f64;

    let dx = vec![dz; n];
    let center_spacing = vec![dz; n - 1];
    let chi_planck = vec![chi; n];
    let chi_ross = vec![chi; n];

    // Cold start as a tiny uniform equilibrium (T = floor, E = aT⁴): no spurious emission, and
    // it converges to the T(0) = 0 benchmark as the floor → 0.
    let tfloor = 1e-4_f64;
    let mut temp = vec![tfloor; n];
    let mut e_rad = vec![a * tfloor.powi(4); n];
    let mut e_mat = vec![0.25 * alpha * tfloor.powi(4); n]; // e_m = (α/4) T⁴

    for _ in 0..nsteps {
        let cv_vol: Vec<f64> = temp.iter().map(|&t| alpha * t.powi(3)).collect();
        let medium = Medium {
            dx: &dx,
            center_spacing: &center_spacing,
            temp: &temp,
            cv_vol: &cv_vol,
            chi_planck: &chi_planck,
            chi_ross: &chi_ross,
            source: None,
        };
        let delta_e = fld_substep(
            &medium,
            &mut e_rad,
            RadBc::Marshak(e_inc),
            RadBc::Reflecting,
            dt,
            consts,
            Limiter::Fick,
        );
        // Energy-based matter update: deposit Δe, then invert e_m = (α/4) T⁴ for the new T.
        for j in 0..n {
            e_mat[j] += delta_e[j];
            temp[j] = (4.0 * e_mat[j] / alpha).powf(0.25);
        }
    }

    // U = E_rad, V = a T⁴ (ener_in = 1). Compare to the published table.
    let v_field: Vec<f64> = temp.iter().map(|&t| a * t.powi(4)).collect();
    let mut max_du = 0.0_f64;
    let mut max_dv = 0.0_f64;
    for ((&xq, &ur), &vr) in XPOS.iter().zip(U_REF.iter()).zip(V_REF.iter()) {
        let un = interp_at_x(&e_rad, dz, xq);
        let vn = interp_at_x(&v_field, dz, xq);
        max_du = max_du.max((un - ur).abs());
        max_dv = max_dv.max((vn - vr).abs());
    }

    // Discretization error at this resolution is ~1.3e-3 (first-order in dt); the published table is
    // itself 5-decimal. 3e-3 is a comfortable, honest acceptance for the pure-diffusion benchmark.
    assert!(
        max_du < 3e-3 && max_dv < 3e-3,
        "Su–Olson Marshak mismatch: L∞(U)={max_du:e}, L∞(V)={max_dv:e}"
    );
}
