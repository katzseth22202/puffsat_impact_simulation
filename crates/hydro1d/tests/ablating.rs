//! Rung E (E1): the quasi-steady **ablating wall** — the ablation mass/energy source (ADR-0014).
//!
//! The ablating wall is the rigid coupled bounce ([`CoupledBounce`]) plus a wall **mass source** —
//! the inverse of the `CondensingBounce` wall sink. Each step the incoming wall flux `q_in`
//! (radiative + conductive) boils off `ṁ = q_in / Q*` of ablator, injected as cold vapor into the
//! wall cell. Acceptance tests, written first:
//!
//!  - **`Q* → ∞` recovers the rigid floor** — with an infinite heat of ablation `ṁ → 0`, so the
//!    ablating bounce reduces *exactly* to the rigid [`CoupledBounce`] `e_eff`. This is the
//!    conservative-floor consistency gate (ADR-0014/0013): every later recovery is a gain measured
//!    against a trusted baseline.
//!  - **the incoming wall flux fully drives ablation** — the quasi-steady surface energy balance
//!    `q_in = ṁ·Q*`, i.e. `loss_ablation == ablated_mass·Q* == loss_radiative_wall + loss_conductive`.
//!  - **the ablation rate is monotone in `Q*`** (and vanishes as `Q* → ∞`).
//!  - **the ablated mass is a small fraction of the cloud** — the quasi-steady regime (~1.5 %,
//!    ADR-0014), so the bounce is a perturbation of the rigid floor, not a different problem.

use hydro1d::conduction::Solid;
use hydro1d::eos::TableEos;
use hydro1d::kernel::{AblatingBounce, Ablation, CoupledBounce, Tube, Viscosity};
use hydro1d::radiation::{Limiter, RadConstants};
use tables::Table;

const GAMMA: f64 = 5.0 / 3.0; // monatomic; e = T, p = (γ−1)ρT, c_s = √(γ(γ−1)T)
const CONSTS: RadConstants = RadConstants { c: 3.0, a: 1.0 }; // normalized; radiation active
const LIMITER: Limiter = Limiter::LevermorePomraning;
const T_VAPOR: f64 = 0.01; // cold injected-vapor temperature (the slug's T₀ — conservative)

/// An ideal-gas table (`p = (γ−1)ρT`, `e = T`, `c_v = 1`) with a **constant, nonzero** opacity, so
/// the gray-FLD radiation step is active and the wall absorbs a real radiative flux (the `q_in` that
/// drives ablation). Power laws in `(ρ, T)` ⇒ the log-log interpolation is exact; `κ = 1` is flat.
fn radiating_ideal_table() -> Table {
    let n = 8;
    let rho_grid: Vec<f64> = (0..n)
        .map(|i| 1e-3 * 1e5_f64.powf(i as f64 / (n - 1) as f64)) // 1e-3 … 100
        .collect();
    let t_grid: Vec<f64> = (0..n)
        .map(|j| 1e-3 * 1e6_f64.powf(j as f64 / (n - 1) as f64)) // 1e-3 … 1000
        .collect();
    let (mut p, mut e, mut cs) = (Vec::new(), Vec::new(), Vec::new());
    for &r in &rho_grid {
        for &t in &t_grid {
            p.push((GAMMA - 1.0) * r * t);
            e.push(t);
            cs.push((GAMMA * (GAMMA - 1.0) * t).sqrt());
        }
    }
    let one = vec![1.0; n * n]; // κ_R = κ_P = 1 m²/kg ⇒ χ = ρ
    let json = serde_json::json!({
        "rho_grid": rho_grid,
        "T_grid": t_grid,
        "shape": [n, n],
        "fields": { "p": p, "e": e, "c_s": cs, "kappa_rosseland": one, "kappa_planck": one },
    });
    Table::from_json(&json.to_string()).unwrap()
}

/// A cold water-like slug (ρ = 1, v = 1, M ≈ 9.5 at T₀ = 0.01) coasting into the wall, on `table`.
fn slug(table: &Table) -> Tube<TableEos> {
    Tube::slug_si(
        200,
        1.0,
        1.0,
        1.0,
        T_VAPOR,
        TableEos::new(table.clone()),
        Viscosity::VON_NEUMANN_RICHTMYER,
    )
}

/// `Q* → ∞` ⇒ `ṁ → 0` ⇒ the ablating bounce *is* the rigid coupled bounce, to round-off. The
/// conservative-floor consistency gate (ADR-0013/0014).
#[test]
fn q_star_infinite_recovers_rigid_floor() {
    let table = radiating_ideal_table();
    let rigid = CoupledBounce::new(slug(&table), None, CONSTS, LIMITER).run();
    let ablating = AblatingBounce::new(
        slug(&table),
        None,
        CONSTS,
        LIMITER,
        Ablation::new(1e30, T_VAPOR),
    )
    .run();
    let d = (ablating.bounce.e_eff - rigid.bounce.e_eff).abs();
    assert!(
        d < 1e-9,
        "Q*→∞ must recover the rigid floor: |Δe_eff| = {d:e}"
    );
}

/// The quasi-steady surface energy balance: the incoming wall flux is exactly the energy spent
/// ablating, and the ablated mass is that energy over `Q*`. With `wall = None` (the realistic high-v
/// config) the only incoming flux is radiative, so all three quantities coincide to round-off.
#[test]
fn incoming_flux_drives_ablation() {
    let table = radiating_ideal_table();
    let q_star = 50.0;
    let r = AblatingBounce::new(
        slug(&table),
        None,
        CONSTS,
        LIMITER,
        Ablation::new(q_star, T_VAPOR),
    )
    .run();
    assert!(
        r.ablated_mass > 0.0,
        "finite Q* must ablate mass: {}",
        r.ablated_mass
    );

    // loss_ablation == ablated_mass · Q* (the ṁ = q_in/Q* definition).
    let rel_qstar = (r.loss_ablation - r.ablated_mass * q_star).abs() / r.loss_ablation;
    assert!(
        rel_qstar < 1e-12,
        "loss_ablation ≠ ablated_mass·Q*: rel = {rel_qstar:e}"
    );

    // loss_ablation == the incoming wall flux (radiative + conductive). wall = None ⇒ conductive 0.
    let q_in = r.loss_radiative_wall + r.loss_conductive;
    let rel_flux = (r.loss_ablation - q_in).abs() / q_in;
    assert!(
        rel_flux < 1e-12,
        "loss_ablation ≠ incoming wall flux: rel = {rel_flux:e}"
    );
}

/// A larger heat of ablation boils off less mass (and `Q* → ∞` boils off ≈ nothing) — the ṁ ∝ 1/Q*
/// direction of the surface energy balance.
#[test]
fn ablation_rate_is_monotone_in_q_star() {
    let table = radiating_ideal_table();
    let run = |q_star: f64| {
        AblatingBounce::new(
            slug(&table),
            None,
            CONSTS,
            LIMITER,
            Ablation::new(q_star, T_VAPOR),
        )
        .run()
        .ablated_mass
    };
    let (m25, m100, m_inf) = (run(25.0), run(100.0), run(1e30));
    assert!(m25 > m100, "more ablation at smaller Q*: {m25} vs {m100}");
    assert!(m100 > m_inf, "monotone in Q*: {m100} vs {m_inf}");
    assert!(
        m_inf < 1e-9 * m25,
        "Q*→∞ ablates ≈ nothing: {m_inf} vs {m25}"
    );
}

/// The ablated mass is a small fraction of the cloud (the cloud mass per area is ρ·L = 1 here), so
/// the quasi-steady picture holds and the ablating bounce is a perturbation of the rigid floor
/// (ADR-0014: ~1.5 %). A physically-sized `Q*` (here 50 in normalized units) is assumed.
#[test]
fn ablated_mass_is_small_quasi_steady() {
    let table = radiating_ideal_table();
    let r = AblatingBounce::new(
        slug(&table),
        None,
        CONSTS,
        LIMITER,
        Ablation::new(50.0, T_VAPOR),
    )
    .run();
    let cloud_mass = 1.0; // ρ·L per unit area
    assert!(
        r.ablated_mass > 0.0 && r.ablated_mass < 0.2 * cloud_mass,
        "ablated mass {} outside the quasi-steady band (0, 0.2·cloud)",
        r.ablated_mass
    );
}

// ---- E2: blowing correction on the conductive flux (verify-and-bound) ----
//
// Blowing — the injected vapor thickening the boundary layer — cuts the conductive flux into the
// plate (ADR-0014). We model it as a monotone factor `φ(B) ∈ (0,1]` attenuating the conducted heat
// delivered to the plate (the vapor curtain convects the intercepted fraction back into the gas),
// with the dimensionless blowing rate `B ∝ ablated_mass` (so `φ → 1` as blowing → 0). The key
// finding it bounds: at the *science anchors* conduction is off (the high-v table carries no
// `k_gas`), so blowing is identically null there — the live recovery is vapor shielding (E3).

/// An ideal-gas radiating table that *also* carries a constant gas conductivity `k_gas`, so the
/// gas-side conduction operator engages (blowing has something to reduce).
fn radiating_conducting_table() -> Table {
    let n = 8;
    let rho_grid: Vec<f64> = (0..n)
        .map(|i| 1e-3 * 1e5_f64.powf(i as f64 / (n - 1) as f64))
        .collect();
    let t_grid: Vec<f64> = (0..n)
        .map(|j| 1e-3 * 1e6_f64.powf(j as f64 / (n - 1) as f64))
        .collect();
    let (mut p, mut e, mut cs) = (Vec::new(), Vec::new(), Vec::new());
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
            "kappa_rosseland": one, "kappa_planck": one, "k_gas": one,
        },
    });
    Table::from_json(&json.to_string()).unwrap()
}

/// A cold, high-effusivity conducting plate behind the wall (the conduction sink). Coarse (40 cells)
/// to keep the conducting-wall bounces cheap — these tests check the blowing *direction*, not a
/// converged number.
fn conducting_wall() -> Solid {
    Solid::new(40, 1.0, T_VAPOR, 0.1, 10.0)
}

/// A small slug (60 cells) for the conducting-wall tests: the coupled conduction over a gas+solid
/// mesh is the expensive path, so these run on a coarse grid (the claims are qualitative).
fn small_slug(table: &Table) -> Tube<TableEos> {
    Tube::slug_si(
        60,
        1.0,
        1.0,
        1.0,
        T_VAPOR,
        TableEos::new(table.clone()),
        Viscosity::VON_NEUMANN_RICHTMYER,
    )
}

fn run_with(
    table: &Table,
    wall: Option<Solid>,
    ablation: Ablation,
) -> hydro1d::kernel::AblatingBounceResult {
    AblatingBounce::new(small_slug(table), wall, CONSTS, LIMITER, ablation).run()
}

/// Blowing strictly cuts the conductive wall loss: at a finite `Q*` and a conducting wall, turning
/// blowing on lowers `loss_conductive` below the unblown ablating bounce.
#[test]
fn blowing_reduces_conductive_loss() {
    let table = radiating_conducting_table();
    let unblown = run_with(
        &table,
        Some(conducting_wall()),
        Ablation::new(50.0, T_VAPOR),
    );
    let blown = run_with(
        &table,
        Some(conducting_wall()),
        Ablation::new(50.0, T_VAPOR).with_blowing(50.0),
    );
    assert!(
        unblown.loss_conductive > 0.0,
        "conduction must be active for the test: {}",
        unblown.loss_conductive
    );
    assert!(
        blown.loss_conductive < unblown.loss_conductive,
        "blowing must cut conduction: {} vs unblown {}",
        blown.loss_conductive,
        unblown.loss_conductive
    );
}

/// In the strong-blowing limit the conductive channel is largely choked off — the vapor curtain
/// insulates the plate (the blowing → ∞ limit `φ → 0`). The residual is the early conduction that
/// occurs *before* the curtain accumulates (φ starts at 1 and falls as ablation builds — blowing
/// cannot retroactively block heat already delivered), so the bar is a >90 % cut, not zero.
#[test]
fn strong_blowing_chokes_conduction() {
    let table = radiating_conducting_table();
    let unblown = run_with(
        &table,
        Some(conducting_wall()),
        Ablation::new(50.0, T_VAPOR),
    );
    let choked = run_with(
        &table,
        Some(conducting_wall()),
        Ablation::new(50.0, T_VAPOR).with_blowing(1e6),
    );
    assert!(
        choked.loss_conductive < 0.1 * unblown.loss_conductive,
        "strong blowing must choke conduction by >90%: {} vs unblown {}",
        choked.loss_conductive,
        unblown.loss_conductive
    );
}

/// The bound: at the realistic high-v config (`wall = None`, conduction off) blowing acts on
/// nothing — `loss_conductive` is identically 0 and `e_eff` is unchanged by the blowing coefficient.
/// So blowing is foreclosed from the headline recovery; that is vapor shielding (E3).
#[test]
fn blowing_is_null_without_conduction() {
    let table = radiating_ideal_table();
    let base = run_with(&table, None, Ablation::new(50.0, T_VAPOR));
    let blown = run_with(&table, None, Ablation::new(50.0, T_VAPOR).with_blowing(1e6));
    assert_eq!(base.loss_conductive, 0.0, "no wall ⇒ no conduction");
    assert_eq!(blown.loss_conductive, 0.0, "no wall ⇒ no conduction, blown");
    assert_eq!(
        base.bounce.e_eff, blown.bounce.e_eff,
        "blowing changes nothing when conduction is off"
    );
}
