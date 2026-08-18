//! Gray discrete-ordinates (Sₙ) transport on a **frozen** slab — the independent radiation model
//! that checks flux-limited diffusion where FLD is weakest (ADR-0012, design §12.1 Q9/Q21).
//!
//! FLD is exact when `τ ≫ 1` and degrades as the slab thins. The quantity this study actually reads
//! off the radiation model is the **escape flux to space**, and FLD reports it through its Marshak
//! surface as
//!
//! ```text
//! F_escape = (c/2)·E_surface,
//! ```
//!
//! which is the *half-range-isotropic* answer: it presumes the emergent intensity is the same in
//! every outward direction. Transport makes no such presumption and computes the first moment of
//! the actual angular distribution,
//!
//! ```text
//! F_escape = 2π ∫₀¹ I⁺(μ) μ dμ.
//! ```
//!
//! The two agree when `I⁺` really is isotropic (a thick emitter) and part company when it is not —
//! limb-brightened at small `τ`, forward-peaked across a transparent gap. **Measuring that gap is
//! the whole purpose of this module.**
//!
//! # Scope, and what this is not
//!
//! This is a **one-way diagnostic**: it reads a state the FLD run produced and reports what
//! transport would have said about the escaping flux. It does not feed back, so it yields a *bias
//! estimate*, not a corrected `e_eff`. Coupling would require replacing the radiation operator (a
//! coupled M1 or full transport solve), which is the escalation this diagnostic exists to decide
//! for or against.
//!
//! # Method
//!
//! **Short characteristics** with a piecewise-constant source. Along an ordinate `μ`, the gray
//! transfer equation on a frozen medium is
//!
//! ```text
//! μ dI/dx = χ (S − I),
//! ```
//!
//! whose exact solution across a cell of optical depth `Δτ = χ dx/|μ|` and constant `S` is
//!
//! ```text
//! I_out = I_in e^{−Δτ} + S (1 − e^{−Δτ}),
//! ```
//!
//! with cell average `Ī = S + (I_in − S)(1 − e^{−Δτ})/Δτ`. Composing that recurrence over the mesh
//! is *exact* for a uniform slab at any resolution, so the scheme carries no attenuation error of
//! its own — its errors are angular (finite ordinate count) and from the piecewise-constant source.
//!
//! Angles use **double Gauss–Legendre** (Sykes 1951): a Gauss–Legendre rule on `μ ∈ (0, 1]` and its
//! mirror on `[−1, 0)`. The half-range split matters — at a vacuum boundary `I(μ)` is discontinuous
//! across `μ = 0` (outgoing hemisphere filled, incoming empty), and a full-range rule straddling
//! that jump converges far more slowly.
//!
//! # Geometry convention
//!
//! Matches the 1D kernel: **cell `0` is at the wall**, the last cell faces **space**, and flux is
//! positive toward `+x` (toward space). So `Solution::escape_flux` is the flux the kernel's
//! `loss_escape_space` channel accumulates, and it is directly comparable to FLD's `(c/2)·E_last`.

use std::f64::consts::PI;

/// Discrete ordinates: direction cosines `μ` and quadrature weights over `μ ∈ [−1, 1]`.
///
/// The weights integrate to 2 (the full range), so an isotropic intensity `I₀` gives
/// `E = (2π/c)·Σ w I₀ = 4π I₀/c` and `F = 2π·Σ w μ I₀ = 0`, as it must.
#[derive(Debug, Clone)]
pub struct Ordinates {
    /// Direction cosines, ascending. Never zero (a grazing ordinate carries no net transport but
    /// would divide by zero in `Δτ = χ dx/|μ|`).
    pub mu: Vec<f64>,
    /// Quadrature weights, summing to 2.
    pub weight: Vec<f64>,
}

impl Ordinates {
    /// **Double Gauss–Legendre** with `half` ordinates per hemisphere (so `S_{2·half}`).
    ///
    /// A Gauss–Legendre rule on `[−1, 1]` is affinely mapped onto `(0, 1]` and mirrored, which
    /// resolves the `μ = 0` discontinuity of a vacuum boundary that a full-range rule smears
    /// across.
    ///
    /// # Panics
    /// Panics if `half` is zero.
    #[must_use]
    pub fn double_gauss(half: usize) -> Self {
        assert!(half > 0, "need at least one ordinate per hemisphere");
        let (nodes, weights) = gauss_legendre(half);
        let mut mu = Vec::with_capacity(2 * half);
        let mut weight = Vec::with_capacity(2 * half);
        // `nodes` descends, so negating the mapped values already ascends on the negative
        // hemisphere while the positive one has to be walked backwards. Net effect: the whole set
        // ascends in μ, −0.98 … −0.02, +0.02 … +0.98.
        for i in 0..half {
            mu.push(-(0.5 * nodes[i] + 0.5));
            weight.push(0.5 * weights[i]);
        }
        for i in (0..half).rev() {
            mu.push(0.5 * nodes[i] + 0.5);
            weight.push(0.5 * weights[i]);
        }
        debug_assert!(
            mu.windows(2).all(|w| w[0] < w[1]),
            "ordinates must be ascending in μ"
        );
        Self { mu, weight }
    }

    /// Number of ordinates (both hemispheres).
    #[must_use]
    pub fn len(&self) -> usize {
        self.mu.len()
    }

    /// Whether the set is empty. Never true for a constructed set; present because clippy asks for
    /// it alongside [`Ordinates::len`].
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.mu.is_empty()
    }
}

/// `n`-point Gauss–Legendre nodes and weights on `[−1, 1]`, by Newton iteration on the Legendre
/// polynomial `P_n` (the standard `gauleg`). Exact for polynomials up to degree `2n − 1`.
///
/// Nodes come out **descending** — the initial guess `cos(π(i + 3/4)/(n + 1/2))` walks the roots
/// from `+1` toward `−1`. [`Ordinates::double_gauss`] undoes that.
fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = vec![0.0; n];
    let mut weights = vec![0.0; n];
    let nf = n as f64;
    for (i, (node, weight)) in nodes.iter_mut().zip(weights.iter_mut()).enumerate() {
        // Chebyshev-like initial guess for the i-th root of P_n; converges in a handful of steps.
        let mut z = (PI * (i as f64 + 0.75) / (nf + 0.5)).cos();
        let mut dp = 0.0;
        for _ in 0..100 {
            // Bonnet recurrence: P_{j+1} = ((2j+1) z P_j − j P_{j−1})/(j+1).
            let (mut p_cur, mut p_prev) = (1.0, 0.0);
            for j in 0..n {
                let p_prev2 = p_prev;
                p_prev = p_cur;
                p_cur =
                    ((2.0 * j as f64 + 1.0) * z * p_prev - j as f64 * p_prev2) / (j as f64 + 1.0);
            }
            dp = nf * (z * p_cur - p_prev) / (z * z - 1.0);
            let step = p_cur / dp;
            z -= step;
            if step.abs() < 1e-15 {
                break;
            }
        }
        *node = z;
        *weight = 2.0 / ((1.0 - z * z) * dp * dp);
    }
    (nodes, weights)
}

/// The radiation entering the slab through one boundary.
#[derive(Debug, Clone, Copy)]
pub enum Incident<'a> {
    /// Nothing enters — a cold black absorber (the kernel's wall) or open space.
    Vacuum,
    /// The same intensity in every inward direction, e.g. a blackbody surface at `I = a c T⁴/(4π)`.
    Isotropic(f64),
    /// An arbitrary angular distribution, indexed by ordinate exactly as [`Ordinates::mu`] is.
    /// Entries for outgoing ordinates are ignored. Used to inject a collimated beam.
    PerOrdinate(&'a [f64]),
}

impl Incident<'_> {
    /// Incident intensity in ordinate `k`.
    fn at(&self, k: usize) -> f64 {
        match *self {
            Incident::Vacuum => 0.0,
            Incident::Isotropic(i) => i,
            Incident::PerOrdinate(v) => v[k],
        }
    }
}

/// A frozen slab for the transport sweep: `N` cells, wall-side first.
///
/// `chi` is the extinction coefficient **per length** (`χ = κ ρ`, matching
/// [`crate::radiation::Medium`]), and `source` is the source function `S` in intensity units. For
/// LTE gray absorption/emission with no scattering `S = B = a c T⁴/(4π)`, which
/// [`planck_source`] computes.
#[derive(Debug, Clone, Copy)]
pub struct Slab<'a> {
    /// Cell widths.
    pub dx: &'a [f64],
    /// Extinction coefficient per length `χ = κ ρ`.
    pub chi: &'a [f64],
    /// Source function `S` (intensity units).
    pub source: &'a [f64],
}

/// The gray LTE source function `S = B = a c T⁴/(4π)` — the frequency-integrated Planck intensity,
/// normalized so that an infinite medium in equilibrium has `E = a T⁴`.
#[must_use]
pub fn planck_source(a: f64, c: f64, temp: f64) -> f64 {
    a * c * temp.powi(4) / (4.0 * PI)
}

/// What one transport sweep produced.
#[derive(Debug, Clone)]
pub struct Solution {
    /// Cell-averaged radiation energy density `E = (2π/c)·Σ w Ī`, one per cell.
    pub e_rad: Vec<f64>,
    /// Net flux `F = 2π·Σ w μ I` at each of the `N+1` faces, positive toward space (`+x`).
    pub face_flux: Vec<f64>,
    /// Direction cosines of the outgoing (`μ > 0`) ordinates at the space-side face, ascending.
    pub mu_out: Vec<f64>,
    /// Emergent intensity `I⁺(μ)` at the space-side face, parallel to [`Solution::mu_out`]. This is
    /// the angular distribution FLD *assumes* to be flat; its shape is the diagnostic's substance.
    pub exit_intensity: Vec<f64>,
}

impl Solution {
    /// Net flux escaping the space-side face — the transport counterpart of the kernel's
    /// `loss_escape_space` integrand.
    ///
    /// # Panics
    /// Panics if the slab had no cells.
    #[must_use]
    pub fn escape_flux(&self) -> f64 {
        *self.face_flux.last().expect("slab has at least one face")
    }

    /// Eddington-like ratio `F/(cE)` at the space-side surface, using the last cell's energy
    /// density. FLD's Marshak boundary fixes this at exactly `1/2` by construction; transport lets
    /// it float, so the departure from `1/2` *is* the model error.
    ///
    /// # Panics
    /// Panics if the slab had no cells.
    #[must_use]
    pub fn surface_flux_ratio(&self, c: f64) -> f64 {
        let e = *self.e_rad.last().expect("slab has at least one cell");
        if e <= 0.0 {
            return 0.0;
        }
        self.escape_flux() / (c * e)
    }
}

/// `(1 − e^{−t})/t`, the short-characteristic cell-average weight, continued to `1` at `t = 0`.
fn cell_average_weight(t: f64) -> f64 {
    if t <= 0.0 { 1.0 } else { -(-t).exp_m1() / t }
}

/// Solve the gray transfer equation on a frozen slab by short characteristics.
///
/// `wall_incident` enters at the `x = 0` (cell 0) face travelling toward `+x`; `space_incident`
/// enters at the far face travelling toward `−x`.
///
/// # Panics
/// Panics unless `dx`, `chi` and `source` have the same non-zero length, and unless every
/// `PerOrdinate` slice is as long as `ordinates`.
#[must_use]
pub fn solve(
    slab: &Slab<'_>,
    ordinates: &Ordinates,
    c: f64,
    wall_incident: Incident<'_>,
    space_incident: Incident<'_>,
) -> Solution {
    let n = slab.dx.len();
    assert!(n > 0, "slab needs at least one cell");
    assert!(
        slab.chi.len() == n && slab.source.len() == n,
        "slab arrays must all have length {n}"
    );

    let mut e_rad = vec![0.0; n];
    let mut face_flux = vec![0.0; n + 1];
    let mut mu_out = Vec::new();
    let mut exit_intensity = Vec::new();

    for (k, (&mu, &w)) in ordinates.mu.iter().zip(&ordinates.weight).enumerate() {
        // The two hemispheres sweep in opposite directions; each visits every face and every cell.
        let forward = mu > 0.0;
        let mut intensity = if forward {
            wall_incident.at(k)
        } else {
            space_incident.at(k)
        };

        let entry_face = if forward { 0 } else { n };
        face_flux[entry_face] += 2.0 * PI * w * mu * intensity;

        for step in 0..n {
            let i = if forward { step } else { n - 1 - step };
            let d_tau = slab.chi[i] * slab.dx[i] / mu.abs();
            let attenuation = (-d_tau).exp();
            let s = slab.source[i];

            // Cell average of the exact in-cell solution, before advancing to the outgoing face.
            e_rad[i] += (2.0 * PI / c) * w * (s + (intensity - s) * cell_average_weight(d_tau));

            intensity = intensity * attenuation + s * (1.0 - attenuation);
            let exit_face = if forward { i + 1 } else { i };
            face_flux[exit_face] += 2.0 * PI * w * mu * intensity;
        }

        if forward {
            mu_out.push(mu);
            exit_intensity.push(intensity);
        }
    }

    Solution {
        e_rad,
        face_flux,
        mu_out,
        exit_intensity,
    }
}

/// The two radiation models' escape flux on one frozen state, side by side — the diagnostic's
/// output record.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EscapeComparison {
    /// Total optical depth `Σ χ dx` of the slab. The whole question is what happens near `1`.
    pub optical_depth: f64,
    /// What FLD reports: `(c/2)·E_surface` from its own evolved radiation field.
    pub fld_flux: f64,
    /// What transport reports: `2π ∫₀¹ I⁺(μ) μ dμ` from the formal solution.
    pub transport_flux: f64,
    /// FLD's surface radiation energy density, as handed in.
    pub fld_e_surface: f64,
    /// Transport's own surface radiation energy density, for separating a closure disagreement
    /// from the two models simply carrying different `E`.
    pub transport_e_surface: f64,
    /// Transport's `F/(cE)` at the surface. FLD pins this at `1/2`; the departure is the anisotropy
    /// FLD cannot see.
    pub surface_flux_ratio: f64,
    /// Optical depth `χ dx` of the **space-facing cell alone**.
    ///
    /// FLD's Marshak boundary reads a cell *average* and turns it into a surface flux, which is
    /// only meaningful while that cell is thin. Once `χ dx ≳ 1` the boundary is converting
    /// sub-photospheric gas into emergent radiation and its error is a discretization artifact, not
    /// a closure one — so any comparison has to report this number or it cannot tell the two apart.
    pub surface_cell_optical_depth: f64,
    /// `π·S` for the space-facing cell: the flux an *opaque* surface at that cell's temperature
    /// would emit.
    ///
    /// A yardstick, **not a ceiling**. It bounds the escape only when the space-facing cell is
    /// itself optically thick; when that cell is thin (as it is throughout the re-expansion here)
    /// the escaping radiation was emitted deeper and hotter and legitimately exceeds it. Read it
    /// together with [`EscapeComparison::surface_cell_optical_depth`]: the pair says how far the
    /// photosphere has retreated below the last cell.
    pub blackbody_flux_surface: f64,
}

impl EscapeComparison {
    /// Signed relative departure `(transport − FLD)/FLD`. Positive means FLD *under*-reports the
    /// escape, which would make the FLD bounce too elastic.
    ///
    /// Returns zero when FLD reports no escape at all, since a relative error on zero is not a
    /// number the caller can act on.
    #[must_use]
    pub fn relative_difference(&self) -> f64 {
        if self.fld_flux.abs() <= f64::MIN_POSITIVE {
            return 0.0;
        }
        (self.transport_flux - self.fld_flux) / self.fld_flux
    }
}

/// Run the transport solve on a frozen slab and pair its escape flux with FLD's.
///
/// `fld_e_surface` is the FLD run's radiation energy density in the space-facing cell; its Marshak
/// boundary turns that into `(c/2)·E` with no incident radiation (see
/// [`crate::radiation::RadBc::Marshak`]).
///
/// # The quasi-static assumption
///
/// Transport here is the *formal* (steady-state) solution on the instantaneous temperature and
/// opacity fields, while FLD's `E` carries time history. That is sound because the two clocks are
/// nowhere near each other: light crosses a ~10 m cloud in ~3·10⁻⁸ s against a bounce lasting
/// ~10⁻⁴ s, so the radiation field is quasi-static to a part in ~10³. It does mean the comparison
/// attributes *all* of any `E` mismatch to the models rather than to lag.
#[must_use]
pub fn compare_escape(
    slab: &Slab<'_>,
    ordinates: &Ordinates,
    c: f64,
    fld_e_surface: f64,
) -> EscapeComparison {
    let solution = solve(slab, ordinates, c, Incident::Vacuum, Incident::Vacuum);
    let optical_depth = slab
        .chi
        .iter()
        .zip(slab.dx)
        .map(|(&chi, &dx)| chi * dx)
        .sum();
    let last = slab.dx.len() - 1;
    EscapeComparison {
        optical_depth,
        fld_flux: 0.5 * c * fld_e_surface,
        transport_flux: solution.escape_flux(),
        fld_e_surface,
        transport_e_surface: *solution.e_rad.last().expect("slab has at least one cell"),
        surface_flux_ratio: solution.surface_flux_ratio(c),
        surface_cell_optical_depth: slab.chi[last] * slab.dx[last],
        blackbody_flux_surface: PI * slab.source[last],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// Below this relative error a quadrature has stopped resolving anything but double-precision
    /// roundoff, so "the next refinement must be better" is no longer a physical claim.
    const ROUNDOFF_FLOOR: f64 = 1e-12;

    /// A stretched mesh whose widths sum to `total`. Non-uniform on purpose: a uniform mesh would
    /// let a per-cell `Δτ` bug cancel itself across the sweep.
    fn stretched_mesh(n: usize, total: f64, ratio: f64) -> Vec<f64> {
        let raw: Vec<f64> = (0..n)
            .map(|i| ratio.powi(i32::try_from(i).unwrap()))
            .collect();
        let sum: f64 = raw.iter().sum();
        raw.iter().map(|w| w * total / sum).collect()
    }

    /// The angular quadrature underpins every moment the solver reports, so it is checked against
    /// its defining property before anything is built on it: an `n`-point Gauss–Legendre rule
    /// integrates polynomials of degree `2n − 1` **exactly**. Split across two hemispheres, the
    /// double-Gauss set of `2·half` ordinates is exact to degree `2·half − 1` on each half, so the
    /// moments `∫₋₁¹ μᵏ dμ` come out right for every `k` the rule claims.
    #[test]
    fn double_gauss_integrates_moments_exactly() {
        for &half in &[2_usize, 4, 8, 16] {
            let ord = Ordinates::double_gauss(half);
            assert_eq!(ord.len(), 2 * half);
            assert!(!ord.is_empty());
            assert!(ord.mu.windows(2).all(|w| w[0] < w[1]), "μ must ascend");

            for k in 0..(2 * half) {
                let numeric: f64 = ord
                    .mu
                    .iter()
                    .zip(&ord.weight)
                    .map(|(&mu, &w)| w * mu.powi(i32::try_from(k).unwrap()))
                    .sum();
                // ∫₋₁¹ μᵏ dμ = 2/(k+1) for even k, 0 for odd.
                let exact = if k % 2 == 0 {
                    2.0 / (k as f64 + 1.0)
                } else {
                    0.0
                };
                assert_relative_eq!(numeric, exact, epsilon = 1e-13);
            }
        }
    }

    /// **Acceptance 1 — the formal solution.** With nothing incident and a constant source
    /// function `S`, the intensity emerging from a slab of optical depth `τ` along direction `μ` is
    ///
    /// ```text
    /// I(μ) = S (1 − e^{−τ/μ}),
    /// ```
    ///
    /// the closed-form solution of `μ dI/dx = χ(S − I)`. The solver must reproduce it on a
    /// *stretched* mesh, which is the real content: it composes 200 separate exponentials and they
    /// must telescope to the single analytic one.
    #[test]
    fn uniform_slab_attenuates_to_the_formal_solution() {
        let s = 2.5;
        let n = 200;
        let dx = stretched_mesh(n, 1.0, 1.01);

        for &tau in &[0.02, 0.5, 2.0, 10.0] {
            let chi = vec![tau; n]; // slab is 1.0 long, so total optical depth is τ
            let source = vec![s; n];
            let slab = Slab {
                dx: &dx,
                chi: &chi,
                source: &source,
            };
            let ord = Ordinates::double_gauss(8);
            let sol = solve(&slab, &ord, 1.0, Incident::Vacuum, Incident::Vacuum);

            for (&mu, &i_out) in sol.mu_out.iter().zip(&sol.exit_intensity) {
                let expected = s * (1.0 - (-tau / mu).exp());
                assert_relative_eq!(i_out, expected, max_relative = 1e-12);
            }
        }
    }

    /// **Acceptance 1b — the sweep is direction-aware.** A uniform slab cannot detect a reversed
    /// sweep: with one constant source the cell visit order does not change the answer, and a
    /// deliberate mutation swapping the backward sweep for a forward one passes acceptance 1
    /// untouched. So the order has to be pinned by an *asymmetric* slab.
    ///
    /// Two layers — an emitting one at the wall (`S`, `τ_a`) and a cold absorbing blanket facing
    /// space (`S = 0`, `τ_b`). The closed-form emergent intensities differ by direction:
    ///
    /// ```text
    /// space side (μ > 0): I = S(1 − e^{−τ_a/μ})·e^{−τ_b/μ}   emitted, then blanketed
    /// wall side  (μ < 0): I = S(1 − e^{−τ_a/|μ|})            the blanket adds nothing
    /// ```
    ///
    /// Face fluxes are compared against these analytic intensities summed on the solver's own
    /// ordinates, which isolates the sweep from the angular discretization: the quadrature is the
    /// same on both sides, so any disagreement is the transport, not the rule.
    #[test]
    fn layered_slab_emerges_differently_at_each_face() {
        let (s, tau_a, tau_b) = (3.0, 0.8, 1.5);
        let half = 24;
        let dx = [
            stretched_mesh(half, 1.0, 1.02),
            stretched_mesh(half, 1.0, 0.98),
        ]
        .concat();
        let chi: Vec<f64> = (0..2 * half)
            .map(|i| if i < half { tau_a } else { tau_b })
            .collect();
        let source: Vec<f64> = (0..2 * half)
            .map(|i| if i < half { s } else { 0.0 })
            .collect();
        let slab = Slab {
            dx: &dx,
            chi: &chi,
            source: &source,
        };
        let ord = Ordinates::double_gauss(8);
        let sol = solve(&slab, &ord, 1.0, Incident::Vacuum, Incident::Vacuum);

        // Space side: emitted by layer A, then attenuated through the cold layer B.
        for (&mu, &i_out) in sol.mu_out.iter().zip(&sol.exit_intensity) {
            let expected = s * (1.0 - (-tau_a / mu).exp()) * (-tau_b / mu).exp();
            assert_relative_eq!(i_out, expected, max_relative = 1e-12);
        }

        // Wall side: layer B is transparent to its own (absent) emission, so only A contributes.
        let expected_wall: f64 = ord
            .mu
            .iter()
            .zip(&ord.weight)
            .filter(|&(&mu, _)| mu < 0.0)
            .map(|(&mu, &w)| 2.0 * PI * w * mu * s * (1.0 - (-tau_a / mu.abs()).exp()))
            .sum();
        assert_relative_eq!(sol.face_flux[0], expected_wall, max_relative = 1e-12);
        assert!(
            sol.face_flux[0] < 0.0,
            "wall-bound flux must point toward −x, got {}",
            sol.face_flux[0]
        );

        // And the blanket really is doing something — otherwise the direction check is vacuous.
        let space_side: f64 = sol.escape_flux();
        assert!(
            space_side < 0.5 * expected_wall.abs(),
            "cold blanket should cut the escaping flux well below the wall-side flux"
        );
    }

    /// **Acceptance 2 — free streaming, and the one thing FLD structurally cannot do.**
    ///
    /// FLD's flux is `F = −D ∇E`: *proportional to the gradient*, whatever the limiter does to `D`.
    /// So a **uniform** radiation field carries **zero** flux in FLD — at any opacity, however
    /// transparent. That is not a tuning problem, it is the form of the closure.
    ///
    /// Transport has no such constraint. A collimated beam crossing a transparent gap has a flat
    /// `E` profile and still carries the maximum flux the field can carry: `F = μ·cE`, which tends
    /// to `cE` as the beam approaches the normal. Both halves are asserted here — the transport
    /// beam, and FLD sitting still on the identical state.
    #[test]
    fn a_transparent_beam_streams_at_ce_where_fld_carries_nothing() {
        use crate::radiation::{Limiter, Medium, RadBc, RadConstants, fld_substep};

        let c = 3.0;
        let n = 40;
        let dx = vec![0.25; n];
        let chi = vec![0.0; n]; // vacuum gap
        let source = vec![0.0; n];
        let slab = Slab {
            dx: &dx,
            chi: &chi,
            source: &source,
        };
        let ord = Ordinates::double_gauss(8);

        // A unit beam in the single most-forward ordinate, nothing in any other direction.
        let mut beam = vec![0.0; ord.len()];
        let k_beam = ord.len() - 1;
        beam[k_beam] = 1.0;
        let mu_beam = ord.mu[k_beam];
        let sol = solve(
            &slab,
            &ord,
            c,
            Incident::PerOrdinate(&beam),
            Incident::Vacuum,
        );

        // Crosses the gap untouched.
        assert_relative_eq!(
            *sol.exit_intensity.last().unwrap(),
            1.0,
            max_relative = 1e-14
        );
        for (&mu, &i_out) in sol.mu_out.iter().zip(&sol.exit_intensity) {
            if mu < mu_beam {
                assert_relative_eq!(i_out, 0.0, epsilon = 1e-14);
            }
        }

        // Flux is the full `μ·cE` at every face, and `E` is flat — a gradient-free maximal flux.
        for i in 0..n {
            assert_relative_eq!(sol.e_rad[i], sol.e_rad[0], max_relative = 1e-14);
            assert_relative_eq!(
                sol.face_flux[i],
                mu_beam * c * sol.e_rad[0],
                max_relative = 1e-14
            );
        }
        // The quadrature's most-forward ordinate approaches normal incidence as S_N refines, so the
        // achievable ratio approaches the free-streaming ceiling F = cE.
        assert!(
            mu_beam > 0.98,
            "S16's outermost ordinate should be near-normal, got μ={mu_beam}"
        );

        // Now FLD on the very same uniform field: whatever the opacity, ∇E = 0 ⇒ no transport.
        for &chi_r in &[0.0, 1e-3, 1.0, 100.0] {
            let mut e_fld = vec![sol.e_rad[0]; n];
            let before = e_fld.clone();
            let medium = Medium {
                dx: &dx,
                center_spacing: &vec![0.25; n - 1],
                temp: &vec![0.0; n],
                cv_vol: &vec![1.0; n],
                chi_planck: &vec![0.0; n], // matter decoupled: isolate the transport term
                chi_ross: &vec![chi_r; n],
                source: None,
            };
            let _ = fld_substep(
                &medium,
                &mut e_fld,
                RadBc::Reflecting,
                RadBc::Reflecting,
                1e-3,
                RadConstants { c, a: 1.0 },
                Limiter::LevermorePomraning,
            );
            for i in 0..n {
                assert_relative_eq!(e_fld[i], before[i], max_relative = 1e-12);
            }
        }
    }

    /// `E₁(x) = ∫₁^∞ t⁻¹ e^{−xt} dt`, by the textbook route: the series when `x < 1`, a
    /// modified-Lentz continued fraction when `x ≥ 1`.
    ///
    /// Deliberately a *different algorithm* from anything in the solver, so the `Eₙ` oracles below
    /// are independent truth rather than a restatement of the code.
    fn exp_integral_1(x: f64) -> f64 {
        assert!(x > 0.0, "E₁ diverges at the origin");
        if x < 1.0 {
            // E₁(x) = −γ − ln x + Σ_{k≥1} (−1)^{k+1} xᵏ/(k·k!)
            const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;
            let mut sum = 0.0;
            let mut term = 1.0;
            for k in 1..60 {
                let kf = f64::from(k);
                term *= -x / kf;
                sum -= term / kf;
            }
            -EULER_GAMMA - x.ln() + sum
        } else {
            // Continued fraction E₁(x) = e^{−x}/(x + 1 − 1²/(x + 3 − 2²/(x + 5 − …))).
            let tiny = 1e-300;
            let (mut b, mut c_l) = (x + 1.0, 1.0 / tiny);
            let mut d = 1.0 / b;
            let mut h = d;
            for i in 1..200 {
                let a = -f64::from(i) * f64::from(i);
                b += 2.0;
                d = 1.0 / (a * d + b);
                c_l = b + a / c_l;
                let del = c_l * d;
                h *= del;
                if (del - 1.0).abs() < 1e-16 {
                    break;
                }
            }
            h * (-x).exp()
        }
    }

    /// `E₂(x) = ∫₀¹ e^{−x/μ} dμ`, via the recurrence `E₂ = e^{−x} − x E₁`. `E₂(0) = 1`.
    fn exp_integral_2(x: f64) -> f64 {
        if x == 0.0 {
            return 1.0;
        }
        (-x).exp() - x * exp_integral_1(x)
    }

    /// `E₃(x) = ∫₀¹ μ e^{−x/μ} dμ`, via `E₃ = (e^{−x} − x E₂)/2`. `E₃(0) = 1/2`.
    fn exp_integral_3(x: f64) -> f64 {
        if x == 0.0 {
            return 0.5;
        }
        ((-x).exp() - x * exp_integral_2(x)) / 2.0
    }

    /// The oracle has to be right before it can judge anything, so it is pinned against published
    /// values of `E₁` (Abramowitz & Stegun table 5.1) carried through the recurrence, plus the two
    /// endpoints `E₃(0) = 1/2` and `E₃(∞) = 0`.
    ///
    /// Tolerance is `1e-8`, not machine epsilon: the table's `E₁` values are quoted to ten
    /// significant figures, and at `x = 5` the recurrence `E₃ = (e^{−x} − x(e^{−x} − x E₁))/2`
    /// amplifies that truncation to ~1e-9 relative. The limit is the published constant, not the
    /// algorithm.
    #[test]
    fn exponential_integral_oracle_matches_published_values() {
        for &(x, e1) in &[
            (0.5_f64, 0.559_773_594_776_16_f64),
            (1.0, 0.219_383_934_395_52),
            (2.0, 0.048_900_510_708_06),
            (5.0, 1.148_295_591_2e-3),
        ] {
            let e2 = (-x).exp() - x * e1;
            let expected = ((-x).exp() - x * e2) / 2.0;
            assert_relative_eq!(exp_integral_3(x), expected, max_relative = 1e-8);
        }
        assert_relative_eq!(exp_integral_3(0.0), 0.5, max_relative = 1e-14);
        assert!(exp_integral_3(40.0) < 1e-18);
        assert_relative_eq!(exp_integral_2(0.0), 1.0, max_relative = 1e-14);
    }

    /// **Acceptance 5 — the energy density, not just the flux.**
    ///
    /// Added because a mutation survived the rest of the ladder: swapping the cell-average weight
    /// `(1 − e^{−Δτ})/Δτ` for a plain `e^{−Δτ}` left every other test green. Nothing pinned
    /// `e_rad` in an optically active medium — the flux tests read face intensities, and the beam
    /// test runs at `Δτ = 0`, where the two forms coincide. `e_rad` is what
    /// [`Solution::surface_flux_ratio`] divides by and what the FLD comparison is calibrated
    /// against, so leaving it unpinned was not acceptable.
    ///
    /// The truth is analytic. A uniform slab of total optical depth `τ_L` with no incident
    /// radiation sees, at optical depth `t` from the wall face,
    ///
    /// ```text
    /// E(t) = (2πS/c)·[2 − E₂(t) − E₂(τ_L − t)],
    /// ```
    ///
    /// one `E₂` per hemisphere. Deep inside a thick slab both vanish and `E → 4πS/c = aT⁴` —
    /// radiative equilibrium, which is the check that the normalization is right.
    #[test]
    fn energy_density_matches_the_analytic_two_hemisphere_profile() {
        let (s, c, tau_l) = (1.3, 2.0, 3.0);
        let n = 400;
        let dx = vec![1.0 / n as f64; n];
        let chi = vec![tau_l; n]; // unit-length slab, so χ = τ_L
        let source = vec![s; n];
        let sol = solve(
            &Slab {
                dx: &dx,
                chi: &chi,
                source: &source,
            },
            &Ordinates::double_gauss(24),
            c,
            Incident::Vacuum,
            Incident::Vacuum,
        );

        for (i, &e) in sol.e_rad.iter().enumerate() {
            let t = tau_l * (i as f64 + 0.5) / n as f64;
            let exact = (2.0 * PI * s / c) * (2.0 - exp_integral_2(t) - exp_integral_2(tau_l - t));
            assert_relative_eq!(e, exact, max_relative = 2e-3);
        }

        // Radiative equilibrium deep inside a thick slab: E = 4πS/c.
        let thick = vec![80.0; n];
        let deep = solve(
            &Slab {
                dx: &dx,
                chi: &thick,
                source: &source,
            },
            &Ordinates::double_gauss(24),
            c,
            Incident::Vacuum,
            Incident::Vacuum,
        );
        assert_relative_eq!(deep.e_rad[n / 2], 4.0 * PI * s / c, max_relative = 1e-9);
    }

    /// **Acceptance 6 — the comparison record itself**, and the thick-limit anchor that has to hold
    /// before any `τ ~ 1` disagreement can be read as a closure error.
    ///
    /// A thick isothermal slab radiates the blackbody flux `σT⁴ = acT⁴/4` — an answer fixed by
    /// thermodynamics, not by either model. Transport must produce it. FLD's Marshak rule
    /// `F = (c/2)E` then reproduces the same number **iff** its surface cell has drained to
    /// `E = aT⁴/2`, half the equilibrium `aT⁴`, which is exactly the known behaviour of that
    /// boundary condition. So the thick limit is where the two models are *supposed* to agree, and
    /// the record must show them agreeing.
    #[test]
    fn thick_isothermal_slab_radiates_the_blackbody_flux_from_both_models() {
        let (a, c, temp) = (2.0, 3.0, 1.7);
        let n = 200;
        let dx = vec![0.2; n];
        let chi = vec![1.0; n]; // τ_total = 40
        let s = planck_source(a, c, temp);
        let source = vec![s; n];
        let slab = Slab {
            dx: &dx,
            chi: &chi,
            source: &source,
        };

        // Stefan–Boltzmann: σT⁴ = acT⁴/4, the flux a black surface emits.
        let sigma_t4 = a * c * temp.powi(4) / 4.0;
        // The surface energy density that makes FLD's (c/2)E equal that flux.
        let drained = 2.0 * sigma_t4 / c;
        assert_relative_eq!(drained, 0.5 * a * temp.powi(4), max_relative = 1e-14);

        let cmp = compare_escape(&slab, &Ordinates::double_gauss(16), c, drained);

        assert_relative_eq!(cmp.optical_depth, 40.0, max_relative = 1e-14);
        assert_relative_eq!(cmp.transport_flux, sigma_t4, max_relative = 1e-6);
        assert_relative_eq!(cmp.fld_flux, 0.5 * c * drained, max_relative = 1e-14);
        assert!(
            cmp.relative_difference().abs() < 1e-6,
            "thick limit must agree, got {:.3e}",
            cmp.relative_difference()
        );

        // Had FLD's surface cell *not* drained — sitting at full equilibrium aT⁴ — its Marshak rule
        // would report twice the blackbody flux, and the record would say so.
        let undrained = compare_escape(&slab, &Ordinates::double_gauss(16), c, a * temp.powi(4));
        assert_relative_eq!(undrained.relative_difference(), -0.5, max_relative = 1e-6);
    }

    /// **Acceptance 7 — the anisotropy FLD cannot see.** FLD's Marshak rule fixes `F/(cE) = 1/2`
    /// identically, for every state it will ever be handed. Transport lets the ratio follow the
    /// actual angular distribution, and in a *thin* slab that distribution is limb-brightened
    /// (`I⁺ ∝ 1/μ`), not flat.
    ///
    /// The ratio is formed on cell-averaged `E`, deliberately: FLD's boundary reads its last
    /// **cell**, so the comparison has to read the same thing to be about the closure rather than
    /// about where each model samples. The price is that the ratio — unlike the flux, which is an
    /// exact face quantity — is **mesh-sensitive**: a surface cell spanning `Δτ = 0.25` averages in
    /// gas well below the photosphere and reads ~0.38 even in the thick limit. Pinned below, so
    /// nobody reads a coarse-mesh ratio as evidence of anisotropy.
    #[test]
    fn surface_flux_ratio_departs_from_one_half_as_the_slab_thins() {
        let (a, c, temp) = (1.0, 1.0, 1.0);
        let s = planck_source(a, c, temp);

        // Resolve the surface cell to Δτ ≤ 0.01 at every τ, so the ratio is reading the
        // photosphere rather than the mesh.
        let run = |tau: f64, n: usize| {
            let dx = vec![1.0 / n as f64; n];
            let chi = vec![tau; n];
            let source = vec![s; n];
            compare_escape(
                &Slab {
                    dx: &dx,
                    chi: &chi,
                    source: &source,
                },
                &Ordinates::double_gauss(24),
                c,
                0.5 * a * temp.powi(4),
            )
        };

        // Cell counts chosen so the surface cell spans Δτ ≈ 0.01 at every τ.
        let mut ratios = Vec::new();
        for &(tau, n) in &[
            (0.01_f64, 400_usize),
            (0.1, 400),
            (1.0, 400),
            (10.0, 1_000),
            (100.0, 10_000),
        ] {
            ratios.push((tau, run(tau, n).surface_flux_ratio));
        }

        // Thick: the emergent intensity really is near-isotropic, so transport approaches FLD's
        // 1/2. Only *approaches*: the residual converges like `τ_cell·ln τ_cell`, so Δτ = 0.01
        // still leaves ~3%. Chasing 1% would take Δτ ~ 5·10⁻⁴ and a 200 000-cell mesh — which is
        // the practical reason the production diagnostic reads the **flux**, an exact face
        // quantity, and treats this ratio as a qualitative anisotropy indicator only.
        let (_, thick) = *ratios.last().unwrap();
        assert_relative_eq!(thick, 0.5, max_relative = 0.05);

        // Thin: it does not. The ratio moves monotonically away as τ falls, and by τ = 0.01 the
        // disagreement with FLD's fixed 1/2 is large enough to matter to `e_eff`.
        assert!(
            ratios.windows(2).all(|w| w[0].1 < w[1].1),
            "ratio should rise monotonically with τ toward 1/2, got {ratios:?}"
        );
        let (_, thin) = ratios[0];
        assert!(
            (thin - 0.5).abs() > 0.1,
            "τ = 0.01 should be far from FLD's fixed 1/2, got {thin:.4}"
        );

        // The mesh sensitivity itself, pinned: the same thick slab read through a coarse surface
        // cell (Δτ = 0.25) reports ~0.38 — an artifact of averaging, not anisotropy. The escaping
        // *flux* is a face quantity and does not move.
        let coarse = run(100.0, 400);
        let fine = run(100.0, 10_000);
        assert!(
            coarse.surface_flux_ratio < 0.42,
            "coarse surface cell should read low, got {:.4}",
            coarse.surface_flux_ratio
        );
        assert_relative_eq!(
            coarse.transport_flux,
            fine.transport_flux,
            max_relative = 1e-9
        );
    }

    /// **Acceptance 3 — angular convergence of the escape flux**, the quantity the diagnostic
    /// exists to report.
    ///
    /// A uniform slab has a closed-form emergent flux: integrating the formal solution over the
    /// outgoing hemisphere,
    ///
    /// ```text
    /// F = 2π ∫₀¹ S(1 − e^{−τ/μ}) μ dμ = 2πS [1/2 − E₃(τ)],
    /// ```
    ///
    /// which tends to `πS = σT⁴` as `τ → ∞` — the blackbody surface flux, the answer we already
    /// know must come out. `Sₙ` must converge to it as the ordinate count rises.
    ///
    /// **Convergence is much slower at small `τ`, and that is physics, not a defect.** The emergent
    /// intensity `S(1 − e^{−τ/μ}) ≈ Sτ/μ` is limb-brightened — it diverges as `μ → 0` — so the
    /// integrand a thin slab hands the quadrature is nearly singular. The required ordinate count
    /// is therefore a function of `τ`, which is exactly why the production diagnostic reports its
    /// angular convergence rather than fixing `N` once.
    #[test]
    fn escape_flux_converges_on_the_exponential_integral() {
        let s = 1.7;
        for &tau in &[0.1, 1.0, 5.0, 30.0] {
            let exact = 2.0 * PI * s * (0.5 - exp_integral_3(tau));

            let mut previous_error = f64::INFINITY;
            for &half in &[2_usize, 4, 8, 16, 32] {
                let n = 64;
                let dx = vec![1.0 / n as f64; n];
                let chi = vec![tau; n];
                let source = vec![s; n];
                let slab = Slab {
                    dx: &dx,
                    chi: &chi,
                    source: &source,
                };
                let sol = solve(
                    &slab,
                    &Ordinates::double_gauss(half),
                    1.0,
                    Incident::Vacuum,
                    Incident::Vacuum,
                );
                let error = (sol.escape_flux() - exact).abs() / exact;
                // Thick slabs converge to roundoff well before S64, where "still improving" stops
                // being a meaningful demand — the floor is double precision, not the quadrature.
                assert!(
                    error < previous_error || error < ROUNDOFF_FLOOR,
                    "τ={tau}: S{} error {error:.3e} did not improve on {previous_error:.3e}",
                    2 * half
                );
                previous_error = error.min(previous_error);
            }
            // S64 resolves even the limb-brightened thin case to well under a percent.
            assert!(
                previous_error < 5e-3,
                "τ={tau}: S64 escape flux still off by {previous_error:.3e}"
            );
        }

        // The τ ≫ 1 limit is the blackbody surface flux πS — an answer known without any solver.
        let n = 64;
        let dx = vec![1.0 / n as f64; n];
        let chi = vec![60.0; n];
        let source = vec![s; n];
        let sol = solve(
            &Slab {
                dx: &dx,
                chi: &chi,
                source: &source,
            },
            &Ordinates::double_gauss(16),
            1.0,
            Incident::Vacuum,
            Incident::Vacuum,
        );
        assert_relative_eq!(sol.escape_flux(), PI * s, max_relative = 1e-6);
    }

    /// **Acceptance 4 — the diffusion limit.** Deep inside an optically thick medium with a source
    /// function linear in `x`, transport theory gives the Fickian flux exactly (the `K = E/3`
    /// Eddington closure is exact for a linear source in an infinite medium):
    ///
    /// ```text
    /// F = −(c/3χ)·dE/dx = −(4π/3χ)·dS/dx,   using E = 4πS/c.
    /// ```
    ///
    /// This is where FLD is supposed to be right, so agreement here is the *precondition* for
    /// trusting the τ ~ 1 comparison — if the two models disagreed in the thick limit the
    /// diagnostic would be measuring a bug, not a closure error.
    ///
    /// Also an order-of-accuracy check: the piecewise-constant source makes the scheme first-order
    /// in cell optical depth, so halving `Δτ` must roughly halve the error.
    #[test]
    fn thick_slab_reproduces_the_fickian_flux_at_first_order() {
        let (chi, slope, s0) = (1.0, 0.5, 4.0);
        let length = 40.0; // τ_total = 40: the boundary layers are e^{−20} away from midplane
        let exact = -(4.0 * PI / (3.0 * chi)) * slope;

        let mut previous_error = f64::INFINITY;
        let mut errors = Vec::new();
        // Refine by cell count rather than by Δτ, so the mesh size is exact rather than rounded
        // off a float. τ_total = 40, so these are Δτ = 0.4, 0.2, 0.1, 0.05.
        for &n in &[100_usize, 200, 400, 800] {
            let d_tau = length * chi / n as f64;
            let dx = vec![length / n as f64; n];
            let chi_v = vec![chi; n];
            // S(x) = s0 + slope·x, sampled at cell centers.
            let source: Vec<f64> = (0..n)
                .map(|i| s0 + slope * (i as f64 + 0.5) * length / n as f64)
                .collect();
            let sol = solve(
                &Slab {
                    dx: &dx,
                    chi: &chi_v,
                    source: &source,
                },
                &Ordinates::double_gauss(16),
                1.0,
                Incident::Vacuum,
                Incident::Vacuum,
            );

            let mid = n / 2;
            let error = (sol.face_flux[mid] - exact).abs() / exact.abs();
            errors.push((d_tau, error));
            assert!(
                error < previous_error,
                "Δτ={d_tau}: error {error:.3e} did not improve on {previous_error:.3e}"
            );
            previous_error = error;
        }

        // First order: each halving of Δτ must cut the error by at least ~1.7x.
        for pair in errors.windows(2) {
            let ratio = pair[0].1 / pair[1].1;
            assert!(
                ratio > 1.7,
                "expected ~first-order convergence, got {ratio:.2} between Δτ={} and {}",
                pair[0].0,
                pair[1].0
            );
        }
        assert!(
            previous_error < 0.02,
            "finest mesh should sit within 2% of the Fickian flux, got {previous_error:.3e}"
        );
    }
}
