//! Grid convergence of the reflected wall spike in the 1-D spherical chamber (methane, γ 1.2).
use euler2d::sphere::Sphere1D;

fn main() {
    let (radius, energy, mass, fill_rho) = (1.684, 0.5 * 2.5 * 75_000.0_f64.powi(2), 6.9, 3.44);
    for cells in [842_usize, 1684, 3368, 6736] {
        let mut s = Sphere1D::new(cells, radius, 1.2, |_| (fill_rho, 0.0, 5.0e5));
        s.deposit(energy, mass, 0.15);
        let mut peak: f64 = 0.0;
        let mut qsp = (0.0, 0_u32);
        s.run_to(0.0, 4.0e-3, |t, p| {
            peak = peak.max(p);
            if t > 3.0e-3 {
                qsp = (qsp.0 + p, qsp.1 + 1);
            }
        });
        let mean = qsp.0 / f64::from(qsp.1);
        println!(
            "cell {:.3} mm: peak {:.3} GPa, QSP {:.1} MPa, peak/QSP {:.1}",
            radius / cells as f64 * 1e3,
            peak / 1e9,
            mean / 1e6,
            peak / mean
        );
    }
}
