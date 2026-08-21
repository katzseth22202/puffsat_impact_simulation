"""Does f ~ 0.818 transfer from the 100 kg / 15 m heavy plate to the paper's 25 kg / 5 m plate?

Runs the repo's OWN contour construction (puffsat.contour) on both plates. The only change is
(mass, plate_radius); every other input -- geometry sweep, e_eff table, c_stag, P_limit, shape
box -- is untouched.
"""
import sys; sys.path.insert(0, "python")
import math
from puffsat import contour, heavyplate
from puffsat.analysis import impact_density

MASS_H, R_H = 100.0, 15.0   # heavy plate (design 12.1)
MASS_P, R_P = 25.0, 5.0     # paper's reference plate (design 2)
LOD_LO, LOD_HI = contour.L_OVER_D_BOX
RF_LO, RF_HI = contour.R_FOOT_BOX

rows = heavyplate.read_heavyplate_sweep()
vels = heavyplate.sweep_velocities(rows)

def best_shape(rho, mass, radius, d_over_d=contour.D_OVER_D_HEADLINE, samples=2048):
    """contour_point's inner loop, parameterized by plate. Returns (eta, L/D, rf/R) or None."""
    best = None
    for i in range(samples):
        rf = RF_LO + i * (RF_HI - RF_LO) / (samples - 1)
        lod = contour.l_over_d_for(rho, rf, mass, radius)
        if not (LOD_LO <= lod <= LOD_HI):
            continue
        eta = contour.eta_at(lod, rf, d_over_d)
        if best is None or eta > best[0]:
            best = (eta, lod, rf)
    return best

def rho_window(mass, radius):
    """Contour densities this plate's shape box can deliver at all."""
    lo = impact_density(LOD_HI, RF_HI, mass, radius)   # longest cloud, widest footprint
    hi = impact_density(LOD_LO, RF_LO, mass, radius)   # shortest cloud, tightest footprint
    return lo, hi

print(f"{'plate':>8}  {'rho_min':>9}  {'rho_max':>9}   (kg/m^3 the shape box can deliver)")
for tag, m, r in (("heavy", MASS_H, R_H), ("paper", MASS_P, R_P)):
    lo, hi = rho_window(m, r)
    print(f"{tag:>8}  {lo:9.4f}  {hi:9.4f}")

print()
hdr = (f"{'v':>5} {'rho':>7} | {'H rf/R':>6} {'H L/D':>6} {'H eta':>6} {'H L[m]':>6} {'H f':>6}"
       f" | {'P rf/R':>6} {'P L/D':>6} {'P eta':>6} {'P L[m]':>6} {'P f':>6} | {'df':>7}")
print(hdr); print("-" * len(hdr))

report = [16000, 22000, 28000, 34000, 40000, 45000, 50000, 55000, 60000, 63000]
for v in vels:
    if v not in report:
        continue
    c_stag = heavyplate.stagnation_coefficient_at_v(rows, v)
    e_at = heavyplate.e_eff_interpolator_at_v(rows, v)
    ceil = contour.rho_ceiling(v, c_stag, contour.P_LIMIT_BASELINE
                               if hasattr(contour, "P_LIMIT_BASELINE") else 400.0e6)
    out = {}
    for tag, m, r in (("H", MASS_H, R_H), ("P", MASS_P, R_P)):
        rho = min(ceil, impact_density(LOD_LO, RF_LO, m, r))   # min(ceiling, box max)
        b = best_shape(rho, m, r)
        if b is None:
            out[tag] = (rho, None)
            continue
        eta, lod, rf = b
        e_eff = e_at(rho)
        f = eta * (1.0 + e_eff) / 2.0
        L = 2.0 * lod * rf * r          # actual cloud length delivered
        out[tag] = (rho, (eta, lod, rf, L, f))
    (rho_h, H), (rho_p, P) = out["H"], out["P"]
    hs = (f"{H[2]:6.3f} {H[1]:6.3f} {H[0]:6.4f} {H[3]:6.2f} {H[4]:6.3f}"
          if H else f"{'--':>6} {'--':>6} {'--':>6} {'--':>6} {'--':>6}")
    ps = (f"{P[2]:6.3f} {P[1]:6.3f} {P[0]:6.4f} {P[3]:6.2f} {P[4]:6.3f}"
          if P else f"{'--':>6} {'--':>6} {'--':>6} {'--':>6} {'INFEAS':>6}")
    df = f"{P[4]-H[4]:+7.3f}" if (H and P) else f"{'--':>7}"
    print(f"{v/1000:5.0f} {rho_h:7.4f} | {hs} | {ps} | {df}")

# The velocity at which the paper's plate falls off its own contour.
lo_p, _ = rho_window(MASS_P, R_P)
c55 = heavyplate.stagnation_coefficient_at_v(rows, 55000.0)
k = 400.0e6 / c55                      # rho_ceiling = k / v^2
print(f"\npaper plate: dilutest deliverable cloud = {lo_p:.4f} kg/m^3")
print(f"             falls off the survivable contour above v = {math.sqrt(k/lo_p)/1000:.1f} km/s")
print(f"growth push band 45.58-56.53 km/s -> "
      f"{'INSIDE' if math.sqrt(k/lo_p) > 56530 else 'NOT fully inside'} the paper plate's envelope")
