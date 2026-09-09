# Rate control and local fixed-energy temperature feedback

## Scope and retained assumptions

Continue Steps 7–8 on their saved 100–250 microsecond **equilibrium** flow.
No new rate coefficients, nonlinear chemistry integration, hydrodynamic
trajectory or wall impulse are produced. The rate set remains ADR-0046's
ten reversible neutral reactions with study-matched reverse constants.
The pressure-dependent single-channel curve of ADR-0047 is not substituted
where it lacked coverage. All previous high-temperature/dense-fluid rate
limitations remain; all artifacts retain `kinetics_validated=false`.

Evaluate both fixed-temperature and **fixed-volume, fixed-total-energy**
local linear response. Ions and electrons remain chemically fixed; their
sensible heat capacity is retained. Neither constraint implies a parcel
actually remains at fixed density during its evolving loading history.
The history only supplies local states, weights and comparison clocks.

## Energy-consistent local temperature response

Let Q span the four-dimensional reactive subspace in the sqrt(n)-scaled
coordinates, and let `y=Q^T delta(n)/sqrt(n)`. B is ADR-0046's symmetric
positive reactive relaxation matrix. For per-neutral **total internal
energies** u (ground energy plus translational, rotational and vibrational
energy), the fixed-energy constraint gives

```
a = Q^T diag(sqrt(n)) u
delta(T) = -a^T y / Cv
H = I + a a^T / (k_B T² Cv)
dy/dt = -B H y.
```

Here Cv is the **volumetric fixed-composition heat capacity**, not the
equilibrium heat capacity. It includes neutral internal modes, sensible
heat of the fixed ions/electrons, and the composition-independent dense
residual `rho*cv_ref(rho)*(1000/T)^2`. Using equilibrium Cv would count the
chemical response twice. The residual does not add a composition affinity
at fixed mass/density. Concentration-form van't Hoff is
`d ln K_r / dT = nu_r^T u / (k_B T²)`, verified against the existing
partition functions. Elemental energy gauges and the common storage shift
do not change the reactive energy constraint.

Use the Step 6 log-density cubic-Hermite continuation for `cv_ref`, with
512 knots over 1e-5–2000 kg/m3 and slopes `-1000*p_res,TT/rho`. Reject
out-of-curve states and nonpositive total frozen Cv. Parent metadata must
specify the 1000 K, T^-2 residual closure. Source-point tests and a
1024-knot history rerun check this interpolation, not the physical EOS.

For the ground/bond-energy observable `b=Q^T diag(sqrt(n)) e0`, the
fixed-energy static covariance is proportional to `H^-1`. Therefore
symmetrize and transform the observable **together**:

```
C = sqrt(H) B sqrt(H)
g = H^(-1/2) b
tau_isoE = g^T C^(-1) g / (g^T g).
```

This is the normalized integrated bond-energy autocorrelation response
under a small conjugate chemical-energy bias, with conserved local total
energy. It does not describe every possible initial composition
perturbation. Reusing the iso-T energy weights with new eigenvalues would
not be the same observable/constraint. An infinite heat-capacity bath
recovers the earlier isothermal diagnostic exactly.

A shorter iso-E time is **not more recombination energy or useful work**:
the available chemical susceptibility changes too, and reaction-induced
temperature changes alter the restoring chemical affinity. This local
relaxation is neither nonlinear thermal tracking of moving equilibrium
nor a finite-rate thrust calculation.

## Which reaction rates control the time?

For a given constraint, write `C=sum q_r d_r d_r^T`, with `d_r=Q^T nu_r/sqrt(n)`
for iso-T or its sqrt(H)-transformed counterpart for iso-E. Then

```
S_r = -d ln(tau) / d ln(q_r)
    = q_r [d_r^T C^(-1) g]^2 / [g^T C^(-1) g].
```

The nonnegative sensitivities sum to one for resolved spectra. They are
local logarithmic **rate-control** weights, not fractions of reaction
heat, bond energy returned or plate impulse. Summaries average them using
the same positive consecutive parcel bond-store decreases as Step 7,
conditioned on resolved rates with coverage stated explicitly.

Evaluate baseline and zero-direct-water-channel sensitivities. Also slow
**all six association channels** (source reactions 1, 2, 12, 13, 14, 15)
together by 1e3 or 1e6, holding all four exchange rates fixed. Each factor
scales both reaction directions; equilibrium composition does not change.
These factors are deliberately parameterized tests, not evaluated rate
uncertainties. The generic, H2-collider and H2O-collider H2-formation
channels remain separate and are never double-counted.

Deleting all associations exactly leaves an extra conserved neutral
particle count: exchange reactions alone span only three reactive
directions. Do not floor this zero eigenvalue into a fast or finite full
equilibration time. This differs from removing direct water association
alone, where the other reactions still have full reactive rank.

## Checks and artifact contract

Retain the Step 7 weighting, exact time windows, parent conservation/mass
gates and min(expansion, sampled store-change) clock. All missing clocks,
tiny changes and unresolved spectra remain in summary denominators.
Conditional sensitivity means are normalized only over resolved rates;
their coverage is reported separately from clock coverage.
Accept a rate variant only when **both** iso-T and iso-E spectra and
sensitivity sums pass, so the paired comparison uses a common subset.
The spectral condition ceiling is 1e12 and the sensitivity-sum error
ceiling is 1e-7. In the million-fold slowdown stress case, marginal
numerical acceptance can change under tiny heat-capacity perturbations;
report that coverage/sampling effect, not just the physical-cv difference.

Acceptance tests cover van't Hoff, frozen-Cv differentiation, iso-T and
infinite-bath limits, finite-difference logarithmic sensitivities, rate
rank loss, elemental energy-gauge invariance, residual refinement, and
an independent finite-difference Jacobian whose temperature is solved
from **conserved total energy** after composition perturbations. No new
hydro/EOS production closure is introduced.

`make water-plate-thermal-kinetics-audit` retains separate CSV/JSON evidence.
It accepts the same parent/stride Makefile overrides as the preceding
audits. For the residual check run `thermal_kinetics_audit` with N=80,
`--strides 2 --residual-knots 1024 --output
data/results/water_plate/thermal_kinetics_refined.csv`.
