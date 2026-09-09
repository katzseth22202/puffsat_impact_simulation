# Nonlinear neutral chemistry on prescribed parcel histories

## Scope and accounting

Integrate the six-neutral network in Rust, with fixed ion/electron populations
per unit mass. Python prepares thermochemistry and saved-flow density histories;
no Python time integration and no changed hydrodynamic trajectory. Start all
parcels at the same 153.5 microsecond switch, report net bond-store return at
193 and 250 microseconds. These are loading-window diagnostics, not a universal
pressure-decoupling time or a thrust attribution.

Use populations per original H2O formula unit, y=n*M_H2O/rho. Two elemental
budgets eliminate atomic H and O; the four molecular populations and temperature
are the independent variables. Chemical energy changes exchange with sensible
energy. Include fixed-ion/electron sensible heat and the existing composition-
independent dense Helmholtz residual. Match the initial parent energy/pressure
using the same constant energy and gas-constant offsets as the frozen restarts.
Report these offsets, not new physical data. The total energy offset includes
the invariant frozen-ion chemical energy and elemental/storage gauge constants
as well as interpolation matching; its magnitude is NOT all interpolation error.

At each prescribed density change, solve chemistry and total internal energy
together using backward Euler and an adaptive step-doubling error estimate.
The energy equation is `e_new-e_old = -p_new*(v_new-v_old) + q_external*dt`.
The two-half-step state is accepted without extrapolation; positivity is
enforced by Newton line search, never by clipping completed populations. A failed
parcel remains unresolved in the denominator. Report work/forcing separately
and verify the integrated energy ledger.

Density varies log-linearly between saved samples. An optional parent forcing
is the sampled residual `delta(e_parent) + mean(p_parent)*delta(v_parent)`
per interval. It contains shock heating **and** sampling/EOS/discretization
defects; it is not an independently measured heat source. Compare replay of
that forcing with reversible-work-only paths. Do not prescribe both parent
temperature and energy, count the inferred forcing as recombination, or claim
that prescribed paths include chemistry's feedback on momentum/expansion.

## Chemistry and range interpretation

Forward rates come from ADR-0046; reverse rates use the study partition functions.
Export the six-species log-equilibrium reference at one temperature, then use
analytic translational/rotational/vibrational partition-function ratios in Rust.
This preserves van't Hoff with the same species internal energies, without an
independent interpolated equilibrium table. Generic and explicit H2/H2O collider
channels are not double-counted. Alternative Li H2 dissociation replaces all
three original H2 channels; reverse association is study-matched.

Documented source alternatives, channel removals, and group slowdowns are
**conditional scenarios**, not a probability distribution or evaluated hot/dense
uncertainty interval. The 1e3/1e6 association slowdowns retain Step 9's meaning.
Frozen chemistry is an exact zero-return control, not a claimed physical bound
on every coupled-flow outcome. No artificial high-pressure rate cap is invented.

Net return is `sum m*(bond_start-bond_end)`, with negative returns retained;
do not add successive positive releases. Fractions use the named start store
of the same parcel ensemble. Fixed ions can strand an additional molecular
store; the result is not a complete ion/molecular reaction model. Radiation,
wall heat transfer, HO2/H2O2, real injection and finite-plate escape are absent.
Results remain conditional on the lossless provisional EOS as well as rates.

## Acceptance gates

Before a flow-history result is quoted, test analytic frozen-gas expansion,
analytic reversible chemistry, first-order temporal convergence, positivity,
elemental and energy conservation, fixed-volume thermal relaxation, detailed
balance, and Rust/Python thermochemistry agreement. Refine solver tolerances and
history sampling; compare N40/N80. Preserve parent conservation/indexing gates.
Report unresolved mass/start-store coverage explicitly; do not renormalize it
away. This gate does not validate the physical extrapolation identified in
Step 10, nor produce an impulse gain.
