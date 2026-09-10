# The N16 grid runs on a cone, and the N9 load case is a snowplow with a clock

Status: **accepted**, 2026-09-10. Answers asks **N16** and **N9 items 1-7** of
`Balloon-Pulse-Propulsion` @ `dd6ad2a`, `docs/nozzle_asks_for_impact_sim.md`.
Extends [ADR-0050](0050-walled-thermal-nozzle-methane-chamber.md); every scope
statement there still holds.

## Context

N16 replaces N14 and N15 and asks for one surface over chamber temperature,
chamber volume and throat area, for methane **and** water, plus the density
dependence of the held store, which it calls "the most valuable thing this ask
can return". N9 items 1-7 ask what the arriving front does to the wall, and the
ask itself names the gate: below about 170 m^3 the front may leave the column
before its cone reaches the liner, which would close items 1-4 rather than
answer them.

Two things about the geometry moved on the paper side between ADR-0050 and this
ask, and they are what force new decisions rather than a wider sweep of the old
machinery: chamber volumes fell to 50-400 m^3, and throat areas fell to
0.05-0.2 m^2 -- area ratios of 140 to 600 against the flown 4.

## Decision 1: the throat is not a third dimension of the grid

A chamber at `(rho_c, T_c)` has exactly one isentrope. The throat area decides
only **where along it** the gas stops, because the exit area is the 3 m bore in
every case. So `surface.py` runs one deep expansion per `(fluid, T_c, V)` out to
the widest area ratio in the grid and reads every throat off it.

That is 32 expansions rather than 256, and -- the reason it is a decision and not
an optimisation -- it makes the throat comparison exactly like-for-like: the same
gas, the same adiabat, one cut point moved. `test_walled_nozzle_surface.py` pins a
point read off the deep isentrope against an expansion run directly to that ratio.

## Decision 2: the nozzle clock is rebuilt on a cone, not inherited from a length

`expansion.cooling_history` opens the area **linearly in `A`** over a length the
caller supplies. That was right for `freeze.py`, where the geometry was a field
profile of fixed length. It is wrong across two decades of area ratio: at fixed
length a larger area ratio means a faster opening, so the Damkoehler number would
fall with area ratio for a reason that is pure bookkeeping.

So the grid puts a **conical nozzle at 15 degrees** from the throat radius out to
the 3 m bore and integrates `dt = dx/u` along it:

    x(A) = (sqrt(A/pi) - r*) / tan(theta),   L = (R_bore - r*) / tan(theta)

Consequences, both reported rather than hidden: across the throats the ask
recommends (0.5 m^2 and below) the length is nearly constant at 9.7-10.7 m, so
those columns differ only in where the cone starts; the flown 7 m^2 throat is a
genuinely *shorter* nozzle at 5.6 m, because a 1.5 m throat radius is already half
the bore. `Da` scales linearly in the length, so the half-angle is a stated
sensitivity and not a hidden one.

## Decision 3: the equilibrium branch is reported with a freeze cap beside it

The deep throats this ask wants **do** freeze -- which is new, and is the opposite
of what W5 found inside the flown bore. The equilibrium conversion past a freeze
station is not a physical number, so every grid row carries a second column,
`conversion_capped`, which charges the equilibrium branch for the chemical energy
it releases *downstream* of the freeze station.

It is a **lower bound**, and deliberately: the frozen gas is also colder, so part
of that energy would not have reached the jet as thermal enthalpy either. The
truth is between the two columns, and the two coincide wherever nothing freezes.
Re-solving a genuinely frozen methane expansion needs a
`pressure_energy_frozen` that `eos_methane` does not have (`eos_water` does), and
that is recorded as deferred rather than approximated.

## Decision 4: the mixture N16 names is assembled, not solved, and says so

N16 asks to price "90% water with 10% liquid hydrogen" first. There is no joint
H/O equilibrium EOS here, so `surface.mixture` composes `eos_water` and
`hydrogen` as **two independent subsystems** at a common temperature and their own
partial densities. That is exact for the thermal terms and each component's
internal chemistry, and wrong about exactly one thing: the free hydrogen does not
shift water's dissociation equilibrium. Le Chatelier puts the real mixture's held
store *above* this, so the mixture rung is an optimistic edge and carries
`estimated=True` on every row -- the same discipline `propellants.ammonia_rung`
already uses.

## Decision 5: the front's wall load is priced on the layer's clock, not its temperature

`wall.py` integrates `coupling.snowplow`'s system with a clock and a general EOS,
pinned in the tests against `front.integrate` on the magnetic geometry so the two
models cannot drift apart.

The load case is then split three ways and only one of them is a real delivery:

- **Radiative** is `held x (1 - exp(-dt/tau_rad))`, with `tau_rad` the binding
  branch of `expansion.radiation_check`, capped at `sigma T^4 dt`. **That cap
  alone settles the item**, which is why the verdict does not depend on the
  opacity model.
- **Convective** is reported as an *incident* stagnation-enthalpy flux times a
  Stanton bracket of 0.001-0.01, because gas arriving at a wall does not deposit
  its energy -- the boundary layer decides what crosses. This is the crudest
  number in the module and it is the one that binds, which is exactly what the
  ask says of its own item 7.
- **The ceiling** -- everything the shocked column holds -- is printed beside both
  so a reader can see how far from it the answers are.

## Decision 6: the opacity is the free-free floor, bracketed, and the floor is the conservative choice

The shocked layer is a partly ionised C/H plasma at tens of thousands of kelvin.
`plasma_opacity` computes Thomson plus Kramers free-free on the **solved**
composition and nothing else. Bound-free and line opacity are omitted, and every
omitted term makes the layer *thicker*; a thicker layer in the diffusion regime
delivers *less*. So the floor is conservative for a question about how hard the
wall is hit, and `OPACITY_BRACKET` carries 0.1x-10x to show the verdict surviving
it. This is the same posture `opacity_bracket.py` takes for the magnetic nozzle
and the opposite of quoting one number from a table this study does not have.

## Decision 7: carbon at the throat is a thermodynamic verdict plus a kinetic ceiling

Item 5 asks for net deposition per pulse. What is owned here is where deposition
is *favoured* -- graphite's Clausius-Clapeyron vapour pressure, anchored at its
measured 1 atm sublimation point, against the solved monatomic-carbon partial
pressure -- and a Hertz-Knudsen **ceiling** on the rate with unit sticking and no
boundary layer.

The ceiling is deliberately useless as a prediction (it is centimetres per pulse)
and useful as a bound: it says the answer is set by the boundary layer and not by
the thermodynamics, so a real number needs a boundary-layer solve this study does
not have. What the thermodynamics does settle is the verdict, and the verdict is
that no wall temperature both survives and self-cleans.

## Consequences

- New module `python/puffsat/walled_nozzle/surface.py` (N16) and
  `python/puffsat/walled_nozzle/wall.py` (N9 items 1-7).
- New targets `make walled-nozzle-surface` and `make walled-nozzle-wall`, both
  added to `make walled-nozzle-test`.
- Six new committed CSVs under `data/results/walled_nozzle/`, un-ignored for the
  cross-repo reason ADR-0025 and the existing walled-nozzle entries already record.
- The answer document is `docs/walled_nozzle_grid_and_wall.md`, numbered W10
  onward so it continues the W series the paper repository already cites.
