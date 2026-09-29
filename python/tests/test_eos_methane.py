"""Acceptance tests for the methane equilibrium EOS (ADR-0050, asks N9/N10).

The pattern is `test_eos_water.py`'s: the unit under test is *the EOS reproducing a state whose
answer is known before the solver exists*, not a micro-function. The known answers here are the
0 K thermochemistry (checkable by hand against NIST-JANAF), the two endpoints (bound methane,
stripped plasma), element conservation, and a hand-solved `H2 <=> 2H` equilibrium.
"""

from __future__ import annotations

import math

import pytest

from puffsat import eos_methane as eos

CHAMBER_DENSITIES = (2.495, 1.2475, 0.7415)  # 499 kg over 200, 400 and 673 m^3 (ADR-0016)


# --- Thermochemistry: hand-checkable against NIST-JANAF 0 K heats of formation -----------------


def test_methane_atomization_matches_janaf() -> None:
    """`C + 4H - CH4` at 0 K is 1641.96 kJ/mol, i.e. 102.35 MJ/kg -- and the mean C-H bond it
    implies (410 kJ/mol) is the textbook value, which is the check that no sign is wrong."""
    assert eos.CH4.d_at * eos.N_A / 1e3 == pytest.approx(1641.96, abs=0.05)
    assert pytest.approx(102.35e6, rel=1e-3) == eos.FULL_ATOMIZATION_ENERGY
    assert eos.CH4.d_at * eos.N_A / 1e3 / 4.0 == pytest.approx(410.5, abs=1.0)


def test_recombination_channels_split_the_store_as_the_adr_says() -> None:
    """ADR-0016 states methane returns 52% of its atomisation as H2 and 43% as condensing carbon.

    That split is derived here from the heats of formation alone, so agreeing with it is an
    independent confirmation of both sides' thermochemistry -- and the 4% that is neither is
    methane's own formation enthalpy, which never returns.
    """
    assert pytest.approx(0.526, abs=0.005) == eos.H2_RETURN_ENERGY / eos.FULL_ATOMIZATION_ENERGY
    assert (
        pytest.approx(0.433, abs=0.005)
        == eos.CARBON_CONDENSATION_ENERGY / eos.FULL_ATOMIZATION_ENERGY
    )
    total = eos.H2_RETURN_ENERGY + eos.CARBON_CONDENSATION_ENERGY + eos.FORMATION_ENTHALPY
    assert total == pytest.approx(eos.FULL_ATOMIZATION_ENERGY, rel=1e-12)


def test_h2_constants_agree_with_the_water_eos() -> None:
    """H2 appears in both species sets. If the two tables ever disagreed, a freeze verdict would
    depend on which module asked -- so they are pinned equal rather than merely both present."""
    from puffsat import eos_water

    assert eos.H2.d_at == pytest.approx(eos_water.H2.d0, rel=1e-12)
    assert eos.H2.theta_rot[0] == pytest.approx(eos_water.H2.theta_rot, rel=1e-12)
    assert eos.H2.theta_vib[0] == pytest.approx(eos_water.H2.theta_vib, rel=1e-12)


# --- Conservation: the solve's own constraints, checked on the reconstructed composition -------


@pytest.mark.parametrize("rho", CHAMBER_DENSITIES)
@pytest.mark.parametrize("temp", [500.0, 1500.0, 3000.0, 5000.0, 8000.0, 10000.0, 15000.0, 30000.0])
def test_element_conservation_and_neutrality(rho: float, temp: float) -> None:
    """C nuclei conserved, H:C = 4:1, and charge neutral -- checked on the *reconstructed*
    densities rather than on the residual, so a bug in `_densities` cannot hide behind a
    converged solve."""
    comp = eos.composition(rho, temp)
    n_f = rho / eos.M_CH4
    assert comp.n_c_nuclei == pytest.approx(n_f, rel=1e-6)
    assert comp.n_h_nuclei == pytest.approx(4.0 * n_f, rel=1e-6)
    charge = comp.n_hp + sum((k + 1) * n for k, n in enumerate(comp.n_c_ions))
    assert comp.n_e == pytest.approx(charge, rel=1e-6)


# --- Endpoints: both are known without a solver ------------------------------------------------


def test_cold_endpoint_is_bound_methane_at_zero_energy() -> None:
    """At 500 K the charge is intact CH4: one particle per formula unit, no free store, and an
    internal energy that is sensible heat only. `e` must still be strictly positive -- the Rust
    table loader interpolates `ln e`."""
    rho = 1.2475
    comp = eos.composition(rho, 500.0)
    n_f = rho / eos.M_CH4
    assert comp.density_of(eos.CH4) / n_f > 0.999
    assert comp.hydrogen_dissociated < 1e-9
    p, e = eos.pressure_energy(rho, 500.0)
    assert 0.0 < e < 2.0e6
    assert p == pytest.approx(n_f * eos.K_B * 500.0, rel=1e-3)
    # The trace C2H2/H2 pool holds a real but negligible store: ~0.7 J/kg of 102 MJ/kg.
    assert eos.bond_energy_held(rho, 500.0) / eos.FULL_ATOMIZATION_ENERGY < 1e-7


def test_hot_endpoint_strips_to_atoms_and_holds_the_whole_store() -> None:
    """At 30 kK nothing molecular survives, the bond store is fully charged, and the ionisation
    energy is on top of it."""
    rho = 1.2475
    comp = eos.composition(rho, 30000.0)
    assert comp.hydrogen_dissociated > 0.999
    assert comp.carbon_free > 0.999
    held = eos.bond_energy_held(rho, 30000.0)
    assert held == pytest.approx(eos.FULL_ATOMIZATION_ENERGY, rel=1e-3)
    _, e = eos.pressure_energy(rho, 30000.0)
    assert e > eos.FULL_ATOMIZATION_ENERGY


# --- The number ask N10.4b turns on -------------------------------------------------------------


def test_hydrogen_dissociation_matches_a_hand_solved_equilibrium() -> None:
    """`H2 <=> 2H` at 10 kK, solved by hand from the module's own `ln K`, must match the full
    nine-species solve to a few percent.

    This is the acceptance test for item 4b. At 10 kK the carbon-bearing molecules are gone, so
    the full solve *should* reduce to the single binary equilibrium; if it does not, the extra
    species are doing something they should not.
    """
    rho, temp = 2.495, 10000.0
    k_h2 = math.exp(eos.ln_k_molecule(eos.H2, temp))  # n_H^2 / n_H2 [m^-3]
    n_nuclei = 4.0 * rho / eos.M_CH4
    # n_H^2 / ((n_nuclei - n_H)/2) = K  ->  2 n_H^2 + K n_H - K n_nuclei = 0
    n_h = (-k_h2 + math.sqrt(k_h2**2 + 8.0 * k_h2 * n_nuclei)) / 4.0
    hand = n_h / n_nuclei

    comp = eos.composition(rho, temp)
    assert comp.hydrogen_dissociated == pytest.approx(hand, rel=0.05)
    assert hand > 0.9  # the standing paper-side estimate of 49-76% is not what equilibrium gives


@pytest.mark.parametrize("rho", CHAMBER_DENSITIES)
def test_dissociation_rises_as_the_chamber_grows(rho: float) -> None:
    """Le Chatelier, and it is the direction ADR-0016 argues from: a bigger (thinner) chamber
    charges more of the store. The finding is that the *size* of the effect is small at 10 kK,
    not that the sign is wrong."""
    lo = eos.composition(CHAMBER_DENSITIES[0], 10000.0).hydrogen_dissociated
    here = eos.composition(rho, 10000.0).hydrogen_dissociated
    assert here >= lo - 1e-12
    assert 0.9 < here < 1.0


# --- Thermodynamic consistency ------------------------------------------------------------------


@pytest.mark.parametrize("temp", [3000.0, 8000.0, 12000.0])
def test_energy_and_pressure_are_monotone_in_temperature(temp: float) -> None:
    """`de/dT > 0` and `dp/dT > 0` at fixed density: c_v > 0 is the minimum a usable EOS owes,
    and a mis-signed chemical term is the way this one could fail."""
    rho = 1.2475
    p_lo, e_lo = eos.pressure_energy(rho, temp)
    p_hi, e_hi = eos.pressure_energy(rho, temp * 1.02)
    assert e_hi > e_lo
    assert p_hi > p_lo


def test_sound_speed_is_finite_and_sane_in_the_chamber() -> None:
    """The sealed-vessel argument of N9 item 0 rests on the chamber sound speed, so it gets a
    test: an equilibrium methane plasma at 10 kK runs a few km/s, well below the frozen-ideal
    value because dissociation absorbs the compression."""
    c_s = eos.sound_speed(1.2475, 10000.0)
    assert 2.0e3 < c_s < 2.0e4


def test_half_dissociation_needs_thirty_times_the_chamber_density() -> None:
    """The quantitative form of why the paper-side 49-76% estimate cannot be right.

    Half dissociation of `H2 <=> 2H` happens where the H-nucleus density equals the equilibrium
    constant itself (`n_H = N/2`, `n_H2 = N/4`, so `n_H^2/n_H2 = N`). At 10 kK that constant is
    1.2e28 m^-3, which is 81.5 kg/m^3 of methane -- thirty times the densest candidate chamber.
    Pinned because the answer document quotes it as the reason the correction overshot.
    """
    k_h2 = math.exp(eos.ln_k_molecule(eos.H2, 10000.0))
    assert k_h2 == pytest.approx(1.22e28, rel=0.02)
    rho_half = k_h2 / 4.0 * eos.M_CH4
    assert rho_half == pytest.approx(81.5, rel=0.02)
    assert rho_half / max(CHAMBER_DENSITIES) > 30.0
    assert eos.composition(rho_half, 10000.0).hydrogen_dissociated == pytest.approx(0.46, abs=0.03)


# --- Feeds other than pure methane (the near-term chamber) --------------------------------------

A1_FEED = eos.mixture_feed("A1", [(68.8, eos.METHANE_GAS), (6.9, eos.POLYETHYLENE)])


def test_default_feed_is_pure_methane() -> None:
    """`METHANE` is the default everywhere and must reproduce the pre-feed EOS exactly -- an
    explicitly built 4:1 methane feed takes the same code path and returns the same numbers."""
    built = eos.mixture_feed("CH4", [(1.0, eos.METHANE_GAS)])
    assert built.hc_ratio == pytest.approx(4.0, rel=1e-14)
    assert built.reference_atomization == pytest.approx(eos.CH4.d_at, rel=1e-14)
    assert eos.pressure_energy(3.4, 7000.0) == eos.pressure_energy(3.4, 7000.0, eos.METHANE)


def test_polyethylene_atomisation_is_its_combustion_heat_restated() -> None:
    """(CH2)n: 1143.3 kJ/mol of atoms less its -27 kJ/mol heat of formation, 83.4 MJ/kg."""
    per_kg = eos.POLYETHYLENE.atomization / eos.POLYETHYLENE.mass
    assert per_kg == pytest.approx(83.43e6, rel=2e-3)


@pytest.mark.parametrize("temp", [1500.0, 3000.0, 5000.0, 7000.0, 10000.0])
def test_mixed_feed_conserves_its_own_stoichiometry(temp: float) -> None:
    """At H:C = 3.77 the general H:C row replaces CH4's divided-out form; nuclei and charge
    must still close on the reconstructed densities."""
    rho = 3.4
    comp = eos.composition(rho, temp, A1_FEED)
    n_f = rho / A1_FEED.formula_mass
    assert comp.n_c_nuclei == pytest.approx(n_f, rel=1e-6)
    assert comp.n_h_nuclei == pytest.approx(A1_FEED.hc_ratio * n_f, rel=1e-6)
    charge = comp.n_hp + sum((k + 1) * n for k, n in enumerate(comp.n_c_ions))
    assert comp.n_e == pytest.approx(charge, rel=1e-6)


@pytest.mark.parametrize(("rho", "temp"), [(1.0, 5500.0), (3.0, 5500.0), (0.2, 10000.0)])
def test_trace_carbon_feed_is_the_pure_hydrogen_eos(rho: float, temp: float) -> None:
    """A hydrogen feed with 1e-5 carbon must be `walled_nozzle.hydrogen`, an independent solver
    on the same H2 constants, to the size of the carbon it carries."""
    from puffsat.walled_nozzle import hydrogen

    feed = eos.mixture_feed("H2", [(1.0, eos.HYDROGEN_GAS), (1e-5, eos.METHANE_GAS)])
    p, e = eos.pressure_energy(rho, temp, feed)
    p_ref, e_ref = hydrogen.pressure_energy(rho, temp)
    assert p == pytest.approx(p_ref, rel=5e-5)
    assert e == pytest.approx(e_ref, rel=5e-5)


# --- The opt-in JANAF partition-function correction --------------------------------------------


@pytest.mark.parametrize("name", ["CH4", "CH", "C2", "C3", "C2H2", "H2"])
@pytest.mark.parametrize("temp", [2000.0, 4000.0, 6000.0])
def test_janaf_thermo_reproduces_janaf_up_to_heat_of_formation(name: str, temp: float) -> None:
    """With `thermo="janaf"` a molecule's atomisation constant is JANAF's, except for the part a
    different 0 K heat of formation explains -- which is deliberately kept (C2H's JANAF value is
    ~92 kJ/mol out of date)."""
    from puffsat import janaf_carbon

    mol = next(m for m in eos.MOLECULES if m.name == name)
    i = janaf_carbon.TEMPERATURES.index(temp)
    row = {k: v[i] for k, v in janaf_carbon.LOG_KF.items()}
    hf_j = 0.0 if name == "H2" else janaf_carbon.HF0[name]
    offset = (
        (
            mol.n_c * (eos._HF0_C - janaf_carbon.HF0["C"])
            + mol.n_h * (eos._HF0_H - janaf_carbon.HF0["H"])
            - (eos._HF0[name] - hf_j)
        )
        * 1e3
        / (eos.K_B * eos.N_A * temp)
    )
    expected = eos._janaf_ln_k(mol, temp, row) - offset
    assert eos.ln_k_molecule_thermo(mol, temp, "janaf") == pytest.approx(expected, abs=0.01)


def test_harmonic_c3_is_over_bound_thirtyfold_at_6000_k() -> None:
    """The reason the correction exists, pinned so a future EOS change that fixes C3 directly
    shows up here: harmonic C3 is ~30x too stable, harmonic C2 ~9x not stable enough."""
    assert eos.janaf_correction(eos.C3, 6000.0)[0] == pytest.approx(-3.45, abs=0.05)
    assert eos.janaf_correction(eos.C2, 6000.0)[0] == pytest.approx(2.21, abs=0.05)


@pytest.mark.parametrize("mol_name", ["C3", "C2", "H2"])
def test_janaf_energy_obeys_van_t_hoff(mol_name: str) -> None:
    """`d ln K_n / dT = (E_atoms - E_molecule) / kT^2` must hold with the corrected `K` *and* the
    corrected molecular energy -- the check that `e` and `K` are one thermodynamics."""
    mol = next(m for m in eos.MOLECULES if m.name == mol_name)
    temp, dt = 5000.0, 1.0
    slope = (
        eos.ln_k_molecule_thermo(mol, temp + dt, "janaf")
        - eos.ln_k_molecule_thermo(mol, temp - dt, "janaf")
    ) / (2 * dt)
    kt = eos.K_B * temp
    e_atoms = (mol.n_c + mol.n_h) * 1.5 * kt
    e_mol = (
        1.5 * kt
        + eos._e_rot(mol, temp)
        + eos._e_vib(mol, temp)
        + kt * temp * eos.janaf_correction(mol, temp)[1]
        - mol.d_at
    )
    assert slope * kt * temp == pytest.approx(e_atoms - e_mol, rel=1e-6)


@pytest.mark.parametrize("hc_ratio_feed", ["hydrogen_rich", "carbon_rich"])
@pytest.mark.parametrize("temp", [150.0, 366.0, 1200.0])
@pytest.mark.parametrize("rho", [1e-5, 0.1, 3.0])
def test_off_stoichiometric_feeds_converge_cold(
    hc_ratio_feed: str, temp: float, rho: float
) -> None:
    """A deep expansion takes the gas to a few hundred kelvin, where a 4:1 atom-pool init is
    decades off for a feed that is not 4:1. The molecular init must carry these."""
    feed = (
        eos.mixture_feed("H2+PE", [(60.0, eos.HYDROGEN_GAS), (6.9, eos.POLYETHYLENE)])
        if hc_ratio_feed == "hydrogen_rich"
        else A1_FEED
    )
    comp = eos.composition(rho, temp, feed)
    n_f = rho / feed.formula_mass
    assert comp.n_h_nuclei == pytest.approx(feed.hc_ratio * n_f, rel=1e-6)


@pytest.mark.parametrize("hc_ratio", [0.05, 0.5, 1.0, 2.0])
@pytest.mark.parametrize("temp", [2500.0, 3900.0])
def test_carbon_rich_feeds_converge(hc_ratio: float, temp: float) -> None:
    """A gas that has dissolved graphite from a wall can be far below H:C = 1; the molecular
    init must not take the log of a negative CH4 share there."""
    feed = eos.Feed("carbon-rich", hc_ratio, eos.CH4.d_at)
    comp = eos.composition(1.0, temp, feed)
    n_f = 1.0 / feed.formula_mass
    assert comp.n_c_nuclei == pytest.approx(n_f, rel=1e-6)
    assert comp.n_h_nuclei == pytest.approx(hc_ratio * n_f, rel=1e-6)
