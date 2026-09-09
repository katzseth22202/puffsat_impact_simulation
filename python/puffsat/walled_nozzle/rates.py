"""Literature three-body recombination coefficients for the walled nozzle (N10.1-3).

**This module computes no rate coefficients of its own.** Every number below is a published
value, named with its source, its bath gas, the temperature window it was fitted over, and the
uncertainty the evaluation itself states. That discipline is inherited from `recombination.py`,
and it is the reason the 1080 s / 793 s fork can be reported as a *margin* rather than as a
point value: an uncertainty stated as a factor converts directly into decades of headroom.

**Where they came from.** All entries were read from the NIST Chemical Kinetics Database
(`kinetics.nist.gov`), which reproduces each evaluation's own fit rather than refitting it.
The squib in each docstring is NIST's record identifier, so any number here can be traced back
to the exact record. Everything is stated in NIST's own form and unit,

    k(T) = A (T / 298 K)^n     [cm^6 molecule^-2 s^-1]

and converted to SI once, by `_si`, because the rest of the codebase works in `m^6/s`.

**Two channels, and they are not equally well known.** Hydrogen's has six independent
evaluations and three shock-tube studies spanning 2500-7000 K. Nitrogen's rests on one
high-temperature shock-tube study and one modern ab initio calculation that disagree by more
than a decade. Carbon's barely exists -- see `CARBON_CHANNEL_NOTE`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

#: `1 cm^6 = 1e-12 m^6`. Named rather than inlined because it is applied to every constant here
#: and a silent factor of `1e12` would be invisible in a verdict.
CM6_IN_M6 = 1.0e-12

#: NIST states every fit against this reference temperature.
T_REFERENCE = 298.0


def _si(a_cm6: float) -> float:
    """Convert a published pre-exponential from `cm^6 molecule^-2 s^-1` to `m^6/s`."""
    return a_cm6 * CM6_IN_M6


@dataclass(frozen=True)
class ThreeBodyRate:
    """One published low-pressure three-body recombination coefficient.

    The fields are deliberately the ones a reader needs in order to *distrust* the number
    appropriately: which third body it was measured with, over what temperature window, how
    the evaluation itself was categorised, and what it costs to be wrong.
    """

    #: `A` in SI [m^6/s], already converted from the published `cm^6 molecule^-2 s^-1`.
    a: float
    #: Temperature exponent `n` in `k = A (T/298)^n`.
    n: float
    #: Third body the fit belongs to. **Not decorative** -- see `THIRD_BODY_EFFICIENCY`.
    bath_gas: str
    #: Fitted window [K]. Evaluating outside it is an extrapolation and `extrapolated` says so.
    t_min: float
    t_max: float
    #: Multiplicative uncertainty the source states (`3.16` means "a factor of 3.16 either way",
    #: i.e. exactly half a decade). `None` where the source states none.
    uncertainty_factor: float | None
    #: NIST record identifier, e.g. `1992BAU/COB411-429:216`.
    squib: str
    #: Short human citation.
    source: str
    #: NIST's own category: `"Review"` for an evaluation, `"Experiment"` for a measurement.
    category: str

    def __call__(self, temp: float) -> float:
        """`k(T)` [m^6/s]. No clamping: extrapolation is allowed but must be declared."""
        return float(self.a * (temp / T_REFERENCE) ** self.n)

    def extrapolated(self, temp: float) -> bool:
        """Whether `temp` lies outside the window this fit was actually established over."""
        return not (self.t_min <= temp <= self.t_max)

    def decades(self) -> float:
        """Stated uncertainty expressed in decades of rate, for comparison with a margin.

        A factor of 3.16 is 0.5 decades. Returned as `0.0` where the source states nothing,
        so a caller that sums exposures does not silently gain headroom from silence.
        """
        if self.uncertainty_factor is None or self.uncertainty_factor <= 1.0:
            return 0.0
        return math.log10(self.uncertainty_factor)


# ---- H + H + M -> H2 + M -----------------------------------------------------------------------
# Methane's recoverable store is the H2 one: 52% of full atomisation (`eos_methane`), against 43%
# locked in carbon. So for the flown fluid this reaction *is* N10.

#: Baulch et al., "Evaluated kinetic data for combustion modelling", J. Phys. Chem. Ref. Data
#: **21**, 411 (1992). The standard combustion evaluation, third body H2.
H_H2 = ThreeBodyRate(
    a=_si(8.85e-33),
    n=-0.60,
    bath_gas="H2",
    t_min=100.0,
    t_max=5000.0,
    uncertainty_factor=3.16,
    squib="1992BAU/COB411-429:216",
    source="Baulch et al. 1992, J. Phys. Chem. Ref. Data 21, 411",
    category="Review",
)

#: Cohen & Westberg, "Chemical kinetic data sheets for high-temperature chemical reactions",
#: J. Phys. Chem. Ref. Data **12**, 531 (1983). An *independent* evaluation of the same third
#: body that lands within 2% of `H_H2` with a tighter stated uncertainty. Two evaluations
#: agreeing that closely is the strongest evidence available that this channel is really known.
H_H2_COHEN = ThreeBodyRate(
    a=_si(9.04e-33),
    n=-0.60,
    bath_gas="H2",
    t_min=50.0,
    t_max=5000.0,
    uncertainty_factor=2.51,
    squib="1983COH/WES531:42",
    source="Cohen & Westberg 1983, J. Phys. Chem. Ref. Data 12, 531",
    category="Review",
)

#: Cohen & Westberg 1983, third body Ar. Carried as the inert reference the efficiency ratios
#: below are stated against.
H_AR = ThreeBodyRate(
    a=_si(6.48e-33),
    n=-1.00,
    bath_gas="Ar",
    t_min=77.0,
    t_max=5000.0,
    uncertainty_factor=2.0,
    squib="1983COH/WES531:44",
    source="Cohen & Westberg 1983, J. Phys. Chem. Ref. Data 12, 531",
    category="Review",
)

#: Hurle, Jones & Rosenfeld, "Shock-wave observations of rate constants for atomic hydrogen
#: recombination from 2500 to 7000 K", Proc. 8th Int. Shock Tube Symp. 253 (1969).
#:
#: **The only direct measurement inside our temperature window**, and the reason this study is
#: not extrapolating. Fitted as temperature-*independent* across 2500-7000 K, which is already
#: a physical statement: the steep negative exponents in the evaluations above are low-
#: temperature behaviour and do not continue upward. Extrapolating Baulch's `T^-0.6` to 5000 K
#: gives 1.6e-45 m^6/s where Hurle measures 4.8e-45 -- so the low-temperature fits are
#: *conservative* here by about a factor of 3, not optimistic.
H_H2_HOT = ThreeBodyRate(
    a=_si(4.83e-33),
    n=0.0,
    bath_gas="H2",
    t_min=2500.0,
    t_max=7000.0,
    uncertainty_factor=None,
    squib="1969HUR/JON253-276:1",
    source="Hurle, Jones & Rosenfeld 1969, Proc. 8th Int. Shock Tube Symp. 253",
    category="Experiment",
)

#: Hurle et al. 1969, third body Ar, same window.
H_AR_HOT = ThreeBodyRate(
    a=_si(1.68e-32),
    n=0.0,
    bath_gas="Ar",
    t_min=2500.0,
    t_max=7000.0,
    uncertainty_factor=None,
    squib="1969HUR/JON253-276:2",
    source="Hurle, Jones & Rosenfeld 1969, Proc. 8th Int. Shock Tube Symp. 253",
    category="Experiment",
)

#: Jacobs, Giedt & Cohen, J. Chem. Phys. **47**, 54 (1967), third body **atomic H**, 2900-4700 K.
#: See `THIRD_BODY_EFFICIENCY` for why this is the entry that matters for the flown charge.
H_ATOMIC_H = ThreeBodyRate(
    a=_si(1.85e-31),
    n=-1.00,
    bath_gas="H",
    t_min=2900.0,
    t_max=4700.0,
    uncertainty_factor=None,
    squib="1967JAC/GIE54-57:2",
    source="Jacobs, Giedt & Cohen 1967, J. Chem. Phys. 47, 54",
    category="Experiment",
)

#: Jacobs, Giedt & Cohen 1967, third body H2, from the *same* experiment as `H_ATOMIC_H`. Paired
#: with it deliberately: an efficiency ratio taken within one study cancels that study's
#: systematics, which a ratio across two studies does not.
H_H2_JACOBS = ThreeBodyRate(
    a=_si(2.31e-32),
    n=-1.00,
    bath_gas="H2",
    t_min=2900.0,
    t_max=4700.0,
    uncertainty_factor=None,
    squib="1967JAC/GIE54-57:3",
    source="Jacobs, Giedt & Cohen 1967, J. Chem. Phys. 47, 54",
    category="Experiment",
)

#: Jacobs, Giedt & Cohen 1967, third body Ar, same experiment.
H_AR_JACOBS = ThreeBodyRate(
    a=_si(9.26e-33),
    n=-1.00,
    bath_gas="Ar",
    t_min=2900.0,
    t_max=4700.0,
    uncertainty_factor=None,
    squib="1967JAC/GIE54-57:1",
    source="Jacobs, Giedt & Cohen 1967, J. Chem. Phys. 47, 54",
    category="Experiment",
)

#: Patch, J. Chem. Phys. **36**, 1919 (1962): the same three third bodies, 2950-5330 K, an
#: independent shock tube. Carried only to bound the *spread* in the efficiency of atomic H.
H_ATOMIC_H_PATCH = ThreeBodyRate(
    a=_si(4.63e-31),
    n=-1.00,
    bath_gas="H",
    t_min=2950.0,
    t_max=5330.0,
    uncertainty_factor=None,
    squib="1962PAT1919-1924:2",
    source="Patch 1962, J. Chem. Phys. 36, 1919",
    category="Experiment",
)

#: Patch 1962, third body Ar.
H_AR_PATCH = ThreeBodyRate(
    a=_si(6.95e-33),
    n=-1.00,
    bath_gas="Ar",
    t_min=2950.0,
    t_max=5330.0,
    uncertainty_factor=None,
    squib="1962PAT1919-1924:1",
    source="Patch 1962, J. Chem. Phys. 36, 1919",
    category="Experiment",
)

#: Rink, J. Chem. Phys. **36**, 262 (1962): third body atomic H, 2800-5000 K, a third
#: independent shock tube. The *low* end of the atomic-H efficiency range.
H_ATOMIC_H_RINK = ThreeBodyRate(
    a=_si(9.26e-32),
    n=-1.00,
    bath_gas="H",
    t_min=2800.0,
    t_max=5000.0,
    uncertainty_factor=None,
    squib="1962RIN262-265:2",
    source="Rink 1962, J. Chem. Phys. 36, 262",
    category="Experiment",
)

#: Rink 1962, third body H2, same experiment.
H_H2_RINK = ThreeBodyRate(
    a=_si(2.78e-32),
    n=-1.00,
    bath_gas="H2",
    t_min=2800.0,
    t_max=5000.0,
    uncertainty_factor=None,
    squib="1962RIN262-265:1",
    source="Rink 1962, J. Chem. Phys. 36, 262",
    category="Experiment",
)

#: **The load-bearing find of the N10 literature search**, as `(study, k_H/k_Ar, k_H2/k_Ar)`.
#:
#: Three independent shock-tube studies each measured `H + H + M` with M = Ar, H2 and *atomic
#: H* over 2800-5330 K, all fitted with the same `n = -1`. Taking the ratio within a single
#: study cancels its systematics, and all three agree on the ordering and roughly on the size:
#: **atomic hydrogen is an order of magnitude better at stabilising the collision than argon,
#: and several times better than H2.**
#:
#: This is what makes the walled-nozzle case stronger than a hand estimate suggests. At the
#: flown 10 kK the methane charge is 93-98% dissociated (W1), so the third body in the throat
#: is overwhelmingly *atomic hydrogen*, not H2 and certainly not an inert. A margin computed on
#: an H2 or generic third body understates the recombination rate by roughly a factor of 5-10.
THIRD_BODY_EFFICIENCY = (
    ("Jacobs, Giedt & Cohen 1967", 1.85e-31 / 9.26e-33, 2.31e-32 / 9.26e-33),
    ("Rink 1962", 9.26e-32 / 1.39e-32, 2.78e-32 / 1.39e-32),
    ("Patch 1962", 4.63e-31 / 6.95e-33, 6.95e-32 / 6.95e-33),
)


# ---- N + N + M -> N2 + M -----------------------------------------------------------------------
# Ammonia's rung only. Methane contains no nitrogen; carrying this reaction against the flown
# fluid would be a category error, which is why `freeze.py` selects the channel by propellant.

#: Byron, "Shock-tube measurement of the rate of dissociation of nitrogen", J. Chem. Phys. **44**,
#: 1378 (1966), third body N2, obtained by detailed balance from the measured dissociation rate.
#: Its 6000-9000 K window *straddles the chamber temperature*, so unlike the hydrogen
#: evaluations it needs no extrapolation upward at all.
N_N2 = ThreeBodyRate(
    a=_si(4.16e-33),
    n=-0.50,
    bath_gas="N2",
    t_min=6000.0,
    t_max=9000.0,
    uncertainty_factor=None,
    squib="1966BYR1378-1388:4",
    source="Byron 1966, J. Chem. Phys. 44, 1378",
    category="Experiment",
)

#: Byron 1966, third body **atomic N** -- the relevant one in a dissociated ammonia plume, and
#: again far more efficient than the molecule, mirroring hydrogen's behaviour.
N_ATOMIC_N = ThreeBodyRate(
    a=_si(1.29e-30),
    n=-1.50,
    bath_gas="N",
    t_min=6000.0,
    t_max=9000.0,
    uncertainty_factor=None,
    squib="1966BYR1378-1388:5",
    source="Byron 1966, J. Chem. Phys. 44, 1378",
    category="Experiment",
)

#: Byron 1966, third body Ar.
N_AR = ThreeBodyRate(
    a=_si(1.60e-33),
    n=-0.50,
    bath_gas="Ar",
    t_min=6000.0,
    t_max=9000.0,
    uncertainty_factor=None,
    squib="1966BYR1378-1388:6",
    source="Byron 1966, J. Chem. Phys. 44, 1378",
    category="Experiment",
)

#: Notey, Jo & Panesi, "Master equation study of three-body recombination of nitrogen and oxygen
#: in non-equilibrium hypersonic flows", J. Chem. Phys. **163**, 194311 (2025); arXiv:2506.17452.
#:
#: A state-to-state master equation on ab initio potential energy surfaces, run in exactly our
#: physical setup: atoms held at 10 000 K and then plunged into a cold bath, which is what a
#: nozzle does to a chamber. Reported as `(temperature [K], k [m^6/s])` for `N + N + N` in the
#: quasi-steady state; the paper's own `kbar` column, which it shows agrees with the value
#: obtained by micro-reversibility from the QSS dissociation rate.
#:
#: **It disagrees with Byron by more than a decade** once Byron's fit is brought down to the
#: same temperatures. That disagreement is not a defect in this study -- it is the honest width
#: of the nitrogen channel, and it is why `freeze.py` reports nitrogen as a band.
N_ABINITIO = ((2500.0, 9.125e-46), (5000.0, 7.60e-46))


# ---- Carbon ------------------------------------------------------------------------------------

#: Slack, "Kinetics and Thermodynamics of the CN Molecule III", J. Chem. Phys. **64**, 228 (1976):
#: `C + C + M -> C2 + M` with M = Ar, 5000-6000 K, `5.46e-31 (T/298)^-1.60` with a stated
#: +/-3.03e-31 -- an uncertainty of 55% on the only measurement there is.
C_AR = ThreeBodyRate(
    a=_si(5.46e-31),
    n=-1.60,
    bath_gas="Ar",
    t_min=5000.0,
    t_max=6000.0,
    uncertainty_factor=None,
    squib="1976SLA228:4",
    source="Slack 1976, J. Chem. Phys. 64, 228",
    category="Experiment",
)

#: **Why methane's carbon store gets no three-body answer, stated as a finding rather than a gap.**
#:
#: The NIST evaluated database holds exactly one measurement of `C + C + M -> C2 + M` (`C_AR`,
#: one third body, a 1000 K window, 55% uncertainty) and **no record at all** for
#: `C + H + M -> CH + M`. Hydrogen's channel, by contrast, has six independent evaluations and
#: three shock-tube studies spanning 2500-7000 K.
#:
#: That asymmetry is itself informative. An association channel gets measured when it controls
#: something measurable, and the carbon channels do not, because carbon in a cooling
#: carbon-rich gas does not go to diatomics -- it goes to clusters and then to soot. So the 43%
#: of methane's atomisation locked in carbon cannot be settled by a rate coefficient the way
#: the 52% held as H2 can. This is independent support for ADR-0050's exclusion of condensed
#: carbon from the EOS, and for N10 item 5 being a *nucleation* problem needing different
#: machinery rather than one more entry in this module.
CARBON_CHANNEL_NOTE = (
    "C + C + M has one measurement (Slack 1976, M = Ar, 5000-6000 K, +/-55%); "
    "C + H + M has none. The carbon store is a nucleation problem, not a three-body one."
)
