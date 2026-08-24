"""Q-P: does more slug mass buy back the frozen dissociation store? Three levers, all weak.

Raised by Seth 2026-08-24 -- "larger mass would lead to more collisions and more energy recovery,
right?" The direction is right and the mechanism is not: at fixed density a bigger blob has the
*same* collision rate per particle. What it has is more time, because it is physically larger and
takes longer to thin out. That makes the scaling a cube root rather than linear:

    Da ~ tau_exp * n^2 ~ (R/u) * (M/R^3)^2 = M^(1/3) rho^(5/3) / u        =>  rho_freeze ~ M^(-1/5)

Three readings of "more mass", plus the lever aimed at the trigger itself. Monkeypatches the
module-level geometry constants, which is why this lives here rather than in `python/`.
"""

import sys
from itertools import pairwise

sys.path.insert(0, "python")

import puffsat.recombination as rc
from puffsat import eos_water, expansion, fireball, plume

COLD = 45.58e3
R0, L0, A0 = expansion.THROAT_RADIUS, expansion.FIELD_LENGTH, expansion.AREA_RATIO_EXIT
BOND = eos_water.D_AT / eos_water.M_H2O  # 50.9 MJ/kg, the whole store


def freeze(temp_0, rho_0):
    """First station where the expansion outruns atomic recombination, on current geometry."""
    fireball._cached_history.cache_clear()
    rows = fireball._history(temp_0, 45.0, rho_0, 320, fireball.EXPANSION_RATIO)
    for prev, row in pairwise(rows[::4]):
        c = eos_water.composition(row.rho, row.temp)
        n3 = c.n_h2o + c.n_h + c.n_o + c.n_hp + sum(c.n_o_ions) + c.n_e
        st = rc.freeze_station(
            time=row.time, temp=row.temp, rho=row.rho, n_e=c.n_e, n_atom=c.n_h, n_third=n3,
            tau_expansion=rc.expansion_time(prev.time, prev.rho, row.time, row.rho),
            dissociated_fraction=row.dissociated_fraction)
        if st.da_dissociation < 1.0:
            return st
    return None


def show(label, value, st):
    """`stranded` is what freezes in; `returned` is what recombination gave back before it did.
    The lever's worth is the *difference* between a row's `stranded` and the baseline row's."""
    held = st.dissociated_fraction
    print(f"  {label:>12} {value:>10} {st.rho:11.4e} {held:11.4f} "
          f"{held * BOND / 1e6:9.1f}M {(1 - held) * BOND / 1e6:9.1f}M")


HEAD = f"  {'':>12} {'':>10} {'rho_freeze':>11} {'f_diss held':>11} {'stranded':>10} {'returned':>10}"

print("(1) BIGGER BAG, SAME DENSITY -- self-similar scale-up, the fixed-k case")
print("    Da ~ M^(1/3): more *time*, not more collisions.")
print(HEAD)
for lam in (1.0, 2.0, 4.0):
    expansion.THROAT_RADIUS, expansion.FIELD_LENGTH = R0 * lam, L0 * lam
    show("mass x", f"{lam**3:.0f}", freeze(14700.0, expansion.BAG_RHO))
expansion.THROAT_RADIUS, expansion.FIELD_LENGTH = R0, L0

print()
print("(2) SAME BAG, MORE MASS IN IT -- higher density. The strongest of the three.")
print("    tau_rec ~ rho^-2 directly, AND a denser start at equal specific energy is a")
print("    lower-entropy start, so the plume runs colder downstream and returns more.")
print(HEAD)
for rho_0 in (expansion.BAG_RHO, 1.0, 3.0):
    show("mass x", f"{rho_0 / expansion.BAG_RHO:.1f}", freeze(plume.plume_state(COLD, rho_0).temp, rho_0))

print()
print("(3) LONGER NOZZLE -- aimed at the trigger (the field letting go), same opening rate.")
print(HEAD)
rate = L0 / (A0 - 1.0)
for a_exit in (4.0, 8.0, 16.0, 32.0):
    expansion.AREA_RATIO_EXIT, expansion.FIELD_LENGTH = a_exit, rate * (a_exit - 1.0)
    show("bore m", f"{expansion.FIELD_LENGTH:.0f}", freeze(14700.0, expansion.BAG_RHO))
expansion.AREA_RATIO_EXIT, expansion.FIELD_LENGTH = A0, L0

print()
print("None of the three is a rescue. Recombination needs rho^2 and the clock only gives R,")
print("so the geometry is against it. What would actually move this is the missing OH")
print("chemistry, not scale -- the scaling laws above are geometric and would survive it.")
