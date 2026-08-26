# The seed window's `v L` is an output of the solved expansion, not a back-solve of the paper

`tab:seed_window`'s magnetic Reynolds number is `Rm = mu0 sigma(T) (v L)`. This repository owns
`sigma`; it did **not** own `v L`. The paper reports the `Rm` column without stating either the
expansion speed or the field length (audit Q-G), so the column cannot be reproduced from the paper
alone, and `conductivity.REF_V_L` was seeded by **back-solving** it: the product that makes the
paper's own `Rm = 361` row at 15 000 K come out against a hand-computed conductivity. That gave
`1.81e4 m^2/s`.

Two things were wrong with it, and they compound.

1. **It encodes a conductivity this module contradicts.** `361 / (mu0 x 1.81e4) = 15 872 S/m`,
   against the `6993 S/m` computed here — the same 2.3x gap that
   `test_blended_conductivity_reproduces_the_audit_hand_calculation` records against the audit's
   corrected hand blend of ~6950 S/m. So at `1.81e4` this module gives `Rm = 159`, not 361: the
   constant and the row it was fitted to cannot both be right, and the module says the row is the
   one that is wrong. Fitting a constant to that row propagated the error into `Rm` instead.
2. **`v L` had since become an output.** `expansion.HistoryRow.v_l` is the local flow speed times
   the local flux-tube radius at every station of the solved cooling history. Once that exists,
   taking `v L` from a back-solved table is choosing a worse source over a better one that is
   already in the repository.

## Decision

`REF_V_L = 7.4e4 m^2/s` — the **nozzle-exit** `v L` of the 56.53 km/s equilibrium leg (7.44e4), the
leg the paper states. `V_L_BAND = (5.5e4, 9.7e4)` is what the four flown legs actually span at the
exit: 5.54e4 at 45.58 km/s on the frozen branch, 9.72e4 at 75 km/s on the equilibrium branch. Both
are reproduced by the last station of `data/results/cooling_history.csv` (`make analysis-expansion`).

`RETIRED_V_L = 1.81e4` is kept as a named constant, and `make analysis-conductivity` prints the
cliff it implied beside the one that holds. Deleting it would let it drift back in.

`magnetic_reynolds` keeps **no default** for `v_l`. The caller must still say what expansion it
means; the change is that the module constant it reaches for is now solved rather than invented.

## Consequence 1 — the cliff moves down, and the band says it barely matters

Bisecting `Rm(T) = 1` on the continuous `sigma` at `rho = 0.32`, `x_K = 0.01`:

| `v L` [m^2/s] | 1.81e4 (retired) | 5.5e4 | **7.4e4 (stated)** | 9.7e4 |
| --- | ---: | ---: | ---: | ---: |
| cliff, solved | 2859 K | 2524 K | **2449 K** | 2386 K |
| cliff, log-interpolated from the six table rows | 2910 K | 2639 K | 2566 K | 2500 K |

The whole flown band spans 138 K of cliff, so nothing rests on the middle of it being exactly
right. The leak limit at 3800 K binds ~1350 K earlier at the stated leg, and 1276–1414 K earlier
across the band — **the field-grip argument gets stronger, not weaker**, because a faster expansion
holds the field to a lower temperature.

## Consequence 2 — the interpolated row is wrong, and is printed anyway

The second row above is kept executable as `interpolated_cliff_temperature`, purely for contrast.
`tab:seed_window` samples every 1000 K while `sigma` climbs 60x between its first two rows, and the
crossing lies **inside that first interval** — so log-interpolating the table and solving on the
interpolant overstates the crossing by +51 to +117 K. The tabulated values are not at fault: they
*are* the `sigma` the solver uses, sampled at the six temperatures the table reports. Asking six
points for a crossing inside their first interval is.

This is why a companion repository printed 2570 K, and printing both columns side by side is what
stops that recurring. **Do not quote the interpolated column.**

## Consequence 3 — the `Rm` and leak columns of `tab:seed_window` move with it

`Rm` is linear in `v L`, so the regenerated table's `Rm` column scales by 7.4e4/1.81e4 = **4.09x**
and its `leak = 1/Rm` column by the reciprocal:

| T [K] | 2000 | 3000 | 5000 | 8000 | 11 000 | 15 000 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `sigma` [S/m] (unchanged) | 1.056 | 63.52 | 492.6 | 734.6 | 2184 | 6993 |
| `Rm` at the retired 1.81e4 | 0.0240 | 1.445 | 11.20 | 16.71 | 49.67 | 159.0 |
| **`Rm` at 7.4e4** | **0.0982** | **5.907** | **45.81** | **68.31** | **203.1** | **650.2** |
| **leak = min(1/`Rm`, 1)** | **1.000** | **0.169** | **0.022** | **0.015** | **0.005** | **0.002** |

Anything downstream that quotes a leak fraction from this table is quoting a **smaller** leak than
before. The 2000 K row is still capped at 1.000: below the cliff the `1/Rm` bracket has broken down,
and reporting the breakdown of an approximation as a measurement is what the cap exists to prevent.

## Consequence 4 — ADR-0038's electron-temperature margin widens

ADR-0038 quotes Q-F's table for the cliff under a crude electron-temperature elevation and concludes
the cliff is not reached on any leg. Those figures are at the retired `v L`; at the stated one the
cliff falls further, so the margin grows and the verdict is unchanged:

| `t_e_offset` | retired 1.81e4 | **stated 7.4e4** | band (9.7e4 → 5.5e4) |
| --- | ---: | ---: | ---: |
| +500 K | 2344 K | **1946 K** | 1883–2019 K |
| +1000 K | 1841 K | **1446 K** | 1383–1518 K |

The +1000 K row is below `cliff_temperature`'s default `t_lo = 1500 K` everywhere except the very
bottom of the band, so calling it at the stated `v L` **raises** rather than extrapolating. That is
the intended behaviour, not a regression: the figures above were obtained by passing `t_lo = 600.0`
explicitly, which is the caller saying it knows it is asking about a colder gas than the seed window
covers.

## Status

Accepted 2026-08-26. Closes the `puffsat_impact_simulation` half of deferred item D9
(`docs/deferred_to_companion_repos.md`). Pinned by
`test_the_reference_v_l_is_the_solved_expansion_exit_not_a_back_solved_table`,
`test_the_cliff_at_the_stated_expansion_is_2450_k_and_the_band_barely_moves_it` and
`test_interpolating_the_six_tabulated_conductivities_overstates_the_cliff`.
