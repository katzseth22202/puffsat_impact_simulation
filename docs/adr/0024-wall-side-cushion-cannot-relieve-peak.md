# Wall-side cushions cannot relieve peak facesheet pressure: a buffer that works is just a longer cloud

A natural reviewer question — the pressure cousin of "why not deepen the dish?" (ADR-0021) — is
*"why not cushion the impact?"*: interpose a thin low-pressure buffer gas (a ~0.01 atm balloon) or a
thicker ablative spray between the incoming cloud and the plate, so the gas decelerates gradually and a
high-`eta_capture` flat disk could survive the survivability frontier (§7, ADR-0010). This ADR records
why that is foreclosed by momentum conservation, so it is logged as considered-and-rejected rather than
relitigated — and so the ablator's role (ADR-0014) is not mistaken for pressure relief.

**The governing principle: you cannot remove the impulse, only spread it in time.** Peak-pressure
relief means stretching the load out in time, because the time-integral of pressure *is* the delivered
momentum — `∫P dt ≈ (1+e_eff)·ρ·v·L` — which is fixed: it is the thrust the plate exists to receive.
The peak `≈ 2·ρv²` (the cloud's ram pressure recompressed at the wall, verified physical from the 1D
`peak_wall_force`, ADR-0010) is reached when the cloud is brought to rest. A cushion can lower that peak
*only* if it makes the stagnation take longer.

**Column density is the criterion, and a 0.01 atm balloon fails it by ~100x.** Deceleration is a
momentum-exchange problem, so what matters is areal mass (column density `ρ·L`), not pressure. The
incoming cloud's column density is the `Σ = m/(π r_foot²)` contract (ADR-0003): `Σ ≈ 1.3 kg/m²` at
`m = 25 kg`, `r_foot/R = 0.5`. A 0.01 atm balloon (water vapor, ~300 K, `ρ ≈ 0.007 kg/m³`) over a 1 m
standoff is `≈ 0.007 kg/m²` — about **180x lighter** than the cloud; even a 5 m-thick balloon is ~40x
short. A cushion ~100x lighter than the projectile is transparent to it: the cloud drives a shock
through the thin gas, sweeps it up, and arrives at the wall at essentially full `v`, stagnating at the
full `≈ 2·ρv²`. The shock-compressed buffer adds ~kPa against hundreds of MPa from the cloud —
negligible.

**A cushion dense enough to matter is just "a longer cloud."** To match the cloud's column density at a
1 m standoff needs `ρ ≈ 1.3 kg/m³` — roughly **1 atm**, ~25 kg of gas (as massive as the projectile
itself), replenished every pulse. And a graded gas column the cloud merges into is *physically
identical* to making the cloud itself longer (higher `L/D`). So a tamper that works recreates, at a
per-pulse mass/logistics cost, the exact effect the geometric lever gives for free — which is why the
real lever is to stretch the actual cloud, not to add a buffer at the wall.

**The ablator is therefore not a pressure device either** (clarifies ADR-0014). The ablative
spray/vapor and the thin balloon are both thin near-wall layers that remove none of the incoming `ρv²`
momentum flux. The ablator's job is wall-loss and oxidation: vapor-shielding to cut radiative/conductive
wall losses (raising `e_eff`) and SiO2 passivation against atomic O — a loss-recovery and survival
device, never peak-pressure relief.

## Consequence

Peak-pressure relief is bought only by spreading the cloud's *own* arrival in time — the geometric
`L/D` / footprint / shallow-concave trade of the survivability frontier (§7, ADR-0010/0021). No
wall-side cushion (buffer gas, thicker ablator) is a survivability lever. This joins the two other
foreclosed peak-relief mechanisms so the set is complete: the **bulk-recoil shock absorber** (wrong
timescale, ~10³x — bulk recoil acts on ms-s, the facesheet stress wave on µs; ADR-0011) and **facesheet
material strength** (low-leverage on `e_eff`; the binding mode is gas-imposed compression no material
reduces; §5, ADR-0010). All three point the same way: **survivability is a cloud-geometry property, not
a wall property.**

## Considered Options

- **Thin buffer-gas balloon (~0.01 atm) in front of the plate.** Rejected: column density ~100x below
  the cloud's, so it is transparent to the cloud and the peak stays at `≈ 2·ρv²`.
- **Dense/thick tamper (~1 atm, comparable column density).** Rejected: physically identical to a longer
  cloud (higher `L/D`), but with a ~25 kg per-pulse gas mass + replenishment cost — strictly worse than
  stretching the real cloud, which the cloud schedule can do for free.
- **Rely on the ablative spray to cushion the peak.** Rejected: the ablator is a wall-loss/oxidation
  device (ADR-0014); it is a thin near-wall layer that removes none of the incoming momentum flux, so it
  cannot lower `≈ 2·ρv²`.
