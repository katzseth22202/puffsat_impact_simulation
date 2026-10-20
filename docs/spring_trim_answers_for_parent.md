# Keeping the lossless spring on its cycle: answers (S48)

Answers 2026-10-20 to the parent's S48 (its `docs/spring_trim_asks_for_impact_sim.md`, commit
94f8808). **Written to be copied into the parent.** Evidence: Result 19 of
[`xray_collision_study.md`](xray_collision_study.md) and
[ADR-0070](adr/0070-spring-held-by-next-charge-k-trim.md); `make xray-collision-spring-trim`
(`spring_trim_*.csv`). **Status: proposed.**

Two levers are taken, as Seth set them (2026-10-20):
- **The next charge's k.** The craft deploys each polyethylene charge, so it can size it from the
  swing it has just measured.
- **The craft's conventional thrusters.** They take up the gap between the push expected and the
  push achieved. The budget is small and is not yet charged.

The answer: k does the per-pulse work and the thrusters clear the remainder. Both stay modest.

## Short answers

1. **The scatter is ~0.85% of J in heave, from tolerances we assumed.** The largest part is how
   well each sheet's mass is known on arrival (1% of mass gives 0.77% of J). Charge mass, standoff,
   aim loss, speed and plume drag add 0.01-0.3% each. The aim matters far more for *rocking*:
   10 cm off the axis kicks the floor's edge at 4.2% of the heave kick.
2. **The drift is faster than your estimate.** On the flown 2 pi/3 spring, a 1% scatter walks the
   plate 4 m past its stroke in a median 85 pulses (40 at 2%, 15 at 5%). You estimated 625-2500.
   The reason is that an error rings at the spring's free frequency, so 1 m/s of error is about
   0.24 m of swing, while the steady stroke is only a 2 pi/3 arc. Result 17's loss barely helps,
   and no seed survives the push.
3. **The trim is the next charge's k, every pulse.** k varies by 0.24 rms (6.5-7.7), 6.3 kg of a
   185 kg charge. That holds the plate to ~1 m past each end of its stroke (p95), so the bags need
   ~10 m, not 8. The controller also holds the summed impulse, so the thrusters see only
   0.16-0.45 m/s at the push's end: 61-173 kg of propellant at 300 s. It holds to 2% scatter. At
   5% the k = 5-9 we ran saturates.
4. **theta = pi is unstable.** It is not more forgiving: the gas column's stiffening puts its cycle
   in Hill's first instability tongue, and an error doubles every pulse. The k trim holds it only
   at <= 1% scatter and a lag of <= 1 pulse. The cycle is stable up to theta = 0.85 pi. **Fly
   5 pi/6 instead**: an 11.9 m column (against 18.6), held better than the flown spring.

## 1. The scatter (item 1)

These are new 2-D runs of the floor-wide 9.5 m sheet on the flown plate at k = 5-9 and standoff
1.5-3 m. Only the impulse inside the floor counts, since the skirt hangs on the frame. At fixed
sheet mass the impulse rises **3.3% per unit of k** (`d ln J / d ln k` = 0.23), and it falls 3.6%
per metre of standoff. The 1-sigma errors are assumed tolerances, not measurements.

| source | 1-sigma error (assumed) | heave, of J | rock, edge kick / heave kick |
|---|---|---|---|
| sheet mass on arrival | 1% | 0.77% | - |
| charge mass | 0.5% | 0.12% | - |
| standoff | 5 cm | 0.18% | - |
| aim offset | 10 cm per axis | 0.30% (and a 0.58% mean loss, taken out by the mean k) | **4.2%** |
| closing speed | 1e-4 | 0.01% | - |
| plume drag (sim ADR-0069, 3-5 m/s) | 5 m/s | 0.005% | - |
| sheet areal density | 5% per m^2 | (its total is the mass line) | 0.46% |
| charge areal density | 2% per m^2 | (its total is the mass line) | 0.06% |
| arrival time | 0.1 ms | a phase error of 1.2 cm | - |
| **total** | | **0.85%** | **4.2%** |

Rocking is `4 d / R` of the heave kick for an impulse centroid `d` off the axis of a uniform
floor. A 10 cm aim error is therefore five times the whole heave budget.

## 2. The drift, solved (items 2 and 4)

This is Result 17's adiabatic nitrogen column, flown for 2000 pulses with 100 seeds. Its 0.14-0.37%
loss is modelled as an equivalent drag. The table gives the median pulses until the plate strays
1 m, then 4 m, past its steady stroke (0.14% loss).

| spring | column | 1% | 2% | 5% | 10% |
|---|---|---|---|---|---|
| 2 pi/3 (flown) | 18.6 m | 26 / 85 | 8 / 40 | 4 / 15 | 2 / 5 |
| 5 pi/6 | 11.9 | 37 / 543 | 12 / 111 | 4 / 21 | 3 / 8 |
| pi | 8.3 | 6 / 7 | 5 / 6 | 3 / 5 | 3 / 4 |
| your estimate | | 2500 | 625 | 100 | 25 |

At theta = pi the one-period map has a real eigenvalue of -2.006. A linear spring at theta = pi
sits on the stability edge. The gas column stiffens as it compresses, and that pushes the
perturbation's phase past pi. The edge is at **theta = 0.85 pi** (an 11.5 m column on an 8.45 m
stroke, peak 1.39x the mean). The bags' loss does not change it.

## 3. The trim (item 3)

A discrete LQR on the one-period map works from three quantities: the swing's position and
velocity, read over the 0.5 s interval to 1 cm and 1 cm/s (assumed), and the summed impulse error.
It sets the next charge's k within the 5-9 the 2-D runs cover. The lag counts extra pulses between
committing a charge and its pulse. All rows are at 0.14% loss.

| spring | scatter | lag | extra excursion, p95 | k rms (1-99%) | thruster dv / propellant |
|---|---|---|---|---|---|
| 2 pi/3 | 0.85% | 0 | 1.02 m | 0.24 (6.54-7.66) | 0.45 m/s / 173 kg |
| 2 pi/3 | 0.85% | 2 | 1.15 | 0.22 | 0.26 / 99 |
| 2 pi/3 | 2% | 0-2 | 2.4-2.6 | 0.51-0.57 (5.8-8.5) | 0.55-1.07 / 210-409 |
| 2 pi/3 | 5% | 0 | lost (k at its limit on 15% of pulses) | | |
| **5 pi/6** | **0.85%** | **0 / 2** | **0.77 / 1.02** | **0.22 (6.58-7.62)** | **0.34 / 130; 0.16 / 61** |
| 5 pi/6 | 2% | 0 | 1.80 | 0.54 | 0.88 / 338 |
| pi | 0.85% | 0-1 | 0.83-0.84 | 0.49 (6.0-8.3) | 0.51-0.54 / 196-209 |
| pi | 0.85% at lag 2, or 2% | | lost | | |

The other options:
- **Thrusters alone cannot hold the spring.** To cancel each pulse's error at 1% they would spend
  117 MN s (40 t at 300 s) over the push. Left to take only the sum, they spend 3.3 MN s (1.1 t),
  but the spring is lost. With the summed impulse in the k controller, they need 61-173 kg.
- **Standoff** is a lever of about the same size (3.6% per metre), but it is not needed. To use
  it, the charge must be placed to the centimetre against a plate that swings +-4 m.
- **Arrival time** is not needed. At 2 pi/3 the next charge reaches the swing's phase as well as
  its amplitude (0.27% of J cancels a 10 cm phase error in two pulses). A 0.1 ms arrival error
  changes nothing.
- **A switched orifice** is not needed. The bags' own loss barely slows the drift, so an orifice
  would have to dissipate far more to help, and it adds a valve to a sealless bag.

The floor of ~1 m per 1% of scatter comes from the one-pulse lag: each pulse's own error rings
for a period before the next charge can answer it.

## 4. Rocking, which k cannot reach

On a uniform floor over bags spread under all of it, rocking obeys the heave's equation. It drifts
the same way, and it is unstable at theta = pi. k is symmetric and cannot touch it. Two answers:

- **A wedged charge.** Its areal density tilts across the disk (up to 30% here), which moves the
  impulse's centroid by `0.23 x wedge x R / 4`. On the same controller, the edge strays (p95, both
  axes) as follows. Every untrimmed case is lost.

| spring | aim 1 cm | 3 cm | 10 cm |
|---|---|---|---|
| 2 pi/3 | 0.84 m (wedge 2.7% rms) | 1.81 m (5.7%) | lost |
| 5 pi/6 | 0.64 | 1.38 | 30% lost |
| pi | 0.62 | 5% lost | lost |

- **Tilt guides.** Your absorber budget already carries "alignment guides". Guides that hold the
  tilt hand each pulse's torque to the craft, which has no spring, so only the sum matters. At
  10 cm aim, thruster couples at the plate's radius clear it with 6.9 MN s (2.4 t) over the push.
  The guides' loads and the craft's attitude are not sized here.

So the aim needs to be ~3 cm per axis with a wedged charge, or the plate needs tilt guides.

## Not modelled

- The pulse's taper as propellant burns. It moves the steady orbit slowly and predictably, and a
  feed-forward on k follows it (Seth: not needed at first order).
- k outside 5-9.
- The loss as thermal hysteresis rather than an equivalent drag.
- The 4-cell 1 MPa layout, whose rocking stiffness differs from its heave stiffness.

## Suggested changes (for the parent to make)

- `sec:charge_pusher_plume`, after the spring paragraph: "A lossless spring never forgets an
  impulse error: at 1% scatter the plate would leave its stroke within about a hundred pulses.
  The craft holds it by sizing each charge from the swing it measures between pulses. At an
  estimated 0.85% scatter the charge's k varies by 0.24 rms about 7.07, and the plate strays about
  1 m past each end of its stroke. The same trim holds the push's total impulse, leaving the
  thrusters about 0.5 m/s, some 170 kg of propellant \cite{Katz_puffsat_impact_sim_2026}."
- Same section, on the column: replace theta = pi as the gas-column fix with theta = 5 pi/6
  (11.9 m). "Stiffer than about 0.85 pi, the gas column's own stiffening makes the cycle unstable,
  doubling any error each pulse \cite{Katz_puffsat_impact_sim_2026}."
- `sec:charge_pusher_open`: add the aim. "Rocking cannot be trimmed by the charge's mass. An
  impulse 10 cm off the plate's axis rocks its edge at 4% of the heave kick. The plate needs ~3 cm
  aim with a charge whose areal density tilts across it, or alignment guides that hand the torque
  to the craft."
- The unpriced list: add the ~2 m of extra bag stroke and the trim's thruster propellant
  (~0.1-0.2 t a push).
