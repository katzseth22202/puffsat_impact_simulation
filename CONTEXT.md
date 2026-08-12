# PuffSat Impact Simulation

Glossary for two studies in this repository:

- the **per-collision study** that computes the paper's fudge factor `f(v)` — the
  momentum-transfer efficiency of a PuffSat gas pulse bouncing off the pusher plate.
  Full design: [`puffsat_impact_sim_design.md`](puffsat_impact_sim_design.md).
- the **tamped-nozzle study** that computes effective Isp for a tamped head-on collision.
  Full requirements: [`puffsat_tamper_isp_prd.md`](puffsat_tamper_isp_prd.md).

They share the vehicle and the pusher plate, not the deliverable.

## Language — per-collision study (`f(v)`)

**Fudge factor (`f`)**:
The axial momentum actually delivered to the plate by one pulse, as a fraction of the
theoretical maximum (full capture + perfect elastic bounce). `f = eta_capture · (1 + e_eff) / 2`,
ranging 0.5 (dead stick, full capture) to 1 (elastic, full capture). The single deliverable,
reported as `f(v)` across the velocity envelope.
_Avoid_: efficiency, coefficient of restitution (that is `e_eff`)

**Wall impulse (`J_wall`)**:
The physical momentum the gas transfers to the plate in one pulse — the time-integrated axial
wall force, `J_wall = ∫ P_wall(t)·A dt`. Equals `p_in · (1 + e_eff)` in 1D: the incoming
momentum plus the rebound. This is what the 1D solver measures directly; `e_eff` is its
normalized form, not a separately-measured quantity.
_Avoid_: momentum delivered (imprecise — say whether you mean `J_wall` or just `p_in`)

**Effective restitution (`e_eff`)**:
The rebound fraction of incident axial momentum that survives radiative, conductive, and
condensation losses: `e_eff = p_rebound / p_in`, measured as `J_wall / p_in − 1`. Ranges
0 (dead stick) up to the **bounce ceiling** (below) — `1` is only the idealized `M → 0` /
specular-reflection limit, not a value a finite gas slug reaches. `(1 − e_eff)` is the
per-pulse momentum loss, decomposed by channel. The thermophysics output (1D rad-hydro track).
_Avoid_: bounce factor, coefficient of restitution `e`

**Bounce ceiling (lossless gas-dynamic `e_eff`)**:
The maximum `e_eff` a re-expanding gas slug can return with **zero** radiative/conductive/
condensation losses — strictly `< 1`, set by `γ` and incident Mach `M`. The rebound is a
rarefaction fan with a velocity spread, so by Cauchy–Schwarz the coherent rebound momentum is
below the incident even in the lossless case. This is the true upper bound on `e_eff`, not 1;
the `f → 1` smoke test is therefore a bookkeeping + `M → 0` limit check, not a target (ADR-0001).
_Avoid_: "elastic limit" as a numeric target of 1.

**Capture efficiency (`eta_capture`)**:
The geometric efficiency of the bounce relative to a perfectly-collimated 1D collision — the
fraction of axial momentum that lands and rebounds usefully rather than escaping sideways. Pure
geometry, set by plate radius, curvature, and cloud footprint; the perfectly-collimated 1D case
is the `eta_capture = 1` ceiling. The geometry output (2D Euler track).
_Avoid_: collection efficiency, catch fraction

**Incident momentum (`p_in`)**:
The axial momentum a pulse carries into the collision in the plate frame, `p_in = m_pulse · v`.
The normalizing reference for both `e_eff` and `f`.

**Pulse shape**:
The geometry of the gas pulse at the moment of impact, at fixed pulse mass and speed: footprint
coverage `r_foot/R`, aspect ratio `L/D`, edge taper, and radial divergence. Shape changes at fixed
mass move `f` through both factors — `eta_capture` (2D geometry) and `e_eff` via the areal density
`Σ` (1D thermophysics).
_Avoid_: plume, cloud shape ("plume" suggests an engine exhaust; this is a delivered pulse)

**Shape box**:
The assumed dispersion region around the nominal pulse shape over which shape sensitivity is
assessed. An *assumption* standing in for real delivery dispersion, which is unquantified until
the deferred cloud-schedule study; claims made over the shape box must say so.

**Normalized shape sensitivity (`S`)**:
Per shape axis `x`, the relative response of the fudge factor: `S_x = (Δf/f)/(Δx/x)`, reported as
a max over the shape box. The quotable form of "slight shape change → slight impulse change";
`|S| ≲ 1` means a 1% shape error costs ≲ 1% impulse. The claim requires bounded `S` *and* no
cliff (no second-difference outlier surviving grid refinement).

**Plate radius (`R`)**:
The radius of the circular pusher plate — always a *radius*, so the plate's width/diameter is
`2R`. The canonical size variable: footprint coverage is the ratio `r_foot/R`, and impact density
scales as `ρ ∝ m_pulse/R²` through the Σ contract. "A 15 m-wide plate" therefore means `2R = 15 m`
(`R = 7.5 m`) — distinct from tripling `R` itself. State whether a plate dimension is `R` or `2R`
whenever it is not a ratio.
_Avoid_: "plate width/size" as a bare number (ambiguous between `R` and `2R`).

## Language — tamped-nozzle study (effective Isp)

**Pusher plate**:
The permanent structure on the vehicle that receives the impulse, with a hole at its vertex
through which the projectile passes. Always "plate", never "dish" — earlier tamper notes used
"dish" and it is the outlier.
_Avoid_: pusher dish, nozzle (the nozzle is the tamper + plate acting together, not the plate)

**Slug ratio (`k`)**:
Slug mass per **projectile** kg, per pulse — the carried mass the projectile buries itself in,
divided by the projectile's own mass. A continuous real, ≈ 7 at the bare-plate Isp optimum, not
a count. Stated in this direction always; the inverse (projectile per carried kg) reads as a
plausible but different quantity and has caused a real collision.
_Avoid_: "impactor-to-rocket-mass ratio" (ambiguous in direction), treating `k` as an integer

**Tamper ratio (`τ_t`)**:
Tamper mass per **slug** kg, `τ_t = m_t/m_s`. The subscript is required because optical depth is
`τ_opt`. This is a mass ratio, not an optical property.
_Avoid_: bare `τ` in the tamped-nozzle study

**Interlayer ratio (`μ`)**:
Interlayer mass per **slug** kg, `μ = m_int/m_s`. It is zero for the bare and vacuum-gap reference
cases, but nonzero for filled Arm D and must not be silently folded into `τ_t`.

**Hydrodynamic carried-mass ratio (`K`)**:
All nonprojectile carried mass represented in the near-field hydrodynamics per projectile kg,
`K = m_hydro/m_i`. When the only regions are slug, tamper, and interlayer,
`K = (m_s+m_t+m_int)/m_i = k(1+τ_t+μ)`; simulated spacers or supports add their own terms.
It equals the charged-mass ratio only when no other expendable carried mass is present.
_Avoid_: `K = k(1+τ_t)` when an interlayer is present; using `k` where `K` is meant

**Charged-mass ratio (`C`)**:
All expended mass carried by the vehicle per projectile kg, `C = m_charged/m_i`. Thus
`C = K + a_abl + a_other`, where `a_abl = m_abl/m_i` and `a_other` covers named expendables
outside the near-field model. Effective Isp uses `C`:
`Isp_eff = βw/(g₀C)`.

**Ejecta-mass ratio (`K_ej`)**:
Nonprojectile mass that exits with the event's ejecta per projectile kg. The thermodynamic
ceiling uses this mass: `j_max = J_max/m_i = w(√(1+K_ej) − 1)`. In a closed no-ablator
near-field calculation, `K_ej = K`; in Pass 2 it also includes whatever ablator mass the wall
model ejects across the system boundary.
This is distinct from `C`, which charges expended carried mass whether or not all of it exits.

**Tamper**:
Carried mass placed beyond the slug, on the far side from the plate, that turns part of the
away-going fireball back toward it. Fully vaporised — it works by **areal density**, which is
conserved through vaporisation, not by surviving. It is a **piston**, not a mirror: its
backward momentum is credited, while its entropy production is one part of the ceiling-shortfall
ledger (ADR-0030).
_Avoid_: reflector, mirror (both imply reflected-speed is the figure of merit — it is not)

**Ballistic capture fraction**:
`max[0, (1 − 1/√k)/2]` — the share of the isotropic, pressure-free ballistic fireball model
that outruns its own centre-of-mass recoil and reaches the plate. It is zero at `k ≤ 1`; this is
a limit of the ballistic model, not a claim that pressure-driven hydrodynamic thrust is exactly
zero there. Pure geometry, distinct from `eta_capture` in the `f(v)` study.
_Avoid_: conflating with `eta_capture` (that is a 2D/1D wall-impulse ratio, this is a
free-flight geometric fraction)

**Realization fraction**:
Net axial vehicle impulse as a share of the thermodynamic ceiling for the same ejecta mass,
`r_real = J/J_max = β/(√(1+K_ej) − 1)`. This is the hydrodynamic comparison metric that feeds the
effective-Isp deliverable and makes "spend the mass as tamper" and "spend it as slug" directly
comparable. The often-quoted 62.9% break-even value applies only to the reference candidate
`k = 7.06, τ_t = 1, μ = a_abl = 0` compared with the bare `k = 7.06` design; swept configurations
use their own charged-mass break-even condition, `β_candidate/C_candidate > β_reference/C_reference`.
Equivalently, the candidate's required realization fraction is
`β_reference C_candidate/[C_reference(√(1+K_ej_candidate)−1)]`.
_Avoid_: tamper multiplier, `Λ` (a ratio of two realized numbers, which hides the ceiling and
miscounts recoil as loss — ADR-0030)

**Projectile economy**:
Net axial vehicle impulse per **projectile** kg, `β·w = J/m_i`, where
`β = J/(m_i·w)` is dimensionless and `J` includes the incoming projectile momentum debit.
`β·w` has velocity units but is an impulse-per-mass normalisation, not the velocity of any
material or body. It is the infrastructure-throughput counterpart to effective Isp's
carried-mass economy. Reported beside Isp, never combined with it: converting projectiles to
payload-equivalent needs program economics this repository does not model.
_Avoid_: two-currency ledger (jargon, and it implied a combined figure of merit that is not
constructed)
