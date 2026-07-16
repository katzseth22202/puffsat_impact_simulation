# PuffSat Impact Simulation

Glossary for the per-collision study that computes the paper's fudge factor `f(v)` —
the momentum-transfer efficiency of a PuffSat gas pulse bouncing off the pusher plate.
Full design: [`puffsat_impact_sim_design.md`](puffsat_impact_sim_design.md).

## Language

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
