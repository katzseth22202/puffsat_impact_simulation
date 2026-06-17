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
0 (stick) to 1 (elastic). `(1 − e_eff)` is the per-pulse momentum loss, decomposed by channel.
The thermophysics output (1D rad-hydro track).
_Avoid_: bounce factor, coefficient of restitution `e`

**Capture efficiency (`eta_capture`)**:
The geometric efficiency of the bounce relative to a perfectly-collimated 1D collision — the
fraction of axial momentum that lands and rebounds usefully rather than escaping sideways. Pure
geometry, set by plate radius, curvature, and cloud footprint; the perfectly-collimated 1D case
is the `eta_capture = 1` ceiling. The geometry output (2D Euler track).
_Avoid_: collection efficiency, catch fraction

**Incident momentum (`p_in`)**:
The axial momentum a pulse carries into the collision in the plate frame, `p_in = m_pulse · v`.
The normalizing reference for both `e_eff` and `f`.
