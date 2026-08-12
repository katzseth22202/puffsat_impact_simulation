"""Cold path for the **tamped-nozzle effective-Isp study** — a *different study* from `f(v)`.

Two studies share this repository (`CONTEXT.md`). The older one computes the paper's fudge factor
`f(v)`, the momentum-transfer efficiency of one gas pulse bouncing off the pusher plate; its
modules are the flat ones in `puffsat/` (`analysis`, `heavyplate`, `shape`, ...). **This package
is the other one**: effective specific impulse for a tamped head-on collision, specified by
[`puffsat_tamper_isp_prd.md`](../../../puffsat_tamper_isp_prd.md).

They share the vehicle, the pusher plate, the kernels, and the validation discipline — not the
deliverable, not the regime, and not the plate-side conventions. Where this study departs from a
decision recorded for the other one, PRD §12 says so and an ADR records it (ADR-0030..0033).
Keeping it in its own package is the point: importing across the boundary should be a deliberate,
visible act rather than an accident of a shared module.

Build entry points live under the `tamper-*` Makefile targets, never the bare `f(v)` ones.
"""
