"""The walled thermal nozzle study (asks N9-N11 of the companion repository).

ADR-0050 is the scope. **Nothing here shares the magnetic nozzle's geometry or the `f(v)`
bag**: 3 m bore, no field, a 200-673 m^3 chamber, 25 kg at 75 km/s into methane. A module that
reads `expansion.BAG_RHO`, `expansion.AREA_RATIO_EXIT` or the graded 20-5 T field is reading the
wrong study.
"""
