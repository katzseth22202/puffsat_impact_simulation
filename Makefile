# PuffSat impact simulation — single build entry point (ADR-0018).
# Delegates to cargo (Rust hot path) and uv (Python cold path); the two meet only in data/.
#
# TWO STUDIES live here (CONTEXT.md). Unprefixed targets belong to the per-collision `f(v)` study
# (puffsat_impact_sim_design.md); `tamper-*` targets belong to the tamped-nozzle effective-Isp
# study (puffsat_tamper_isp_prd.md), whose block is at the bottom of this file.

PY := uv run python

.PHONY: tamper-ledger tamper-test
.PHONY: all smoke build test lint fmt clean tables sweep analysis sensitivity sweep-geometry-m40 analysis-lte tables-lowv sweep-lowv analysis-lowv sweep-transitional analysis-transitional sweep-geometry analysis-geometry analysis-survivability analysis-margin sweep-ablating analysis-ablating sweep-frozen-probe tables-frozen sweep-frozen analysis-frozen tables-jupiter sweep-jupiter analysis-jupiter sweep-frozen-probe-jupiter tables-frozen-jupiter sweep-frozen-jupiter analysis-frozen-jupiter fetch-tops sweep-heavyplate analysis-heavyplate analysis-structure-heavyplate sweep-frozen-probe-heavyplate tables-frozen-heavyplate sweep-frozen-heavyplate analysis-frozen-heavyplate sweep-shape analysis-shape sweep-frozen-probe-shape tables-frozen-shape sweep-frozen-shape analysis-frozen-shape

all: smoke

## smoke: boundary round-trip plumbing test (Python -> JSON -> Rust -> JSONL -> Python)
smoke: build
	@mkdir -p data/tables data/results
	@rm -f data/results/smoke.jsonl
	$(PY) python/puffsat/smoke.py write
	cargo run --quiet -p smoke
	$(PY) python/puffsat/smoke.py check

## build: compile the Rust workspace
build:
	cargo build

## test: run all tests (cargo + pytest)
test:
	cargo test
	uv run pytest

## lint: ruff + mypy + clippy + fmt checks (CI gate)
lint:
	uv run ruff check python
	uv run ruff format --check python
	uv run mypy
	cargo clippy --all-targets --all-features
	cargo fmt --all -- --check

## fmt: auto-format Python and Rust
fmt:
	uv run ruff format python
	uv run ruff check --fix python
	cargo fmt --all

## clean: remove build artifacts and smoke outputs
clean:
	cargo clean
	rm -f data/tables/smoke.json data/results/smoke.jsonl

# --- Physics pipeline (stubs filled in at the corresponding build rungs) ---
## tables: generate the water EOS/opacity table (rung B) -> data/tables/water.json
tables: data/tables/water.json

data/tables/water.json: python/puffsat/eos_water.py python/puffsat/tables.py
	@mkdir -p data/tables
	PYTHONPATH=python $(PY) -m puffsat.tables

## tables-lowv: generate the Rung C cool-gas two-phase table (CoolProp) -> data/tables/water_lowv.json
tables-lowv: data/tables/water_lowv.json

data/tables/water_lowv.json: python/puffsat/eos_cool.py python/puffsat/tables.py
	@mkdir -p data/tables
	PYTHONPATH=python uv run --extra sci python -m puffsat.tables --lowv

## sweep: run the 16 km/s e_eff(rho) sweep (rung B) -> data/results/sweep.jsonl; depends on tables
sweep: data/results/sweep.jsonl

data/results/sweep.jsonl: data/tables/water.json $(wildcard crates/sweep/src/*.rs) $(wildcard crates/hydro1d/src/*.rs)
	@mkdir -p data/results
	cargo run --release -p sweep

## analysis: frontier extraction + plots (rung B) -> data/results/frontier.csv + figures; depends on sweep
analysis: data/results/frontier.csv

data/results/frontier.csv: data/results/sweep.jsonl python/puffsat/analysis.py
	PYTHONPATH=python uv run --extra sci python -m puffsat.analysis

## sweep-lowv: 3.2 km/s condensing e_eff(rho) anchor (rung C) -> data/results/sweep_lowv.jsonl
sweep-lowv: data/results/sweep_lowv.jsonl

data/results/sweep_lowv.jsonl: data/tables/water_lowv.json $(wildcard crates/sweep/src/*.rs) $(wildcard crates/hydro1d/src/*.rs)
	@mkdir -p data/results
	cargo run --release -p sweep -- --lowv

## analysis-lowv: frontier + figures for the low-v anchor -> data/results/frontier_lowv.csv; depends on sweep-lowv
analysis-lowv: data/results/frontier_lowv.csv

data/results/frontier_lowv.csv: data/results/sweep_lowv.jsonl python/puffsat/analysis.py
	PYTHONPATH=python uv run --extra sci python -m puffsat.analysis \
		--sweep data/results/sweep_lowv.jsonl --summary data/results/frontier_lowv.csv --tag lowv_

## sweep-transitional: transitional-anchor e_eff(v) sweep (ADR-0012) over V_GRID x RHO_GRID with the
## high-v table; emits the EOS-only and radiation-on curves into two files in one run. Depends on tables.
sweep-transitional: data/results/sweep_transitional_eos.jsonl

data/results/sweep_transitional_eos.jsonl: data/tables/water.json $(wildcard crates/sweep/src/*.rs) $(wildcard crates/hydro1d/src/*.rs)
	@mkdir -p data/results
	cargo run --release -p sweep -- --transitional

## analysis-transitional: e_eff(v) frontier + EOS-vs-rad overlay + dip locator (ADR-0012) ->
## data/results/frontier_transitional.csv + figure; depends on sweep-transitional
analysis-transitional: data/results/frontier_transitional.csv

data/results/frontier_transitional.csv: data/results/sweep_transitional_eos.jsonl data/results/sweep_transitional_rad.jsonl python/puffsat/analysis.py
	PYTHONPATH=python uv run --extra sci python -m puffsat.analysis --axis v

## sweep-geometry: 2D eta_capture(curvature x L/D x r_foot/R) sweep (Rung D follow-on) ->
## data/results/sweep_geometry.jsonl. Radiation-free (euler2d, effective-gamma), so no table needed.
sweep-geometry: data/results/sweep_geometry.jsonl

data/results/sweep_geometry.jsonl: $(wildcard crates/sweep/src/*.rs) $(wildcard crates/euler2d/src/*.rs)
	@mkdir -p data/results
	cargo run --release -p sweep -- --geometry

## analysis-lte: McWhirter LTE-validity check at the probe turnaround states (Q5) ->
## data/results/lte_validity.csv. Both the Saha EOS and the LTE TOPS opacity assume collisionally
## controlled level populations; this reports where the probed states sit relative to that.
analysis-lte: data/results/lte_validity.csv

data/results/lte_validity.csv: data/results/frozen_probe_heavyplate.jsonl data/results/frozen_probe_jupiter.jsonl python/puffsat/lte.py python/puffsat/eos_water.py
	PYTHONPATH=python $(PY) -m puffsat.lte

## sweep-geometry-m40: the same eta_capture sweep at M = 40 (the high-v scenarios' Mach anchor) ->
## data/results/sweep_geometry_m40.jsonl. Five frontier CSVs consume it (jupiter, frozen-jupiter,
## heavyplate, structure-heavyplate, frozen-heavyplate) as their eta_capture source.
sweep-geometry-m40: data/results/sweep_geometry_m40.jsonl

data/results/sweep_geometry_m40.jsonl: $(wildcard crates/sweep/src/*.rs) $(wildcard crates/euler2d/src/*.rs)
	@mkdir -p data/results
	cargo run --release -p sweep -- --geometry-m40

## analysis-geometry: f = eta_capture*(1+e_eff)/2 reconciliation + eta/f figures (Rung D follow-on)
## -> data/results/frontier_geometry.csv; depends on sweep-geometry.
analysis-geometry: data/results/frontier_geometry.csv

data/results/frontier_geometry.csv: data/results/sweep_geometry.jsonl python/puffsat/analysis.py
	PYTHONPATH=python uv run --extra sci python -m puffsat.analysis --axis geometry

## analysis-survivability: peak facesheet pressure vs P_limit + the survivability-resolved f frontier
## (Rung S) -> data/results/frontier_survivability.csv. Resolves each geometry case to a peak
## stagnation pressure via the Sigma contract (c_stag from the 1D sweeps) and classifies it against
## the SiC+Ti limits (ADR-0010/0011); reuses existing results, no new sweep.
analysis-survivability: data/results/frontier_survivability.csv

data/results/frontier_survivability.csv: data/results/sweep_geometry.jsonl data/results/sweep.jsonl data/results/sweep_transitional_eos.jsonl python/puffsat/analysis.py
	PYTHONPATH=python uv run --extra sci python -m puffsat.analysis --axis survivability

## analysis-margin: closed-form f-margin exploration (design §7, ADR-0010 amendment) ->
## data/results/frontier_margin.csv. Rescales the survivability frontier over the (plate radius R,
## pulse mass m) grid (peak ∝ m/R³, eta_capture scale-invariant) to map how much survivable f a
## wider plate / smaller pulse buys above the passing baseline. Reuses existing results, no sweep.
analysis-margin: data/results/frontier_margin.csv

data/results/frontier_margin.csv: data/results/sweep_geometry.jsonl data/results/sweep.jsonl data/results/sweep_transitional_eos.jsonl python/puffsat/analysis.py
	PYTHONPATH=python uv run --extra sci python -m puffsat.analysis --axis margin

## sweep-ablating: ablating-wall recovery sweep (Rung E, ADR-0014) over (v x rho x opacity-scale x
## Q*) -> data/results/sweep_ablating.jsonl. Rigid floor vs shielding+injection ablating wall;
## depends on the high-v table (opacity-scaled in-process).
sweep-ablating: data/results/sweep_ablating.jsonl

data/results/sweep_ablating.jsonl: data/tables/water.json $(wildcard crates/sweep/src/*.rs) $(wildcard crates/hydro1d/src/*.rs)
	cargo build --release -p sweep
	cargo run --release -p sweep -- --ablating

## analysis-ablating: tau-bracketed e_eff recovery + the 16 km/s f>=0.8-at-a-survivable-shape call
## (Rung E, ADR-0014/0009) -> data/results/frontier_ablating.csv + figure; folds in the geometry +
## survivability results (no new sweep). Depends on sweep-ablating.
analysis-ablating: data/results/frontier_ablating.csv

data/results/frontier_ablating.csv: data/results/sweep_ablating.jsonl data/results/sweep_geometry.jsonl data/results/sweep.jsonl python/puffsat/analysis.py
	PYTHONPATH=python uv run --extra sci python -m puffsat.analysis --axis ablating

## sweep-frozen-probe: turnaround-state probe for the frozen-recombination check (audit finding 3)
## -> data/results/frozen_probe.jsonl. EOS-only transitional grid; records each case's mass-weighted
## (rho*, T*) at global momentum zero.
sweep-frozen-probe: data/results/frozen_probe.jsonl

data/results/frozen_probe.jsonl: data/tables/water.json $(wildcard crates/sweep/src/*.rs) $(wildcard crates/hydro1d/src/*.rs)
	@mkdir -p data/results
	cargo run --release -p sweep -- --frozen-probe

## tables-frozen: per-case frozen-composition tables (sudden-freeze rebound EOS) + the pure-H2O
## no-chemistry bracket -> data/tables/frozen/. Depends on the probe.
tables-frozen: data/tables/frozen/h2o.json

data/tables/frozen/h2o.json: data/results/frozen_probe.jsonl python/puffsat/eos_water.py python/puffsat/tables.py
	PYTHONPATH=python $(PY) -m puffsat.tables --frozen-from-probe data/results/frozen_probe.jsonl

## sweep-frozen: the three-curve frozen-recombination bounding sweep (equilibrium vs
## freeze-after-the-plate vs freeze-before-the-plate) -> data/results/sweep_frozen.jsonl
sweep-frozen: data/results/sweep_frozen.jsonl

data/results/sweep_frozen.jsonl: data/tables/frozen/h2o.json data/tables/water.json $(wildcard crates/sweep/src/*.rs) $(wildcard crates/hydro1d/src/*.rs)
	cargo run --release -p sweep -- --frozen

## analysis-frozen: e_eff(v) freeze-timing bracket overlay + dip impact on f ->
## data/results/frontier_frozen.csv + figure; depends on sweep-frozen.
analysis-frozen: data/results/frontier_frozen.csv

data/results/frontier_frozen.csv: data/results/sweep_frozen.jsonl python/puffsat/analysis.py
	PYTHONPATH=python uv run --extra sci python -m puffsat.analysis --axis frozen

## tables-jupiter: extended-grid table for the 69 km/s scenario (multi-stage O Saha ladder,
## T to 1.2e6 K) -> data/tables/water_jupiter.json. Overlays the real TOPS/OPLIB opacity when
## the pull (data/tables/tops/, see fetch-tops) is present; interim Kramers bracket otherwise.
TOPS_PULL := data/tables/tops/tops_water_gray.html
tables-jupiter: data/tables/water_jupiter.json

data/tables/water_jupiter.json: python/puffsat/eos_water.py python/puffsat/tables.py python/puffsat/tops.py $(wildcard $(TOPS_PULL))
	@mkdir -p data/tables
	PYTHONPATH=python $(PY) -m puffsat.tables --jupiter $(if $(wildcard $(TOPS_PULL)),--tops $(TOPS_PULL),)

## fetch-tops: re-pull the TOPS/OPLIB water gray opacities (network; two-stage web form) ->
## data/tables/tops/tops_water_gray.html. The saved HTML is the citable provenance artifact.
fetch-tops:
	PYTHONPATH=python uv run --extra fetch python -m puffsat.fetch_tops

## sweep-jupiter: 69 km/s (rho x length x opacity-scale) coupled-bounce grid ->
## data/results/sweep_jupiter.jsonl; depends on tables-jupiter
sweep-jupiter: data/results/sweep_jupiter.jsonl

data/results/sweep_jupiter.jsonl: data/tables/water_jupiter.json $(wildcard crates/sweep/src/*.rs) $(wildcard crates/hydro1d/src/*.rs)
	@mkdir -p data/results
	cargo run --release -p sweep -- --jupiter

## analysis-jupiter: plate sizing + survivable-f frontier for the 69 km/s scenario ->
## data/results/frontier_jupiter.csv; depends on sweep-jupiter (+ the M=40 geometry sweep)
analysis-jupiter: data/results/frontier_jupiter.csv

data/results/frontier_jupiter.csv: data/results/sweep_jupiter.jsonl data/results/sweep_geometry_m40.jsonl python/puffsat/jupiter.py
	PYTHONPATH=python uv run --extra sci python -m puffsat.jupiter

## sweep-frozen-probe-jupiter: turnaround-state probe for the 69 km/s freeze-timing bracket
## (ADR-0026 instrument at the L=12 m realistic-cloud anchor) -> data/results/frozen_probe_jupiter.jsonl.
## EOS-only on the extended-grid Jupiter table; records each JUP_RHO case's (rho*, T*) at turnaround.
sweep-frozen-probe-jupiter: data/results/frozen_probe_jupiter.jsonl

data/results/frozen_probe_jupiter.jsonl: data/tables/water_jupiter.json $(wildcard crates/sweep/src/*.rs) $(wildcard crates/hydro1d/src/*.rs)
	@mkdir -p data/results
	cargo run --release -p sweep -- --frozen-probe-jupiter

## tables-frozen-jupiter: per-case frozen-composition tables (+ pure-H2O bracket) on the extended
## Jupiter grid (T to 1.2e6 K) -> data/tables/frozen_jupiter/. Depends on the Jupiter probe.
tables-frozen-jupiter: data/tables/frozen_jupiter/h2o.json

data/tables/frozen_jupiter/h2o.json: data/results/frozen_probe_jupiter.jsonl python/puffsat/eos_water.py python/puffsat/tables.py
	PYTHONPATH=python $(PY) -m puffsat.tables --frozen-from-probe data/results/frozen_probe_jupiter.jsonl --jupiter

## sweep-frozen-jupiter: three-curve freeze-timing bracket at 69 km/s (equilibrium vs
## freeze-after-the-plate vs freeze-before-the-plate) -> data/results/sweep_frozen_jupiter.jsonl
sweep-frozen-jupiter: data/results/sweep_frozen_jupiter.jsonl

data/results/sweep_frozen_jupiter.jsonl: data/tables/frozen_jupiter/h2o.json data/tables/water_jupiter.json $(wildcard crates/sweep/src/*.rs) $(wildcard crates/hydro1d/src/*.rs)
	cargo run --release -p sweep -- --frozen-jupiter

## analysis-frozen-jupiter: translate the 69 km/s EOS-only e_eff freeze bracket onto the survivable
## f -> data/results/frontier_frozen_jupiter.csv; depends on sweep-frozen-jupiter (+ the coupled
## sweep and M=40 geometry for the headline design point).
analysis-frozen-jupiter: data/results/frontier_frozen_jupiter.csv

data/results/frontier_frozen_jupiter.csv: data/results/sweep_frozen_jupiter.jsonl data/results/sweep_jupiter.jsonl data/results/sweep_geometry_m40.jsonl python/puffsat/jupiter.py
	PYTHONPATH=python uv run --extra sci python -m puffsat.jupiter --frozen

# --- Heavy-plate 16-28 km/s special scenario (design §12.1, ADR-0027): reuses the Jupiter table ---
## sweep-heavyplate: heavy-plate coupled-bounce grid (v x rho headline at fixed L + L-sensitivity
## spot rows + opacity tau-check) on the reused Jupiter extended-grid table ->
## data/results/sweep_heavyplate.jsonl; depends on tables-jupiter.
sweep-heavyplate: data/results/sweep_heavyplate.jsonl

data/results/sweep_heavyplate.jsonl: data/tables/water_jupiter.json $(wildcard crates/sweep/src/*.rs) $(wildcard crates/hydro1d/src/*.rs)
	@mkdir -p data/results
	cargo run --release -p sweep -- --heavyplate

## analysis-heavyplate: f(v) + facesheet-survivability frontier at the pinned 30 m / <=40 t plate ->
## data/results/frontier_heavyplate.csv + f(v) figure; depends on sweep-heavyplate (+ M=40 geometry).
analysis-heavyplate: data/results/frontier_heavyplate.csv

data/results/frontier_heavyplate.csv: data/results/sweep_heavyplate.jsonl data/results/sweep_geometry_m40.jsonl python/puffsat/heavyplate.py
	PYTHONPATH=python uv run --extra sci python -m puffsat.heavyplate

## analysis-structure-heavyplate: ADR-0027 closed-form whole-plate structural bound (rigid-during-
## pulse / f-validity, areal-impulse membrane, SiC-Ti spall) at the survivable design points ->
## data/results/frontier_structure_heavyplate.csv; depends on sweep-heavyplate (+ M=40 geometry).
analysis-structure-heavyplate: data/results/frontier_structure_heavyplate.csv

data/results/frontier_structure_heavyplate.csv: data/results/sweep_heavyplate.jsonl data/results/sweep_geometry_m40.jsonl python/puffsat/structure.py python/puffsat/heavyplate.py
	PYTHONPATH=python $(PY) -m puffsat.structure

## sweep-frozen-probe-heavyplate: turnaround-state probe for the heavy-plate freeze-timing bracket
## (ADR-0026, at the 16/22/28 km/s anchors) -> data/results/frozen_probe_heavyplate.jsonl.
sweep-frozen-probe-heavyplate: data/results/frozen_probe_heavyplate.jsonl

data/results/frozen_probe_heavyplate.jsonl: data/tables/water_jupiter.json $(wildcard crates/sweep/src/*.rs) $(wildcard crates/hydro1d/src/*.rs)
	@mkdir -p data/results
	cargo run --release -p sweep -- --frozen-probe-heavyplate

## tables-frozen-heavyplate: per-case frozen-composition tables (+ pure-H2O bracket) on the extended
## Jupiter grid -> data/tables/frozen_heavyplate/. Depends on the heavy-plate probe.
tables-frozen-heavyplate: data/tables/frozen_heavyplate/h2o.json

data/tables/frozen_heavyplate/h2o.json: data/results/frozen_probe_heavyplate.jsonl python/puffsat/eos_water.py python/puffsat/tables.py
	PYTHONPATH=python $(PY) -m puffsat.tables --frozen-from-probe data/results/frozen_probe_heavyplate.jsonl --outdir data/tables/frozen_heavyplate --jupiter

## sweep-frozen-heavyplate: three-curve freeze-timing bracket at 16/22/28 km/s (equilibrium vs
## freeze-after-the-plate vs freeze-before-the-plate) -> data/results/sweep_frozen_heavyplate.jsonl
sweep-frozen-heavyplate: data/results/sweep_frozen_heavyplate.jsonl

data/results/sweep_frozen_heavyplate.jsonl: data/tables/frozen_heavyplate/h2o.json data/tables/water_jupiter.json $(wildcard crates/sweep/src/*.rs) $(wildcard crates/hydro1d/src/*.rs)
	cargo run --release -p sweep -- --frozen-heavyplate

## analysis-frozen-heavyplate: translate the 16-28 km/s EOS-only e_eff freeze bracket onto the
## survivable f -> data/results/frontier_frozen_heavyplate.csv; depends on sweep-frozen-heavyplate
## (+ the coupled sweep and M=40 geometry for the per-anchor design points).
analysis-frozen-heavyplate: data/results/frontier_frozen_heavyplate.csv

data/results/frontier_frozen_heavyplate.csv: data/results/sweep_frozen_heavyplate.jsonl data/results/sweep_heavyplate.jsonl data/results/sweep_geometry_m40.jsonl python/puffsat/heavyplate.py
	PYTHONPATH=python uv run --extra sci python -m puffsat.heavyplate --frozen

# --- Pulse-shape sensitivity study (design §13, ADR-0028) ---
## sweep-shape: raw f(shape) inputs at the fixed baseline design — the fixed-grid 2D shape box
## (+ refined noise-floor repeats) and the fresh Sigma-contract 1D e_eff runs -> two JSONLs;
## depends on tables.
sweep-shape: data/results/sweep_shape_geometry.jsonl

data/results/sweep_shape_geometry.jsonl: data/tables/water.json $(wildcard crates/sweep/src/*.rs) $(wildcard crates/euler2d/src/*.rs) $(wildcard crates/hydro1d/src/*.rs)
	@mkdir -p data/results
	cargo run --release -p sweep -- --shape

## analysis-shape: assemble f over the shape box, compute S per axis, run the cliff detector, the
## Sigma-profile bound, and the survivability margin -> data/results/shape_sensitivity.csv + .png;
## depends on sweep-shape.
analysis-shape: data/results/shape_sensitivity.csv

data/results/shape_sensitivity.csv: data/results/sweep_shape_geometry.jsonl python/puffsat/shape.py
	PYTHONPATH=python uv run --extra sci python -m puffsat.shape

## sweep-frozen-probe-shape: turnaround-state probe for the three-point dip-anchor frozen
## spot-check (design §13, ADR-0026 instrument) -> data/results/frozen_probe_shape.jsonl.
sweep-frozen-probe-shape: data/results/frozen_probe_shape.jsonl

data/results/frozen_probe_shape.jsonl: data/tables/water.json $(wildcard crates/sweep/src/*.rs) $(wildcard crates/hydro1d/src/*.rs)
	@mkdir -p data/results
	cargo run --release -p sweep -- --frozen-probe-shape

## tables-frozen-shape: per-case frozen-composition tables (+ pure-H2O bracket) for the shape
## spot-check -> data/tables/frozen_shape/. Depends on the shape probe.
tables-frozen-shape: data/tables/frozen_shape/h2o.json

data/tables/frozen_shape/h2o.json: data/results/frozen_probe_shape.jsonl python/puffsat/eos_water.py python/puffsat/tables.py
	PYTHONPATH=python $(PY) -m puffsat.tables --frozen-from-probe data/results/frozen_probe_shape.jsonl --outdir data/tables/frozen_shape

## sweep-frozen-shape: three-curve freeze-timing spot-check at the dip anchor's three Sigma points
## -> data/results/sweep_frozen_shape.jsonl
sweep-frozen-shape: data/results/sweep_frozen_shape.jsonl

data/results/sweep_frozen_shape.jsonl: data/tables/frozen_shape/h2o.json data/tables/water.json $(wildcard crates/sweep/src/*.rs) $(wildcard crates/hydro1d/src/*.rs)
	cargo run --release -p sweep -- --frozen-shape

## analysis-frozen-shape: the frozen-vs-equilibrium slope comparison across the Sigma box
## (design §13 exit criterion) -> data/results/shape_frozen_spotcheck.csv
analysis-frozen-shape: data/results/shape_frozen_spotcheck.csv

data/results/shape_frozen_spotcheck.csv: data/results/sweep_frozen_shape.jsonl python/puffsat/shape.py
	PYTHONPATH=python uv run --extra sci python -m puffsat.shape --frozen

## sensitivity: opacity-insensitivity scan (rung B, B5d-3) — sweep at 0.1x/1x/10x opacity, show
## e_eff barely moves. Builds the release sweep first; writes data/results/opacity_scan/.
sensitivity:
	cargo build --release -p sweep
	PYTHONPATH=python uv run --extra sci python -m puffsat.sensitivity

# =================================================================================================
# TAMPED-NOZZLE EFFECTIVE-Isp STUDY  —  puffsat_tamper_isp_prd.md
#
# A SEPARATE STUDY from everything above. Everything above computes the paper's fudge factor `f(v)`
# (puffsat_impact_sim_design.md); the targets below compute effective specific impulse for a tamped
# head-on collision. They share the vehicle, the pusher plate, the kernels, and the validation
# discipline — not the deliverable, not the regime, and not the plate-side conventions (PRD §12
# lists every departure; ADR-0030..0033 record them).
#
# Conventions for this block, so the two studies cannot be run into each other by accident:
#   * every target is prefixed `tamper-`;
#   * every artifact lands under data/results/tamper/ or data/tables/tamper/;
#   * the Python lives in its own package, python/puffsat/tamper/, never the flat f(v) modules.
# Rungs are PRD §10. Later rungs reuse crates/hydro1d and crates/euler2d, but always through
# tamper-specific sweep modes and tamper-specific tables.
# =================================================================================================

## tamper-ledger: Rung 0 — the analytic reference ledger (PRD §10). The single cold-path calculator
## that owns every closed-form number the PRD quotes, so the PRD, the ADRs, and the analysis cannot
## drift apart; downstream rungs quote closed-form figures only from here. Pure algebra, stdlib
## only (no `sci` extra, no tables, no kernel) -> data/results/tamper/ledger_*.csv
tamper-ledger: data/results/tamper/ledger_anchors.csv

data/results/tamper/ledger_anchors.csv: python/puffsat/tamper/ledger.py
	@mkdir -p data/results/tamper
	PYTHONPATH=python $(PY) -m puffsat.tamper.ledger

## tamper-test: the tamped-nozzle study's tests alone (analytic anchors + invariants, PRD §8)
tamper-test:
	uv run pytest python/tests/test_tamper_ledger.py
