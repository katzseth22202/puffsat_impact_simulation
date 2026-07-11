# PuffSat impact simulation — single build entry point (ADR-0018).
# Delegates to cargo (Rust hot path) and uv (Python cold path); the two meet only in data/.

PY := uv run python

.PHONY: all smoke build test lint fmt clean tables sweep analysis sensitivity tables-lowv sweep-lowv analysis-lowv sweep-transitional analysis-transitional sweep-geometry analysis-geometry analysis-survivability analysis-margin sweep-ablating analysis-ablating sweep-frozen-probe tables-frozen sweep-frozen analysis-frozen tables-jupiter sweep-jupiter analysis-jupiter sweep-frozen-probe-jupiter tables-frozen-jupiter sweep-frozen-jupiter analysis-frozen-jupiter fetch-tops

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

## sensitivity: opacity-insensitivity scan (rung B, B5d-3) — sweep at 0.1x/1x/10x opacity, show
## e_eff barely moves. Builds the release sweep first; writes data/results/opacity_scan/.
sensitivity:
	cargo build --release -p sweep
	PYTHONPATH=python uv run --extra sci python -m puffsat.sensitivity
