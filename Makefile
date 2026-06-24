# PuffSat impact simulation — single build entry point (ADR-0018).
# Delegates to cargo (Rust hot path) and uv (Python cold path); the two meet only in data/.

PY := uv run python

.PHONY: all smoke build test lint fmt clean tables sweep analysis sensitivity tables-lowv sweep-lowv analysis-lowv sweep-transitional analysis-transitional sweep-geometry analysis-geometry analysis-survivability

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

## sensitivity: opacity-insensitivity scan (rung B, B5d-3) — sweep at 0.1x/1x/10x opacity, show
## e_eff barely moves. Builds the release sweep first; writes data/results/opacity_scan/.
sensitivity:
	cargo build --release -p sweep
	PYTHONPATH=python uv run --extra sci python -m puffsat.sensitivity
