# PuffSat impact simulation — single build entry point (ADR-0018).
# Delegates to cargo (Rust hot path) and uv (Python cold path); the two meet only in data/.

PY := uv run python

.PHONY: all smoke build test lint fmt clean tables sweep analysis sensitivity

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

## sweep: run the 16 km/s e_eff(rho) sweep (rung B) -> data/results/sweep.jsonl; depends on tables
sweep: data/results/sweep.jsonl

data/results/sweep.jsonl: data/tables/water.json $(wildcard crates/sweep/src/*.rs) $(wildcard crates/hydro1d/src/*.rs)
	@mkdir -p data/results
	cargo run --release -p sweep

## analysis: frontier extraction + plots (rung B) -> data/results/frontier.csv + figures; depends on sweep
analysis: data/results/frontier.csv

data/results/frontier.csv: data/results/sweep.jsonl python/puffsat/analysis.py
	PYTHONPATH=python uv run --extra sci python -m puffsat.analysis

## sensitivity: opacity-insensitivity scan (rung B, B5d-3) — sweep at 0.1x/1x/10x opacity, show
## e_eff barely moves. Builds the release sweep first; writes data/results/opacity_scan/.
sensitivity:
	cargo build --release -p sweep
	PYTHONPATH=python uv run --extra sci python -m puffsat.sensitivity
