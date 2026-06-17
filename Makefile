# PuffSat impact simulation — single build entry point (ADR-0018).
# Delegates to cargo (Rust hot path) and uv (Python cold path); the two meet only in data/.

PY := uv run python

.PHONY: all smoke build test lint fmt clean tables sweep analysis

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

## test: run all tests (cargo; pytest once Python tests exist)
test:
	cargo test

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
## tables: generate EOS/opacity tables (rung B+)
tables:
	@echo "TODO (rung B): uv run python -m puffsat.tables  ->  data/tables/*.json"

## sweep: run the parameter sweep (rung B+); depends on tables
sweep:
	@echo "TODO (rung B): cargo run -p sweep --release  ->  data/results/*.jsonl"

## analysis: frontier extraction + plots (rung B+); depends on sweep
analysis:
	@echo "TODO (rung B): uv run python -m puffsat.analysis"
