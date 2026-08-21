# `audit/` — reproduction scripts for the cross-repo calculation audit

Standalone reductions backing the findings in
[`../companion_repo_calculations_2026-08-21.md`](../companion_repo_calculations_2026-08-21.md).
Each is cited by name from the item it supports. Run from the repo root:

```
uv run python audit/plate_transfer_check.py    # does f ~ 0.818 transfer to the paper's 25 kg / 5 m plate?
uv run python audit/wide_footprint_check.py    # design §7's remaining r_foot/R nodes (needs `make sweep-geometry-wide`)
uv run python audit/focusing_check.py          # Q-D: the contour omits the focusing factor the frontier applies
```

**These are exploratory reductions, not library code.** They read committed results and the
`puffsat` package, monkeypatch a cache or two, and print a table. They are deliberately outside
`python/` — the `make lint` gate (`ruff check python`, `mypy files = ["python"]`) covers the
package, and holding one-off audit scripts to mypy-strict would be ceremony without payoff.

When a finding here graduates into a decision, it moves *into* `python/puffsat/` with tests and an
ADR, and the script here becomes redundant. Q-D's focusing fix is the first candidate; it is
pending sign-off because it moves a published deliverable.
