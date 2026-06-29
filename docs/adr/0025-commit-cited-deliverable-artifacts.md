# The cited deliverable artifacts are committed in place, despite the `data/` boundary being gitignored

The `data/` directory is the generated Python↔Rust boundary (ADR-0002/0007/0019) and is
gitignored: tables JSON and sweep JSONL are produced by `make tables` / `make sweep` and were
declared "never checked in." But `CONCLUSION.md` — the study's headline deliverable, linked from
the README — cites specific files in that directory as its evidence: `data/results/frontier_margin.csv`
(the margin sweep), `data/results/frontier.csv` (the loss budget), and generically
`data/results/frontier*.csv` and "the `*.png` plots." On a fresh clone those citations dangle. This
ADR records why we resolved that by committing the cited artifacts *in place* rather than de-citing
them or copying them elsewhere — because committing files inside a directory the repo's own rules
call "never checked in" is surprising without context.

**The cited artifacts are the study's output, not a regenerable intermediate.** The data in `data/`
is not homogeneous. The input tables (`tables/*.json`) and raw sweep dumps (`results/sweep*.jsonl`,
~1.5 MB) are bulk regenerable intermediates. The frontier CSVs are ~27 KB of text — the actual
quoted numbers, diffable, and the *result* of the analysis stage, not its input. A headline document
whose credibility rests on inspectable numbers should ship those numbers. So the carve-out is along
the input/intermediate-vs-deliverable seam, not the whole directory.

**"Just run `make`" is not a guaranteed reproduction.** The study is wrapped, and the real-opacity
table was firewall-blocked (default-deny network policy; CONCLUSION.md §Status). A fresh clone cannot
re-derive these artifacts on demand — which removes the usual argument for leaving generated outputs
uncommitted and makes the committed copy the only available evidence.

**Committed in place, not copied.** The files are un-ignored where `make analysis` writes them, via a
negation block in `.gitignore`:

```gitignore
/data/*
!/data/results/
/data/results/*
!/data/results/frontier*.csv
!/data/results/*.png
```

This keeps a **single source of truth** — the file the pipeline produces *is* the file git tracks, so
there is no second copy to drift — and leaves `CONCLUSION.md`'s existing citations resolving without
edits. The cost is that `data/` is no longer "purely generated" and the gitignore carries a four-line
negation; the explanatory comment there points back to this ADR.

**Scope: headline CSVs + figures only.** The seven `frontier*.csv` and the eight top-level
`results/*.png` are committed (~630 KB total, the figures being the bulk). The `results/opacity_scan/`
subdirectory (a supporting non-blocking refinement, not a headline-cited figure) and all bulk
tables/JSONL stay ignored.

## Consequence

`CONCLUSION.md`'s evidence is inspectable on a fresh clone with no broken links, while the regenerable
boundary stays out of git. The repo's "data is gitignored" statements (`.gitignore`, CLAUDE.md
architecture note) are amended to name this exception so a future contributor does not "helpfully"
re-ignore or delete the committed artifacts. The ~600 KB of binary PNGs entering git history is
accepted as a one-time cost for a wrapped study; it is not a license to commit further generated
binaries — new bulk outputs remain ignored by default.

## Considered Options

- **Commit the cited deliverable artifacts in place (chosen).** Single source of truth, citations
  unchanged, evidence viewable on clone. Cost: `data/` no longer purely generated; gitignore negation
  block; ~600 KB of binaries in history.
- **De-cite entirely — make `CONCLUSION.md` stand alone with inline numbers.** Rejected: honors the
  "never checked in" rule literally but discards inspectable provenance for the headline claim, and the
  firewall-blocked table means a reader cannot regenerate it themselves.
- **Copy the artifacts to a checked-in location (e.g. `docs/results/`).** Rejected: keeps `data/`
  conceptually pure but introduces a duplicate that drifts from the canonical `make`-written file every
  re-run, and forces rewriting every citation.
- **Commit all of `data/` (drop the ignore).** Rejected: 2.1 MB including regenerable intermediates,
  against the boundary's intent (ADR-0019); solves the symptom by abandoning the rule rather than
  scoping an exception to it.
