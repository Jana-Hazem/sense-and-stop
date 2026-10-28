# Contributing

1. Never commit to `main` directly. Branch as `member<N>/<short-topic>`, for example `member3/solver-baseline`.
2. Open a pull request. At least one other member reviews before merging.
3. Each PR states: what changed, which interface or data contract it touches, and any units or grid assumptions.
4. Never commit raw data, large arrays, credentials, or restricted files.
5. Every dataset used must have an entry in `data/MANIFEST.md` first.
6. Interface changes (grids, units, H/y formats, output metadata) need an issue and a note in `docs/decision_log.md` before merging.
7. Don't label EEG as CSD, and make no clinical claims in code, docs, or commit messages.
8. Keep notebooks reproducible: run top to bottom, with no hidden state.
