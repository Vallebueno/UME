# Changelog

## v0.1.0 - 2026-04-21

First publication-ready public release candidate.

### Added

- Root-level repository documentation
- user-facing `bin/ume` wrapper
- `config/`, `docs/`, `examples/`, and `tests/` scaffolding
- citation metadata
- BSD-3-Clause license
- example environment configuration
- toy workflow example

### Changed

- stabilized active workflow around `code/UMCAL/src/UME_RCALL_V2.1.sh`
- reduced hard-coded path assumptions in the core path
- documented `.ll` as the discovery product and `.tlone`/`.hh` as merge helpers
- added local fallback execution for production when SLURM is unavailable
- improved production merge consistency and helper generation

### Notes

- Historical scripts remain in the repository for backward compatibility and
  traceability.
- Upstream pre-merge helpers are still more environment-specific than the core
  active pipeline.
