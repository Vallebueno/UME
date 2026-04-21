# Repository Notes

## Active engine

The active engine remains in:

- `code/UMCAL/src`

This tree is preserved for backward compatibility and to minimize risk of
changing scientific behavior during publication hardening.

## Public-facing layer

The publication-oriented repository layout adds:

- `bin/`
- `config/`
- `docs/`
- `examples/`
- `tests/`

This layer improves usability without forcing a rewrite of the scientific core.

## Legacy code

Historical variants such as older `UME_RCALL_*`, `DB_Union_*`, and production
helper versions are still present. They have not yet been removed in order to
avoid accidental loss of scientific traceability or backward compatibility.
