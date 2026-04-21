# Test Notes

This directory contains lightweight validation helpers for the publication-ready
repository layout.

The current repository does not yet include a fully automated end-to-end test
suite because the workflow depends on Perl, domain-specific compressed inputs,
and optional HPC infrastructure.

The main current smoke test is:

- confirm that the public wrapper exists
- confirm that documentation/config/example files are present
- confirm that the active engine path is documented
