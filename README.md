# UME: Unified Multi-caller Ensemble

UME (Unified Multi-caller Ensemble) is a variant-calling integration framework
for heterogeneous sequencing datasets. It combines information from multiple
variant callers and then projects the resulting discovery site set into
production samples, with the goal of reducing coverage-related bias and making
comparisons across uneven datasets more robust.

This repository is being prepared as a publication-ready software package for
long-term use on Linux servers and HPC systems. The active pipeline is still
implemented primarily in Bash and Perl, and this repository preserves that
implementation while adding portability, clearer structure, and user-facing
documentation.

## Repository layout

- [bin](D:\MV_Postdoc_GMI_KS\UME_MER_REPO\Caronte\bin)
  Public entry points and convenience wrappers.
- [code/UMCAL](D:\MV_Postdoc_GMI_KS\UME_MER_REPO\Caronte\code\UMCAL)
  Active UME engine and legacy source tree.
- [config](D:\MV_Postdoc_GMI_KS\UME_MER_REPO\Caronte\config)
  Example environment/configuration files.
- [docs](D:\MV_Postdoc_GMI_KS\UME_MER_REPO\Caronte\docs)
  Detailed usage, workflow, and publication notes.
- [examples](D:\MV_Postdoc_GMI_KS\UME_MER_REPO\Caronte\examples)
  Small toy inputs showing expected file structure and command usage.
- [tests](D:\MV_Postdoc_GMI_KS\UME_MER_REPO\Caronte\tests)
  Lightweight smoke tests and validation helpers.

## What UME does

The active workflow currently centers on three stages:

1. `clist`
   Builds caller metadata and empirical quality lookup tables from a merged DB
   and the list of source files/callers.
2. `discovery`
   Computes the ensemble union, applies a phred threshold, optionally enriches
   the site list with GBS diversity, and writes the discovery `.ll` file.
3. `production`
   Uses the discovery `.ll` file to call/project those sites into production
   samples and merge all per-sample temporary outputs into a multi-sample
   VCF-like output.

The main scientific discovery product is the `.ll` site list.

The `.tlone` and `.hh` files are merge/header helper artifacts used during
assembly of the final production output.

## Main entry point

The recommended user-facing entry point is:

- [bin/ume](D:\MV_Postdoc_GMI_KS\UME_MER_REPO\Caronte\bin\ume)

This wrapper delegates to:

- [code/UMCAL/src/UME_RCALL_V2.1.sh](D:\MV_Postdoc_GMI_KS\UME_MER_REPO\Caronte\code\UMCAL\src\UME_RCALL_V2.1.sh)

## Software requirements

Core requirements:

- Linux
- `bash`
- `perl`
- `gzip` or `pigz`

Optional HPC requirements:

- `sbatch` for SLURM execution
- `ml`/environment modules when available


The core active workflow can now run without those optional helpers, but some
downstream convenience outputs will be skipped.

## Installation

UME currently uses a lightweight installation model:

1. Clone the repository.
2. Ensure `bash`, `perl`, and `gzip` are available.
3. Copy and edit:
   - [config/ume.env.example](D:\MV_Postdoc_GMI_KS\UME_MER_REPO\Caronte\config\ume.env.example)
4. Source your environment file before running UME.

Example:

```bash
git clone <repo-url> ume
cd ume
cp config/ume.env.example my-ume.env
# edit my-ume.env
source my-ume.env
bin/ume --help
```

## Minimal usage

See:

- [docs/USAGE.md](D:\MV_Postdoc_GMI_KS\UME_MER_REPO\Caronte\docs\USAGE.md)
- [examples/toy-workflow](D:\MV_Postdoc_GMI_KS\UME_MER_REPO\Caronte\examples\toy-workflow)

## Inputs and outputs

Expected inputs and outputs for the active pipeline are documented in:

- [docs/USAGE.md](D:\MV_Postdoc_GMI_KS\UME_MER_REPO\Caronte\docs\USAGE.md)
- [code/UMCAL/README.md](D:\MV_Postdoc_GMI_KS\UME_MER_REPO\Caronte\code\UMCAL\README.md)

## Computational notes

- `production` is designed to parallelize across samples and works best on HPC
  systems with SLURM, but a local sequential fallback is available.
- The upstream pre-merge stage remains more environment-specific than the core
  `clist -> discovery -> production` path.
- Memory and runtime depend strongly on the number of samples, number of sites,
  and compression/decompression throughput.

## Citation

Citation metadata for software release is provided in:

- [CITATION.cff](D:\MV_Postdoc_GMI_KS\UME_MER_REPO\Caronte\CITATION.cff)

If you use UME in academic work, please cite:

1. The manuscript describing UME.
2. The archived software release associated with the repository DOI.

## Versioning and release

This repository is being prepared for an initial stable public release.

Recommended first release tag:

- `v0.1.0`

Rationale:

- first cleaned and publication-ready public release
- suitable for Zenodo archiving
- semantically indicates a stable but early public version

See also:

- [CHANGELOG.md](D:\MV_Postdoc_GMI_KS\UME_MER_REPO\Caronte\CHANGELOG.md)

## License

This repository includes:

- [LICENSE](D:\MV_Postdoc_GMI_KS\UME_MER_REPO\Caronte\LICENSE)

## Current status

The active workflow has been stabilized for portability and documentation, but
some historical scripts and legacy helper variants are still preserved in the
source tree for backward compatibility and manuscript traceability.
