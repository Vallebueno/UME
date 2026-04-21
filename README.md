# UME: Unified Multi-caller Ensemble

UME (Unified Multi-caller Ensemble) is a variant-calling integration framework
for heterogeneous sequencing datasets. The public workflow preserved in this
repository centers on `UME_RCALL_V2.1.sh`, which builds caller metadata,
generates a discovery `.ll` site list, and projects that site list into
production samples.

The repository has been cleaned for public release with these goals:

- keep the core scientific workflow intact
- remove local cluster paths and infrastructure details
- separate active workflow code from historical material
- make the runnable path clearer for Linux and HPC users

## Active workflow

The recommended entry point is [`bin/ume`](bin/ume), which wraps
[`code/UMCAL/src/UME_RCALL_V2.1.sh`](code/UMCAL/src/UME_RCALL_V2.1.sh).

The active stages are:

1. `clist`
   Creates caller metadata and empirical quality lookup tables from an existing
   merged DB and the source input list.
2. `discovery`
   Builds the ensemble union and writes the discovery `.ll` file.
3. `production`
   Projects the discovery `.ll` into production samples and assembles the final
   merged production output.

Scientific and helper products:

- `.ll` is the discovery site-definition product.
- `.tlone` and `.hh` are helper files used to assemble the final production
  merge and header.

## Repository layout

- [`bin/`](bin)
  Public entry points.
- [`code/UMCAL/`](code/UMCAL)
  Active UME engine and workflow-specific documentation.
- [`config/`](config)
  Example environment settings.
- [`docs/`](docs)
  Usage and repository notes.
- [`examples/`](examples)
  Minimal example inputs and file-layout guidance.
- [`tests/`](tests)
  Lightweight smoke checks.
- [`legacy/`](legacy)
  Historical, superseded, or non-core material retained outside the main
  execution surface.

## Requirements

Core runtime:

- Linux
- `bash`
- `perl`
- `gzip` or `pigz`

Optional HPC/runtime helpers:

- `sbatch` for SLURM array execution
- `ml` environment modules when available

Optional postprocessing helpers:

- `Modules/VCF/VCF_RANDOM_SAMPLE.pl`
- `Modules/Tassel/TasselTRG_D_Matrix_downsample_HETS.sh`

The core workflow runs without those optional helper modules.

## Quick start

```bash
git clone https://github.com/Vallebueno/UME.git
cd UME
cp config/ume.env.example my-ume.env
# edit my-ume.env for your environment
source my-ume.env
bin/ume --help
```

Minimal commands:

```bash
bin/ume clist \
  -f /path/to/All.tlone.Merge.db.gz \
  -k /path/to/files.mergein.lst \
  -r /path/to/reference.fa.OL

bin/ume discovery \
  -f /path/to/All.tlone.Merge.db.gz \
  -k /path/to/files.mergein.lst \
  -c 1 \
  -l 5 \
  -r /path/to/reference.fa.OL

bin/ume production \
  -f /path/to/in_prod_mpileup_files.lst \
  -i /path/to/mpileup_dir \
  -k /path/to/discovery.ll \
  -r /path/to/reference.fa.OL \
  -d /path/to/output_dir
```

## Documentation

- Workflow usage: [`docs/USAGE.md`](docs/USAGE.md)
- Engine details: [`code/UMCAL/README.md`](code/UMCAL/README.md)
- Cleanup notes: [`docs/REPOSITORY_NOTES.md`](docs/REPOSITORY_NOTES.md)
- Toy example: [`examples/toy-workflow`](examples/toy-workflow)

## Release metadata

- Citation metadata: [`CITATION.cff`](CITATION.cff)
- Changelog: [`CHANGELOG.md`](CHANGELOG.md)
- License: [`LICENSE`](LICENSE)
- Current version: [`VERSION`](VERSION)
