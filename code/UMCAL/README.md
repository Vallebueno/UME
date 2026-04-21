# UMCAL Core Workflow

This directory contains the active UME/UMCAL pipeline code used to build a
discovery site list and project those sites into production calls.

The current workflow is centered on:

- [`src/UME_RCALL_V2.1.sh`](D:\MV_Postdoc_GMI_KS\UME_MER_REPO\Caronte\code\UMCAL\src\UME_RCALL_V2.1.sh)

This README documents only the core runnable workflow around that script. It
does not try to describe every historical helper or manuscript-support file in
the repository.

## Workflow Overview

The active workflow has three main stages:

1. `clist`
   Builds caller metadata and empirical quality lookup tables from an existing
   merged DB plus the list of source files/callers.
2. `discovery`
   Computes the ensemble union, applies a phred cutoff, optionally enriches the
   site list with GBS diversity, and writes a discovery `.ll` file.
3. `production`
   Uses the discovery `.ll` file to create one temporary genotype column per
   production input and then merges those columns into a VCF-like output.

There is also an upstream pre-merge stage in `src/Pre/` that prepares the
merged DB consumed by `clist` and `discovery`, but that stage is not yet fully
portable and is documented here only as context.

## Expected Stage Inputs

### `clist`

Required arguments:

- `-x clist`
- `-f <merged_db.gz>`
- `-k <input_list.txt>`

Optional but recommended:

- `-r <reference.fa.OL>`

What the script expects:

- `-f` points to the merged DB used to estimate empirical caller quality
  distributions.
- `-k` points to the original list of per-sample inputs. This list is used to
  derive:
  - `input_list.txt.Disc.lst`
  - `input_list.txt.tlone`
  - `input_list.txt.hh`

Key outputs:

- `<input_list>.Disc.lst`
- `<input_list>.tlone`
- `<input_list>.hh`
- `<merged_db>.pval2f`

Artifact roles:

- `.Disc.lst` is caller metadata used during discovery.
- `.tlone` is an auxiliary merge instruction file.
- `.hh` is an auxiliary header file used when assembling the merged output.
- `.tlone` stores the shell `paste` recipe used to combine all sample columns.

### `discovery`

Required arguments:

- `-x discovery`
- `-f <merged_db.gz>`
- `-k <input_list.txt>`
- `-c <chromosome>`
- `-l <phred_cutoff>`
- `-r <reference.fa.OL>`

What the script expects:

- `<input_list>.Disc.lst` already exists, usually from the `clist` stage.
- `<merged_db>.pval2f` already exists, usually from the `clist` stage.

Key outputs:

- `<merged_db>.qual0.Union.max.dbx`
- `<merged_db>.qual0.Union.max.dbx.cof<cutoff>.lst`
- `<merged_db>.qual0.Union.max.dbx.cof<cutoff>.lst2` when GBS enrichment is
  available
- final discovery `.ll` file derived from the enriched list when possible, or
  from the cutoff list otherwise

Artifact role:

- `.ll` is the main scientific product of the discovery stage used by
  `production`.

### `production`

Required arguments:

- `-x production`
- `-f <production_keyfile.txt>`
- `-i <production_input_dir>`
- `-k <discovery.ll>`
- `-r <reference.fa.OL>`
- `-d <output_dir>`

What the script expects:

- `-f` is a key file containing one production input per line. The first column
  must match the filename found under `-i`.
- `-k` is the discovery `.ll` site-definition file.
- `<keyfile>.tlone` and `<keyfile>.hh` are helper files used to merge all
  production samples and build the final header. The workflow now creates them
  automatically when missing.

Key outputs:

- one `*.tmpt.gz` file per production sample in `-d`
- merged production VCF-like output from the `Merge_columns2DB_production.sh`
  step

## File Contracts

Important intermediate files are connected by naming convention:

- `<list>.Disc.lst`
- `<list>.tlone`
- `<list>.hh`
- `<merged_db>.pval2f`
- `<merged_db>.qual0.Union.max.dbx`
- `<union_db>.cof<cutoff>.lst`
- `<cutoff_list>.lst2`
- `<cutoff_or_enriched_list>.ll`
- `<sample>.tmpt.gz`

The current pipeline still relies on these suffixes. Renaming outputs manually
will usually break downstream steps.

Semantic distinction:

- `.ll` is the discovery output that defines which sites will be projected into
  production.
- `.tlone` and `.hh` are merge/header helper artifacts used to assemble the
  final multi-sample production output.

## Environment Variables

The portability layer currently supports these environment variables:

- `UME_SRC_DIR`
  Overrides the source directory containing `Pre/`, `Discovery/`, and
  `Production/`.
- `UME_USE_SLURM`
  Set to `0` to force local sequential execution instead of `sbatch`.
- `UME_GBS_SUMMARY`
  Path to the GBS summary file used by
  `src/Discovery/Discovery_enrichList_GBS.pl`.
- `UME_CARONTE_ROOT`
  Root directory used to look for optional external `Modules/` helpers during
  production merge/postprocessing.
- `UME_SUM_PVAL_TABLE`
  Optional lookup table for the `sum` aggregation mode in
  `src/Discovery/DB_Union_V10_FAST.pl`. The current default workflow uses
  `max`, but this environment variable preserves portability if `sum` is used.

## Minimal Examples

### `clist`

```bash
bash src/UME_RCALL_V2.1.sh \
  -x clist \
  -f /path/to/All.tlone.Merge.db.gz \
  -k /path/to/files.mergein.lst \
  -r /path/to/reference.fa.OL
```

### `discovery`

```bash
bash src/UME_RCALL_V2.1.sh \
  -x discovery \
  -f /path/to/All.tlone.Merge.db.gz \
  -k /path/to/files.mergein.lst \
  -c 1 \
  -l 5 \
  -r /path/to/reference.fa.OL
```

### `production` with local fallback

```bash
export UME_USE_SLURM=0

bash src/UME_RCALL_V2.1.sh \
  -x production \
  -f /path/to/in_prod_mpileup_files.lst \
  -i /path/to/mpileup_dir \
  -k /path/to/Discovery_output.ll \
  -r /path/to/reference.fa.OL \
  -d /path/to/output_dir
```

## Required Software

Core runtime:

- `bash`
- `perl`
- `gzip` or `pigz`

Cluster/HPC execution:

- `sbatch` for SLURM array execution
- environment modules such as `ml` are optional after the portability cleanup,
  but still used when available

Optional external helpers:

- `Modules/VCF/VCF_RANDOM_SAMPLE.pl`
- `Modules/Tassel/TasselTRG_D_Matrix_downsample_HETS.sh`

The core production merge now skips these optional helpers when they are not
available, but users should expect fewer downstream convenience artifacts in
that case.

## Current Limits

- The upstream pre-merge stage is still only partially portable.
- The repository still contains many historical script variants not yet
  curated.
- The workflow depends on custom input conventions such as `.fa.OL` reference
  files and specific production keyfile formatting.
- This README documents the active central path only.
