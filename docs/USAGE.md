# UME Usage Guide

This guide documents the active public workflow around
[`bin/ume`](../bin/ume) and
[`code/UMCAL/src/UME_RCALL_V2.1.sh`](../code/UMCAL/src/UME_RCALL_V2.1.sh).

## Workflow summary

### 1. `clist`

Purpose:

- derive caller metadata from the source list
- create helper files for later merge steps
- estimate empirical caller-specific quality distributions

Command:

```bash
bin/ume clist \
  -f /path/to/All.tlone.Merge.db.gz \
  -k /path/to/files.mergein.lst \
  -r /path/to/reference.fa.OL
```

Outputs:

- `<list>.Disc.lst`
- `<list>.tlone`
- `<list>.hh`
- `<merged_db>.pval2f`

### 2. `discovery`

Purpose:

- combine evidence across callers
- apply the phred cutoff
- produce the final discovery `.ll`

Command:

```bash
bin/ume discovery \
  -f /path/to/All.tlone.Merge.db.gz \
  -k /path/to/files.mergein.lst \
  -c 1 \
  -l 5 \
  -r /path/to/reference.fa.OL
```

Outputs:

- `<merged_db>.qual0.Union.max.dbx`
- `<merged_db>.qual0.Union.max.dbx.cof<cutoff>.lst`
- `<merged_db>.qual0.Union.max.dbx.cof<cutoff>.lst2` when GBS enrichment is enabled
- final discovery `.ll`

### 3. `production`

Purpose:

- project the discovery `.ll` into production samples
- write one temporary genotype column per sample
- merge those columns into the final multi-sample production output

Command:

```bash
bin/ume production \
  -f /path/to/in_prod_mpileup_files.lst \
  -i /path/to/mpileup_dir \
  -k /path/to/discovery.ll \
  -r /path/to/reference.fa.OL \
  -d /path/to/output_dir
```

Outputs:

- one `*.tmpt.gz` file per production sample
- merged production VCF-like output in the output directory

## File roles

- `.ll`
  Discovery site-definition product used by `production`.
- `.tlone`
  Merge recipe used to combine all per-sample columns.
- `.hh`
  Header file used during final output assembly.

## Environment variables

- `UME_SRC_DIR`
  Overrides the engine source directory.
- `UME_USE_SLURM`
  Set to `0` to force local sequential execution.
- `UME_GBS_SUMMARY`
  Optional path to the summary table used for discovery enrichment. If unset,
  discovery continues without the optional enrichment step.
- `UME_REPO_ROOT`
  Optional repository root containing `Modules/` for extra postprocessing.
- `UME_SUM_PVAL_TABLE`
  Required only if the discovery aggregation mode is changed to `sum`.
- `UME_PRODUCTION_PREFIX`
  Prefix for merged production output names.

## Local vs SLURM execution

To force local execution:

```bash
export UME_USE_SLURM=0
```

To allow SLURM execution when `sbatch` is available:

```bash
export UME_USE_SLURM=1
```

The workflow automatically falls back to local execution if `sbatch` is not
available.
