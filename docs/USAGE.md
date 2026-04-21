# UME Usage Guide

This document focuses on the active workflow:

- `clist`
- `discovery`
- `production`

The recommended entry point is:

- [bin/ume](D:\MV_Postdoc_GMI_KS\UME_MER_REPO\Caronte\bin\ume)

## Workflow summary

### 1. `clist`

Purpose:

- derive caller metadata from the source list
- create helper files for future merge steps
- estimate empirical caller-specific quality distributions

Typical command:

```bash
bin/ume clist \
  -f /path/to/All.tlone.Merge.db.gz \
  -k /path/to/files.mergein.lst \
  -r /path/to/reference.fa.OL
```

Key outputs:

- `<list>.Disc.lst`
- `<list>.tlone`
- `<list>.hh`
- `<merged_db>.pval2f`

### 2. `discovery`

Purpose:

- combine evidence from callers
- threshold known positions
- produce the final discovery `.ll`

Typical command:

```bash
bin/ume discovery \
  -f /path/to/All.tlone.Merge.db.gz \
  -k /path/to/files.mergein.lst \
  -c 1 \
  -l 5 \
  -r /path/to/reference.fa.OL
```

Key scientific output:

- discovery `.ll`

### 3. `production`

Purpose:

- use the discovery `.ll` to project/call all production samples
- merge all sample columns into a final multi-sample output

Typical command:

```bash
bin/ume production \
  -f /path/to/in_prod_mpileup_files.lst \
  -i /path/to/mpileup_dir \
  -k /path/to/discovery.ll \
  -r /path/to/reference.fa.OL \
  -d /path/to/output_dir
```

Key helper artifacts:

- `<keyfile>.tlone`
- `<keyfile>.hh`

These are created automatically if missing.

## Expected file roles

- `.ll`
  Discovery site-definition product used scientifically by `production`.
- `.tlone`
  Shell recipe used to merge many sample columns.
- `.hh`
  VCF-like header used during final assembly.

## SLURM vs local execution

To run locally:

```bash
export UME_USE_SLURM=0
```

To run with SLURM:

```bash
export UME_USE_SLURM=1
```

The active workflow will use `sbatch` when it is available and SLURM is not
disabled.
