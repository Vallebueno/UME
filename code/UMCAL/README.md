# UMCAL Core Workflow

This directory contains the active UME engine retained for the public release.
The workflow is centered on
[`src/UME_RCALL_V2.1.sh`](src/UME_RCALL_V2.1.sh).

## Core execution path

`UME_RCALL_V2.1.sh` drives three public stages:

1. `clist`
   Uses:
   `src/Pre/MK_List2IN.pl`,
   `src/Pre/Merge_list2paste3.pl`,
   `src/Pre/Final_header_make.pl`,
   `src/Discovery/EmpDist2PvalV3.pl`
2. `discovery`
   Uses:
   `src/Discovery/DB_Union_V10_FAST.pl`,
   `src/Discovery/Discovery_mklst_DB.pl`,
   `src/Discovery/Discovery_enrichList_GBS.pl`,
   `src/Discovery/Discovery_sort_alts.pl`
3. `production`
   Uses:
   `src/Production/Production_Parallel_mpileup_lst_CLIP.sh`,
   `src/Production/Production_Call_mpileup_list_OC_V9.pl`,
   `src/Production/CHECK_tmptOLS.sh`,
   `src/Production/Merge_columns2DB_production.sh`

## Workflow contracts

Important files exchanged between stages:

- `<list>.Disc.lst`
- `<list>.tlone`
- `<list>.hh`
- `<merged_db>.pval2f`
- `<merged_db>.qual0.Union.max.dbx`
- `<union_db>.cof<cutoff>.lst`
- `<cutoff_or_enriched_list>.ll`
- `<sample>.tmpt.gz`

Semantic roles:

- `.ll` is the scientific discovery product.
- `.tlone` and `.hh` are merge/header helper files for production assembly.

## Environment variables

- `UME_SRC_DIR`
  Override the location of the `src/` tree.
- `UME_USE_SLURM`
  Set to `0` to force local execution.
- `UME_GBS_SUMMARY`
  Optional summary table for discovery enrichment.
- `UME_REPO_ROOT`
  Optional root containing `Modules/` for extra postprocessing helpers.
- `UME_CARONTE_ROOT`
  Legacy alias still accepted for backward compatibility.
- `UME_SUM_PVAL_TABLE`
  Optional lookup table used only with `sum` aggregation mode.

## Notes

- The public repository surface is intentionally focused on the active
  `clist -> discovery -> production` path.
- Historical or superseded material has been moved under `legacy/`.
