# Repository Notes

## Public workflow surface

The public repository now keeps the active execution path deliberately small.
The runnable engine is:

- `code/UMCAL/src/UME_RCALL_V2.1.sh`

Its retained direct dependencies are:

- `code/UMCAL/src/Pre/MK_List2IN.pl`
- `code/UMCAL/src/Pre/Merge_list2paste3.pl`
- `code/UMCAL/src/Pre/Final_header_make.pl`
- `code/UMCAL/src/Discovery/EmpDist2PvalV3.pl`
- `code/UMCAL/src/Discovery/DB_Union_V10_FAST.pl`
- `code/UMCAL/src/Discovery/Discovery_mklst_DB.pl`
- `code/UMCAL/src/Discovery/Discovery_enrichList_GBS.pl`
- `code/UMCAL/src/Discovery/Discovery_sort_alts.pl`
- `code/UMCAL/src/Production/Production_Parallel_mpileup_lst_CLIP.sh`
- `code/UMCAL/src/Production/Production_Call_mpileup_list_OC_V9.pl`
- `code/UMCAL/src/Production/CHECK_tmptOLS.sh`
- `code/UMCAL/src/Production/Merge_columns2DB_production.sh`

## Legacy material

Historical, superseded, exploratory, figure-generation, and non-core helper
files were moved under:

- `legacy/`

This keeps them available for traceability while reducing the surface area of
the active public workflow.

## Path sanitization

The active workflow and public-facing documentation were cleaned to avoid
revealing local workstation or cluster directory layouts. Legacy text files
retained in `legacy/` were also sanitized to remove cluster-specific absolute
paths and usernames.
