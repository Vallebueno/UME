#!/bin/bash

# Master UME Launcher in VBC CLIP
#
# Author: Miguel Vallebueno,
# Date: 2022-07-29


#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p c
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --qos=medium
#SBATCH --array=5-5

ml build-env/.f2021
ml build-env/f2021
ml pigz/2.6-gcccore-10.2.0

CHR=$SLURM_ARRAY_TASK_ID

#$caronte/UMCAL/src/UME_RCALL_V1.0.sh -f /path/to/workdir/DB/1.ALL.Merge.db.gz -c $CHR -f 1.ALL.Merge.db.gz -x union -l 5


### first run in a single core
##$caronte/UMCAL/src/UME_RCALL_V2.1.sh clist



### then run discovery
#$caronte/UMCAL/src/UME_RCALL_V2.1.sh -x discovery -f /path/to/workdir/MERGE/All.tlone.Merge.db.gz -c $CHR -k /path/to/workdir/MERGE/files.mergein.lst



### then run production
#srcDIR="/path/to/legacy/Caronte/code/UMCAL/src"
# perl ${srcDIR}/Pre/Merge_list2paste3.pl /path/to/mpileup/in_prod_mpileup_files.lst


$caronte/UMCAL/src/UME_RCALL_V2.1.sh -f /path/to/mpileup/in_prod_mpileup_files.lst -x production -i /path/to/mpileup -k /path/to/maizedb/Production/Discovery_06_2023_cof5.ll -r /path/to/reference/Zm-B73-REFERENCE-NAM-5.0.fa.OL -d /path/to/workdir/Production


#$caronte/UMCAL/src/UME_RCALL_V2.1.sh -f /path/to/mpileup/in_prod_mpileup_files.lst -x production -i /path/to/mpileup -k /path/to/workdir/Production/GBS2/MasterLandraceTeoInbredGBS_collapseDist0.02.ll -r /path/to/reference/Zm-B73-REFERENCE-NAM-5.0.fa.OL -d /path/to/workdir/Production/GBS2
#$caronte/UMCAL/src/UME_RCALL_V2.1.sh -f /path/to/mpileup/in_prod_mpileup_files.lst -x production -i /path/to/mpileup -k /path/to/workdir/Production/cof${CHR}/Discovery_LL.ALLC.cof${CHR}.ll -r /path/to/reference/Zm-B73-REFERENCE-NAM-5.0.fa.OL -d /path/to/workdir/Production/cof${CHR}
#$caronte/UMCAL/src/UME_RCALL_V2.1.sh -f /path/to/mpileup/in_prod_mpileup_files.lst -x production -i /path/to/mpileup -k /path/to/workdir/Production/HM3/Discovery_LL.ALLC.HM3.ll -r /path/to/reference/Zm-B73-REFERENCE-NAM-5.0.fa.OL -d /path/to/workdir/Production/HM3
#$caronte/UMCAL/src/UME_RCALL_V2.1.sh -f /path/to/mpileup/in_prod_mpileup_files.lst -x production -i /path/to/mpileup -k /path/to/workdir/Production/GBS/ZeaGBSv27_publicSamples_raw_AGPv5_crossmap.ll -r /path/to/reference/Zm-B73-REFERENCE-NAM-5.0.fa.OL -d /path/to/workdir/Production/GBS
#$caronte/UMCAL/src/UME_RCALL_V2.1.sh -f /path/to/mpileup/in_prod_mpileup_files.lst -x production -i /path/to/mpileup -k /path/to/workdir/Production/GBSK/Maize_GBS_Kriztian_build.ll -r /path/to/reference/Zm-B73-REFERENCE-NAM-5.0.fa.OL -d /path/to/workdir/Production/GBSK
