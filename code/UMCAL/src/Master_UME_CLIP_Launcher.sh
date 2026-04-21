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

#$caronte/UMCAL/src/UME_RCALL_V1.0.sh -f /scratch-cbe/users/miguel.vallebueno/DB/1.ALL.Merge.db.gz -c $CHR -f 1.ALL.Merge.db.gz -x union -l 5


### first run in a single core
##$caronte/UMCAL/src/UME_RCALL_V2.1.sh clist



### then run discovery
#$caronte/UMCAL/src/UME_RCALL_V2.1.sh -x discovery -f /scratch-cbe/users/miguel.vallebueno/MERGE/All.tlone.Merge.db.gz -c $CHR -k /scratch-cbe/users/miguel.vallebueno/MERGE/files.mergein.lst



### then run production
#srcDIR="/groups/swarts/lab/Resources/Scripts/Caronte/code/UMCAL/src"
# perl ${srcDIR}/Pre/Merge_list2paste3.pl /groups/swarts/lab/Resources/Maize/WGS/mpileup/in_prod_mpileup_files.lst


$caronte/UMCAL/src/UME_RCALL_V2.1.sh -f /groups/swarts/lab/Resources/Maize/WGS/mpileup/in_prod_mpileup_files.lst -x production -i /groups/swarts/lab/Resources/Maize/WGS/mpileup -k /groups/swarts/lab/MAVE/MAIZEDB/Production/Discovery_06_2023_cof5.ll -r /groups/swarts/lab/Resources/ReferenceGenomes/Zm-B73-REFERENCE-NAM-5.0.fa.OL -d /scratch-cbe/users/miguel.vallebueno/Production


#$caronte/UMCAL/src/UME_RCALL_V2.1.sh -f /groups/swarts/lab/Resources/Maize/WGS/mpileup/in_prod_mpileup_files.lst -x production -i /groups/swarts/lab/Resources/Maize/WGS/mpileup -k /scratch-cbe/users/miguel.vallebueno/Production/GBS2/MasterLandraceTeoInbredGBS_collapseDist0.02.ll -r /groups/swarts/lab/Resources/ReferenceGenomes/Zm-B73-REFERENCE-NAM-5.0.fa.OL -d /scratch-cbe/users/miguel.vallebueno/Production/GBS2
#$caronte/UMCAL/src/UME_RCALL_V2.1.sh -f /groups/swarts/lab/Resources/Maize/WGS/mpileup/in_prod_mpileup_files.lst -x production -i /groups/swarts/lab/Resources/Maize/WGS/mpileup -k /scratch-cbe/users/miguel.vallebueno/Production/cof${CHR}/Discovery_LL.ALLC.cof${CHR}.ll -r /groups/swarts/lab/Resources/ReferenceGenomes/Zm-B73-REFERENCE-NAM-5.0.fa.OL -d /scratch-cbe/users/miguel.vallebueno/Production/cof${CHR}
#$caronte/UMCAL/src/UME_RCALL_V2.1.sh -f /groups/swarts/lab/Resources/Maize/WGS/mpileup/in_prod_mpileup_files.lst -x production -i /groups/swarts/lab/Resources/Maize/WGS/mpileup -k /scratch-cbe/users/miguel.vallebueno/Production/HM3/Discovery_LL.ALLC.HM3.ll -r /groups/swarts/lab/Resources/ReferenceGenomes/Zm-B73-REFERENCE-NAM-5.0.fa.OL -d /scratch-cbe/users/miguel.vallebueno/Production/HM3
#$caronte/UMCAL/src/UME_RCALL_V2.1.sh -f /groups/swarts/lab/Resources/Maize/WGS/mpileup/in_prod_mpileup_files.lst -x production -i /groups/swarts/lab/Resources/Maize/WGS/mpileup -k /scratch-cbe/users/miguel.vallebueno/Production/GBS/ZeaGBSv27_publicSamples_raw_AGPv5_crossmap.ll -r /groups/swarts/lab/Resources/ReferenceGenomes/Zm-B73-REFERENCE-NAM-5.0.fa.OL -d /scratch-cbe/users/miguel.vallebueno/Production/GBS
#$caronte/UMCAL/src/UME_RCALL_V2.1.sh -f /groups/swarts/lab/Resources/Maize/WGS/mpileup/in_prod_mpileup_files.lst -x production -i /groups/swarts/lab/Resources/Maize/WGS/mpileup -k /scratch-cbe/users/miguel.vallebueno/Production/GBSK/Maize_GBS_Kriztian_build.ll -r /groups/swarts/lab/Resources/ReferenceGenomes/Zm-B73-REFERENCE-NAM-5.0.fa.OL -d /scratch-cbe/users/miguel.vallebueno/Production/GBSK