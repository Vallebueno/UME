#!/bin/bash

# Master UME Launcher in VBC CLIP
#
# Author: Miguel Vallebueno,
# Date: 2022-07-29

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p c
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --qos=medium
#SBATCH --array=1-10

ml build-env/.f2021
ml build-env/f2021
ml pigz/2.6-gcccore-10.2.0

KEY="/path/to/workdir/MERGE/SPLIT/segmentac"
line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${KEY})
my_array=($(echo $line | tr "\t" "\n"))

echo "<EasyTaskarr><$SLURM_ARRAY_TASK_ID>"


CHR=$SLURM_ARRAY_TASK_ID

#$caronte/UMCAL/src/UME_RCALL_V1.0.sh -f /path/to/workdir/MERGE/SPLIT/$my_array -c $CHR -f 1.ALL.Merge.db.gz -x union -l 5


### first run in a single core
##$caronte/UMCAL/src/UME_RCALL_V2.1.sh clist



### then run discovery

#zcat All.tlone.Merge.db.${CHR}.poli.db.part*.gz.qual0.Union.max.dbx.gz > cat_All.tlone.Merge.db.${CHR}.poli.db.qual0.Union.max.dbx.gz


$caronte/UMCAL/src/UME_RCALL_V2.1_LST.sh -x discovery -f /path/to/workdir/MERGE/SPLIT/cat_All.tlone.Merge.db.${CHR}.poli.db.qual0.Union.max.dbx.gz -c $CHR -k /path/to/workdir/MERGE/SPLIT/files.mergein.lst -l 30



### then run production
#$caronte/UMCAL/src/UME_RCALL_V2.1.sh -f /path/to/workdir/MERGE/All.tlone.Merge.db.gz -c $CHR -x production -l 5
