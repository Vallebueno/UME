#!/usr/bin/env bash
# Transpose a gvcf into merged hapmaplike
# Author: Miguel Vallebueno
# Date: 2021-07-23

# === SLURM Config ===

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p c                 #Node type(c1,c2,m2)
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --qos=long      #-time +priority
#SBATCH --array=1-1

# === Load Modules ===

# === INPUT Config ===

REF="/path/to/reference/Zm-B73-REFERENCE-NAM-5.0.fa";
KEY="/path/to/workdir/TEST/Merge/list.txt";
BED="/path/to/reference/Zm-B73-REFERENCE-NAM-5.0.bed";

line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${KEY})
my_array=($(echo $line | tr "\t" "\n"))

echo "$my_array"




perl /path/to/legacy/Caronte/code/UMCAL/src/VCF_Merge_HidraV20.pl $REF $my_array $BED
