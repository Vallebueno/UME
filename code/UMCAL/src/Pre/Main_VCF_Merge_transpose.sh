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

REF="/groups/swarts/lab/Resources/ReferenceGenomes/Zm-B73-REFERENCE-NAM-5.0.fa";
KEY="/scratch-cbe/users/miguel.vallebueno/TEST/Merge/list.txt";
BED="/groups/swarts/lab/Resources/ReferenceGenomes/Zm-B73-REFERENCE-NAM-5.0.bed";

line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${KEY})
my_array=($(echo $line | tr "\t" "\n"))

echo "$my_array"




perl /groups/swarts/lab/Resources/Scripts/Caronte/code/UMCAL/src/VCF_Merge_HidraV20.pl $REF $my_array $BED