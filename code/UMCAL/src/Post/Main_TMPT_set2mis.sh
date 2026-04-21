#!/bin/bash

# mv files SLURM job script
#
# Author: Miguel Vallebueno,
# Date: 2021-02-22

# === SLURM Config ===

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p c
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --qos=short
#SBATCH --array=1-884

ml build-env/f2021
ml pigz/2.6-gcccore-10.2.0

# === INSTRUCTIONS ===

KEY="/scratch-cbe/users/miguel.vallebueno/Production/cof5/tmpt.lst"

line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${KEY})
my_array=($(echo $line | tr "\t" "\n"))
file=${my_array[0]}
echo "<EasyTaskarr><$SLURM_ARRAY_TASK_ID><$file>"

perl $caronte/UMCAL/src/Post/tmpt_set2mis.pl $file