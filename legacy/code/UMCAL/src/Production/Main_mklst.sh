#!/bin/bash

# mv files SLURM job script
#
# Author: Miguel Vallebueno,
# Date: 2021-02-22

############  Designed to run with C1 Nodes as they have High clock speeds


# === SLURM Config === for C1 Nodes(hi-clock speed),  qos short(hi-priority)

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p c                 #Node type(c1,c2,m2)
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --qos=medium        #-time +priority
#SBATCH --array=1-10

# === INSTRUCTIONS ===
KEY="/path/to/maizedb/polis2.lst"


line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${KEY})
my_array=($(echo $line | tr "\t" "\n"))

echo "<EasyTaskarr><$my_array>"

perl $caronte/UMCAL/src/Production/Production_call_mklst.pl $my_array
