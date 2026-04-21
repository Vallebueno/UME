#!/bin/bash

# mv files SLURM job script
#
# Author: Miguel Vallebueno,
# Date: 2021-02-22

############  Designed to run with C1 Nodes as they have High clock speeds


# === SLURM Config === for C1 Nodes(hi-clock speed),  qos short(hi-priority)

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p c                #Node type(c1,c2,m2)
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --qos=medium



# === INSTRUCTIONS ===
ml build-env/f2021
ml pigz/2.6-gcccore-10.2.0

##KEY="/scratch-cbe/users/miguel.vallebueno/TEST/file.parts.lst"
#KEY="/groups/swarts/lab/MAVE/Experiments/Hydra_coverage/OL/list.txt"

#line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${KEY})
#my_array=($(echo $line | tr "\t" "\n"))

echo "<Main_union_parallel><$my_array>"

#perl $caronte/Modules/VCF/DB_Union_V8_FAST.pl max $my_array
perl $caronte/Modules/VCF/DB_Union_V7_FAST.pl
#sbatch $caronte/Modules/RAND/Phred_score.sh $my_array
