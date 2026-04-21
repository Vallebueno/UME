#!/bin/bash

# Author: Miguel Vallebueno,
# Date: 2021-02-22

############  Designed to run with C1 Nodes as they have High clock speeds


# === SLURM Config === for C1 Nodes(hi-clock speed),  qos short(hi-priority)

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p c                 #Node type(c1,c2,m2)
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --qos=long
#SBATCH --array=1-1
#SBATCH --constraint=c1

ml build-env/f2022
ml pigz/2.6-gcccore-11.2.0
# === INSTRUCTIONS ===

K=$1
KEY=${K}.tlone
WDIR=$2

line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${KEY})

cd $WDIR

ulimit -n 3000
#getconf OPEN_MAX

echo "<$line>";

eval $line
