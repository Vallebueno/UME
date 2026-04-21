#!/bin/bash

# Author: Miguel Vallebueno,
# Date: 2021-02-22

############  Designed to run with C1 Nodes as they have High clock speeds


# === SLURM Config === for C1 Nodes(hi-clock speed),  qos short(hi-priority)

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p c                 #Node type(c1,c2,m2)
#SBATCH --cpus-per-task=20
#SBATCH --mem=10G
#SBATCH --qos=long
#SBATCH --time=7-10:00:00
#SBATCH --array=1-1
#SBATCH --constraint=c1

ml build-env/f2021
ml pigz/2.6-gcccore-10.2.0
# === INSTRUCTIONS ===
#KEY=/scratch-cbe/users/miguel.vallebueno/ALL/prt.tmp.lst.tl
#KEY=/scratch-cbe/users/miguel.vallebueno/ALL2/vcf_tmptx.lst.tlone
KEY=/groups/swarts/lab/MAVE/Experiments/Hydra_coverage/files.mergein.lst.tlone
line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${KEY})

ulimit -n 10000
getconf OPEN_MAX

echo "<$line>";

eval $line



