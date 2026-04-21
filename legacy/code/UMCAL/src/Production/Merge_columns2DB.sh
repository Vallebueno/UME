#!/bin/bash

# Author: Miguel Vallebueno,
# Date: 2021-02-22

############  Designed to run with C1 Nodes as they have High clock speeds


# === SLURM Config === for C1 Nodes(hi-clock speed),  qos short(hi-priority)

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p c                 #Node type(c1,c2,m2)
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --qos=long

ml build-env/f2021
ml pigz/2.6-gcccore-10.2.0
# === INSTRUCTIONS ===

#KEY="/path/to/workdir/Downsampling2/mpileup.lst.tlone"

KEY=${1}.tlone

line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${KEY})

ulimit -n 10000
getconf OPEN_MAX

echo "<$line>";
#cd /path/to/workdir/Production
eval $line

head=${1}.hh.gz
oname="UME_production_2023_08_24.vcf.gz"

zcat ${head} files_mp.lst.Merge.db.gz | pigz > ${oname}

perl $caronte/Modules/VCF/VCF_RANDOM_SAMPLE.pl ${oname}

sbatch $caronte/Modules/Tassel/TasselTRG_D_Matrix_downsample_HETS.sh ${oname}.smp.vcf X VCF

