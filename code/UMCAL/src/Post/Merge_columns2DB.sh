#!/bin/bash

# Author: Miguel Vallebueno,
# Date: 2021-02-22

############  Designed to run with C1 Nodes as they have High clock speeds


# === SLURM Config === for C1 Nodes(hi-clock speed),  qos short(hi-priority)

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p c                 #Node type(c1,c2,m2)
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --qos=medium
#SBATCH --array=1-1
#SBATCH --constraint=c1

ml build-env/f2021
ml pigz/2.6-gcccore-10.2.0
# === INSTRUCTIONS ===
#KEY=/scratch-cbe/users/miguel.vallebueno/ALL/prt.tmp.lst.tl
#KEY=/scratch-cbe/users/miguel.vallebueno/ALL2/vcf_tmptx.lst.tlone
KEY=/scratch-cbe/users/miguel.vallebueno/Production/files_mp.lst.tlone
line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${KEY})
K=$2
K="Discovery_2023_08_cof10_ALLC.ll"
#Head="/scratch-cbe/users/miguel.vallebueno/Prod_cofs_test/in_prod_mpileup_files.lst.hh.gz"
Head="/scratch-cbe/users/miguel.vallebueno/Production/files_mp.lst.hh.gz"

ulimit -n 10000
getconf OPEN_MAX
echo "<$1><$K>";

#cd $1

#rm UME_production.vcf.gz
#first="/scratch-cbe/users/miguel.vallebueno/MERGE/build06_2023.sort.ll"
#second="/scratch-cbe/users/miguel.vallebueno/Production/${K}"
#linex=${line/paste/$second}
#echo "<$linex> ";
eval $line

#second=in_prod_mpileup_files.lst.Merge.db.gz

oname=UME_PROD_2023_0818.vcf.gz

zcat ${Head} files_mp.lst.Merge.db.gz | pigz > ${oname}




#sbatch $caronte/Modules/Tassel/TasselTRG_D_Matrix_downsample_HETS.sh ${oname} X VCF

