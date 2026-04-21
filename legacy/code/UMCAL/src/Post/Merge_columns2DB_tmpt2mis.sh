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
#SBATCH --qos=short
#SBATCH --constraint=c1
#SBATCH --array=1-1
#SBATCH --job-name=MERGE


ml build-env/f2021
ml pigz/2.6-gcccore-10.2.0
# === INSTRUCTIONS ===
#KEY=/path/to/workdir/ALL/prt.tmp.lst.tl
#KEY=/path/to/workdir/ALL2/vcf_tmptx.lst.tlone
#KEY=/path/to/project_inputs/files.mergein.lst.
DST="HM3"

F=/path/to/mpileup/in_prod_mpileup_files.lst
KEY=${F}.tlone
Head=${F}.hh.gz
WDIR=/path/to/workdir/Production/${DST}
K=/path/to/workdir/Production/${DST}/Discovery_LL.ALLC.${DST}.ll
#K=/path/to/workdir/Production/${DST}/ZeaGBSv27_publicSamples_raw_AGPv5_crossmap.ll

echo "<$F> <$KEY> <$Head> <$WDIR> <$K>"


line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${KEY})



second="paste <(cat ${K})"
linex=${line/paste/$second}



K=${K/.ll/}

my_array=($(echo $K | tr "/" "\n"))
second="${my_array[-1]}.UME_pres2m.vcf.gz"

linex=${linex/All.tlone.Merge.db.gz/$second}
linex=${linex//tmpt.gz/tmpt.gz.s2m.gz}

cd $WDIR

ulimit -n 10000
getconf OPEN_MAX
#echo "$linex";
echo "<$linex>";

eval $linex


oname=${second/UME_pres2m/UME_productions2m}

#gzip ${Head}

zcat ${Head} ${second} | pigz > ${oname}




