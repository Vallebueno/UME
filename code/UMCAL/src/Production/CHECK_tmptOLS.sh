#!/bin/bash

# mv files SLURM job script
#
# Author: Miguel Vallebueno,
# Date: 2021-02-22

############  Designed to run with C1 Nodes as they have High clock speeds


# === SLURM Config === for C1 Nodes(hi-clock speed),  qos short(hi-priority)

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p c
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --qos=medium


set -euo pipefail

if command -v ml >/dev/null 2>&1; then
    ml build-env/f2022
    ml pigz/2.6-gcccore-11.2.0
fi

# === INSTRUCTIONS ===


WDIR=$1
LLST=$2
RUNS=$3


#cd /scratch-cbe/users/miguel.vallebueno/Production/cof5
#cd /scratch-cbe/users/miguel.vallebueno/Production/GBS
#cd /scratch-cbe/users/miguel.vallebueno/Production/HM3

cd "$WDIR"

shopt -s nullglob
tmpt_files=( *.tmpt.gz )
shopt -u nullglob

: > CCU.lst
for x in "${tmpt_files[@]}"; do
    if command -v unpigz >/dev/null 2>&1; then
        unpigz -c "$x" | wc -l >> CCU.lst
    else
        gzip -dc "$x" | wc -l >> CCU.lst
    fi
    echo "$x" >> CCU.lst
done


LKEY=$(wc -l < "${LLST}")
LCCU=$(wc -l CCU.lst)

QQ=$(grep -c "${LKEY}" CCU.lst)





if [ "${#tmpt_files[@]}" -eq 0 ]; then
    echo "FATAL ERROR: no tmpt.gz files were found in <$WDIR>"
elif [ "$QQ" = "$RUNS" ]; then
    echo "ALL seems fine!!!!; check if you got either a 0/0 or 1/1 s at the files"
else
    echo "FATAL ERROR: some tmpt.gz files have diferent lenghts check CCU.lst for quantification "
fi


