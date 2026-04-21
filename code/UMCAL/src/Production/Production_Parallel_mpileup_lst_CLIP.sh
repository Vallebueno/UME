#!/bin/bash

# Parallel Production from independant mpileup files Launcher in VBC CLIP

# Author: Miguel Vallebueno,
# Date: 2023-01-25

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p c
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --qos=short


###list to array handling
set -euo pipefail

KEY=$1
TASK_ID="${SLURM_ARRAY_TASK_ID:-${UME_TASK_ID:-1}}"

line=$(sed -n "${TASK_ID}"p "${KEY}")
my_array=($(echo $line | tr "\t" "\n"))



SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
srcDIR="${UME_SRC_DIR:-$(cd "${SCRIPT_DIR}/.." && pwd)}"


file="$2/${my_array[0]}"
list=$3
ref=$4
Odir=$5

echo "info  file<$file> list<$list> ref<$ref> odir<$Odir>"

perl "$srcDIR/Production/Production_Call_mpileup_list_OC_V9.pl" "$2" "${my_array[0]}" "$ref" "$list" "$Odir"

