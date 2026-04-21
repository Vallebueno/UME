#!/bin/bash

# Author: Miguel Vallebueno,
# Date: 2021-02-22

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p c
#SBATCH --cpus-per-task=20
#SBATCH --mem=10G
#SBATCH --qos=long
#SBATCH --array=1-1

set -euo pipefail

if command -v ml >/dev/null 2>&1; then
    ml build-env/f2021
    ml pigz/2.6-gcccore-10.2.0
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_ROOT="${UME_REPO_ROOT:-${UME_CARONTE_ROOT:-$(cd "${SCRIPT_DIR}/../../.." && pwd)}}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-${UME_TASK_ID:-1}}"
MERGE_PREFIX="${UME_PRODUCTION_PREFIX:-MaizeDB_UME}"
F=$1
KEY=${F}.tlone
HeadGz=${F}.hh.gz
HeadPlain=${F}.hh
WDIR=$2
K=$3

if [ -f "${HeadGz}" ]; then
    Head="${HeadGz}"
    HEAD_CMD="gzip -dc \"${Head}\""
elif [ -f "${HeadPlain}" ]; then
    Head="${HeadPlain}"
    HEAD_CMD="cat \"${Head}\""
else
    echo "Missing header file. Expected ${HeadGz} or ${HeadPlain}" 1>&2
    exit 1
fi

echo "<$F> <$KEY> <$Head> <$WDIR> <$K>"

line=$(sed -n "${TASK_ID}"p "${KEY}")

second="paste <(cat ${K})"
linex=${line/paste/$second}
MERGED_TMP="${MERGE_PREFIX}_pre.vcf.gz"
MERGED_OUT="${MERGE_PREFIX}_production.vcf.gz"
SOURCE_MERGE="${F}.Merge.db.gz"
MERGE_CMD="gzip"

#linex=${linex/All.tlone.Merge.db.gz/$second}
linex=$(echo "$linex" | sed 's/mpileup.gz/mpileup.gz.tmpt.gz/g')
linex=${linex/${SOURCE_MERGE}/${MERGED_TMP}}
if command -v pigz >/dev/null 2>&1; then
    linex=${linex/| gzip >/| pigz >}
    MERGE_CMD="pigz"
fi

cd "$WDIR"

ulimit -n 10000
getconf OPEN_MAX
#echo "$linex";
echo "<$linex>";

eval "$linex"

oname=${MERGED_OUT}

eval "${HEAD_CMD} | cat - <(gzip -dc \"${MERGED_TMP}\") | ${MERGE_CMD} > \"${oname}\""

rm -f "${MERGED_TMP}"

if [ -f "${CODE_ROOT}/Modules/VCF/VCF_RANDOM_SAMPLE.pl" ]; then
    perl "${CODE_ROOT}/Modules/VCF/VCF_RANDOM_SAMPLE.pl" "${oname}"
fi

if [ -f "${CODE_ROOT}/Modules/Tassel/TasselTRG_D_Matrix_downsample_HETS.sh" ] && command -v sbatch >/dev/null 2>&1; then
    sbatch "${CODE_ROOT}/Modules/Tassel/TasselTRG_D_Matrix_downsample_HETS.sh" "${oname}.smp.vcf" X VCF
fi
