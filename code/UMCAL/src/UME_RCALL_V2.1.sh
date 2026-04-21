#!/bin/bash

# UME Main
# Author: Miguel Vallebueno
# Date: 2023-01-16
#https://stackoverflow.com/questions/16483119/an-example-of-how-to-use-getopts-in-bash

set -euo pipefail

usage() { echo "Usage: $0 [-f <merged_db_or_keyfile>] [-c <chromosome>] [-l <Phredscore cutoff>] [-k <list_of_vcf_files.txt>] [-i <IN_DIR_MPILEUPS>] [-r <OL_REFERENCE>] [-d <Output_DIR>] [-x <clist | merge | discovery | production>]" 1>&2; exit 1; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
srcDIR="${UME_SRC_DIR:-$SCRIPT_DIR}"
f=""
c=""
i=""
k=""
d=""
r=""
l=""
x=""

require_file() {
    local path="$1"
    local label="$2"
    if [ ! -f "$path" ]; then
        echo "Missing ${label}: <$path>" 1>&2
        exit 1
    fi
}

run_production_locally() {
    local key_file="$1"
    local input_dir="$2"
    local discovery_ll="$3"
    local reference="$4"
    local output_dir="$5"
    local sample

    while IFS=$'\t' read -r sample _; do
        if [ -z "${sample}" ]; then
            continue
        fi
        perl "$srcDIR/Production/Production_Call_mpileup_list_OC_V9.pl" "$input_dir" "$sample" "$reference" "$discovery_ll" "$output_dir"
    done < "$key_file"
}

prepare_production_merge_inputs() {
    local key_file="$1"
    local reference="$2"

    # These are assembly helpers for the final merged production VCF-like file.
    # They are not scientific discovery outputs like the .ll site list.
    if [ ! -f "${key_file}.tlone" ]; then
        perl "${srcDIR}/Pre/Merge_list2paste3.pl" "${key_file}"
    fi

    if [ ! -f "${key_file}.hh" ] && [ ! -f "${key_file}.hh.gz" ]; then
        perl "${srcDIR}/Pre/Final_header_make.pl" "${key_file}" "${reference}"
    fi
}
while getopts ":f:i:c:l:x:k:r:d:" o; do
    case "${o}" in
        x)
            x=${OPTARG}
            case "${x}" in
                clist|merge|discovery|production) ;;
                *) usage ;;
            esac
            ;;
        f)
            f=${OPTARG}
            #((s == 45 || s == 90)) || usage
            ;;
        c)
            c=${OPTARG}
            ;;
        i)
            i=${OPTARG}
            ;;
        k)
            k=${OPTARG}
            ;;
       d)
            d=${OPTARG}
            ;;
      r)
            r=${OPTARG}
            ;;
        l)
            l=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

a=max
q=0
echo "###############      UME main      ###############"
echo "f = <${f}>"
echo "k = <${k}>"
echo "c = ${c}"
echo "l = ${l}"
echo "x = <${x}>"
echo "a = ${a}"
echo "q= ${q}"
echo "r= <${r}>"
echo "d= <${d}>"
echo "i= <${i}>"

echo "file to analize = <${f}>"

if [  "${x}" != "clist" ] && [ "${x}" != "merge" ] && [ "${x}" != "discovery" ] && [ "${x}" != "production" ] ; then
 usage
fi

if [  "${x}" == "clist" ] ; then

  if [ -z "${f:-}" ] || [ -z "${k:-}" ]; then
      usage
  fi

  require_file "${f}" "merged DB for empirical distributions"
  require_file "${k}" "input list"

  DISC_LIST="${k}.Disc.lst"

  echo "UME creating list ....."

  perl "${srcDIR}/Pre/MK_List2IN.pl" "${k}"
  perl "${srcDIR}/Pre/Merge_list2paste3.pl" "${k}"

   echo "UME creating header ....."

  perl "${srcDIR}/Pre/Final_header_make.pl" "${k}" "${r:-}"

 echo "UME creating empirical db ....."

   perl "$srcDIR/Discovery/EmpDist2PvalV3.pl" "${DISC_LIST}" "${f}"
#perl ${srcDIR}/Pre/create_input_list_from_merge.pl ${k}
fi


if [  "${x}" == "merge" ] ; then

  #if [ -z "${f}" ] || [ -z "${c}" ] || [ -z "${l}" ] || [ -z "${x}" ] || [ -z "${k}" ]; then
  #    usage
  #fi

  echo "UME creating list ....."

fi


if [  "${x}" == "discovery" ] ; then

   echo "UME Discovery..."

  if [ -z "${f:-}" ] || [ -z "${c:-}" ] || [ -z "${k:-}" ] || [ -z "${l:-}" ] || [ -z "${r:-}" ]; then
     usage
  fi

  require_file "${f}" "merged DB"
  require_file "${k}" "input list"
  require_file "${k}.Disc.lst" "discovery caller list (.Disc.lst)"
  require_file "${f}.pval2f" "empirical p-value table (.pval2f)"
  require_file "${r}" "reference"

  DISC_LIST="${k}.Disc.lst"
  UNION_DB="${f%.db.gz}.qual${q}.Union.${a}.dbx"
  if [ "${UNION_DB}" = "${f}.qual${q}.Union.${a}.dbx" ]; then
      UNION_DB="${f}.qual${q}.Union.${a}.dbx"
  fi
  KNOWN_LIST="${UNION_DB}.cof${l}.lst"
  ENRICHED_LIST="${KNOWN_LIST%.lst}.lst2"
  ALT_INPUT="${ENRICHED_LIST}"

  echo "UME UNION Preparing input dataset..."
   #perl $srcDIR/Discovery/DB2mono_poli_V3.pl ${f}  #split chromosome and mono and poli files
 echo "Done!!!"


  echo "UME UNION Discovery..."
   perl "$srcDIR/Discovery/DB_Union_V10_FAST.pl" "${f}" "${c}" "${a}" "${k}" "${q}"
  echo "Done!!!"


  echo "UME creating DB of known positions..."
    perl "$srcDIR/Discovery/Discovery_mklst_DB.pl" "${f}" "${c}" "${l}" "${a}" "${q}"
  echo "Done!!!"


  echo "UME enriching DB with known GBS diversity..."
    perl "$srcDIR/Discovery/Discovery_enrichList_GBS.pl" "${KNOWN_LIST}" "${c}"
  echo "Done!!!"

  if [ ! -f "${ALT_INPUT}" ]; then
      ALT_INPUT="${KNOWN_LIST}"
  fi

  echo "UME creating Discovery.ll file..."
    perl "$srcDIR/Discovery/Discovery_sort_alts.pl" "${ALT_INPUT}" "${r}"
  LL_FILE="${ALT_INPUT}"
  if [[ "${LL_FILE}" == *.lst2 ]]; then
      LL_FILE="${LL_FILE%.lst2}.ll"
  elif [[ "${LL_FILE}" == *.lst ]]; then
      LL_FILE="${LL_FILE%.lst}.ll"
  fi
  echo "UME discovery product (.ll) = <${LL_FILE}>"
  echo "Done!!!"


fi


if [  "${x}" == "production" ] ; then

  echo "UME Production..."
  if [ -z "${f:-}" ] || [ -z "${k:-}" ] || [ -z "${d:-}" ] || [ -z "${x:-}" ] || [ -z "${r:-}" ] || [ -z "${i:-}" ]; then
      usage
  fi

  require_file "${f}" "production key file"
  require_file "${k}" "discovery .ll site list"
  require_file "${r}" "reference"

  if [ ! -d "${i}" ]; then
      echo "Missing production input directory: <${i}>" 1>&2
      exit 1
  fi

  mkdir -p "${d}"

  echo "UME production site list (.ll) = <${k}>"
  echo "UME production merge helpers = <${f}.tlone> and <${f}.hh>"

  prepare_production_merge_inputs "${f}" "${r}"

LKEY=$(wc -l < "${f}") #number of files in keyfile.


if command -v sbatch >/dev/null 2>&1 && [ "${UME_USE_SLURM:-1}" != "0" ]; then
  sbatch -W --array=1-"${LKEY}" "$srcDIR/Production/Production_Parallel_mpileup_lst_CLIP.sh" "${f}" "${i}" "${k}" "${r}" "${d}"
else
  run_production_locally "${f}" "${i}" "${k}" "${r}" "${d}"
fi

echo "Done!!!"
echo "Checking outs....."

if command -v sbatch >/dev/null 2>&1 && [ "${UME_USE_SLURM:-1}" != "0" ]; then
  sbatch -W "$srcDIR/Production/CHECK_tmptOLS.sh" "${d}" "${k}" "${LKEY}"
  echo "Done!!!"
  echo "Merge..."
  sbatch -W "$srcDIR/Production/Merge_columns2DB_production.sh" "${f}" "${d}" "${k}"
else
  bash "$srcDIR/Production/CHECK_tmptOLS.sh" "${d}" "${k}" "${LKEY}"
  echo "Done!!!"
  echo "Merge..."
  UME_TASK_ID=1 bash "$srcDIR/Production/Merge_columns2DB_production.sh" "${f}" "${d}" "${k}"
fi

   echo "Done!!!"
   echo "DB2VCF"

#this needs its own sbatch louncher due ram efiency
#cd ${d}

    # perl $srcDIR/Production/DB2VCF_V8.pl ${f} ${k} ${r} ${d}


 #  echo "UME Gzipping..."
 #  namex="${f}.${c}.poli.cof${l}.Union.V8.${a}.production.db"
 #  name=${namex/.db.gz/}
 #  gzip ${name}
 #  echo "Done 1!!!"

 #  namex="${f}.${c}.poli.qual${q}.Union.V8.${a}.dbx"
   #name=${namex/.db.gz/}
  # gzip ${name}
  # echo "Done!!!"

fi

 echo "UME Success !!!"
