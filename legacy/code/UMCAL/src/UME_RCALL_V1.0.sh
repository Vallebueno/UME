#!/bin/bash

# UME Main
# Author: Miguel Vallebueno
# Date: 2022-07-29
#https://stackoverflow.com/questions/16483119/an-example-of-how-to-use-getopts-in-bash

usage() { echo "Usage: $0 [-f <Merged.gvcf>] [-c <chromosome>] [-l <Phredscore cutoff>] [-x <clist | emp | union>]" 1>&2; exit 1; }
while getopts ":f:c:l:x:" o; do
    case "${o}" in
        x)
            x=${OPTARG}
            ((${x} == "clist" || ${x} == "emp" || ${x} == "union")) || usage
            ;;
        f)
            f=${OPTARG}
            #((s == 45 || s == 90)) || usage
            ;;
        c)
            c=${OPTARG}
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
if [ -z "${f}" ] || [ -z "${c}" ] || [ -z "${l}" ] || [ -z "${x}" ] ; then
    usage
fi
a=max
q=0
echo "###############      UME main      ###############"
echo "f = <${f}>"
echo "c = ${c}"
echo "l = ${l}"
echo "x = ${x}"
echo "a = ${a}"
echo "q= ${q}"

srcDIR="/path/to/legacy/Caronte/code/UMCAL/src"

if [  ${x} != "clist" ] && [ ${x} != "emp" ] && [ ${x} != "union" ] ; then
 usage
fi

if [  ${x} == "clist" ] ; then

  echo "UME creating list ....."
  perl ${srcDIR}/Discovery/create_input_list_from_merge.pl ${f}

fi


if [  ${x} == "emp" ] ; then

  echo "UME creating empirical db ....."
   #$srcDIR/Discovery/EmpDist2PvalV2.pl

fi


if [  ${x} == "union" ] ; then

  echo "UME UNION Preparing input dataset..."
   #perl $srcDIR/Discovery/DB2mono_poli_V3.pl ${f}   #split chromosome and mono and poli files
 echo "Done!!!"

  echo "UME UNION Discovery..."
   #perl $srcDIR/Discovery/DB_Union_V8_FAST.pl ${f} ${c} ${a} /path/to/phred/Emp_Dist_CALLERS.pval2f ${q}
  echo "Done!!!"

  echo "UME creating DB of known positions..."
    perl $srcDIR/Discovery/Discovery_mklst_DB.pl ${f} ${c} ${l} ${a} ${q}
  echo "Done!!!"

  echo "UME enriching DB with known GBS diversity..."
    perl $srcDIR/Discovery/Discovery_enrichList_GBS.pl ${f} ${c} ${l} ${a} ${q}
  echo "Done!!!"

  echo "UME Production..."
   perl $srcDIR/Production/Production_call_from_DB.pl ${f} ${c} ${l} ${a} ${q}
  echo "Done!!!"
  echo "UME Gzipping..."
  namex="${f}.${c}.poli.cof${l}.Union.V8.${a}.production.db"
  name=${namex/.db.gz/}
  gzip ${name}
  echo "Done 1!!!"

  namex="${f}.${c}.poli.qual${q}.Union.V8.${a}.dbx"
  name=${namex/.db.gz/}
  gzip ${name}
 echo "Done 2!!!"
fi

 echo "UME Success !!!"
