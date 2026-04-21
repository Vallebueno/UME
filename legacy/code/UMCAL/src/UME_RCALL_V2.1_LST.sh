#!/bin/bash

# UME Main
# Author: Miguel Vallebueno
# Date: 2023-01-16
#https://stackoverflow.com/questions/16483119/an-example-of-how-to-use-getopts-in-bash

usage() { echo "Usage: $0 [-f <Merged.gvcf>] [-c <chromosome>] [-l <Phredscore cutoff>] [-k <list_of_vcf_files.txt>] [-x <clist | emp | discovery | production>]" 1>&2; exit 1; }
while getopts ":f:c:l:x:k:" o; do
    case "${o}" in
        x)
            x=${OPTARG}
            ((${x} == "clist" || ${x} == "merge" || ${x} == "discovery" || ${x} == "production"))  || usage
            ;;
        f)
            f=${OPTARG}
            #((s == 45 || s == 90)) || usage
            ;;
        c)
            c=${OPTARG}
            ;;
        k)
            k=${OPTARG}
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
echo "k = ${k}"
echo "c = ${c}"
echo "l = ${l}"
echo "x = ${x}"
echo "a = ${a}"
echo "q= ${q}"

echo "file to analize = <${fa}>"

srcDIR="/path/to/legacy/Caronte/code/UMCAL/src"

if [  ${x} != "clist" ] && [ ${x} != "merge" ] && [ ${x} != "discovery" ] && [ ${x} != "production" ] ; then
 usage
fi

if [  ${x} == "clist" ] ; then

  if [ -z "${k}" ]; then
      usage
  fi

  echo "UME creating list ....."

  perl ${srcDIR}/Pre/MK_List2IN.pl ${k}
  perl ${srcDIR}/Pre/Merge_list2paste3.pl ${k}

 echo "UME creating empirical db ....."

   perl $srcDIR/Discovery/EmpDist2PvalV3.pl ${k} ${f}
#perl ${srcDIR}/Pre/create_input_list_from_merge.pl ${k}
fi


if [  ${x} == "merge" ] ; then

  #if [ -z "${f}" ] || [ -z "${c}" ] || [ -z "${l}" ] || [ -z "${x}" ] || [ -z "${k}" ]; then
  #    usage
  #fi

  echo "UME creating list ....."

fi


if [  ${x} == "discovery" ] ; then

   echo "UME Discovery..."

  if [ -z "${f}" ] || [ -z "${c}" ] ||  [ -z "${k}" ]; then
     usage
  fi

  echo "UME UNION Preparing input dataset..."
   #perl $srcDIR/Discovery/DB2mono_poli_V3.pl ${f}  #split chromosome and mono and poli files
 echo "Done!!!"


  echo "UME UNION Discovery..."
  #perl $srcDIR/Discovery/DB_Union_V10_FAST_LST.pl ${f} ${c} ${a} ${k} ${q}
  echo "Done!!!"


  echo "UME creating DB of known positions..."
    perl $srcDIR/Discovery/Discovery_mklst_DB_V2.pl ${f} ${c} ${l} ${a} ${q}
  echo "Done!!!"



  echo "UME enriching DB with known GBS diversity..."
    #perl $srcDIR/Discovery/Discovery_enrichList_GBS.pl ${f} ${c} ${l} ${a} ${q}
  echo "Done!!!"


 echo "Done 2!!!"
fi


if [  ${x} == "production" ] ; then

  echo "UME Production..."
  #if [ -z "${f}" ] || [ -z "${c}" ] || [ -z "${l}" ] || [ -z "${x}" ] || [ -z "${k}" ]; then
  #    usage
  #fi


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
   echo "Done!!!"

fi

 echo "UME Success !!!"
