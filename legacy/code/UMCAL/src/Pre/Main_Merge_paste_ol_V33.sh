#!/bin/bash

# SCRIPT to MERGE VCFs into a multi VCF file
#it uses the paste strategy with is the fastest


# Author: Miguel Vallebueno
# Date: 2022-03-21

# === SLURM Config ===

#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH -p c                 #Node type(c1,c2,m2)
#SBATCH --mem=10G
#SBATCH --qos=long         #-time +priority

# === Load Modules ===
ml build-env/f2021
ml pigz/2.6-gcccore-10.2.0

IDIR="/path/to/gvcf_inputs"
WDIR="/path/to/workdir/MERGE"
KEY="/path/to/gvcf_inputs/files.mergein_Discovery_2023_05.lst"

###make list of VCF files ordered by sample name
#find ${WDIR} -type f -name *.vcf > files.mergein.lst

###paralel run for each sample into a one column sample
#
#sbatch --array=1-2457 $caronte/UMCAL/src/Pre/Merge_OL_parallel.sh $KEY $IDIR $WDIR
#sbatch --array=422-2457 $caronte/UMCAL/src/Pre/Merge_OL_parallel.sh $KEY $IDIR $WDIR
#sbatch --array=344-421 $caronte/UMCAL/src/Pre/Merge_OL_parallel.sh $KEY $IDIR $WDIR
#exit

perl $caronte/UMCAL/src/Pre/CHECK_ll.pl

##make innput line
perl $caronte/UMCAL/src/Pre/Merge_list2paste3.pl $KEY

#### merge using paste
sbatch -W $caronte/UMCAL/src/Pre/Merge_columns2DB.sh $KEY $WDIR

### put back the coordinates and split by cromosome

perl $caronte/Modules/Hidra-Heracles/DB2mono_poli_V3.pl /path/to/reference/Zm-B73-REFERENCE-NAM-5.0.fa.OL ${WDIR}/All.tlone.Merge.db /path/to/reference/Zm-B73-REFERENCE-NAM-5.0.bed $KEY.tlone
perl $srcDIR/Discovery/DB2mono_poli_V3.pl ${f}  #split chromosome and mono and poli files




gzip --fast All.tlone.Merge.db.1.poli.db
gzip --fast All.tlone.Merge.db.2.poli.db
gzip --fast All.tlone.Merge.db.3.poli.db
gzip --fast All.tlone.Merge.db.4.poli.db
gzip --fast All.tlone.Merge.db.5.poli.db
gzip --fast All.tlone.Merge.db.6.poli.db
gzip --fast All.tlone.Merge.db.7.poli.db
gzip --fast All.tlone.Merge.db.8.poli.db
gzip --fast All.tlone.Merge.db.9.poli.db
gzip --fast All.tlone.Merge.db.10.poli.db

#sbatch $caronte/Modules/VCF/Main_union_parallel.sh





