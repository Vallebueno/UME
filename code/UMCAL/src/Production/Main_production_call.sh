#!/usr/bin/env bash
# UMCAL Production SLURM job script

# uses a list of genotypes per position and extract the genotypes from the reads (mpileup) data for each individual creating a VCF

# Author: Miguel Vallebueno
# Date: 2022-04-29

# === SLURM Config ===

#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -p c
#SBATCH --mem=10G
#SBATCH --qos=medium
#SBATCH --array=1-120

# === Load Modules ===

GP="/groups/swarts/lab/Resources/Scripts/Caronte/code"
PACK="UMCAL/src/Production/"

# === INPUT Config ===

IDIR="/scratch-cbe/users/miguel.vallebueno/mpilp2db"
KEY="/scratch-cbe/users/miguel.vallebueno/mpilp2db/mpileups.lst"
REF="/groups/swarts/lab/Resources/ReferenceGenomes/Zm-B73-REFERENCE-NAM-5.0.fa.OL"

ODIR="/scratch-cbe/users/miguel.vallebueno/mpilp2db"
#LST="/groups/swarts/lab/MAVE/Experiments/Hydra_coverage/OL/files.mergein3.lst.tlone.Merge.db.7.poli.qual0.Union.V7-1.max.db.lst"
LST="/scratch-cbe/users/miguel.vallebueno/mpilp2db/Maize_GBS_Kriztian_build.lst"

###      >>>> LIST CONTROL ARRAY<<<<

line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${KEY})
my_array=($(echo $line | tr "\t" "\n"))

acc=${my_array[1]}
acc="${acc%%[[:cntrl:]]}" #### remove weird EOF stuff
name=$acc

file=${my_array[0]}


# === INSTRUCTIONS ===

### One column per file

### PASTE

#### prod
echo "<<Caronte><HYDRA>>      <$file>   <$name>   <$REF> <$IDIR> <$ODIR>"





#perl ${GP}/${PACK}/Production_call_OCflist_V3.pl $file $REF $LST
#perl ${GP}/${PACK}/Production_Call_mpileup_list_OC_V4.pl $file $REF $LST $ODIR
#gzip ${file}.ol

### gzip new OL
### make list


#perl ${GP}/${PACK}/Merge_list2paste3.pl files.txt


#excec
sbatch {GP}/${PACK}/Merge_columns2DB.sh files.txt

DB2VCF








