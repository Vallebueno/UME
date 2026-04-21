#!/bin/bash

# Author: Miguel Vallebueno,
# Date: 2021-02-22

############  Designed to run with C1 Nodes as they have High clock speeds


# === SLURM Config === for C1 Nodes(hi-clock speed),  qos short(hi-priority)

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p c                #Node type(c1,c2,m2)
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --qos=medium

# === INSTRUCTIONS ===
KEY=$1
IDIR=$2
WDIR=$3

line=$(sed -n "$SLURM_ARRAY_TASK_ID"p ${KEY})
my_array=($(echo $line | tr "\t" "\n"))

echo "<Merge OL><$my_array><${my_array[0]}> <$IDIR> <$WDIR>"

cd ${IDIR}

#perl $caronte/UMCAL/src/Pre/VCF_Merge_HidraV25.pl /path/to/reference/Zm-B73-REFERENCE-NAM-5.0.fa.OL ${my_array[0]} /path/to/reference/Zm-B73-REFERENCE-NAM-5.0.bed $WDIR

cd ${WDIR}

gzip --fast ${my_array[0]}.tmpt

#filename="${filename%.gz}"
