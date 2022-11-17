#!/bin/bash -l
###Script to run Augustus on a bunch of genomes. This script submit each job independently to the HPC-UGent (see: https://www.vscentrum.be/)

LOWER=$1 #e.g. LOWER=160
UPPER=$2 #e.g. UPPER=319
SPECIES_LIST=all_euk_genomes_list.csv #Modify this accordingly! List of species-associated genome assemblies to be processed using AUGUSTUS
###Cast list of job submission scripts into array form
readarray -t all_subjobs < ${SPECIES_LIST}
for (( i = $LOWER; i < $UPPER; i++ )); do
    SPECIES=${all_subjobs[i]}
    python augustus_orf_predict.py ${SPECIES}
done
