#!/usr/bin/env python

"""Script to annotate the TARA SMAGs protein DB against the UniRef90 DB using diamond

   Author: jayson.gutierrez@vliz.be
"""

import numpy as np
from subprocess import Popen, call, STDOUT, PIPE

def runShell(cmd):
    process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = process.communicate()   
    return out, err

shell_script_template = '''#!/bin/bash
TARA_FASTA=/kyukon/data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/Annotation/Data/TARA_SMAGs/
UNIREF_DB=/kyukon/data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/HUMAnN3/data/databases/uniref/uniref90_201901b_full.dmnd
QUERY_FASTA=${TARA_FASTA}"%(query_fast)s"
SEQ_ANNOT_HITS=${TARA_FASTA}"%(matches_uniref_annots)s"
SEQ_ALIGMENT_HITS=${TARA_FASTA}"%(matches_annot_seqs)s"

diamond blastp -d ${UNIREF_DB} -q ${QUERY_FASTA} -o ${SEQ_ANNOT_HITS} --evalue 1e-3 --query-cover 80 --id 50 --max-target-seqs 1 --block-size 4 --index-chunks 2 --threads 16 --header --al ${SEQ_ALIGMENT_HITS} --alfmt fasta --quiet'''

submit_job_script_template = '''#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=60gb
#PBS -l walltime=72:00:00

cd /data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/Annotation/Scripts/Diamond_runs
module --force purge
source activate /data/gent/vo/001/gvo00125/vsc43582/Software/dmnd_annotation

bash %(shell_script_name)s
'''

def main():
    #Run diamond over the different seq chunks: annotating seqs against UniRef90
    num_chunks = np.arange(0,10)
    for n in num_chunks: 

        #Name different input/output files
        query_fast="SMAGs_filtered_chunk_{}.faa".format(str(n))
        matches_uniref_annots="SMAGs_filtered_chunk_annots_{}.tsv".format(str(n))
        matches_annot_seqs="SMAGs_filtered_chunk_annot_seqs_{}.faa".format(str(n))
        shell_script_name = "run_filtered_smags_chunk_{}.sh".format(str(n))
        submit_job_script_name = "submit_dmnd_annot_filtered_smags_chunk_{}.sh".format(str(n))

        #Create shell script for running diamnond on a given chunk of sequences from the TARA SMAGs
        shell_text = shell_script_template % {'query_fast':query_fast,
                                              'matches_uniref_annots':matches_uniref_annots,
                                              'matches_annot_seqs':matches_annot_seqs}

        with open(shell_script_name,"w") as ss:
            ss.write(shell_text)

        #Create submit job script and place job in the queue
        submit_job_script = submit_job_script_template % {'shell_script_name':shell_script_name}

        with open(submit_job_script_name,"w") as ss:
                ss.write(submit_job_script)

        runShell("qsub {}".format(submit_job_script_name))
    
if __name__ == "__main__":
    
    ###NOTE: Run command to split up original TARA SMAGs DB into smaller chunks of sequence batches
    # awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%1000==0){file=sprintf("SMAGs_filtered_chunk%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < ${TARA_SMAGs}
    main()