"""
Here we invoke the qc_samples.sh script to run the BBTools suite, a fast, multithreaded bioinformatics tools set designed for analysis of DNA and RNA sequence data.
NOTE: Make sure that bbtools is installed on your system, see: https://jgi.doe.gov/data-and-tools/bbtools/

Author: jayson.gutierrez@vliz.be
"""

import glob
from subprocess import Popen, call, STDOUT, PIPE
from utilities import run_shell, seq_len_filter

def run_preproc_pipeline():
    '''Run entire bbtools pipeline to preprocess a bunch of fastq files and obtain QC seq reads. See predefined parameters in qc_samples.sh file'''
    
    ###List fastq files
    sample_cont = []
    for i in glob.glob("*fastq.gz"):
        sample_id = i.split("_")[0]
        if(sample_id not in sample_cont):
            sample_cont.append(sample_id)
            #print(sample_id)

    ###Run bbtools over each sample (.fastq.gz file) fetched from a public repository, e.g. ENA
    for s in sample_cont:
        fastq1, fastq2  = sorted(glob.glob("{}*.fastq.gz".format(s)),key=lambda s: int(s.split(".fastq")[0].split("_")[-1]))
        print("Initialize running bbtols on sample: {}".format(s))
        cmd = "bash qc_samples.sh {} {}".format(fastq1, fastq2)
        print(cmd)
        run_shell(cmd, stdout=True)
        #Now filter seq reads by lenght
        qc_reads_infile = '{}_Qual_Filt.fastq'.format(s)
        qc_reads_outfile = qc_reads_infile.replace("Filt","Filt2")
        seq_len_filter(qc_reads_infile,qc_reads_outfile,min_seq_len=100)
        #Remove unfiltered file and rename filtered one
        cmd2 = "rm {0}; mv {1} {0}".format(qc_reads_infile,qc_reads_outfile)
        run_shell(cmd2, stdout=False)
        print("Finished running bbtols on sample: {}".format(s))
        print("---------------------------------------------------------")

if __name__ == "__main__":
    run_preproc_pipeline()