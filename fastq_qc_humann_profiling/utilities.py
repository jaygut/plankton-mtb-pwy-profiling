"""Utility functions are being invoked within other scripts

Author: jayson.gutierrez@vliz.be
"""

from subprocess import Popen, call, STDOUT, PIPE
import pandas as pd
import os
import sys
import Bio.SeqIO as bioseqio


def run_shell(cmd, stdout=True):
    '''Simple function to invoke and run command-line tools
    
       --cmd: command line args in string format
    '''
    
    if(stdout):
        return Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE).communicate()[0].decode("utf-8")
    else:
        Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE).communicate()[0]


def fetch_fastq_files(fid):
    '''Function to fetch fastq files from an FTP server, e.g. ftp://ftp.sra.ebi.ac.uk

       --fid: name of file with a list of ftp addresses to the fastq files bundleded  
              in a given bioproj, e.g. PRJEB37134. Such list can be obtained from:
              https://sra-explorer.info/
    '''
    
    #Load in list of ftp addresses to fastq files
    fastq_list = pd.read_csv(fid,header=None).values.flatten()

    #Fetch all the fastq files in the list
    for fsq in fastq_list:
        print("Fetching file from: {}\n".format(fsq))
        cmd = "curl -L {} -o {}".format(fsq, os.path.basename(fsq))
        run_shell(cmd)


def seq_len_filter(fid, outf, min_seq_len=100):
    """Function to filter reads from a QC based on seq lenght"""

    with open(fid,"r") as input_handle:
        for r in bioseqio.parse(input_handle, "fastq"):
            if len(r) >= min_seq_len:
                with open(outf, "a") as output_handle:
                    bioseqio.write(r, output_handle, "fastq")


if __name__ == "__main__":
    fid = sys.argv[1] #e.g. PRJEB37134_mtx_samples.txt
    fetch_fastq_files(fid)
