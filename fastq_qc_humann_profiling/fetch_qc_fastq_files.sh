#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00

###-------------------------------------------------------------------------------------------------------------------------------------
###Run script to fetch from a public repository a bunch of fastq files from a bioproj and then run bbtools on them to obtain QC seq reads
###-------------------------------------------------------------------------------------------------------------------------------------

###Run script on already QC reads files: removing very short reads; i.e. < 100 nucleotides
#cd /data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/HUMAnN3/data/reads/Elferink_etal_2020_mtx_data

###Activate appropriate conda env
#source activate /home/VLIZ2000/jayson.gutierrez/anaconda3/envs/biobakery3

###Fetch samples: these must be stored in a .txt file, e.g. the PRJEB37134_mtx_samples.txt contains:
###ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR398/004/ERR3980584/ERR3980584_1.fastq.gz
###ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR398/004/ERR3980584/ERR3980584_2.fastq.gz
###ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR398/009/ERR3980589/ERR3980589_1.fastq.gz
###ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR398/009/ERR3980589/ERR3980589_2.fastq.gz
python utilities.py PRJEB37134_mtx_samples.txt

###Run bbtools: this python script invokes qc_samples.sh to automatize the QC process
python run_bbtools_preprocess.py
