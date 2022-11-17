#!/usr/bin/env python
# oding: utf-8

"""Collection of functions developed to create a wrapper to run BUSCO on a set of genomes.  
   These functions are not only suitable for the identification of BUSCO genes, but also for extracting hints necessary to run Augustus in a post-BUSCO run,
   which can help us refine the initial ab initio structural gene prediction.

   Author: jayson.gutierrez@vliz.be
"""

from subprocess import Popen, call, STDOUT, PIPE
import time
import sys
import pandas as pd
import os
import glob
import re
import argparse

from keep_busco_outdirs_lightweight import *

def set_submit_busco_job_script(GENOME_ID = "/Genomes/genbank/GCA_900660405.1/GCA_900660405.1_ASM90066040v1_genomic.fna",
                                SPECIES_NAME = "Pseudo-nitzschia multistriata",
                                LINEAGE_DB = "eukaryota_odb10",
                                GENE_PREDICTOR='augustus',
                                NUM_PROCESSORS_PERNODE = "10",
                                WALL_TIME="24", SUBMIT_JOB=False):
    
    '''Function to create two scripts:
        1) a script, shell_script_name, which sets the parameters to run BUSCO on a given genome assembly
        2) a script, , used to submit the job to the HPC cluster'''
    
    SPECIES_NAME = SPECIES_NAME.replace(" ","_")
        
    shell_script_template = '''#!/bin/bash  

###script name: run_busco.sh 

#source activate /data/gent/vo/001/gvo00125/vsc43582/Software/annotation 

###$1 can be set to e.g. "stramenopiles_odb10" 
###$2 can be set to e.g. "augustus" or "metaeuk" 

###Testing busco with genome assembly
INPUT="%(GENOME_ID)s" 

MODE="genome" 

OUTPUT="%(SPECIES_NAME)s" 

OUTPUT_PATH="./Results_Default" 

CONFIG="./config/config.ini" 

if [ "$1" == "auto" ]; then 

    busco --in ${INPUT} --mode ${MODE} --out ${OUTPUT} --out_path ${OUTPUT_PATH} --auto-lineage-euk --config ${CONFIG} --quiet --force --cpu %(NUM_PROCESSORS_PERNODE)s 

else ###Passing in name of lineage dataset 
   
   ###Run BUSCO with default gene predictor: MetaEuk 
   if [ "$2" == "metaeuk"  ]; then 
      
      busco --in ${INPUT} --mode ${MODE} --out ${OUTPUT} --out_path ${OUTPUT_PATH} --lineage_dataset $1 --config ${CONFIG} --quiet --force --cpu %(NUM_PROCESSORS_PERNODE)s 
   
   ###Else run Busco with Augustus 
   else 
     OUTPUT_PATH="./Results_Augustus" 
     CONFIG="./config/config_augustus.ini" 
     export PATH="/data/gent/vo/001/gvo00125/vsc43582/Software/annotation/bin:$PATH" 
     export PATH="/data/gent/vo/001/gvo00125/vsc43582/Software/annotation/bin:$PATH" 
     export AUGUSTUS_CONFIG_PATH="/data/gent/vo/001/gvo00125/vsc43582/Software/annotation/config/" 
     busco --in ${INPUT} --mode ${MODE} --out ${OUTPUT} --out_path ${OUTPUT_PATH} --lineage_dataset $1 --config ${CONFIG} --augustus --long --quiet --force --cpu %(NUM_PROCESSORS_PERNODE)s 
   fi
   
fi'''
    
    shell_text =  shell_script_template % {'GENOME_ID':GENOME_ID,\
                                           'SPECIES_NAME':SPECIES_NAME,\
                                           'NUM_PROCESSORS_PERNODE':NUM_PROCESSORS_PERNODE}
    
    shell_script_name = "run_busco_{}.sh".format(SPECIES_NAME.replace("-",'_').lower())
    
    with open(shell_script_name,"w") as ss:
        ss.write(shell_text)
        
    submit_job_script_template = '''#!/bin/bash -l

#PBS -l nodes=1:ppn=%(NUM_PROCESSORS_PERNODE)s 
#PBS -l walltime=%(WALL_TIME)s:00:00

cd $PBS_O_WORKDIR

module --force purge
module load cluster/swalot

source activate /data/gent/vo/001/gvo00125/vsc43582/Software/annotation

###Set name of lineage dataset. options: 1) "auto";  or 2) any odb10 e.g. "eukaryota_odb10"
LD="%(LINEAGE_DB)s"

###Set gene predictor algorithm. Options: 1) "metaeuk" (default); or 2) "augustus"
GP="%(GENE_PREDICTOR)s" 

bash %(shell_script_name)s ${LD} ${GP}'''


    submit_job_script = submit_job_script_template % {'NUM_PROCESSORS_PERNODE':NUM_PROCESSORS_PERNODE,\
                                                      'WALL_TIME':WALL_TIME,\
                                                      'LINEAGE_DB':LINEAGE_DB,\
                                                      'GENE_PREDICTOR':GENE_PREDICTOR,\
                                                      'shell_script_name':shell_script_name}

    submit_job_script_name = shell_script_name.replace("run_","submit_")
    
    with open(submit_job_script_name,"w") as ss:
        ss.write(submit_job_script)
            

    if(SUBMIT_JOB):
        #print("##########################################")
        print("Submitting job via: {}".format(submit_job_script_name))
        #Submit job to the HPC cluster
        runShell("qsub {}".format(submit_job_script_name))
    

def submit_busco_job(sps_name,ftp_path,lineage_db,download_genome=True,submit_job=False):
    '''Function to submit a BUSCO job based on the scripts created by set_submit_busco_job_script.
       Params:
           - sps_name: name of species whose genome will be processed using BUSCO
           - ftp_path: NCBI ftp to fetch genome assembly from
           - lineage_db: set busco lineage database to run the genome for BUSCO
    '''
    
    #Set path to dir to store genome
    genome_dir = "./Genomes/genbank/" + os.path.basename(os.path.dirname(ftp_path))
    #unzip file
    gunzip_fn = os.path.join(genome_dir,os.path.basename(ftp_path))
    
    if(download_genome):
        
        #Create dir if genome dir doesn't exist
        if(~os.path.exists(genome_dir)):
            runShell("mkdir -p {}".format(genome_dir))
        
        #Keep trying downloading genome file until valid format is fetched
        print("Fetching genome from: {}".format(ftp_path))
        err = "invalid"
        while(err!=''):
            #First remove anything from this dir
            runShell("rm {}".format(os.path.join(os.path.dirname(gunzip_fn),"*")))
            #Fetch genome assembly from ftp    
            runShell("wget {} -P {}".format(ftp_path, genome_dir))
            #Unzip genome
            out, err = runShell("gunzip {}".format(gunzip_fn))
            #If faulty .gz file exists then remove it!
            if(os.path.exists(gunzip_fn)):
                runShell("rm {}".format(gunzip_fn))

        #print("##########################################")
        print("Genome fetched!")
    #Run BUSCO for downloaded genome assembly and obtain hints for running Augustus afterwards
    unzipped_gf = gunzip_fn.replace(".gz","")
    set_submit_busco_job_script(GENOME_ID = unzipped_gf,
                                SPECIES_NAME = sps_name,
                                LINEAGE_DB = lineage_db,
                                GENE_PREDICTOR='augustus',
                                NUM_PROCESSORS_PERNODE = "10",
                                WALL_TIME="24",
                                SUBMIT_JOB=submit_job)
    

def generate_pending_run_idxs(consolidated_eukaryotic_genomes_df, dir_main = "Results_Augustus"):
    """Generate a list of idxs in metadata df (consolidated_eukaryotic_genomes_df) that are still pending for
       processing using BUSCO. This list is dumped to pending_run_idxs.csv"""

    busco_out_dirs = get_all_busco_out_dirs(dir_main = dir_main)

    covered_idxs = []
    #Run over each BUSCO folder
    for sps_name in busco_out_dirs:
        #Fetch str to extract patt for querying species dir name list
        str2search = sps_name.replace("_", " ").replace("[", "").replace("]", "")

        fn_patts = np.array(re.split("\W+|_|\.\W+",str2search))

        search_fn = [fn_checker(fn_patts,target_str) for target_str in consolidated_eukaryotic_genomes_df["SpeciesName"]]

        get_entry_in_df = consolidated_eukaryotic_genomes_df[search_fn]

        if (get_entry_in_df.shape[0]>0):
            covered_idxs.append(get_entry_in_df.index.values[0])    

    #Get IDs in consolidated_eukaryotic_genomes_df that are still pending for running to completion
    pridxs = sorted(set(consolidated_eukaryotic_genomes_df.index.values).difference(set(covered_idxs)))
    pending_run_idxs = pd.Series(pridxs)
    pending_run_idxs.to_csv("pending_run_idxs.csv",index=False) 

    
def deploy_subset_busco_runs(consolidated_eukaryotic_genomes_df, dir_main = "Results_Augustus"):
    """Function to submit to an HPC system a subset of an entire list of species/genomes to be processed via BUSCO.
       The starting point is a DF with metadata (consolidated_eukaryotic_genomes_df) associated to 
       individual species. """

    #Clean already generated dirs and remove the faulty ones!
    #Func imported from module: keep_busco_outdirs_lightweight
    #Get species dirs generated by BUSCO runs
    busco_out_dirs = get_all_busco_out_dirs(dir_main = dir_main)
    if(len(busco_out_dirs)>0):

        for sps_name in busco_out_dirs:
            print("Cleaning dir: {}".format(sps_name))
            clean_busco_run_outdirs(sps_name, dir_main = dir_main)
            print("\n")    
        print("###############################################")
        #Get total number of dirs after cleaning
        busco_out_dirs = get_all_busco_out_dirs(dir_main = dir_main)    
        n = str(len(busco_out_dirs))
        print("Total number of BUSCO output directories kept for running Augustus with hint files is: {} ".format(n))

        #Get the idxs list for the pending runs
        generate_pending_run_idxs(consolidated_eukaryotic_genomes_df, dir_main = dir_main)
        #Import list
        pending_run_idxs = pd.read_csv("pending_run_idxs.csv").values.flatten()
        print("Total number of BUSCO pending: {} ".format(len(pending_run_idxs)))
        
        return pending_run_idxs
    
    else:
        return []
    

def main():
    
    """Main function to execute BUSCO on a bunch of genomes, the feats of which are encoded in a metadata df"""

    parser = argparse.ArgumentParser(description = "Wrapper to run BUSCO on a set of genomes, not only to identify BUSCO genes, but mainly for extracting hints for running Augustus and refine the ab initio structural gene prediction in a post-BUSCO run.\n",
                                     formatter_class = argparse.RawTextHelpFormatter,
                                     prog = "euk_genome_annotate.py")    

    parser.add_argument("--dry_run",
                        help = "A dry-run of this program yields the number of genomes that have been already processed through previous BUSCO runs.\n[OPTIONAL]",
                        metavar = "<bool: True/False>",
                        type = bool,
                        action = 'store',
                        nargs = '?',
                        default = True,
                        required = False)     

    parser.add_argument("--run_mode",
                        help = "Set whether to run the pipeline in default (by indicating a list of idxs) or fallback mode (by letting the program randomly choose a set of idxs pending for running).\n[REQUIRED!]",
                        type = str,
                        action = 'store',
                        nargs = '?',
                        default = "default",
                        required = False)    
    
    parser.add_argument("--output_dir",
                        help = "Directory where BUSCO run should put all output files.\n[REQUIRED!]",
                        type = str,
                        action = 'store',
                        nargs = '?',
                        default = "./Results_Augustus",
                        required = False)        
    
    parser.add_argument("--md_fn",
                        help = "Metadata file name.\n[REQUIRED!]",
                        metavar = "<file.csv>",
                        type = str,
                        action = 'store',
                        nargs = '?',
                        default = "NCBI_EukaryoticGenomicResources.csv",
                        required = False)        

    parser.add_argument("--md_path",
                        help = "Path to metadata file.\n[REQUIRED!]",
                        type = str,
                        action = 'store',
                        nargs = '?',
                        default = ".",
                        required = False)     
    
    parser.add_argument("--idxl_fn",
                        help = "File name of list containing a set of idx in the metadata file (set via --md_fn). If this arg is not provided the program will then run in fallback mode by checking which idx in metadata are pending for running BUSCO.\n[OPTIONAL]",
                        metavar = "<file.csv>",
                        type = str,
                        action = 'store',
                        nargs = '?',
                        default = "pending_run_idxs_subset1.csv",
                        required = False)  

    parser.add_argument("--idxl_path",
                        help = "Path to list containing a set of idx in the metadata file.\n[OPTIONAL]",
                        metavar="<./>",
                        type = str,
                        action = 'store',
                        nargs = '?',
                        default = ".",
                        required = False)     
    
    parser.add_argument("--max_num_jobs",
                        help = "Set number of jobs (distinct genomes) to submit by default when program is set to run in fallback mode\n[OPTIONAL]",
                        type = int,
                        action = 'store',
                        nargs = '?',
                        default = 50,
                        required = False)        
    

    # Execute parse_args()
    args = parser.parse_args()   
    
    #These fields are mandatory in the metadata file (DF)!
    FIELDS = ["Species","FtpPathGenBank","BUSCO_Lineage"]  
    
    #BUSCO output dir
    DIR_MAIN = args.output_dir

    #Print message when program is run in dry-run mode
    if(args.dry_run is None):
        subdirs_list = glob.glob(os.path.join(DIR_MAIN,"*"))    
        if(len(subdirs_list)>0):
            busco_out_dirs = [os.path.basename(f) for f in subdirs_list if os.path.isdir(f)]
            print("DRY-RUN: number of already processed files in previous BUSCO runs is {}".format(len(busco_out_dirs)))
        else:
            print("DRY-RUN: BUSCO has not been previously run on sample genomes!")
            
    #Else, run program according to  command line args provided!
    else:
        #Metadata df
        md_fn = os.path.join(args.md_path,args.md_fn)
        consolidated_eukaryotic_genomes_df = pd.read_csv(md_fn)    
        
        try:
        
            if(args.run_mode == "default"):
                idxl_fn = os.path.join(args.idxl_path,args.idxl_fn)
                set_idxs = pd.read_csv(idxl_fn).values.flatten();
                print("Running in default mode!\n")
                print(set_idxs)
                
            elif(args.run_mode == "fallback"):
                MAX_NUM_JOBS = int(args.max_num_jobs)    
                #Fallback running mode: clean already generated dirs and remove the faulty ones!
                pending_run_idxs = deploy_subset_busco_runs(consolidated_eukaryotic_genomes_df, dir_main = DIR_MAIN)    
                #Set a limit of jobs to be submitted to prevent jobs from crashing due to space disk quota
                if(len(pending_run_idxs)<=MAX_NUM_JOBS):
                    set_idxs = pending_run_idxs
                else:
                    set_idxs = pending_run_idxs[:MAX_NUM_JOBS]
                print("Running in fallback mode!\n")
                print(set_idxs)

            #Submit jobs!
            for i in set_idxs:
                sps_name, ftp_path, lineage_db = consolidated_eukaryotic_genomes_df[FIELDS].iloc[i].values
                print("##########################################")
                print("Preparing BUSCO run for: {}".format(sps_name))

                submit_busco_job(sps_name,ftp_path,lineage_db,download_genome=True,submit_job=True) 
                print("Job submitted!\n")
                #Delay fetching genomes to ensure integrity of file
                time.sleep(5)

            print("ALL JOBS SUBMITTED!")
        except:
            pass    

if __name__ == "__main__":
    
    main()

