#!/usr/bin/env python

"""Code developed to run Augustus (a gene prediction tool for eukaryotic genomic sequences: https://bioinf.uni-greifswald.de/augustus/) 
   Importantly, Augustus is run using a hints file obtained by previously running the BUSCO (https://busco.ezlab.org/) software on a given genome sequence.

   Author: jayson.gutierrez@vliz.be

"""

from subprocess import Popen, call, STDOUT, PIPE
import time
import sys
import pandas as pd
import os
import glob
import shutil
import numpy as np
from collections import OrderedDict
#deepest is a tool to determine the maximum depth and path length within the current (or a specified) directory tree. To install do: pip install deepest
import deepest
import argparse

from keep_busco_outdirs_lightweight import *

#Template for the augustus job script
shell_script_template = """#!/bin/bash

###Running Augustus using hints obtained from a previous BUSCO run

###To run this script first activate the appropriate conda env: 
###e.g. source activate /data/gent/vo/001/gvo00125/vsc43582/Software/annotation

###QUERY: path to genome assembly, it should be a *.fna file
QUERY="%(QUERY_SEQ)s" 

###OUTPUT_PATH: path to output Augustus results, e.g. "./Results_Augustus"
OUTPUT_PATH="%(AUGUSTUS_OUTPUT_PATH)s" 

###SPECIES: name of the species the genome assembly belongs to, e.g. "Aspergillus_puulaauensis"
SPECIES="%(SPECIES_NAME)s" 

echo "Running Augustus to conduct structural annotation for the genome of species ${SPECIES}"

AUGUS_OUTPUT="$OUTPUT_PATH/Predicted_ORFs_${SPECIES}.gff"
### Check for dir, if not found, create it
if [[ ! -d "$OUTPUT_PATH" ]]; then
  mkdir -p $OUTPUT_PATH
fi
echo "Augustus output is directed to:"
echo ${AUGUS_OUTPUT}
#Run Augustus
augustus --species=${SPECIES} ${QUERY} > ${AUGUS_OUTPUT} 

#Put predicted gene sequences into a fasta files
getAnnoFasta.pl $AUGUS_OUTPUT

echo "Augustus output is on: ${AUGUS_OUTPUT_FASTA}"
"""

#Template for submitting to the HPC the augustus job
submit_job_script_template = '''#!/bin/bash -l

#PBS -l nodes=1:ppn=%(NUM_PROCESSORS_PERNODE)s 
#PBS -l walltime=%(WALL_TIME)s:00:00

cd $PBS_O_WORKDIR

module --force purge
module load cluster/swalot

source activate /data/gent/vo/001/gvo00125/vsc43582/Software/annotation

bash %(shell_script_name)s'''


def fetch_seq_data_from_NCBI(ftp_path, dest_dir = "./Genomes/genbank/", fetch_genome = True):
    """Function to fetch a genome/proteome from NCBI, or any other public repository"""
    #Set path to dir to store genome
    genome_dir = dest_dir + os.path.basename(os.path.dirname(ftp_path))
    #unzip file
    gunzip_fn = os.path.join(genome_dir,os.path.basename(ftp_path))

    if(fetch_genome):
        #Create dir if genome dir doesn't exist
        if(~os.path.exists(genome_dir)):
            runShell("mkdir -p {}".format(genome_dir))

        #Keep trying downloading genome file until valid format is fetched
        print("Fetching genome/proteome from: {}".format(ftp_path))
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
            
    #Path to unzipped genome assembly, which can be passed in to BUSCO/Augustus or any other software
    unzipped_gf = gunzip_fn.replace(".gz","")
    
    return os.path.join(genome_dir,unzipped_gf)
    

def species_name2NCBI_ftpath(dir_main = "Results_Augustus", path2meta_df = "NCBI_EukaryoticGenomicResources.csv"):
    """Generate a dict mapping species dir name to their corresponding NCBI ftp path to fetch genome assembly"""
    
    #Load in metadata, a df containing fileds like SpeciesName & FtpPathGenBank
    consolidated_eukaryotic_genomes_df = pd.read_csv(path2meta_df)
    
    #List all busco folders generated in previous runs
    busco_out_dirs = get_all_busco_out_dirs(dir_main = dir_main)
    
    try:
        #Let's map species folder name to NCBI ftpath, and cast data into a dict
        sps_fns2ftpath_dict = OrderedDict()
        #Run over each BUSCO folder
        for sps_name in busco_out_dirs:
            #Fetch str to extract patt for querying species dir name list
            str2search = sps_name.replace("_", " ").replace("[", "").replace("]", "")

            fn_patts = np.array(re.split("\W+|_|\.\W+",str2search))

            search_fn = [fn_checker(fn_patts,target_str) for target_str in consolidated_eukaryotic_genomes_df["SpeciesName"]]

            get_entry_in_df = consolidated_eukaryotic_genomes_df[search_fn]

            if (get_entry_in_df.shape[0]>0):
                sps_fns2ftpath_dict[sps_name] = get_entry_in_df["FtpPathGenBank"].values[0] 

        return sps_fns2ftpath_dict
    
    except:
        print("No previous BUSCO runs have been deployed!")


def rename_move_all_augustus_hints_files(sps_name = "Calanus_finmarchicus",
                                         dir_main = "Results_Augustus",
                                         path2meta_df = "NCBI_EukaryoticGenomicResources.csv",
                                         assmbl_dest_dir = "./Genomes/genbank/",
                                         augustus_config_path = "/data/gent/vo/001/gvo00125/vsc43582/Software/annotation/config/",
                                         fetch_genome = False):
    """Moving all the Augustus hints files generated during a BUSCO run to $AUGUSTUS_CONFIG_PATH, 
       and renaming those appropriately. 
       NOTE: Use this newly created hints folder to further refine the gene prediction task using Augustus, following 
       the initial BUSCO run.
    """
    
    #Get dict mapping species dir name to NCBI ftpath to assembly
    sps_fns2ftpath_dict = species_name2NCBI_ftpath(dir_main = dir_main, path2meta_df = path2meta_df)
   
    try:
        #Fetch genome assembly from NCBI ftp
        ftp_path = sps_fns2ftpath_dict[sps_name]

        unzipped_gf = fetch_seq_data_from_NCBI(ftp_path, dest_dir = assmbl_dest_dir, fetch_genome = False)

        #Check if retraining_parameters subdir exists for the current genome
        aug_retrain_params_dir = glob.glob(os.path.join(dir_main,sps_name,"*_odb10","augustus_output","retraining_parameters"))

        #Set new dir to store hints where Augustus can find them!
        set_new_dir2hints = os.path.join(augustus_config_path,"species",sps_name)

        if(len(aug_retrain_params_dir)>0):
            #Get deepest, nested subdir
            deepest_augustus_path = deepest.get_depth(aug_retrain_params_dir[0])[0]
            #Get the prefix for the hints
            hints_prefix = os.path.split(deepest_augustus_path)[-1]
            #List all hint files and their location
            original_path2all_hints = glob.glob(os.path.join(deepest_augustus_path,"*"))
            #Get only original hint file names
            original_hints_fns = [os.path.split(f)[-1] for f in original_path2all_hints] 

            #Make dir if it doesn't exist for hints where Augustus can find them!
            if(os.path.exists(set_new_dir2hints)==False):
                runShell("mkdir -p {}".format(set_new_dir2hints))

            #Create a dict mapping original hint file name to new name: create this mapping just in case
            ori2new_hints_fns = OrderedDict()
            for f in original_hints_fns:
                if(hints_prefix in f):
                    orifpath = os.path.join(deepest_augustus_path,f)
                    ori2new_hints_fns[orifpath] = os.path.join(set_new_dir2hints, f.replace(hints_prefix,sps_name))

            #Now use dict ori2new_hints_fns to copy original hint files to dir where Augustus can find them!
            if(len(ori2new_hints_fns)>0):
                
                for (src, dst) in ori2new_hints_fns.items():
                    #Copy original to $AUGUSTUS_CONFIG_PATH anyway!
                    runShell("cp {} {}".format(src, os.path.dirname(dst)))
                    #No prefix BUSCO appended to new file names
                    #Remove dst first anyway!
                    runShell("rm {}".format(dst))
                    runShell("cp {} {}".format(src, dst))
                    #shutil.copy(src, dst) 
                    #Naming with prefix BUSCO
                    dst = os.path.join(os.path.dirname(dst), "BUSCO_" + os.path.basename(dst))  
                    #Remove dst first anyway!
                    runShell("rm {}".format(dst))
                    runShell("cp {} {}".format(src, dst))
                    #shutil.copy(src, dst)                       

                print("Copying and moving files to AUGUSTUS_CONFIG_PATH")
                return set_new_dir2hints, unzipped_gf

        else:
            print("Copying and moving files operations failed!")

    except:
        print("It seems like no previous BUSCO runs have been deployed yet!")
        
        
def create_augustus_job_script(sps_name,
                               dir_main = "/data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/Annotation/Scripts/Busco_runs/Results_Augustus",
                               augustus_out_dir = "/data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/Annotation/Scripts/Busco_runs/Results_Augustus_Post_BUSCO",
                               path2meta_df = "/data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/Annotation/Scripts/Busco_runs/NCBI_EukaryoticGenomicResources.csv",
                               assmbl_dest_dir = "/data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/Annotation/Scripts/Busco_runs/Genomes/genbank/",
                               augustus_config_path = "/data/gent/vo/001/gvo00125/vsc43582/Software/annotation/config",
                               num_proc_pernode = "1",
                               walltime="48",
                               create_submit_job_script=False,
                               fetch_genome = True):
    """Function that creates a script to invoke augustus and run it on a given genome (associated to sps_name, which is
       linked to the input genome assembly via a metadata df indicated by path2meta_df). 
       In addition to creating the script for invoking augustus, this function can optionally create a script to
       submit the augustus job on the HPC (currently the UGent HPC system) and send it to the queue
    
    params:
        - sps_name: name of the species (e.g. Calanus_finmarchicus) whose genome assembly is to be processed via augustus
    """
    
    #Get hint files obtained via previous BUSCO run, move to $AUGUSTUS_CONFIG_PATH and rename those hint files to match
    #the name of the corresponding species/organism
    set_new_dir2hints, query_seq = rename_move_all_augustus_hints_files(sps_name,
                                                                        dir_main,
                                                                        path2meta_df,
                                                                        assmbl_dest_dir,
                                                                        augustus_config_path,
                                                                        fetch_genome)

    #Formatting shell script used to invoke Augustus: shell_script_template defined above!
    shell_text = shell_script_template % {'QUERY_SEQ':query_seq,
                                          'AUGUSTUS_OUTPUT_PATH':augustus_out_dir,
                                          'SPECIES_NAME':sps_name}

    #Set path to shell script used to invoke Augustus
    shell_script_name = os.path.join(os.path.dirname(dir_main), "run_augustus_{}.sh".format(sps_name.lower()))
    #Save shell script to file
    with open(shell_script_name,"w") as ss:
        ss.write(shell_text)

    #Submission job script to HPC
    submit_job_script_name = shell_script_name.replace("run_","submit_")
    #Set template for job submission to HPC, fill it in and save to file with a given label/name
    if(create_submit_job_script):

        #Fill in job submission script with respective fields for a given species/organism: submit_job_script_template
        #defined above!
        submit_job_script = submit_job_script_template % {'NUM_PROCESSORS_PERNODE':num_proc_pernode,\
                                                           'WALL_TIME':walltime,
                                                           'shell_script_name':shell_script_name}

        #Save job submission script to file
        with open(submit_job_script_name,"w") as ss:
            ss.write(submit_job_script)

        print("Submitting job via: {}".format(submit_job_script_name))
        #Submit job to the HPC cluster
        runShell("qsub {}".format(submit_job_script_name))

    return shell_script_name,submit_job_script_name 

        
def run_augustus_local_test(sps_name):
    """Running a local test (rbi@vliz) with augustus based on hint files obtained in previous BUSCO runs"""
    
    dir_main = "/data/home/VLIZ2000/jayson.gutierrez/meta_omics_data/annotation/Busco_runs/Results_Augustus"
    augustus_out_dir = "/data/home/VLIZ2000/jayson.gutierrez/meta_omics_data/annotation/Busco_runs/Results_Augustus_Post_BUSCO"
    path2meta_df = "/data/home/VLIZ2000/jayson.gutierrez/meta_omics_data/Struo/PlanktonSeqDB_draft1/Files/NCBI_EukaryoticGenomicResources.csv"
    assmbl_dest_dir = "/data/home/VLIZ2000/jayson.gutierrez/meta_omics_data/annotation/Busco_runs/Genomes/genbank/"
    augustus_config_path = "/home/VLIZ2000/jayson.gutierrez/anaconda3/envs/bioinformatics/config"
    num_proc_pernode = "1"
    walltime="48"
    create_submit_job_script=False
    fetch_genome = False

    shell_script_name,_ = create_augustus_job_script(sps_name,
                                                     dir_main,
                                                     augustus_out_dir,
                                                     path2meta_df,
                                                     assmbl_dest_dir,
                                                     augustus_config_path,
                                                     num_proc_pernode,
                                                     walltime,create_submit_job_script,
                                                     fetch_genome)            

    print("Running augustus on the genome of species: {}, using hints obtained via a previous BUSCO run".format(sps_name))
    print("###############################################")   
    print(shell_script_name)
    out, err = runShell("bash {}".format(shell_script_name))
    
    if(len(err)>0):
        print("Script: {} exited with error: {}\n".format(shell_script_name,err))


def submit_augustus_job2hpc(sps_name = "Acanthamoeba_castellanii",
                            dir_main = "/data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/Annotation/Scripts/Busco_runs/Results_Augustus",
                            augustus_out_dir = "/data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/Annotation/Scripts/Busco_runs/Results_Augustus_Post_BUSCO",
                            path2meta_df = "/data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/Annotation/Scripts/Busco_runs/BUSCO_Processed_PlanktEukary_Genomes.csv",
                            assmbl_dest_dir = "/data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/Annotation/Scripts/Busco_runs/Genomes/genbank/",
                            augustus_config_path = "/data/gent/vo/001/gvo00125/vsc43582/Software/annotation/config",
                            num_proc_pernode = "1", walltime="48"):
    
    """Function to submit a job to the HPC system to run Augustus on a previously BUSCO-processed genome assembly,
       by which a set of hint files were obtained"""

    #Read in metadata DF
    consolidated_eukaryotic_genomes_df_final = pd.read_csv(path2meta_df)

    #Get list of all the organisms whose genomes were previously processed via BUSCO
    all_busco_preprocessed_dirs = consolidated_eukaryotic_genomes_df_final["BUSCOutputDir"].apply(lambda s: s.split("/")[1]).values
    #Get path to subdir containing the genome assembly
    genomes_subdirs = consolidated_eukaryotic_genomes_df_final["FtpPathGenBank"].apply(lambda s: "/".join(s.split("/")[-2:]).replace("fna.gz","fna"))
    #Cast data into a dict
    sps2ftpath_subdir_dict = dict(zip(all_busco_preprocessed_dirs, genomes_subdirs))

    #Path to genome assembly
    query_seq = os.path.join(assmbl_dest_dir,sps2ftpath_subdir_dict[sps_name])

    #Formatting shell script used to invoke Augustus: shell_script_template defined above!
    shell_text = shell_script_template % {'QUERY_SEQ':query_seq,
                                          'AUGUSTUS_OUTPUT_PATH':augustus_out_dir,
                                          'SPECIES_NAME':sps_name}

    #Set path to shell script used to invoke Augustus
    shell_script_name = os.path.join(os.path.dirname(dir_main), "run_augustus_{}.sh".format(sps_name.lower()))

    #Save shell script to file
    with open(shell_script_name,"w") as ss:
        ss.write(shell_text)

    #Submission job script to HPC
    submit_job_script_name = shell_script_name.replace("run_","submit_")

    #Fill in job submission script with respective fields for a given species/organism: submit_job_script_template
    #defined above!
    submit_job_script = submit_job_script_template % {'NUM_PROCESSORS_PERNODE':num_proc_pernode,\
                                                      'WALL_TIME':walltime,
                                                      'shell_script_name':shell_script_name}

    #Save job submission script to file
    with open(submit_job_script_name,"w") as ss:
        ss.write(submit_job_script)

    print("Submitting job via: {}".format(submit_job_script_name))
    #Submit job to the HPC cluster
    runShell("qsub {}".format(submit_job_script_name))    


if __name__ == "__main__":

    sps_name = sys.argv[1] #e.g. Acanthamoeba_castellanii
    dir_main = "/data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/Annotation/Scripts/Busco_runs/Results_Augustus"
    augustus_out_dir = "/data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/Annotation/Scripts/Busco_runs/Results_Augustus_Post_BUSCO"
    path2meta_df = "/data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/Annotation/Scripts/Busco_runs/BUSCO_Processed_PlanktEukary_Genomes.csv"
    assmbl_dest_dir = "/data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/Annotation/Scripts/Busco_runs/Genomes/genbank/"
    augustus_config_path = "/data/gent/vo/001/gvo00125/vsc43582/Software/annotation/config"
    num_proc_pernode = "1"
    walltime="48"

    submit_augustus_job2hpc(sps_name, dir_main, augustus_out_dir, path2meta_df, assmbl_dest_dir, augustus_config_path, num_proc_pernode, walltime)
 
