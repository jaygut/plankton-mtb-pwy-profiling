"""Functions to run humann3 on a series of samples: submitting jobs to the HPC
  
   Humann3 can be deployed on an HPC system using either:
    - the default DBs (chocophlan + UniRef90)
    - a custom DB assembled from a bunch of genome assemblies fetched from e.g. NCBI  
Author: jayson.gutierrez@vliz.be
"""

from subprocess import Popen, call, STDOUT, PIPE
import time
import glob
import os
import sys

from utilities import run_shell


#Use template to create a job script for the VSC HPC system and run humann3 on a given sample
job_script_template = '''#!/bin/bash -l

#PBS -l nodes=1:ppn=%(num_processors_pernode)s
#PBS -l mem=%(memory)sgb
#PBS -l walltime=%(wall_time)s:00:00

###Running humann3 on a particular sample using only translated search via diamond on uniref90 2020_1
###This uniref version contains 130661074 representative sequences

WORK_DIR=%(work_dir)s

INPUT=${WORK_DIR}%(input_file)s

OUTPUT=${WORK_DIR}%(output_file)s

#Folder containing both full UniRef90 DB and our custom DB created using Struo2 + BUSCO + Augustus
PROT_DB=${WORK_DIR}/data/databases/%(database_name)s

cd ${WORK_DIR}
module swap cluster/swalot
source activate biobakery

RUNNING_MODE=%(running_mode)s

NUM_THREADS=%(num_threads)s

if  [ "$RUNNING_MODE" = "default_db" ];
then
  humann3 --input ${INPUT} --output ${OUTPUT} --threads ${NUM_THREADS}
else
  humann3 --input ${INPUT} --output ${OUTPUT} --protein-database ${PROT_DB} --bypass-nucleotide-search --threads ${NUM_THREADS}
fi

'''

def create_submit_humann_job_script(num_processors_pernode = '1', memory = '60', wall_time = '72',
                                    work_dir = '/data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/HUMAnN3/',
                                    input_file = 'data/reads/Elferink_etal_2020_mtx_data/ERR3980595_Qual_Filt2.fastq',
                                    output_file = 'results/metatransc_dataset_PRJEB37134_0',
                                    submit_job_script_name = 'run_humann3_uniref_prot_db_only_sample_ERR3980595.sh',
                                    database_name = 'uniref',
                                    running_mode = 'default',
                                    num_threads = '24',
                                    submit_job = True):
   
    '''Function to run humann3 on a given input sample. This generates a script to submit a job to an HPC system. NOTE: set function's args appropriately.'''
    
    #Create script using the template above
    submit_job_script = job_script_template % {'num_processors_pernode':num_processors_pernode,
                                               'memory':memory,
                                               'wall_time':wall_time,
                                               'work_dir':work_dir,
                                               'input_file':input_file,
                                               'output_file':output_file,
                                               'database_name':database_name,
                                               'running_mode':running_mode,
                                               'num_threads':num_threads}
    
    #Save to a .sh file
    with open(submit_job_script_name,"w") as ss:
        ss.write(submit_job_script)
    
    #Check if job is to be submitted to the HPC system
    if(submit_job):
        time.sleep(1)    
        run_shell("qsub {}".format(submit_job_script_name))


def run_hummann_on_samples(samples_list, running_mode = 'custom_db', num_threads = '24', 
                           work_dir = '/data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/HUMAnN3/',
                           submit_job = True):
    
    '''Automatize job submission: running humann3 on a bunch of samples, i.e. QC fastq files. NOTE: set function's args appropriately'''

    if (running_mode=='default_db'):
        #Run humann3 on a set of samples in both nucleotide and translated mode using the default Chocophlan and UniRef90 DB (diamond-formatted)
        outf = 'results/metatransc_dataset_PRJEB37134_0'
        for sfn in samples_list:
            sid = os.path.basename(sfn).split('_')[0]
            in_file = 'data/reads/Elferink_etal_2020_mtx_data/{}_Qual_Filt.fastq'.format(sid)
            job_sn = 'run_humann3_prof_with_default_dbs_sample_{}.sh'.format(sid)
            create_submit_humann_job_script(num_processors_pernode = '1', memory = '60', wall_time = '72',
                                            work_dir = work_dir,
                                            input_file = in_file,
                                            output_file = outf,
                                            submit_job_script_name = job_sn,
                                            database_name = 'uniref',
                                            running_mode = running_mode,
                                            num_threads = num_threads,
                                            submit_job = submit_job)    
    elif (running_mode=='custom_db'):
        #Run humann3 on a set of samples in translated mode only (with nucleotide bypassed) using a custom DB annotated against UniRef90 (diamond-formatted)
        outf = 'results/metatransc_dataset_PRJEB37134_1'
        for sfn in samples_list:
            sid = os.path.basename(sfn).split('_')[0]
            in_file = 'results/metatransc_dataset_PRJEB37134_0/{0}_Qual_Filt_humann_temp/{0}_Qual_Filt_diamond_unaligned.fa'.format(sid)
            job_sn = 'run_humann3_prof_with_custom_prot_db_sample_{}.sh'.format(sid)
            create_submit_humann_job_script(num_processors_pernode = '1', memory = '60', wall_time = '72',
                                            work_dir = work_dir,
                                            input_file = in_file,
                                            output_file = outf,
                                            submit_job_script_name = job_sn,
                                            database_name = 'custom_db',
                                            running_mode = running_mode,
                                            num_threads = num_threads,                                          
                                            submit_job = submit_job)
    else:
        print("Running mode: {} is not valid! Choose either default_db or custom_db".format(str(running_mode)))

if __name__ == "__main__":
    #Set to run either with humann3 default DBs (chocophlan + UniRef90) or with custom DB
    running_mode = sys.argv[1] #Set to: 'default_db' | 'custom_db'
    #Set data_dir appropriately
    data_dir = sys.argv[1] #In my case: '/data/gent/vo/001/gvo00125/vsc43582/Bioinformatics/HUMAnN3/data/reads/Elferink_etal_2020_mtx_data'
    #Make sure .fastq files are stored in data_dir
    samples_list = glob.glob(os.path.join(data_dir,'*.fastq'))
    #Create sample-specific job scripts and submit to VSC-HPC system 
    run_hummann_on_samples(samples_list, running_mode = running_mode, submit_job = True)


