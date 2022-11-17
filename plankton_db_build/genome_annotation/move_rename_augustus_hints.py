#!/usr/bin/env python

"""Script containing functions used to move and rename Augustus hints files generated during a BUSCO run 

   Author: jayson.gutierrez@vliz.be

"""

import os
import glob
import shutil
import sys
import re

def rename_augustus_hints_file(fn,fmt,AUGUSTUS_CONFIG_PATH,OUTPUT):
    """Moving and renaming a given hint file"""
    
    #Format the name of the folder where Augustus seeks for hints
    new_folder = os.path.join(AUGUSTUS_CONFIG_PATH,re.sub('[-,.,:,;]', '_', OUTPUT.lower()))

    if not os.path.exists(new_folder):
        os.makedirs(new_folder)  

    #Set path to file
    mv_file = shutil.copy(fn,new_folder)

    #To save hint files with different formatting
    if(fmt=="modified"):
        #Formatted file to make it compatible with the Augustus nomenclature and rename destination
        prefix4files = os.path.basename(mv_file).replace("BUSCO_","").lower()
        prefix4files = re.sub('[-,:,;]', '_', prefix4files)
        rename_file = os.path.join(new_folder,prefix4files)
    
    else:
        rename_file = os.path.join(new_folder,os.path.basename(mv_file))
    
    if not os.path.exists(rename_file):
        shutil.move(mv_file, rename_file) 
        
        
def rename_move_all_augustus_hints_files(OUTPUT_PATH,OUTPUT,LINEAGE_DATASET,AUGUSTUS_CONFIG_PATH):
    """Moving all the Augustus hints files generated during a BUSCO run to $AUGUSTUS_CONFIG_PATH, 
       and renaming those appropriately. 
       NOTE: Use this newly created hints folder to further refine the gene prediction task using Augustus, following 
       the initial BUSCO run"""
    
    #Path to Augustus hints file generated during the BUSCO run
    augustus_hints_path="{0}/{1}/run_{2}/augustus_output/retraining_parameters/BUSCO_{1}".format(OUTPUT_PATH,OUTPUT,LINEAGE_DATASET)
    
    #Folder where Augustus seeks for hints
    new_folder = os.path.join(AUGUSTUS_CONFIG_PATH,re.sub('[-,.,:,;]', '_', OUTPUT.lower()))
    
    #Always remove new folder before processing and moving original files to this new location
#     if os.path.exists(new_folder):
#         shutil.rmtree(new_folder)
    
    #Get file names in the originally generated folder
    augustus_hints_from_busco_run = glob.glob(os.path.join(augustus_hints_path,"*"))
    
    #Run over each original file and move to new destination (AUGUSTUS_CONFIG_PATH)
    for f in augustus_hints_from_busco_run:
        rename_augustus_hints_file(f,"original",AUGUSTUS_CONFIG_PATH,OUTPUT)

    for f in augustus_hints_from_busco_run:
        rename_augustus_hints_file(f,"modified",AUGUSTUS_CONFIG_PATH,OUTPUT)

if __name__ == "__main__":

    ###Example set of args to run this script
    # OUTPUT_PATH="./Results_Augustus"
    # OUTPUT="Pseudo-nitzschia_multistriata"
    # LINEAGE_DATASET="stramenopiles_odb10"
    # AUGUSTUS_CONFIG_PATH="../config/species"    
    OUTPUT_PATH,OUTPUT,LINEAGE_DATASET,AUGUSTUS_CONFIG_PATH = sys.argv[1:]
    
    rename_move_all_augustus_hints_files(OUTPUT_PATH,OUTPUT,LINEAGE_DATASET,AUGUSTUS_CONFIG_PATH)
    
    
    
    
    
