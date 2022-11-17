#!/usr/bin/env python

'''Script to scoop out from the Eukzoo DB (https://zenodo.org/record/1476236#.YfpQwv7MJPY) the species for which no genome assembly is available via NCBI. 
   For these species, the corresponding assemblied transcriptome could be annotated against UniRef as well. This will increase massively the diversity of 
   our custom DB! However, given that those seqs have already been annotated in the EukZoo DB, we can simply filter those sequences and format the header 
   appropriately! 
   NOTE: here I filter only those species for which an NCBI taxID is available, down to the species level, because we want high taxonomic resolution.

Author: jayson.gutierrez@vliz.be
'''

from subprocess import Popen, call, STDOUT, PIPE
import os
import pandas as pd
import numpy as np

import Bio.SeqIO as bioseqio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import Entrez

import urllib
from collections import OrderedDict 

from utilities import *

def filter_eukzoo_db(in_fasta, filtered_MMETSP, out_fasta):
    #List of MMETSP ID, the annotated protein sequences of which will be pulled out
    source_ids2filter = filtered_MMETSP["Source_ID"].values

    with open(in_fasta,"r") as handle:

        for (i,record) in enumerate(bioseqio.parse(handle, "fasta")):
            #Fetch source MMETSP ID
            source_id = record.description.split("_")[0]
            #If current record in list then relabeled and dump into new .faa file
            if(source_id in source_ids2filter):        
                #Changing entry header's label
                tax_id = filtered_MMETSP.query("Source_ID=='{}'".format(source_id))["ncbi_species_taxid"].values[0]
                full_ncbi_tax_tree = filtered_MMETSP.query("Source_ID=='{}'".format(source_id))["full_ncbi_tax_tree"].values[0]
                ncbi_taxonomy = ";".join(full_ncbi_tax_tree.split(";")[-2:])

                #Format the entry in a way that can be modified properly when the seq has been re-annotated against latest UniRef version
                record.description = "{}:UniRefID|{}|{}__taxID{}".format(source_id,len(record.seq),ncbi_taxonomy,tax_id)
                record.id = record.description
                record.name = record.description        

                ##Save as .fasta        
                with open(out_fasta, "a") as output_handle:
                    bioseqio.write(record, output_handle, "fasta")         

if __name__ == "__main__":
    
    ###Set path to EukZoo DB
    fn = '../marine_databases/TARA_Oceans/EukZoo/{}'

    #Input/Output .faa
    in_fasta = fn.format("EukZoo_v_0.2.faa")
    out_fasta = fn.format("EukZoo_filtered_relabeled.faa")

    #Load in initial DB. It should be noticed that struo_metadata_df contains info on mostly bacterioplankton species, while
    #leftout_struo_metadata_df contains info on eukaryotic species
    db_fn = "/data/home/VLIZ2000/jayson.gutierrez/meta_omics_data/Struo/PlanktonSeqDB_draft1/"
    #Metadata file
    filtered_MMETSP = pd.read_csv(os.path.join(db_fn,'Files/filtered_MMETSP_transcriptomes_assemblies_metadata.txt'),sep="\t",index_col="Unnamed: 0")

    filter_eukzoo_db(in_fasta, filtered_MMETSP, out_fasta)