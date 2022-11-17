#!/home/VLIZ2000/jayson.gutierrez/anaconda3/envs/bioinformatics/bin/python

"""Processing the SMAGs DB published by Delmont et.al. 2020 (https://www.biorxiv.org/content/10.1101/2020.10.15.341214v2):  

Excerpt from abstract:
We used 280 billion Tara Oceans metagenomic reads from polar, temperate, and tropical sunlit oceans to reconstruct and manually curate more than 700 abundant 
and widespread eukaryotic environmental genomes ranging from 10 Mbp to 1.3 Gbp. This genomic resource covers a wide range of poorly characterized eukaryotic 
lineages that complement long-standing contributions from culture collections while better representing plankton in the upper layer of the oceans.  

The DB is accessible from TARA Oceans repository (https://www.genoscope.cns.fr/tara/).  

The goal here is:
To properly format this protein seq DB in a way that one can use it with humann to profile environmental samples at the metabolic pathway level. To achieve this, 
we have to: 

    1. Scoop out sequences that can be taxonomically associated with an NCBI taxID down to the genus level, or at least the family level (using BioPython here!).
    2. Make sure that a given sequence is sufficiently large, e.g. seq-len > 80 AA.
    3. Create sequence header placeholders for rapid assignment of UniRef90 (2020-2) IDs to each sequence match upon running Diamond.
    
Jayson Gutierrez: 1/2/2021
"""

from subprocess import Popen, call, STDOUT, PIPE
import os
import shutil
import pandas as pd
import json
import glob
import re
import gzip
import sys
import csv
import time
import io
import pathlib

import Bio.SeqIO as bioseqio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC
from Bio import Entrez

from ete3 import NCBITaxa
from taxonomy_ranks import TaxonomyRanks

from utilities import *

def entry_filter(tax_ranks):
    '''Function to filter taxonomic ranks in the SMAGs metadata'''
    gr, fr, orr, cr = tax_ranks
    if( (isinstance(gr,str)) and ('_' not in gr) and ('-' not in gr)):
        return gr
    elif( (isinstance(fr,str)) and ('_' not in fr) and ('-' not in fr)):
        return fr
    elif( (isinstance(orr,str)) and ('_' not in orr) and ('-' not in orr)):
        return orr
    elif( (isinstance(cr,str)) and ('_' not in cr) and ('-' not in cr)):
        return cr
    else:
        return 'NA'

def orgn2taxranks(tax_ranks,best_tax_g2taxlvls):
    '''Utility to map the computed tax ranks corresponding to a given entry in the SMAGs metadata'''
    gr, fr, orr, cr = tax_ranks
    if(gr in best_tax_g2taxlvls.keys()):
        return best_tax_g2taxlvls[gr]
    elif(fr in best_tax_g2taxlvls.keys()):
        return best_tax_g2taxlvls[fr]
    elif(orr in best_tax_g2taxlvls.keys()):
        return best_tax_g2taxlvls[orr]
    elif(cr in best_tax_g2taxlvls.keys()):
        return best_tax_g2taxlvls[cr]
    else:
        return 'NA'

def modify_TOSAG_labels(s):
    '''Utility to modify entries with TOSAG substring'''
    if('TOSAG' in s): 
        return s.split('TOSAG')[0] + 'TOSAG' + s.split('TOSAG')[1].replace("_","-")
    else:
        return s    

def make_fasta_header_templ(orgs_list):
    '''Function that takes in a list of organisms and query the ncbi taxid in order to fetch the 
       taxonomic range and tax id down to the genus level. 
       The function receives a list of organism names (orgs_list) and returns a dict (best_tax_g2taxlvls)
    '''
    best_tax_g2taxlvls = {}
    for g in orgs_list:
        if(g not in best_tax_g2taxlvls.keys()):
            #Format tax ranks for fasta headers
            try:
                rank_taxon = TaxonomyRanks(g)
                rank_taxon.get_lineage_taxids_and_taxanames()
                rank_dict = list(rank_taxon.lineages.values())[0]

                #Fetch tax ranks in a given order for formatting fasta headers according to humann
                gfo_ranks = []
                for rrt in ['genus','species','class','order','family']:
                    txn = rank_dict[rrt][0]
                    if(txn=='NA'):
                        gfo_ranks.append('unclassified')
                    else:
                        gfo_ranks.append(txn)

                #Check for NCBI tax IDs (favoring lower ranks)
                if(rank_dict['genus'][1]!='NA'):
                    ncbi_txid = rank_dict['genus'][1]
                elif(rank_dict['family'][1]!='NA'):
                    ncbi_txid = rank_dict['family'][1]
                elif(rank_dict['order'][1]!='NA'):
                    ncbi_txid = rank_dict['order'][1]
                elif(rank_dict['class'][1]!='NA'):
                    ncbi_txid = rank_dict['class'][1]
                else:
                    ncbi_txid = 'NA'

                #Add ncbi taxid
                gfo_ranks.insert(2,ncbi_txid)
                #Update dict
                best_tax_g2taxlvls[g] =  "UniRef90|SeqLen|g__{}.s__{}|{}|c__{}.o__{}.f__{}".format(*gfo_ranks)    

            except:
                pass
           
    #Then filter out those cases for which taxID are NA
    new_best_tax_g2taxlvls = dict()
    # Iterate over all the items in dictionary and filter items which has even keys
    for (key, value) in best_tax_g2taxlvls.items():
       # Check if key is even then add pair to new dictionary
       if ('|NA|' not in value):
        new_best_tax_g2taxlvls[key] = value
        
    return new_best_tax_g2taxlvls


def reformat_SMAGs_db(input_table,fasta_fid,new_fasta_fid,seq_len_thr=90):
    '''Main function: reformatting SMAGs DB to make seq headers easily compatible with UniRef90 re-annotation via diamond, 
       which can be used downstream as a protein DB for profiling environmental samples using HumanN3.
    '''
    smags_meta_df = pd.read_csv(input_table,sep='\t')
    #Filtering out entries that correspond to METDB, which lack any representation in the original fasta file SMAGs_v1_concat.faa.tar.gz 
    smags_meta_df = smags_meta_df[~smags_meta_df['Genome_Id final names'].str.contains("METDB")]

    #Fetch list of genre reported in the study
    cols2search = ['Best_taxonomy_GENRE','Best_taxonomy_FAMILY','Best_taxonomy_ORDER','Best_taxonomy_CLASS']
    orgs_list = smags_meta_df[cols2search].apply(entry_filter,axis=1).value_counts().index.values
    orgs_list = list(filter(lambda s: s!='NA',orgs_list))

    #Make dict mapping genus name to tax ranks
    best_tax_g2taxlvls = make_fasta_header_templ(orgs_list)

    #Here we filter entries for which NCBI taxID was found at least down to the class level
    smags_tax_ranks_df = smags_meta_df[cols2search]
    filtered_smags_meta_df = smags_meta_df[smags_tax_ranks_df.apply(lambda l: any([i in best_tax_g2taxlvls.keys() for i in list(map(str,l))]),axis=1)]

    #Next we add the header template corresponding to each entry, which will be used to correctly label each seq entry in SMAGs_v1_concat.faa.tar.gz
    pd.options.mode.chained_assignment = None
    header_templates = filtered_smags_meta_df[cols2search].apply(lambda l: orgn2taxranks(l, best_tax_g2taxlvls), axis=1)
    filtered_smags_meta_df["Humann_header_template"] = header_templates
    #Modify entries labelled with TOSAG substring to properly query fasta file
    filtered_smags_meta_df["Genome_Id final names"] = filtered_smags_meta_df["Genome_Id final names"].apply(modify_TOSAG_labels)

    #Query original SMAGs DB and scoop out sequences for particular samples, which are dumped into a another fasta file
    with gzip.open(fasta_fid, "rt") as handle:
        #Loop over seq DB
        for record in bioseqio.parse(handle, "fasta"):
            fasta_id = record.description
            query_tag = "_".join(fasta_id.split("_")[:-1])
            #Get corresponding sample metadata
            sample_metadata = filtered_smags_meta_df[filtered_smags_meta_df["Genome_Id final names"].str.contains(query_tag)]
            #Replace current header with a new one for diamond annotation against UniRef
            if(len(sample_metadata)==1):
                record_len = len(record.seq)
                #Check whether the query seq meets the following criteria
                if(record_len>=seq_len_thr):
                    #Replace fasta header and modify accordingly
                    #NOTE: this provides a template for the diamond annotation against UniRef90 DB
                    seq_header_temp = sample_metadata['Humann_header_template'].iloc[0].replace('UniRef90',fasta_id)
                    new_seq_header = seq_header_temp.replace("SeqLen",str(record_len*3))
                    record.id = ''
                    record.name = ''
                    record.description = ''
                    record.id = new_seq_header
                    #Dump record into file
                    with open(new_fasta_fid, "a") as output_handle:
                        bioseqio.write(record, output_handle, "fasta")       

                        
if __name__ == "__main__":

    #For searches over NCBI, let them know who you are in advance!
    Entrez.email = 'jayson.gutierrez@vliz.be'

    #Loading metadata on SMAGs
    input_table = '/home/VLIZ2000/jayson.gutierrez/meta_omics_data/marine_databases/TARA_Oceans/TARA_SMAGs/Delmont_etal_2020_Table-3.tsv'
    #Set output folder and input/output fasta files
    data_dir = '/home/VLIZ2000/jayson.gutierrez/meta_omics_data/marine_databases/TARA_Oceans/TARA_SMAGs/'
    fasta_fid = os.path.join(data_dir,"SMAGs_v1_concat.faa.tar.gz")
    new_fasta_fid = os.path.join(data_dir,"SMAGs_v1_concat_filtered.faa")
    #Set minimal seq len for a peptide seq to be included in the new/filtered DB
    seq_len_thr=90

    #Run main function
    reformat_SMAGs_db(input_table,fasta_fid,new_fasta_fid,seq_len_thr)
    