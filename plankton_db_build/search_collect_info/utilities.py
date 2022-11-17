#!/home/VLIZ2000/jayson.gutierrez/anaconda3/envs/bioinformatics/bin/python

"""
This module contains a bunch of functions to conduct several multiple tasks linked to fetching particular genomes/protein sequences, taxonomic assignments, 
for mapping pangenomes and and for structural annotation of prokaryotic genomes.

Jayson Gutierrez: 1/2/2021

"""

import Bio.SeqIO as bioseqio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC
from Bio import Entrez
from collections import OrderedDict
from subprocess import Popen, call, STDOUT, PIPE
import urllib
import json
import urllib.request
import requests
from ete3 import NCBITaxa
import os
import numpy as np
import pandas as pd
import glob
import re


def get_tax_id(species):
    """to get data from ncbi taxomomy, we need to have the taxid. we can
    get that by passing the species name to esearch, which will return
    the tax id"""
    species = species.replace(' ', "+").strip()
    try:
        search = Entrez.esearch(term = species, db = "taxonomy", retmode = "xml")
        record = Entrez.read(search)
        return record['IdList'][0]
    except:
        pass


def format_tax4struo(row_data):
    
    """Function to generate a properly formatted NCBI taxonomic for Struo2"""
    
    #Select taxonomic levels to join into one single string
    tax_lvls = ["d__","p__","c__","o__","f__","g__","s__"]
    #formatter function, drop .sp at the end if exists
    tax_str = ";".join(map(lambda l: "".join(l), zip(tax_lvls, row_data))).replace(" sp.","")
    return tax_str


def get_full_tax_lineage(init_taxid=1519565):
    "This function relies on ete3, and is more accurate than the func get_formatted_tax_tree described below"
    #By default, we set taxlevels to unclassified
    ncbi_tax_dict = OrderedDict()
    ncbi_tax_dict['kingdom'] = 'unclassified'
    ncbi_tax_dict['phylum'] = 'unclassified'
    ncbi_tax_dict['class'] = 'unclassified'
    ncbi_tax_dict['order'] = 'unclassified'
    ncbi_tax_dict['family'] = 'unclassified'
    ncbi_tax_dict['genus'] = 'unclassified'
    ncbi_tax_dict['species'] = 'unclassified'      
    
    ncbi = NCBITaxa()    
    lineage = ncbi.get_lineage(init_taxid) 
    names = ncbi.get_taxid_translator(lineage)

    for taxid in lineage:
        rank = list(ncbi.get_rank([taxid]).values())[0]
        if(rank in ncbi_tax_dict.keys()):
            ncbi_tax_dict[rank] = names[taxid]
            
    taxlevels = []
    for k,v in ncbi_tax_dict.items():
        txl = '{}__'.format(k[0]) + v
        taxlevels.append(txl)

    return ';'.join(taxlevels)


def get_super_kingdom(taxid = "253245"):
    ncbi = NCBITaxa()    
    lineage = ncbi.get_lineage(taxid) 
    names = ncbi.get_taxid_translator(lineage)
    sk = list(names.values())[1]
    return sk


def get_formatted_tax_tree(taxid=1519565):
    "This function fetches data from Entrez and return a properly formatted tax tree suitable for Struo/HumanN"
    
    #By default, we set taxlevels to unclassified
    ncbi_tax_dict = OrderedDict()
    ncbi_tax_dict['kingdom'] = 'unclassified'
    ncbi_tax_dict['phylum'] = 'unclassified'
    ncbi_tax_dict['class'] = 'unclassified'
    ncbi_tax_dict['order'] = 'unclassified'
    ncbi_tax_dict['family'] = 'unclassified'
    ncbi_tax_dict['genus'] = 'unclassified'
    ncbi_tax_dict['species'] = 'unclassified'  

    #Fetch the Entrez taxonomy-containing object
    search = Entrez.efetch(id = taxid, db = "Taxonomy", retmode = "xml")
    tax_dict_obj = Entrez.read(search)

    #Grab the antire lineage available for this organism
    lineage = {d['Rank']:d['ScientificName'] for d in tax_dict_obj[0]['LineageEx']}

    #Redefine originally set ncbi_tax_dict 
    for k in lineage.keys():
        if(k in ncbi_tax_dict.keys()):
            ncbi_tax_dict[k] = lineage[k]

    species = tax_dict_obj[0]['ScientificName'] 
    #Species name is this field        
    if(len(species.split(" "))>1):
        ncbi_tax_dict['species'] = species        

    taxlevels = []
    for k,v in ncbi_tax_dict.items():
        txl = '{}__'.format(k[0]) + v
        taxlevels.append(txl)

    return '.'.join(taxlevels)


def get_fields4genome_assmb(db="assembly", term = "Aureococcus anophagefferens", num_entry = "32", 
                            genomes_dir = "/home/VLIZ2000/jayson.gutierrez/meta_omics_data/Struo/PlanktonSeqDB_draft1/Genomes/genbank/",
                            download = True):
    
    """Function that returns the fields required for building a metadata df used as starting point for Struo"""
    assmb_fields = OrderedDict()
   
    handle = Entrez.esearch(db=db, term=term, retmax='200')
    record = Entrez.read(handle)

    id_list = record["IdList"]
    
    if(len(id_list)>0):
        
        if(db == "assembly"):
            
            ref_id = id_list[0]
            esummary_handle = Entrez.esummary(db=db, id=ref_id, report="full")
            esummary_record = Entrez.read(esummary_handle,validate=False)

            #Fetch all fields of interest
            doc_summary = esummary_record["DocumentSummarySet"]["DocumentSummary"][0]

            ncbi_organism_name = doc_summary["SpeciesName"].upper().split(" ")[0] + "_" + num_entry

            ncbi_genbank_assembly_accession = doc_summary["AssemblyAccession"]

            ncbi_species_taxid = int(doc_summary["Taxid"])

            ftp_path_refseq = doc_summary["FtpPath_GenBank"] #doc_summary["FtpPath_RefSeq"] 
            ftp_path = os.path.join(ftp_path_refseq,os.path.basename(ftp_path_refseq) + '_genomic.fna.gz')

            full_ncbi_tax_tree = get_formatted_tax_tree(taxid=ncbi_species_taxid)
            ncbi_taxonomy = "g__" + full_ncbi_tax_tree.split("g__")[-1]

            # assmb_seq = os.path.basename(ftp_refseq) + '_genomic.fna.gz'
            fasta_file_path = os.path.join(genomes_dir,ncbi_genbank_assembly_accession)

            #Grab fields into OrderedDict
            assmb_fields["ncbi_organism_name"] = ncbi_organism_name
            assmb_fields["ncbi_genbank_assembly_accession"] = ncbi_genbank_assembly_accession
            assmb_fields["ncbi_species_taxid"] = ncbi_species_taxid
            assmb_fields["ncbi_taxonomy"] = ncbi_taxonomy
            assmb_fields["ftp_path"] = ftp_path	
            assmb_fields["fasta_file_path"] = fasta_file_path
            

            #If genome assembly is to be donwloaded
            if(download):
                cmd = "wget {} -P {}".format(ftp_path, fasta_file_path)
                runShell(cmd)

            assmb_fpp = os.path.join(fasta_file_path,"*fna.gz")
            assmb_fn = glob.glob(assmb_fpp)

            if(len(assmb_fn)>0):
                #Compute simple genome seq stats using seqkit commandline tool
                seq_stats = get_read_stats4sample(assmb_fn[0])
                assmb_fields["num_contigs"] = int(seq_stats["num_seqs"].replace(",",""))
                assmb_fields["tot_seq_len"] = int(seq_stats["sum_len"].replace(",",""))
            else:
                assmb_fields["num_contigs"] = np.nan
                assmb_fields["tot_seq_len"] = np.nan
            
            assmb_fields["full_ncbi_tax_tree"] = full_ncbi_tax_tree                
                                
            return assmb_fields
        
    else:
        return record


def get_assembly_summary(id):
    """Get esummary for an entrez id"""
    esummary_handle = Entrez.esummary(db="assembly", id=id, report="full")
    esummary_record = Entrez.read(esummary_handle, validate=False)
    return esummary_record


def get_tax_data(taxid):
    """once we have the taxid, we can fetch the record"""
    search = Entrez.efetch(id = taxid, db = "Taxonomy", retmode = "xml")
    return Entrez.read(search)


def format_taxonomy(lineage_ext):
    '''This function receives a dict obj returned by Entrez.esummary and returns the taxonomy formatted in a way 
       compatible with humann2, example: 
       'k__Metazoa;p__Chordata;c__Mammalia;o__Primates;f__Hominidae;g__Homo;s__unclassified'
    '''

    #By default, we set taxlevels to unclassified
    ncbi_tax_dict = OrderedDict()
    ncbi_tax_dict['kingdom'] = 'unclassified'
    ncbi_tax_dict['phylum'] = 'unclassified'
    ncbi_tax_dict['class'] = 'unclassified'
    ncbi_tax_dict['order'] = 'unclassified'
    ncbi_tax_dict['family'] = 'unclassified'
    ncbi_tax_dict['genus'] = 'unclassified'
    ncbi_tax_dict['species'] = 'unclassified'


    #Check which taxlevel exists in the filed LineageEx
    for l in lineage_ext:
        rank = l['Rank']
        if(rank in ncbi_tax_dict.keys()):
            ncbi_tax_dict[rank] = l['ScientificName']

    #Formatting tax levels: Any taxonomy lacking genus and/or species levels will be labeled:
    #g__unclassified (if no genus)
    #s__unclassified (if no species)
    #Required for running humann2
    taxlevels = []
    for k,v in ncbi_tax_dict.items():
        txl = '{}__'.format(k[0]) + v
        taxlevels.append(txl)

    return ';'.join(taxlevels)


def genbank_assembly_accessions_feats(ncbi_organism_name = 'gamma proteobacterium SCGC AAA076-D02'):
    '''Use function to fetch for a given species assembly accession feats:
        - ncbi_organism_name
        - ncbi_genbank_assembly_accession
        - ncbi_species_taxid
        - ncbi_taxonomy
        - ftp_path
        '''
    #Fetch IdList for downstream processing
    handle = Entrez.esearch(db='assembly',term=ncbi_organism_name)
    record = Entrez.read(handle)
    # print(record)
    handle.close()   

    if(int(record['Count'])>0):
    
        ###Fetch id to be passed as arg to esummary
        idl = record['IdList'][0]

        esummary_record = get_assembly_summary(idl)

        doc_summary = esummary_record['DocumentSummarySet']['DocumentSummary'][0]
        # print(doc_summary)

        ncbi_genbank_assembly_accession = doc_summary['AssemblyAccession']

        ###FTP to download refseq
        ftp_path_refseq = doc_summary['FtpPath_GenBank'] 

        ###Get tax_id and use it to format lineage
        tax_id = doc_summary['Taxid']
        tax_data = get_tax_data(tax_id)[0]
        
        #Format lineage
        lineage = format_taxonomy(tax_data['LineageEx'])
        
        #Cast data into a dict
        dict2df = OrderedDict()
        dict2df['ncbi_organism_name'] = ncbi_organism_name
        dict2df['ncbi_genbank_assembly_accession'] = ncbi_genbank_assembly_accession
        dict2df['ncbi_species_taxid'] = tax_id
        dict2df['ncbi_taxonomy'] = lineage

        #get the fasta link - change this to get other formats
        label = os.path.basename(ftp_path_refseq)
        link = os.path.join(ftp_path_refseq,label+'_genomic.fna.gz')
        #NCBI FTP to access fna.gz file
        dict2df['ftp_path'] = link

        return dict2df
    
    else:
        return []


def convert_bytes(num):
    """
    this function will convert bytes to MB.... GB... etc
    """
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
        if num < 1024.0:
            return "%3.1f %s" % (num, x)
        num /= 1024.0


def runShell(cmd, stdout=True):
    if(stdout):
        return Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE).communicate()[0].decode("utf-8")
    else:
        Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE).communicate()[0]


def runShell2(cmd):
    process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = process.communicate()   
    return out, err        
        
        
def get_read_stats4sample(sample_fn):
    '''Get reads stats. Using the command line tool seqkit 
       Alternative for patt='*_singletons_non_rRNA_filtered_reads.fastq'
    '''
    
    out = Popen('seqkit stats {}'.format(sample_fn),
                shell=True,stdout=PIPE,stderr=PIPE).communicate()[0].decode("utf-8")

    #Extract fields
    fields = [h for h in out.split('\n')[0].split('  ') if ''!=h]
    fields = [h if len(h.split(' '))==1 else h.split(' ')[1] for h in fields]
    #Extract atts
    atts = [n for n in out.split('\n')[1].split('  ') if n!='']

    return dict(zip(fields,atts))


def get_metadata4pangenomes(query_orgs_list):
    """Function to create a df containing information on a set of small-sized, simple genomes for a group of
       organisms, from which one can create pangenomes for each group of related species """
    #To collect entries and concatenate all of them in one single df
    meta_data_collect = []

    for org in query_orgs_list:

        dict_collect = []

        search_term = "{}[Orgn]".format(org)
        record = Entrez.read(Entrez.esearch(db="assembly", term=search_term, prop="complete[prop]"))
        ids_list = record["IdList"]
        if(len(ids_list)>1):
            for sid in ids_list:
                summary_id = Entrez.read(Entrez.esummary(db="assembly",id=sid),validate=False)
                id_doc_sum = summary_id["DocumentSummarySet"]["DocumentSummary"][0]

                org_name = id_doc_sum["Organism"]
                ftpath_gbk = id_doc_sum["FtpPath_GenBank"]

                if((('Complete' in id_doc_sum["AssemblyStatus"]) or ('Contig' in id_doc_sum["AssemblyStatus"])) 
                                                                 and (ftpath_gbk!='')):
                    assmb_fn =  os.path.basename(ftpath_gbk) + '_genomic.fna.gz'
                    full_ftpath_gbk = os.path.join(ftpath_gbk,assmb_fn)
                    tax_id = id_doc_sum["Taxid"]
                    tax_tree = get_formatted_tax_tree(taxid=tax_id)
                    entries = [org, org_name, tax_id, tax_tree, full_ftpath_gbk]
                    field_names = ["Organism","Sps_name","Tax_id","Tax_tree","FTP"]
                    dict_collect.append(dict(zip(field_names,entries)))

        if(len(dict_collect)>0):
            #Get metadata for a given organism/genus/family            
            meta_data = pd.DataFrame.from_dict(dict_collect) 
            meta_data_collect.append(meta_data)

    #Now concatenate all the entries into a single df
    orgs4pangenome_df = pd.concat(meta_data_collect)
    orgs4pangenome_df.reset_index(drop=True,inplace=True)
    #Remove duplicates for Tax_id
    orgs4pangenome_df.drop_duplicates(subset='Tax_id',inplace=True)
    
    return orgs4pangenome_df


def run_prokka(orgs4pangenome_df, target = "Prochlorococcus marinus",
               path2folder = "../annotation/Prokka_runs/data/genome_assemblies/",
               path2annot = "../annotation/Prokka_runs/output/gene_predict_annots/"):
    
    """Run prokka on a list of input genomes, the details of which are being stored in df orgs4pangenome_df"""
    
    entries_list = orgs4pangenome_df[orgs4pangenome_df["Organism"]==target][["Tax_id","FTP"]].values
    
    target = target.replace(" ","_")
    dir2store_asmbl = os.path.join(path2folder,target)
    dir2store_annot = os.path.join(path2annot,target)
    

    for (taxid,ftp) in entries_list:

        f2gunzip = os.path.join(dir2store_asmbl,os.path.basename(ftp))

        #Fetch genome assembly from ftp server
        if(~os.path.exists(f2gunzip)):          

            dwl_cmd = "wget {} -P {}".format(ftp,dir2store_asmbl)
            runShell(dwl_cmd)

            #Unzip fasta file
            gunzip_cmd = "gunzip {}".format(f2gunzip)
            runShell(gunzip_cmd)

            print("Assembly: {} gunzipped and ready for annotation via prokka!".format(os.path.basename(f2gunzip)))

        #Set assembly to be processed via prokka
        input_assmbl = f2gunzip.replace(".gz","")

        #Run prokka
        outdir = dir2store_annot + "/NCBI_taxid_" + taxid
        prokka_cmd = "prokka --prefix {0} --locustag {0} --outdir {1} {2}".format(taxid,outdir,input_assmbl)
        runShell(prokka_cmd)

        #Remove downloaded genome assembly
        #runShell("rm {}".format(input_assmbl))

        print("prokka job finished with assembly: {}\n".format(os.path.basename(input_assmbl)))

        
def run_prokka2(pangenome_entries_grby_df, target = "Owenweeksia",
               path2folder = "../annotation/Prokka_runs/data/genome_assemblies/",
               path2annot = "../annotation/Prokka_runs/output/gene_predict_annots/",
               remove_genome = True):
    
    """Run prokka on a list of input genomes, the details of which are being stored in df pangenome_entries_grby_df"""
    
    entries_list = pangenome_entries_grby_df.get_group(target)[["NCBITaxid","FtpPathGenBank"]].values
    
    target = target.replace(" ","_")
    dir2store_asmbl = os.path.join(path2folder,target)
    dir2store_annot = os.path.join(path2annot,target)
    
    for (taxid,ftp) in entries_list:

        f2gunzip = os.path.join(dir2store_asmbl,os.path.basename(ftp))

        #Fetch genome assembly from ftp server
        if(~os.path.exists(f2gunzip)):          

            dwl_cmd = "wget {} -P {}".format(ftp,dir2store_asmbl)
            runShell(dwl_cmd)

            #Unzip fasta file
            gunzip_cmd = "gunzip {}".format(f2gunzip)
            runShell(gunzip_cmd)

            print("Assembly: {} gunzipped and ready for annotation via prokka!".format(os.path.basename(f2gunzip)))

        #Set assembly to be processed via prokka
        input_assmbl = f2gunzip.replace(".gz","")

        #Run prokka
        outdir = dir2store_annot + "/NCBI_taxid_" + str(taxid)
        prokka_cmd = "prokka --prefix {0} --locustag {0} --outdir {1} {2}".format(taxid,outdir,input_assmbl)
        runShell(prokka_cmd)
        
        if(remove_genome):
            #Remove downloaded genome assembly
            runShell("rm {}".format(input_assmbl))

        print("prokka job finished with assembly: {}\n".format(os.path.basename(input_assmbl)))
        

def run_roary_prok_pangenomes(target = "Prochlorococcus marinus",
                              path2inputs = "../annotation/Prokka_runs/output/gene_predict_annots/",
                              path2pang = "../annotation/Prokka_runs/output/pangenomes/"):
    
    """Runing roary on a set of gffs files previously generated via prokka, to create a pangenome 
       (non-redundant gene set of all the species/strains for a given group of organisms)
    """
    
    target = target.replace(" ","_")

    #Set this to be tmp
    tmp_folder4gffs = os.path.join(path2inputs,target)

    #List all the gff files required by roary
    gffs_paths = glob.glob(tmp_folder4gffs + "/NCBI_*/*.gff")

    #Copy all the gff files to this tmp
    cp_cmd = "cp {} {}".format(" ".join(gffs_paths), tmp_folder4gffs)
    runShell(cp_cmd)

    #Run roary on all the gffs
    roary_cmd = "roary -f {} -e -n -v {}/*.gff".format(os.path.join(path2pang,target), tmp_folder4gffs)
    runShell(roary_cmd)

    #Remove all gff files from tmp
    rm_cmd = "rm {}/*.gff".format(tmp_folder4gffs)
    runShell(rm_cmd)
    

def run_prokaryotic_pangenome_pipeline(pangenome_entries_grby_df, target = "Prochlorococcus marinus",
                                       path2folder = "../annotation/Prokka_runs/data/genome_assemblies/",
                                       path2annot = "../annotation/Prokka_runs/output/gene_predict_annots/",
                                       path2pang = "../annotation/Prokka_runs/output/pangenomes/"):
    
    """Runing pangenome pipeline for groups of prokaryotic genomes, using prokka + roary"""


    run_prokka2(pangenome_entries_grby_df, 
                target = target,
                path2folder = path2folder,
                path2annot = path2annot,
                remove_genome = True)

    run_roary_prok_pangenomes(target = target,
                              path2inputs = path2annot,
                              path2pang = path2pang)    


def transl_fna2faa(input_fna_file='../annotation/Prokka_runs/output/pangenomes/Microcystis_flos-aquae/pan_genome_reference.fa',
                   output_faa_file='../annotation/Prokka_runs/output/pangenomes/Microcystis_flos-aquae/pan_genome_reference.faa',
                   transl_table = 11, 
                   min_aa_len = 50):
    """Function that takes every entry in a .fna file and get it translated and saved it to a .faa file"""
    
    with open(input_fna_file,"r") as handle:

            for (i,record) in enumerate(bioseqio.parse(handle, "fasta")):
                #The Bacterial, Archaeal and Plant Plastid Code: table = 11
                #The Chlorophycean Mitochondrial Code: table = 16
                transl_record = record.translate(table=transl_table,to_stop=True,
                                                 name=record.name,
                                                 description=record.description,
                                                 id=record.id)

                if(len(transl_record)>min_aa_len):
                    with open(output_faa_file, "a") as output_handle:
                        bioseqio.write(transl_record, output_handle, "fasta")                     
                        
                        
def queryAphiaRecords(set_id, query_mode="WoRMS2NCBI"):
    '''Function to query WoRMS for mapping data to NCBI or the other way around'''
    
    if(query_mode=="NCBI2WoRMS"): #NCBI to WoRMS
        url = "http://www.marinespecies.org/rest/AphiaRecordByExternalID/{}?type=ncbi".format(set_id)
    
    elif(query_mode=="WoRMS2NCBI"): #WoRMS to NCBI,
        url = "http://www.marinespecies.org/rest/AphiaExternalIDByAphiaID/{}?type=ncbi".format(set_id)
    else: #Otherwise check for attributes of the AphiaID being queried
        url = "http://www.marinespecies.org/rest/AphiaAttributesByAphiaID/{}?include_inherited=true".format(set_id)
    try:
        response = urllib.request.urlopen(url)
        data = json.loads(response.read())
        return data
    except:
        return None


def check4KeywordInAphiaRecord(set_id, kw="Plankton"):
    """Function to check for the presence of a given keyword in an Aphia record given by its ID"""
    
    url = "http://www.marinespecies.org/rest/AphiaAttributesByAphiaID/{}?include_inherited=true".format(set_id)
    data = runShell("curl -X GET {}".format(url))
    
    try: 
        data_dict = json.loads(data)
        kw_inst = []
        for d in data_dict: 
            concat_vals = "_".join(map(str,d.values()))
            if re.search('Meiobenthos', concat_vals, re.IGNORECASE):
                kw_inst.append(1)
        if len(kw_inst)>0:
            return 1
        else:
            return 0
    except:
        return None    
    
    
def queryNCBI4genomic_resources(tax_id = "2316423", db = "assembly"):
    """Function to search NCBI for genomic resources for a given taxID"""
    
    Entrez.email = "jayson.gutierrez@vliz.be"
    
    esummary_handle = Entrez.esummary(db=db, id=tax_id, report="full")
    
    try:
        handle = Entrez.esearch(db=db,term=tax_id)
        record = Entrez.read(handle)
        if(len(record['IdList'])>0):
            return record['IdList']
    except:
        return None

    
def get_genomic_resources_from_NCBI(tax_id = "2316423", db = "assembly"):
    """Check if any genomic resource is available for the tax_id"""
    genomic_resource_ids_collect = []
    genomic_resource_ids = queryNCBI4genomic_resources(tax_id = tax_id, db = db)
    if(genomic_resource_ids!=None):
        for gr_id in genomic_resource_ids:
            genomic_resource_ids_collect.append(gr_id)
    
    return genomic_resource_ids_collect


def search_ncbi4genomic_resources(kw_query_list):
    """Function to query NCBI for genomic resources (genome assembly/annotated proteome) based on a list of
       keywords, which can be e.g. family/genus name, etc."""

    ref_tax_ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

    entry_sum_collect = []
    ncbi_taxid_collect = []

    for search_kw in kw_query_list:

        record = Entrez.read(Entrez.esearch(db="assembly", term="{}[Orgn]".format(search_kw), prop="complete[prop]"))
        ids_list = record["IdList"]

        if(len(ids_list)>0):

            for sid in ids_list:

                summary_id = Entrez.read(Entrez.esummary(db="assembly",id=sid),validate=False)
                id_doc_sum = summary_id["DocumentSummarySet"]["DocumentSummary"][0]             
                entry_sum = OrderedDict()

                #If entry has a genome assembly append to annot_gs_collect
                if ((id_doc_sum["FtpPath_GenBank"] != '') or (id_doc_sum["FtpPath_RefSeq"] != '')):
                    
                    try: 
                        entry_sum["NCBITaxid"] = id_doc_sum["Taxid"]
                        entry_sum["SpeciesName"] = id_doc_sum["SpeciesName"]
                        entry_sum["Organism"] = id_doc_sum["Organism"]
                        entry_sum["SuperKingdom"] = get_super_kingdom(id_doc_sum["Taxid"])

                        #Fetch tax ranks individually
        #                 entry_sum["FullLineage"] = get_full_tax_lineage(id_doc_sum["Taxid"])
                        conc_tax_lineage = get_full_tax_lineage(id_doc_sum["Taxid"])
                        calc_tax_ranks = list(map(lambda s: s.split("__")[1], conc_tax_lineage.split(";")))
                        #Cast into dict the tax ranks computed
                        calc_tax_ranks_dict = dict(zip(ref_tax_ranks, calc_tax_ranks))                
                        for taxr in ref_tax_ranks:
                            entry_sum[taxr] = calc_tax_ranks_dict[taxr]

                        #Check which genomic resources: only unannotated genome assembly vs annotated proteome
                        if (id_doc_sum["FtpPath_RefSeq"] != ''):
                            ftp_path = os.path.join(id_doc_sum["FtpPath_RefSeq"], 
                                                    "{}_protein.faa.gz".format(os.path.basename(id_doc_sum["FtpPath_RefSeq"])))

                            entry_sum["FtpPathGenBank"] = ftp_path
                            entry_sum["MolecularResourceType"] = "Annotated Proteome"

                        elif (id_doc_sum["FtpPath_GenBank"] != ''):
                            ftp_path = os.path.join(id_doc_sum["FtpPath_GenBank"],
                                                    "{}_genomic.fna.gz".format(os.path.basename(id_doc_sum["FtpPath_GenBank"])))
                            entry_sum["FtpPathGenBank"] = ftp_path
                            entry_sum["MolecularResourceType"] = "Genome Assembly"

                        entry_sum["LastMajorReleaseAccession"] = id_doc_sum["LastMajorReleaseAccession"]
                        entry_sum["AssemblyName"] = id_doc_sum["AssemblyName"]
                        entry_sum["AssemblyType"] = id_doc_sum["AssemblyType"]
                        entry_sum["AssemblyStatus"] = id_doc_sum["AssemblyStatus"]
        #                 entry_sum["BioprojectAccn"] = id_doc_sum["GB_BioProjects"][0]["BioprojectAccn"]
                        entry_sum["Coverage"] = id_doc_sum["Coverage"]
                        entry_sum["AsmReleaseDate_GenBank"] = id_doc_sum["AsmReleaseDate_GenBank"]
                        entry_sum["LastUpdateDate"] = id_doc_sum["LastUpdateDate"]
                        entry_sum["SubmitterOrganization"] = id_doc_sum["SubmitterOrganization"]
                        entry_sum["RefSeq_category"] = id_doc_sum["RefSeq_category"]
                        entry_sum["ContigN50"] = id_doc_sum["ContigN50"]
                        entry_sum["ScaffoldN50"] = id_doc_sum["ScaffoldN50"]
                        #Retrieve infor from WORMS
                        aphia_record = queryAphiaRecords(id_doc_sum["Taxid"], query_mode="NCBI2WoRMS")
                        entry_sum["AphiaTaxid"] = np.nan
                        entry_sum["IsMarine"] = np.nan
                        entry_sum["IsBrackish"] = np.nan
                        entry_sum["IsFreshwater"] = np.nan
                        if(aphia_record != None):
                            entry_sum["AphiaTaxid"] = aphia_record["AphiaID"]
                            entry_sum["IsMarine"] = aphia_record["isMarine"]
                            entry_sum["IsBrackish"] = aphia_record["isBrackish"]
                            entry_sum["IsFreshwater"] = aphia_record["isFreshwater"]


                        if(entry_sum["NCBITaxid"] not in ncbi_taxid_collect):
                            entry_sum_collect.append(entry_sum)
                    
                    except:
                        
                        pass

    #Cast full list as DF
    entry_sum_collect_df = pd.DataFrame.from_dict(entry_sum_collect)
    entry_sum_collect_df.drop_duplicates(subset=["NCBITaxid"], inplace=True)

    return entry_sum_collect_df



def fetch_assemblies(ftp, download=True, kingdom = "prokaryotes", 
                     path2assembl = "/data/home/VLIZ2000/jayson.gutierrez/meta_omics_data/Struo2/data/custom_db/genomes"):
    
    """Function to fetch genome assemblies and place them on a given destination"""
    
    dir2store_asmbl = os.path.join(path2assembl,kingdom)

    f2gunzip = os.path.join(dir2store_asmbl,os.path.basename(ftp))

    if(download):
        #Fetch genome assembly from ftp server
        if(~os.path.exists(f2gunzip)):          

            dwl_cmd = "wget {} -P {}".format(ftp,dir2store_asmbl)
            runShell(dwl_cmd)
    
    return f2gunzip