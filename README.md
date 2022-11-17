Metabolic pathway analysis of planktonic ecosystems based on metagenomes (MTGs)/metatranscriptomes (MTXs) derived from marine environments 
======

# Introduction
Adequate trait parameterization of microbial community structure and biogeochemical potential under hypothesized climate change scenarios is critical for prediction purposes via e.g. next-generation ecosystem biology models able to integrate features derived from omics analysis. Systems biology-inspired profiling tools such as [HUMAnN](https://elifesciences.org/articles/65088) can deliver this information by profiling microbial activity at the level of metabolic pathways (as opposed to the individual transcript level). However, current application of these tools is limited to well studied prokaryotic assemblages. Here, by combining information extracted from global eukaryotic reference sequence databases and a suite of state-of-the-art bioinformatics tools, we were able to profile the metabolic pathways of poorly characterized natural eukaryotic communities.

<br>
> WARNING:  Many of the tools in this repo may undergo substantial changes over time without warning in order to meet increasing demand for computational efficiency and/or data processing/analysis capabilities!

# Directory structure

#### /plankton_db_build
  * /search_collect_info: 
     1. Code organized as a series of scripts and Jupyter notebooks implemented to search, fetch and process plankton genome assemblies available through NCBI and other types of sequence resources, including the [EukZoo](https://zenodo.org/record/1476236#.YfpQwv7MJPY) and the [TARA SMAGs](https://www.genoscope.cns.fr/tara/) databases.
     2. This subfolder also contains a series of utilities for dealing with taxonomic assignments of genomic/gene sequences, as well as for parsing different sorts of information fetched from NCBI and the WoRMs catalogue (Aphia Information) via Biopython and their APIs.  
     
  * /genome_annotation:
    1. This subfolder contains a series of scripts developed to make gene prediction/calling on eukaryotic genomes using [BUSCO](https://busco.ezlab.org/) + [Augustus](https://bioinf.uni-greifswald.de/augustus/).
    2. Here we can also find utilities for dealing with the different intermediate files generated while the pipeline is being deployed on a HPC system.  
    
  * /funct_profil_data_analysis:
     1. This subfolder contains a series of jupyter notebooks, scripts and modules used to analyze three case studies (metagenomics datasets), where we demonstrate the potential of our systems biology-inspired approach to map the metabolic activity of complex planktonic communities (see description within each notebook).

#### /fastq_qc_humann_profiling 
  1. Scripts to conduct QC of (short) sequence reads using the [bbtools suite](https://jgi.doe.gov/data-and-tools/bbtools/).
  2. Scripts to conduct functional profiling of QC MTG/MTX datasets using humann3 (https://github.com/biobakery/humann).

#### Protein database
A plankton-specific database, with a format compatible with the profiling tool [HUMAnN](https://elifesciences.org/articles/65088) can be accessed [here](https://doi.org/10.14284/573): 

#### NOTE: 
Each file provides documentation on what the code does and how to utilize it. 
