## Code to create tissue-specific sets of cCREs, pooling across multiple samples to minimize noise
## Script includes pulling data from the ENCODE API 

import os
import json
import urllib
from urllib import request, parse, error
import requests
import sys
from tqdm import tqdm, trange
from multiprocessing import Pool
import pybigtools
import pandas as pd
import numpy as np
import math
import statistics
from plotnine import *
from plotnine.data import diamonds
from mizani.formatters import percent_format

## DNase API pull
def Extract_DNase_Experiment_Metatdata(exp, genome):
    #dataDir="/data/projects/encode/json/exps/"+exp
    dataDir="https://www.encodeproject.org/experiments/"+exp+"/?format=json"
    json_data=open(dataDir+".json").read()
    data = json.loads(json_data) 
    ## Noting here that this is the name of the biosample as it would appear in SCREEN
    biosample_summary = data['biosample_summary']) 

    for entry in data["files"]:
        if entry['file_format'] == 'bigWig':
            if entry['analyses'][0]['pipeline_award_rfas'] == ['ENCODE4']:
                if entry['biological_replicates'] == [1]:
                    acc = entry['accession']
                    #file_path = f'/data/projects/encode/data/{exp}/{acc}.bigWig'
                    file_path = f'https://www.encodeproject.org/files/{acc}/@@download/{acc}.bigWig'  ## You can use this path to download bigwigs from ENCODE
                    return(file_path)

def pull_DNase_files(tissue):
    genome = "GRCh38"
    species = "Homo+sapiens"
    
    ## Build query
    urlMain = "https://www.encodeproject.org/search/?type=Experiment&status=released" + \
        "&perturbed=false&perturbed=true&assay_title=DNase-seq" + \
        "&replicates.library.biosample.donor.organism.scientific_name=" + species + \
        "&format=json&limit=all&biosample_ontology.organ_slims=" + tissue + \
        "&replicates.library.biosample.life_stage=adult&replicates.library.biosample.life_stage=unknown" + \ ## for my purposes, I only wanted adult data, but you can remove this line if you also want to download fetal/child data
        "&biosample_ontology.classification=tissue&biosample_ontology.classification=primary+cell" ## I also limited to tissues and primary cells because I didn't want data from cell lines
    
    response = urllib.request.urlopen(urlMain)
    data = json.loads(response.read())

    paths = []
    for entry in data["@graph"]: #loops through experiments
        paths.append(Extract_DNase_Experiment_Metatdata(entry["accession"], genome))
    return(paths)

## ATAC API pull 
def Extract_ATAC_Experiment_Metatdata(exp, genome):
    #dataDir="/data/projects/encode/json/exps/"+exp
    dataDir="https://www.encodeproject.org/experiments/"+exp+"/?format=json"
    json_data=open(dataDir+".json").read()
    data = json.loads(json_data)
    ## Noting here that this is the name of the biosample as it would appear in SCREEN
    biosample_summary = data['biosample_summary']) 

    for entry in data["files"]:
        if entry['file_format'] == 'bigWig':
            if entry['output_type'] == 'fold change over control':
                if entry['analyses'][0]['pipeline_award_rfas'] == ['ENCODE4']:
                    if entry['biological_replicates'] == [1]:
                        acc = entry['accession']
                        #file_path = f'/data/projects/encode/data/{exp}/{acc}.bigWig'
                        file_path = f'https://www.encodeproject.org/files/{acc}/@@download/{acc}.bigWig' ## You can use this path to download bigwigs from ENCODE
                        return(file_path)

def pull_ATAC_files(tissue):
    genome = "GRCh38"
    species = "Homo+sapiens"
    
    ## Build query
    urlMain = "https://www.encodeproject.org/search/?type=Experiment&status=released" + \
        "&perturbed=false&perturbed=true&assay_title=ATAC-seq" + \
        "&replicates.library.biosample.donor.organism.scientific_name=" + species + \
        "&format=json&limit=all&biosample_ontology.organ_slims=" + tissue + \
        "&replicates.library.biosample.life_stage=adult&replicates.library.biosample.life_stage=unknown" + \ ## for my purposes, I only wanted adult data, but you can remove this line if you also want to download fetal/child data
        "&biosample_ontology.classification=tissue&biosample_ontology.classification=primary+cell" ## I also limited to tissues and primary cells because I didn't want data from cell lines
    
    response = urllib.request.urlopen(urlMain)
    data = json.loads(response.read())

    paths = []
    for entry in data["@graph"]: #loops through experiments
        paths.append(Extract_ATAC_Experiment_Metatdata(entry["accession"], genome))
    return(paths)

## Notes on pulling other file types:
## For ChIP-seq, you can modify the ATAC-seq query (both search for fold-change over control bigwigs),
## replacing "&assay_title=ATAC-seq" in the url
## with "&assay_title=Histone+ChIP-seq&target.label=H3K27ac" or 
## "&assay_title=Histone+ChIP-seq&target.label=H3K4me3" or 
## "&assay_title=TF+ChIP-seq&target.label=CTCF"

## Open bigwig and compute signal and z-scores (but just return zscores based on alphabetical order of rDHSs)
def compute_zscores(bigWig):
    with pybigtools.open(bigWig, "r") as bw:
        x = list(bw.average_over_bed('/data/projects/encode/Registry/V4/GRCh38/GRCh38-Anchors.bed', names=4, stats='mean'))
    df = pd.DataFrame(x, columns=["rDHS", "signal"])
    
    signalmean = statistics.mean([ math.log(x) for x in df["signal"] if x > 0.0 ])
    signalstd = statistics.stdev([ math.log(x) for x in df["signal"] if x > 0.0 ])
    
    df["zscore"] = [ (math.log(x) - signalmean) / signalstd if x > 0.0 else -10.0 for x in df["signal"] ]
    
    return(df.sort_values(by='rDHS')['zscore'].to_list())

## Parallelize any function
def run_imap_multiprocessing(func, argument_list, num_processes):

    pool = Pool(processes=num_processes)

    result_list_tqdm = []
    for result in tqdm(pool.imap(func=func, iterable=argument_list), total=len(argument_list)):
        result_list_tqdm.append(result)

    return result_list_tqdm

## Create histogram visualizing active cCREs
def plot_set(tmp, tissue, assay):
    exps = tmp.columns
    tmp['counts'] = tmp[exps].ge(1.64).sum(axis=1)
    tmp['Freq5'] = tmp[exps].ge(1.64).sum(axis=1) >= 5
    tmp['BioZ'] = tmp[exps].ge(2.32).sum(axis=1) >= 1
    num_active = len(tmp)
    num_freq = len(tmp[tmp.Freq5])
    num_bioZ = len(tmp[~tmp.Freq5][tmp.BioZ])
    tmp['color'] = 'not included'
    tmp.loc[tmp.BioZ, 'color'] = 'Z >= 2.32'
    tmp.loc[tmp.Freq5, 'color'] = '>= 5 biosamples'

    p = (
        ggplot(tmp, aes(x="counts", fill="color"))
        + geom_histogram(binwidth=1)
        + theme_classic() 
        + ggtitle(tissue)
        + xlab('Number of Supporting Experiments')
        + ylab('Number of cCREs')
        + labs(subtitle=f'{num_active:,} cCREs active in any biosample,\n{num_freq:,} cCREs active in >= 5 biosamples,\n{num_bioZ:,} additional cCREs with Z >= 2.32')
    )

    p.save(f'/data/zusers/sheddn/TCGA/healthy_tissue_sets/plots/{tissue}_{assay}.pdf')

## Compute consensus active cCREs any matrix of zscores
def compute_set(df, tissue, assay):
    tissue = tissue.replace('+','-')

    plot_set(df[df.ge(1.64).sum(axis=1) >= 1], tissue, assay)
    
    freq = df.ge(1.64).sum(axis=1) >= 5
    bioz = df.ge(2.32).sum(axis=1) >= 1

    ccre_set = ccres[ccres.cCRE.isin(df[freq | bioz].index)]
    ccre_set.to_csv(f'/data/zusers/sheddn/TCGA/healthy_tissue_sets/{tissue}_{assay}_cCREs.bed',
                    sep='\t', header=False, index=False)
    

anchors = pd.read_csv('/data/projects/encode/Registry/V4/GRCh38/GRCh38-Anchors.bed',
                      sep='\t', header=None, usecols=[3])[3].to_list()
ccres = pd.read_csv('/data/projects/encode/Registry/V4/GRCh38/GRCh38-cCREs.bed',
                      sep='\t', header=None, names=['chr','start','end','rDHS','cCRE','class'])

DNase_tissues = ['adrenal+gland','breast','liver','colon','esophagus',
                 'kidney','lung','prostate+gland','skin+of+body','stomach',
                 'testis','thyroid+gland','uterus']

for DNase_tissue in DNase_tissues:
    print(DNase_tissue)
    file_paths = pull_DNase_files(DNase_tissue)
    zscores = run_imap_multiprocessing(compute_zscores, file_paths, 12)

    zscores_df = pd.DataFrame(zscores).transpose()
    zscores_df.index = anchors
    zscores_df.columns = [x.replace('/data/projects/encode/data/','').replace('/','-').replace('.bigWig','') for x in file_paths]
    zscores_df = zscores_df[zscores_df.index.isin(ccres.rDHS)]
    
    compute_set(zscores_df, DNase_tissue, 'DNase')


ATAC_tissues = ['adrenal+gland','urinary+bladder','breast','liver','colon',
                'esophagus','kidney','lung','prostate+gland','stomach',
                'testis','thyroid+gland','uterus', 'brain']

for ATAC_tissue in ATAC_tissues:
    print(ATAC_tissue)
    file_paths = pull_ATAC_files(ATAC_tissue)
    zscores = run_imap_multiprocessing(compute_zscores, file_paths, 12)

    zscores_df = pd.DataFrame(zscores).transpose()
    zscores_df.index = anchors
    zscores_df.columns = [x.replace('/data/projects/encode/data/','').replace('/','-').replace('.bigWig','') for x in file_paths]
    zscores_df = zscores_df[zscores_df.index.isin(ccres.rDHS)]
    
    compute_set(zscores_df, ATAC_tissue, 'ATAC')
