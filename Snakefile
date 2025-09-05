## How to run this?

## module load R/4.4.1-gfbf-2023b 
## module load snakemake/7.23.1-foss-2022a or snakemake/8.27.0-foss-2024a 
## snakemake --cores 1 -np ### For dry run

## Libraries
import os 
import sys
import glob
from pathlib import Path
import time

## Config

global CONFIG_PATH
CONFIG_PATH = "config.yaml"
configfile: CONFIG_PATH

## Container

## Directories
DATA_DIR = config["data_dir"]
BASE_OUTPUT_DIR = config["output_dir"]
FILTERED_COUNTMATRIX_OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, "countmatrix_filtered/")
SUBSAMPLES_BASE_DIR = os.path.join(BASE_OUTPUT_DIR, "subsampled/")

## Input files ##
INPUT_SEURAT = os.path.join(DATA_DIR, "singlecell.RData")

## Other inputs ##
CELLTYPES = config["celltypes"]
RUN_NR = config["subsample_run"]


## Output files ##
FILTERED_COUNTMATRIX_RDATA = os.path.join(FILTERED_COUNTMATRIX_OUTPUT_DIR, "countmatrix_filtered.RData")
FILTERED_METADATA_TSV = os.path.join(FILTERED_COUNTMATRIX_OUTPUT_DIR, "metadata_single_cells_filtered.tsv")
LIST_OF_SUBSAMPLES_COUNT_RDATA = os.path.join(SUBSAMPLES_BASE_DIR, "{subsample_run}", "subsamples.RData")
PSEUDOBULK_TSV = os.path.join(SUBSAMPLES_BASE_DIR, "{subsample_run}", "{celltypes}", "pseudobulk.tsv")
PSEUDOBULK_METADATA_TSV = os.path.join(SUBSAMPLES_BASE_DIR, "{subsample_run}", "{celltypes}", "pseudobulk_metadata.tsv")

## Rule ALL ##
rule all:
    input:
        FILTERED_COUNTMATRIX_RDATA, \
        FILTERED_METADATA_TSV, \
        expand(LIST_OF_SUBSAMPLES_COUNT_RDATA, subsample_run = RUN_NR), \
        expand(PSEUDOBULK_TSV, subsample_run = RUN_NR, celltypes = CELLTYPES), \
        expand(PSEUDOBULK_METADATA_TSV, subsample_run = RUN_NR, celltypes = CELLTYPES)

## Rules ##

rule create_count_matrix:
    """ 
    This rule takes as input a seurat object and filters the counts based on celltypes of interest and a threshold of genes expressed in a minimal amount of cells (1%).
    And outputs the filtered counts, and in the future hopefully also metadata.
    
    Inputs
    ------
    INPUT_SEURAT
        Path to file containing Seurat in RData format
    Outputs
    -------
    SUBSAMPLE_COUNTMATRIX_RDATA
        File containing a filtered count matrix, and hopefully in the future also subsampled count_matrices (as a list or other format) in RData format   
    """
    input:
        seurat = INPUT_SEURAT
    output:
        FILTERED_COUNTMATRIX_RDATA, \
        FILTERED_METADATA_TSV
    message:
        "; Running count matrix filtering on seurat object."
    params:
        bin = os.path.join(config["bin"]), \
        output_dir = FILTERED_COUNTMATRIX_OUTPUT_DIR, \
        celltypes = CELLTYPES
    shell:
        """
        celltype_list=$(echo "{params.celltypes}" | tr ' ' ',')
        echo $celltype_list 
        Rscript {params.bin}/create_count_matrix.R \
            -s {input.seurat} \
            -c $celltype_list \
            -o {params.output_dir}
        """

rule run_subsampling:
    """
    This rule makes sure to run one batch of subsampling. This rule is aimed to run 10 times, so I have 10 varieties of subsamples (random). As input it takes a count_matrix (dgCmatrix is fine) and metadata (tab separated), both filtered in the create_count_matrix rule.
    """
    input:
        count_matrix = FILTERED_COUNTMATRIX_RDATA, \
        metadata = FILTERED_METADATA_TSV
    output:
        list_of_subsamples = LIST_OF_SUBSAMPLES_COUNT_RDATA
    message:
        "; running subsampling."
    params:
        bin = os.path.join(config["bin"]), \
        output_dir = os.path.join(SUBSAMPLES_BASE_DIR, "{subsample_run}")
    shell:
        """
        Rscript {params.bin}/run_subsampling.R \
            -c {input.count_matrix} \
            -m {input.metadata} \
            -o {params.output_dir}
        """

rule create_pseudobulk:
    """
    This rule creates pseudobulk from single cell count matrices, prepped in the run_subsampling rule. It takes all the subsamples of one subsampling run as input and creates pseudobulk for every patient and aggregates them by cell type, the celltype as defined by the wildcard. A pseudobulk per celltype is stored in the RData file. 
    """
    input:
        list_of_subsamples = LIST_OF_SUBSAMPLES_COUNT_RDATA
    output:
        PSEUDOBULK_TSV, \ 
        PSEUDOBULK_METADATA_TSV
    message:
        "; pseudobulking"
    params:
        bin = os.path.join(config["bin"]), \
        output_dir = os.path.join(SUBSAMPLES_BASE_DIR, "{subsample_run}", "{celltypes}")
    shell:
        """
        echo "; This file is input for pseudobulk" ;
        echo {input.list_of_subsamples}
        Rscript {params.bin}/pseudobulk.R \
            -s {input.list_of_subsamples} \
            -o {params.output_dir}
        """


# rule construct_networks:
#     """
#     This rule should take as input the pseudobulk dataframe per cell type from the pseudobulk rule. Then networks should be built for each pseudobulk using sisana in this rule.
#     """
#     input:
#         pseudobulk = PSEUDOBULKS_TSV
#     output:
#         networks
#     shell:
#         """
#         conda activate sisana_sept2025
        
#         """