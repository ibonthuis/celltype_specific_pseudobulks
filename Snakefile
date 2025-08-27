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


## Output files ##
FILTERED_COUNTMATRIX_RDATA = os.path.join(FILTERED_COUNTMATRIX_OUTPUT_DIR, "countmatrix_filtered.RData")
FILTERED_METADATA_TSV = os.path.join(FILTERED_COUNTMATRIX_OUTPUT_DIR, "metadata_single_cells_filtered.tsv")
LIST_OF_SUBSAMPLES_COUNT_RDATA = os.path.join(SUBSAMPLES_BASE_DIR, "subsamples.RData")

## Rule ALL ##
rule all:
    input:
        FILTERED_COUNTMATRIX_RDATA, \
        FILTERED_METADATA_TSV, \
        LIST_OF_SUBSAMPLES_COUNT_RDATA


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
    This rule makes sure to run one batch of subsampling. This rule is aimed to run 10x, so I have 10 varieties of subsamples (random). As input it takes count_matrix (dgCmatrix is fine) and metadata, both filtered in the create_count_matrix rule.
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
        output_dir = SUBSAMPLES_BASE_DIR
    shell:
        """
        Rscript {params.bin}/run_subsampling.R \
            -c {input.count_matrix} \
            -m {input.metadata} \
            -o {params.output_dir}
        """

# rule pseudobulk:
#     """
#     This rule creates pseudobulks from single cell count matrices, prepped in the run_subsampling rule. It takes all the subsamples of one subsampling run as input and creates pseudobulks for every patients and aggregates them by cell type. The RData file contains a list per cell type??? Or this rule outputs per celltype?
#     """
#     input:
#         list_of_subsamples = LIST_OF_SUBSAMPLES_COUNT_RDATA
#     output:
#         pseudobulk = PSEUDOBULK_RDATA
#         pseudo_metadata = PSEUDOBULK_METADATA_RDATA
#     message:
#         "; pseudobulking"
#     params:
#         bin = os.path.join(config["bin"]), \
#         output_dir = SUBSAMPLES_BASE_DIR
#     shell:
#         """
#         echo "; This file is input for pseudobulks" ;
#         echo {input.list_of_subsamples}
#         Rscript {params.bin}/pseudobulk.R \
#             -s {input.list_of_subsamples} \
#             -o {params.output_dir}
#         """