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
PYTHON_CONTAINER = config["python_container"]


## Directories
DATA_DIR = config["data_dir"]
YAML_DIR = config["yaml_dir"]
BASE_OUTPUT_DIR = config["output_dir"]
FILTERED_COUNTMATRIX_OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, "countmatrix_filtered/")
SUBSAMPLES_BASE_DIR = os.path.join(BASE_OUTPUT_DIR, "subsampled/")
PSEUDOBULK_DIR = os.path.join(SUBSAMPLES_BASE_DIR, "{subsample_run}", "{celltypes}")
NETWORK_DIR = os.path.join(SUBSAMPLES_BASE_DIR, "{subsample_run}", "{celltypes}", "network")

## Input files ##
INPUT_SEURAT = os.path.join(DATA_DIR, "singlecell.RData")
INPUT_YAML = os.path.join(YAML_DIR, "params.yml")

## Other inputs ##
CELLTYPES = config["celltypes"]
RUN_NR = config["subsample_run"]

# # SiSaNA outputs
# EXPRESSION_FILTERED = os.path.join(SISANA_DIR, "preprocess", EXP + "_filtered.txt")
# MOTIF_PRIOR_FILTERED = os.path.join(SISANA_DIR, "preprocess", MOTIF_PRIOR_TAG + "_filtered.txt")
# PPI_PRIOR_FILTERED = os.path.join(SISANA_DIR, "preprocess", PPI_PRIOR_TAG + "_filtered.txt")
# STATS = os.path.join(SISANA_DIR, "preprocess", EXP + "_filtering_statistics.txt")
# PANDA_NET = os.path.join(SISANA_DIR, "network", "panda_network.txt")
# MOTIF_PRIOR_FILTERED = os.path.join(SISANA_DIR, "preprocess", MOTIF_PRIOR_TAG + "_filtered.txt")

## Output files ##
FILTERED_COUNTMATRIX_RDATA = os.path.join(FILTERED_COUNTMATRIX_OUTPUT_DIR, "countmatrix_filtered.RData")
FILTERED_METADATA_TSV = os.path.join(FILTERED_COUNTMATRIX_OUTPUT_DIR, "metadata_single_cells_filtered.tsv")
LIST_OF_SUBSAMPLES_COUNT_RDATA = os.path.join(SUBSAMPLES_BASE_DIR, "{subsample_run}", "subsamples.RData")
PSEUDOBULK_TSV = os.path.join(SUBSAMPLES_BASE_DIR, "{subsample_run}", "{celltypes}", "pseudobulk.tsv")
PSEUDOBULK_METADATA_TSV = os.path.join(SUBSAMPLES_BASE_DIR, "{subsample_run}", "{celltypes}", "pseudobulk_metadata.tsv")
#LIONESS_NPY = os.path.join(SUBSAMPLES_BASE_DIR, "{subsample_run}", "{celltypes}", "network", "lioness.npy")
COPIED_PARAMS = os.path.join(SUBSAMPLES_BASE_DIR, "{subsample_run}", "{celltypes}", "params.yml")


## Rule ALL ##
rule all:
    input:
        FILTERED_COUNTMATRIX_RDATA, \
        FILTERED_METADATA_TSV, \
        expand(LIST_OF_SUBSAMPLES_COUNT_RDATA, subsample_run = RUN_NR), \
        expand(PSEUDOBULK_TSV, subsample_run = RUN_NR, celltypes = CELLTYPES), \
        expand(PSEUDOBULK_METADATA_TSV, subsample_run = RUN_NR, celltypes = CELLTYPES), \
        #expand(LIONESS_NPY, subsample_run = RUN_NR, celltypes = CELLTYPES),\
        expand(COPIED_PARAMS, subsample_run = RUN_NR, celltypes = CELLTYPES)
        

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


# rule preprocess_networks:
#     """
#     This rule should take as input the pseudobulk dataframe per cell type from the pseudobulk rule.  
#     """
#     input:
#         yaml_file = "params.yml"
#     output:
#         # I have to check whether sisana automatically 
#         preprocessed_file = "pseudobulk_preprocessed.txt"

#     params:
#         output_dir = os.path.join(SUBSAMPLES_BASE_DIR, "{subsample_run}", "{celltypes}")
#     shell:
#         """
#         conda activate sisana_sept2025
#         sisana preprocess config.yaml
#         conda deactivate
#         """

# rule build_networks:
#     """
#     Then networks should be built for each pseudobulk using sisana in this rule.
#     """
#     input:
#         PSEUDOBULK_TSV
#         # preprocessed_file = "pseudobulk_preprocessed.txt", \
#         # preprocessed_motif_prior = "motif_prior_names_2024.tsv", \
#         # preprocessed_ppi = "ppi_prior_2024.tsv"
#     output:
#         LIONESS_NPY#, \
#         # lioness_pickle = "lioness.pickle", \
#         # lioness_indegree = "lioness_indegree.csv", \
#         # lioness_outdegree = "lioness_outdegree.csv", \
#         # panda_output = "panda_network.txt"
#     params:
#         "config.yaml"
#         #output_dir = os.path.join(SUBSAMPLES_BASE_DIR, "{subsample_run}", "{celltypes}")
#     # conda:
#     # #     #"sisana_env.yaml"
#     #     "/storage/kuijjerarea/ine/projects/BRCA_SINGLECELL_TO_PSEUDOBULK/celltype_specific_pseudobulks/sisana_sept2025"
#     shell:
#         """
#         #conda run -p /storage/kuijjerarea/ine/projects/BRCA_SINGLECELL_TO_PSEUDOBULK/celltype_specific_pseudobulks/sisana_sept2025 sisana generate config.yaml
#         echo {PSEUDOBULK_TSV}
#         """
#         #/storage/kuijjerarea/ine/projects/BRCA_SINGLECELL_TO_PSEUDOBULK/celltype_specific_pseudobulks/ # path for sisana env
#         #sisana generate config.yaml


# rule list_inputs:
#     """
#     This script should echo the filepaths I need for my params.yml files, as PSEUDOBULK_TSV is made up of wild cards
#     """
#     input: 
#         PSEUDOBULK_TSV
#         params_file = "params.yml"
#     shell:
#         """
#         echo PSEUDOBULK_TSV
#         """

# rule place_yaml_file:
#     """
#     This script should copy a basic yaml file into the folders specified by wildcards
#     """
#     input: 
#         input_params = INPUT_YAML
#     output:
#         output_params = COPIED_PARAMS
#     params:
#        # output_dir = PSEUDOBULK_DIR, \
#         subsampleruns_dir = SUBSAMPLES_BASE_DIR
#     shell:
#         """
#         echo ${subsampleruns_dir}
#         echo ${input_params}
#         echo ${output_params}

#         """

        # for dir in ${subsampleruns_dir}
        # do
        #     # Get the directory name
        #     directory=${dir}
        #     echo "$directory"
            
        #     # Copy the params file to the directory
        #     cp ${input_params} $directory

        #     # Run sisana
        #     path_to_params=$directory/params.yml
        #     echo "$path_to_params"
        #     # sisana generate $path_to_params
        # done

rule place_yaml_file:
    """
    This script should copy a basic yaml file into the folders specified by wildcards
    """
    input: 
        params_yaml = INPUT_YAML
    output:
      #  output_params = f"{params.subsampleruns_dir}/subsampled/{wildcards.subsample_run}/{wildcards.celltypes}/params.yml"
        output_params = COPIED_PARAMS
    # params:
    #     subsampleruns_dir = SUBSAMPLES_BASE_DIR
    shell:
        """
        cp {input.params_yaml} {output.output_params}
        """

# rule run_sisana:
#     """
#     With the params.yml files copied into every directory together with the pseudobulk.tsv, it should be able to run sisana, if the environment can correctly be identified.
#     Outputs
#     -------
#     EXPRESSION_FILTERED:
#         A TXT file with filtered expression.
#     MOTIF_PRIOR_FILTERED:
#         A TXT file with filtered motif prior.
#     PPI_PRIOR_FILTERED:
#         A TXT file with filtered PPI prior.
#     STATS:
#         A file containing information about genes filtered.
#     PANDA_NET:
#         A TXT file with the PANDA network.
#     """
#     input:
#         paramsfiles = COPIED_PARAMS
#     output:
#         EXPRESSION_FILTERED, \
#         MOTIF_PRIOR_FILTERED, \
#         PPI_PRIOR_FILTERED, \
#         STATS, \
#         PANDA_NET
#     container:
#         PYTHON_CONTAINER
#     shell:
#         """
#         echo ${input.paramsfiles}
#         sisana generate ${input.paramsfiles}
#         """