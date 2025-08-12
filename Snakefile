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

