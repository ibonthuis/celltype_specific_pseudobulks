### Loading libraries
required_libraries <- c(
    "data.table",    
    "dplyr",
    "Seurat",
#   "limma",
 #  "rlang",
  # "purrr"
    "optparse"
  )

for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
}

### Options
options(stringsAsFactors = FALSE)

### Command line options
option_list <- list(

    optparse::make_option(
        c("-s", "--seuratobject"),
        type = "character",
        default = NULL,
        help = "Path to the seurat file in RData format.",
        metavar = "character"),
    # optparse::make_option(
    #     c("-m", "--metadata"),
    #     type = "character",
    #     default = NULL,
    #     help = "Path to the metadata file.",
    #     metavar = "character"),
    optparse::make_option(
        c("-c", "--celltypes"),
        type = "character",
        default = NULL,
        help = "A list of cell types. In the format of",
        metavar = "character"),
    optparse::make_option(
        c("-o", "--output_dir"),
        type = "character",
        default = NULL,
        help = "Path to the output directory.",
        metavar = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

## Initialize variable
SEURAT <- opt$seuratobject
OUTPUT_DIR <- opt$output_dir
CELLTYPES <- opt$celltypes

## Debug

## Functions
source("bin/create_count_matrix_fn.R")

## Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


## Computing
ALL <- get(load(SEURAT))
count_matrix <- LayerData(ALL, assay = "RNA", layer = "counts")
metadata <- ALL@meta.data
list_of_celltypes <- as.list(CELLTYPES)

# Filter count_matrix for cell types of interest
filtered_counts <- create_cell_counts_of_cell_types(metadata, list_of_celltypes, count_matrix)

# Filter count_matrix for genes expressed in more than 1% of cells

counts_perTE_filtered <- filter_counts_by_nonzero_geneentries(filtered_counts, 1)

# Maybe subsampling now
