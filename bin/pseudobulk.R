### Loading libraries
required_libraries <- c(
    "data.table",    
    "dplyr",
    "tidyr",
    "stringr",
   # "ggplot2",
   # "Seurat",
#   "limma",
 #   "rlang",
   "purrr",
    "optparse"#,
  #  "RColorBrewer"
  )

for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
}

### Options
options(stringsAsFactors = FALSE)

### Command line options
option_list <- list(

    optparse::make_option(
        c("-s", "--subsamples_list"),
        type = "character",
        default = NULL,
        help = "Path to the seurat file in RData format, containing the subsamples.",
        metavar = "character"),
    # optparse::make_option(
    #     c("-m", "--metadata"),
    #     type = "character",
    #     default = NULL,
    #     help = "Path to the metadata file.",
    #     metavar = "character"),
    # optparse::make_option(
    #     c("-c", "--celltypes"),
    #     type = "character",
    #     default = NULL,
    #     help = "A list of cell types.",
    #     metavar = "character"),
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
SUBSAMPLES_LIST <- opt$subsamples_list
METADATA <- opt$metadata
#CELLTYPES <- opt$celltypes
OUTPUT_PATH <- opt$output_dir


## Debug
#setwd("/storage/kuijjerarea/ine/projects/BRCA_SINGLECELL_TO_PSEUDOBULK/celltype_specific_pseudobulks")
#SUBSAMPLES_LIST <- "/storage/kuijjerarea/ine/projects/BRCA_SINGLECELL_TO_PSEUDOBULK/celltype_specific_pseudobulks/snakemake_results/subsampled/subsamples.RData"

#OUTPUT_PATH

## Functions
source("bin/pseudobulk_fn.R")

subsamples <- get(load(SUBSAMPLES_LIST))
#metadata <- fread(METADATA)
CT <- str_split_i(OUTPUT_PATH, "/", 4)
print(paste0(CT))
length(subsamples)

new_metadata <- data.frame(patient_name = c(names(subsamples)))
new_metadata <- new_metadata %>%
  tidyr::separate(patient_name, into = c("patientID", "tumor_type", "celltype"), sep = "_", remove = FALSE)

# celltype_specific_count_matrices <- gather_celltypes(celltype_ind = "T", count_matrices_list = subsampled_cm_, metadata = new_metadata)


# different_celltypes <- unique(new_metadata$celltype)
# pseudo_per_celltype <- purrr::map(different_celltypes, ~ gather_and_sum(.x, subsamples, new_metadata))
# pseudo_per_celltype_dfs <- purrr::map(pseudo_per_celltype, function(inner_list){
#    reduce(inner_list, function(x, y) full_join(x, y, by = "gene"))
# })

# pseudo_per_celltype_dfs <- purrr::map(pseudo_per_celltype_dfs, function(df) {
#   df <- df[, c(2, 1, 3:ncol(df))]
# })

# # For a single celltype:
pseudo_per_celltype <- gather_and_sum(CT, subsamples, new_metadata)
pseudo_per_celltype_df <- reduce(pseudo_per_celltype, function(x, y) full_join(x, y, by = "gene"))
pseudo_per_celltype_df <- pseudo_per_celltype_df[, c(2, 1, 3:ncol(pseudo_per_celltype_df))]

## End for a single celltype

# save(
#     pseudo_per_celltype_df,
#     file = file.path(OUTPUT_PATH, "pseudobulk.RData") # It should be pseudobulks per cell type of a single subsampling run in one file
# )

fwrite(
  pseudo_per_celltype_df, 
  file = file.path(OUTPUT_PATH, "pseudobulk.tsv"),
  sep = "\t",
  col.names = TRUE#,
#  row.names = TRUE
)

fwrite(
  new_metadata, 
  file = file.path(OUTPUT_PATH, "pseudobulk_metadata.tsv"),
  sep = "\t",
  col.names = TRUE#,
#  row.names = TRUE
)

# Then for the separate subsample runs, I have separate files, which can function be run in parallel
