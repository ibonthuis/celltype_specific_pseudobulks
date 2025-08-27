### Loading libraries
required_libraries <- c(
    "data.table",    
    "dplyr",
    "tidyr",
  #  "stringr",
   # "ggplot2",
   # "Seurat",
#   "limma",
 #   "rlang",
   "purrr",
    "optparse",
    "RColorBrewer"
  )

for (lib in required_libraries) {
  suppressPackageStartupMessages(library(lib, character.only = TRUE, quietly = TRUE))
}

### Options
options(stringsAsFactors = FALSE)

### Command line options
option_list <- list(

    optparse::make_option(
        c("-c", "--countmatrix"),
        type = "character",
        default = NULL,
        help = "Path to the seurat file in RData format.",
        metavar = "character"),
    optparse::make_option(
        c("-m", "--metadata"),
        type = "character",
        default = NULL,
        help = "Path to the metadata file.",
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
COUNTMATRIX <- opt$countmatrix
METADATA <- opt$metadata
OUTPUT_PATH <- opt$output_dir


## Debug
#COUNTMATRIX <- "/storage/kuijjerarea/ine/projects/BRCA_SINGLECELL_TO_PSEUDOBULK/celltype_specific_pseudobulks/test/countmatrix_subsample.RData"
COUNTMATRIX <- "/storage/kuijjerarea/ine/projects/BRCA_SINGLECELL_TO_PSEUDOBULK/celltype_specific_pseudobulks/snakemake_results/countmatrix_filtered/countmatrix_filtered.RData"
METADATA <- "/storage/kuijjerarea/ine/projects/BRCA_SINGLECELL_TO_PSEUDOBULK/celltype_specific_pseudobulks/snakemake_results/countmatrix_filtered/metadata_single_cells_filtered.tsv"
#METADATA <- "/storage/kuijjerarea/ine/projects/BRCA_SINGLECELL_TO_PSEUDOBULK/celltype_specific_pseudobulks/test/metadata_single_cells.tsv"
setwd("/storage/kuijjerarea/ine/projects/BRCA_SINGLECELL_TO_PSEUDOBULK/celltype_specific_pseudobulks")
OUTPUT_DIR <- "subsamples"

## Functions
#setwd("/storage/kuijjerarea/ine/projects/BRCA_SINGLECELL_TO_PSEUDOBULK/celltype_specific_pseudobulks")
source("bin/run_subsampling_fn.R")

count_matrix <- get(load(COUNTMATRIX))

# Do I want filtered metadata?
metadata <- fread(METADATA)

# head(metadata)
# head(colnames(count_matrix))
# nrow(metadata)
# Important to have the correct list structure for celltypes_list:
celltypes_list <- as.list(unique(metadata$celltype)) 

metadata_per_ct <- map(celltypes_list, ~ meta_per_ct(metadata = metadata, celltype_list = .x))

counts_per_ct <- filter_count_matrix(count_matrix = count_matrix, metadata_per_ct = metadata_per_ct, celltype_list = celltypes_list)
# ncol(counts_per_ct) # 97894
# counts_per_ct[1:5, 1:5]


patients <- unique(metadata_per_ct[[1]]$patient_ID)
tumor_type <- unique(metadata_per_ct[[1]]$GroupID)

data <- list(
  patients = patients,
  tumor_type = tumor_type,
  metadata_per_ct = metadata_per_ct
)

# Create a combination of all patients and tumor types
combos <- tidyr::expand_grid(
  patient = patients,
  tumor_type = tumor_type,
  cell_type_index = seq_along(metadata_per_ct) # This will allow me to reference the correct metadata
)


# Apply the function using pmap
results <- pmap(combos, ~ meta_per_patient_and_tt(..1, ..2, ..3))


metadata_filtered <- map(results, ~ filter_metadata_100(.x))
metadata_filtered_list <- Filter(Negate(is.null), metadata_filtered)
subsampled_cm <- map(metadata_filtered_list, ~ make_subsample_count_matrix(filtered_metadata_ct_pp = .x, counts = counts_per_ct, nr_of_cells_to_sample = 100, OUTPUT_DIR = OUTPUT_DIR))

#subsampled_cm_ <- unlist(subsampled_cm)
flat_subsampled_cm <- purrr::flatten(subsampled_cm)
length(flat_subsampled_cm)

save(
  flat_subsampled_cm,
  file = file.path(OUTPUT_PATH, "subsamples.RData")
)


# Quick check if now orig.ident are kind of correct in the colnames of the subsamples
# col <- colnames(subsampled_cm[[1]])
# class(col)

# for(i in 1:length(col)) {
#   print(paste0(stringr::str_split_1(col[i], "_" )))
# }

# unique(stringr::str_split_1(, "_" ))
#unique(metadata_ct_pp[[3]]$celltype)




# When doing analysis on just one celltype: 
#metadata_per_ct <- meta_per_ct(metadata = metadata, celltype_list = celltypes_list)
#metadata_ct_pp <- map(patients, ~ meta_per_patient(.x, metadata_ct = metadata_per_ct))



# TO DO I need to check the order of the column names and metadata names
# Maybe they don't need to be in the same order, but they should have overlap, so that the right columns are selected for the count_matrix

