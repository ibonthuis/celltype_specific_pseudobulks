# Function from create_count_matrix_fn.R
create_cell_counts_of_cell_types <- function(metadata, celltypes_list, counts_per_cell) {
  # Filter the metadata for T and Epithelial cell types, excluding unpaired patients
  metadata <- metadata %>%
    filter(celltype %in% celltypes_list) %>%
    filter(patient_ID != "unpaired")
  
  # Get the identifiers for T and Epithelial cells
  idents_to_keep <- rownames(metadata)
  
  # Print the count of identified cells
  cat("; Number of cell identifiers for ", paste0(celltypes_list), length(idents_to_keep), "\n")
  
  # Filter counts_per_cell based on the identifiers
  counts_per_cell <- counts_per_cell[, colnames(counts_per_cell) %in% idents_to_keep]
  
  # Return the filtered counts
  return(counts_per_cell)
}


meta_per_ct <- function(metadata, celltype_list){
    metadata_ct <- metadata %>%
      dplyr::filter(celltype %in% celltype_list) %>%
      dplyr::filter(patient_ID != "unpaired")
    return(metadata_ct)
}

# meta_per_patient_and_tt <- function(patient, tumor_type, metadata_ct) { # counts_per_cell
#   filtered_metadata_ct_pp <- metadata_ct %>%
#     filter(patient_ID %in% patient) %>%
#     filter(GroupID %in% tumor_type)

#   # idents_to_keep <- rownames(metadata)

#   # counts_ct_pp <- counts_per_cell[, colnames(counts_per_cell) %in% idents_to_keep]
#   return(filtered_metadata_ct_pp)
# }


#' Extract Patient-Specific Metadata for a Given Tumor Type
#'
#' This function retrieves metadata for a specific patient and tumor type (being primary or metastatic (Lymph Node
#' in this case)), based on the provided cell type index. The metadata is filtered to include only the records
#' corresponding to the specified patient ID and tumor type.
#'
#' It is easiest to accompany this function with the expand_grid function from tidyr package.
#' So that one automatically makes all the patient, tumor type, cell type combinations.
#' @param patient A character string representing the unique identifier of the patient.
#' @param tumor_type A character string representing the tumor type of interest.
#' @param cell_type_index An integer representing the index for the cell type in the metadata list.
#'
#' @return A data frame containing the filtered metadata for the specified patient and tumor type.
#'         The data frame will include all relevant metadata fields specific to the treatment conditions.
#'
#' @details The function accesses a global metadata list, `metadata_per_ct`, which should
#'          be defined in the environment prior to calling this function. Each index in this list
#'          corresponds to a different cell type's metadata. The function ensures that the output
#'          contains only records that match both the provided patient ID and tumor type.
#'
#' @examples
#' combos <- tidyr::expand_grid(
#'  patient = patients,
#'  tumor_type = tumor_type,
#'  cell_type_index = seq_along(metadata_per_ct) 
#' )
#' Apply it with a function from purrr package:
#' results <- pmap(combos, ~ meta_per_patient_and_tt(..1, ..2, ..3))
#' 
#' Or with just one combination of patient and celltype and tumor type:
#' Assuming patient ID is "P001", tumor type is "Primary", and the cell type index is 1
#' filtered_metadata <- meta_per_patient_and_tt("P001", "Primary", 1)
#'
#' @export
meta_per_patient_and_tt <- function(patient, tumor_type, cell_type_index) {
  metadata_ct <- metadata_per_ct[[cell_type_index]]  # Get the appropriate metadata for this index

  filtered_metadata_ct_pp <- metadata_ct %>%
    filter(patient_ID == patient) %>%
    filter(GroupID == tumor_type)

  return(filtered_metadata_ct_pp)
}

filter_metadata_100 <- function(metadata_per_subset) {
  if (nrow(metadata_per_subset) < 100) {
     metadata_per_subset <- NULL
    print(paste0(unique(metadata_per_subset$patient_ID), " has less than 100 ", unique(metadata_per_subset$celltype), " cells" ))

   } 
   #else {
  #    print(paste0(unique(metadata_per_subset$patient_ID), " has more than 100 ", unique(metadata_per_subset$celltype), " cells" ))
  # }
  return(metadata_per_subset)
}


# I want to filter the count matrix for the celltypes as well, because then in the final function it doesn't load all the unnecessary cells like Fibroblasts
filter_count_matrix <- function(count_matrix, metadata_per_ct, celltype_list) {
  # first make the list concatenated again. In that way, you get counts of all cell types of interest in one matrix
  # Then later you need to split the count matrix per celltype again
  metadata_cts <- purrr::list_rbind(metadata_per_ct)
  idents_to_keep <- metadata_cts$V1
  count_matrix <- as.matrix(count_matrix)
  counts <- count_matrix[, colnames(count_matrix) %in% idents_to_keep] 
  return(counts)
}

# TO DO: check whether it is 100 cells and not 101 or something, because in the past I've seen stuff having 101 columns, but maybe the first column is just the gene column
make_subsample_count_matrix <- function(filtered_metadata_ct_pp, counts, nr_of_cells_to_sample, OUTPUT_DIR){
  # metadata must be already filtered for celltype
  idents_to_keep <- filtered_metadata_ct_pp$V1
  patient_id <- paste(unique(filtered_metadata_ct_pp$patient_ID), unique(filtered_metadata_ct_pp$GroupID), unique(filtered_metadata_ct_pp$celltype), sep = "_")
  print(patient_id)

 # file_path <- file.path(OUTPUT_DIR, paste0(name, ".RData"))
 # print(file_path)

 # print(paste0(idents_to_keep[1:5]))
  counts_ct_pp <- counts[, colnames(counts) %in% idents_to_keep] 

  set.seed(42)
  idx <- sample(c(1:ncol(counts_ct_pp)), nr_of_cells_to_sample)
#  print(paste0(length(idx)))

  sampled_df <- counts_ct_pp[, idx]
 # save(sampled_df, file = file_path)

  return(setNames(list(sampled_df), patient_id))
}

