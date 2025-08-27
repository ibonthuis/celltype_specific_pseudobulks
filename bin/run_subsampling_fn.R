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

# Modified the above function to take the index of the metadata
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

