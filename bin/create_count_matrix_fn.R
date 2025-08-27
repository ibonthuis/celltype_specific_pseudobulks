create_cell_counts_of_cell_types <- function(metadata, celltypes_list, counts_per_cell) {
  # Filter the metadata for T and Epithelial cell types, excluding unpaired patients
  # celltypes_list can be of class character, but the celltypes have to be separate
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


filter_counts_by_nonzero_geneentries <- function(count_matrix, percent_of_cells) {
    # count_matrix can be of datatype dataframe or matrix
    nr_of_cells <- ncol(count_matrix)*(percent_of_cells/100)
    min_nr_of_cells <- round(nr_of_cells)
    count_matrix <- count_matrix[rowSums(count_matrix != 0) > min_nr_of_cells,]
    cat( ";", paste0(nrow(count_matrix)), " number of genes expressed in ", paste0(percent_of_cells), "% of cells out of a total cell number ", paste0(ncol(count_matrix)), "\n")

    return(count_matrix)
}

meta_per_ct <- function(metadata, celltype_list){
    metadata_ct <- metadata %>%
      dplyr::filter(celltype %in% celltype_list) %>%
      dplyr::filter(patient_ID != "unpaired")
    return(metadata_ct)
}