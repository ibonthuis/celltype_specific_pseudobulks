# make_pseudobulk <- function(count_matrix){
#     samples <- unique(colnames(count_matrix))
#     list_of_pseudo_counts <- map(samples, ~sum_gene_expression(.x, count_matrix)) %>%
#         set_names(samples[samples %in% colnames(count_matrix)]) %>%
#         discard(is.null)
#     return(list_of_pseudo_counts)
# }

sum_gene_expression <- function( count_matrix){ 
    if (is.null(ncol(count_matrix))) {
       print("had less than 2 columns so not calculated")
       return(NULL)
    } else {
        pseudo_counts <- rowSums(count_matrix)
        return(pseudo_counts)
    }
}


gather_celltypes <- function(celltype_ind, count_matrices_list, metadata) {
    # This function gather_celltypes filters metadata, and then takes the cell identifiers of the individual celltype specified by celltype_ind 
    metadata <- metadata %>%
        dplyr::filter(celltype == celltype_ind) 
    print(paste0(unique(metadata$celltype)))
    
    count_matrices_ct <- count_matrices_list[names(count_matrices_list) %in% metadata$patient_name]
    return(count_matrices_ct)
}


gather_and_sum <- function(celltype_ind, count_matrices_list, metadata) {
    celltype_specific_count_matrices <- gather_celltypes(celltype_ind, count_matrices_list, metadata)
    # Put the rownames, which are the genenames, into an extra column so that they can later function as column to be merged with other pseudo_counts
    ct_spec_pseudos <- purrr::map(celltype_specific_count_matrices, ~ sum_gene_expression(.x)) 

ct_spec_pseudos <- purrr::imap(ct_spec_pseudos, function(x, name) {
    x <- as.data.frame(x)
    x$gene <- rownames(x)
    colnames(x) <- c(name, "gene")  # add the name as a new column
    return(x)
})
    return(ct_spec_pseudos)

}