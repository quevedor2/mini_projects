## Create a mapping between the cellranger_aggr file and 
## pre-existing seurat merged files

mapSmergeToCRaggr <- function(seu_smerge, path_to_craggr, reduction='umap', 
                              return_reduction=FALSE, return_metadata=FALSE){
  cat("\tThis function is designed to map the merge of multiple seurat objects
      from multiple samples to the merge of the same samples using cellranger
      aggr. The intention for this is to be able to map the barcodes between
      the two in order to transfer the umap and metadata over correctly
      for Loupe-visualization and inspection\n")
  
  ##---- Read in the cellranger aggr files
  mtx <- Read10X(data.dir = path_to_craggr, #file.path(PDIR, "data", 'cellranger_aggr', 'cd45', 'outs', 'count', 'filtered_feature_bc_matrix'), 
                 strip.suffix=TRUE)
  seu_craggr <- CreateSeuratObject(counts = mtx, project = 'aggr')
  
  ##---- Find the best matching barcode identifiers between the two seurat objects
  merge_barcodes <- split(gsub("^.*_", "", Cells(seu_smerge)), 
                          gsub("_[ACGT]*$", "", Cells(seu_smerge)))
  aggr_barcodes <- split(gsub("-[0-9]$", "", colnames(seu_craggr)), 
                         gsub("^.*-", "", colnames(seu_craggr)))
  merge_aggr_map_mat <- sapply(merge_barcodes, function(i){
    sum_map <- sapply(aggr_barcodes, function(j){
      sum(i %in% j)
    })
    if(max(sum_map)/min(sum_map) <= 20){ warning("Low matching-detected; please manually inspect the mapping matrix: obj$mapping_mat")}
    return(sum_map)
  }) 
  merge_aggr_map <-  apply(merge_aggr_map_mat, 2, which.max)
  
  
  ##---- Creating a mapping between the barcodes using the best-approximation
  # assumes seurat merged IDs are:   [ID]_[barcode]
  # assumes cellranger aggr IDs are:   [barcode]-[sample_idx]
  barcode_map <- Cells(seu_smerge)
  for(idx in seq_along(merge_aggr_map)){
    barcode_map <- gsub(paste0(names(merge_aggr_map)[idx], "_([ACGT]*)$"), 
                         paste0("\\1-", merge_aggr_map[idx]),
                        barcode_map)
  }
  names(barcode_map) <- Cells(seu_smerge)
  
  umap <- NULL
  metadata <- NULL
  if(return_reduction){
    umap <- cbind("Barcode" = barcode_map[Cells(seu_smerge)],
                  Embeddings(object = seu_smerge, 
                             reduction = reduction)) %>%
      as.data.frame
  } 
  
  if(return_metadata){
    metadata <- seu_smerge@meta.data %>%
      mutate(Barcode=barcode_map[rownames(.)]) %>%
      relocate(Barcode) 
  }

  
  return(list("mapping_mat"=merge_aggr_map_mat,
              "mapping"=barcode_map,
              "umap"=umap,
              "metadata"=metadata))
}
