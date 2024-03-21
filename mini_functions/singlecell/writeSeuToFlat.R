
writeSeuToFlat <- function(seurat_obj, 
                           reduction='umap', 
                           pca_reduction='pca',
                           assay='RNA', layer='counts',
                           out_metaf='metadata.csv',
                           out_cntsf='counts.mtx',
                           out_pcaf='pca.csv',
                           out_featuref='genes.csv'
                           ){
  ## Save metadata table
  red <- Embeddings(seurat_obj, reduction=reduction)
  seurat_obj$barcode <- colnames(seurat_obj)
  seurat_obj$UMAP_1 <- red[,1]
  seurat_obj$UMAP_2 <- red[,2]
  write.csv(seurat_obj@meta.data, file=out_metaf, quote=F, row.names=F)
  
  ## Write expression counts matrix
  counts_matrix <- GetAssayData(seurat_obj, assay=assay, layer=layer)
  Matrix::writeMM(counts_matrix, file=out_cntsf)
  
  ## Write dimensional reduction matrix
  write.csv(Embeddings(seurat_obj, reduction=pca_reduction), 
            file=out_pcaf, quote=F, row.names=F)

  ## Write feature names
  write.table(
    data.frame('gene'=rownames(counts_matrix)),
    file=out_featuref, sep=',',
    quote=F,row.names=F,col.names=F
  )
}
