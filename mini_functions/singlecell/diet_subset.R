
## Runs the seurat Subset function, as well as the essentials of DietSeurat
## in one step.  Less prone to error than subset since it looks for
## all cells matching to the given ident, and subsets based on that alone
diet_subset <- function(seu, ident, assay='RNA', slot='counts'){
  cells_to_keep <- WhichCells(seu, idents = ident)
  rna_counts_subset <- GetAssayData(seu, assay = assay, slot = slot)[, cells_to_keep]
  meta_data_subset <- seu@meta.data[cells_to_keep, ]
  seusub <- CreateSeuratObject(
    counts = rna_counts_subset,
    meta.data = meta_data_subset
  )
  return(seusub)
}