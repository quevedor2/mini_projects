## Used in conjunction with WGCNA R package
## https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html#define-a-minimum-counts-cutoff
## Expanded function for more flexibility

make_module_heatmap <- function(module_name,
                                expression_mat = normalized_counts,
                                metadata_df = metadata,
                                gene_module_key_df = gene_module_key,
                                module_eigengenes_df = module_eigengenes,
                                col_list=NULL, meta_arrange_ord=NULL, 
                                genekey='gene', gene.cex=8) {
  # Create a summary heatmap of a given module.
  #
  # Args:
  # module_name: a character indicating what module should be plotted, e.g. "ME19"
  # expression_mat: The full gene expression matrix. Default is `normalized_counts`.
  # metadata_df: a data frame with Sample and time_point
  #              as columns. Default is `metadata`.
  # gene_module_key: a data.frame indicating what genes are a part of what modules. Default is `gene_module_key`.
  # module_eigengenes: a sample x eigengene data.frame with samples as row names. Default is `module_eigengenes`.
  # col_list:  a list object of group names in metadata and colors for each grouping:
  #            col_list <- list(celltype = c("C" = "#f1a340", "G" = "#998ec3"),
  #                             group= c('DAB'='#7fc97f', 'DMSO'='#beaed4'))
  # meta_arrange_ord: a vector of the ordering that the metadata should take for plotting
  #
  # Returns:
  # A heatmap of expression matrix for a module's genes, with a barplot of the
  # eigengene expression for that module.
  
  # Set up the module eigengene with its sample_ids (rownames)
  module_eigengene <- module_eigengenes_df %>%
    dplyr::select(all_of(module_name)) %>%
    tibble::rownames_to_column("Sample")
  
  # Set up column annotation from metadata
  col_annot_df <- metadata_df %>%
    dplyr::inner_join(module_eigengene, by = "Sample") %>%
    dplyr::arrange(list(meta_arrange_ord)) %>%
    tibble::column_to_rownames("Sample")
  
  # Create the ComplexHeatmap column annotation object
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    df=col_annot_df[,names(col_list)],
    module_eigengene = ComplexHeatmap::anno_barplot(dplyr::select(col_annot_df, module_name)),
    col = col_list
  )
  
  # Get a vector of the Ensembl gene IDs that correspond to this module
  module_genes <- gene_module_key_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(genekey)
  
  # Set up the gene expression data frame
  mod_mat <- expression_mat %>%
    t() %>%
    as.data.frame() %>%
    dplyr::filter(rownames(.) %in% module_genes) %>%
    dplyr::select(rownames(col_annot_df)) %>%
    as.matrix()
  
  # Normalize the gene expression values
  mod_mat <- mod_mat %>%
    t() %>%
    scale() %>%
    t()
  
  # Create a color function based on standardized scale
  color_func <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#67a9cf", "#f7f7f7", "#ef8a62")
  )
  
  # Plot on a heatmap
  heatmap <- ComplexHeatmap::Heatmap(mod_mat,
                                     name = module_name,
                                     col = color_func,
                                     bottom_annotation = col_annot,
                                     cluster_columns = FALSE,
                                     show_row_names = TRUE,
                                     show_column_names = FALSE,
                                     row_names_gp = grid::gpar(fontsize = gene.cex)
  )
  
  return(heatmap)
}


