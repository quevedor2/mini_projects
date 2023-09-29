# Used to load in the gprofiler database stored in a gmt flatfile
# and parse it into a list of ontology-databases that contain
# a lit of gene sets in a named vector format
# Params: gene_map is expecting a named vector of gene names, where the names
# of the vector correpond to the Ensembl IDs
loadGprofilerDb <- function(gprofiler_dir='/cluster/projects/mcgahalab/ref/gprofiler',
                            species='human', gene_map=NULL){
  species_gmt <- switch(species,
                        human='gprofiler_full_hsapiens.ENSG.gmt',
                        mouse='gprofiler_full_mmusculus.ENSG.gmt')
  
  gmt <- GSA::GSA.read.gmt(file.path(gprofiler_dir, species_gmt))
  gprof_ds <-setNames(gmt$genesets, 
                      paste0(unlist(gmt$geneset.names), "_", unlist(gmt$geneset.descriptions)))
  gprof_ds <- gprof_ds[grep("^CORUM", names(gprof_ds), invert=T)]
  
  
  msig_ds <- lapply(names(gprof_ds), function(sublvl){
    data.frame("gs_name"=sublvl,
               "entrez_gene"=gprof_ds[[sublvl]])
  }) %>% do.call(rbind,.)
  
  msig_ds$classification <- gsub(":.*", "", msig_ds$gs_name)
  if(!is.null(gene_map)){
    cat("Converting ENSEMBL IDs and parsing list")
    msig_l <- lapply(split(gene_map[msig_ds$entrez_gene], msig_ds$gs_name), function(i){
      i[!is.na(i)]
    })
  } else {
    cat("Parsing genesets containing ENSEMBL ids")
    msig_l <- lapply(split(msig_ds$entrez_gene, msig_ds$gs_name), function(i){
      i[!is.na(i)]
    })
  }
  
  msig_l <- split(msig_l, f=gsub(":.*", "", names(msig_l)))
  return(msig_l)
}