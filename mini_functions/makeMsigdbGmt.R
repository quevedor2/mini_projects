## Mini script meant to take the R package msigdbr and create the 
# expected .gmt files used for enrichmentMap in cytoscape. There
# requires a little post-processing from this output to remove the
# NA values from the .gmt file.
#   Bash code: 
#   for i in $(ls *gmt); do sed -E "s/\tNA//g" ${i} > tmp; mv tmp ${i}; done
library(msigdbr)
outdir <- "~/gsea_home/db/msigdbr"
msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                  'C2'=list('CP:REACTOME'),             # curated gene sets
                  'C5'=list('GO:BP', 'GO:CC', 'GO:MF'), # ontology gene sets
                  'C8'=list(NULL))                      # cell type signature gene sets
species <- 'Homo sapiens'
columns <- c('entrez_gene', 'gene_symbol', 'ensembl_gene')
colid <- 'gene_symbol'

dir.create(file.path(outdir, gsub(" ", "_", species), colid), showWarnings = F)

.mapGsToExactSource <- function(msig_lvls, species){
  if(missing(msig_lvls)){
    warning("'msig_lvls' undefined: Defaulting to H, C2, C5, and C8")
    msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                      'C2'=list('CP:REACTOME'),             # curated gene sets
                      'C5'=list('GO:BP', 'GO:CC', 'GO:MF'), # ontology gene sets
                      'C8'=list(NULL)) 
  }
  if(missing(species)){
    warning("'species' undefined: Defaulting to Homo sapiens")
    species <- 'Homo sapiens'
  }
  
  main_obj <- lapply(names(msig_lvls), function(mlvl){
    sub_obj <- lapply(msig_lvls[[mlvl]], function(sublvl){
      print(paste0(">", mlvl, ":", sublvl, "..."))
      msig_ds <- msigdbr(species = species, category = mlvl, subcategory = sublvl) %>% 
        as.data.frame %>%
        select(gs_name, gs_exact_source) %>% unique
      with(msig_ds, setNames(gs_exact_source, gs_name))
    })
  }) %>% unlist %>% 
    gsub("[ :-?]", "_", .)
}
gs_map <- .mapGsToExactSource(species=species, msig_lvls=msig_lvls)


main_obj <- lapply(names(msig_lvls), function(mlvl){
  sub_obj <- lapply(msig_lvls[[mlvl]], function(sublvl){
    print(paste0(">", mlvl, ":", sublvl, "..."))
    msig_ds <- msigdbr(species = species, category = mlvl, subcategory = sublvl) %>%
      dplyr::select(gs_name, !!rlang::sym(colid)) %>%
      as.data.frame
    msig_df <- lapply(split(msig_ds, msig_ds[,1]), function(m_i){
      t(m_i[,2]) %>% 
        as.data.frame %>%
        mutate(ID=m_i[1,1,drop=T],
               DESCRIPTION=gs_map[m_i[1,1,drop=T]]) %>%
        relocate(ID,DESCRIPTION)
    }) %>% plyr::rbind.fill(.)
    return(msig_df)
  })
  return(sub_obj)
})
main_obj <- unlist(main_obj, recursive = F)
ids <- sapply(names(msig_lvls), function(i){
  paste0(i, ".", msig_lvls[[i]])
})  %>% unlist %>% as.character
names(main_obj) <- ids


for(id in names(main_obj)){
  write.table(main_obj[[id]], 
              file=file.path(outdir, gsub(" ", "_", species), colid, paste0(id, ".gmt")),
              sep="\t", col.names = F, row.names = F, quote = F)
}
