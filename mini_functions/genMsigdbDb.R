genMsigdbDb <- function(msig_lvls, species){
  require(msigdbr)
  
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
  
  main_obj <- lapply(setNames(names(msig_lvls),names(msig_lvls)), function(mlvl){
    sub_obj <- lapply(setNames(msig_lvls[[mlvl]],msig_lvls[[mlvl]]), function(sublvl){
      print(paste0(">", mlvl, ":", sublvl, "..."))
      msig_ds <- msigdbr(species = species, category = mlvl, subcategory = sublvl) %>% 
        as.data.frame %>%
        dplyr::select(gs_name, gene_symbol) %>% 
        unique
      msig_tbl <- msig_ds %>%
        group_by(gs_name) %>%
        dplyr::summarise(Gene=paste0("", paste(gene_symbol, collapse=","))) %>% 
        mutate(gs_name2=gs_name) %>%
        relocate(gs_name2, .after=gs_name) %>%
        tibble::column_to_rownames('gs_name')
        
      return(msig_tbl)
    })
    return(sub_obj)
  })
  return(main_obj)
}

# msig <- genMsigdbDb(species='Mus musculus')
# msig <- unlist(msig, recursive=F)
# for(i in names(msig)){
#   dir <- '/Users/rquevedo/gsea_home/db/msigdbr/Mus_musculus/gene_symbol'
#   write.table(msig[[i]], file=file.path(dir, paste0(i, ".gmt")),
#               sep=",", col.names = F, row.names = T, quote = F)
# }
# 
# 
