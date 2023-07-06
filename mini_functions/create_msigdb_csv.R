library(msigdbr)
library(SCPA)
msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                  'C2'=list('CP:REACTOME', 'CP:KEGG'),  # curated gene sets
                  'C5'=list('GO:BP', 'GO:CC', 'GO:MF')) #, # ontology gene sets

species <- 'Mus musculus'
maxncol <- 5000
for(mainlvl in names(msig_lvls)){
  for(sublvl in msig_lvls[[mainlvl]]){
    mdb <- msigdbr(species = species, category = mainlvl, subcategory=sublvl) %>%
      format_pathways()
    names(mdb) <- sapply(mdb, function(i) i$Pathway[1])
    mdb_df <- lapply(mdb, function(i) as.data.frame(t(i$Genes))) %>%
      do.call(plyr::rbind.fill, .) %>% 
      as.data.frame %>% 
      mutate(Pathway=names(mdb)) %>%
      dplyr::relocate(., Pathway)
    if(ncol(mdb_df) > maxncol) stop("maxncol is smaller than the largest geneset, increase maxncol")
    mdb_df <- plyr::rbind.fill(as.data.frame(matrix(nrow=0, ncol=maxncol)) %>%
                                 mutate(Pathway="") %>%
                                 dplyr::relocate(., Pathway),
                               mdb_df)
    # mdb_df <- data.frame("Pathway"=names(mdb),
    #                      "Genes"=sapply(mdb, function(i) paste(i$Genes, collapse=",")))
    
    write.table(mdb_df, file=
                  paste0(
                    paste(c(tolower(mainlvl), gsub(":", "", tolower(sublvl))), collapse="_"),
                    ".", gsub(" ", "_", tolower(species)), 
                    ".msigdb.csv"), 
                sep=",", col.names = F, row.names = F, quote = F, na='')
  }
}
