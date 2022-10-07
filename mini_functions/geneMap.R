# Returns a mapping from key to value (e.g. SYMBOL -> ENSEMBL == map$SYMBOL$ENSEMBL)
geneMap <- function(species='Homo sapiens'){
  genome_gse <- switch(species,
                       "Homo sapiens"={
                         library(org.Hs.eg.db)
                         org.Hs.eg.db
                       },
                       "Mus musculus"={
                         library(org.Mm.eg.db)
                         org.Mm.eg.db
                       })
  db <- c('SYMBOL', 'ENSEMBL', 'ENTREZID')
  db <- setNames(db,db)
  
  gene_map <- lapply(db, function(db_i){
    lapply(db[-match(db_i, db)], function(db_j){
      txby <- keys(genome_gse, db_i)
      i2j_ids <- mapIds(genome_gse, keys=txby, column=db_j,
                        keytype=db_i, multiVals="first")
      i2j_ids
    })
  })
  
  return(gene_map)
}