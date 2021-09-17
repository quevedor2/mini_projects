library(dplyr)
lapply(c('unmapped', 'mapped'), function(m){
  df <- read.csv(paste0(m, ".blast_results.txt"), header=FALSE)
  df$species <- gsub("_chromosom.*", "", df$V1) %>% gsub("^.*?_", "", .)
  cnt.df <- as.matrix(sort(table(df$species), decreasing=TRUE))
  write.table(cnt.df, file=paste0(m,".counts.txt"), col.names=FALSE, row.names=TRUE, quote=FALSE)
})
