library(tidyverse)
analysis <- 'local' # 'hpc', 'local'
if(analysis  == 'teresa'){
  pdir <- '/cluster/projects/mcgahalab/data/mcgahalab/teresa_metagenomics/shotgun/results'
} else {
  pdir <- '/Users/rquevedo/Projects/mcgaha/benchmark/metagenomics'
}
kdir <- file.path(pdir, '/kraken/bracken')
m3dir <- file.path(pdir, '/metaphlan/main')
m4dir <- file.path(pdir, '/metaphlan4/main')

readData <- function(dir, pattern, labels, lvl_pattern, 
                     read.header=F, use.coverage=F){
  taxas <- lapply(list.files(dir, pattern=pattern), function(f){
    read.table(file.path(dir, f), sep="\t", header = read.header, 
               check.names = F, stringsAsFactors = F) %>%
      magrittr::set_colnames(., labels)
  })
  names(taxas) <- gsub(pattern, "", list.files(dir, pattern=pattern))
  
  taxa_df <- lapply(taxas, function(m){
    if(use.coverage){
      m[grep(lvl_pattern, m$taxa),] %>% 
        dplyr::select(., c(taxa, coverage))
    } else {
      m[grep(lvl_pattern, m$taxa),] %>% 
        dplyr::select(., c(taxa, relabund))
    }
  }) %>% 
    purrr::reduce(., full_join, by=c('taxa')) %>%
    tibble::column_to_rownames(., "taxa") %>%
    magrittr::set_colnames(names(taxas)) %>%
    round(., 3)
  return(taxa_df)
}
.getTaxa <- function(df, taxa_lvl){
  taxa_letter <- tolower(substr(taxa_lvl, 1, 1))
  taxa_pattern <- paste0('\\|', taxa_letter, '__[a-zA-Z0-9]*(_[a-zA-Z0-9]*)?$')
  df <- df %>% 
    filter(grepl(taxa_pattern, rownames(.)))
  rownames(df) <- tolower(gsub("^.*__", "", rownames(df))) %>%
    gsub("_", " ", .)
  df
}
.mapIds <- function(df, map.ids=NULL){
  if(is.null(map.ids)){
    map.ids <- c('bacillota'='firmicutes',
                 'bacteroidota'='bacteroidetes')
  }
  rownames(df) <- tolower(rownames(df))
  newids <- map.ids[rownames(df)]
  if(any(!is.na(newids))){
    rownames(df)[which(!is.na(newids))] <- na.omit(newids)
  }
  return(df)
}


m3df <- readData(m3dir, '.s1.tsv$', 
                 c('taxa', 'taxa_id', 'relabund', 'relcoverage', 
                   'coverage'),
                 ".*$",
                 use.coverage=T)
m4df <- readData(m4dir, '.s1.tsv$', 
                 c('taxa', 'taxa_id', 'relabund', 'relcoverage', 
                   'coverage'),
                 ".*$",
                 use.coverage=T)

taxa_lvl <- 'phylum'
top_taxa <- 10
k2df <- readData(file.path(kdir, taxa_lvl), ".bracken$", 
                 c('taxa', 'taxa_id', 'taxa_lvl', 'prevcov', 'addedcov', 
                   'coverage', 'relabund'), ".*",
                 use.coverage=T,
                 read.header=T) %>%
  .mapIds

m3df_lvl <- .getTaxa(m3df, taxa_lvl)
m4df_lvl <- .getTaxa(m4df, taxa_lvl)

cnt_l <- lapply(setNames(colnames(k2df),colnames(k2df)), function(id){
  .norm <- function(df){
    df[,1] <- round(df[,1] / sum(df[,1], na.rm=T), 3)
    return(df)
  }
  list(tibble::rownames_to_column(.norm(k2df[,id,drop=F])),
       tibble::rownames_to_column(.norm(m3df_lvl[,id,drop=F])),
       tibble::rownames_to_column(.norm(m4df_lvl[,id,drop=F]))) %>%
    purrr::reduce(., full_join, by='rowname') %>%
    tibble::column_to_rownames(., 'rowname') %>%
    magrittr::set_colnames(., c('kraken2', 'metaphlan3', 'metaphlan4'))
})

ggps <- lapply(names(cnt_l), function(id){
  cnt_i <- cnt_l[[id]]
  
  cnt_i <- cnt_i %>%
    filter(rownames(cnt_i) %in% (names(tail(sort(rowSums(cnt_i, na.rm=T)), top_taxa))))
  reshape2::melt(t(cnt_i)) %>%
    arrange(value) %>% 
    mutate(Var2=factor(Var2, levels=rev(unique(Var2))),
           ID=gsub("^(.*?)_.*?_(.*?\\..*)$", "\\1.\\2", id))
}) %>% do.call(rbind, .)

ggplot2::ggplot(ggps %>% 
                  filter(grepl("2500K", ID)),
                aes(x=Var2, y=value, col=Var1, group=Var1)) +
  geom_point() + 
  geom_line() + 
  facet_wrap(~ID, ncol=2) + 
  cowplot::theme_cowplot() +
  coord_flip() +
  theme(axis.text.x = element_text(angle=90)) + 
  xlab(taxa_lvl)
