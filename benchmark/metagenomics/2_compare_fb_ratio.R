library(tidyverse)
analysis <- 'teresa' # 'benchmark'
if(analysis  == 'teresa'){
  pdir <- '/cluster/projects/mcgahalab/data/mcgahalab/teresa_metagenomics/shotgun/results'
  dir <- file.path(pdir, '/kraken/bracken/phylum')
  dir3 <- file.path(pdir, '/metaphlan/main')
  dir4 <- file.path(pdir, '/metaphlan4/main')
} else {
  pdir <- '/cluster/projects/mcgahalab/data/mcgahalab/benchmark/metagenomics/results'
  dir <- file.path(pdir, '/kraken/bracken/phylum')
  dir3 <- file.path(pdir, '/metaphlan/main')
  dir4 <- file.path(pdir, '/metaphlan4/main')
}


#### Metaphlan3 ####
metas <- lapply(list.files(dir3, pattern=".s1.tsv$"), function(f){
  read.table(file.path(dir3, f), sep="\t", header = F, 
             check.names = F, stringsAsFactors = F) %>%
    magrittr::set_colnames(., c('taxa', 'taxa_id', 'relabund', 'coverage', 'estimated_number_of_reads_from_the_clade'))
})
names(metas) <- gsub(".s1.tsv$", "", list.files(dir3, pattern=".s1.tsv$"))

meta3_df <- lapply(metas, function(m){
  m[grep("\\|p__[a-zA-Z0-9]*$", m$taxa),] %>% 
    dplyr::select(., c(taxa, relabund))
}) %>% 
  purrr::reduce(., full_join, by=c('taxa')) %>%
  tibble::column_to_rownames(., "taxa") %>%
  magrittr::set_colnames(names(metas)) %>%
  round(., 3)
m3ratio <- meta3_df['k__Bacteria|p__Firmicutes',] / meta3_df['k__Bacteria|p__Bacteroidetes',]
t(m3ratio)

#### Metaphlan4 ####
metas <- lapply(list.files(dir4, pattern=".s1.tsv$"), function(f){
  read.table(file.path(dir4, f), sep="\t", header = F, 
             check.names = F, stringsAsFactors = F) %>%
    magrittr::set_colnames(., c('taxa', 'taxa_id', 'relabund', 'coverage', 'estimated_number_of_reads_from_the_clade'))
})
names(metas) <- gsub(".s1.tsv$", "", list.files(dir4, pattern=".s1.tsv$"))

meta4_df <- lapply(metas, function(m){
  m[grep("\\|p__[a-zA-Z0-9]*$", m$taxa),] %>% 
    dplyr::select(., c(taxa, relabund))
}) %>% 
  purrr::reduce(., full_join, by=c('taxa')) %>%
  tibble::column_to_rownames(., "taxa") %>%
  magrittr::set_colnames(names(metas)) %>%
  round(., 3)

m4ratio <- meta4_df['k__Bacteria|p__Firmicutes',] / meta4_df['k__Bacteria|p__Bacteroidetes',]



#### Kraken ####
brackens <- lapply(list.files(dir, pattern=".bracken$"), function(f){
  read.table(file.path(dir, f), sep="\t", header = T, 
             check.names = F, stringsAsFactors = F)
})
names(brackens) <- gsub(".bracken$", "", list.files(dir, pattern=".bracken$"))

bracken_df <- lapply(brackens, function(b){
  b %>%
    select(name, fraction_total_reads)
}) %>% 
  purrr::reduce(., full_join, by=c('name')) %>%
  tibble::column_to_rownames(., "name") %>%
  magrittr::set_colnames(names(brackens)) %>%
  round(., 3)

ratio <- bracken_df['Bacillota',] / bracken_df['Bacteroidota',]


if(analysis  == 'teresa'){
  cbind(t(m3ratio), t(m4ratio), t(ratio)) %>%
    round(., 3) %>%
    magrittr::set_colnames(., c('meta3', 'meta4', 'kraken')) %>%
    as.data.frame %>%
    filter(grepl("PreOp", rownames(.))) %>%
    group_split(., gsub("_[0-9]*$", "", rownames(.)))
} else {
  cbind(t(m3ratio), t(m4ratio), t(ratio)) %>%
    round(., 3) %>%
    magrittr::set_colnames(., c('meta3', 'meta4', 'kraken')) %>%
    as.data.frame %>% 
    mutate("ID"=gsub("\\..*$", "", rownames(.)),
           "subset"=as.integer(gsub("^.*\\.(.*)K$", "\\1", rownames(.)))) %>%
    dplyr::arrange(subset) %>%
    group_split(., ID)
}