# This script is meant to reduce the Uniref90 mapping in the biobakery3-humann
# directory to all other pathways (e.g. kegg, go, eggnog, ...)

library(dplyr)
library(data.table)
PDIR <- '/cluster/projects/mcgahalab/ref/metagenomics/biobakery_workflows_databases/humann/utility_mapping/rds'
setwd(PDIR)

method <- 'names'
method <- 'pathways'

files <- switch(method,
                names=list.files(PDIR, pattern="name.txt$"),
                pathways=list.files(PDIR, pattern="uniref90.txt$"))

util_map <- lapply(files, function(f){
  print(f)
  # maxcnt <- max(count.fields(f, sep = '\t'))
  dat = readLines(f)
  dat2 <- strsplit(dat, split="\\\t")
  names(dat2) <- sapply(dat2, function(i) i[[1]])
  return(dat2)
})
names(util_map) <- gsub("^.*_(.*)_.*", "\\1", files)
saveRDS(util_map, file=paste0(method, ".rds"))

