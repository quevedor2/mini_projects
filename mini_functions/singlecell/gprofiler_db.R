## Create a gene-set signature list from gProfiler database
gprofiler_dir <- '/cluster/projects/mcgahalab/ref/gprofiler'
gmt <- GSA::GSA.read.gmt(file.path(gprofiler_dir, 'gprofiler_full_mmusculus.ENSG.gmt'))
gprof_ds <-setNames(gmt$genesets, 
                    paste0(unlist(gmt$geneset.names), "_", unlist(gmt$geneset.descriptions)))
gprof_ds <- gprof_ds[grep("^CORUM", names(gprof_ds), invert=T)]


msig_ds <- lapply(names(gprof_ds), function(sublvl){
  data.frame("gs_name"=sublvl,
             "entrez_gene"=gprof_ds[[sublvl]])
}) %>% do.call(rbind,.)
msig_ds$classification <- gsub(":.*", "", msig_ds$gs_name)
