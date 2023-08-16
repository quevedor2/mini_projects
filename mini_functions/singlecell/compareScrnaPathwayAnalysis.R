## OBJECTIVE: Compares gene-set enrichment workflows for scRNA ##
library(AUCell)
library(Seurat)
library(ggplot2)
library(cowplot)
library(scales)
library(msigdbr)
# library(scRNAseq)
library(SingleR)
library(SCPA)
# library(testSctpa)
library(tidyverse)
# library(EnsDb.Hsapiens.v86)

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/ipmn_pdac'
datadir <- file.path(PDIR, "data")
outdir <- file.path(PDIR, "results")

###################
#### Functions ####
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")
source("~/git/mini_projects/mini_functions/geneMap.R")
source("~/.Rprofile")
gm <- geneMap('Homo sapiens')

aucellFun <- function(msig_ds, expr_mat, gm, mapfrom='SYMBOL', mapto='ENTREZID'){
  mapped_id <- if(mapfrom==mapto){
    msig_ds$entrez_gene
  } else {
    gm[[mapfrom]][[mapto]][as.character(msig_ds$entrez_gene)]
  }
  
  msig_l <- lapply(split(mapped_id, msig_ds$gs_name), function(i){
    i[!is.na(i)]
  })
  auc <- tryCatch({
    AUCell_run(expr_mat, msig_l)
  }, error=function(e){NULL})
  return(auc)
}

wilcox_auc <- function(auc, seu, grps=c("CD4+ T cells", "B cells"), colid='manual_anno'){
  if(class(auc) != 'Matrix'){
    auc <- auc
  } else {
    auc <- assay(auc)
  }
  spl_auc <- split(as.data.frame(t(auc[,colnames(seu)])),
                   seu@meta.data[,colid])
  spl_auc <- spl_auc[grps]
  
  wilcoxp <- sapply(1:ncol(spl_auc[[1]]), function(idx){
    tryCatch({
      wilcox.test(spl_auc[[1]][,idx], spl_auc[[2]][,idx])$p.value
    }, error=function(e){NA})
  })
  fc <- sapply(1:ncol(spl_auc[[1]]), function(idx){
    pc <- 0.0001
    tryCatch({
      log2(mean(spl_auc[[1]][,idx])+pc) - log2(mean(spl_auc[[2]][,idx])+pc)
    }, error=function(e){NA})
  })
  
  auc_df <- data.frame("Geneset"=colnames(spl_auc[[1]]),
                       'FC'=fc,
                       "p"=wilcoxp,
                       "padj"=p.adjust(wilcoxp, method='BH')) %>%
    mutate(Dataset=gsub("_.*", "", Geneset)) %>%
    arrange(padj)
  return(auc_df)
}

##############
#### Main ####
#---------------------------------
#--- Load and Parameter setup ----
for(i in c('scpa', 'aucell', 'ssgsea')) dir.create(file.path(outdir, "gse_benchmark", i), recursive = F, showWarnings = F)
scpa_f <- file.path(outdir, "gse_benchmark", "scpa", "scpa.cd4t_b.rds")
auc_f <- file.path(outdir, "gse_benchmark", "aucell", paste0("msigdb.", "LEVEL", ".rds"))
auc_manual_f <- file.path(outdir, "gse_benchmark", "aucell", paste0("msigdb_manual.", "LEVEL", ".rds"))
ssgsea_f <- file.path(outdir, "gse_benchmark", "ssgsea", paste0("ssgsea.", "LEVEL", ".rds"))

seu <- readRDS(file.path(outdir, "seurat_obj", "3_annotated.rds"))
set.seed(1234)
sidx <- sample(x = 1:30000, size=5000, replace=F)
seu <- subset(seu, cells=Cells(seu)[sidx])
DefaultAssay(seu) <- 'RNA'

# Annotate
bed.se <- readRDS("/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets/Monaco.rds") 
singler_anno <- SingleR(test=GetAssayData(seu, assay = "RNA", slot = "data"), 
                        ref=bed.se, 
                        assay.type.test=1, labels=bed.se[['label.main']],
                        de.method="wilcox", genes='sd',
                        clusters=seu$seurat_clusters)
cluster_ids <- setNames(singler_anno$pruned.labels, 
                        as.character(rownames(singler_anno)))
seu@meta.data[,'anno'] <- cluster_ids[as.character(seu$seurat_clusters)]


species <- 'Homo sapiens'
msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                  'C2'=list('CP:REACTOME', 'CP:KEGG'),  # curated gene sets
                  'C5'=list('GO:BP', 'GO:CC', 'GO:MF')) #, # ontology gene sets

# Read in Msigdbr Pathwys
pathways <- lapply(names(msig_lvls), function(mlvl){
  lapply(msig_lvls[[mlvl]], function(sublvl){
    msigdbr(species = species, category = mlvl, subcategory=sublvl) %>%
      format_pathways()
  })
}) %>%  unlist(., recursive=F) %>% do.call(c, .)

#----------=-----------------------------------------
#--- 0. a) Specify pathway data for benchmarking ----
DefaultAssay(seu) <- 'RNA'
expr <- GetAssayData(seu, slot='counts')
idx <- which(rowSums(expr==0) < (ncol(seu)*0.95))
mlvl <- 'C2'

stopifnot(mlvl == 'C2')
msig_ds <- msigdbr(species = 'Homo sapiens', category = mlvl, subcategory = 'REACTOME') %>%
  dplyr::select(gs_name, entrez_gene) %>%
  as.data.frame()
rand_gs <- unique(msig_ds$gs_name)[sample(1:length(unique(msig_ds$gs_name)), size=50,replace=F)]
# msig_ds <- msig_ds %>% 
#   filter(gs_name %in% c(scpa_top$Pathway,
#                         auc_top$Geneset,
#                         rand_gs))
gm <- geneMap('Homo sapiens')
mapped_id <- gm[['ENTREZID']][['SYMBOL']][as.character(msig_ds$entrez_gene)]
msig_l <- lapply(split(mapped_id, msig_ds$gs_name), function(i){
  i[!is.na(i)]
})

#--- 0. b) Run SCPA-baseline ----
scpa_out <- compare_seurat(seu,
                           group1 = "anno", 
                           group1_population = c("CD4+ T cells", "B cells"),
                           pathways = pathways, parallel=F)
scpa_out <- scpa_out %>%
  mutate(Dataset=gsub("_.*", "", Pathway))
scpa_out %>% filter(Dataset=='REACTOME') %>% filter(qval > quantile(qval, 0.95))
saveRDS(scpa_out, file=scpa_f)

#-------------------------#
#--- 0. c) Run AUCell ----
overwrite=TRUE
DefaultAssay(seu) <- 'RNA'
expr <- GetAssayData(seu, slot='counts')
idx <- which(rowSums(expr==0) < (ncol(seu)*0.95)) # gene expressed in >=5% of all cells
auc_ls <- lapply(names(msig_lvls[c(1:2)]), function(mlvl){
  print(paste0("Msigdb lvl: ", mlvl, "..."))
  outf <- gsub("LEVEL", mlvl, auc_f)
  if(!file.exists(outf) | overwrite){
    auc_l <- iterateMsigdb(species='Homo sapiens', msig_lvls=msig_lvls[mlvl], 
                           fun=aucellFun, expr_mat=expr[idx,], gm=gm,
                           mapfrom='ENTREZID', mapto='SYMBOL') %>% 
      unlist(., recursive=F)
    saveRDS(auc_l, file=outf)
  } else {
    auc_l <- readRDS(outf)
  }
  return(auc_l)
})

#--- 0. d) Run ssGSEA ----
ssgsea_score = GSVA::gsva(expr = expr[idx,], msig_l, 
                          method = "ssgsea", 
                          parallel.sz = 1, verbose = T)
saveRDS(ssgsea_score, file=gsub("LVL", mlvl, ssgsea_f))
# ssgsea_df <- wilcox_auc(ssgsea_score, seu, grps=c("CD4+ T cells", "B cells"))

#--- 0. e) Run Barkley-Enrichment ----
msig_l

#--- 0. f) VISION ----
#--- 0. g) PAGODA2 ----

#-------------------------------------------------
#--- 1. a) Load in precomputed pathway scores ----
# Load in SCPA 
scpa_out <- readRDS(file=scpa_f)
scpa_out <- scpa_out %>%
  mutate(Dataset=gsub("_.*", "", Pathway))
# poi <- head((scpa_out %>% filter(Dataset=='REACTOME'))$Pathway, 20)

# Load in AUCell
aucell_scores <- lapply(c('C2'), function(mlvl){
  unlist(readRDS(file=gsub("LEVEL", mlvl, auc_f)), recursive=F)
}) %>% 
  unlist(., recursive=F) %>%
  do.call(rbind, .)
# seu[['AUCell']] <- CreateAssayObject(data = assay(aucell_scores))

# Load in ssGSEA
ssgsea_score <- readRDS(file=gsub("LVL", mlvl, ssgsea_f))

res = GeneToEnrichment(srt=seu, db=msig_l[1:500], method='rand_distribution', addtoseu=FALSE)

barkley_df <- wilcox_auc(auc=as.data.frame(t(res$l2r)), seu, grps=c("CD4+ T cells", "B cells"), colid='anno') 
ssgsea_df <- wilcox_auc(auc=as.data.frame(ssgsea_score), seu, grps=c("CD4+ T cells", "B cells"), colid='anno') 
aucell_df <- wilcox_auc(auc=as.data.frame(assay(aucell_scores)), seu, grps=c("CD4+ T cells", "B cells"), colid='anno') 
head(barkley_df, 20)



ds <- 'REACTOME'
pcutoff <- 0.01
mkSigOrder <- function(x, adjp, geneset){
  x %>% 
    mutate(rank=order(x[,adjp])) %>% 
    filter(x[,adjp] < pcutoff) %>% 
    dplyr::select(!!rlang::sym(geneset), rank)
}
scpa_sig <- scpa_out %>%
  filter(Dataset==ds) %>%
  mkSigOrder(., adjp='adjPval', geneset='Pathway')


scpa_sig <- scpa_out %>%
  filter(Dataset==ds) %>%
  mkSigOrder(., adjp='adjPval', geneset='Pathway')

aggregate_df <- aucell_df %>%
  filter(Dataset=='REACTOME') %>% 
  mkSigOrder(., adjp='padj', geneset='Geneset') %>%
  dplyr::rename_with(., ~ paste0("aucell.", .), .cols = -Geneset) %>%
  dplyr::full_join(., ssgsea_df %>% mkSigOrder(., adjp='padj', geneset='Geneset') %>%
                     rename_with(., ~paste0("ssGSEA.", .), .cols = -Geneset), by='Geneset') %>%
  dplyr::full_join(., barkley_df %>% mkSigOrder(., adjp='padj', geneset='Geneset') %>%
                     rename_with(., ~paste0("Barkley.", .), .cols = -Geneset), by='Geneset') %>%
  dplyr::full_join(., scpa_sig %>% rename_with(., ~paste0("SCPA.", .), .cols = -Pathway), by=c('Geneset'='Pathway')) %>%
  tibble::column_to_rownames(., "Geneset")

head(aggregate_df, 30)
metric_df <- apply(aggregate_df, 2, function(i){
  apply(aggregate_df, 2, function(j){
    cor(i,j, method='spearman', use='complete.obs')
    # sum((!is.na(i)) & (!is.na(j))) / sum((!is.na(i)) | (!is.na(j)))
  })
})
metric_df[metric_df>0.5] <- 0.5




#--- 1. c) BEYOND ----
# #--- AUCell Manual
# ac_rankings = AUCell::AUCell_buildRankings(expr[idx,], nCores = 1, 
#                                            plotStats = FALSE, verbose = F)
# sc_AUC = AUCell::AUCell_calcAUC(msig_l, ac_rankings, normAUC = T, 
#                                 aucMaxRank = ceiling(0.05 * nrow(ac_rankings)), verbose = F)
# manual_auc = AUCell::getAUC(sc_AUC)
# saveRDS(manual_auc, file="~/xfer/manualauc_score.rds")
# manual_auc <- readRDS(file="~/xfer/manualauc_score.rds")
# manualauc_df <- wilcox_auc(manual_auc, seu, grps=c("CD4+ T cells", "B cells"))

#--- ssGSEA from GSVA package
ssgsea_score = GSVA::gsva(expr = expr[idx,], msig_l, 
                          method = "ssgsea", 
                          parallel.sz = 1, verbose = T)
saveRDS(ssgsea_score, file=gsub("LVL", mlvl, ssgsea_f))

#---- GSVA pathway enrichment
gsva_score = GSVA::gsva(expr[idx,], msig_l, 
                        method = "gsva", parallel.sz = 1, 
                        verbose = T)
saveRDS(gsva_score, file="~/xfer/gsva_score.rds")
gsva_score <- readRDS(file="~/xfer/gsva_score.rds")
gsva_df <- wilcox_auc(gsva_score, seu, grps=c("CD4+ T cells", "B cells"))


# Calculate differential of AUCell using wilcoxon
aucell_scores_sub <- split(as.data.frame(t(assay(aucell_scores)[,colnames(seu)])), 
                           seu$monaco.main.cluster)
aucell_scores_sub <- aucell_scores_sub[c("CD4+ T cells", "B cells")]
wilcoxp <- sapply(1:ncol(aucell_scores_sub[[1]]), function(idx){
  wilcox.test(aucell_scores_sub[[1]][,idx], aucell_scores_sub[[2]][,idx])$p.value
})
fc <- sapply(1:ncol(aucell_scores_sub[[1]]), function(idx){
  pc <- 0.0001
  log2(mean(aucell_scores_sub[[1]][,idx])+pc) - log2(mean(aucell_scores_sub[[2]][,idx])+pc)
})

auc_df <- data.frame("Geneset"=colnames(aucell_scores_sub[[1]]),
                     'FC'=fc,
                     "p"=wilcoxp,
                     "padj"=p.adjust(wilcoxp, method='BH')) %>%
  mutate(Dataset=gsub("_.*", "", Geneset)) %>%
  arrange(padj)



ds <- 'REACTOME'
n<-30
top_scpa <- scpa_out %>% filter(Dataset==ds) %>%
  left_join(auc_df %>% filter(Dataset==ds),
            by=c("Pathway"="Geneset", "Dataset")) %>%
  arrange(desc(qval)) %>% head(., n)
top_scpa
print("")
top_aucell <- scpa_out %>% filter(Dataset==ds) %>%
  left_join(auc_df %>% filter(Dataset==ds),
            by=c("Pathway"="Geneset", "Dataset")) %>%
  arrange(padj) %>% head(., n)
top_aucell

auc_top <- auc_df %>% filter(Dataset==ds) %>% head(., n) %>% dplyr::select(Geneset)
scpa_top <- scpa_out %>% filter(Dataset==ds) %>% head(., n) %>% dplyr::select(Pathway)
list("scpa"=setdiff(scpa_top$Pathway, auc_top$Geneset),
     "both"=intersect(scpa_top$Pathway, auc_top$Geneset),
     "aucell"=setdiff(auc_top$Geneset, scpa_top$Pathway))




msig_ds <- msigdbr(species = 'Homo sapiens', category = 'C2', subcategory = 'REACTOME') %>%
  dplyr::select(gs_name, entrez_gene) %>%
  as.data.frame()
rand_gs <- unique(msig_ds$gs_name)[sample(1:length(unique(msig_ds$gs_name)), size=50,replace=F)]
msig_ds <- msig_ds %>% 
  filter(gs_name %in% c(scpa_top$Pathway,
                        auc_top$Geneset,
                        rand_gs))
gm <- geneMap('Homo sapiens')
mapped_id <- gm[['ENTREZID']][['SYMBOL']][as.character(msig_ds$entrez_gene)]
msig_l <- lapply(split(mapped_id, msig_ds$gs_name), function(i){
  i[!is.na(i)]
})

DefaultAssay(seu) <- 'RNA'
expr <- GetAssayData(seu, slot='counts')
idx <- which(rowSums(expr==0) < (ncol(seu)*0.95))
seqs <- seq(0.01, 0.21, by=0.04)
#--- AUCell Wrapper
auc_dfs <- lapply(setNames(seqs, paste0("frac", seqs)), function(maxrank){
  print(maxrank)
  auc <- tryCatch({
    AUCell_run(exprMat=expr[idx,], 
               geneSets=msig_l, 
               aucMaxRank=(maxrank*nrow(expr[idx,])))
  }, error=function(e){NULL})
  auc_df <- wilcox_auc(auc, seu, grps=c("CD4+ T cells", "B cells"))
  return(auc_df)
})
saveRDS(auc_dfs, file="~/xfer/auc_dfs.rds")
auc_dfs <- readRDS(file="~/xfer/auc_dfs.rds")

#--- AUCell Manual
ac_rankings = AUCell::AUCell_buildRankings(expr[idx,], nCores = 1, 
                                           plotStats = FALSE, verbose = F)
sc_AUC = AUCell::AUCell_calcAUC(msig_l, ac_rankings, normAUC = T, 
                                aucMaxRank = ceiling(0.05 * nrow(ac_rankings)), verbose = F)
manual_auc = AUCell::getAUC(sc_AUC)
saveRDS(manual_auc, file="~/xfer/manualauc_score.rds")
manual_auc <- readRDS(file="~/xfer/manualauc_score.rds")
manualauc_df <- wilcox_auc(manual_auc, seu, grps=c("CD4+ T cells", "B cells"))

#--- ssGSEA from GSVA package
ssgsea_score = GSVA::gsva(expr = expr[idx,], msig_l, 
                          method = "ssgsea", 
                          parallel.sz = 1, verbose = T)
saveRDS(ssgsea_score, file="~/xfer/ssgsea_score.rds")
ssgsea_score <- readRDS(file="~/xfer/ssgsea_score.rds")
ssgsea_df <- wilcox_auc(ssgsea_score, seu, grps=c("CD4+ T cells", "B cells"))

#---- GSVA pathway enrichment
gsva_score = GSVA::gsva(expr[idx,], msig_l, 
                        method = "gsva", parallel.sz = 1, 
                        verbose = T)
saveRDS(gsva_score, file="~/xfer/gsva_score.rds")
gsva_score <- readRDS(file="~/xfer/gsva_score.rds")
gsva_df <- wilcox_auc(gsva_score, seu, grps=c("CD4+ T cells", "B cells"))


# testSctpa:::cal_vision
vis = VISION::Vision(counts, signatures = gSets_path, 
                     projection_method = "UMAP", sig_gene_threshold = 0)
print(gc())
options(mc.cores = n_cores)
vis = VISION::analyze(vis)


ssgsea_dfs <- lapply(setNames(seqs, paste0("frac", seqs)), function(maxrank){
  x = cal_PAS(seurat_object = seu,
              tool = 'AUCell',   ## GSVA, ssGSEA, plage, zscore or Vision
              normalize = 'log',
              species = 'human', 
              pathway='reactome')
})
wilcox_auc <- function(auc, seu, grps=c("CD4+ T cells", "B cells")){
  if(class(auc) != 'Matrix'){
    auc <- auc
  } else {
    auc <- assay(auc)
  }
  spl_auc <- split(as.data.frame(t(auc[,colnames(seu)])),
                   seu$monaco.main.cluster)
  spl_auc <- spl_auc[grps]
  
  wilcoxp <- sapply(1:ncol(spl_auc[[1]]), function(idx){
    wilcox.test(spl_auc[[1]][,idx], spl_auc[[2]][,idx])$p.value
  })
  fc <- sapply(1:ncol(spl_auc[[1]]), function(idx){
    pc <- 0.0001
    log2(mean(spl_auc[[1]][,idx])+pc) - log2(mean(spl_auc[[2]][,idx])+pc)
  })
  
  auc_df <- data.frame("Geneset"=colnames(spl_auc[[1]]),
                       'FC'=fc,
                       "p"=wilcoxp,
                       "padj"=p.adjust(wilcoxp, method='BH')) %>%
    mutate(Dataset=gsub("_.*", "", Geneset)) %>%
    arrange(padj)
  return(auc_df)
}

pcutoff <- 0.001
mkSigOrder <- function(x, adjp, geneset){
  x %>% 
    mutate(rank=order(x[,adjp])) %>% 
    filter(x[,adjp] < pcutoff) %>% 
    dplyr::select(!!rlang::sym(geneset), rank)
}
scpa_sig <- scpa_out %>%
  filter(Dataset==ds) %>%
  mkSigOrder(., adjp='adjPval', geneset='Pathway')

aggregate_df <- lapply(auc_dfs, function(i){
  i%>% filter(Dataset==ds) %>%
    mutate(rank=order(padj)) %>% 
    filter(padj < pcutoff) %>% 
    dplyr::select(Geneset,rank)
}) %>% 
  purrr::reduce(full_join, by='Geneset') %>%
  dplyr::rename_with(., ~ paste0("aucell.", names(auc_dfs)), .cols = -Geneset) %>%
  dplyr::full_join(., gsva_df %>% mkSigOrder(., adjp='padj', geneset='Geneset') %>%
                     rename_with(., ~paste0("GSVA.", .), .cols = -Geneset), by='Geneset') %>% 
  dplyr::full_join(., ssgsea_df %>% mkSigOrder(., adjp='padj', geneset='Geneset') %>%
                     rename_with(., ~paste0("ssGSEA.", .), .cols = -Geneset), by='Geneset') %>%
  dplyr::full_join(., manualauc_df %>% mkSigOrder(., adjp='padj', geneset='Geneset') %>%
                     rename_with(., ~paste0("mAUCell.", .), .cols = -Geneset), by='Geneset') %>%
  dplyr::full_join(., scpa_sig %>% rename_with(., ~paste0("SCPA.", .), .cols = -Pathway), by=c('Geneset'='Pathway')) %>%
  tibble::column_to_rownames(., "Geneset")

head(aggregate_df, 30)
metric_df <- apply(aggregate_df, 2, function(i){
  apply(aggregate_df, 2, function(j){
    cor(i,j, method='spearman', use='complete.obs')
    # sum((!is.na(i)) & (!is.na(j))) / sum((!is.na(i)) | (!is.na(j)))
  })
})
metric_df[metric_df>0.5] <- 0.5


pdf("~/xfer/pathway.pdf")
ggplot(reshape2::melt(metric_df), aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  viridis::scale_fill_viridis(discrete=FALSE) +
  theme_cowplot() + 
  theme(axis.text.x=element_text(angle = 90))
dev.off()



by=c("Pathway"="Geneset", "Dataset")

DefaultAssay(seu) <- 'AUCell'
pdf("~/xfer/x.pdf", width = 14, height = 14)
FeaturePlot(seu, feature=gsub("_", "-", head(top_scpa$Pathway)), reduction='umap', raster=T)
FeaturePlot(seu, feature=gsub("_", "-", head(top_aucell$Pathway)), reduction='umap', raster=T)
dev.off()
