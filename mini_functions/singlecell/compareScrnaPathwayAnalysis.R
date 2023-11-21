renv::load("/cluster/home/quever/downloads/renvs/")
## OBJECTIVE: Compares gene-set enrichment workflows for scRNA ##
library(AUCell)
library(Seurat)
library(ggplot2)
library(cowplot)
library(scales)
library(msigdbr)
# library(scRNAseq)
library(pagoda2)
library(SingleR)
library(SCPA)
# library(testSctpa)
library(tidyverse)
# library(EnsDb.Hsapiens.v86)

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/benchmark/pathway_activity'
outdir <- file.path(PDIR, "results")
datadir <- file.path(PDIR, "data")

dir.create(outdir, showWarnings = F, recursive = T)
# IPMN_dir <- '/cluster/projects/mcgahalab/data/mcgahalab/ipmn_pdac'


###################
#### Functions ####
source("~/.Rprofile"); wideScreen()
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")
source("~/git/mini_projects/mini_functions/geneMap.R")
source("~/git/mini_projects/mini_functions/singlecell/pagoda2.R")
gml <- list("Homo sapiens"=geneMap('Homo sapiens'),
            "Mus musculus"=geneMap('Mus musculus'))

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
    mean
    tryCatch({
      log2((mean(spl_auc[[1]][,idx], na.rm=T)) / (mean(spl_auc[[2]][,idx], na.rm=T)+pc))
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

wilcoxDifferentialPathway <- function(seu, group1, group1_population, assay){
  samples <- list()
  for (i in group1_population) {
    samples[[i]] <- SCPA::seurat_extract(seu, assay = assay, 
                                         meta1 = group1, value_meta1 = i)
  }
  avg_expression <- lapply(samples, function(x) data.frame(rowMeans(x)))
  samp_combined <- cbind(avg_expression[[1]], avg_expression[[2]]) %>% 
    magrittr::set_colnames(c("Pop1", "Pop2")) %>%
    mutate(log2FC = log2(Pop1 / Pop2))
  samp_combined$p <- sapply(seq_along(rownames(samples[[1]])), function(idx){
    wilcox.test(samples[[1]][idx,], samples[[2]][idx,])$p.value
  })
  samp_combined$p.adj <- p.adjust(samp_combined$p, method='BH')
  samp_combined <- samp_combined %>% arrange(p.adj)
  return(samp_combined)
}


#############>> Main <<###############
#### 0. Load and prepare data ####
#--- a) Load and Parameter setup ----
dir.create(file.path(outdir, "gse_benchmark"), showWarnings = F)
scpa_f <- file.path(outdir, "gse_benchmark", paste0("scpa.", "LEVEL", ".rds"))
auc_f <- file.path(outdir, "gse_benchmark", paste0("aucell.", "LEVEL", ".rds"))
ssgsea_f <- file.path(outdir, "gse_benchmark", paste0("ssgsea.", "LEVEL", ".rds"))
pagoda2_f <- file.path(outdir, "gse_benchmark", paste0("pagoda2.", "LEVEL", ".rds"))

seu_files <- list('ST2'=list('file'='st2_mouse_seu.rds', 
                          'species'='Mus musculus',
                          'group1'='manual_anno',
                          'group1_population'=c('B', 'CD4_Tcell')),
                  'IPMN'=list('file'='ipmnPDAC_human_seu.rds',
                           'species'='Homo sapiens',
                           'group1'='hpca.main.cluster',
                           'group1_population'=c("T_cells", "B_cell")),
                  'P14'=list('file'='P14_seurs_with_geneset_enrichment_scores.rds',
                          'species'='Mus musculus',
                          'group1'='cell_type',
                          'group1_population'=c('Proliferative P14', 'Terminal Tex')))
seul <- lapply(seu_files, function(seuf){
  seu <- readRDS(file.path(datadir, seuf$file))
  if(is.list(seu)) seu <- seu$ST2
  set.seed(1234)
  sidx <- sample(x = 1:ncol(seu), size=5000, replace=F)
  seu <- subset(seu, cells=Cells(seu)[sidx])
  DefaultAssay(seu) <- 'RNA'
  return(seu)
})
seulids <- setNames(names(seul), names(seul))

species <- c('Homo sapiens', 'Mus musculus') %>%
  setNames(., .)
msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                  'C2'=list('CP:KEGG'),  # curated gene sets
                  'C5'=list('GO:BP')) #, # ontology gene sets
mlvl <- 'C5'
msig_lvls <- msig_lvls[mlvl]

# Read in Msigdbr Pathways
pathways <- lapply(species, function(species_x){
  print(species_x)
  pathway2 <- lapply(names(msig_lvls), function(mlvl){
    lapply(msig_lvls[[mlvl]], function(sublvl){
      x <- msigdbr(species = species_x, category = mlvl, subcategory=sublvl)
      
      msig_ds <- x %>% dplyr::select(gs_name, entrez_gene) %>%
        as.data.frame
      msig_l <- split(gml[[species_x]]$ENTREZ$SYMBOL[msig_ds$entrez_gene],
                      msig_ds$gs_name)
      pathway <- x %>% format_pathways()
      dat <- list("msig_ds"=msig_ds,
                  "msig_l"=msig_l,
                  "pathway"=pathway)
      names(dat$pathway) <- sapply(dat$pathway, function(i) unique(i$Pathway))
      
      return(dat)
    })
  }) 

  .get <- function(pw, pwid){
    lapply(pw, function(i){
      lapply(i, function(j) j[[pwid]])
    }) %>% unlist(., recursive=F)
  }
  
  return(list("msig_ds"=.get(pathway2, 'msig_ds')  %>% do.call(rbind,.),
              "msig_l"=.get(pathway2, 'msig_l')  %>% do.call(c,.),
              "pathway"=.get(pathway2, 'pathway')  %>% do.call(c,.)))
})

#--- b) Specify pathway data for benchmarking ----
run_old <- FALSE
overwrite <- FALSE

seul_dat <- lapply(seul, function(seu){
  DefaultAssay(seu) <- 'RNA'
  expr <- GetAssayData(seu, slot='counts')
  idx <- which(rowSums(expr==0) < (ncol(seu)*0.95))
  return(list('seu'=seu,
              'expr'=expr,
              'idx'=idx))
})

# Create a subset to work on
set.seed(1234)
pathway_samples <- sapply(pathways, function(pathway){
  names(pathway$pathway)
}) %>% 
  unlist %>% unique %>%
  sample(., size=1000, replace=F)


pathways_sub <- lapply(pathways, function(pathway){
  list("pathway"=pathway$pathway[pathway_samples],
       "msig_l"=pathway$msig_l[pathway_samples],
       "msig_ds"=pathway$msig_ds %>% filter(gs_name %in% pathway_samples))
})

if(run_old){
  mlvl <- 'C5'
  
  # stopifnot(mlvl == 'C2')
  msig_ds <- msigdbr(species = 'Homo sapiens', category = mlvl, subcategory = 'GO:BP') %>%
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
}


###########################################
#### 1. Run pathway activity scorers  #####
#--- a) Run SCPA-baseline ----
outf <- gsub("LEVEL", paste0("sub", mlvl), scpa_f)
if(file.exists(outf)) {
  scpas <- readRDS(file=outf)
} else {
  scpas <- lapply(seulids, function(id){
    print(paste0("SCPA: ", id, "..."))
    seu <- seul_dat[[id]]$seu
    
    seu[['RNA2']] <- CreateAssayObject(data=seu@assays$RNA$data)
    DefaultAssay(seu) <- 'RNA2'
    scpa_out <- compare_seurat(seu,
                               assay='RNA2',
                               group1 = seu_files[[id]]$group1, 
                               group1_population = seu_files[[id]]$group1_population,
                               pathways = pathways_sub[[seu_files[[id]]$species]]$pathway, 
                               parallel=F)
    scpa_out <- scpa_out %>%
      mutate(Dataset=gsub("_.*", "", Pathway))
    return(scpa_out)
  })
  saveRDS(scpas, file=outf)
}




#--- b) Run AUCell ----
outf <- gsub("LEVEL", paste0("sub", mlvl), auc_f)

if(file.exists(outf)) {
  auc_res <- readRDS(file=outf)
} else {
  auc_res <- lapply(seulids, function(id){
    print(paste0("AUCell: ", id, "..."))
    seu <- seul_dat[[id]]$seu
    expr <- seul_dat[[id]]$expr
    idx <- seul_dat[[id]]$idx
    
    aucl <- aucellFun(msig_ds = pathways_sub[[seu_files[[id]]$species]]$msig_ds, 
                      expr_mat=expr[idx,], gml[[seu_files[[id]]$species]],
                      mapfrom='ENTREZID', mapto='SYMBOL')
    
    seu[['aucell']] <- CreateAssayObject(data=as(assay(aucl), 'dgCMatrix'))
    auc_res <- wilcoxDifferentialPathway(seu, group1=seu_files[[id]]$group1, 
                                         group1_population = seu_files[[id]]$group1_population,
                                         assay='aucell')
    list("score"=aucl,
         "res"=auc_res)
  })  
  saveRDS(auc_res, file=outf)
}

#--- c) Run ssGSEA ----
outf <- gsub("LEVEL", paste0("sub", mlvl), ssgsea_f)

if(file.exists(outf)) {
  ssgsea_res <- readRDS(file=outf)
} else {
  ssgsea_res <- lapply(seulids, function(id){
    print(paste0("ssGSEA: ", id, "..."))
    seu <- seul_dat[[id]]$seu
    expr <- seul_dat[[id]]$expr
    idx <- seul_dat[[id]]$idx
    
    ssgsea_score = GSVA::gsva(expr = expr[idx,], 
                              pathways[[seu_files[[id]]$species]]$msig_l, 
                              method = "ssgsea", 
                              parallel.sz = 4, verbose = T)
    seu[['ssgsea']] <- CreateAssayObject(data=as(ssgsea_score, 'dgCMatrix'))
    ssgsea_res <- wilcoxDifferentialPathway(seu, group1=seu_files[[id]]$group1, 
                                            group1_population = seu_files[[id]]$group1_population,
                                            assay='ssgsea')
    list("score"=ssgsea_score,
         "res"=ssgsea_res)
  })
  saveRDS(ssgsea_res, file=outf)
}

#--- e) PAGODA2 ----
outf <- gsub("LEVEL", paste0("sub", mlvl), pagoda2_f)

if(file.exists(outf)) {
  pagoda2_res <- readRDS(file=outf)
} else {
  pagoda2_res <- lapply(seulids, function(id){
    print(paste0("Pagoda2: ", id, "..."))
    seu <- seul_dat[[id]]$seu
    expr <- seul_dat[[id]]$expr
    idx <- seul_dat[[id]]$idx
    
    msig_l <- pathways_sub[[seu_files[[id]]$species]]$msig_l
    msig_l <- lapply(msig_l, na.omit)
    
    while(is.character(scores)){
      scores <- tryCatch({
        cal_pagoda2(expr[idx,], gSets=msig_l,
                    z.score=0) #pathways[[seu_files[[id]]$species]]$msig_l)
      }, error=function(e){
        return(NULL)
      })
    }
    print(paste0("CLASS ------------------------- ", class(scores)))
    if(!is.character(scores)){
      seu[['pagoda2']] <- CreateAssayObject(data=as(scores, 'dgCMatrix'))
      res <- wilcoxDifferentialPathway(seu, group1=seu_files[[id]]$group1, 
                                       group1_population = seu_files[[id]]$group1_population,
                                       assay='pagoda2')
    } else {
      res <- NULL
    }
    
    
    list("score"=scores,
         "res"=res)
  })
  saveRDS(pagoda2_res, file=outf)
}




#--- f) VISION ----

lapply(scpas, head, n=20)
lapply(auc_res, function(i) head(i$res, 20))
lapply(ssgsea_res, function(i) head(i$res, 20))
lapply(pagoda2_res, function(i) head(i$res, 20))


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

# Load in Pagoda2
pagoda2_score <- readRDS(file=gsub("LVL", mlvl, pagoda2_f))

res = GeneToEnrichment(srt=seu, db=msig_l[1:500], method='rand_distribution', addtoseu=FALSE)

barkley_df <- wilcox_auc(auc=as.data.frame(t(res$l2r)), seu, grps=c("CD4+ T cells", "B cells"), colid='anno') 
ssgsea_df <- wilcox_auc(auc=as.data.frame(ssgsea_score), seu, grps=c("CD4+ T cells", "B cells"), colid='anno') 
aucell_df <- wilcox_auc(auc=as.data.frame(assay(aucell_scores)), seu, grps=c("CD4+ T cells", "B cells"), colid='anno')
pagoda2_df <- wilcox_auc(auc=as.data.frame(pagoda2_score), seu, grps=c("CD4+ T cells", "B cells"), colid='anno') 
head(pagoda2_df, 20)



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
  dplyr::full_join(., pagoda2_df %>% mkSigOrder(., adjp='padj', geneset='Geneset') %>%
                     rename_with(., ~paste0("Pagoda2.", .), .cols = -Geneset), by='Geneset') %>%
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

#--- 1. c) BEYOND - PAGODA2 ----
library(pagoda2)
cm <- countMatrix <- p2data::sample_BM1
p2.processed <- basicP2proc(countMatrix, n.cores=1, min.cells.per.gene=10, 
                            n.odgenes=2e3, get.largevis=FALSE, make.geneknn=FALSE)
ext.res <- extendedP2proc(p2.processed, organism = 'hs')



counts <- gene.vs.molecule.cell.filter(cm, min.cell.size=500)
counts <- counts[rowSums(counts)>=10, ]
#create pagoda2 object
rownames(counts) <- make.unique(rownames(counts))
r <- Pagoda2$new(counts, log.scale=TRUE, n.cores=1)
r$adjustVariance(plot=TRUE, gam.k=10)
r$calculatePcaReduction(nPcs=50, n.odgenes=3e3)

r$makeKnnGraph(k=40, type='PCA', center=TRUE, distance='cosine')
r$getKnnClusters(method=infomap.community, type='PCA')
M <- 30
r$getEmbedding(type='PCA', embeddingType = 'largeVis', M=M, perplexity=30, gamma=1/M)
r$getEmbedding(type='PCA', embeddingType='tSNE', perplexity=50, verbose=FALSE)

r$getKnnClusters(method=multilevel.community, type='PCA', name='multilevel')
r$getKnnClusters(method=walktrap.community, type='PCA', name='walktrap')

r$getDifferentialGenes(type='PCA', verbose=TRUE, clusterType='community')

suppressMessages(library(org.Hs.eg.db))
ids <- unlist(lapply(mget(colnames(r$counts), org.Hs.egALIAS2EG, ifnotfound=NA), function(x) x[1]))
rids <- names(ids)
names(rids) <- ids
# list all the ids per GO category
go.env <- list2env(eapply(org.Hs.egGO2ALLEGS,function(x) as.character(na.omit(rids[x]))))

r$testPathwayOverdispersion(go.env, verbose=TRUE, 
                            recalculate.pca=F, min.pathway.size=1)

hdea <- r$getHierarchicalDiffExpressionAspects(type='PCA', clusterName='community', z.threshold=3)
genesets <- hierDiffToGenesets(hdea)
termDescriptions <- Term(GOTERM[names(go.env)]) # saves a good minute or so compared to individual lookups

sn <- function(x) { names(x) <- x; x}  ## utility function
genesets.go <- lapply(sn(names(go.env)),function(x) {
  list(properties=list(locked=TRUE, genesetname=x, shortdescription=as.character(termDescriptions[x])), genes=c(go.env[[x]]))
})
## concatenate
genesets <- c(genesets, genesets.go)

deSets <- get.de.geneset(r, groups = r$clusters$PCA[['community']], prefix = 'de_')
## concatenate
genesets <- c(genesets, deSets)
r$makeGeneKnnGraph(n.cores = 1)




r <- Pagoda2$new(counts, log.scale=TRUE, n.cores=1)
r$adjustVariance(plot=TRUE, gam.k=10)
r$calculatePcaReduction(nPcs=50, n.odgenes=3e3)
go.list <- eapply(org.Hs.egGO2ALLEGS,function(x) as.character(na.omit(rids[x])))
go.env <- list2env(go.list[1:100])
pdf("~/xfer/pagoda2.pdf")
r$testPathwayOverdispersion(setenv = go.env, verbose = T, recalculate.pca = T,
                             plot = T, min.pathway.size = 3)
dev.off()

path_names = names(r$misc$pwpca)
score = matrix(NA,nrow=length(path_names),ncol=ncol(counts))
rownames(score) = path_names
colnames(score) = colnames(counts)
for(i in 1:length(r$misc$pwpca)){
  if(!is.null(r$misc$pwpca[[i]]$xp$score)){
    score[i,] = as.numeric(r$misc$pwpca[[i]]$xp$scores)
  }
}

return(score)
},error = function(e){
  print(e)
  return("error")
})


