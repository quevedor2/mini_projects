##-- Preprocessing functions ----
preprocessAzimuth <- function(seu){
  seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")
  seu <- subset(seu, subset = nFeature_RNA > 200 & 
                  nFeature_RNA < 6000 & 
                  percent.mt < 25 & 
                  nCount_RNA > 400)
  DefaultAssay(seu) <- 'RNA'
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000,
                              verbose = FALSE)
  
  seu <- SCTransform(seu, assay = 'RNA', new.assay.name = 'SCT',
                     vars.to.regress = 'percent.mt')
  DefaultAssay(seu) <- 'SCT'
  seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
  seu
}

preprocessSeu <- function(seu, ncount_min=400, ncount_max=Inf,
                          nfeature_min=200, nfeature_max=6000,
                          mt_max=25, org='mouse', numpcs=50, getPCs=FALSE,
                          variable_features=2000, res=1.2){
  require(glmGamPoi)
  
  if(org=='mouse'){
    mt_pattern <- "^mt-"
    s.features <- stringr::str_to_title(cc.genes$s.genes)
    g2m.features <- stringr::str_to_title(cc.genes$g2m.genes)
  } else {
    mt_pattern <- "^MT-"
    s.features <- cc.genes$s.genes
    g2m.features <- cc.genes$g2m.genes
  }
  
  seu <- PercentageFeatureSet(seu, pattern = mt_pattern, col.name = "percent.mt")
  
  # Filter
  seu <- subset(seu, subset = 
                  nCount_RNA > ncount_min &
                  nCount_RNA < ncount_max &
                  nFeature_RNA > nfeature_min & 
                  nFeature_RNA < nfeature_max & 
                  percent.mt < mt_max)
  
  # Normalization and Scaling
  DefaultAssay(seu) <- 'RNA'
  seu <- NormalizeData(seu,
                       normalization.method = "LogNormalize") %>%
    FindVariableFeatures(., selection.method = "vst",
                         nfeatures = variable_features, verbose = FALSE) %>% 
    ScaleData(.)
  
  # CellCycle Scoring
  seu <- CellCycleScoring(seu, s.features = s.features, 
                          g2m.features = g2m.features,
                          assay='RNA', set.ident = TRUE)
  seu$CC.Difference <- seu$S.Score - seu$G2M.Score
  
  # Normalization using SCTransform-v2
  seu <- SCTransform(seu, assay = 'RNA', new.assay.name = 'SCT', vst.flavor='v2',
                     vars.to.regress = c('percent.mt', "CC.Difference"), 
                     conserve.memory = T) %>%
    RunPCA(., npcs = max(c(100, numpcs)), verbose = FALSE) 
  
  if(getPCs){
    DefaultAssay(seu) <- 'SCT'
    numpcs <- getNumOfPcs(seu)
  }
  
  # Generate SCT-based UMAP
  DefaultAssay(seu) <- 'SCT'
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:numpcs, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:numpcs, verbose = FALSE) %>%
    FindClusters(resolution = res, verbose = FALSE)
  
  seu <- RenameCells(seu, new.names=paste0(colnames(seu), "_", 
                                           unique(seu$orig.ident)))
  return(seu)
}

getNumOfPcs <- function(seu, min_var=0.9){
  sdev <- seu@reductions$pca@stdev
  percent_var <- round(sdev^2/sum(sdev^2),5)
  max_pcs <- which(cumsum(percent_var) <= min_var) %>% max
  return(max_pcs)
}

##-- Doublet functions ----
runBcds <- function(seu, bcds.ntop=3000, bcds.srat=1, n=25){
  # Binary classification based doublet scoring (bcds)
  # Creates artificial doublets and trains a classifier to predict them
  # based on the most variable genes features
  
  sce <- as.SingleCellExperiment(seu)
  # bcds.ntop <- 3000 # number of variable gene features to consider
  # bcds.srat <- 1    # ratio of nCells to doublets
  # retRes <- TRUE    # return the classifier
  # varImp <- TRUE    # return variable importance
  # nmax <- 25        # number of iterations to train model
  sce = scds::bcds(sce, verb = TRUE, ntop = bcds.ntop,
                   srat = bcds.srat, retRes=TRUE, varImp=TRUE)
  bcds_score <- as.data.frame(sce$bcds_score)
  top_n <- getTopScdsPairs(sce, method='bcds', n=n)
  return(list("bcds"=bcds_score, "top_pairs"=top_n))
}

runCxds <- function(seu, cxds.ntop=3000, cxds.binThresh=0, n=5){
  # Co-expression based doublet scoring
  # Binarizes gene-expression per cell (present/absent), then finds the rate of
  # co-expression between genes. Attributes a score to each gene based on how
  # often they coexpress (LESS often has a hgiher score). Then for each cell,
  # sum the co-expressed gene pair scores (high scores denotes doublets as they
  # will continue both high-scoring single-cell genes)
  
  sce <- as.SingleCellExperiment(seu)
  # cxds.ntop <- 3000   # number of variable gene features to consider 
  # cxds.binThresh <- 0 # minimum counts to consider a gene "present"
  # retRes <- TRUE      # whether to return gene pair scores & top-scoring gene pairs
  
  sce = scds::cxds(sce, retRes = TRUE, ntop =cxds.ntop, 
                   binThresh = cxds.binThresh)
  cxds_score <- data.frame("cxds"=sce$cxds_score,
                           "z_cxds"=scale(sce$cxds_score)[,1])
  top_n <- getTopScdsPairs(sce, method='cxds', n=n)
  return(list("cxds"=cxds_score, "top_pairs"=top_n))
}

getTopScdsPairs <- function(sce, method, n=5){
  rs   = rownames(sce)
  if(method=='bcds'){
    top_n = metadata(sce)$bcds_vimp[1:n,] %>%
      as.data.frame %>%
      mutate(gene=rs[gene_index])
  } else if(method=='cxds'){
    ho   = rowData(sce)$cxds_hvg_ordr[which(rowData(sce)$cxds_hvg_bool)]
    hgs  = rs[ho]
    top_n = metadata(sce)$cxds$topPairs[1:n,] %>%
      as.data.frame %>%
      mutate(V1=hgs[V1], V2=hgs[V2])
  }
  
  return(top_n)
}

getExpMultiplet <- function(cnt){
  multiplet_chart <- data.frame(
    "rate"=(c(0.4, 0.8, 1.6, 2.3, 3.1, 3.9, 4.6, 5.4, 6.1, 6.9, 7.6)/100),
    "cells_loaded"=(c(0.8, 1.6, 3.2, 4.8, 6.4, 8, 9.6, 11.2, 12.8, 14.4, 16)*1000),
    "cells_recovered"=(c(0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10) * 1000)
  )
  multiplet_lm <- lm(rate~cells_recovered, data=multiplet_chart)
  predict.lm(multiplet_lm, data.frame("cells_recovered"=cnt)) * cnt %>%
    round(., 4)
}


runDoubletFinderv3 <- function(seu, annotation_col, k=500){
  # pK identification
  pcs <- 30
  
  # how many cells to compare against for doublets
  pk <- round(k/ncol(seu),5)
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(seu@meta.data[,annotation_col])
  nExp_poi <- getExpMultiplet(ncol(seu))  # multiplet chart from https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/76
  nExp_poi.adj <- nExp_poi*(1-homotypic.prop)
  
  gen_sct <- ("SCT" %in% names(seu@assays))
  DefaultAssay(seu) <- ifelse(gen_sct, 'SCT', 'RNA')
  seu <- doubletFinder_v3(seu, PCs = 1:pcs, pN = 0.25, pK = pk,  # first bump under 0.1 from bcmvn_seu
                          nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = gen_sct)
  seu$df_doublets <- seu@meta.data[,grep('DF.classifications', colnames(seu@meta.data), value=T)]
  seu$df_pANN <- seu@meta.data[,grep('^pANN', colnames(seu@meta.data), value=T)]
  return(seu)
}

# Calculates the delta between the doublet score for each cell compared to the
# median doublet score across all cells, then calculates the mean MAD score 
# for each seurat cluster to identify potential clusters with "high" doublets
getDoubletDeltaZ <- function(x, median_bcds=NULL, median_cxds=NULL, 
                             median_df=NULL, mad_bcds=NULL, mad_cxds=NULL, 
                             mad_df=NULL){
  
  if(is.null(median_bcds)) median_bcds <- median(x$bcds_score)
  if(is.null(median_cxds)) median_cxds <- median(x$cxds_z)
  if(is.null(median_df)) median_df <- median(x$df_pANN)
  if(is.null(mad_bcds)) mad_bcds <- mad(x$bcds_score)
  if(is.null(mad_cxds)) mad_cxds <- mad(x$cxds_z)
  if(is.null(mad_df)) mad_df <- mad(x$df_pANN)
  
  # pdf("~/xfer/t.pdf")
  # VlnPlot(object = seu, features = c('bcds_score', 'cxds_z', 'df_pANN'),
  #         group.by = 'seurat_clusters', stack=T, combine=FALSE)
  # dev.off()
  x %>% 
    group_by(seurat_clusters) %>%
    dplyr::summarise(delta_bcds=median(bcds_score-median_bcds),
                     delta_cxds=median(cxds_z-median_cxds),
                     delta_df=median(df_pANN-median_df)) %>%
    mutate(z_bcds=(delta_bcds-median_bcds) / mad_bcds,
           z_cxds=(delta_cxds-median_cxds) / mad_cxds,
           z_df=(delta_df-median_df) / mad_df) %>%
    tibble::column_to_rownames(.,'seurat_clusters') %>%
    as.data.frame
}

# Uses the output of getDoubletDeltaZ() to select the top clusters that 
# have higher doublet scores than a given Z-score (MAD score)
getTopDub <- function(ds, zmin=2){
  zmat <- ds[,grep("^z", colnames(ds),value=T)]
  rownames(ds)[which(zmat>zmin, arr.ind = T)[,1]] %>% 
    unique
}

# This function will look at the doublet scores (cxds, bcds, doubletFinder) across
# all cells and find lcusters that are enriched in high-probability doublets. Then
# it will subcluster just that cluster to see whether the doublets are a sub-sample
# of that cluster and only flag high-confident doublets.
refineDoubletClusters <- function(seu, z1=2, z2=5){
  features <- c('bcds_score', 'cxds_z', 'df_pANN')
  
  # Get doublet score across all cells
  seu$dub_aggregate <- 0
  doublet_scores <- seu@meta.data[,c('seurat_clusters', features)]
  ds <- doublet_scores %>%
    getDoubletDeltaZ
  dub_clus <- getTopDub(ds, z1) # High (z>2) doublet probability clusters
  if(length(dub_clus) == 0){
    return(list("seu"=seu, "vlns"=NULL))
  }
  
  # Subset for high-doublet clusters, sub-cluster and detect doublet scores  
  Idents(seu) <- 'seurat_clusters'
  dub_subclus <- lapply(setNames(dub_clus,dub_clus), function(dub_clus_i){
    # dub_clus_i <- dub_clus[1]
    
    # Subset for individual clusters, and preprocess
    seu_dub_clus <- subset(x = seu, idents = dub_clus_i)   
    npc_val <- if(ncol(seu_dub_clus) > 50){
      seu_dub_clus <- RunPCA(seu_dub_clus, npcs = 50, verbose = FALSE)
      max_pcs <- getNumOfPcs(seu_dub_clus)
      seu_dub_clus <- FindNeighbors(object = seu_dub_clus, dims = 1:(max_pcs+1))
      seu2 <- FindClusters(object = seu_dub_clus, resolution = 1.4)
      
      # Get doublet score across all subsetted cells
      doublet_scores2 <- seu2@meta.data[,c('seurat_clusters', features)]
      ds2 <- doublet_scores2 %>%
        getDoubletDeltaZ(., 
                         median_bcds = median(seu$bcds_score), 
                         median_cxds = median(seu$cxds_z), 
                         median_df = median(seu$df_pANN),
                         mad_bcds = mad(seu$bcds_score), 
                         mad_cxds = mad(seu$cxds_z), 
                         mad_df = mad(seu$df_pANN))
      if(length(getTopDub(ds2,5)) > 0){
        barcodes <- colnames(seu2[,which(seu2$seurat_clusters %in% getTopDub(ds2,z2))])
        
        sub_vln <- VlnPlot(object = seu2, features = c('bcds_score', 'cxds_z', 'df_pANN'),
                           group.by = 'seurat_clusters', stack=T, combine=FALSE) +
          ggtitle(dub_clus_i)
      } else {
        barcodes = NULL; sub_vln=NULL
      }
    } else {
      barcodes=colnames(seu_dub_clus); sub_vln=NULL
    }
    
    return(list("barcodes"=barcodes, "sub_vln"=sub_vln))
  })
  barcodes <- as.character(unlist(sapply(dub_subclus, function(i) i$barcodes)))
  seu@meta.data[barcodes, 'dub_aggregate'] <- 1
  return(list("seu"=seu, "vlns"=lapply(dub_subclus, function(i) i$sub_vln)))
}


