# 
# 
# seu=seu_i
# PCs = 1:PCNum; pN = 0.25; pK=pk
# nExp = nExp_poi.adj; sct = gen_sct
# annotations=as.character(seu_i$immgen_fine_cell); delim=','; grp.size=25

dfhelp_preprocess <- function(data_wdoublets, PCs, orig.commands, sct=TRUE){
  orig.commands <- seu@commands
  if (sct == FALSE) {
    print("Creating Seurat object...")
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
    print("Normalizing Seurat object...")
    seu_wdoublets <- NormalizeData(seu_wdoublets, 
                                   normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method, 
                                   scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor, 
                                   margin = orig.commands$NormalizeData.RNA@params$margin)
    print("Finding variable genes...")
    seu_wdoublets <- FindVariableFeatures(seu_wdoublets, 
                                          selection.method = orig.commands$FindVariableFeatures.RNA$selection.method, 
                                          loess.span = orig.commands$FindVariableFeatures.RNA$loess.span, 
                                          clip.max = orig.commands$FindVariableFeatures.RNA$clip.max, 
                                          mean.function = orig.commands$FindVariableFeatures.RNA$mean.function, 
                                          dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function, 
                                          num.bin = orig.commands$FindVariableFeatures.RNA$num.bin, 
                                          binning.method = orig.commands$FindVariableFeatures.RNA$binning.method, 
                                          nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures, 
                                          mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff, 
                                          dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
    print("Scaling data...")
    seu_wdoublets <- ScaleData(seu_wdoublets, 
                               features = orig.commands$ScaleData.RNA$features, 
                               model.use = orig.commands$ScaleData.RNA$model.use, 
                               do.scale = orig.commands$ScaleData.RNA$do.scale, 
                               do.center = orig.commands$ScaleData.RNA$do.center, 
                               scale.max = orig.commands$ScaleData.RNA$scale.max, 
                               block.size = orig.commands$ScaleData.RNA$block.size, 
                               min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
    print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                            npcs = length(PCs), 
                            rev.pca = orig.commands$RunPCA.RNA$rev.pca, 
                            weight.by.var = orig.commands$RunPCA.RNA$weight.by.var, 
                            verbose = FALSE)
    pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[,  PCs]
    
  }
  if (sct == TRUE) {
    require(sctransform)
    print("Creating Seurat object...")
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
    print("Running SCTransform...")
    seu_wdoublets <- SCTransform(seu_wdoublets)
    print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
    pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[,PCs]
  }
  cell.names <- rownames(seu_wdoublets@meta.data)
  nCells <- length(cell.names)
  rm(seu_wdoublets)
  gc()
  # seu_wdoublets <- FindNeighbors(object = seu_wdoublets, dims = 1:length(PCs))
  # seu_wdoublets <- FindClusters(object = seu_wdoublets, resolution = 1.2)
  # seu_wdoublets <- RunUMAP(seu_wdoublets, dims = 1:length(PCs))
  # pdf("~/xfer/test_wdub.pdf")
  # seu_wdoublets$dubs <- 'singlet'
  # seu_wdoublets@meta.data[c((n_real.cells+1):(n_real.cells+ncol(doublets))),]$dubs <- doublet_anno
  # DimPlot(seu_wdoublets, label = TRUE, repel = TRUE, reduction = "umap",
  #         group.by='dubs', pt.size=0.5, shuffle=TRUE)
  # dev.off()
  return(pca.coord)
}
dfhelp_annotations <- function(seu, annotations, delim, grp.size=NULL){
  stopifnot(typeof(annotations) == "character")
  stopifnot(length(annotations) == length(Cells(seu)))
  stopifnot(!any(is.na(annotations)))
  
  # Set up annotations
  real.cells <- Cells(seu)
  annotations <- factor(annotations)
  annotations_chr <- as.character(annotations)
  names(annotations_chr) <- names(annotations) <- real.cells
  
  # Find all unique celltype-celltype annotation [2-col data.frame]
  doublet_combn <- utils::combn(unique(annotations_chr), 2) %>%
    cbind(matrix(c(unique(annotations_chr),unique(annotations_chr)),
                 byrow = T, nrow=2)) %>%
    t %>% as.data.frame
  
  # Subset an equal amount of cells for each doublet
  real.cells12 <- apply(doublet_combn, 1, function(comb_i, grp.size){
    .SampleCt <- function(x, grp.size) sample(x, grp.size, replace=T)
    data.frame("cells1"=.SampleCt(names(which(annotations_chr == comb_i[1])), grp.size),
               "cells2"=.SampleCt(names(which(annotations_chr == comb_i[2])), grp.size))
  }, grp.size=grp.size) %>% do.call(rbind, .)
  
  # Track the selected cell-IDs
  real.cells1 <- real.cells12$cells1  # cell-id for doublet-cell1
  real.cells2 <- real.cells12$cells2  # cell-id for doublet-cell2
  doublet_types1 <- annotations[real.cells1] # actual cell-type for doublet-cell1
  doublet_types2 <- annotations[real.cells2] # actual cell-type for doublet-cell2
  doublet_anno <- paste0(doublet_types1, delim, doublet_types2) # joined doublet-id
  list("real1"=real.cells1, "real2"=real.cells2,
       "ct1"=doublet_types1, "ct2"=doublet_types2,
       "ct12"=doublet_anno, "annotations"=annotations,
       "uniq_ct"=doublet_combn)
}

.pANN <- function(neighbors, dubs){
  round(length(which(dubs))/k,4)
}
.ES <- function(neighbors,k, dubs){
  x <- (as.integer(dubs)*2) - 1
  sumx <- cumsum(x)
  sumx[sumx<0] <- 0
  round(max(sumx)/k,4)
}
.gsea <- function(neighbors, dist_i, dubs, alpha=1){
  # test gsea
  dubs <- neighbors[which(dubs)]
  Ra <- matrix(dist_i[neighbors]^alpha, ncol=1, dimnames = list(neighbors,"j"))
  if(length(dubs)==0){
    return(0)
  } else {
    GSVA:::.fastRndWalk(gSetIdx=as.character(dubs), 
                        geneRanking=as.character(neighbors), 
                        j='j', Ra)
  }
}

dfhelp_getpANN <- function(dist.mat, pK_i, cell.names, n_real.cells,
                           annotations=NULL, normalize=F, metric='GSEA'){
  k <- round(length(cell.names) * pK_i)
  cell.names.2 <- cell.names
  
  if(!is.null(annotations)){
    names(cell.names.2)[c(1:n_real.cells)] <- 'NA'
    dub_idx <- split(seq_along(cell.names.2), f=names(cell.names.2))
    dub_idx <- dub_idx[unique(names(cell.names.2))]
    dub_idx <- c(dub_idx, list("all"=c(n_real.cells+1):length(cell.names)))
  } else {
    dub_idx <- list("all"=c(n_real.cells+1):length(cell.names))
  }
  
  pANN_dub <- lapply(dub_idx, function(dub_idx_i){
    pANN <- apply(dist.mat[,1:n_real.cells], 2, function(dist_i){
      neighbors <- order(dist_i)
      neighbors <- neighbors[2:(k + 1)]
      dubs_boolean <- suppressWarnings(neighbors == dub_idx_i)
      
      switch(metric,
           'GSEA'=.gsea(neighbors, dist_i, dubs_boolean),
           'pANN'=.pANN(neighbors, dubs_boolean),
           'ES'=.ES(neighbors, k, dubs_boolean))
    })
    names(pANN) <- cell.names[c(1:n_real.cells)]
    pANN
  })
  if(normalize) {
    pANN_rng <- unlist(pANN_dub) %>%
      range(., na.rm=T)
    pANN_dub <- lapply(pANN_dub, function(i) i / (pANN_rng[2] - pANN_rng[1]))
    pANN_dub
  }
  pANN_dub
}

doubletFinder_v4 <- function (seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, 
                              sct = FALSE, annotations = NULL, delim=",", grp.size=10,
                              metric='GSEA'){
  print(grp.size)
  require(Seurat)
  require(fields)
  require(KernSmooth)
  if (reuse.pANN != FALSE) {
    pANN.old <- seu@meta.data[, reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    return(seu)
  }
  if (reuse.pANN == FALSE) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts[, real.cells]
    n_real.cells <- length(real.cells)
    
    if(is.null(annotations)){
      n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
      print(paste("Creating", n_doublets, "artificial doublets...", 
                  sep = " "))
      dbs_map <- list()
      dbs_map$real1 <- sample(real.cells, n_doublets, replace = TRUE)
      dbs_map$real2 <- sample(real.cells, n_doublets, replace = TRUE)
      dbs_map$ct12 <- paste("X", 1:n_doublets, sep = "")  #set doublet-IDs
    }  else {
      dbs_map <- dfhelp_annotations(seu = seu, annotations = annotations, 
                                    delim=delim, grp.size=grp.size)
      # dbs_map$real1/2 -> real.cells1/2
      # dbs_map$ct1/2 -> doublet_types1/2
      # dbs_map$ct12 -> doublet_anno
      # anno -> annotations
    }
    
    # Set up doublet data
    doublets <- (data[, dbs_map$real1] + data[, dbs_map$real2])/2  # simulate doublets
    colnames(doublets) <- paste("X", 1:ncol(doublets), sep = "")  #set doublet-IDs
    data_wdoublets <- cbind(data, doublets) # combine real and doublet
    cell.names <- setNames(colnames(data_wdoublets),
                           c(colnames(data), dbs_map$ct12)) # real+doublet IDs
    nCells <- length(cell.names)            # number of real+doublet cells
      
    orig.commands <- seu@commands
    # Preprocess and dimensional reduction of the data with doublets
    pca.coord <- dfhelp_preprocess(data_wdoublets, PCs = PCs, 
                                   orig.commands = orig.commands, sct = sct)
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    
    print("Computing pANN...")
    pk_l <- c(0.01, seq(0.05, 0.30, by=0.05), pK)
    
    pANN_l <- lapply(setNames(pk_l, pk_l), function(pK_i){
      dfhelp_getpANN(dist.mat, pK_i, cell.names, n_real.cells, 
                     annotations, normalize = T, metric=metric)
    })
    
    print("Getting all score for proportion of doublets around real cells...")
    pANN_df <- sapply(setNames(names(pANN_l),names(pANN_l)), function(pK_i){
      pANN_l[[pK_i]]$all
    }) %>% 
      as.data.frame
    pANN_df[is.na(pANN_df)] <- 0
    colnames(pANN_df) <- paste0("all_pK", colnames(pANN_df))
    colnames(pANN_df)[ncol(pANN_df)] <- 'pANN'
    
    if (!is.null(annotations)) {
      print(paste0("Computing neighborhood composition based on ", metric, "..."))
      # Find all unique celltype-celltype annotation
      neighbor_types <- do.call(cbind, pANN_l[[as.character(pK)]]) %>%
        as.data.frame
      neighbor_types[is.na(neighbor_types)] <- 0
    }
    
    print("Classifying doublets..")
    classifications <- rep("Singlet", n_real.cells)
    classifications[order(pANN_df$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, "DF.classifications"] <- classifications
    seu@meta.data <- cbind(seu@meta.data, pANN_df)
    
    if (!is.null(annotations)) {
        seu@meta.data <- cbind(seu@meta.data,
                               neighbor_types)
    }
    
    return(seu)
  }
}