#### Test Data ####
# seu=seu_i
# PCs = 1:PCNum; pN = 0.25; pK=pk
# nExp = nExp_poi.adj; sct = gen_sct
# annotations=as.character(seu_i$immgen_fine_cell); delim=','; grp.size=25

#### Functions for Seurat and expected doublets ####
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

getResForStableSmallClus <- function(seu_j, max.res=5, min.res=0.1, topn=5){
  # Identify the size of the clusters at each resolution
  res_seq <- seq(min.res, max.res, by=0.1)
  size_per_clus <- sapply(res_seq, function(res){
    x <- FindClusters(object = seu_j, resolution = res, verbose=F)
    table(x$seurat_clusters)
  })
  
  # Identify Pielou evenness index for resolution-based cluster_sizes
  # 1 = even [c(1,1,1)],   0 = uneven [c(1,0,0)]
  sdx <- sapply(size_per_clus, function(i) vegan::diversity(i) / log(sum(i))) %>% round(.,4)
  
  # Use the first-derivative of the evenness~res to idetnify when there are
  # big changes in the evennness of the clusters
  sdx_diff <- diff(sdx) %>% abs
  
  # Find peaks in the 1st derivative based on a 66% diff (1xSD)
  pks <- pracma::findpeaks(sdx_diff, npeaks=topn, threshold=sd(sdx_diff)*1)
  
  # pdf("~/xfer/test.pdf")
  # plot(sdx, type='l', ylim=c(0,0.5))
  # lines(x=c(2:length(sdx))-0.5, y=sdx_diff, col='red')
  # points(x=pks[,2]+0.5, y=pks[,1], col='black')
  # dev.off()
  
  # Return the resolution that it was calculated at
  res <- res_seq[pks[which.max(pks[,2]),2]+1]
  return(res)
}


#### DoubletFinder_v4 ####
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
                              metric='GSEA', annotation_col='immgen_fine_cell'){
  print(grp.size)
  require(Seurat)
  require(fields)
  require(KernSmooth)
  if (reuse.pANN) {
    pANN.old <- seu@meta.data[, reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    return(seu)
  } else {
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

#### DoubletFinder_v4 Viz ####

getCellCols <- function(x=NULL, overwrite=F){
  cell_cols <- list('B cells'=colorRampPalette(c("firebrick1","firebrick4")),
                    "Basophils"=colorRampPalette(c('bisque')),
                    "Eosinophils"=colorRampPalette(c('darkgoldenrod1')),
                    "Monocytes"=colorRampPalette(c('darkolivegreen1', 'darkolivegreen4')),
                    "NK cells"=colorRampPalette(c('darkseagreen1', 'darkseagreen4')),
                    "NKT"=colorRampPalette(c('aquamarine', 'aquamarine4')),
                    "Neutrophils"=colorRampPalette(c('pink', 'pink4')),
                    "Tgd"=colorRampPalette(c('cyan3')),
                    'T cells'=colorRampPalette(c("deepskyblue","dodgerblue4")), 
                    'DC'=colorRampPalette(c("darkorchid","darkorchid4")),
                    'Unknown'=colorRampPalette('grey'))
  if(!is.null(x)) {
    stopifnot(is.function(x[[1]]))
    cell_cols <- c(cell_cols, x)
  }
  return(cell_cols)
}

# Contour plot for the doublet enrichment score
#'
#' @param seu seurat object from doubletFinder_v4 with annotations
#' @param delim Delimiter used to separate doublet IDs in meta.data (default=,)
#' @param top_n The top heterotypic doublets to plot (default=23)
#' @param b The breakpoints for contour plotting (default=0.5)
#' @param sample_col (default=orig.ident)
#' @param cluster_col (default=seurat_clusters)
#' @param doublet_col (default=DF.classification)
#'
#' @return A list to be used for cowplot::plot_grid
#' @export
doubletFinder_contour <- function(seu, delim=",", top_n=23, b=0.5,
                                  sample_col='orig.ident', 
                                  cluster_col='seurat_clusters',
                                  doublet_col='DF.classifications'){
  # data validation
  stopifnot(any(grepl("umap", names(seu@reductions))))  # ensure UMAP was run
  stopifnot(any(grepl(delim, colnames(seu@meta.data))))
  
  # Calculate the mean ES for each cluster for each doublet celltype
  cluster_proportion <- .getClusterProp(seu, delim, cluster_col)
  
  # Isolate, order and select the top heterotypic doublets
  ct_pairings <- do.call(rbind, strsplit(colnames(cluster_proportion), split=delim))
  hetero_idx <- apply(ct_pairings,1, function(i) !any(duplicated(i))) # remove homotypic
  heteros <- apply(ct_pairings[hetero_idx,], 1, paste, collapse=",")
  top_heteros <- apply(cluster_proportion[,heteros], 2, max) %>% 
    sort %>% 
    tail(., top_n)
  top_heteros <- top_heteros[!grepl("Unknown", names(top_heteros))]
  
  # Create a blank UMAP grid for contour plotting
  umap_val <- as.data.frame(seu@reductions$umap@cell.embeddings) 
  umap_val <- round(umap_val /b) * b
  umap_grid <- expand.grid(seq(min(umap_val$UMAP_1), max(umap_val$UMAP_1), by=b),
                           seq(min(umap_val$UMAP_2), max(umap_val$UMAP_2), by=b)) %>%
    as.data.frame %>%
    setNames(., c('UMAP_1', 'UMAP_2'))
  
  # Seurat cluster dimplot, also showing doublet labels
  doublet_umap <- seu@reductions$umap@cell.embeddings %>%
    cbind(., data.frame('doublets'=seu@meta.data[,doublet_col])) %>%
    filter(doublets=='Doublet')
  dp_clus <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                     group.by=cluster_col, pt.size=0.5, shuffle=TRUE) 
  dp_clus_dub <- dp_clus +
    geom_point(data=doublet_umap, aes(x=UMAP_1, y=UMAP_2), 
               pch=21, col='black', fill='black') +
    ggtitle(unique(seu@meta.data[,sample_col]))
  
  
  # Iterate through every doublet and plot the contour map for these doublets
  mycolors<- function(x) {
    colors<-colorRampPalette(c("black", "yellow", "white"))( x )
    colors[1:x]
  }
  dp_contour <- lapply(names(top_heteros), function(hetero_i){
    umap_val$z <- seu@meta.data[,hetero_i]
    umap_val <- umap_val[which(umap_val$z != 0),]
    umap_grid <- left_join(umap_grid, umap_val, by=c('UMAP_1', 'UMAP_2'))
    umap_grid$z[is.na(umap_grid$z)] <- 0    # Remove NA, too low to compare
    umap_grid$z[with(umap_grid, z<0)] <- 0  # Negative means enrichment far away, so ignore

    mybreaks <- seq(0, 2.5, by=0.1)
    dp_clus + 
      ggtitle(hetero_i) +
      theme(legend.position='none',
            plot.title = element_text(size = 10, face = "bold")) +
      stat_contour_filled(data=umap_grid, aes(x=UMAP_1, y=UMAP_2,z=z), 
                          alpha=0.5, breaks=mybreaks) +
      scale_fill_manual(palette=mycolors, drop=FALSE)
  })
  
  plotlist <- c(list(dp_clus_dub + theme(legend.position='none')), dp_contour)
  return(plotlist)
}

.getClusterProp <- function(seu, delim=',', 
                            cluster_col='seurat_clusters'){
  # Calculate the mean ES for each cluster for each doublet celltype
  cts <- grep(delim, colnames(seu@meta.data), value=T) # all doublet subtypes
  cluster_proportion <- seu@meta.data %>%
    as.data.frame %>%
    group_by_at(cluster_col) %>%
    select(cts) %>%
    summarise_all(., mean, na.rm=T) %>%
    tibble::column_to_rownames(., cluster_col)
  return(cluster_proportion)
}

doubletFinder_heatmap <- function(seu, cell_cols_lbl, delim=',',
                                  cluster_col='seurat_clusters', top_n=23){
  # Calculate the mean ES for each cluster for each doublet celltype
  cluster_proportion <- .getClusterProp(seu, delim, cluster_col)
  annotation_colors <- list("V1"=cell_cols_lbl,
                            "V2"=cell_cols_lbl)
  
  # Isolate top X doublets and report mean ES per doublet-subtype and cluster
  cluster_proportion[cluster_proportion<0] <- 0
  cluster_proportion <- apply(cluster_proportion, 1, scales::rescale,to=c(0,1))  %>% t
  cluster_proportion <- cluster_proportion[,order(colSums(cluster_proportion), decreasing = T)]
  cluster_proportion <- cluster_proportion[,1:min(ncol(cluster_proportion),top_n)]
  cell_comb <- strsplit(colnames(cluster_proportion), split = ",") %>%
    do.call(rbind, .) %>%
    as.data.frame
  colnames(cluster_proportion) <- as.character(1:ncol(cluster_proportion))
  rownames(cell_comb) <- colnames(cluster_proportion)
  
  # Visualize the dimplot and the heatmap
  dp_clus <- DimPlot(seu, label = TRUE, repel = TRUE, reduction = "umap",
                     group.by=cluster_col, pt.size=0.5, shuffle=TRUE) 
  prop_ch <- ComplexHeatmap::pheatmap(cluster_proportion,
                                      annotation_col = cell_comb, cluster_rows = TRUE,
                                      annotation_colors=annotation_colors)
  
  heatmap_pl <- plot_grid(dp_clus, ggplotify::as.grob(prop_ch),
                          ncol=2, rel_widths = c(1,2.5))
  heatmap_pl
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











