#### Cluster level proportion ####
.assignClusterToPartition <- function(seu, clusid, partid, annoid=NULL){
  # clusid <- 'seurat_clusters'
  # partid <- 'monocle3_partitions'
  tbl <- table(seu@meta.data[,clusid], seu@meta.data[,partid])
  partitions <- apply(tbl / rowSums(tbl), 1, which.max)
  partitions <- setNames(colnames(tbl)[partitions], names(partitions))
  
  prt_clus_df <- as.data.frame(partitions) %>% 
    tibble::rownames_to_column(., "clusters") %>%
    arrange(partitions)
  
  if(!is.null(annoid)){
    tbl <- table(seu@meta.data[,clusid], seu@meta.data[,annoid])
    annos <- apply(tbl / rowSums(tbl), 1, which.max)
    annos <- setNames(colnames(tbl)[annos], names(annos))
    
    prt_clus_df <- prt_clus_df %>% 
      left_join(as.data.frame(annos) %>% 
                  tibble::rownames_to_column(., "clusters"),
                by='clusters')
  }
  return(prt_clus_df)
} 

# Gets the proportion of cells in cluster [clusid] that belong to either 
# smaple groupings of [sampleid] based on the conditions set by
# [ident_meta] and [dvar].  For example, proportion of cells
# per seurat_clusters belonging to samples A,B or C,D, where ident_meta
# indicates that A,B are GroupX and C,D are group Y
.getDeltaTbl <- function(seu, clusid, sampleid, dvar, ident_meta, method='stdres'){
  sort_ord <- unique(c(dvar, colnames(ident_meta), 'sampleID'))
  ident_spl <- ident_meta %>% 
    dplyr::arrange(!!! rlang::syms(sort_ord)) %>% 
    group_split(!! rlang::sym(dvar))
  
  # sub_ident_meta <- ident_spl[[1]][,sort_ord[c(2,3,1)]]
  sub_samples <- sapply(ident_spl, function(i) i%>% pull('sampleID'))
  
  cnt_tbl <- table(list(seu@meta.data[,clusid], seu@meta.data[,sampleid]))
  chi_tbl <- chisq.test(cnt_tbl)
  tbl <- apply(cnt_tbl, 2, function(i) i/sum(i)) %>%
    round(., 3)
  
  
  cluster_hc <- hclust(dist(tbl))$order
  cluster_hc <- rownames(tbl)[cluster_hc]
  dtbl <- switch(method,
                 rr=log2(tbl[,sub_samples[1]] / tbl[,sub_samples[2]]),
                 dist=tbl[,sub_samples[1]] - tbl[,sub_samples[2]],
                 stdres=chi_tbl$stdres[,1])
  
  list("dtbl"=dtbl, "cluster_ord"=cluster_hc, 'tbl'=tbl,
       'cnt_tbl'=cnt_tbl, 'chi_tbl'=chi_tbl)
}

# Wrapper function for annotation and getting the proportion of cells
# in cluster [clusid] that belong to either smaple groupings of [sampleid] based 
# on the conditions set by [ident_meta] and [dvar].  For example, 
# proportion of cells per seurat_clusters belonging to samples A,B or C,D, 
# where ident_meta indicates that A,B are GroupX and C,D are group Y
getDiffProportionsPerClusters <- function(
  obj, clusid='seurat_clusters', sampleid='orig.ident', 
  delta, ident_meta, annoid='functional.cluster', partid='functional.cluster'
){
  # Creates a mapping of cluster to annotations
  cluster_meta <- .assignClusterToPartition(seu, clusid,partid,annoid)
  
  # Get the proportion of GroupA to GroupB per cluster and
  # format to add in the annotations
  dtbl_l <- .getDeltaTbl(obj, clusid, sampleid, delta, ident_meta, method='stdres')
  
  dtbl_melt <- melt(dtbl_l$dtbl) %>%
    rename_with(., ~c('Delta')) %>%
    tibble::rownames_to_column(., "Cluster") %>%
    mutate(Cluster=factor(as.character(Cluster)), 
           Direction=(Delta>0) %>%
             if_else(., "Pos","Neg")) %>%
    left_join(., cluster_meta, by=c('Cluster'='clusters'))
  dtbl_melt$Cluster <- factor(dtbl_melt$Cluster, levels=dtbl_l$cluster_ord)
  return(dtbl_melt)
}

conditions <- c('Condition')
dtbl_melt <- getDiffProportionsPerClusters(seu, clusid='seurat_clusters', 
                                           sampleid='orig.ident', delta='Condition',
                                           ident_meta, annoid=annoid,
                                           partid=partid) %>%
  select(-c(partitions, annos))
  



#### Nearest-neighbour Proportions ####
reduction <- 'mnn'
category_label <- 'orig.ident'

seu <- seul$LN
## Calculate the nearest neighbour for k for the latent-space
# reductions
labels <- factor(seu@meta.data[,category_label])
lspace_obj <- seu@reductions[[reduction]]@cell.embeddings
ks <- seq(5, 50, by=5)
k_wproportions <- lapply(setNames(ks, ks), function(k){
  nn_obj <- RANN::nn2(data = lspace_obj, query = lspace_obj, 
                      k=k, treetype='kd')
  wdist <- round(1-(nn_obj$nn.dist/max(nn_obj$nn.dist)), 2) * 100
  ulbl <- unique(labels)
  wproportions <- sapply(c(1:nrow(lspace_obj)), function(idx){
    wcnt <- split(wdist[idx,-1], labels[nn_obj$nn.idx[idx,-1]]) %>%
      sapply(., sum)
    wcnt / sum(wcnt)
  })
  t(wproportions)
})
w <- (ks/max(ks))
wproportion <- sapply(1:ncol(k_wproportions[[1]]), function(colidx){
  sapply(k_wproportions, function(i){
    i[,colidx]
  }) %>%
    apply(., 1, weighted.mean, w=w)
}) %>%
  magrittr::set_colnames(paste0("lprop_", colnames(k_wproportions[[1]])))

## weigh the 
wproportions
seu$localProportion <- wproportion