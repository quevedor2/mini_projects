# # Read in gtf file to map ENS->Biotype
# gtf_file <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
# GTF <- rtracklayer::import(gtf_file)
# ens2biotype_ids <- with(GTF, setNames(gene_biotype, gene_id))
# sym2biotype_ids <- with(GTF, setNames(gene_biotype, gene_name))

#' .clusterwiseFindMarkers
#' @description Employs the simple FindMarkers() function from seurat
#' to do a pairwise wilcoxon test on RNA::data between two groups (group_id). 
#' It will iterate this test over all clusters (idents_id).
#' For example, iterate over all seurat_clusters (idents_id), and do a test
#' between both samples in [group_id] (Sample-2[comp.id] vs Sample-1[ref.id])
#' 
#' @param obj  seurat object
#' @param idents_id Group ids to use for differential testing (e.g. 'seurat_clusters')
#' @param group_id Identity id to use two groups (e.g. 'orig.ident')
#' @param comp.id IdentityID[1] from orig.ident to use as comparison group
#' @param ref.id IdentityID[2] from orig.ident to use as reference group
#'
#' @examples 
#' .clusterwiseFindMarkers(seu, id, 'orig.ident', 'MDSC_DAB', 'MDSC_DMSO')
.clusterwiseFindMarkers <- function(obj, idents_id, group_id, comp.id, ref.id){
  Idents(obj) <- idents_id
  cluster_ids <- as.character(unique(Idents(obj)))
  
  # Iterate through all clusters and find DEG between each group
  all_markers <- lapply(setNames(cluster_ids,cluster_ids), function(cl_id){
    print(paste0(cl_id, "..."))
    cl_markers <- FindMarkers(obj, assay='RNA', slot = "data",
                              test.use = "wilcox", base=2,
                              ident.1=comp.id, ident.2=ref.id,
                              subset.ident=cl_id, group.by=group_id)
    cl_markers <- cl_markers %>% 
      tibble::rownames_to_column('gene') %>%
      mutate(cluster=cl_id,
             #biotype=sym2biotype_ids[gene],
             comp_group=comp.id,
             ref_group=ref.id) %>%
      relocate(., c('gene', 'cluster', 'comp_group', 'ref_group')) #'biotype'
    return(cl_markers)
  })
  return(all_markers)
}

#' .betweenClustersFindMarkers
#' @description Employs the simple FindMarkers() function from seurat
#' to do a pairwise wilcoxon test on RNA::data between clusters as specified
#' in the [idents_id] parameter. For instance, if idents_id is 'seurat_clusters',
#' it will person pairwise tests between every unique pairing of clusters. 
#' With the resulting DEGs, it will assign a rank to each gene based on 
#' how discriminative that gene is across all clusters. This is done based
#' on the aggregate Log2FC across each comparison, aggregate -log10(padj), or 
#' a combination of the two (|LFC| * -log10(padj)). 
#' 
#' @param obj  seurat object
#' @param idents_id Group ids to use for differential testing (e.g. 'seurat_clusters')
#' @param score A score to rank discriminative genes one: 'lfc', 'padj', 'score', 'gsea'
#'
#' @examples 
#' .betweenClustersFindMarkers(seu, 'seurat_clusters', 'score')
.betweenClustersFindMarkers <- function(obj, idents_id, score='score'){
  require(dplyr)
  require(tibble)
  require(purrr)
  
  # Create all unique combinations of cluster-pairings
  Idents(obj) <- idents_id
  cluster_ids <- as.character(unique(Idents(obj)))
  cluster_ids <- setNames(cluster_ids,cluster_ids)
  cluster_ids_pair <- combn(cluster_ids, m=2)
  
  # Iterate through and compare every cluster to every other cluster
  all_markers <- apply(cluster_ids_pair, 2,function(cl_ids){
    print(paste0(cl_ids[1], " - ", cl_ids[2], "..."))
    cl_markers <- FindMarkers(obj, assay='RNA', slot = "data",
                              test.use = "wilcox", base=2,
                              ident.1=cl_ids[1], ident.2=cl_ids[2],
                              group.by=idents_id)
    return(cl_markers)
  })
  names(all_markers) <- apply(cluster_ids_pair, 2, paste, collapse="_")
  
  # Merge the LFC and padj across all cluster-wise comparisons into a single matrix
  .cleanit <- function(i, colid){
    i %>% 
      tibble::rownames_to_column(., "gene") %>%
      dplyr::select(gene, !!rlang::sym(colid))
  }
  .mergeit <- function(x, colid){
    lapply(x, .cleanit, colid=colid) %>% 
      purrr::reduce(full_join, by='gene') %>%
      tibble::column_to_rownames(., "gene") %>%
      rename_with(., ~paste(names(x), colid, sep="."))
  }
  pseudop <- 1*10^-5
  lfc_mat <- .mergeit(all_markers, colid='avg_log2FC')
  padj_mat <- .mergeit(all_markers, colid='p_val_adj')
  score_mat <- abs(lfc_mat) * -log10(padj_mat+pseudop)
  gsea_mat <- sign(lfc_mat) * -log10(padj_mat+pseudop)
  
  # Rank the genes based on their aggregate score to find most discriminative
  score_mat <- switch(score,
                      'score'=score_mat,
                      'gsea'=gsea_mat,
                      'lfc'=abs(lfc_mat),
                      'padj'=-log10(padj_mat+pseudop))
  score_rank_mat <- apply(score_mat, 2, order, 
                          decreasing=T, na.last=T) # low-value = discriminative
  gene_rank <- order(rowSums(score_rank_mat, na.rm = T), 
                     decreasing = F, na.last=T) # low-value = discriminative
  
  return(list("lfc"=lfc_mat[gene_rank,,drop=F],
            "padj"=padj_mat[gene_rank,,drop=F],
            "score"=score_mat[gene_rank,,drop=F]))
}

#' .identifyClusterSpecificMarkers
#' @description Given a matrix of values (set up for padj as of this moment), 
#' it will isolate cluster specific comparisons from all other comparisons and 
#' find features that are significant to that cluster. For example, in a 
#' pariwise comparison between clusters 1, 2, 3 and focusing cluster 1, it will
#' find features that are signficant in cluster 1 comparisons (1-2, 1-3) and compare
#' them to nfeatures that are not signifcant  in non-cluster 1 comparisons (2-3). It
#' finds these features based on the quantile (q[1]) of cluster 1 padj values compared 
#' to the quantile_2 (q[2]) of non-cluster 1 padj values per feature. As this function 
#' computes the -1*log10(padj), q[1] should be low and q[2] should be high to maximize 
#' feature specificity
#' 
#' @param bcobj object returned from .betweenClustersFindMarkers() function
#' @param slot character: only works with 'padj' & 'score' for now
#' @param pseudoval pseudoval to add [default=0.0001]
#' @param n number of top features to return [n=100]
#' @param q quantiles to use as threshold [q=c(0.05, 0.95)]
#' 
#' @examples
#' markers <- .betweenClustersFindMarkers(seu, 'seurat_clusters', 'score')
#' .identifyClusterSpecificMarkers(markers)
.identifyClusterSpecificMarkers <- function(bcobj, slot='padj', q=c(0.05, 0.95),
                                            pseudoval=0.0001, n=100){
  if(!(slot %in% c('padj', 'score'))) warning("slot was not tested with anything other than padj/score")
  dat <- bcobj[[slot]]
  if(slot == 'padj'){
    dat[is.na(dat)] <- 1
    dat <- -1*log10(dat + pseudoval)
  } else if (slot=='score'){
    dat[is.na(dat) | dat < 0] <- 0
  }
  
  if(ncol(dat) == 1){
    uniq_markers <- list(
      sort(setNames(dat[,1], rownames(dat)), decreasing=T) %>%
        head(., n) %>% 
        t %>% melt %>%
        mutate(Var1=uc_i)
    )
  } else {
    meta <- gsub("\\..*", "", colnames(dat)) %>%
      strsplit(., split="_") %>%
      do.call(rbind, .)
    uniq_clus <- unique(as.character(meta))
    
    # Isolate individual clusters
    uniq_markers <- lapply(uniq_clus, function(uc_i){
      i_idx <- which(rowSums(meta == uc_i) > 0)
      j_idx <- c(1:nrow(meta))[-i_idx]
      if((length(i_idx) > 4) & (length(j_idx) > 5)){
        i_summ <- apply(dat[,i_idx], 1, quantile, na.rm = T, probs=q[1])
        j_summ <- apply(dat[,j_idx], 1, quantile, na.rm = T, probs=q[2])
      } else {
        i_summ <- apply(dat[,i_idx,drop=F], 1, median, na.rm = T)
        j_summ <- apply(dat[,j_idx,drop=F], 1, median, na.rm = T)
      }
      ij_ratio <- (i_summ / (j_summ + pseudoval))
      sort(ij_ratio, decreasing=T) %>%
        head(., n) %>% 
        t %>% melt %>%
        mutate(Var1=uc_i)
    })
  }
  
  return(uniq_markers)
}