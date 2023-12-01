## This function is designed to extract the top 5% of genes per cell that were used in the 
# AUCell function.  It is pieced together from the AUC function, but instead of 
# reporting the auc (in the .auc() function), it now reports the genes that are
# included in the AUC metric (.auc_report())
# Afterwards, this function will tally up how many times gene_X was used in the 
# auc score across all cells and return a dataframe of gene by cell-count
getGenesUsedInAucell <- function(exprMat,  geneSets,
                     featureType='genes', keepZeroesAsNA = FALSE,
                     nGenesDetected=numeric(0)){
  aucMaxRank = ceiling(0.05 * nrow(exprMat))

  ########## BuildRanks #### 
  ## Calculate expression stats?
  ## Rank genes
  rowNames <- rownames(exprMat)
  colNames <- colnames(exprMat)
  exprMat <- -exprMat # ro rank decreasingly
  exprMat2 <- do.call(cbind,
                     DelayedArray::blockApply(DelayedArray::DelayedArray(exprMat),
                                              FUN=DelayedMatrixStats::colRanks,
                                              ties.method="random", preserveShape=TRUE,
                                              BPPARAM=NULL,
                                              grid=DelayedArray::colAutoGrid(exprMat)))
  rownames(exprMat2) <- rowNames
  colnames(exprMat2) <- colNames
  
  if(keepZeroesAsNA){
    exprMat2[which(zeroesCoords==0, arr.ind=TRUE)] <- NA
  }
  
  ## Format & return
  names(dimnames(exprMat2)) <- c(featureType, "cells")
  rankings <- new("aucellResults",
                  SummarizedExperiment::SummarizedExperiment(assays=list(ranking=exprMat2)),
                  nGenesDetected=nGenesDetected)
  
  ########## CalcAUC #### 
  rankings2 <- AUCell:::getRanking(rankings)
  
  .auc_report <- function(oneRanking, aucThreshold, maxAUC) {
    x <- unlist(oneRanking)
    
    x <- sort(x[x<aucThreshold])
    data.frame("genes"=names(x),
               "rank"=as.integer(x))
  }
  .AUC.geneSet_norm <- function(geneSet, rankings, aucMaxRank, gSetName="") {
    geneSet <- unique(geneSet)
    nGenes <- length(geneSet)
    geneSet <- geneSet[which(geneSet %in% rownames(rankings))]
    missing <- nGenes-length(geneSet)
    
    gSetRanks <- rankings[which(rownames(rankings) %in% geneSet),,drop=FALSE]
    rm(rankings)
    
    aucThreshold <- round(aucMaxRank)
    ########### NEW version:  #######################
    x_th <- 1:nrow(gSetRanks)
    x_th <- sort(x_th[x_th<aucThreshold])
    y_th <- seq_along(x_th)
    maxAUC <- sum(diff(c(x_th, aucThreshold)) * y_th) 
    ############################################
    
    # Apply by columns (i.e. to each ranking)
    auc <- apply(gSetRanks, 2, .auc_report, aucThreshold, maxAUC)
    auc <- lapply(auc, function(i) i$genes) %>% 
      unlist %>% 
      table %>% sort
    # x <- auc %>% 
    #   purrr:::reduce(left_join, by='genes')
    
    return(c(auc, missing=missing, nGenes=nGenes))
  }
  
  aucMatrix <- lapply(names(geneSets), function(gSetName){
    # gSetName <- names(geneSets)[1]
    .AUC.geneSet_norm(geneSet=geneSets[[gSetName]], 
                 rankings=rankings2,
                 aucMaxRank=aucMaxRank, 
                 gSetName=gSetName)
  })
  cnt_mat <- lapply(aucMatrix, function(i){
    data.frame('gene'=names(i), 'count'=as.integer(i))
  }) %>% 
    purrr::reduce(full_join, by='gene') %>%
    magrittr::set_colnames(c('gene', names(geneSets)))
  return(cnt_mat)
}