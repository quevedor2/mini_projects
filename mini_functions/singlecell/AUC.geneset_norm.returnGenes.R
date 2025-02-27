## Modified AUCell run function to return the ranks of the genes meeting the treshold in each geneset
#' .AUC.geneset_norm.returnGenes
#' @description A slight modification of the AUCell:::AUCell_calcAUC() function. Instead of 
#' just returning the aucell per score, now you have an added parameter: returntype, 
#' which allows you return either the cumulative AUC score across all top-scoring
#' genes in the geneset, or a list of all the genes that made the leading edge
#' 
#' @param geneSet 
#' @param rankings 
#' @param aucMaxRank 
#' @param gSetName 
#' @param returntype Either 'cumulative_auc', 'genes', or 'auc' where 'auc' is the 
#' modified AUCell:::AUCell_calcAUC() function.  and 'cumulative_auc' breaks down the
#' auc score to the cumulative-sum (cumulative auc) across each gene in the top leading
#' genes. while the 'genes' param will just return the genes in the leading edge and their
#' unscaled ranks
.AUC.geneset_norm.returnGenes <- function(geneSet, rankings, aucMaxRank, gSetName="", returntype='genes'){
  message("returntype  can be either 'genes', 'cumulative_auc', or 'auc'")
  
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
  .aucrank <- function(oneRanking, aucThreshold, maxAUC, returntype=returntype, verbose=F) {
    x <- unlist(oneRanking)
    x <- sort(x[x < aucThreshold])
    y <- seq_along(x)                             # normal function for AUCell:::.auc
    scaled_diff <- diff(c(x, aucThreshold)) * y   # normal function for AUCell:::.auc
    auc_final <- sum(scaled_diff)/maxAUC          # normal function for AUCell:::.auc
    auc_cumulative <- cumsum(scaled_diff) / maxAUC
    names(auc_cumulative) <- names(x)
    
    ranks <- c(x, aucThreshold)
    names(ranks)[length(ranks)] <- 'n'
    
    if(returntype == 'genes'){
      if(verbose) message("Returning the top genes that were used for AUC scoring")
      returnval <- as.data.frame(t(ranks))
    } else if(returntype == 'cumulative_auc'){
      if(verbose) message("Returning the cumulative AUC score for the top genes")
      returnval <-  as.data.frame(t(auc_cumulative))
    } else if(returntype == 'auc'){
      if(verbose) message("Running in standard AUCell mode, returning AUC value")
      returnval <- auc_final
    } else {
      stop("Return value must be either 'genes', 'cumulative_auc', or'auc'")
    }
    
    return(returnval)
  }
  
  # Apply by columns (i.e. to each ranking)
  auc <- apply(gSetRanks, 2, .aucrank, aucThreshold, maxAUC, returntype=returntype) %>%
    plyr::rbind.fill(.)
  return(auc)
}

#' plotCumulativeAUC
#' @description  X
#' 
#' @param X 
#' @param summary_data 
#' @param groupid 
#' @param mk.plot 
#'
#' @return
#' @export
#'
#' @examples
plotCumulativeAUC <- function(X=NULL, summary_data=NULL, n=NULL, 
                              title=NULL,  groupid = NULL, mk.plot=TRUE){
  if(is.null(summary_data)){
    # Transform data into a simple triplet matrix  with x indices added
    df <- reshape2::melt(t(X)) %>%
      dplyr::filter(!is.na(value)) %>%
      dplyr::arrange(Var2, value) %>%
      group_split(Var2) %>%
      lapply(., function(i) {i$xidx <- seq_len(nrow(i)); i}) %>%
      do.call(rbind, .)
    if(!is.null(groupid)) df$groupid <- groupid
    
    # Calculate mean and confidence interval around each index
    summary_data <- df %>%
      dplyr::group_by(xidx) %>%
      dplyr::mutate(n=n()) %>%
      dplyr::summarize(
        mean_y = mean(value),
        sd_y = sd(value),
        se_y = sd_y / sqrt(n), # Standard error
        ci_lower = mean_y - 1.96 * se_y, # 95% CI lower bound
        ci_upper = mean_y + 1.96 * se_y  # 95% CI upper bound
      ) %>% unique
    if(!is.null(groupid)) summary_data$groupid <- groupid
  }
  returnval <- summary_data
  
  if(mk.plot){
    if(!'groupid' %in% colnames(summary_data)){
      summary_data$groupid <- 'XX'
    }
    max_data <- summary_data %>%
      dplyr::group_by(groupid) %>%
      dplyr::filter(xidx == max(xidx))
    
    # pdf("~/xfer/cumauc.pdf")
    p <- ggplot(summary_data, aes(x = xidx, y = mean_y, color=groupid, fill=groupid)) +
      geom_line(linewidth = 1) +
      geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2) +
      geom_hline(data=max_data, aes(yintercept=mean_y, color=groupid), linetype = "dashed") +
      labs(title = paste0("Cumulative AUC: ", title),
           x = "Geneset Rank",
           y = "AUC") + 
      cowplot::theme_cowplot() +
      xlim(1, n) +
      scale_y_log10(
        breaks = c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 1),
        labels = c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 1), #scales::trans_format("log10", scales::math_format(10^.x)),
        limits = c(min(summary_data$mean_y[summary_data$mean_y > 0]), 1) # Ensure no 0s, set limits
      ) +
      annotation_logticks(sides = "l")
    # p
    # dev.off()
    returnval <- p
  }
  return(returnval)
}


run_demo = FALSE
if(run_demo){
  returntype <- 'cumulative_auc' # 'genes', 'auc'
  
  DefaultAssay(cd8seu) <- 'RNA'
  Idents(cd8seu) <- 'group'
  
  expr <- GetAssayData(subset(cd8seu, ident=id_i), slot='counts')
  msig_l <- list("GenesetX"="Ackr4","Alpk1","Apobec3","Apol9b","Apol7e","Apol8",
                 "Apol6","Batf2","Bst2","C1rb","C1s2","C3","Casp1","Ccdc68","Ccl12",
                 "Cd274","Cd40","Cd74","Ceacam2","Cfb","Cfh","Cfhr1","Cmpk2","Ctss",
                 "Cxcl10","Cxcl11","Cxcl9","Cyp1b1","Cyp2j7","Ddx60","Dhx58","Dtx3l",
                 "Eif2ak2","Epsti1","Fam117b","Fam20a","Fbxo6","Flt3l","Gabbr1","Gbp2",
                 "Gbp3","Gbp5","Gmpr","Hapln3","Helz2","Herc6","H2-T22","H2-DMa","H2-DMb2",
                 "H2-Ab1","H2-Ea","Icam1","Ido1","Ifi207","Ifi27","Ifi30","Ifi35","Ifi44",
                 "Ifi44l","Ifih1","Ifit1bl1","Ifit2","Ifit3b","Ifitm1","Ifitm2","Ifitm3",
                 "Igflr1","Il15","Il15ra","Il18bp","Irf1","Irf7","Irf9","Isg15","Isg20",
                 "Itpkc","Jak2","Lamp3","Lap3","Lgals3bp","Lmo2","Mlkl","Muc1","Mx2","Myd88",
                 "Nlrc5","Nmi","Oas1e","Oas2","Oas3","Oasl1","P2rx7","Parp10","Parp12","Parp14",
                 "Parp9","Pdcd1lg2","Phf11c","Pik3r2","Plscr2","Plscr4","Pmaip1","Pml","Pnpt1",
                 "Psmb10","Psmb8","Psmb9","Psme2","Rarres1","Rasgrp3","Ripk2","Rnf19b","Rnf213",
                 "Rsad2","Rtp4","Samd9l","Samhd1","Sectm1b","Serping1","Slc15a3","Socs1","Socs3",
                 "Sp100","Sp110","Stat1","Stat2","Tap1","Tap2","Tapbpl","Tdrd7","Tesk2","Tlr3",
                 "Tmem140","Tnfsf10","Tnfsf13b","Trafd1","Trank1","Trim21","Trim25","Trim5",
                 "Trim69","Uba7","Ube2l6","Unc93b1","Usp18","Xaf1","Zc3hav1","Znfx1")
    
  # Breaking down the  [AUCell::AUCell_run(expr, msig_l)] function into its constituents:
  ranked <- AUCell:::AUCell_buildRankings(expr, splitByBlocks=FALSE, 
                                          featureType='genes',
                                          keepZeroesAsNA=FALSE, 
                                          plotStats=FALSE, verbose=FALSE)
  rankings <- getRanking(ranked)
  geneSets <- msig_l
  aucMaxRank <- ceiling(0.05*nrow(ranked))
  
  
  # Breaking down the [AUCell:::AUCell_calcAUC(geneSets, ranked, normAUC=T, 
  #                              aucMaxRank=aucMaxRank, verbose=FALSE)]
  # function into a similar function but with options to extract different bits of information
  aucMatrices <- lapply(names(geneSets), function(gSetName)
    # AUCell:::.AUC.geneSet_norm(geneSet=geneSets[[gSetName]], rankings=rankings,
    #              aucMaxRank=aucMaxRank, gSetName=gSetName)
    .AUC.geneset_norm.returnGenes(geneSet=geneSets[[gSetName]], rankings=rankings,
                                  aucMaxRank=aucMaxRank, gSetName=gSetName, returntype=returntype)
  ) %>% setNames(., names(geneSets))
  aucMatrices
}
