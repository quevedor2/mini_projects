#' Linear Model Eigenvalue matrix to a design matrix
#' @description in an attempt to identify the principle components (or eigengene)
#' module associated with a given metadata component (Deisgn matrix), this 
#' function will fit a linear model to a one-hot-encoded design matrix 
#' (multifactorial) and smoooth for standard errors to report the most
#' likely component associated with the metadata
#'
#' @param eigen nxm eigenvalue matrix, where n=components and m=samples
#' @param design_mat a mxi model.matrix(), where m=samples and i is the 
#' factors encoded
#' 
#' @examples
#' bwnet <- WGCNA::blockwiseModules(...)
#' module_eigengenes <- bwnet$MEs
#' group <- 'celltype' 
#' des_mat <- model.matrix(as.formula(paste0("~ coldata$", group)))
#' linearModelEigen(eigen=t(module_eigengenes), des_mat)
linearModelEigen <- function(eigen, design_mat){
  ## Git a linear model to the design matrix given the module eigengenes
  fit <- limma::lmFit(eigen, design = design_mat)
  fit <- limma::eBayes(fit) # Apply empirical Bayes to smooth standard errors
  
  ## Apply multiple testing correction and obtain stats
  stats_df <- limma::topTable(fit, number = nrow(eigen)) %>%
    tibble::rownames_to_column("module")
  return(stats_df)
}