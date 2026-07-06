
imbalance_score_enrichment <- function(rd, conditions, k = 10, smooth = k){
  message("Based on the 2-condition chisq std.res version of condiments:::.imbalance_score")
  if (length(conditions) != nrow(rd)) {
    stop("The conditions and reduced dimensions do not contain the same cells")
  }
  if (is.factor(conditions)) {
    groups <- levels(conditions)
  } else {
    groups <- unique(conditions)
  }
  groupidx <- match(as.character(conditions), groups) ## Match index of conditions to groups
  props <- as.vector(table(conditions)[groups]/length(conditions))
  if (length(groups) == 1) 
    stop("conditions should have at least 2 classes")
  tmp <- RANN::nn2(rd, rd, k + 1, searchtype = "standard")
  neighborMatrix <- tmp[[1]]
  cdMatrix <- matrix(factor(conditions)[neighborMatrix], ncol = k + 1)
  
  scores <- .chisq.stdres(cdMatrix, groups, props, groupidx)
  scores <- unlist(scores)
  names(scores) <- rownames(rd)
  formula <- paste0("scores ~ s(", 
                    paste0("rd[, ", seq_len(ncol(rd)), 
                           "], ", collapse = ""), "k = smooth)")
  mm <- mgcv::gam(stats::as.formula(formula))
  scaled_scores <- mgcv::predict.gam(mm, type = "response")
  return(list(scores = scores, scaled_scores = scaled_scores))
}

.chisq.stdres <- function(cdMatrix, groups, props, groupidx) {
  .tbl2df <- function(x) {
    as.data.frame(t(as.matrix(x)))
  }
  
  res <- apply(cdMatrix, 1, function(conds) {
    real <- as.vector(table(factor(conds, levels = groups)))
    .tbl2df(chisq.test(real, p = props)$stdres)
  })
  # Extract the STD.Res for the group matching the cell of interest
  res2 <- do.call(rbind, res)
  res <- res2[cbind(seq_len(nrow(res2)), groupidx)]  
  return(res)
}
