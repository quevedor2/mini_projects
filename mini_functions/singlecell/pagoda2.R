# counts=expr[idx,]
# gSets=msig_l
# trim = 5
# n_cores=1
# min_gset_size=5
# max_gset_size=1000

cal_pagoda2 = function(counts,
                       gSets,
                       trim = 5,
                       n_cores=1,
                       min_gset_size=5,
                       max_gset_size=1000, ...){
  
  # counts <- expr
  # gSets <- msig_l
  # max_gset_size <- 1500
  # trim = 5
  # n_cores=1
  # min_gset_size=5
  ### must be counts matrix !!!!!
  
  ### other parameters for knn.error.models
  # min.nonfailed = 5, min.count.threshold = 1,
  # max.model.plots = 50,
  # min.size.entries = 2000, min.fpm = 0, cor.method = "pearson",
  # verbose = 0, fpm.estimate.trim = 0.25, linear.fit = TRUE,
  # local.theta.fit = linear.fit, theta.fit.range = c(0.01, 100),
  # alpha.weight.power = 1/2
  
  ### other parameters for pagoda.varnorm
  # batch = NULL, prior = NULL,
  # fit.genes = NULL, minimize.underdispersion = FALSE,
  # n.cores = detectCores(), n.randomizations = 100, weight.k = 0.9,
  # verbose = 0, weight.df.power = 1, smooth.df = -1,
  # theta.range = c(0.01, 100), gene.length = NULL
  
  nPcs = min(round(ncol(counts)/5),5)
  #counts = apply(counts,2,function(x) {storage.mode(x) = 'integer'; x})
  tryCatch({
    p2 = Pagoda2$new(as(counts, "dgCMatrix"), n.cores = n_cores,log.scale=F, modelType="raw")
    p2$adjustVariance(plot=F)
    p2$calculatePcaReduction(nPcs = nPcs,use.odgenes=F,fastpath=F)
    
    proper.gene.names <- rownames(counts)
    sum_cnt <- sapply(gSets, function(sn) sum(proper.gene.names %in% sn))
    rm_idx <- which(sum_cnt <= min_gset_size)
    if(length(rm_idx)>0) {
      warning(paste("Removing genesets for being too small: ", names(gSets)[rm_idx], sep="\n"))
      gSets <- gSets[-rm_idx]
    }
    
    path_names = c()
    env <- list2env(gSets)
    
    # p2$calculatePcaReduction(nPcs=50,n.odgenes=3e3)
    # p2$makeKnnGraph(k=40,type='PCA',center=T,distance='cosine');
    # p2$getKnnClusters(method=infomap.community,type='PCA')
    # p2$getHierarchicalDiffExpressionAspects(type='PCA',clusterName='community',z.threshold=3)

    p2 <- testPathwayOverdispersion.relaxed(self=p2, setenv = env, verbose = T,
                                            recalculate.pca = F,
                                            min.pathway.size =min_gset_size,
                                            max.pathway.size=max_gset_size, ...)
    
    path_names = names(p2$misc$pwpca)
    score = matrix(NA,nrow=length(path_names),ncol=ncol(counts))
    rownames(score) = path_names
    colnames(score) = colnames(counts)
    for(i in 1:length(p2$misc$pwpca)){
      if(!is.null(p2$misc$pwpca[[i]]$xp$score)){
        score[i,] = as.numeric(p2$misc$pwpca[[i]]$xp$scores)
      }
    }
    
    return(score)
  },error = function(e){
    print(e)
    return("error")
  })
}


# self <- p2
# type='counts'
# max.pathway.size=1e3
# min.pathway.size=10
# n.randomizations=5
# verbose=FALSE
# n.cores=self$n.cores
# score.alpha=0.05
# plot=FALSE
# cells=NULL
# adjusted.pvalues=TRUE
# z.score = qnorm(0.05/2, lower.tail = FALSE)
# use.oe.scale = FALSE
# return.table=FALSE
# name='pathwayPCA'
# correlation.distance.threshold=0.2
# loading.distance.threshold=0.01
# top.aspects=Inf
# recalculate.pca=FALSE
# save.pca=TRUE
# setenv = env
# verbose = T
# recalculate.pca = F
# min.pathway.size =1

testPathwayOverdispersion.relaxed=function(self, setenv, type='counts', max.pathway.size=1e3, min.pathway.size=10, 
                                           n.randomizations=5, verbose=FALSE, n.cores=self$n.cores, score.alpha=0.05, plot=FALSE, cells=NULL, adjusted.pvalues=TRUE,
                                           z.score = qnorm(0.05/2, lower.tail = FALSE), use.oe.scale = FALSE, return.table=FALSE, name='pathwayPCA',
                                           correlation.distance.threshold=0.2, loading.distance.threshold=0.01, top.aspects=Inf, recalculate.pca=FALSE, save.pca=TRUE) {
  
  if (!requireNamespace("scde", quietly=TRUE)){
    stop("You need to install package 'scde' to be able to use testPathwayOverdispersion().")
  }
  
  nPcs <- 1
  if (type=='counts') {
    x <- self$counts
    # apply scaling if using raw counts
    x@x <- x@x*rep(self$misc[['varinfo']][colnames(x),'gsf'],diff(x@p))
  } else {
    if (!type %in% names(self$reductions)) { stop("Reduction ",type,' not found')}
    x <- self$reductions[[type]]
  }
  if (!is.null(cells)) {
    x <- x[cells,]
  }
  
  proper.gene.names <- colnames(x)
  
  if (is.null(self$misc[['pwpca']]) || recalculate.pca) {
    if (verbose) {
      message("determining valid pathways")
    }
    
    # determine valid pathways
    gsl <- ls(envir = setenv)
    gsl.ng <- unlist(parallel::mclapply(pagoda2:::sn(gsl), function(go) sum(unique(get(go, envir = setenv)) %in% proper.gene.names),mc.cores=n.cores,mc.preschedule=TRUE))
    gsl <- gsl[gsl.ng >= min.pathway.size & gsl.ng<= max.pathway.size]
    names(gsl) <- gsl
    
    if (verbose) {
      message("processing ", length(gsl), " valid pathways")
    }
    
    cm <- Matrix::colMeans(x)
    
    pwpca <- pagoda2:::papply(gsl, function(sn) {
      lab <- proper.gene.names %in% get(sn, envir = setenv)
      if (sum(lab)<1) { 
        return(NULL)
      }
      pcs <- irlba::irlba(x[,lab], nv=nPcs, nu=0, center=cm[lab])
      pcs$d <- pcs$d/sqrt(nrow(x))
      pcs$rotation <- pcs$v
      pcs$v <- NULL
      
      # get standard deviations for the random samples
      ngenes <- sum(lab)
      z <- do.call(rbind,lapply(seq_len(n.randomizations), function(i) {
        si <- sample(ncol(x), ngenes)
        pcs <- irlba::irlba(x[,si], nv=nPcs, nu=0, center=cm[si])$d
      }))
      z <- z/sqrt(nrow(x))
      
      # local normalization of each component relative to sampled PC1 sd
      avar <- pmax(0, (pcs$d^2 - mean(z[, 1]^2))/sd(z[, 1]^2))
      print(paste0(sn, " - ", avar))
      
      if (TRUE) { #(avar>=0) {
        # flip orientations to roughly correspond with the means
        pcs$scores <- as.matrix(t(x[,lab] %*% pcs$rotation) - as.numeric((cm[lab] %*% pcs$rotation)))
        cs <- unlist(lapply(seq_len(nrow(pcs$scores)), function(i) sign(cor(pcs$scores[i,], colMeans(t(x[, lab, drop = FALSE])*abs(pcs$rotation[, i]))))))
        pcs$scores <- pcs$scores*cs
        pcs$rotation <- pcs$rotation*cs
        rownames(pcs$rotation) <- colnames(x)[lab]
      } # don't bother otherwise - it's not significant
      return(list(xp=pcs,z=z,n=ngenes))
    }, n.cores = n.cores,mc.preschedule=TRUE)
    if (save.pca) {
      self$misc[['pwpca']] <- pwpca
    }
  } else {
    if (verbose) {
      message("reusing previous overdispersion calculations")
      pwpca <- self$misc[['pwpca']]
    }
  }
  
  if (verbose) {
    message("scoring pathway od signifcance")
  }
  
  # score overdispersion
  true.n.cells <- nrow(x)
  
  pagoda.effective.cells <- function(pwpca, start = NULL) {
    n.genes <- unlist(lapply(pwpca, function(x) rep(x$n, nrow(x$z))))
    var <- unlist(lapply(pwpca, function(x) x$z[, 1]))
    if (is.null(start)) { start <- true.n.cells*2 } # start with a high value
    of <- function(p, v, sp) {
      sn <- p[1]
      vfit <- (sn+sp)^2/(sn*sn+1/2) -1.2065335745820*(sn+sp)*((1/sn + 1/sp)^(1/3))/(sn*sn+1/2)
      residuals <- (v-vfit)^2
      return(sum(residuals))
    }
    x <- nlminb(objective = of, start = c(start), v = var, sp = sqrt(n.genes-1/2), lower = c(1), upper = c(true.n.cells))
    return((x$par)^2+1/2)
  }
  n.cells <- pagoda.effective.cells(pwpca)
  
  vdf <- data.frame(do.call(rbind, lapply(seq_along(pwpca), function(i) {
    vars <- as.numeric((pwpca[[i]]$xp$d))
    cbind(i = i, var = vars, n = pwpca[[i]]$n, npc = seq(1:ncol(pwpca[[i]]$xp$rotation)))
  })))
  
  # fix p-to-q mistake in qWishartSpike
  qWishartSpikeFixed <- function (q, spike, ndf = NA, pdim = NA, var = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)  {
    params <- RMTstat::WishartSpikePar(spike, ndf, pdim, var, beta)
    qnorm(q, mean = params$centering, sd = params$scaling, lower.tail, log.p)
  }
  
  # add right tail approximation to ptw, which gives up quite early
  pWishartMaxFixed <- function (q, ndf, pdim, var = 1, beta = 1, lower.tail = TRUE) {
    params <- RMTstat::WishartMaxPar(ndf, pdim, var, beta)
    q.tw <- (q - params$centering)/(params$scaling)
    p <- RMTstat::ptw(q.tw, beta, lower.tail, log.p = TRUE)
    p[p == -Inf] <- pgamma((2/3)*q.tw[p == -Inf]^(3/2), 2/3, lower.tail = FALSE, log.p = TRUE) + lgamma(2/3) + log((2/3)^(1/3))
    p
  }
  
  vshift <- 0
  ev <- 0
  
  vdf$var <- vdf$var-(vshift-ev)*vdf$n
  basevar <- 1
  vdf$exp <- RMTstat::qWishartMax(0.5, n.cells, vdf$n, var = basevar, lower.tail = FALSE)
  #vdf$z <- qnorm(pWishartMax(vdf$var, n.cells, vdf$n, log.p = TRUE, lower.tail = FALSE, var = basevar), lower.tail = FALSE, log.p = TRUE)
  vdf$z <- qnorm(pWishartMaxFixed(vdf$var, n.cells, vdf$n, lower.tail = FALSE, var = basevar), lower.tail = FALSE, log.p = TRUE)
  vdf$cz <- qnorm(pagoda2:::bh.adjust(pnorm(as.numeric(vdf$z), lower.tail = FALSE, log.p = TRUE), log = TRUE), lower.tail = FALSE, log.p = TRUE)
  vdf$ub <- RMTstat::qWishartMax(score.alpha/2, n.cells, vdf$n, var = basevar, lower.tail = FALSE)
  vdf$ub.stringent <- RMTstat::qWishartMax(score.alpha/nrow(vdf)/2, n.cells, vdf$n, var = basevar, lower.tail = FALSE)
  
  if (plot) {
    test_pathway_par <- par(mfrow = c(1, 1), mar = c(3.5, 3.5, 1.0, 1.0), mgp = c(2, 0.65, 0))
    on.exit(par(test_pathway_par))
    un <- sort(unique(vdf$n))
    on <- order(vdf$n, decreasing = FALSE)
    pccol <- colorRampPalette(c("black", "grey70"), space = "Lab")(max(vdf$npc))
    plot(vdf$n, vdf$var/vdf$n, xlab = "gene set size", ylab = "PC1 var/n", ylim = c(0, max(vdf$var/vdf$n)), col = adjustcolor(pccol[vdf$npc],alpha=0.1),pch=19)
    lines(vdf$n[on], (vdf$exp/vdf$n)[on], col = 2, lty = 1)
    lines(vdf$n[on], (vdf$ub.stringent/vdf$n)[on], col = 2, lty = 2)
  }
  
  rs <- (vshift-ev)*vdf$n
  vdf$oe <- (vdf$var+rs)/(vdf$exp+rs)
  vdf$oec <- (vdf$var+rs)/(vdf$ub+rs)
  
  df <- data.frame(name = names(pwpca)[vdf$i], npc = vdf$npc, n = vdf$n, score = vdf$oe, z = vdf$z, adj.z = vdf$cz, stringsAsFactors = FALSE)
  if (adjusted.pvalues) {
    vdf$valid <- vdf$cz  >=  z.score
  } else {
    vdf$valid <- vdf$z  >=  z.score
  }
  
  if (!any(vdf$valid)) { 
    stop("No significantly overdispersed pathways found at z.score threshold of ",z.score) 
  }
  
  # apply additional filtering based on >0.5 sd above the local random estimate
  vdf$valid <- vdf$valid & unlist(lapply(pwpca,function(x) !is.null(x$xp$scores)))
  vdf$name <- names(pwpca)[vdf$i]
  
  if (return.table) {
    df <- df[vdf$valid, ]
    df <- df[order(df$score, decreasing = TRUE), ]
    return(df)
  }
  if (verbose) {
    message("compiling pathway reduction")
  }
  # calculate pathway reduction matrix
  
  # return scaled patterns
  xmv <- do.call(rbind, lapply(pwpca[vdf$valid], function(x) {
    xm <- x$xp$scores
  }))
  
  if (use.oe.scale) {
    xmv <- (xmv -rowMeans(xmv))* (as.numeric(vdf$oe[vdf$valid])/sqrt(apply(xmv, 1, var)))
    vdf$sd <- as.numeric(vdf$oe)
  } else {
    # chi-squared
    xmv <- (xmv-rowMeans(xmv)) * sqrt((qchisq(pnorm(vdf$z[vdf$valid], lower.tail = FALSE, log.p = TRUE), n.cells, lower.tail = FALSE, log.p = TRUE)/n.cells)/apply(xmv, 1, var))
    vdf$sd <- sqrt((qchisq(pnorm(vdf$z, lower.tail = FALSE, log.p = TRUE), n.cells, lower.tail = FALSE, log.p = TRUE)/n.cells))
    
  }
  rownames(xmv) <- paste("#PC", vdf$npc[vdf$valid], "# ", names(pwpca)[vdf$i[vdf$valid]], sep = "")
  rownames(vdf) <- paste("#PC", vdf$npc, "# ", vdf$name, sep = "")
  self$misc[['pathwayODInfo']] <- vdf
  
  # collapse gene loading
  if (verbose) {
    message("clustering aspects based on gene loading ... ",appendLF=FALSE)
  }
  tam2 <- pagoda.reduce.loading.redundancy(list(xv=xmv,xvw=matrix(1,ncol=ncol(xmv),nrow=nrow(xmv))),
                                           pwpca,NULL,plot=FALSE,
                                           distance.threshold=loading.distance.threshold,n.cores=n.cores)
  if (verbose) {
    message(nrow(tam2$xv)," aspects remaining")
  }
  if (verbose) {
    message("clustering aspects based on pattern similarity ... ",appendLF=FALSE)
  }
  tam3 <- tryCatch({
    pagoda.reduce.redundancy(tam2, distance.threshold=correlation.distance.threshold,top=top.aspects)
  }, error=function(e){
    tam2
  })
  if (verbose) {
    message(nrow(tam3$xv)," aspects remaining\n")
  }
  tam2$xvw <- tam3$xvw <- NULL # to save space
  tam3$env <- setenv
  
  # clean up aspect names, as GO ids are meaningless
  names(tam3$cnam) <- rownames(tam3$xv) <- paste0('aspect',1:nrow(tam3$xv))
  
  self$misc[['pathwayOD']] <- tam3
  self$reductions[[name]] <- tam3$xv
  # invisible(tam3)
  return(self)
}
