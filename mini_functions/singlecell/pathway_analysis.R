
cal_pagoda2 = function(counts,
                       gSets,
                       trim = 5,
                       n_cores=1){
  
  
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
    p2 = Pagoda2$new(counts, n.cores = n_cores,log.scale=F)
    p2$adjustVariance(plot=F)
    p2$calculatePcaReduction(nPcs = nPcs,use.odgenes=F,fastpath=F)
    
    path_names = c()
    env <- list2env(gSets)
    
    p2$testPathwayOverdispersion(setenv = env, verbose = T,
                                 recalculate.pca = T,
                                 min.pathway.size = 1)
    
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


#### Run Pagoda2 ####
## Set up the gene-sets:
msig_ds <- msigdbr(species = 'Homo sapiens', category = 'C2', subcategory = 'REACTOME') %>%
  dplyr::select(gs_name, entrez_gene) %>%
  as.data.frame()
rand_gs <- unique(msig_ds$gs_name)[sample(1:length(unique(msig_ds$gs_name)), size=50,replace=F)]

gm <- geneMap('Homo sapiens')
mapped_id <- gm[['ENTREZID']][['SYMBOL']][as.character(msig_ds$entrez_gene)]
msig_l <- lapply(split(mapped_id, msig_ds$gs_name), function(i){
  i[!is.na(i)]
})

## Run Pagoda2:
expr <- GetAssayData(seu, slot='counts')
idx <- which(rowSums(expr==0) < (ncol(seu)*0.95))
pagoda2_scores<- cal_pagoda2(expr[idx,], msig_l)