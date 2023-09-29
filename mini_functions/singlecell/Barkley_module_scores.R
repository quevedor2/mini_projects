run_GSVA_ssGSEA <- function (seur=NULL, msigdbr_specie="Homo sapiens", assay='RNA', 
                             slot='counts', msigdbr_category='H', 
                             msigdbr_subcategory=NULL, distribution="Poisson", 
                             method='gsva', genesets_ls=NULL, ncores=10, 
                             expr_data=NULL) {
  
  if (is.null(seur) & is.null(expr_data)) {
    print('Error: a Seurat object or a gene expression matrix is needed')
    return(NULL)    
  }
  
  if (!is.null(seur)) {
    meta_data = seur@meta.data
    expr_data = as.matrix(LayerData(seur, assay=assay, slot=slot)) # can be applied on slot data or counts
  }
  
  if (is.null(genesets_ls)) {
    genesets_df = msigdbr(species=msigdbr_specie, category=msigdbr_category, subcategory=msigdbr_subcategory) 
    genesets_df = msigdbr_df %>% 
      dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
    genesets_ls = geneset_df %>% split(x=.$gene_symbol, f=.$gs_name)
  }
  
  # run GSVA on genesets
  res = GSVA::gsva(expr=expr_data, gset.idx.list=genesets_ls, kcdf=distribution, 
                   method=method, min.sz=5, max.sz=500, parallel.sz=ncores)
  
  # if only one geneset is tested return a vector instead of a dataframe
  if (nrow(res)==1) { 
    res = structure(c(res), names=colnames(res))
  }
  res
}  
# library(GSVA)
# library(dplyr)
# library(msigdbr)




########### From Barkley et al ############
########### https://www.nature.com/articles/s41588-022-01141-9 #######
# Modules to cells
GeneToEnrichment = function(
    srt,
    type = 'GO',
    db = NULL,
    method = 'rand',
    genes = NULL,
    assay = NULL,
    do.rescale = FALSE,
    min.cells = 0,
    min.genes = 0,
    min.var = 0,
    min.var.rescaled = 0,
    auc_percentile = 0.05,
    db_rand = NULL,
    nrand = 4,
    nbin = 25,
    rand_method='rapid_original',
    addtoseu=TRUE,
    ...
){
  if (is.null(assay)){
    assay = DefaultAssay(srt)
  }
  if (is.null(db)){
    db = FindMSigDB(type)
  }
  
  counts = as.matrix(GetData(srt, assay = assay, slot = 'counts'))
  genes = rownames(counts)
  genes.expr = rownames(counts)[rowSums(counts) > min.cells]
  
  if (method == 'metagene'){
    
    data = as.matrix(GetAssayData(srt, assay = assay, slot = 'scale.data'))
    
    db = lapply(db, intersect, genes.expr)
    
    enrichment.profile = t(sapply(names(db), function(m){
      colMeans(data[db[[m]], ], na.rm = TRUE)
    }))
    
    enrichment.profile = enrichment.profile[sapply(names(db), function(x){
      v = var(enrichment.profile[x, ])
      l = length(db[[x]])
      return(l > min.genes
             && v > min.var
             && v*l^2 > min.var.rescaled)
    }), ]
    
    if (do.rescale){
      mn = apply(enrichment.profile, 1, mean)
      v = apply(enrichment.profile, 1, var)
      enrichment.profile = (enrichment.profile - mn) / sqrt(v)
    }
    
    srt = AddMetaData(srt,  t(enrichment.profile), col.name = rownames(enrichment.profile))
  }
  
  if (method == 'auc'){
    
    data = as.matrix(GetData(srt, assay = assay, slot = 'data'))
    
    cells_rankings = AUCell_buildRankings(data)
    cells_AUC = AUCell_calcAUC(db, cells_rankings, aucMaxRank=nrow(cells_rankings)*auc_percentile)
    enrichment.profile = getAUC(cells_AUC)
    
    if (do.rescale){
      mn = apply(enrichment.profile, 1, mean)
      v = apply(enrichment.profile, 1, var)
      enrichment.profile = (enrichment.profile - mn) / sqrt(v)
    }
    
    srt = AddMetaData(srt,  t(enrichment.profile), col.name = rownames(enrichment.profile))
  }
  
  if (method == 'score'){
    
    temp = AddModuleScore(srt, features = db, assay = assay, name = names(db), nbin = nbin, ...)
    
    enrichment.profile = t(temp@meta.data[, names(db)])
    
    if (do.rescale){
      mn = apply(enrichment.profile, 1, mean)
      v = apply(enrichment.profile, 1, var)
      enrichment.profile = (enrichment.profile - mn) / sqrt(v)
    }
    
    srt = AddMetaData(srt,  t(enrichment.profile), col.name = rownames(enrichment.profile))
  }
  
  if (method == 'rand_distribution'){
    data = as.matrix(GetData(srt, assay = assay, slot = 'scale.data'))
    
    db = lapply(db, intersect, genes)
    
    data.avg = sort(rowMeans(x = data))
    data.cut = cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                          n = nbin, labels = FALSE, right = FALSE)
    names(data.cut) = names(data.avg)
    binned = split(data.avg, data.cut)
    ecdf_binned <- sapply(binned, ecdf)
    mean_binned <- sapply(binned, mean)
    
    
    pseudoval <- 1*10^-5
    enrichment.profile = t(sapply(names(db), function(m){
      print(m)
      pools <- data.cut[db[[m]]]
      pools_ecdf <- lapply(pools, function(pool_i){
        ecdf_binned[[pool_i]]
      })
      re <- t(sapply(seq_along(pools_ecdf), function(g_idx){
        pools_ecdf[[g_idx]](data[db[[m]][g_idx],])
      })) %>%
        apply(., 2, function(i) {
          metap::sumlog(1-na.omit(i))$p
          # mean(1-i)
          })
      print("re-generated")
      p = -log10(re + pseudoval)
      return(p)
    }))
    
    
    l2r.profile = t(sapply(names(db), function(m){
      pools <- data.cut[db[[m]]]
      pools_mean <- lapply(pools, function(pool_i){
        mean_binned[[pool_i]]
      })
      l2r <- colMeans(t(sapply(seq_along(pools_mean), function(g_idx){
        log2((data[db[[m]][g_idx],] / (pools_mean[[g_idx]] + pseudoval)) + 1)
      })))
      return(l2r)
    }))
    
    if(addtoseu){
      srt = AddMetaData(srt,  t(enrichment.profile), col.name = rownames(enrichment.profile))
    } else {
      srt = list("p"=as.data.frame(t(enrichment.profile)) %>%
                   magrittr::set_rownames(., Cells(srt)),
                 "l2r"=t(l2r.profile))
    }
  }
  
  if (method == 'rand'){
    
    data = as.matrix(GetData(srt, assay = assay, slot = 'scale.data'))
    
    db = lapply(db, intersect, genes)
    
    if (is.null(db_rand)){
      db_rand = MakeRand(srt, db, nrand = nrand, nbin = nbin)
    } else {
      nrand = log10(length(db_rand[[1]]))
    }
    
    enrichment.profile = t(sapply(names(db), function(m){
      print(m)
      if(rand_method == 'rapid_original'){
          # d1 <- Sys.time()
          genes_pop <- unique(unlist(db_rand[[m]]))
          data_pop <- round(data[genes_pop,],3)
          idx <- sapply(db_rand[[m]], match, rownames(data_pop))
          ra = apply(idx, 2, function(i){
            colMeans(data[i, ], na.rm = TRUE)
          })
          # print(Sys.time() - d1)
        } else if(rand_method == 'original'){
          # d1 <- Sys.time()
          ra = sapply(db_rand[[m]], function(i){
            colMeans(data[i, ], na.rm = TRUE)
          })
          # print(Sys.time() - d1)
        }
      
      
      print("ra-generated")
      re = colMeans(data[db[[m]], ], na.rm = TRUE)
      print("re-generated")
      p = rowMeans(ra >= re)
      p = -log10(p)
      return(p)
    }))
    enrichment.profile[is.infinite(enrichment.profile)] = nrand
    enrichment.profile = enrichment.profile/nrand
    
    srt = AddMetaData(srt,  t(enrichment.profile), col.name = rownames(enrichment.profile))
  }
  
  return(srt)
}


# Modified from Seurat
AddModuleScore = function (object, features, pool = NULL, nbin = 24, ctrl = 100, 
                           k = FALSE, assay = NULL, name = "Cluster", seed = 1) 
{
  set.seed(seed = seed)
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay
  assay.data <- GetData(object)
  features.old <- features
  if (k) {
    .NotYetUsed(arg = "k")
    features <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      features[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == 
                                         i))
    }
    cluster.length <- length(x = features)
  } else {
    if (is.null(x = features)) {
      stop("Missing input feature list")
    }
    features <- lapply(X = features, FUN = function(x) {
      return(intersect(x = x, y = rownames(x = object)))
    })
    cluster.length <- length(x = features)
  }
  if (!all(Seurat:::LengthCheck(values = features))) {
    warning(paste("Could not find enough features in the object from the following feature lists:", 
                  paste(names(x = which(x = !Seurat:::LengthCheck(values = features)))), 
                  "Attempting to match case..."))
    features <- lapply(X = features.old, FUN = CaseMatch, 
                       match = rownames(x = object))
  }
  if (!all(Seurat:::LengthCheck(values = features))) {
    stop(paste("The following feature lists do not have enough features present in the object:", 
               paste(names(x = which(x = !Seurat:::LengthCheck(values = features)))), 
               "exiting..."))
  }
  pool <- pool %||% rownames(x = object)
  data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                         n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(ctrl.use[[i]], 
                         names(x = sample(data.cut[which(data.cut == data.cut[features.use[j]])], size=ctrl, replace=FALSE)))
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), 
                        ncol = ncol(x = object))
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, ])
  }
  features.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length, 
                            ncol = ncol(x = object))
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- assay.data[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = object)
  object[[colnames(x = features.scores.use)]] <- features.scores.use
  Seurat:::CheckGC()
  DefaultAssay(object = object) <- assay.old
  return(object)
}

# To get full expression matrix
GetData = function(
    srt,
    genes = NULL,
    slot = 'scale.data',
    assay = NULL
){
  if (is.null(assay)){
    assay = DefaultAssay(srt)
  }
  if (is.null(genes)){
    if ('RNA' %in% names(srt@assays)){
      genes = rownames(GetAssayData(srt, assay = 'RNA', slot='counts'))
    } else if  ('Spatial' %in% names(srt@assays)){
      genes = rownames(LayerData(srt, assay = 'Spatial', layer = 'counts'))
    } else {
      genes = rownames(LayerData(srt, assay = assay, layer = 'counts'))
    }
  }
  data = GetAssayData(srt, assay = assay, layer = slot)
  missing = setdiff(genes, rownames(data))
  add = matrix(0, nrow = length(missing), ncol = ncol(data))
  rownames(add) = missing
  data = rbind(data, add)
  data = data[genes, ]
  return(data)
}

# Make equivalent random modules
MakeRand = function(
    srt,
    db,
    assay = NULL,
    nrand = 3, #nrand is the exponent of 10, number of random sets (nrand 3 -> 10^3 random sets)
    nbin = 25
){
  if (is.null(assay)){
    assay = DefaultAssay(srt)
  }
  data = GetData(srt, slot='data')
  db = lapply(db, intersect, rownames(data))
  data.avg = sort(rowMeans(x = data))
  for(nbin_i in c(2:100)){
    data.cut = cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
                          n = nbin_i, labels = FALSE, right = FALSE)
    names(data.cut) = names(data.avg)
    binned = split(names(data.cut), data.cut)
    date()
    db_rand = lapply(names(db), function(m){
      dup_check <- TRUE
      pool <- data.cut[db[[m]]] # get the pool for each gene
      rand_sel_pool <- sapply(pool, function(p) sample(x=binned[[p]], size=10^nrand, replace=T))
      
      counter <- 1
      max_iterations=50
      while(dup_check){
        dup_idx <- apply(rand_sel_pool, 1, function(i) any(duplicated(i)))
        if(counter == max_iterations){
          print(paste0("Max-iterations hit, forcing non-duplicates for ", m))
          # If no solution could be reached after x number of iterations, force a solution by 
          # going through each random vector and replace only the duplicated entries
          dup_check <- FALSE
          used <- lapply(seq_along(dup_idx), function(dup_i){
            if(!dup_idx[[dup_i]]){
              nonrepeat_rand <- rand_sel_pool[dup_i,]
            } else {
              xrand <- rand_sel_pool[dup_i,]
              force_dup_check <- TRUE
              force_counter <- 1
              while(force_dup_check){
                # Replace only the duplicated enrtices
                dup_xrand <- duplicated(xrand)
                newval <- sapply(pool[which(dup_xrand)], function(p) sample(x=binned[[p]], size=1, replace=T))
                xrand[which(dup_xrand)] <- as.character(newval)
                
                force_counter <- force_counter + 1
                if(force_counter == max_iterations  |!any(dup_xrand)){
                  force_dup_check <- FALSE
                }
              }
              nonrepeat_rand <- xrand
            }
            return(nonrepeat_rand)
          })
        }
        if(!any(dup_idx)) {
          dup_check <- FALSE
          rand_sel_pool <- rand_sel_pool[which(!dup_idx),]
          used <- lapply(1:nrow(rand_sel_pool), function(i) as.character(rand_sel_pool[i,]))
        } else {
          counter <- counter + 1
          rand_sel_pool <- rand_sel_pool[-which(dup_idx),]
          rand_sel_pool <- rbind(rand_sel_pool,
                                 sapply(pool, function(p) sample(x=binned[[p]], size=sum(dup_idx), replace=T)))
        }
      }
      
      lapply(used, function(i){
        
      })
      return(used)
    })
    date()
    names(db_rand) = names(db)
  }
  return(db_rand)
}


visualize_scores = function(score_colnames, seurs_ls, pdf_width, pdf_height, suffix, out_dir, limits=NULL) {
  # score_colnames (character) names of the columns in the meta.data that contain the scores
  # suffix (character) choose a suffix for the pdf name to indicate which score was visualized
  #                    e.g. 'ssGSEA', 'UCell', 'Barkley'
  # limits (numeric) vector of 2 numbers delimiting the range of color values for the color legend
  dot_size = rep(2, sum(sapply(seurs_ls, function(seur) length(unique(seur$sample)))))
  names(dot_size) = unlist(sapply(seurs_ls, function(seur) unique(seur$sample)))
  dot_size['WholeLiver_Guilliams_H35_S2'] = dot_size['WholeLiver_Guilliams_H38'] = 3
  dot_size['WholeLiver_Guilliams_H37'] = 3
  dot_size['Lung_Kadur_D1'] = 3
  dot_size['cSCC_Ji4.rep1'] = dot_size['cSCC_Ji4.rep2'] = dot_size['cSCC_Ji6.rep1'] = dot_size['cSCC_Ji6.rep2'] = 7
  dot_size['HealthyLiver_Guilliams_Capsule_M2_S1'] = dot_size['HealthyLiver_Guilliams_M1_S1'] = 3
  dot_size['HealthyLiver_Guilliams_M1_S2'] = dot_size['HealthyLiver_Guilliams_M1_S3'] = 2.5
  dot_size['Liver_Std_diet_Guilliams_M1_S1'] = 3
  
  SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
  
  for (i in 1:length(seurs_ls)) { 
    seurs_i = seurs_ls[[i]]
    samples = unique(seurs_i$sample)

    pdf(file.path(out_dir, paste0(suffix, '_', names(seurs_ls)[i], '_spatial_expression.pdf')), 
        width=pdf_width, height=pdf_height) 
    for (s in 1:length(samples)) {
      sample_s = samples[s]
      p = SpatialFeaturePlot(subset(seurs_i, sample==sample_s), 
                             features=c(score_colnames, 'nFeature_Spatial'), 
                             pt.size.factor=dot_size[sample_s], alpha=.7, image.alpha=0.5, 
                             images=paste0('image_', sample_s), combine = F) 
      for (i in 1:length(p)) {
        p[[i]] = p[[i]] +
          theme(legend.key.height= unit(.15, 'cm'), # .3
                legend.key.width= unit(.27, 'cm'), # .5
                legend.title = element_text(size=6), #8
                legend.text = element_text(size=5),
                legend.margin = unit(c(5.5, 5.5, 0, 5.5), "pt"),
                plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt"))
        if (!is.null(limits) & i != length(p)) {
          p[[i]] = p[[i]] +
            scale_fill_gradientn(limits=limits, colours=SpatialColors(n = 100))
        }
      }
      
      p = ggarrange(plotlist = p, nrow=2, ncol=4) +
        plot_annotation(title=sample_s, theme=theme(plot.title=element_text(size=7)))
      print(p)
    }
    dev.off()
  }
}




#### 03 calculate and visualize module score as in Barkley et al. #####
genesets_rand = MakeRand(seu, db=msig_l, nrand=2, nbin=25) 
for (g in names(msig_l)) {
  print(paste0(g, ": ", match(g, names(msig_l)), "/", length(msig_l[poi])))
  seu = GeneToEnrichment(seu, db=msig_l[g], method='rand', db_rand=genesets_rand[g], addtoseu=T)
}
cell_l <- split(Cells(seu), seu$manual_anno)

res = GeneToEnrichment(seu, db=msig_l, method='rand_distribution', addtoseu=FALSE)
wilcox_auc <- function(auc, seu, grps=c("CD4+ T cells", "B cells"), colid='manual_anno'){
  if(class(auc) != 'Matrix'){
    auc <- auc
  } else {
    auc <- assay(auc)
  }
  spl_auc <- split(as.data.frame(t(auc[,colnames(seu)])),
                   seu@meta.data[,colid])
  spl_auc <- spl_auc[grps]
  
  wilcoxp <- sapply(1:ncol(spl_auc[[1]]), function(idx){
    wilcox.test(spl_auc[[1]][,idx], spl_auc[[2]][,idx])$p.value
  })
  fc <- sapply(1:ncol(spl_auc[[1]]), function(idx){
    pc <- 0.0001
    log2(mean(spl_auc[[1]][,idx])+pc) - log2(mean(spl_auc[[2]][,idx])+pc)
  })
  
  auc_df <- data.frame("Geneset"=colnames(spl_auc[[1]]),
                       'FC'=fc,
                       "p"=wilcoxp,
                       "padj"=p.adjust(wilcoxp, method='BH')) %>%
    mutate(Dataset=gsub("_.*", "", Geneset)) %>%
    arrange(padj)
  return(auc_df)
}
wilcox_auc(auc=t(res$p), seu, grps=c("CD4+ T cells", "B cells"), colid='anno') %>% head(., 20)
wilcox_auc(auc=t(seu@meta.data[,grep("REACTOME", colnames(seu@meta.data))]), seu, grps=c("CD4+ T cells", "B cells"), colid='anno') %>% head(., 20)





res = seu@meta.data[, names(msig_l)]

for (i in 1:length(seurs_ls)) {
  seurs   = seurs_ls[[i]]
  seur_ls = SplitObject(seu, split.by = 'orig.ident')
  specie  = strsplit(names(seurs_ls)[i], '_')[[1]][1]
  genesets_ls = IFN_sign_ls[[specie]]
  
  mod_scores = do.call(rbind, lapply(seur_ls, function(seur) {
    is.cSSC_data = length(grep('cSCC_Ji', seur$sample[1]) > 0)
    genesets_rand = MakeRand(seu, db=msig_l, nrand=3, nbin=25) 
    for (g in names(msig_l)) {
      print(paste0("g: ", g))
      res = GeneToEnrichment(seu, db=msig_l[g], method='rand_distribution', addtoseu=FALSE)
    }
    res = seu@meta.data[, names(msig_l)]
  }))
  seurs_ls[[i]]@meta.data = cbind(seurs_ls[[i]]@meta.data, mod_scores)
  
  print(paste(names(seurs_ls)[i], 'Done!'))
}

## visulatize scores in pdf
score_colnames = names(IFN_sign_ls$Human)
visualize_scores(score_colnames, seurs_ls, pdf_width=10, pdf_height=6, suffix='Barkley_score', 
                 out_dir, limits=c(0, 1))

