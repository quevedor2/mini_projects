## This function was editted from the standard Run.STACAS() function add in
## (...) to the main function, allowing it to be passeded into the 
## STACAS:::FindAnchors.STACAS() function. The poitn of this is to allow the
## user to pass in the future.maxSize parameter, allowing the user to 
## increase the maximum memory size from 16gb to anything they wish
## options(future.maxSize)
Run.STACAS <- function (object.list = NULL, assay = NULL, reference = NULL, 
          anchor.features = 1000, genesBlockList = "default", dims = 1:30, 
          normalization.method = c("LogNormalize", "SCT"), k.anchor = 5, 
          k.score = 30, k.weight = 100, alpha = 0.8, anchor.coverage = 0.5, 
          correction.scale = 2, cell.labels = NULL, label.confidence = 1, 
          seed = 123, verbose = FALSE, ...) {
  stacas_anchors <- FindAnchors.STACAS(object.list, assay = assay, 
                                       reference = reference, anchor.features = anchor.features, 
                                       genesBlockList = genesBlockList, dims = dims, normalization.method = normalization.method, 
                                       k.anchor = k.anchor, k.score = k.score, alpha = alpha, 
                                       anchor.coverage = anchor.coverage, correction.scale = correction.scale, 
                                       cell.labels = cell.labels, seed = seed, label.confidence = label.confidence, 
                                       verbose = verbose, ...)
  if (is.null(cell.labels)) {
    semisupervised <- FALSE
  } else {
    semisupervised <- TRUE
  }
  tree <- SampleTree.STACAS(anchorset = stacas_anchors, semisupervised = semisupervised, 
                            plot = FALSE)
  integrated <- IntegrateData.STACAS(stacas_anchors, dims = dims, 
                                     sample.tree = tree, k.weight = k.weight, semisupervised = semisupervised, 
                                     normalization.method = normalization.method, features.to.integrate = stacas_anchors@anchor.features)
  normalization.method <- match.arg(arg = normalization.method)
  if (normalization.method == "LogNormalize") {
    integrated <- ScaleData(integrated)
  }
  integrated <- RunPCA(integrated, npcs = max(dims))
  return(integrated)
}

## This function was editted from the standard ProjecTILs make.reference() 
## function to change line 66. Instead of using the ProjecTILs:::prcomp.seurat()
## function, which may through memory errors, it will use the Seurat
## ScaleData() %>% RunPCA() function to achieve virtually the same goal of 
## running PCA on the integrated assay.
make.reference <- function (ref, assay = NULL, atlas.name = "custom_reference", 
                            annotation.column = "functional.cluster", recalculate.umap = FALSE, 
                            umap.method = c("umap", "uwot"), metric = "cosine", min_dist = 0.3, 
                            n_neighbors = 30, ndim = 20, dimred = "umap", nfeatures = 1000, 
                            seed = 123){
  if (is.null(assay)) {
    assay = DefaultAssay(ref)
  }
  
  if (is.null(ref@assays[[assay]])) {
    stop(sprintf("Assay %s not found in reference object. Select a different assay", 
                 assay))
  }
  
  if ("var.features" %in% slotNames(ref@assays[[assay]]) & 
      !is.null(ref@assays[[assay]]@var.features)) {
    varfeat <- ref@assays[[assay]]@var.features
  } else {
    ref <- FindVariableFeatures(ref, assay = assay, nfeatures = nfeatures, 
                                verbose = FALSE)
    varfeat <- ref@assays[[assay]]@var.features
  }
  
  ref <- ScaleData(ref, features=varfeat) %>% 
    RunPCA(., features=varfeat, npcs = ndim, assay = assay, verbose = FALSE)
  ref <- .makeprcomp_obj(ref, varfeat)
  pca.obj <- ref@misc$pca_object
  table(names(pca.obj$center) %in% rownames(pca.obj$rotation))
  
  # ref <- ProjecTILs:::prcomp.seurat(ref, ndim = ndim, assay = assay)
  if (!recalculate.umap) {
    if (dimred %in% names(ref@reductions)) {
      ref.pca <- ref@misc$pca_object
      cell.order = rownames(ref.pca$x)
      low.emb <- ref@reductions[[dimred]]@cell.embeddings[cell.order, ]
      colnames(low.emb) <- c("UMAP_1", "UMAP_2")
      ref@misc$umap_object <- list()
      ref@misc$umap_object$data <- ref.pca$x
      ref@misc$umap_object$layout <- low.emb
    } else {
      stop(sprintf("Dimred %s not found in reference object. Select a different dimensionality reduction, or set recalculate.umap=TRUE to compute UMAP coordinates", 
                   dimred))
    }
  } else {
    umap.method = umap.method[1]
    ref.pca <- ref@misc$pca_object
    if (umap.method == "umap") {
      ref.umap <- run.umap.2(ref.pca, ndim = ndim, seed = seed, 
                             n.neighbors = n_neighbors, min.dist = min_dist, 
                             metric = metric)
      ref@misc$umap_object <- ref.umap
      ref@reductions$umap@cell.embeddings <- ref.umap$layout
    } else if (umap.method == "uwot") {
      warning("There are known issues with saving a loading uwot models. If you plan to save your reference as an .rds file, please use umap.method='umap'")
      ref.umap <- run.umap.uwot(ref.pca, ndim = ndim, seed = seed, 
                                n.neighbors = n_neighbors, min.dist = min_dist, 
                                metric = metric)
      ref.umap$data <- ref.pca$x
      ref@misc$umap_object <- ref.umap
      ref@reductions$umap@cell.embeddings <- ref.umap$embedding
    } else {
      stop("Unsupported UMAP method.")
    }
  }
  
  ## ref <- RenameAssays(object = ref, mnn.reconstructed = 'integrated')
  # names(ref@assays)[names(ref@assays) == assay] = "integrated"
  # DefaultAssay(ref) <- "integrated"
  if (!annotation.column == "functional.cluster") {
    ref$functional.cluster <- ref@meta.data[, annotation.column]
  }
  ref$functional.cluster <- factor(ref$functional.cluster)
  Idents(ref) <- "functional.cluster"
  ref@misc$projecTILs = atlas.name
  return(ref)
}


## The STACAS object requires a prcomp class object in the obj@misc$pca_object 
## slot for ProjecTILs integration. This object is not generated in the 
## RunPCA() function of seurat, but it can be recreated manually. Based
## on the prcomp() code https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/prcomp.R,
## I recreated  a prcomp object using the values extracted from the 
## RunPCA() seurat function, and then use the base::stats scale() function
## to estimate the scale and cen values
.makeprcomp_obj <- function(obj, varfeat){
  assay <- DefaultAssay(obj)
  # varfeat <- VariableFeatures(obj)
  scx <- data.frame(t(as.matrix(obj@assays[[assay]]@data[varfeat, ])),
                    check.names = F) %>%
    scale(., center=T, scale=T)
  # refdata <- refdata[, sort(colnames(refdata))]
  # scx <- scale(refdata, center = TRUE, scale = TRUE)
  cen <- attr(scx, "scaled:center")
  sc <- attr(scx, "scaled:scale")
  
  r <- list(x=obj@reductions$pca@cell.embeddings,
            sdev=obj@reductions$pca@stdev, #sdev = s$d, 
            rotation=obj@reductions$pca@feature.loadings, #rotation = s$v,
            center = cen, #if(is.null(cen)) FALSE else cen,
            scale = sc) #if(is.null(sc)) FALSE else sc)
  class(r) <- "prcomp"
  obj@misc$pca_object <- r
  return(obj)
}

# prcomp.seurat <- function (obj, assay = NULL, ndim = 10, scale = TRUE)  { 
#   if (is.null(assay)) {
#     assay <- DefaultAssay(obj)
#   }
#   varfeat <- VariableFeatures(obj)
#   refdata <- data.frame(t(as.matrix(obj@assays[[assay]]@data[varfeat, ])))
#   refdata <- refdata[, sort(colnames(refdata))]
#   ref.pca <- prcomp(refdata, rank. = ndim, scale. = scale, 
#                     center = TRUE, retx = TRUE)
#   obj@misc$pca_object <- ref.pca
#   obj@reductions$pca@cell.embeddings <- ref.pca$x
#   obj@reductions$pca@feature.loadings <- ref.pca$rotation
#   colnames(obj@reductions$pca@cell.embeddings) <- gsub("PC(\\d+)", "PC_\\1", colnames(ref.pca$x), perl = TRUE)
#   colnames(obj@reductions$pca@feature.loadings) <- gsub("PC(\\d+)", "PC_\\1", colnames(ref.pca$rotation), perl = TRUE)
#   return(obj)
# }