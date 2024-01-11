# # DietSeurat fix to remove layers
# object = seusub
# layers = grep("count", Layers(object), value=T)
# features = NULL
# assays = NULL
# dimreducs = NULL
# graphs = NULL
# misc = TRUE
# counts = lifecycle:::deprecated()
# data = lifecycle:::deprecated()
# scale.data = lifecycle:::deprecated()

DietSeurat <- function(
  object,
  layers = NULL,
  features = NULL,
  assays = NULL,
  dimreducs = NULL,
  graphs = NULL,
  misc = TRUE,
  counts = lifecycle:::deprecated(),
  data = lifecycle:::deprecated(),
  scale.data = lifecycle:::deprecated(),
  ...
) {
  CheckDots(...)
  dep.args <- c(counts = counts, data = data, scale.data = scale.data)
  for (lyr in names(x = dep.args)) {
    if (lifecycle:::is_present(arg = dep.args[[lyr]])) {
      if (is.null(x = layers)) {
        layers <- unique(x = unlist(x = lapply(
          X = Seurat::Assays(object = object),
          FUN = function(x) {
            return(Layers(object = object[[x]]))
          }
        )))
      }
      deprecate_soft(
        when = '5.0.0',
        what = paste0('DietSeurat(', lyr, ' = )'),
        with = 'DietSeurat(layers = )'
      )
      layers <- if (isTRUE(x = dep.args[[lyr]])) {
        c(layers, lyr)
      } else {
        Filter(f = function(x) x != lyr, x = layers)
      }
    }
  }
  object <- UpdateSlots(object = object)
  assays <- assays %||% Seurat::Assays(object = object)
  assays <- intersect(x = assays, y = Seurat::Assays(object = object))
  if (!length(x = assays)) {
    abort(message = "No assays provided were found in the Seurat object")
  }
  if (!DefaultAssay(object = object) %in% assays) {
    abort(
      message = "The default assay is slated to be removed, please change the default assay"
    )
  }
  layers <- layers %||% assays
  layers <- .PropagateList(x = layers, names = assays)
  for (assay in names(x = layers)) {
    layers[[assay]] <- tryCatch(
      expr = Layers(object = object[[assay]], search = layers[[assay]]),
      error = function(...) {
        return(character(length = 0L))
      }
    )
  }
  layers <- Filter(f = length, x = layers)
  if (!length(x = layers)) {
    abort(message = "None of the requested layers found")
  }
  for (assay in Seurat::Assays(object = object)) {
    if (!(assay %in% assays)) {
      object[[assay]] <- NULL
      next
    }
    layers.rm <- setdiff(
      x = Layers(object = object[[assay]]),
      y = layers[[assay]]
    )
    if (length(x = layers.rm)) {
      if (inherits(x = object[[assay]], what = 'Assay') && all(c('counts', 'data') %in% layers.rm)) {
        abort(message = "Cannot remove both 'counts' and 'data' from v3 Assays")
      }
      for (lyr in layers.rm) {
        suppressWarnings(object <- tryCatch(expr = {
          object[[assay]][lyr] <- NULL
          object
        }, error = function(e) {
          if (lyr == "data"){
            object[[assay]][lyr] <- Matrix::sparseMatrix(i = 1, j = 1, x = 1,
                                                   dims = dim(object[[assay]][lyr]),
                                                   dimnames = dimnames(object[[assay]][lyr]))
          } else{
            slot(object = object[[assay]], name = lyr) <- new(Class = "dgCMatrix")
          }
          message("Converting layer ", lyr, " in assay ",
                  assay, " to empty dgCMatrix")
          object
        }))
      }
    }
    if (!is.null(x = features)) {
      features.assay <- intersect(
        x = features,
        y = rownames(x = object[[assay]])
      )
      if (!length(x = features.assay)) {
        warn(message = paste0(
          'No features found in assay ',
          sQuote(x = assay),
          ', removing...'
        ))
        object[[assay]] <- NULL
        next
      }
      suppressWarnings(object[[assay]] <- subset(x = object[[assay]], features = features.assay))
    }
  }
  # remove misc when desired
  if (!isTRUE(x = misc)) {
    slot(object = object, name = "misc") <- list()
  }
  # remove unspecified DimReducs and Graphs
  all.objects <- .FilterObjects(
    object = object,
    classes.keep = c('DimReduc', 'Graph')
  )
  objects.to.remove <- all.objects[!all.objects %in% c(dimreducs, graphs)]
  for (ob in objects.to.remove) {
    object[[ob]] <- NULL
  }
  cells.keep <- list()
  for (assay in  Seurat::Assays(object = object)) {
    cells.keep[[assay]] <- colnames(x = object[[assay]] )
  }
  cells.keep <- intersect(colnames(x = object), unlist(cells.keep))
  if (length(cells.keep) <- ncol(x = object)) {
    object <- subset(object, cells = cells.keep)
  }
  return(object)
}