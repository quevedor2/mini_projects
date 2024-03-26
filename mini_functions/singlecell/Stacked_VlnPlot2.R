# The same as Stacked_VlnPlot from scCustomize, but will return the ggplot object
# from VlnPlot instead of patchworking things together
Stacked_VlnPlot2 <- function(
  seurat_object,
  features,
  group.by = NULL,
  split.by = NULL,
  idents = NULL,
  x_lab_rotate = FALSE,
  plot_legend = FALSE,
  colors_use = NULL,
  color_seed = 123,
  ggplot_default_colors = FALSE,
  plot_spacing = 0.15,
  spacing_unit = "cm",
  vln_linewidth = NULL,
  pt.size = NULL,
  raster = NULL,
  add.noise = TRUE,
  ...
) {
  # Check Seurat
  scCustomize:::Is_Seurat(seurat_object = seurat_object)
  
  # Check features and meta to determine which features present
  all_found_features <- scCustomize:::Feature_PreCheck(object = seurat_object, features = features)
  
  # set pt.size (default is no points)
  if (is.null(x = pt.size)) {
    pt.size <- 0
  }
  
  # Set rasterization
  num_cells <- unlist(x = CellsByIdentities(object = seurat_object, idents = idents))
  
  if (length(x = num_cells) * length(x = all_found_features) > 100000 && is.null(x = raster) && pt.size != 0) {
    raster <- TRUE
    cli_inform(message = c("NOTE: Rasterizing points since total number of points across all plots exceeds 100,000.",
                           "i" = "To plot in vector form set {.code raster=FALSE}")
    )
  }
  
  # Set default color palette based on number of levels being plotted
  if (is.null(x = group.by)) {
    group_by_length <- length(x = unique(x = seurat_object@active.ident))
  } else {
    group_by_length <- length(x = unique(x = seurat_object@meta.data[[group.by]]))
  }
  
  # Check colors use vs. ggplot2 color scale
  if (!is.null(x = colors_use) && ggplot_default_colors) {
    cli_abort(message = "Cannot provide both custom palette to {.code colors_use} and specify {.code ggplot_default_colors = TRUE}.")
  }
  if (is.null(x = colors_use)) {
    # set default plot colors
    if (is.null(x = colors_use)) {
      colors_use <- scCustomize_Palette(num_groups = group_by_length, 
                                        ggplot_default_colors = ggplot_default_colors, 
                                        color_seed = color_seed)
    }
  }
  
  # Setup plot spacing/margin parameter
  plot_margin <- margin(t = plot_spacing, r = 0, b = plot_spacing, l = 0, unit = spacing_unit)
  
  # Create plots
  plot_list <- purrr::map(all_found_features, function(x) scCustomize:::Modify_VlnPlot(seurat_object = seurat_object, 
                                                                                features = x, cols = colors_use, 
                                                                                group.by = group.by, 
                                                                                split.by = split.by, 
                                                                                idents = idents, 
                                                                                plot_margin = plot_margin, 
                                                                                pt.size = pt.size, 
                                                                                raster = raster, 
                                                                                add.noise = add.noise, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  # Add ability to rotate the X axis labels to the function call
  if (isTRUE(x = x_lab_rotate) || x_lab_rotate == 45) {
    plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.ticks.x = element_line())
  }
  
  if (x_lab_rotate == 90) {
    plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.ticks.x = element_line())
  }
  
  if (!is.logical(x = x_lab_rotate) && !x_lab_rotate %in% c(45, 90)) {
    cli_abort(message = "{.code x_lab_rotate} must be either a logical or a numeric value of 45 or 90.")
  }
  
  plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(), axis.ticks.x = element_line())
  return(plot_list)
  
}