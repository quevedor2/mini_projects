publicationDimPlot <- function(seu, grp, simplify_labels=TRUE, colseed=231, 
                               colors_use=NULL, pt.size=NULL, aspect.ratio=NULL,
                               repel=FALSE, reduction='umap', return='plot', ...){
  if(simplify_labels){
    seu@meta.data[,grp] <- gsub("^T_", "", seu@meta.data[,grp])
    grp_ids <- na.omit(unique(seu@meta.data[,grp]))
    grp_map <- setNames(seq_along(grp_ids), grp_ids)
    legend_map <- setNames(paste(grp_map, names(grp_map), sep="_"), grp_map)
    
    seu@meta.data[,grp] <- factor(grp_map[seu@meta.data[,grp]])
    labels <- levels(seu@meta.data[,grp])
  }
  
  if(is.null(pt.size)){
    pt.size <- scCustomize:::AutoPointSize_scCustom(data = ncol(seu), raster = T)
  }
  if(is.null(colors_use)){
    print("Grabbing random colors...")
    colors_use <- scCustomize::scCustomize_Palette(num_groups = length(unique(seu@meta.data[,grp]))+1, 
                                                   ggplot_default_colors = FALSE, 
                                                   color_seed = colseed)
    if('#E4E1E3FF' %in% colors_use){
      colors_use <- colors_use[-which(colors_use %in% '#E4E1E3FF')]
    } else {
      colors_use <- colors_use[c(1:length(unique(seu@meta.data[,grp])))]
    }
  }
  plot <- DimPlot(seu, group.by=grp, raster=T, label=T, reduction=reduction, 
                  repel=repel, cols = colors_use, pt.size=pt.size) +
    theme(...)
  if(simplify_labels){
    plot <- plot + 
      scale_color_manual(labels=legend_map[labels], values=colors_use) + 
      guides(color=guide_legend(ncol =1))
  }
  plot <- plot & NoAxes()
  
  axis_plot <- ggplot(data.frame(x= 100, y = 100), aes(x = .data[["x"]], y = .data[["y"]])) +
    geom_point() +
    xlim(c(0, 10)) + ylim(c(0, 10)) +
    theme_classic() +
    ylab(plot$labels$y) + xlab(plot$labels$x) +
    theme(plot.background = element_rect(fill = "transparent", colour = NA),
          panel.background = element_rect(fill = "transparent"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_line(
            arrow = arrow(angle = 15, length = unit(.5, "cm"), type = "closed"))
    )
  
  figure_layout <- c(
    patchwork::area(t = 1, l = 2, b = 11, r = 11),
    patchwork::area(t = 10, l = 1, b = 12, r = 2))
  
  if(return=='plot'){
    plot_figure <- plot + axis_plot +
      patchwork::plot_layout(design = figure_layout) & 
      theme(aspect.ratio=aspect.ratio)
  } else {
    plot_figure = list('plot'=plot,
                       'axis_plot'=axis_plot,
                       'figure_layout'=figure_layout,
                       'aspect.ratio'=aspect.ratio)
  }
  return(plot_figure)
}
