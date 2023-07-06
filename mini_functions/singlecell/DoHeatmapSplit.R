# Similar to DoHeatmap() command, but added a split.by parameter to allow for splitting
DoHeatmapSplit <- function(seu, split.by='orig.ident', main_title='NULL', verbose=F, ...){
  uids <- unique(as.character(seu@meta.data[,split.by]))
  ref_hm <- DoHeatmap(seu, combine=FALSE, ...)[[1]] + 
    theme(legend.box.margin = margin(0, 0, 0, 12))
  legend <- get_legend(ref_hm)
  yaxis <- get_y_axis(ref_hm, position='left')
  gg_main_title <- ggdraw() +
    draw_label(label = main_title)
  
  dps <- lapply(uids, function(id, ...){
    cl <- setNames(as.character(seu@meta.data[,split.by]), Cells(seu))
    cells_sub <- names(which(cl==id))
    if(verbose) {
      print(paste0("Length: ", length(cells_sub), 
                 "; [", paste(head(cells_sub, 3), collapse=", "), "...]"))
    }
    seusub <- subset(seu, cells=cells_sub)
    hm <- DoHeatmap(seusub, combine=FALSE, label=FALSE, ...)[[1]] + 
      ggtitle(id) +
      NoLegend() + 
      theme(axis.text.y.left=element_blank()) 
    return(hm)
  }, ...)
  pgrid <- cowplot::plot_grid(plotlist=dps, nrow=1)
  plot_gg <- plot_grid(yaxis, pgrid, legend, rel_widths = c(1, length(dps), 1), nrow=1)
  plot_grid(gg_main_title, plot_gg, ncol=1, rel_heights=c(0.1, 1))
}

# pdf("~/xfer/test.pdf", width = 15)
# DoHeatmapSplit(seu_j, split.by='manual_clusters',
#                features=head(VariableFeatures(seu_j)), group.by='orig.ident')
# dev.off()
# pdf("~/xfer/test.pdf")
# j= DimPlot(seu_j, group.by='seurat_clusters', reduction='umap', label=T, combine=F)[[1]]
# l <- get_legend(j)
# l
# dev.off()

# hm_obj[[1]] + facet_grid(.~Sample, space = "fixed", scales="free")
# dev.off()