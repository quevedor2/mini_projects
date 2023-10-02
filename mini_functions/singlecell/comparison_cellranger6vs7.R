## Comparison between cellrangerv6 (intronRemoved) and cellrangerv7 (intronIncluded) analysis
renv::load("/cluster/home/quever/downloads/renvs/")
# Core
library(monocle3)
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
# Visualization
library(cowplot)
library(ggrastr)
library(ggplot2)


options(Seurat.object.assay.version = "v5")
PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/ipmn_pdac/data/rna/intron_testing'
setwd(PDIR)

###################
#### Functions ####

##############
#### Main ####
samples <- list.files(PDIR, pattern="intron")
seul <- lapply(samples, function(s){
  seu.data <- Read10X(data.dir = file.path(PDIR, s))
  seu <- CreateSeuratObject(counts = seu.data, project = s)
  return(seu)
})



seu <- merge(seul[[1]], seul[-1], add.cell.ids=samples)
DefaultAssay(seu) <- 'RNA'
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)


seu <- NormalizeData(seu) %>%
  FindVariableFeatures(.)  %>%
  ScaleData(.)  %>%
  RunPCA(.) %>% 
  FindNeighbors(., dims = 1:30, reduction = "pca")  %>%
  FindClusters(., resolution = 1.2, cluster.name = "unintegrated_clusters") %>%
  RunUMAP(., dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

DefaultAssay(seu) <- 'RNA'
seu_integ <- IntegrateLayers(
  object = seu, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = T
)
seu_integ <- seu_integ %>% 
  FindNeighbors(., reduction = "integrated.mnn", dims = 1:30) %>%
  FindClusters(., resolution = 0.9, cluster.name = "mnn_clusters") %>% 
  RunUMAP(., reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")

seu_integ <- JoinLayers(seu_integ)
seu <- JoinLayers(seu)


dir.create(file.path(PDIR, "seuratobj"), showWarnings = F)
saveRDS(seu_integ, file=file.path(PDIR, "seuratobj", "seu_integrated.rds"))
saveRDS(seu, file=file.path(PDIR, "seuratobj", "seu_nonintegrated.rds"))


dir.create(file.path(PDIR, "results"))
pdf(file.path(PDIR, "results", "compare.pdf"), width = 13)
dp1a <- DimPlot(seu_integ, group.by='orig.ident', raster=T, reduction='umap.mnn', pt.size=1.5)
dp1b <- DimPlot(seu_integ, group.by='seurat_clusters', raster=T, reduction='umap.mnn', pt.size=1.5, label=T)
dp2a <- DimPlot(seu, group.by='orig.ident', raster=T, reduction='umap.unintegrated', pt.size=1.5)
dp2b <- DimPlot(seu, group.by='seurat_clusters', raster=T, reduction='umap.unintegrated', pt.size=1.5, label=T)
cowplot::plot_grid(dp1a, dp1b, ncol=2)
cowplot::plot_grid(dp2a, dp2b, ncol=2)
dev.off()



DefaultAssay(seu_integ) <- 'RNA'
Idents(seu_integ) <- 'seurat_clusters'
clusters <- c('3', '6')
degs_l <- lapply(setNames(clusters,clusters), function(clid){
  seusub <- subset(seu_integ, idents=clid)
  
  Idents(seusub) <- 'orig.ident'
  degs1 <- FindMarkers(seusub, slot='data',
                       ident.1='IPMN_HPB197_intronPos', 
                       ident.2='IPMN_HPB197_intronNeg',
                       logfc.threshold=0)
  degs2 <- FindMarkers(seusub, slot='data',
                       ident.1='PDAC_HPB342_intronPos', 
                       ident.2='PDAC_HPB342_intronNeg',
                       logfc.threshold=0)
  return(list("IPMN"=degs1, "PDAC"=degs2))
})
lapply(unlist(degs_l, recursive = F), head)
saveRDS(degs_l, file.path(PDIR, "results", "degs.uninteg_clusters.rds"))
              
DefaultAssay(seu_integ) <- 'RNA'
Idents(seu_integ) <- 'seurat_clusters'
cell_barcodes <- split(gsub("^.*_([ACGT]*)-.*$", "\\1", Cells(seu_integ)), 
                       f=gsub("^(.*)_[ACGT]*-.*$", "\\1", Cells(seu_integ)))

samples <- c('IPMN_HPB197', 'PDAC_HPB342')
deg_samecells_l <- lapply(setNames(samples,samples), function(s){
  print(s)
  cells <- sample(cell_barcodes[[paste0(s, "_intronNeg")]], size=500)
  cells_keep <- sapply(cells, function(i) grep(i, Cells(seu_integ), value=T)) %>% 
    unlist() %>% as.character
  
  Idents(seu_integ) <- 'orig.ident'
  FindMarkers(seu_integ, 
              ident.1=paste0(s, "_intronPos"),
              ident.2=paste0(s, "_intronNeg"), 
              logfc.threshold=0)
})
saveRDS(deg_samecells_l, file.path(PDIR, "results", "degs.500matching_barcodes.rds"))
lapply(deg_samecells_l, head)

pcutoff <- 1*10^(-c(1:10))
sapply(setNames(pcutoff, pcutoff), function(cutoff){
  genes <- lapply(deg_samecells_l, function(df){
    df %>% as.data.frame %>%
      filter(p_val_adj < cutoff) %>%
      tibble::rownames_to_column("gene") %>%
      pull(gene)
  })
  c(length(intersect(genes[[1]], genes[[2]])), length(unique(unlist(genes))))
  
}) %>%
  t %>%  as.data.frame %>%
  magrittr::set_colnames(c('n_overlap', 'n_uniq_sig_genes')) %>%
  mutate(frac=round(n_overlap/n_uniq_sig_genes, 2))
