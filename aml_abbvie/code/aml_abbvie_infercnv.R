# renv::load("/cluster/projects/mcgahalab/envs/renvs/seuratv5_v2")
library(Seurat)
library(infercnv)

load("/cluster/home/quever/git/infercnv/data/infercnv_genes_example.rda")
seu <- readRDS("/cluster/projects/mcgahalab/data/mcgahalab/aml_abbvie/results/seurat_obj/seu_anno.rds")
seu@meta.data$batch = seu$orig.ident
outcnvs <- "/cluster/projects/mcgahalab/data/mcgahalab/aml_abbvie/results/infercnv"
samples <- unique(seu$batch)
organism = 'homo sapiens'
setwd(outcnvs)

anno <- read.csv("/cluster/projects/mcgahalab/data/mcgahalab/aml_abbvie/ref/Sara_anno.csv")
annomap <- setNames(anno$Sara_anno, anno$Barcode)
seu$anno <- annomap[Cells(seu)]


infercnvdat <- lapply(rev(setNames(samples, samples)), function(s){
  message(paste0(s, "....(", match(s, samples), "/", length(samples), ")\n"))
  outf=file.path(outcnvs, paste0(s, "_infercnv.rds"))
  if(file.exists(outf)){
    infercnv_obj2 <- readRDS(file=outf)
  } else {
    Idents(seu) <- 'batch'
    seusub <- subset(seu, ident=s)
    Idents(seusub) <- 'seurat_clusters'
    # seusub <- subset(seusub, downsample = 150)
    
    
    exprs <- as.matrix(seusub@assays$RNA$counts)
    idx <- which(rownames(exprs) %in% rownames(infercnv_genes_example))
    annos <- data.frame(seusub$anno)
    annos[grep("[BT] cells", annos[,1]),1] <- 'normal'
    annos[annos[,1]=='',1] <- 'misc'
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=exprs,
                                        annotations_file=data.frame(annos),
                                        gene_order_file=infercnv_genes_example,
                                        ref_group_names='normal')
    
    infercnv_obj2 = infercnv::run(infercnv_obj,
                                  cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                  out_dir=file.path(outcnvs, s),  # dir is auto-created for storing outputs
                                  cluster_by_groups=T,   # cluster
                                  denoise=T,
                                  HMM=T,
                                  resume_mode=FALSE)
    
    saveRDS(infercnv_obj2, file=outf)
  }
  return(infercnv_obj2)
})

icnvdirs <- list.files(outcnvs, pattern = "Sample[0-9]*$")
Idents(seu) <- 'orig.ident'
propcnvs <- lapply(icnvdirs, function(i){
  print(i)
  icnvdir=file.path(outcnvs, i)
  if(file.exists(file.path(icnvdir, 'run.final.infercnv_obj'))){
    infercnv_obj <- readRDS(file.path(icnvdir, 'run.final.infercnv_obj'))
    seu2 <- subset(seu, cells=colnames(infercnv_obj@expr.data))
    
    seu2 = infercnv::add_to_seurat(infercnv_output_path=icnvdir,
                                   seurat_obj=seu2, # optional
                                   top_n=10)
    colids <- colnames(seu2@meta.data) %>%
      grep("proportion_cnv", ., value=T)
    prop_cnv <- seu2@meta.data[,colids] %>% rowSums 
    prop_cnv <- cbind(seu2@meta.data[,colids], data.frame("total_cnv"=prop_cnv))
  } else {
    prop_cnv <- NULL
  }
  
  return(prop_cnv)
})
prop_cnv_map <- do.call(rbind, propcnvs)
seu@meta.data <- cbind(seu@meta.data, prop_cnv_map[Cells(seu),])


n_colors <- 100 # Adjust for smoother gradient
color_palette <- colorRampPalette(c("lightgrey", "red"))(n_colors)#(n_colors %/% 2 + 1)
# neg_cols <- colorRampPalette(c("blue", "lightgrey"))(n_colors %/% 2)
# color_palette <- c((neg_cols), pos_cols[-1]) # Combine and remove duplicate grey

gg <- scCustomize::FeaturePlot_scCustom(seu, features='total_cnv',
                                        reduction='umap.harmony', raster=T,
                                        colors_use = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
                                        na_cutoff = NA, pt.size=2) +
  labs(title='InferCNV prop_cnv', subtitle=paste0("red = high-proportions of CNVs"))
gg_chrs <- scCustomize::FeaturePlot_scCustom(seu, features=grep("proportion", colnames(prop_cnv_map), value=T),
                                             reduction='umap.harmony', raster=T,
                                             colors_use = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
                                             na_cutoff = NA, pt.size=2) 
pdf(file.path(outcnvs, paste0("infercnv_propcnv.pdf")), height = 13, width = 13)
gg_chrs + NoLegend()
gg
dev.off()
file.copy(file.path(outcnvs, paste0("infercnv_propcnv.pdf")), to="~/xfer", overwrite = T)