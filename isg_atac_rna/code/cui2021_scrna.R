library(Seurat)
library(ggplot2)
library(cowplot)
library(harmony)
library(dplyr)
# library(EnsDb.Hsapiens.v86)

samples <- setNames(c('GSM3701181_GP33_day30', 'GSM3701180_GP33_day8', 'day8-30'), c('day30', 'day8', 'day8-30'))
PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/nirmin_atac_ifn/'
mtxdir <- file.path(PDIR, 'data/Cui_2021/data/scRNA')
outdir <- file.path(PDIR, 'results/Cui_2021/scrna')
setwd(PDIR)

markers <- setNames(c('Cx3cr1', 'Tcf7', 'Il21', 'Slamf6', 'Pdcd1', 'Cd4', 'Cd8a', 'Gzmb', 
                      'Prdm1', 'Ikzf2', 'Havcr2', 'Cd274', 'Entpd1'), 
                    c("CX3CR1", "TCF-1", "IL-21", "Ly108", 'PD-1', 'CD4', 'CD8', 'GZMB', 
                      'PRDM1', 'IKZF2', 'TIM3', 'PDL1', 'CD39'))
# Table 1: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3274382/
isg_markers <- c('Adar', 'Apobec3', 'Bst2', 'Mb21d1', 'Cd74', 'Ddit4', 'Ddx58',
                 'Ddx60', 'Eif2ak2', 'Gbp2', 'Hpse', 'Ifi44l',
                 'Ifih1', 'Ifit1', 'Ifit2', 'Ifit3', 'Ifitm1',
                 'Ifitm2', 'Ifitm3', 'Irf1', 'Irf7', 'Isg15', 'Isg20', 'Map3k14',
                 'Mov10', 'Ms4a4a', 'Mx1', 'Nampt', 'Nt5c3',
                 'Oas2', 'Oas3', 'P2ry6', 'Phf1', 'Pml', 'Rsad2',
                 'Rtp4', 'Slc15a3', 'Slc25a28', 'Ssbp3', 'Trex1', 'Trim5', 
                 'Trim25', 'Sun2', 'Zc3hav1', 'Irf2')
''
# C6orf150 = CGAS = Cgas = Mb21d1
# GBP1 = no mouse ortholog
# IFI6 = no mouse ortholog
# G1P3 = synonym for IFI6
# IFIT5 = no mouse ortholog
# MX2 = Mx2 = not found
# OAS1 = no mouse ortholog
# OASL = no mouse ortholog
# PHF-15 = JADE2 = Phf15 = Phf1


###################
#### Functions ####
getHClustLabels <- function(seu, marker){
  hc <- seu[marker,]@assays$RNA@scale.data
  # hc <- seu_coi[markers,]@assays$RNA@scale.data
  hc <- hclust(dist(hc))
  return(hc$labels[hc$order])
}

# Compare comparable tables to each other
compareCountsToOriginal <- function(orig, comp){
  data.frame("groups"=names(orig),
             "orig"=as.numeric(orig),
             "filtered"=as.numeric(comp),
             "delta"=as.numeric(comp - orig),
             "pct_delta"=round(as.numeric(comp - orig) / as.numeric(orig), 3))
}

# capitalizes first letter (e.g. mtor to Mtor)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

getMarkers <- function(database='panglaodb', 
                       dbpath='/cluster/projects/mcgahalab/ref/scrna/panglaodb/PanglaoDB_markers_27_Mar_2020.tsv',
                       n=2, species='homo_sapiens'){
  db <- read.table(dbpath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if(database=='panglaodb'){
    db_spl <- split(db, db$organ)
    db_spl <- split(db_spl[['Immune system']], db_spl[['Immune system']]$cell.type)
    db_markers <- lapply(db_spl, function(db){
      if(species=='homo_sapiens'){
        ord <- order(with(db, (1-ubiquitousness.index) + sensitivity_human + specificity_human), 
                     decreasing = TRUE)
        head(db[ord,]$official.gene.symbol,n)
      } else if(species=='mus_musculus'){
        ord <- order(with(db, (1-ubiquitousness.index) + sensitivity_mouse + specificity_mouse), 
                     decreasing = TRUE)
        firstup(tolower(head(db[ord,]$official.gene.symbol,n)))
      } else {
        stop("species must be either 'homo_sapiens' or 'mus_musculus'")
      }
      
    })
    names(db_markers) <- gsub(" ", "_", names(db_markers))
  }
  unlist(db_markers)
}


##############################
#### 0. Read in the GEO data ####
setwd(mtxdir)

seus <- lapply(samples[1:2], function(sample){
  expression_matrix <- ReadMtx(
    mtx = paste0(sample, "_matrix.mtx.gz"),
    features = paste0(sample, "_genes.tsv.gz"),
    cells = paste0(sample, "_barcodes.tsv.gz")
  )
  seu <- CreateSeuratObject(counts = expression_matrix,
                            project= sample)
  SeuratDisk::SaveH5Seurat(seu, filename = paste0(sample, ".h5seurat"),
                           overwrite = TRUE)
  return(seu)
})
seus[[3]] <- merge(seus[[1]], y = seus[[2]], 
                   add.cell.ids = c("day30", "day8"), 
                   project = names(samples)[3])
SeuratDisk::SaveH5Seurat(seus[[3]], filename = paste0(names(samples)[3], ".h5seurat"),
                         overwrite = TRUE)
names(seus) <- names(samples)

#########################
#### 1. QC Filtering ####
dir.create(file.path(outdir, 'qc'), showWarnings = F, recursive = T)

seu_filts <- lapply(samples, function(s){
  seu <- SeuratDisk::LoadH5Seurat(file.path(mtxdir, paste0(s, ".h5seurat")))
  seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")
  
  pdf(file.path(outdir, "qc", paste0(s, ".featureScatter.pdf")), width = 14, height = 7)
  plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", 
                          feature2 = "percent.mt", shuffle = TRUE)
  plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", 
                          feature2 = "nFeature_RNA", shuffle = TRUE)
  CombinePlots(plots = list(plot1, plot2))
  dev.off()
  
  png(file.path(outdir, "qc", paste0(s, ".featureViolin.png")), width=1500, height=400)
  VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  dev.off()
  
  seu_filt <- subset(seu, subset = nFeature_RNA > 200 & 
                     nFeature_RNA < quantile(seu$nFeature_RNA, 0.99) & 
                     percent.mt < quantile(seu$percent.mt, 0.99)+0.01)
  # 99th nFeature_RNA = 2762.01
  # 99th percent.mt = 0
  
  write.table(x=compareCountsToOriginal(table(seu$orig.ident), table(seu_filt$orig.ident)), 
              file=file.path(outdir, "qc", paste0(s, ".sample_cnts.tsv")), sep="\t", 
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  return(seu_filt)
})

#################################
#### 2. Sample preprocessing ####
# Perform standard analysis of each modality independently RNA analysis
seu_harmonies <- lapply(seu_filts, function(seu_filt){
  dir.create(file.path(outdir, "harmony"), showWarnings = F, recursive = T)
  seu_filt <- NormalizeData(seu_filt, normalization.method = "LogNormalize", scale.factor = 10000)
  seu_filt <- FindVariableFeatures(seu_filt, mean.function = ExpMean, 
                                   dispersion.function = LogVMR, nfeatures = 2000)
  seu_filt <- CellCycleScoring(seu_filt, s.features = firstup(tolower(cc.genes$s.genes)), 
                               g2m.features = firstup(tolower(cc.genes$g2m.genes)), set.ident = TRUE)
  seu_filt$CC.Difference <- seu_filt$S.Score - seu_filt$G2M.Score
  seu_filt <- ScaleData(object = seu_filt, vars.to.regress = c("percent.mt", "CC.Difference"))
  seu_filt <- RunPCA(seu_filt, npcs = 30, features = VariableFeatures(object = seu_filt), verbose = FALSE)
  
  #Confirm #PC's determined explain > 90% of variance
  stdev <- seu_filt@reductions$pca@stdev
  var <- stdev^2
  PCNum <-  min(which(cumsum(var)/sum(var) >= 0.9))
  
  seu_filt <- FindNeighbors(seu_filt, reduction = "pca", dims = 1:PCNum)
  seu_filt <- FindClusters(seu_filt, resolution = 1.2)
  seu_filt <- RunUMAP(seu_filt, reduction = "pca", dims = 1:PCNum)
  
  ## Find Neighbours and Cluster with HARMONY
  pdf(file.path(outdir, "harmony", paste0(seu_filt@project.name, ".harmony_convergence2.pdf")))
  seu_harmony <- seu_filt %>% 
    RunHarmony("orig.ident", plot_convergence = TRUE)
  dev.off()
  
  seu_harmony <- FindNeighbors(object = seu_harmony, dims = 1:PCNum, reduction ="harmony")
  seu_harmony <- FindClusters(object = seu_harmony, resolution = 1.2, reduction ="harmony")
  # Tweaked the UMAP parameters here
  seu_harmony <- RunUMAP(object = seu_harmony, dims = 1:PCNum, reduction = "harmony",
                         n.neighbors = 10L, n.components = 2L, n.epochs=400L, min.dist=0.1)
  SeuratDisk::SaveH5Seurat(seu_harmony, 
                           filename = file.path(outdir, "harmony", 
                                                paste0(seu_harmony@project.name, ".harmony.h5seurat")), 
                           overwrite = TRUE)
  return(seu_harmony)
})

seu_harmonies <- lapply(seu_harmonies, function(seu_harmony){
  ids <- c('Ly108+ progenitor', 'CX3CR1-Ly108- exhausted', 'CX3CR1+ cytolytic', 'Day8 CX3CR1+')
  if(seu_harmony@project.name == "day8-30"){
    idmap <- setNames(c(rep(ids[1], 3),  rep(ids[2], 3),  rep(ids[3], 2), rep(ids[4], 4)),
                      c('1', '8', '10',  '7', '3', '9',   '0', '4',       '2', '5', '6', '11'))
  } else if (seu_harmony@project.name == 'GSM3701181_GP33_day30'){
    idmap <- setNames(c(rep(ids[1], 2), rep(ids[2], 3), rep(ids[3], 4)),
                      c('2', '7', '1', '6', '8', '0', '3', '5', '4'))
  } else if (seu_harmony@project.name == 'GSM3701180_GP33_day8'){
    idmap <- setNames(c(rep(ids[1], 2), rep(ids[4], 7)),
                      c('2', '8',       '0','1', '3', '4', '5', '6', '7'))
  }
  seu_harmony$manual_clusters <- idmap[as.character(seu_harmony$seurat_clusters)]
  
  if (seu_harmony@project.name == 'GSM3701181_GP33_day30' & FALSE){
    # this step needs to run manually
    cell_specific <- FALSE
    if(cell_specific){
      seu_harmony$manual_clusters <- 'TCF1+GrB+'
      seu_harmony$manual_clusters[which(x & !y)] <- 'TCF1+GrB-'
      seu_harmony$manual_clusters[which(!x & y)] <- 'TCF1-GrB+'
      seu_harmony$manual_clusters[which(!x & !y)] <- 'TCF1-GrB-'
    } else {
      ids2 <- c("TCF1+GrB-", "TCF1-GrB+")
      idmap <- setNames(c(rep(ids2[1], 2), rep(ids2[2], 7)),
                        c('2', '7',       '1', '6', '8', '0', '3', '5', '4'))
      seu_harmony$manual_clusters <- idmap[as.character(seu_harmony$seurat_clusters)]
    }
    seu_harmonies[[1]] <- seu_harmony
  }
  SeuratDisk::SaveH5Seurat(seu_harmony, 
                           filename = file.path(outdir, "harmony", 
                                                paste0(seu_harmony@project.name, ".harmony.h5seurat")), 
                           overwrite = TRUE)
  return(seu_harmony)
})

####################################
#### 3. Marker Detection & Load ####
clusters <- 'manual'  # 'manual' or 'seurat'
overwrite <- FALSE

seu_markers <- lapply(samples, function(s){
  if(!file.exists(file.path(outdir, "harmony", paste0(s, ".", clusters, ".markers.rda"))) | overwrite){
    print(paste0("Finding markers for: ", s))
    seu_harmony <- SeuratDisk::LoadH5Seurat(file.path(outdir, "harmony", 
                                                      paste0(s, ".harmony.h5seurat")), 
    )
    Idents(seu_harmony) <- if(clusters=='seurat') {
      seu_harmony$seurat_clusters
    } else if(clusters == 'manual'){
      seu_harmony$manual_clusters
    } else {
      stop("clusters must be 'manual' or 'seurat'")
    }
    seu_harmony_markers <- FindAllMarkers(object = seu_harmony, only.pos = TRUE, min.pct = 0.25, 
                                          thresh.use = 0.25)
    # seu.harmony.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
    col_idx <- which(colnames(seu_harmony_markers) == 'avg_log2FC')
    colnames(seu_harmony_markers)[col_idx] <- 'avg_logFC'
    write.table(seu_harmony_markers, 
                file=file.path(outdir, "harmony", paste0(s, ".", clusters, ".seurat_markers.csv")),
                sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
    save(seu_harmony_markers, file=file.path(outdir, "harmony", paste0(s, ".", clusters, ".markers.rda")))
  } else {
    load(file.path(outdir, "harmony", paste0(s, ".", clusters, ".markers.rda")))
  }
  return(seu_harmony_markers)
})
seu_harmonies <- lapply(samples, function(s){
  SeuratDisk::LoadH5Seurat(file.path(outdir, "harmony", paste0(s, ".harmony.h5seurat")))
})

######################################################
#### 4. Feature/Dim plots for markers of interest ####
dir.create(file.path(outdir, "featureplot"), showWarnings = FALSE)

db_markers <- getMarkers(n=1, species = 'mus_musculus')
clusters <- 'manual'  # 'manual' or 'seurat'

# Feature plots of the sorted Marker Genes
lapply(seu_harmonies, function(seu_harmony){
  dp1 <- DimPlot(object = seu_harmony, reduction = "umap", 
                 group.by=if(clusters=='seurat') 'seurat_clusters' else 'manual_clusters', 
                label = TRUE, pt.size = 1)
  if(seu_harmony@project.name == names(samples[3])){
    dp2 <- DimPlot(object = seu_harmony, reduction = "umap", group.by='orig.ident', 
                   label = TRUE, pt.size = 1)
  }
  
  fp1 <- FeaturePlot(seu_harmony, features = db_markers, pt.size = 1, ncol = 4,
              cols=c("#f7f4f9", "#91003f"), order=TRUE, keep.scale=NULL, 
              label=FALSE, raster=TRUE)
  fp2 <- FeaturePlot(seu_harmony, features = markers, pt.size = 1, ncol = 4,
              cols=c("#f7f4f9", "#91003f"), order=TRUE, keep.scale=NULL, 
              label=FALSE, raster=TRUE)
  fp3 <- FeaturePlot(seu_harmony, features = isg_markers, pt.size = 1, ncol = 6,
                     cols=c("#f7f4f9", "#91003f"), order=TRUE, keep.scale=NULL, 
                     label=FALSE, raster=TRUE)
  
  pdf(file.path(outdir, "featureplot", paste0(seu_harmony@project.name, ".featplots.pdf")), 
      height=23, width=23)
  if(seu_harmony@project.name == names(samples[3])){ 
    print(plot_grid(dp2, fp2, nrow=2, rel_widths = c(1, 2)))
  }
  print(plot_grid(dp1, fp2, nrow=2, rel_widths = c(1, 2)))
  # print(fp1) #dbamarkers
  print(fp3)
  dev.off()
})

#######################################
#### 5. Heatmap of cluster markers ####
clusters <- 'manual'  # 'manual' or 'seurat'
n <- 20

lapply(names(seu_markers), function(s){
  smarker <- seu_markers[[s]]
  seu_harmony <- seu_harmonies[[s]]
  
  smarker %>% 
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    top_n(n=n, wt=avg_logFC) -> topn
  
  Idents(seu_harmony) <- if(clusters == 'manual') {
    seu_harmony$manual_clusters
  } else if(clusters == 'seurat'){
    seu_harmony$seurat_clusters
  }
  
  getHClustLabels <- function(seu, marker){
    hc <- seu[marker,]@assays$RNA@scale.data
    # hc <- seu_harmony[markers,]@assays$RNA@scale.data
    hc <- hclust(dist(hc))
    return(hc$labels[hc$order])
  }
  
  hp1 <- DoHeatmap(seu_harmony, features=topn$gene, raster=F, angle = 90,
                   size=3, slot='scale.data')
  hp2 <- DoHeatmap(seu_harmony, features=getHClustLabels(seu_harmony, isg_markers), 
                   raster=F, angle = 90, size=3, slot='scale.data')
  hp3 <- DoHeatmap(seu_harmony, features=getHClustLabels(seu_harmony, markers), 
                   raster=F, angle = 90, size=3, slot='scale.data')
  pdf(file.path(outdir, "heatmap", paste0(seu_harmony@project.name, ".heatmaps.pdf")), 
      height=12)
  print(plot_grid(hp1, hp3, nrow=2))
  print(hp2)
  dev.off()
})

######################################################
#### 6.a Subset for Tcf-1/GranzymeB and recluster ####
seu_harmony <- seu_harmonies[['day30']]
ids2 <- c("TCF1+GrB-", "TCF1-GrB+")
idmap <- setNames(c(rep(ids2[1], 2), rep(ids2[2], 7)),
                  c('2', '7',       '1', '6', '8', '0', '3', '5', '4'))
seu_harmony$manual_clusters <- idmap[as.character(seu_harmony$seurat_clusters)]
Idents(seu_harmony) <- seu_harmony$manual_clusters
seu_coi <- subset(x = seu_harmony, idents = "TCF1-GrB+")

#########################################
##  6.b renormalize just the GnB cells ##
seu_coi <- RunPCA(seu_coi, npcs = 30, features = VariableFeatures(object = seu_coi), verbose = FALSE)

#Confirm #PC's determined explain > 90% of variance
stdev <- seu_coi@reductions$pca@stdev
var <- stdev^2
PCNum <-  min(which(cumsum(var)/sum(var) >= 0.9))

seu_coi <- FindNeighbors(object = seu_coi, dims = 1:PCNum, reduction ="harmony")
seu_coi <- FindClusters(object = seu_coi, resolution = 1.2, reduction ="harmony")
# Tweaked the UMAP parameters here
seu_coi <- RunUMAP(object = seu_coi, dims = 1:PCNum, reduction = "harmony",
                   n.neighbors = 10L, n.components = 2L, n.epochs=400L, min.dist=0.1)

###################################
##  6.c Manually assign clusters ##
ids2 <- c("CX3CR1_hi", "CX3CR1_mid", "CX3CR1_null")
idmap <- setNames(c(rep(ids2[1], 2), rep(ids2[2], 3), rep(ids2[3], 4)),
                  c('0', '2',   '4', '5', '6',   '7', '1', '3', '8'))
seu_coi$manual_clusters <- idmap[as.character(seu_coi$seurat_clusters)]


####################################
## 6.d DEG of the subsetted cells ##
Idents(seu_coi) <- seu_coi$manual_clusters
seu_harmony_markers <- FindAllMarkers(object = seu_coi, only.pos = TRUE, min.pct = 0.25, 
                                      thresh.use = 0.25)
# seu.harmony.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
col_idx <- which(colnames(seu_harmony_markers) == 'avg_log2FC')
colnames(seu_harmony_markers)[col_idx] <- 'avg_logFC'
write.table(seu_harmony_markers, 
            file=file.path(outdir, "harmony", paste0("GnB-subset.seurat_markers.csv")),
            sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
save(seu_harmony_markers, file=file.path(outdir, "harmony", paste0("GnB-subset.seurat_markers.rda")))

############################################################
## 6.e Visualize the GnB specific clusters and projection ##
clusters <- 'seurat'
dp1 <- DimPlot(object = seu_coi, reduction = "umap", 
               group.by='seurat_clusters', 
               label = TRUE, pt.size = 1)
dp2 <- DimPlot(object = seu_coi, reduction = "umap", 
               group.by='manual_clusters', 
               label = TRUE, pt.size = 1)

fp1 <- FeaturePlot(seu_coi, features = markers, pt.size = 1, ncol = 4,
                   cols=c("#f7f4f9", "#91003f"), order=TRUE, keep.scale=NULL, 
                   label=FALSE, raster=TRUE)
fp2 <- FeaturePlot(seu_coi, features = isg_markers, pt.size = 1, ncol = 6,
                   cols=c("#f7f4f9", "#91003f"), order=TRUE, keep.scale=NULL, 
                   label=FALSE, raster=TRUE)

pdf(file.path(outdir, "featureplot", paste0("GnB-subset.featplots.pdf")), 
    height=23, width=23)
plot_grid(dp2, fp1, nrow=2, rel_widths = c(1, 2))
plot_grid(dp1, fp1, nrow=2, rel_widths = c(1, 2))
fp2
dev.off()

#################################
## 6.f Heatmap of DEG and ISGs ##
Idents(seu_coi) <- seu_coi$manual_clusters
n <- 20

seu_harmony_markers %>% 
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(n=n, wt=avg_logFC) -> topn

hp1 <- DoHeatmap(seu_coi, features=topn$gene, raster=F, angle = 90,
                 size=3, slot='scale.data')
hp2 <- DoHeatmap(seu_coi, features=getHClustLabels(seu_coi, isg_markers), 
                 raster=F, angle = 90, size=3, slot='scale.data')
hp3 <- DoHeatmap(seu_coi, features=getHClustLabels(seu_coi, markers), 
                 raster=F, angle = 90, size=3, slot='scale.data')
pdf(file.path(outdir, "heatmap", paste0("GnB-subset.heatmaps.pdf")), 
    height=12)
print(plot_grid(hp1, hp3, nrow=2))
print(hp2)
dev.off()
