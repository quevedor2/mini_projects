library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(harmony)
library(gridExtra)
library(vegan)
library(reshape2)
library(infercnv)
library(RColorBrewer)
library(scales)

##############################################
#### 1. Set up environment and parameters ####
PDIR <- "/cluster/projects/mcgahalab/data/mcgahalab/pdac_ahr"
setwd(file.path(PDIR, "datasets"))
outdir <- file.path(PDIR, "results")
outdir_plus <- file.path(outdir, "features2")
annotation_method <- 'scsa' # natcan or scsa

datasets <- list('PRJCA001063'=c('PRJCA001063_CRC_besca2.annotated', 
                                 'steps_nonorm',
                                 'AHR'),
                 'GSE155698'=c('GSE155698',
                               'steps_norm',
                               'AHR'))

# Genes of interest (Myeloid)
goi <- c("CD14", "ITGAM", "CD68", "MARCO", "FCGR3B", "CD163", "SIGLEC1", "SIRPA", "FCGR3A")
# goi <- c("CD14", "CD163", "MARCO", "SIGLEC1", "FCGR3A", "FCGR3B", "TNFSF15", "TREM1", "IL10", "IRF8", "RAET1E")
gene <- 'AHR'
markers <- c(gene, goi)

# Cell type markers
cell_markers <- list("NatureCancer_auto"=list(
    "Epithelial"=c("KRT7","KRT8","KRT18","KRT19","EPCAM","CDH1"),
    "T-cell"=c("CD3E","CD3G","CD3D","CD4","IL7R","CD8A","LEF1"),
    "Myeloid"=c("CD14","ITGAM","MNDA","MPEG1","ITGAX"),
    "B-cell"=c("CD79A","MS4A1","CD19"),
    "Fibroblasts"=c("CDH11","PDGFRA","PDGFRB","ACTA2"),
    "RBC"=c("HBA1","HBB","HBA2"),
    "NK"=c("NCR3","FCGR3A","NCAM1","KLRF1","KLRC1","CD38","KLRC1"),
    "Endothelial"=c("CDH5","PECAM1"),
    "Acinar"=c("TRY4","SPINK1","AMY2A")
  ),
  "NatureCancer_manual"=list(
    Bcells = c('CD79A', 'MS4A1','SDC1','IGJ','IGLL5','CXCR4','KIT','CD27','HLA-DRA'),
    DCS = c('ITGAE','LYZ','CLEC9A','BATF3','IRF8','IDO1','CD207','CD1A','CD1C', 'HLA-DRA','CCL22','LAMP3','IL22RA2','CD101'),
    Endothelial = c('VWF','CDH5'),
    Epithelial = c('PRSS1', 'CTRB2','REG1A','CLU','MKI67','KRT8','SPINK1','KRT19','KRT18', 'KRT18','TFF1','MUC1'),
    Fibroblasts = c('ACTA2', 'CDH11', 'PDGFRB', 'COL1A1', 'COL3A1', 'RGS5', 'IGFBP7', 'PDPN','DCN','MCAM','IL6','APOE','GLI1','GLI2','GLI3','PDGFA'),
    Mast_cells = c('TPSAB1','CPA3'),
    Myeloid = c('CD14','ITGAM','FCGR3A','FCGR3B','APOE','C1QA','MARCO','LYZ','HLA-DRA'),
    Tcells_Nk = c('CD2', 'CD3D','CD3E','NCAM1','NKG7','CD4','CD8A','PRF1','IFNG', 'GZMB','CD69','FOXP3','TIGIT','TOP2A','FCGR3A')
  ),
  "CellResearch"=list(
    acinar = c('PRSS1', 'CTRB1', 'CTRB2', 'REG1B'),
    ductal_cell_1 = c('AMBP', 'CFTR', 'MMP7'),
    ductal_cell_2 = c('KRT19', 'KRT7', 'TSPAN8', 'SLPI'),
    endocrine_cell = c('CHGB', 'CHGA', 'INS', 'IAPP'),
    endothelial_cell = c('CDH5', 'PLVAP', 'VWF', 'CLDN5'),
    fibroblast = c('LUM', 'DCN', 'COL1A1'),
    stellate_cell = c('RGS5', 'ACTA2', 'PDGFRB', 'ADIRF'),
    macrophage = c('AIF1', 'CD64', 'CD14', 'CD68'),
    b_cell = c('MS4A1', 'CD79A', 'CD79B', 'CD52'),
    t_cell = c('CD3D', 'CD3E', 'CD4', 'CD8')
  )
)

# Compare comparable tables to each other
compareCountsToOriginal <- function(orig, comp){
  data.frame("groups"=names(orig),
             "orig"=as.numeric(orig),
             "filtered"=as.numeric(comp),
             "delta"=as.numeric(comp - orig),
             "pct_delta"=round(as.numeric(comp - orig) / as.numeric(orig), 3))
}

#############
#### 2. GSE ####
ds_idx <- 2
dataset   <- names(datasets)[ds_idx]
file      <- datasets[[ds_idx]][1]
analysis  <- datasets[[ds_idx]][2]
gene      <- datasets[[ds_idx]][3]

###################################
#### 2a. Create Seurat object  ####
PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/pdac_ahr/datasets/GSE155698/raw_data'
samples <- list.files(PDIR, pattern="TISSUE|PBMC")
setwd(PDIR)

seus <- lapply(samples, function(s){
  mtx <- Read10X(data.dir = file.path(s, "filtered_feature_bc_matrix"), strip.suffix=TRUE)
  seu <- CreateSeuratObject(counts = mtx, project = s)
})

seu <- merge(seus[[1]], y = seus[-1], 
             add.cell.ids = samples, 
             project = dataset)
SeuratDisk::SaveH5Seurat(seu, filename = file.path(dataset, paste0(file, ".h5seurat")))


#####################################
#### 2b. Load raw Seurat object  ####
## Load data and use reference PBMC to map PBMC labels
seu <- SeuratDisk::LoadH5Seurat(file.path(dataset, paste0(file, ".h5seurat")))
seu$sample.type <- gsub("_[0-9a-zA-Z]*$", "", seu$orig.ident)
seu$group.ident <- gsub("_[a-zA-Z0-9]+$", "", seu$orig.ident)
seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")

group_cells_orig <- table(seu$group.ident)
sample_cells_orig <- table(seu$orig.ident)

###################################
#### 3a. QC of the cell counts ####
dir.create(file.path(outdir, "qc"), showWarnings = FALSE)
pdf(file.path(outdir, "qc", "qc_metrics.pdf"), width = 16)
Idents(seu) <- seu$group.ident
mt_by_cnt   <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt", shuffle=TRUE)
feat_by_cnt <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", shuffle=TRUE)
mt_by_cnt + feat_by_cnt
dev.off()

# Remove outliers
seu <- subset(seu, subset = nFeature_RNA > quantile(seu$nCount_RNA, 0.05) & 
                 nFeature_RNA < quantile(seu$nFeature_RNA, 0.95) & 
                 percent.mt < 10)
# 95th nCount_RNA = 4836
# 5th nFeature_RNA = 581
# 10 percent.mt = 85th percentile

write.table(x=compareCountsToOriginal(group_cells_orig, table(seu$group.ident)), 
            file=file.path(outdir, "qc", "group_cnts.tsv"), sep="\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE)
write.table(x=compareCountsToOriginal(sample_cells_orig, table(seu$orig.ident)), 
            file=file.path(outdir, "qc", "sample_cnts.tsv"), sep="\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE)

#######################################
#### 3b. Published Method (Harmony) ####
seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- FindVariableFeatures(object = seu, mean.function = ExpMean, 
                            dispersion.function = LogVMR, nfeatures = 2000)
seu <- CellCycleScoring(seu, s.features = cc.genes$s.genes, 
                        g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
seu$CC.Difference <- seu$S.Score - seu$G2M.Score
seu <- ScaleData(object = seu, vars.to.regress = c("percent.mt", "CC.Difference"))
save(seu, file=file.path(dataset, paste0(file, "_scale-mt-cc.rda")))

seu <- RunPCA(seu, npcs = 30, features = VariableFeatures(object = seu), verbose = FALSE)
stdev <- seu@reductions$pca@stdev
var <- stdev^2

EndVar = 0

for(i in 1:length(var)){
  total <- sum(var)
  numerator <- sum(var[1:i])
  expvar <- numerator/total
  if(EndVar == 0){
    if(expvar > 0.9){
      EndVar <- EndVar + 1
      PCNum <- i
    }
  }
}
#Confirm #PC's determined explain > 90% of variance
sum(var[1:PCNum])/ sum(var)

seu <- FindNeighbors(seu, reduction = "pca", dims = 1:PCNum)
seu <- FindClusters(seu, resolution = 1.2)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:PCNum)

## Find Neighbours and Cluster with HARMONY
pdf(file.path(outdir, "harmony_convergence2.pdf"))
seu <- seu %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
dev.off()

seu.harmony <- FindNeighbors(object = seu, dims = 1:PCNum, reduction ="harmony")
seu.harmony <- FindClusters(object = seu, resolution = 1.2, reduction ="harmony")
# Tweaked the UMAP parameters here
seu.harmony <- RunUMAP(object = seu.harmony, dims = 1:PCNum, reduction = "harmony",
                       n.neighbors = 10L, n.components = 2L, n.epochs=400L, min.dist=0.1)
SaveH5Seurat(seu.harmony, filename = file.path(dataset, paste0(file, "_harmony2.h5seurat")), overwrite = TRUE)

#################################
#### 4. Load Annotated Clusters ####
harmony_file <- 'harmony2' # harmony or harmony2
seu.harmony <- SeuratDisk::LoadH5Seurat(file.path(dataset, paste0(file, "_", harmony_file, ".h5seurat")))

if(!file.exists(file.path(dataset, paste0(file, "_", harmony_file, "-markers.rda")))){
  Idents(seu.harmony) <- seu.harmony$seurat_clusters
  seu.harmony.markers <- FindAllMarkers(object = seu.harmony, only.pos = TRUE, min.pct = 0.25, 
                                        thresh.use = 0.25)
  # seu.harmony.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
  write.table(seu.harmony.markers, file=file.path(outdir_plus, "seurat_markers.csv"),
              sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)
  save(seu.harmony.markers, AutomatedClusterMarkerTable, assignClusterId,
       file=file.path(dataset, paste0(file, "_", harmony_file, "-markers.rda")))
  
} else {
  #load(file=file.path(dataset, paste0(file, "_", harmony_file, "-markers.rda")))
  # x <- lapply(seq_along(seu.harmony.markers), function(idx){
  #   seu.harmony.markers[[idx]]$cluster <- as.character(idx - 1)
  #   seu.harmony.markers[[idx]]$gene <- rownames(seu.harmony.markers[[idx]])
  #   seu.harmony.markers[[idx]]
  # })
  # seu.harmony.markers <- do.call(rbind, x)
  # colnames(seu.harmony.markers)[2] <- 'avg_logFC'
  
  print("Reading harmony markers")
  load(file=file.path(dataset, paste0(file, "_", harmony_file, "-markers.rda")))
}

# Annotate the cell types for each cluster
# This function is found in AutomatedClusterMarkerTable.R
Idents(seu.harmony) <- seu.harmony$seurat_clusters
if(annotation_method=='natcan'){
  anno <- lapply(setNames(names(cell_markers), names(cell_markers)), function(cm){
    seu_anno <- AutomatedClusterMarkerTable(Seurat_Object=seu.harmony, 
                                            ClusterList=split(seu.harmony.markers, seu.harmony.markers$cluster),
                                            cell_markers=cell_markers[[cm]],
                                            exact.match=TRUE, cyc=FALSE, rp=FALSE)
    seu_anno[[2]]
  })
  anno <- do.call(cbind, anno)
  anno[grep("^[0-9]+$", anno)] <- NA
  seu_anno <- apply(anno[,c('NatureCancer_manual', 'NatureCancer_auto', 'CellResearch')], 1, function(i) unique(na.omit(i))[1])
  seu_anno[which(is.na(seu_anno))] <- '1234'
  
  cluster_ids <- seu_anno # seu_anno[[2]]
  unk_idx <- grep("^.*[0-9]+$", cluster_ids)
  cluster_ids[unk_idx] <- 'Unknown'
  cluster_ids <- setNames(cluster_ids, seq_along(cluster_ids)-1)
} else if(annotation_method=='scsa' & file.exists(file.path(outdir_plus, 'annotate.txt'))){
  cluster_keep <- length(table(seu.harmony.markers$cluster))
  nL <- R.utils::countLines(file.path(outdir_plus, 'annotate.txt'))
  df <- read.csv(file.path(outdir_plus, 'annotate.txt'), header=FALSE, skip=nL-cluster_keep,
                 check.names = FALSE, stringsAsFactors = FALSE)
  
  df$V1 <-  gsub("\\[", "", df$V1)
  for (coli in seq_along(df)) {
    df[,coli] <- gsub("\\'", "", df[,coli])
  }
  df$sortedID <- sapply(strsplit(df$V3, split="\\|"), function(i) paste(sort(gsub("^ ", "", i)), collapse="/"))
  cluster_ids <- setNames(df$sortedID, df$V1)
}


# Add the annotated cell types to Seurat metadata
seu.harmony$anno_clusters <- cluster_ids[as.character(seu.harmony$seurat_clusters)]

######################################
#### 5.a Feature ordering and prep ####
dir.create(file.path(outdir_plus), recursive = TRUE, showWarnings = FALSE)

# Report the annotation of Harmony clusters
clusters <- data.frame('harmony'=seu.harmony$seurat_clusters,
                       'anno'=seu.harmony$anno_clusters)
clusters_harmony  <- unique(clusters[order(clusters[,1]),])
clusters_anno     <- unique(clusters[order(clusters[,2]),])

# Compute the spearman correlation between mean-expression per cluster and percent-expressed 
# per cluster of the AHR gene compared to a set of marker genes
# Once computer, calculate the geom_mean of the two values and sort
Idents(seu.harmony) <- seu.harmony$seurat_clusters
p <- DotPlot(seu.harmony, features = markers, # unique(unlist(natcan_markers)),
             cols = c('grey', 'red'), dot.scale = 10)

expmat <- dcast(p$data, features.plot ~ id, value.var='avg.exp')
pctmat <- dcast(p$data, features.plot ~ id, value.var='pct.exp')
expcor <- cor(t(expmat[,-1]), method='spearman')
pctcor <- cor(t(pctmat[,-1]), method='spearman')
colnames(expcor) <- rownames(expcor) <- colnames(pctcor) <- rownames(pctcor) <- expmat$features.plot
geomcor <- (expcor + pctcor) / 2 
sorted_markers <- names(sort(geomcor[,'AHR'], decreasing = TRUE))

# Reorder harmony clusters so anno-clusters group
ord <- factor(as.character(Idents(seu.harmony)), levels=as.character(clusters_anno$harmony))
seu.harmony$seurat_clusters_ord <- setNames(ord, names(Idents(seu.harmony)))

######################
#### 5.b Dotplots ####
# Dotplots showing the mean expression + percent-expressed-in-cluster for marker genes
# Harmony Clusters
pdf(file.path(outdir_plus, "nc-dotplot_clusters.pdf"), width=15, height = 12)
Idents(seu.harmony) <- seu.harmony$seurat_clusters_ord
DotPlot(seu.harmony, features = sorted_markers, # unique(unlist(natcan_markers)),
        cols = c('grey', 'darkred'), dot.scale = 10) + 
  RotatedAxis() + FontSize(x.text = 17, y.text = 17)
dev.off()

# Dotplots showing the mean expression + percent-expressed-in-cluster for marker genes
# Annotated clusters
pdf(file.path(outdir_plus, "nc-dotplot_anno.pdf"), width=15, height = 12)
Idents(seu.harmony) <- seu.harmony$anno_clusters
DotPlot(seu.harmony, features = sorted_markers, # unique(unlist(natcan_markers)),
        cols = c('grey', 'darkred'), dot.scale = 10) + 
  RotatedAxis() + FontSize(x.text = 17, y.text = 17)
dev.off()

# Annotated clusters, separated by tissue type
ps <- lapply(unique(seu.harmony$group.ident), function(ident){
  Idents(seu.harmony) <- seu.harmony$group.ident
  sub.seu <- subset(x=seu.harmony, idents=ident)
  Idents(sub.seu) <- sub.seu$anno_clusters
  
  DotPlot(sub.seu, features = sorted_markers, # unique(unlist(natcan_markers)),
          cols = c('grey', 'darkred'), dot.scale = 10) + 
    RotatedAxis() + FontSize(x.text = 17, y.text = 17)
})
names(ps) <- unique(seu.harmony$group.ident)
pdf(file.path(outdir_plus, "nc-dotplot_anno-groups.pdf"), width=15, height=12)
ps[['AdjNorm_TISSUE']]
ps[['Healthy_PBMC']]
ps[['PDAC_PBMC']]
ps[['PDAC_TISSUE']]
dev.off()
Idents(seu.harmony) <- seu.harmony$anno_clusters

group_exps <- lapply(unique(seu.harmony$group.ident), function(ident){
  print(ident)
  Idents(seu.harmony) <- seu.harmony$group.ident
  sub.seu <- subset(x=seu.harmony, idents=ident)
  
  uid_dots <- lapply(unique(sub.seu$orig.ident), function(uid){
    Idents(sub.seu) <- sub.seu$orig.ident
    sub.uid.seu <- subset(x=sub.seu, idents=uid)
    
    Idents(sub.uid.seu) <- sub.uid.seu$anno_clusters
    dat <- DotPlot(sub.uid.seu, features = sorted_markers)$data[,c(1:4)]
    dat$UID <- paste0(dat[,3], "_", dat[,4])
    return(dat)
  })
  
  avg_exp <- Reduce(function(x,y) merge(x,y, by=c('UID')), 
                    lapply(uid_dots, function(i) i[,c('avg.exp', 'UID')]))
  pct_exp <- Reduce(function(x,y) merge(x,y, by=c('UID')), 
                    lapply(uid_dots, function(i) i[,c('pct.exp', 'UID')]))
  
  return(list("exp"=avg_exp,"pct"=pct_exp))
})
names(group_exps) <- unique(seu.harmony$group.ident)

common_uid <- unique(sort(unlist(sapply(group_exps, function(i) i[[1]]$UID))))
uid_df <- data.frame("UID"=common_uid)
t_ret <- 'stat'

# Go by avg_expression or percent_expressed
grp_delta <- lapply(setNames(c("exp", "pct"), c("exp", "pct")), function(grp){
  # Set the comparison group 
  lapply(group_exps, function(groupa){
    ga <- groupa[[grp]]
    ga <- merge(uid_df, ga, by='UID', all.x=TRUE)
    sapply(group_exps, function(groupb){
      gb <- groupb[[grp]]
      gb <- merge(uid_df, gb, by='UID', all.x=TRUE)
      delta <- sapply(seq_along(rownames(ga)), function(rowidx){
        tryCatch({
          tres <- t.test(ga[rowidx,-1], gb[rowidx,-1])
          if(t_ret=='stat') tres$statistic else tres$p.value
        }, error=function(e){NA})
      })
      setNames(delta, ga$UID)
    })
  })
})



Idents(seu.harmony) <- seu.harmony$anno_clusters


#############################
#### 5.c Dimension Plots ####
# Umap visualization, separated by clusters
pdf(file.path(outdir_plus, "nc-dimplot_clusters.pdf"), width=17)
p1 <- DimPlot(object = seu.harmony, reduction = "umap", group.by='seurat_clusters', 
              label = TRUE, pt.size = 0.5)
p2 <- DimPlot(object = seu.harmony, reduction = "umap", group.by='anno_clusters', 
              label = TRUE, pt.size = 0.5)
p1 + p2
p2
dev.off()

# Umap visualization of individual tissue types
ps <- lapply(unique(seu.harmony$group.ident), function(ident){
  Idents(seu.harmony) <- seu.harmony$group.ident
  sub.seu <- subset(x=seu.harmony, idents=ident)
  
  DimPlot(object = sub.seu, reduction = "umap", group.by='anno_clusters', 
          label = TRUE, pt.size = 0.5,
          cols=alpha(colorRampPalette(brewer.pal(9,"Set1"))(length(unique(sub.seu$anno_clusters))),
                    0.1))
  
})
names(ps) <- unique(seu.harmony$group.ident)
pdf(file.path(outdir_plus, "nc-dimplot_groups.pdf"), width=17)
ps[['AdjNorm_TISSUE']]
ps[['Healthy_PBMC']]
ps[['PDAC_PBMC']]
ps[['PDAC_TISSUE']]
dev.off()
Idents(seu.harmony) <- seu.harmony$anno_clusters


# Feature plots of the sorted Marker Genes
pdf(file.path(outdir_plus, "nc-featplots.pdf"), height=14, width=14)
FeaturePlot(seu.harmony, features = sorted_markers, pt.size = 0.5, ncol = 4,
            cols=c("#f7f4f9", "#91003f"), order=TRUE, keep.scale=NULL, label=FALSE)
dev.off()


#######################################
#### 5.d Finding top similar genes ####
### Summarize mean-exp and pct in each cell type/cluster
seu_cl <- lapply(levels(seu.harmony$seurat_clusters), function(cl){
  cl_idx <- which(seu.harmony$seurat_clusters == cl)
  gene_means <- rowMeans(GetAssayData(seu.harmony)[,cl_idx])
  gene_pcts <- rowSums(GetAssayData(seu.harmony)[,cl_idx] > 0) / length(cl_idx)
  
  # Returns the mean gene-expr and percent-expressed in cluster 'cl'
  data.frame("mean"=gene_means, "pct"=gene_pcts)
})

# combine all clusters together and separated into 'expr' and 'percent' dataframes
cl_means  <- as.data.frame(do.call(cbind, lapply(seu_cl, function(i) i[[1]])))
cl_pcts   <- as.data.frame(do.call(cbind, lapply(seu_cl, function(i) i[[2]])))
rownames(cl_means) <- rownames(cl_pcts) <- rownames(seu_cl[[1]])
colnames(cl_means) <- colnames(cl_pcts) <- levels(seu.harmony$seurat_clusters)

# Calculate euclidean distance between AHR and every gene listed for cluster expr/percent-expressed
expr_sim  <- apply(cl_means, 1, function(i) dist(rbind(as.numeric(cl_means['AHR',]), as.numeric(i))))
pct_sim   <- apply(cl_pcts, 1, function(i) dist(rbind(as.numeric(cl_pcts['AHR',]), as.numeric(i))))
comb_sim  <- (expr_sim + pct_sim)/2
ord_idx   <- order(comb_sim, decreasing = FALSE)
sim_df    <- data.frame("gene"=names(comb_sim),
                        "dist_mean_expr"=expr_sim, 
                        "dist_percent_expr"=pct_sim, 
                        "dist_combined"=comb_sim)[ord_idx,]
write.table(sim_df, file = file.path(outdir_plus, "ahr_dist.tsv"), 
            sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

# Calculate a weighted distance between genes similar to AhR
cl_weights    <- round(table(seu.harmony$seurat_clusters) / max(table(seu.harmony$seurat_clusters)),3)
weight_means  <- sweep(cl_means, 2, cl_weights, function(x,y) x * sqrt(y))
w_expr_sim    <- apply(weight_means, 1, function(i) dist(rbind(as.numeric(weight_means['AHR',]), as.numeric(i))))
weight_pcts   <- sweep(cl_pcts, 2, cl_weights, function(x,y) x * sqrt(y))
w_pct_sim     <- apply(weight_pcts, 1, function(i) dist(rbind(as.numeric(weight_pcts['AHR',]), as.numeric(i))))
w_comb_sim    <- (w_expr_sim + w_pct_sim)/2
ord_idx       <- order(w_comb_sim, decreasing = FALSE)
w_sim_df      <- data.frame("gene"=names(w_comb_sim),
                            "dist_mean_expr"=w_expr_sim, 
                            "dist_percent_expr"=w_pct_sim, 
                            "dist_combined"=w_comb_sim)[ord_idx,]
write.table(w_sim_df, file = file.path(outdir, "features", "weighted_ahr_dist.tsv"), 
            sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

# Plot the top 25 most similar features to AHR
pdf(file.path(outdir_plus, "nc-featplots_top25.pdf"), height=15, width=15)
FeaturePlot(seu.harmony, features = names(head(sort(w_comb_sim, decreasing = FALSE), 25)), 
            pt.size = 0.5, ncol = 5, cols=c("#f7f4f9", "#91003f"), order=TRUE, 
            keep.scale='all', label=FALSE)
dev.off()

pdf(file.path(outdir_plus, "nc-featplots_top25.pdf"), height=15, width=15)
FeaturePlot(seu.harmony, features = names(head(sort(w_comb_sim, decreasing = FALSE), 25)), 
            pt.size = 0.5, ncol = 5, cols=c("#f7f4f9", "#91003f"), order=TRUE, 
            keep.scale='all', label=FALSE)
dev.off()

pdf(file.path(outdir_plus, "nc-dotplot_top25.pdf"), height=15, width=20)
Idents(seu.harmony) <- seu.harmony$seurat_clusters_ord
DotPlot(seu.harmony, features =  names(head(sort(w_comb_sim, decreasing = FALSE), 25)), # unique(unlist(natcan_markers)),
        cols = c('grey', 'darkred'), dot.scale = 10) + 
  RotatedAxis() + FontSize(x.text = 17, y.text = 17)
dev.off()

##################
#### 6. inferCNV ####

## **NOTE**: Refer to tn-matched_infercnv.R for script on running the main infercnv pipeline
infercnv_outdir <- '/cluster/projects/mcgahalab/data/mcgahalab/pdac_ahr/results/infercnv/matched_pbmc-miniclust'
infercnv_obj <- "21_denoiseHMMi6.NF_NA.SD_1.5.NL_FALSE.infercnv_obj"
infercnv_cnv <- "expr.infercnv.18_HMM_pred.Bayes_Net.Pnorm_0.5.dat"
ranges <- seq(0, 1, by=0.01)
qvals <- c('5%', '95%')   # Quantiles to use as threshold of CNV calling
norm_limits <- TRUE       # Whether to use just normal cells to established denoised limits

dir.create(file.path(infercnv_outdir,  "pga"))
samples <- grep(list.files(infercnv_outdir), pattern="pga", invert=TRUE, value = TRUE)
pdf(file.path(infercnv_outdir, "pga", "pga_boxplots.pdf"), width = 9, height = 12)
pga_tissues <- lapply(samples, function(tissue){
  #tissue <- list.files(infercnv_outdir)[1]
  
  obj   <- readRDS(file.path(infercnv_outdir, tissue, infercnv_obj))
  cnv_gene <- read.table(file.path(infercnv_outdir, tissue, infercnv_cnv), sep="\t", 
                         header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  
  ## Establsh the limits for loss/normal/gain gene expr
  if(norm_limits){
    message("Using all NORMAL cells to establish CNV limits")
    normref <- (sapply(strsplit(colnames(obj@expr.data), split="_"), function(i) i[[2]]) != 'TISSUE')
    quantile_vals <- round(quantile(as.numeric(obj@expr.data[,which(normref)]), ranges),4)
  } else {
    message("Using all TISSUE and NORMAL cells to establish CNV limits")
    quantile_vals <- round(quantile(as.numeric(obj@expr.data), ranges),4)
  }
  quantile_thresh <- c(0, quantile_vals[qvals], 2)
  lvls <- levels(factor(cut(quantile_thresh, quantile_thresh)))
  
  ## Split gene expr into loss (-1), normal (0), or gain (1)
  thresh_expr <- cut(obj@expr.data, quantile_thresh) %>%
    as.integer(factor(., levels=lvls)) - median(seq_along(lvls))
  thresh_expr <- matrix(thresh_expr, ncol=ncol(obj@expr.data))
  rownames(thresh_expr) <- rownames(obj@expr.data)
  colnames(thresh_expr) <- colnames(obj@expr.data)
  
  ## Calculate percent-genome-altered (PGA)
  genes       <- obj@gene_order
  gene_width  <- with(genes, stop-start)
  
  pga <- apply(thresh_expr, 2, function(cellexpr){
    sum(gene_width[cellexpr != 0])
  })
  pga_thresh <- round(pga/sum(gene_width), 4)
  
  ## Calculate PGA from outputted Bayes GenexSample matrix
  pga <- apply(cnv_gene, 2, function(cellexpr){
    sum(gene_width[cellexpr != 3])
  })
  pga_region <- round(pga/sum(gene_width), 4)
  
  ## Calculate PGA thresholded to CNV regions only
  pga <- sapply(seq_along(colnames(cnv_gene)), function(sidx){
    sum(gene_width[which((cnv_gene[,sidx] != 3) & (thresh_expr[,sidx] != 0))])
  })
  pga_capped <- round(pga/sum(gene_width), 4)
  
  ## Report PGA by group
  par(mfrow=c(3,1))
  sapply(c('pga_capped', 'pga_thresh', 'pga_region'), function(pga_type){
    pga <- switch(pga_type,
                  pga_capped=pga_capped,
                  pga_thresh=pga_thresh,
                  pga_region=pga_region)
    
    pga_grps <- lapply(obj@tumor_subclusters$subclusters, function(idx){
      pga[idx[[1]]]
    })
    tn_grps <- split(pga,normref)
    pga_grps <- c(pga_grps, tn_grps)
    
    boxplot(pga_grps, ylim=c(0,0.6), las=2, cex.axis=0.6, 
            ylab='PGA', xlab='Clusters', main=paste0(pga_type, "_", tissue))
  })
  
  return(setNames(pga_capped, colnames(obj@expr.data)))
})
dev.off()

pga_tissues <- setNames(pga_tissues, samples)
saveRDS(pga_tissues, file = file.path(infercnv_outdir, "pga_tisses.rds"))
pga_tissues <- readRDS(file.path(infercnv_outdir, "pga_tisses.rds"))

# Remove duplicate cells that stem from same ehalthy PBMCs
pga_tissues <- unlist(pga_tissues)
dup_idx <- which(duplicated(gsub("^.*\\.", "", names(pga_tissues))))
pga_tissues <- pga_tissues[-dup_idx]
pga_tissues <- data.frame("pga"=pga_tissues, 
                          "cellids"=gsub("^.*\\.", "", names(pga_tissues)))

# Instantiate the PGA data to be added to seurat object meta, then fill in PGA values
pga_meta <- data.frame("vals"=rep(0, length(seu.harmony$seurat_clusters)),
                       "cellids"=names(seu.harmony$seurat_clusters))
pga_meta <- merge(pga_meta, pga_tissues, all.x=TRUE, by='cellids')
pga_meta <- setNames(rowSums(pga_meta[,2:3], na.rm = TRUE), pga_meta$cellids)
seu.harmony$PGA <- pga_meta

pga_spl <- split(seu.harmony$PGA, seu.harmony$seurat_clusters_ord)
cl_cnts <- sapply(split(clusters_anno, clusters_anno$anno),nrow)
cols <- brewer.pal(n = length(cl_cnts), name = "Set3")
pdf(file.path(outdir, "features", "pga_boxplots.pdf"), width = 9, height = 5)
par(mfrow=c(2,1))
par(mar=c(0, 4.1, 10.1, 2.1))
barplot(sapply(pga_spl, length), col=rep(cols, cl_cnts), 
        xlab='', xaxt='n', ylab='n', las=2, cex.axis=0.6)

par(mar=c(5.1, 4.1, 0.1, 2.1))
box.x <- boxplot(pga_spl, ylim=c(0,0.6), las=2, cex.axis=0.6, pch=16,
        ylab='PGA', xlab='Clusters', border="black",
        col="lightgrey", outline=FALSE)

# +/-1.58 IQR/sqrt(n)
pga_spl_hiiqr <- lapply(seq_along(pga_spl), function(idx){
  pga_x <- pga_spl[[idx]]
  pga_x[which(pga_x > box.x$stats[5,idx])]
})
stripchart(pga_spl_hiiqr,              # Data
           method = "jitter", jitter=0.25,
           pch = 20,          # Pch symbols
           border=rep(cols, cl_cnts),
           col = alpha(rep(cols, cl_cnts), 0.4),           # Color of the symbol
           xlab='Clusters', ylab='PGA',
           vertical = TRUE, add=TRUE)

dev.off()

# Feature plots of AhR expression and PGA blended together
pdf(file.path(outdir, "features", "nc-featplots_pga.pdf"), width = 12, height = 4)
FeaturePlot(seu.harmony, features = c('AHR', 'PGA'), pt.size = 0.5, ncol = 1,
            cols=c("#f7f7f7", "#c51b7d", "#4d9221"), order=TRUE, keep.scale=NULL, label=FALSE,
            blend=TRUE)
dev.off()

# Dot plot of AhR expression and PGA using Harmony Clusters
pdf(file.path(outdir, "features", "nc-dotplot_pga.pdf"), width=6, height = 12)
Idents(seu.harmony) <- seu.harmony$seurat_clusters_ord
DotPlot(seu.harmony, features = c('AHR', 'PGA'), # unique(unlist(natcan_markers)),
        cols = c('grey', 'darkred'), dot.scale = 10) + 
  RotatedAxis() + FontSize(x.text = 17, y.text = 17)
dev.off()


x <- table(seu.harmony$orig.ident)
ids <- unique(gsub("_[a-zA-Z0-9]+$", "", names(x)))
sapply(ids, function(id) sum(x[grep(id, names(x))]))