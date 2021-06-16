# AhR Analysis in PDACs
# Rene Quevedo
# June 16-2021
#
## "Question- Given the wide availability of single cell sequencing 
## datasets of human pancreatic cancer, the authors should address the 
## cell type specific expression of AhR in human tumors."
##
## "dig up some public scRNA seq data (preferably from high impact papers 
## if possible). And look for AhR expression data for all cells and provide 
## a cell type ranked assessment of relative AhR expression"

library(Seurat)
library(SeuratObject)
library(SeuratDisk)
###############################################
#### Functions: Pre-assemble counts matrix ####
preAssembleMatrices <- function(dataset, file){
  switch(dataset,
         "PRJCA001063"=.magliano(file),
         "GSE155698"=.wu(file))
}

.magliano <- function(file='PRJCA001063_CRC_besca2.annotated'){
  ## PRJCA001063
  ## Marina Pasca di Magliano: (Nature Cancer - PRJCA001063) Multimodal mapping of the tumor and peripheral blood immune landscape in human pancreatic cancer
  # 16 treatment-naive PDA samples (6 surgical, 10 fine-needle biopsy)
  # 3 adjacent normal non-malignant pancreas samples (1196, 1258, 19732; duodenal adenoma, ampullary carcinoma, tissue adjacent to PDA)
  # 16 PBMCs from PDA samples
  # 4 PBMCs from 4 healthy donors
  
  Convert(paste0(file, 'h5ad'), dest = "h5seurat", overwrite = TRUE)
}

.wu <- function(dataset='GSE155698'){
  ## GSE155698
  # Wenming Wu: (Cell Research – GSE155698) Single-cell RNA-seq highlights intra-tumoral heterogeneity and malignant progression in pancreatic ductal adenocarcinoma
  # 24 PDAC tumor samples
  # Control adjacent normal pancreas without visible inflammation: 11 (3 non-pancreatic tumor [bile duct tumors, duodenal tumors], and 8 non-malignant pancreatic tumors [pancreatic cyst]).
  dirs <- grep(list.files(dataset), pattern="\\.tar", invert=TRUE, value=TRUE)
  seurobjs <- lapply(dirs, function(d){
    dir <- file.path(dataset, d, "filtered_feature_bc_matrix")
    expression_matrix <- ReadMtx(
      mtx = file.path(dir, "matrix.mtx.gz"), 
      features = file.path(dir, "features.tsv.gz"),
      cells = file.path(dir, "barcodes.tsv.gz")
    )
    CreateSeuratObject(counts=expression_matrix, project=d)
  })
  seurobj <- merge(seurobjs[[1]], y=seurobjs[-1], 
                   add.cell.ids=dirs, project=dataset)
  SeuratDisk::SaveH5Seurat(seurobj, filename=file.path(dataset, paste0(dataset, ".h5seurat")), overwrite = TRUE)
}

###############################################
#### Functions: Pre-assemble counts matrix ####
mapToPMBC <- function(seu, analysis, spotcheck_clusters=TRUE){
  if(analysis == 'sct'){
    seu <- SCTransform(seu, variable.features.n = 2000, verbose = FALSE)
  } else if(analysis == "steps_nonorm") {
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
    # seu <- ScaleData(seu, verbose = FALSE)
    id_map <- c('T'='Tumor', 'N'='Normal')
    seu$sample.type <- id_map[seu$CONDITION]
  } else if(analysis == "steps_norm") {
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
    seu <- ScaleData(seu, verbose = FALSE)
    seu$sample.type <- gsub("_[0-9a-zA-Z]*$", "", seu$orig.ident)
    
    if(spotcheck_clusters){
      seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
      seu <- RunUMAP(seu, reduction = "pca", dims = 1:30)
      seu <- FindNeighbors(seu, reduction = "pca", dims = 1:30)
      seu <- FindClusters(seu, resolution = 0.5)
      
      # pdf("~/xfer/gse_sampletypes.pdf")
      # p1 <- DimPlot(seu, reduction = "umap", group.by = "sample.type")
      # p2 <- DimPlot(seu, reduction = "umap", label = TRUE, repel = TRUE)
      # p1
      # dev.off()
    }
  } else if(analysis == 'integrate'){
    seu <- lapply(seus, function(seu){
      seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
    })
    features <- SelectIntegrationFeatures(object.list = seu)
    seu.anchors <- FindIntegrationAnchors(object.list = seu, anchor.features = features)
    
    # this command creates an 'integrated' data assay
    seu.combined <- IntegrateData(anchorset = seu.anchors)
    seu.combined <- ScaleData(seu.combined, verbose = FALSE)
  }
  # split(seu, seu$CONDITION)                                                                                                   
  
  anchors <- FindTransferAnchors(reference = reference, query = seu,
                                 normalization.method = "SCT", 
                                 reference.reduction = "spca",
                                 dims = 1:50)
  seu <- MapQuery(anchorset = anchors, query = seu,
                  reference = reference,
                  refdata = list(celltype.l1 = "celltype.l1",
                                 celltype.l2 = "celltype.l2",
                                 predicted_ADT = "ADT"),
                  reference.reduction = "spca", 
                  reduction.model = "wnn.umap")
  
  ## Assign cell type labels and tumor/normal default for 'other'
  other_idx <- seu$predicted.celltype.l1 == 'other'
  seu$l1_l2 <- seu$predicted.celltype.l2
  seu$l1 <- seu$predicted.celltype.l1
  seu$l1_l2[which(other_idx)] <- seu$sample.type[which(other_idx)]
  seu$l1[which(other_idx)]    <- seu$sample.type[which(other_idx)]
  
  
  ## Rearrange based on AHR expression
  x <- split(GetAssayData(seu)['AHR',], seu$l1_l2)
  median_ord <- order(sapply(unique(seu$l1_l2), function(i) median(x[[i]][which(x[[i]]!=0)],na.rm=TRUE)), na.last=TRUE)
  seu$l1_l2 <- factor(seu$l1_l2, levels=unique(seu$l1_l2)[median_ord])
  x <- split(GetAssayData(seu)['AHR',], seu$l1)
  median_ord <- order(sapply(unique(seu$l1), function(i) median(x[[i]][which(x[[i]]!=0)],na.rm=TRUE)), na.last=TRUE)
  seu$l1 <- factor(seu$l1, levels=unique(seu$l1)[median_ord])
  
  return(seu)
}

vizGene <- function(seu, gene='AHR', gene_list=c("CD14", "ITGAM", "CD68")){
  ## Dimension plot
  p1 = DimPlot(seu, reduction = "ref.umap", group.by = "l1", label = TRUE, 
               label.size = 3, repel = TRUE) # + NoLegend()
  p1
  p2 = DimPlot(seu, reduction = "ref.umap", group.by = "l1_l2", label = TRUE, 
               label.size = 3, repel = TRUE) # + NoLegend()
  p2
  
  CombinePlots(plots = list(p1 + NoLegend(), p2  + NoLegend()))
  
  
  ## Ridge plot
  p1 <- RidgePlot(seu, features = gene, group.by = "l1", same.y.lims=TRUE) +
    theme(legend.position = 'none')
  p2 <- RidgePlot(seu, features = gene, group.by = "l1_l2", same.y.lims=TRUE) +
    theme(legend.position = 'none')
  CombinePlots(plots = list(p1, p2))
  
  ## Violin plot
  Idents(seu) <- 'l1'
  p1 <- VlnPlot(seu, features = c(gene), sort = TRUE) + NoLegend()
  Idents(seu) <- 'l1_l2'
  p2 <- VlnPlot(seu, features = c(gene), sort = TRUE) + NoLegend()
  CombinePlots(plots = list(p1, p2))
  
  ## Feature plot
  p1 <- FeaturePlot(seu, features = c(gene, gene_list), reduction = 'ref.umap',
                    min.cutoff = 'q10', max.cutoff = 'q99', cols = c("lightgrey", "darkgreen") ,
                    ncol = 3)
  p1
}

##############
#### Main ####
# 1. Convert datasets to Seurat h5 data
setwd("/cluster/projects/pughlab/projects/cancer_cell_lines/tmp/scrna/pancrea_ahr/datasets")
if(!file.exists(file.path('PRJCA001063', 'PRJCA001063_CRC_besca2.annotated.h5seurat'))){
  preAssembleMatrices('PRJCA001063', 'PRJCA001063_CRC_besca2.annotated')
}
if(!.file.exists(file.path('GSE155698', 'GSE155698.h5seurat'))){
  preAssembleMatrices('GSE155698', 'GSE155698')
}

# 2. Read in the reference PBMC CITE-seq data: 162,000 cells, 228 antibodies
##  “Integrated analysis of multimodal single-cell data” (https://pubmed.ncbi.nlm.nih.gov/34062119/),
## https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat
dataset <- 'pbmc_multimodal'
reference <- LoadH5Seurat(file.path(dataset, paste0(dataset, ".h5seurat")))

# 3. Read in and process each query dataset
outdir   <- "~/xfer"
datasets <- list('PRJCA001063'=c('PRJCA001063_CRC_besca2.annotated', 
                                 'steps_nonorm',
                                 'AHR'),
                 'GSE155698'=c('GSE155698',
                               'steps_norm',
                               'AHR'))

seu_maps <- lapply(seq_along(datasets), function(ds_idx){
  dataset   <- names(datasets)[ds_idx]
  file      <- datasets[[ds_idx]][1]
  analysis  <- datasets[[ds_idx]][2]
  gene      <- datasets[[ds_idx]][3]
  
  ## Load data and use reference PBMC to map PBMC labels
  cnts <- SeuratDisk::LoadH5Seurat(file.path(dataset, paste0(file, ".h5seurat")))
  seu <- mapToPMBC(cnts, analysis)
  
  ## Use wilcox test to test gene vs all expression
  wilcox <- sapply(levels(seu$l1_l2), function(i){
    unlist(FindMarkers(seu, features=gene, test.use='wilcox',
                       ident.1 = i, ident.2 = NULL, group.by = 'l1_l2'))
  })
  wilcox <- do.call(rbind, wilcox)
  
  ## Visualize the genes using the mapped query dataset
  pdf(file.path(outdir, paste0(gsub("_.*", "", dataset), ".pdf")), width = 12)
  vizGene(seu, gene)
  dev.off()
  
  ## Validation check if dataset is PRJCA001063 due to the presence of cell labels
  if(dataset=='PRJCA001063'){
    ## Visualize reported cell types to predicted cell types:
    pdf(file.path(outdir, paste0(gsub("_.*", "", dataset), "_celltype.pdf")))
    p1 = DimPlot(seu, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
    p2 = DimPlot(seu, reduction = "ref.umap", group.by='celltype1', label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
    p1 + p2
    dev.off()
    
    ## Precision-recall check of predicted cell labels to actual labels
    hardcodemap <- list('B cell' = 'B',
                        'blood vessel endothelial cell' = 'other',
                        'CD4-positive, alpha-beta T cell' = 'CD4 T',
                        'CD8-positive, alpha-beta T cell' = 'CD8 T',
                        'fibroblast' = 'other',
                        'macrophage' = 'Mono', 
                        'myeloid dendritic cell' = 'DC',
                        'neural cell' = 'other',
                        'pancreatic acinar cell' = 'other',
                        'pancreatic ductal cell' = 'other',
                        'pancreatic epsilon cell' = 'other',
                        'pancreatic stellate cell' = 'other',
                        'plasma cell' = 'B',
                        'type B pancreatic cell' = 'other')
    seu_pred <- factor(as.character(seu$predicted.celltype.l1))
    seu_actual <- as.character(seu$celltype2)
    seu_actual_map <- factor(unlist(hardcodemap[seu_actual]), levels=levels(seu_pred))
    cm <- caret::confusionMatrix(seu_pred, reference = seu_actual_map)
    
    write.table(cm$byClass[,'F1'], file=file.path(outdir, "f1_metrics.csv"),
                sep=",", col.names = TRUE, row.names = TRUE, quote = FALSE)
    write.table(table(seu_pred, seu_actual), file=file.path(outdir, "actual_contigency.csv"),
                sep=",", col.names = TRUE, row.names = TRUE, quote = FALSE)
    write.table(table(seu_pred, seu_actual_map), file=file.path(outdir, "mapped_contigency.csv"),
                sep=",", col.names = TRUE, row.names = TRUE, quote = FALSE)
  }
  
  return(list('seurat'=seu,
              'tests'=wilcox))
})

  


  