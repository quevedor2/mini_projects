library(clusterProfiler)
library(msigdbr)
library(AUCell)
# library(SoupX)
library(monocle3)
library(flexmix)
library(scater)
library(miQC)
library(vegan)
library(DoubletFinder)
# library(cellpypes)
library(Seurat)
library(BPCells)
library(ggplot2)
library(cowplot)
library(harmony)
library(scales)
library(scRNAseq)
library(SingleR)
library(scuttle)
# library(celldex)
library(immunarch)
library(vegan)
library(plyr)
library(circlize)
library(SeuratWrappers)
library(SCPA)
# library(testSctpa)
library(ProjecTILs)
# library(EnsDb.Hsapiens.v86)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1e9)


PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/aml_abbvie'
datadir <- file.path(PDIR, "data")
outdir <- file.path(PDIR, "results")
proj <- 'aml_abbvie'
setwd(PDIR)

samples <- list.files(file.path(datadir, 'cellranger'), pattern="^Sample")

###################
#### Functions 
##-- Recurring functions ----
source("~/git/mini_projects/mini_functions/geneMap.R")
source("~/git/mini_projects/mini_functions/singlecell/scrna.R")
source("~/git/mini_projects/mini_functions/singlecell/pagoda2.R")
source("~/git/mini_projects/mini_functions/singlecell/as.SingleCellExperiment.Seurat5.R")
gm <- geneMap('Homo sapiens')

###################################
#### 1. Create Seurat object  ####
dir.create(file.path("results", "seurat_obj"), 
           recursive = TRUE, showWarnings = FALSE)
seu_rds_f <- file.path(outdir, "seurat_obj", "1_unfiltered.rds")
overwrite <- T
seus <- if(file.exists(seu_rds_f) & !overwrite) readRDS(file=seu_rds_f) else NULL

# Read in and create the seurat objects
seus <- lapply(samples, function(s){
  if(s %in% names(seus)){
    print(paste0(s, ": Returning existing seurat object from a previous iteration"))
    seu <- seus[[s]]
  } else {
    print(paste0(s, ": Processing..."))
    print(file.path(datadir, 'rna', 'filt', s))
    mtx <- Read10X(data.dir = file.path(datadir, 'rna', 'filt', s), strip.suffix=TRUE)
    # raw_mtx <- Read10X(data.dir = file.path(datadir, 'rna', 'raw', s), strip.suffix=TRUE)
    seu <- CreateSeuratObject(counts = mtx, project = s)
  }
  return(seu)
})
names(seus) <- samples
seu <- merge(seus[[1]], seus[-1])

###########################
#### SCT Normalization ####
# seu <- split(DietSeurat(seu, layers='counts'), seu$orig.ident)
percent.mt.cutoff <- 15
ncount.min <- 500
ncount.max <- Inf
nfeature.min <- 0
nfeature.max <- 6500

# normalize data with SCTransform()
seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")
seu.raw <- seu
seu <- subset(seu.raw, subset = 
                nCount_RNA > ncount.min &
                nCount_RNA < ncount.max &
                nFeature_RNA > nfeature.min & 
                nFeature_RNA < nfeature.max & 
                percent.mt < percent.mt.cutoff)

seu <- SCTransform(seu,
                   assay = 'RNA',
                   new.assay.name = 'SCT',
                   vst.flavor='v2',
                   vars.to.regress = c('percent.mt'))


# seu <- CellCycleScoring(seu, 
#                         s.features = stringr::str_to_title(cc.genes$s.genes), 
#                         g2m.features = stringr::str_to_title(cc.genes$g2m.genes), 
#                         set.ident = TRUE,
#                         assay='SCT')
# 
# # seu <- split(seu, seu$sample)
# seu <- SCTransform(seu, vst.flavor = "v2", verbose = TRUE,
#                    assay='RNA', 
#                    new.assay.name='SCT',
#                    # variable.features.n = 1000,
#                    vars.to.regress = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 
#                                        "S.Score", "G2M.Score"))
# seu <- FindVariableFeatures(seu)

seu <- seu %>% 
  ScaleData(.)  %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE, 
          reduction.name='umap.unintegrated')
seu$orig.ident <- factor(seu$orig.ident)

seu <- IntegrateLayers(
  object = seu, 
  method = HarmonyIntegration,
  normalization.method = "SCT",
  new.reduction = "integrated.harmony",
  verbose = T,
  group.by='orig.ident'
) 

set.seed(1234)
DefaultAssay(seu) <- 'SCT'
seu <- seu %>% 
  FindNeighbors(., reduction = "integrated.harmony", dims = 1:30) %>%
  FindClusters(., resolution = 0.8, cluster.name = "harmony_clusters") %>% 
  RunUMAP(., reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony")

seu <- PrepSCTFindMarkers(seu)


#--- a) SingleR: Monaco ----
# RNA preprocess and cell-cycle score
DefaultAssay(seu) <- 'SCT'
references <- list("hpca"="HumanPrimaryCellAtlasData.rds",
                   "encode"="BlueprintEncodeData.rds",
                   "monaco"="Monaco.rds")

for(database_id in names(references)[3]){
  bed.se <- readRDS(file.path("/cluster/projects/mcgahalab/ref/scrna/scRNAseq_datasets", 
                              references[[database_id]])) 
  
  for(lbl in c('label.fine', 'label.main')){
    for(clus in c(TRUE)){
      rds_file <- paste0("celldex_", database_id, ".", gsub("label.", "", lbl),
                         ".", if(clus) 'cluster' else 'cell', ".rds")
      id <- gsub(paste0("^.*", database_id, "(.*).rds"), 
                 paste0(database_id, "\\1"), 
                 rds_file)
      if(!file.exists(file.path(outdir, "annotation", rds_file)) | overwrite){
        print(paste0(" Annotating: ", rds_file))
        singler_anno <- SingleR(test=GetAssayData(seu, assay = "SCT", slot = "counts"), 
                                ref=bed.se, 
                                assay.type.test=1, labels=bed.se[[lbl]],
                                de.method="wilcox", genes='sd',
                                clusters=if(clus) seu$seurat_clusters else NULL)
        saveRDS(singler_anno, file=file.path(outdir, "annotation", rds_file))
      } else {
        print(paste0(" Reading in annotation: ", rds_file))
        singler_anno <- readRDS(file.path(outdir, "annotation", rds_file))
        # pdf("~/xfer/test3.pdf", width=15)
        # plotScoreHeatmap(singler_anno, show_colnames=TRUE)
        # dev.off()
      }
      
      if(clus){
        cluster_ids <- setNames(singler_anno$pruned.labels, 
                                as.character(rownames(singler_anno)))
        seu@meta.data[,id] <- cluster_ids[as.character(seu$seurat_clusters)]
      } else {
        seu@meta.data[,id] <- singler_anno$pruned.labels
      }
    }
  }
}


DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu)
seu <- JoinLayers(seu)


saveRDS(seu, file=file.path(outdir, "seurat_obj", "seu_anno.rds"))
seu <- readRDS(file=file.path(outdir, "seurat_obj", "seu_anno.rds"))


## make loupe
cols <- c('orig.ident', 'seurat_clusters', 
          'monaco.fine.cluster', 'monaco.main.cluster')
seu@meta.data <- seu@meta.data[,cols]
seu[['RNA']] <- CreateAssayObject(counts=seu@assays$SCT$counts)
DefaultAssay(seu) <- 'RNA'
loupeR::create_loupe_from_seurat(seu, output_dir=file.path(outdir, "loupe"), 
                                 output_name='aml_abbvie', 
                                 force=T, tmpdir=PDIR)
file.copy(file.path(outdir, 'loupe', "aml_abbvie.cloupe"),
          to = "~/xfer", overwrite = T)
cat("xfer pbmc_multimodal.cloupe\n")


# Differential expression to find DEGs per cluster
seu <- PrepSCTFindMarkers(seu)
sctallmarkers  <- FindAllMarkers(seu, assay = "SCT", verbose = FALSE)
degs <- lapply(as.character(unique(seu$harmony_clusters)), function(uniq.clus){
  seusub <- subset(seu, cells=Cells(seu)[(seu$harmony_clusters == uniq.clus)])
  Idents(seusub) <- 'orig.ident'
  
  Seurat::FindMarkers(seusub, assay='SCT', ident.1 = "HIF1a_FF", ident.2='HIF1a_FF_LMC', 
                      group.by = 'orig.ident', recorrect_umi=F) %>%
    mutate("Cluster"=uniq.clus,
           'ident.1'='HIF1a_FF',
           'ident.2'='HIF1a_FF_LMC')
})
degs_df <- do.call(rbind, degs)
write.table(degs_df, file=file.path(outdir, "deg", "origident_by_cluster.csv"),
            col.names = T, row.names = F, quote = F)
write.table(sctallmarkers, file=file.path(outdir, "deg", "allmarkers_harmonycluster.csv"),
            col.names = T, row.names = F, quote = F)
file.copy(file.path(outdir, "deg", "allmarkers_harmonycluster.csv"), to="~/xfer", overwrite = T)
file.copy(file.path(outdir, "deg", "origident_by_cluster.csv"), to="~/xfer", overwrite = T)

