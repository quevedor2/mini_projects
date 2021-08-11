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
library(ggplot2)

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

################
#### 2. GSE ####
ds_idx <- 2
dataset   <- names(datasets)[ds_idx]
file      <- datasets[[ds_idx]][1]
analysis  <- datasets[[ds_idx]][2]
gene      <- datasets[[ds_idx]][3]


################################
#### 4a. Load in Seurat data ####
harmony_file <- 'harmony2' # harmony or harmony2
seu.harmony <- SeuratDisk::LoadH5Seurat(file.path(dataset, paste0(file, "_", harmony_file, ".h5seurat")))

####################################
#### 4b. Load Annotated Clusters ####
harmony_file <- 'harmony2' # harmony or harmony2
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
  print("Reading harmony markers")
  load(file=file.path(dataset, paste0(file, "_", harmony_file, "-markers.rda")))
}

###########################################
#### 4c. Annotate the harmony clusters ####
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

#################################################
#### 4d. Refine automatic cluster annotation ####
if(annotation_method=='scsa' & file.exists(file.path(outdir_plus, 'annotate.txt'))){
  if(!file.exists(file.path(outdir_plus, "seurat_markers.manual.csv"))){
    # WARNING: This sections is heavily manually curated, as such, there are hardcoded segments
    # manually annotate some clusters after inspecting UMAP and SCSA annotations
    manual_ids <- setNames(c(1,3,2,1,3,1,2,14,2,3,3,2,3,14,7,15,7,7,
                             14,10,14,12,14,6,14,7,14,8,8,8,11,1,16,
                             14,8,5,3,13,1,4), df$V1)
    seu.harmony$manual_clusters <- manual_ids[as.character(seu.harmony$seurat_clusters)]
    
    # Find differential markers between clusters that I was unable to disambiguate between
    Idents(seu.harmony) <- seu.harmony$manual_clusters
    grp_comp <- list("monocytes-macrophage-neutrophil"=c('14', '1'),
                     'nk-tcell'=c('2', '3'),
                     'plasma-bcell'=c('4', '15'),
                     'b-plasma_b-plasma'=c('4', '12'),
                     'b-plasma_bcell'=c('12', '15'))
    grp_markers <- lapply(grp_comp, function(grps){
      print(grps)
      markers <- FindMarkers(object = seu.harmony, only.pos = TRUE, min.pct = 0.25, 
                             thresh.use = 0.25, ident.1 = grps[1], ident.2 = grps[2])
      return(markers)
    })
    
    # Reduce the list to a dataframe and format for SCSA input
    g_ids <- sapply(grp_comp, paste, collapse="_")
    for (idx in seq_along(grp_comp)) {
      grp_markers[[idx]]$cluster <- g_ids[idx]
      grp_markers[[idx]]$gene <- rownames(grp_markers[[idx]])
      colnames(grp_markers[[idx]])[2] <- 'avg_logFC'
    }
    grpmark <- do.call(rbind, grp_markers)
    write.table(grpmark, file=file.path(outdir_plus, "seurat_markers.manual.csv"), 
                sep=",", col.names = TRUE, row.names=FALSE, quote = FALSE)
  }
  
  ## After running SCSA on the updated seurat markers and referencing the Naturecancer Paper
  cluster_ids <- setNames(c(
    "Neutrophil", "T-Cell", "Natural killer cell", "Neutrophil", "T-Cell", "Neutrophil", 
    "Natural killer cell", "Macrophage", "Natural killer cell", "T-Cell", "T-Cell", 
    "T-Cell", "T-Cell", "Macrophage", "Epithelial cell", "B-cell", "Epithelial cell", 
    "Epithelial cell", "Macrophage", "Mast cell", "Macrophage", "Plasma cell", 
    "Macrophage", "Acinar cell", "Macrophage", "Epithelial cell", "Dendritic cells", 
    "Mesenchymal stem cell", "Mesenchymal stem cell", "Mesenchymal stem cell", 
    "NK T-cell", "Neutrophil", "Endothelial cell", "Macrophage", "Mesenchymal stem cell", 
    "Plasmacytoid dendritic cell", "Natural killer cell", "Megakaryocyte", "Neutrophil", 
    "Hematopoietic stem cell"), df$V1)
}

# Add the annotated cell types to Seurat metadata
seu.harmony$anno_clusters <- cluster_ids[as.character(seu.harmony$seurat_clusters)]


#### 6. InferCNV ####
### inferCNV
load("/cluster/home/quever/downloads/infercnv/data/infercnv_genes_example.rda")

patients <- seu.harmony$orig.ident
idents    <- unique(patients)
idents_t  <- idents[grep("Healthy|Norm|PBMC", idents, invert = TRUE)]
seed <- 1234

## Get PDAC-Tissues
id <- idents_t[as.integer(args[1])]
message(paste0("Looking at ", id, " (", match(id, idents_t), "/", length(idents_t), ")"))
patient_idx <- which(patients %in% id)

## Generate a standardized random set of normal cells
set.seed(seed)
norm_style <- args[2]
if(norm_style == 'random'){
  message("Using frankenstein of healthy PBMCs and adjacent normals")
  norm_idx        <- grep("Healthy|Norm", patients)
  ident_norm_idx  <- grep("Healthy|Norm", idents)
  norm_cells	  <- sapply(idents, function(id) sum(patients %in% id))[ident_norm_idx]
  mean_norm_cells <- ceiling(mean(norm_cells))
  norm_rand_idx   <- sort(sample(norm_idx, size = mean(norm_cells), replace = FALSE))
  norm_ids        <- names(patients)[norm_rand_idx]
  norm_idx        <- norm_rand_idx
} else if(norm_style == 'matched_pbmc'){
  pbmc_normal   <- gsub("TISSUE", "PBMC", id)
  ident_norm_idx <- any(grepl(paste0("^", pbmc_normal, "$"), idents))
  if(ident_norm_idx){
    # used matched normal if it exists
    message("Using matched PBMCs")
    norm_idx        <- which(patients == pbmc_normal)
  } else {
    # default to healthy pbmcs
    message("Using random subset of healthy PBMCs")
    norm_idx <- grep("Healthy_PBMC_", patients)
    ident_norm_idx <- grepl("Healthy_PBMC_", idents)
    norm_idx <- sort(sample(norm_idx, size=ceiling(length(norm_idx)/sum(ident_norm_idx))))
  }
} else if(norm_style == 'healthy_pbmc'){
  message("Using random subset of healthy PBMCs")
  norm_idx <- grep("Healthy_PBMC_", patients)
  ident_norm_idx <- grepl("Healthy_PBMC_", idents)
  norm_idx <- sort(sample(norm_idx, size=ceiling(length(norm_idx)/sum(ident_norm_idx))))
}


# Remove normal clusters under a certain size
min_clust_size <- 5
small_clusts <- table(seu.harmony$anno_clusters[c(norm_idx)]) < min_clust_size
if(any(small_clusts)){
  rm_clust_id <- names(which(small_clusts))
  rm_clust_idx <- !is.na(match(seu.harmony$anno_clusters[c(norm_idx)], rm_clust_id))
  norm_idx <- norm_idx[-which(rm_clust_idx)]
}

small_clusts <- table(seu.harmony$seurat_clusters[c(patient_idx)]) < min_clust_size
if(any(small_clusts)){
  rm_clust_id <- names(which(small_clusts))
  rm_clust_idx <- !is.na(match(seu.harmony$seurat_clusters[c(patient_idx)], rm_clust_id))
  patient_idx <- patient_idx[-which(rm_clust_idx)]
}

# Append normal random cells to tumor cells
exprs <- seu.harmony@assays$RNA@counts[,c(patient_idx, norm_idx)]
annos_tumr <- seu.harmony$seurat_clusters[c(patient_idx)]
annos_norm <- seu.harmony$anno_clusters[c(norm_idx)]
# annos <- seu.harmony$anno_clusters[c(patient_idx, norm_idx)]
annos <- setNames(c(as.character(annos_tumr), as.character(annos_norm)),
                  c(names(annos_tumr), names(annos_norm)))
print(table(gsub("_[ACGT]*-[0-9]$", "", names(annos))))
normal_clusters <- paste0('Normal_', (annos_norm))
annos[-c(1:length(patient_idx))] <- normal_clusters

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=exprs,
                                    annotations_file=data.frame(annos),
                                    gene_order_file=infercnv_genes_example,
                                    ref_group_names=normal_clusters)

outdir <- file.path("output_dir", norm_style, id)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                              out_dir=outdir,  # dir is auto-created for storing outputs
                              cluster_by_groups=T,   # cluster
                              denoise=T,
                              HMM=T,
                              resume_mode=FALSE)

