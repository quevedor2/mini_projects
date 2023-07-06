library(velocyto)
library(Seurat)

##############################
#### 9. Velocity analysis ####
outrds <- file.path(datadir, "seurat_obj", "seu_velocyto.rds")
dir.create(file.path(outdir, "scVelo"), showWarnings = F)

## This analysis assumes that velocyto/0.17.17 was run on your 10x data already:
# gtf='/cluster/tools/data/commondata/cellranger/refdata-gex-mm10-2020-A/genes/genes.gtf'
# mask_gtf='/cluster/projects/mcgahalab/ref/ucsc/repeat_masker/mm10_rmsk.gtf'
# 
# velocyto run10x \
# -m ${mask_gtf} \
# --samtools-threads 1 \
# --samtools-memory 20000 \
# ${sampledir} \
# ${gtf}
# 
# velo_obj <- readRDS(file=outrds)


velocyto_dir <- '/path/to/loom_files'

## Read in the velocyto loom files
sample_ids <- list.files(velocyto_dir)
seu_velo_raw <- lapply(setNames(sample_ids,sample_ids), function(sid){
  f <- list.files(file.path(velocyto_dir, sid))
  loom.data <- ReadVelocity(file = file.path(velocyto_dir, sid, f)) %>%
    as.Seurat(.)
  loom.data$orig.ident <- sid
  loom.data <- RenameCells(loom.data, 
                           new.names=gsub("^McGaha__(MDSC_DAB|MDSC_DMSO)\\:([ACGT]*).*", "\\2_\\1", Cells(loom.data)))
  return(loom.data)
})

## Preprocess the seurat velocyto files
seu_velo <- merge(seu_velo_raw[[1]], seu_velo_raw[-1], 
                  project = 'DMSO_DAB')
seu_velo_l <- list("Merged"=seu_velo,
                   "DAB"=seu_velo_raw$MDSC_DAB,
                   "DMSO"=seu_velo_raw$MDSC_DMSO)
rm(seu_velo, seu_velo_raw); gc()
# seu_velo <- seu_velo_raw[[1]]

make.small <- FALSE
x <- lapply(names(seu_velo_l), function(seu_velo_id){
  print(seu_velo_id)
  seu_velo <- seu_velo_l[[seu_velo_id]]
  DefaultAssay(object = seu_velo) <- "spliced"
  if(make.small){
    seu_velo <- subsetSeu(seu_velo,colid='orig.ident', n=3000) %>%
      subset(., nFeature_spliced>100)
  }
  seu_small <- subset(seu, cells=Cells(seu_velo))
  seu_velo <- subset(seu_velo, cells=Cells(seu_small))
  
  # Normalize and cluster
  set.seed(seed)
  seu_velo[["RNA"]] <- seu_velo[["spliced"]]
  seu_velo <- seu_velo %>%
    SCTransform(.) %>%
    RunPCA %>%
    RunUMAP(., dims=1:30) %>%
    FindNeighbors(., dims = 1:30) %>%
    FindClusters(.)
  DefaultAssay(seu_velo) <- "RNA"
  seu_velo[['umap_orig']] <- CreateDimReducObject(embeddings = Embeddings(seu_small, reduction='umap'),
                                                  key = 'UMAP_', assay = 'RNA')
  seu_velo$functional_cluster <- seu_small$mnn.functional.cluster
  
  SeuratDisk::SaveH5Seurat(seu_velo, 
                           filename = file.path(outdir, "scVelo", paste0("seu_velo.", seu_velo_id, ".h5Seurat")))
  cwd <- getwd()
  setwd(file.path(outdir, "scVelo"))
  SeuratDisk::Convert(source = paste0("seu_velo.", seu_velo_id, ".h5Seurat"), 
                      dest = "h5ad")
  setwd(cwd)
  return(NULL)
})