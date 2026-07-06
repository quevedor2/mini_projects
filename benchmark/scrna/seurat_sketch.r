# renv::load("/cluster/home/quever/downloads/renvs/")
# renv::load("/cluster/projects/mcgahalab/envs/renvs/seuratv5_v2")
library(clusterProfiler)
library(Seurat)
library(ggplot2)
library(cowplot)
library(harmony)
library(scales)
library(scRNAseq)
library(SingleR)
library(tidyverse)
library(SeuratWrappers)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1e9)


PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/benchmark/seurat_sketch'
datadir <- file.path(PDIR, "data", "cellranger")
xeniumdir <- file.path(PDIR, "data", "xenium")
outdir <- file.path(PDIR, "results")
proj <- 'ipmn_pdac_sketch'
setwd(PDIR)

sketch_range <- c(100, 500, 1000, 2000, 3000, 4000, 5000)
variable_features_n <- 3000
n_samples <- 3


###################
#### Functions 
##-- Recurring functions ----
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")
source("~/git/mini_projects/mini_functions/gsea2CytoscapeFormat.R")
source("~/git/mini_projects/mini_functions/geneMap.R")
source("~/git/mini_projects/mini_functions/singlecell/scrna.R")
source("~/git/mini_projects/mini_functions/singlecell/pagoda2.R")
source("~/git/mini_projects/mini_functions/singlecell/as.SingleCellExperiment.Seurat5.R")
source("~/git/mini_projects/mini_functions/singlecell/gptcelltype.R")
gm <- geneMap(species='Homo sapiens')

##-- GSE functions ----
ProjectDataSCT <- function(object, refdata, full_paths,
                           sketched.reduction='pca',
                           sketched.assay = 'sketch',
                           assay = 'SCT', k=5,
                           dims=1:30){
  pca_rotation  <- Loadings(object, reduction = sketched.reduction)
  genes_use     <- rownames(pca_rotation)
  sketch_cells  <- Cells(object[[sketched.assay]])
  sketch_emb    <- Embeddings(object, sketched.reduction)[sketch_cells, dims]
  sketch_labels <- as.character(object@meta.data[sketch_cells,refdata[[1]]])
  names(sketch_labels) <- sketch_cells
  
  # Projected cluster ids
  all_labels <- setNames(rep(NA, ncol(object)), colnames(object))
  all_labels[sketch_cells] <- sketch_labels  # sketch cells labeled directly
  
  for(sid in names(full_paths)){
    message("Projecting: ", sid)
    seu_full  <- qs2::qs_read(full_paths[[sid]])
    DefaultAssay(seu_full) <- assay
    
    # Scale full cells using PCA genes
    genes_present <- intersect(genes_use, rownames(seu_full))
    seu_full <- ScaleData(seu_full, features = genes_present, assay = assay, verbose = FALSE)
    scaled   <- GetAssayData(seu_full, assay = assay, layer = "scale.data")
    
    # Project into sketch PCA space via matrix multiply
    shared_genes <- intersect(rownames(pca_rotation), rownames(scaled))
    emb <- t(scaled[shared_genes, ]) %*% pca_rotation[shared_genes, dims]
    
    # kNN label transfer
    knn_res <- FNN::get.knnx(data = sketch_emb, query = emb, k = k)
    transferred <- apply(knn_res$nn.index, 1, function(idx){
      names(sort(table(sketch_labels[idx]), decreasing = TRUE))[1]
    })
    names(transferred) <- colnames(seu_full)
    
    # Only fill non-sketch cells
    fill_cells <- intersect(setdiff(colnames(seu_full), sketch_cells), names(all_labels))
    all_labels[fill_cells] <- transferred[fill_cells]
    
    rm(seu_full); gc()
  }
  object@meta.data[,names(refdata)] <- all_labels
  return(object)
}
##-- Other functions ----
paste_matrix <- function(...,sep = " ",collapse = NULL){
  n <- max(sapply(list(...),nrow))
  p <- max(sapply(list(...),ncol))
  
  matrix(paste(...,sep = sep,collapse = collapse),n,p)
}

##################
#### Read 10x ####
#.... a) Read in all data ####
samples <- list.files(datadir)
.filtered_dir <- function(d, s){
  d2 <- file.path(d, 'SAMPLEID/outs/per_sample_outs') %>%
    gsub("SAMPLEID", s, .)
  
  sapply(list.files(d2), function(di){
    file.path(d2, di, 'count', 'sample_filtered_feature_bc_matrix')
  })
}

seus <- lapply(samples, function(s){
  dirs <- .filtered_dir(datadir, s)
  seus2 <- lapply(setNames(seq_along(dirs), names(dirs)), function(dir_idx){
    dir_i <- dirs[dir_idx]
    print(dir_i)
    seu2 <- tryCatch({
      mtx <- Read10X(data.dir = dir_i, strip.suffix=TRUE)
      seu <- CreateSeuratObject(counts = mtx, project=s)
      seu$orig.ident <- names(dir_i)
      seu$sample <- s
      seu$batch <- ifelse(all(seu$sample == seu$orig.ident), s, paste0(seu$sample, "_", seu$orig.ident))
      # seu[["RNA"]]$`counts.Antibody Capture` <- NULL
      seu
    }, error=function(e){NULL})
    # mtx=Read10X_h5(file.path(basedir, paste0(s, ".h5")))
    # meta <- read.table(file.path(basedir, paste0(s, "_metadata.txt")), header=T)
    return(seu2)
  })
  return(seus2)
})
names(seus) <- samples
seus <- unlist(seus, recursive=F)
null_samples <- names(seus)[sapply(seus, is.null)]
seus <- Filter(Negate(is.null), unlist(seus, recursive = FALSE))

#... b) Preprocess each seurat object ####
regress_vars <- c('pct_counts_mt')
variable_features_n <- 3000
outf <- file.path(outdir, "objs", "seus_sct_norm.q2")

if(file.exists(outf)){
  seus_norm <- qs2::qs_read(file=outf)
} else {
  seus_norm <- lapply(seus, function(seu){
    seu[[regress_vars]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    
    seu <- SCTransform(seu,
                       assay = 'RNA',
                       new.assay.name = 'SCT',
                       variable.features.n =  variable_features_n,
                       vst.flavor='v2',
                       vars.to.regress = regress_vars)
    return(seu)
  })
  qs2::qs_save(seus_norm, file=outf)
}


#... c) Process the full SCT normalized dataset as reference ####
outf <- file.path(outdir, "objs", "seus_full_harmony_integrated.q2")

if(file.exists(outf)){
  seu_full_merged <- qs2::qs_read(file=outf)
} else {
  seu_full_merged <- merge(
    x    = seus_norm[[1]],
    y    = seus_norm[-1],
    merge.dr = FALSE   # don't try to merge per-sample reductions
  )
  
  seu_full_merged <- FindVariableFeatures(seu_full_merged, nfeatures = variable_features_n) %>%
    ScaleData(.) %>%
    RunPCA(., npcs=50, verbose=F)
  
  params <- list(sigma = 0.1, theta = 2, lambda = 1)
  seu_full_merged <- IntegrateLayers(
    object               = seu_full_merged,
    method               = HarmonyIntegration,
    normalization.method = "SCT",
    orig.reduction       = "pca",
    new.reduction        = "integrated.harmony",
    npcs                 = 50L,
    verbose              = TRUE,
    sigma                = params$sigma,
    theta                = params$theta,
    lambda               = params$lambda,
    group.by             = "batch"
  )
  
  
  seu_full_merged <- FindVariableFeatures(seu_full_merged, assay='SCT') %>%
    RunPCA(., assay = "SCT") %>%
    FindNeighbors(., dims = 1:30, reduction = "integrated.harmony") %>%
    FindClusters(., resolution = 0.5) %>%
    RunUMAP(., dims = 1:30, reduction = "integrated.harmony")
  
  seu_full_merged <- PrepSCTFindMarkers(seu_full_merged)
  markers <- FindAllMarkers(seu_full_merged, logfc.threshold = 0)
  markers <- markers %>%
    mutate(biotype=gm$SYMBOL$gene_biotype[gene]) %>%
    scCustomize::Add_Pct_Diff() 
  markers2 <- markers %>% 
    dplyr::filter(!grepl("^ENSG", gene),
                  biotype=='protein_coding',
                  # pct_diff >= pctdiff,
                  p_val_adj <= 0.05) %>%
    arrange(p_val_adj, desc(round(pct_diff,1)), desc(avg_log2FC))
  
  gptcelltype(markers2, model = 'gpt-4', topgenenumber=40) %>% 
    paste0("Assume that these markers are differentiating between immune cells in an PDAC PBMC population \n", ., 
           "\nPlease put it in a c('cluster' = 'cell-type',...) format. one entry per line\n\n") %>%
    cat
  annomap <- c(
    '0'  = 'Naive/Central Memory CD4+ T cells',
    '1'  = 'Regulatory T cells (Tregs)',
    '2'  = 'Effector Memory CD8+ T cells',
    '3'  = 'Tissue-Resident Memory CD8+ T cells (Trm)',
    '4'  = 'NK cells',
    '5'  = 'Classical Monocytes',
    '6'  = 'Non-Classical/Patrolling Monocytes',
    '7'  = 'T cells (mitochondrial high/low quality)',
    '8'  = 'Classical Dendritic Cells (cDC2)',
    '9'  = 'NKT / Innate-like T cells',
    '10' = 'Neutrophils',
    '11' = 'Proliferating cells',
    '12' = 'Mature/Cytotoxic NK cells',
    '13' = 'Regulatory T cells (Tregs)',
    '14' = 'B cells',
    '15' = 'Non-Classical Monocytes / Monocyte-derived cells',
    '16' = 'Mast cells',
    '17' = 'Epithelial cells / Tumour cells (mitochondrial high)',
    '18' = 'Epithelial / Tumour cells (PDAC)',
    '19' = 'Plasma cells',
    '20' = 'Plasmacytoid Dendritic Cells (pDC)'
  )
  seu_full_merged@meta.data$ai_anno <- annomap[as.character(seu_full_merged@meta.data$seurat_clusters)]
  
  qs2::qs_save(seu_full_merged, file=outf)
}


#... d) Sketch each normalized seurat object ####
get_n_sketch <- function(seu, f) round(ncol(seu) * f)
use.frac= FALSE
n_sketch <- 3000
sketch_range <- c(100, 500, 1000, 2000, 3000, 4000, 5000)
fraction <- 0.1

for(n_sketch in sketch_range){
  outf <- file.path(outdir, "objs", paste0("seus_sketch.", n_sketch, ".q2"))
  
  if(file.exists(outf)){
    seus_sketch <- qs2::qs_read(file=outf)
  } else {
    seus_sketch <- lapply(seq_along(seus_norm), function(i){
      seu <- seus_norm[[i]]
      sid <- names(seus_norm)[i]      
      
      # Save full SCT object BEFORE stripping
      # Only need to do this once, not per sketch_range iteration
      full_outf <- file.path(outdir, "objs", "full", paste0(sid, "_full_SCT.q2"))
      if(!file.exists(full_outf)){
        dir.create(dirname(full_outf), showWarnings = FALSE)
        qs2::qs_save(seu, file = full_outf)
      }
      
      DefaultAssay(seu) <- "SCT"   # ensure SCT is active for SketchData
      n_sketch <- ifelse(use.frac, get_n_sketch(seu, fraction), n_sketch)
      seu <- SketchData(
        object = seu,
        ncells = n_sketch,
        method = "LeverageScore",
        sketched.assay = "sketch"
      )
      DefaultAssay(seu) <- "sketch"
      seu[["RNA"]] <- NULL  # drop raw if memory is a concern
      # seu[["SCT"]] <- NULL  # drop raw if memory is a concern
      return(seu)
    })
    names(seus_sketch) <- names(seus_norm)
    qs2::qs_save(seus_sketch, file=outf)
  }
  
}


#--- e) Integration of sketch ####
for(n_sketch in sketch_range){
  infile <- file.path(outdir, "objs", paste0("seus_sketch.", n_sketch, ".q2"))
  outf <- file.path(outdir, "objs", paste0("nsample_", n_samples), 
                    paste0("seus_sketch_integrated.", n_sketch, ".q2"))
  dir.create(dirname(outf), recursive = T, showWarnings = F)
  
  if(file.exists(outf)){
    message(paste0("Reading in: ", n_sketch))
    seu_merged <- qs2::qs_read(file=outf)
  } else {
    seus_sketch <- qs2::qs_read(infile)
    seus_sketch <- seus_sketch[sample(names(seus_sketch), size = n_samples)]
    seu_merged <- merge(
      x    = seus_sketch[[1]],
      y    = seus_sketch[-1],
      merge.dr = FALSE   # don't try to merge per-sample reductions
    )
    
    
    DefaultAssay(seu_merged) <- "sketch"
    seu_merged <- SCTransform(seu_merged, assay='sketch', vst.flavor = "v2", verbose = FALSE)
    
    DefaultAssay(seu_merged) <- "sketch"
    seu_merged <- FindVariableFeatures(seu_merged, nfeatures = variable_features_n) %>%
      ScaleData(.) %>%
      RunPCA(., npcs=50, verbose=F)
    
    params <- list(sigma = 0.1, theta = 2, lambda = 1)
    seu_merged <- IntegrateLayers(
      object               = seu_merged,
      method               = HarmonyIntegration,
      normalization.method = "SCT",
      orig.reduction       = "pca",
      new.reduction        = "integrated.harmony",
      npcs                 = 50L,
      verbose              = TRUE,
      sigma                = params$sigma,
      theta                = params$theta,
      lambda               = params$lambda,
      group.by             = "batch"
    )
    
    
    seu_merged <- FindVariableFeatures(seu_merged, assay='SCT') %>%
      RunPCA(., assay = "SCT") %>%
      FindNeighbors(., dims = 1:30, reduction = "integrated.harmony") %>%
      FindClusters(., resolution = 0.5) %>%
      RunUMAP(., dims = 1:30, reduction = "integrated.harmony")
  
    
    #--- e2) Projection back to full data ####
    # seu_merged <- ProjectIntegration(
    #   object      = seu_merged,
    #   sketched.assay = "sketch",
    #   assay       = "SCT",
    #   reduction   = "integrated.harmony"
    # )
    
    seu_merged <- ProjectData(
      object      = seu_merged,
      assay       = "SCT",
      sketched.assay = "sketch",
      sketched.reduction = "pca",
      full.reduction     = "full.pca",
      dims               = 1:30,
      refdata    = list(cluster = "seurat_clusters")
    )
    
    # Project labels to full cells using saved full objects
    if(any(is.na(seu_merged$cluster))){
      message("Fixing cluster projections")
      full_paths <- setNames(file.path(outdir, "objs", "full", paste0(names(seus_sketch), "_full_SCT.q2")),
                             names(seus_sketch))
      seu_merged <- ProjectDataSCT(
        object      = seu_merged,
        full_paths = full_paths,
        assay       = "SCT",
        sketched.assay = "sketch",
        sketched.reduction = "pca",
        dims               = 1:30,
        k          = 5,
        refdata    = list(cluster = "seurat_clusters")
      )
    }
    qs2::qs_save(seu_merged, file=outf)
  }
}


#--- f) comparison across sketches ####
sketch_f <- list.files(file.path(outdir, "objs", paste0("nsample_", n_samples)), pattern="sketch.*integrated")
sketch_sizes <- gsub("^.*\\.([0-9]*)\\..*", "\\1", sketch_f)
sketch_f <- setNames(file.path(outdir, "objs", paste0("nsample_", n_samples), sketch_f), sketch_sizes)
seu_sketches <- lapply(sketch_f, qs2::qs_read)

fullf <- file.path(outdir, "objs", "seus_full_harmony_integrated.q2")
seu_full_merged <- qs2::qs_read(file=fullf)

.plotit <- function(seu, group=c('seurat_clusters', 'batch', 'ai_anno')){
  red <- 'umap'
  gg <- lapply(group, function(grp){
    tryCatch({
      scCustomize::DimPlot_scCustom(seu, reduction=red, group.by = grp, 
                                    label = T, label.size = 3) +
        ggtitle(grp)
    }, error=function(e){NULL})
  }) %>%
    cowplot::plot_grid(plotlist = ., nrow=1)
  return(gg)
}

# pdf("~/xfer/sketch.pdf", width=23)
# .plotit(seu_full_merged)
# .plotit(seu_merged_withSCT)
# .plotit(seu_merged_withoutSCT)
# dev.off()




# seu_sketches <- list("w_sct3000"=seu_merged_withSCT, 
#                      "wo_sct3000"=seu_merged_withoutSCT)
for(sketchid in names(seu_sketches)){
  seu_merged <- seu_sketches[[sketchid]]
  shared_cells <- intersect(colnames(seu_full_merged), colnames(seu_merged))
  col_sketchid <- paste0("sketch_", sketchid)
  seu_full_merged@meta.data[,col_sketchid] <- '-1'
  seu_full_merged@meta.data[shared_cells,col_sketchid] <- seu_merged@meta.data[shared_cells,'cluster']
  na.idx <- is.na(seu_full_merged@meta.data[shared_cells,col_sketchid])
  seu_full_merged@meta.data[which(na.idx),col_sketchid] <- -1
}

pdf("~/xfer/sketch_full.pdf", width=23)
.plotit(seu_full_merged, group=c('ai_anno', 'seurat_clusters'))
for(sketchid in sort(as.integer(names(seu_sketches)))){
  .plotit(seu_full_merged, group=c('ai_anno', paste0('sketch_', sketchid))) %>%
    plot
}
dev.off()


#--- g) Single score metric: ARI and NMI similarity scoring ####
library(aricode)   # ARI + NMI
global_metrics <- lapply(names(seu_sketches), function(sketchid) {
  seu_sk     <- seu_sketches[[sketchid]]
  shared     <- intersect(colnames(seu_full_merged), colnames(seu_sk))
  if (!length(shared)) return(NULL)
  
  full_labels   <- seu_full_merged@meta.data[shared, "ai_anno"]
  sketch_labels <- seu_sk@meta.data[shared, "cluster"]
  
  # Drop NAs
  keep          <- !is.na(full_labels) & !is.na(sketch_labels)
  full_labels   <- full_labels[keep]
  sketch_labels <- sketch_labels[keep]
  
  data.frame(
    sketchid = sketchid,
    n_cells  = sum(keep),
    ARI      = aricode::ARI(full_labels,  sketch_labels),
    NMI      = aricode::NMI(full_labels,  sketch_labels),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()


# Quick plot of ARI/NMI across sketches
pdf("~/xfer/ari.pdf")
global_metrics %>%
  pivot_longer(c(ARI, NMI), names_to = "metric", values_to = "value") %>%
  ggplot(aes(x = sketchid, y = value, colour = metric, group = metric)) +
  geom_line() +
  geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Sketch vs full annotation similarity",
       x = "Sketch ID", y = "Score", colour = NULL) +
  theme_minimal(base_size = 12)
dev.off()

#--- h) Alluvial between sketch and refernece ####
library(ggalluvial)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(cowplot)


.build_comparison_df <- function(seu_full, seu_sketch, sketchid) {
  shared_cells <- intersect(colnames(seu_full), colnames(seu_sketch))
  if (!length(shared_cells)) {
    warning("No shared cells for sketch: ", sketchid)
    return(NULL)
  }
  
  data.frame(
    cell          = shared_cells,
    full_anno     = as.character(seu_full@meta.data[shared_cells, "ai_anno"]),
    full_cluster  = as.character(seu_full@meta.data[shared_cells, "seurat_clusters"]),
    sketch_cluster = as.character(seu_sketch@meta.data[shared_cells, "cluster"]),
    sketchid      = sketchid,
    stringsAsFactors = FALSE
  )
}

.plot_sankey <- function(df, title = "") {
  # Summarise cell counts per sketch_cluster → full_anno transition
  alluv <- df %>%
    count(sketch_cluster, full_anno) %>%
    # Order axes by frequency for readability
    mutate(
      sketch_cluster = factor(sketch_cluster,
                              levels = names(sort(table(df$sketch_cluster), decreasing = TRUE))),
      full_anno = factor(full_anno,
                         levels = names(sort(table(df$full_anno),      decreasing = TRUE)))
    )
  
  ggplot(alluv,
         aes(axis1 = sketch_cluster, axis2 = full_anno, y = n)) +
    geom_alluvium(aes(fill = sketch_cluster), alpha = 0.7, width = 0.3) +
    geom_stratum(width = 0.3, fill = "grey90", color = "grey40") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(
      limits = c("Sketch Cluster", "Full Annotation"),
      expand = c(0.1, 0.1)
    ) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = title, y = "Cells", fill = "Sketch cluster") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          panel.grid      = element_blank())
}

.plot_confusion <- function(df, title = "") {
  mat <- df %>%
    count(sketch_cluster, full_anno) %>%
    group_by(sketch_cluster) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    dplyr::select(sketch_cluster, full_anno, prop) %>%
    pivot_wider(names_from = full_anno,
                values_from = prop,
                values_fill = 0) %>%
    tibble::column_to_rownames("sketch_cluster") %>%
    as.matrix()
  
  # Sort rows/cols by hierarchical clustering for cleaner display
  pheatmap(
    mat,
    color            = colorRampPalette(c("white", "#2166ac"))(100),
    border_color     = "grey90",
    cluster_rows     = TRUE,
    cluster_cols     = TRUE,
    display_numbers  = TRUE,
    number_format    = "%.2f",
    fontsize_number  = 7,
    fontsize_row     = 9,
    fontsize_col     = 9,
    main             = title,
    angle_col        = 45,
    silent           = TRUE   # return grob without printing
  )$gtable
}

.plot_stacked_bar <- function(df, title = "") {
  df %>%
    count(sketch_cluster, full_anno) %>%
    group_by(sketch_cluster) %>%
    mutate(prop = n / sum(n),
           sketch_cluster = factor(sketch_cluster,
                                   levels = names(sort(table(df$sketch_cluster), decreasing = TRUE)))) %>%
    ggplot(aes(x = sketch_cluster, y = prop, fill = full_anno)) +
    geom_col(width = 0.8, color = "white", linewidth = 0.2) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = title,
         x     = "Sketch cluster",
         y     = "% of cells",
         fill  = "Full annotation") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1),
          legend.key.size = unit(0.4, "cm"))
}

pdf("~/xfer/sketch_vs_full_comparison.pdf", width = 20, height = 14)
for (sketchid in sort(as.integer(names(seu_sketches)))) {
  print(sketchid)
  sketchid <- as.character(sketchid)
  seu_sk   <- seu_sketches[[sketchid]]
  
  df <- .build_comparison_df(seu_full_merged, seu_sk, sketchid)
  if (is.null(df)) next
  
  n_shared <- nrow(df)
  title_str <- sprintf("Sketch %s  |  %s shared cells", sketchid,
                       formatC(n_shared, big.mark = ",", format = "d"))
  
  p_sankey  <- .plot_sankey(df,      paste("Sankey —",      title_str))
  p_bar     <- .plot_stacked_bar(df, paste("Stacked bar —", title_str))
  p_heatmap <- .plot_confusion(df,   paste("Confusion —",   title_str))
  
  # Combine: sankey left, bar top-right, heatmap bottom-right
  right_col <- plot_grid(p_bar, patchwork::wrap_elements(p_heatmap), ncol = 1)
  page      <- plot_grid(p_sankey, right_col, ncol = 2, rel_widths = c(1.2, 1))
  print(page)
}

dev.off()
message("Saved: sketch_vs_full_comparison.pdf")

#--- h) Alluvial across all sketches ####
library(ggalluvial)
library(dplyr)
library(tidyr)

# Start with full annotation for all cells
df_all <- data.frame(
  cell      = colnames(seu_full_merged),
  full_anno = as.character(seu_full_merged@meta.data[, "ai_anno"]),
  stringsAsFactors = FALSE
)

# Join each sketch's cluster assignment
for (sketchid in sort(as.integer(names(seu_sketches)))) {
  sketchid  <- as.character(sketchid)
  seu_sk    <- seu_sketches[[sketchid]]
  col_name  <- paste0("sketch_", sketchid)
  
  shared    <- intersect(colnames(seu_full_merged), colnames(seu_sk))
  
  df_all[[col_name]] <- NA_character_
  df_all[df_all$cell %in% shared, col_name] <-
    seu_sk@meta.data[shared, "cluster"]
}

# Only keep cells that appear in ALL sketches
# (or relax to: cells in at least one sketch)
sketch_cols <- grep("^sketch_", colnames(df_all), value = TRUE)

df_shared <- df_all %>%
  dplyr::filter(if_all(all_of(sketch_cols), ~ !is.na(.))) %>%
  # Prefix sketch cluster labels to avoid axis ambiguity
  # e.g. sketch_1 cluster "3" becomes "sk1_3"
  mutate(across(all_of(sketch_cols), ~ {
    sk <- sub("sketch_", "sk", cur_column())
    paste0(sk, "_", .x)
  }))

df_counts <- df_shared %>%
  count(across(all_of(c("full_anno", sketch_cols)))) %>%
  # Convert to long/lode form for ggalluvial
  to_lodes_form(axes = c("full_anno", sketch_cols),
                key   = "axis",
                value = "stratum",
                id    = "alluvium") %>%
  mutate(axis = factor(axis, levels = c("full_anno", sketch_cols)))

n_anno <- length(unique(df_shared$full_anno))
cols   <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_anno),
  sort(unique(df_shared$full_anno))
)

# Map alluvium colour to full_anno
# (need to join full_anno back onto the lodes form)
full_anno_map <- df_counts %>%
  dplyr::filter(axis == "full_anno") %>%
  dplyr::select(alluvium, full_anno = stratum)

df_counts <- df_counts %>%
  left_join(full_anno_map, by = "alluvium")

p <- ggplot(df_counts,
            aes(x      = axis,
                y      = n,
                stratum = stratum,
                alluvium = alluvium,
                fill   = full_anno,
                label  = stratum)) +
  geom_flow(alpha = 0.4, width = 0.3, aes(fill = full_anno)) +
  geom_stratum(width = 0.3, color = "grey30", linewidth = 0.2) +
  geom_text(stat = "stratum", size = 2.5, min.y = 50) +  # min.y suppresses labels on tiny strata
  scale_fill_manual(values = cols, na.value = "grey80") +
  scale_x_discrete(
    labels = c("full_anno" = "Full\nAnnotation",
               setNames(paste0("Sketch ", sub("sketch_", "", sketch_cols)), sketch_cols))
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Full annotation vs sketch cluster assignments",
    x     = NULL,
    y     = "Cells",
    fill  = "Full annotation"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid    = element_blank(),
    axis.text.x   = element_text(size = 11, face = "bold"),
    legend.key.size = unit(0.4, "cm")
  )

# Save — width scales with number of sketches
n_sketches <- length(sketch_cols)
ggsave("~/xfer/sankey_full_vs_all_sketches.pdf", p,
       width  = 4 + (n_sketches * 2.5),
       height = 12,
       limitsize = FALSE)


#####################
#### Read Xenium ####
sketch_range <- c(1000, 300, 5000, 10000, 25000, 50000)
variable_features_n <- 3000

samples <- sample(list.files(xeniumdir), size=n_samples)
.checkCoords <- function(obj, fov='fov'){
  print("Centroids")
  print(head(GetTissueCoordinates(obj[[fov]][["centroids"]])))    # Spatial coordinates of the centroid for each cell
  print("Segmentation")
  print(head(GetTissueCoordinates(obj[[fov]][["segmentation"]]))) # Spatial coordinates of polygon describing each cell
  print("Molecules")
  print(head(FetchData(obj[[fov]][["molecules"]], vars = "CD4"))) # Spatial coordinates of each molecule position
}

#.... a) Create xenium obj ####
outf <- file.path(outdir, "objs", "xen_raw.q2")

xens <- lapply(samples, function(s){
  # s = samples[2]
  xenium.obj <- LoadXenium(file.path(xeniumdir, s), 
                           fov = "fov", segmentations = "cell") # requires cell_feature_matrix, transcripts.parquet
  .checkCoords(xenium.obj)
  # remove cells with 0 counts
  xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)
  return(xenium.obj)
})
names(xens) <- samples
for(i in seq_along(samples)){
  xens[[i]]$orig.ident <- gsub("^.*HPB", "HPB", samples[i]) %>%
    gsub("__.*", "", .)
}
qs2::qs_save(xens, file=outf)

#... b) Preprocess each xenium object ####
infile <- file.path(outdir, "objs", "xen_raw.q2")
outfile <- file.path(outdir, "objs", "xen_sct.q2")

# https://satijalab.org/seurat/articles/spatial_vignette.html#working-with-multiple-slices-in-seurat
# Working with multiple slices enables joint dimensional reduction
# and clustering on the underlying RNA expression data.
#
# The data can be jointly visualized in a single UMAP plot.
dir.create(file.path(outdir, "qs2", "sct"), showWarnings = F)
if(!exists("xens")) xens <- qs2::qs_read(infile)
xens_sct <- lapply(names(xens), function(xenid){
  xenium.obj <- xens[[xenid]]
  sampleid <- unique(xenium.obj$orig.ident)
  xenium.obj <- SCTransform(xenium.obj, assay = "Xenium", vst.flavor = "v2") #%>%
  xenium.obj
})
names(xens_sct) <- names(xens)

qs2::qs_save(xens_sct, file=outfile)
xens <- xens_sct
rm(xens_sct); gc()

#... c) Sketch each normalized seurat object ####
infile <- file.path(outdir, "objs", "xen_sct.q2")
if(!exists("xens")) xens <- qs2::qs_read(infile)


for(n_sketch in sketch_range){
  outf <- file.path(outdir, "objs", paste0("xens_sketch.", n_sketch, ".q2"))
  
  if(file.exists(outf)){
    xens_sketch <- qs2::qs_read(file=outf)
  } else {
    xens_sketch <- lapply(seq_along(xens), function(i){
      xen <- xens[[i]]
      sid <- names(xens)[i]      
      
      # Save full SCT object BEFORE stripping
      # Only need to do this once, not per sketch_range iteration
      full_outf <- file.path(outdir, "objs", "full", paste0(sid, "_xen_full_SCT.q2"))
      if(!file.exists(full_outf)){
        dir.create(dirname(full_outf), showWarnings = FALSE)
        qs2::qs_save(xen, file = full_outf)
      }
      
      DefaultAssay(xen) <- "SCT"   # ensure SCT is active for SketchData
      xen <- SketchData(
        object = xen,
        ncells = n_sketch,
        method = "LeverageScore",
        sketched.assay = "sketch"
      )
      DefaultAssay(xen) <- "sketch"
      xen[["RNA"]] <- NULL  # drop raw if memory is a concern
      # xen[["SCT"]] <- NULL  # drop raw if memory is a concern
      return(xen)
    })
    names(xens_sketch) <- names(xens)
    qs2::qs_save(xens_sketch, file=outf)
  }
  
}



#----------------------------
#### Merging SCT objects ####
if(!exists("xens_sct")) xens_sct <- qs2::qs_read(file=file.path(outdir, "qs2", "xens_sct.qs2"))
xen_merge <- merge(xens_sct[[1]], xens_sct[-1])
DefaultAssay(xen_merge) <- "SCT"
merge_varfeats <- sapply(xens_sct, VariableFeatures)  %>%
  unlist %>% as.character %>% sort %>% unique
rm(xens_sct); gc()
VariableFeatures(xen_merge) <- merge_varfeats

xen_merge <- RunPCA(xen_merge, verbose = FALSE)
options(future.globals.maxSize = 2000 * 1024^2) # Sets the limit to 2 GB
xen_merge <- FindNeighbors(xen_merge, dims = 1:30)
xen_merge <- FindClusters(xen_merge, resolution=0.4, verbose = FALSE)
xen_merge <- RunUMAP(xen_merge, dims = 1:30)
qs2::qs_save(xen_merge, file=file.path(outdir, "qs2", "xen_merge.qs2"))
xen_merge <- qs2::qs_read(file=file.path(outdir, "qs2", "xen_merge.qs2"))

pdf("~/xfer/xen_dim.pdf", width=15)
scCustomize::DimPlot_scCustom(xen_merge, reduction='umap',
                              group.by=c('seurat_clusters', 'orig.ident'))
scCustomize::DimPlot_scCustom(xen_merge, reduction='umap',
                              group.by=c('monaco.main.cluster', 'monaco.fine.cluster'))
scCustomize::DimPlot_scCustom(xen_merge, reduction='umap',
                              group.by=c('gpt_anno'))
dev.off()
set.seed(1234)

x <- GetAssayData(xen_merge, layers='SCT')

X <- split(as.data.frame(x[grep("TCR", Features(xen_merge), value=T),]), f=xen_merge$gpt_anno)

#---------------------------