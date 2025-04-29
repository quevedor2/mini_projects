# renv::load("/cluster/home/quever/downloads/renvs/")

library(tidyverse)
library(ComplexHeatmap)
library(ggrastr)
library(DESeq2)
library(ggplot2)
library(reshape2)
library(umap)
library(cowplot)
library(ggrepel)
library(GSVA)
library(org.Hs.eg.db)
library(icellnet)
library("org.Hs.eg.db")
library("hgu133plus2.db")
library(gridExtra)
library(jetset)
library(RColorBrewer)
library(clusterProfiler)


dataset <- 'peritoneal' # 'bonemarrow' 'peritoneal'
colrs <- brewer.pal.info[brewer.pal.info$colorblind == TRUE, ]
col_vec = unlist(mapply(brewer.pal, colrs$maxcolors, rownames(colrs)))
species <- 'Mus musculus'
pdir <- file.path('/cluster/projects/mcgahalab/data/mcgahalab/xin_GCN2', dataset, 'results')

gprofiler_dir <- '/cluster/projects/mcgahalab/ref/gprofiler'
gprofiler_f <- file.path(gprofiler_dir, 'gprofiler_full_mmusculus.ENSG.gmt')
# 
barcodes_f="~/git/mini_projects/ref/barcodes.tsv"

#### Functions ####
# code cloned from https://github.com/quevedor2/mini_projects
source("~/git/mini_projects/mini_functions/geneMap.R")
source("~/git/mini_projects/mini_functions/makeLoupe.R")
gm <- geneMap(species=species)

###################
#### Functions ####
ssGseaFun <- function(msig_ds, lfc_v, ss_method='ssgsea'){
  require(GSVA)
  ssgsea <- tryCatch({
    sig_ens_gs <- split(setNames(msig_ds$entrez_gene, msig_ds$entrez_gene), 
                        f=msig_ds$gs_name)
    gsva(lfc_v, sig_ens_gs, verbose=FALSE, method=ss_method)
  }, error=function(e){NULL})
  return(ssgsea)
}

getDEGedgeRwt <- function(cts, meta, group){
  idx <- match(rownames(meta), colnames(cts))
  if(any(is.na(idx))){
    na_idx <- which(is.na(idx))
    # print(paste0("NA idx: ", paste(na_idx, collapse=",")))
    idx <- idx[-na_idx]
    meta <- meta[-na_idx,]
  }
  
  ## Get EdgeR Results
  se <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(cts[,idx])), 
                             colData = meta)
  
  
  edgey <- tryCatch({
    edgey <- DGEList(counts=cts[,idx],
                     samples=rownames(meta),
                     group=meta[,group])
    # design <- with(meta, model.matrix(as.formula(formula)))
    design <- model.matrix(~ meta[,group])
    keep <- filterByExpr(edgey, design)
    edgey <- edgey[keep, , keep.lib.sizes=FALSE]
    edgey <- calcNormFactors(edgey, method = "TMM")
    edgey <- estimateDisp(edgey)
    edgey
  }, error=function(e){NULL})
  
  # Calculate CPM/TMM values
  edge_tmm <- cpm(edgey)
  edge_tmm_spl <- edge_tmm %>% t %>% as.data.frame %>%
    split(., meta[,group]) 
  
  # Differential testing
  er_dat <- exactTest(edgey)
  et_res <- er_dat$table  %>% 
    as.data.frame %>% 
    tibble::rownames_to_column("ensemble") %>% 
    mutate(padj=p.adjust(PValue, method='BH')) %>%
    mutate(sig=ifelse(padj < 0.05, "er_sig", "er_ns"))
  
  glvl <- levels(factor(meta[,group]))
  wilcox_res <- sapply(seq_along(edge_tmm_spl[[1]]), function(idx){  
    wt <- wilcox.test(edge_tmm_spl[[glvl[1]]][,idx], edge_tmm_spl[[glvl[2]]][,idx])
    fc <- mean(edge_tmm_spl[[glvl[1]]][,idx],na.rm=T) /  mean(edge_tmm_spl[[glvl[2]]][,idx], na.rm=T)
    c('W'=wt$statistic, 'pval'=wt$p.value, "FC"=fc, 'Log2FC'=log2(fc+1),
      'ensemble'=colnames(edge_tmm_spl[[1]][idx]))
  }) %>% 
    t %>% as.data.frame %>%
    mutate(padj=p.adjust(pval, method='BH'),
           FC=as.numeric(FC),
           Log2FC=as.numeric(Log2FC)) %>%
    mutate(sig=ifelse(padj < 0.05, "wt_sig", "wt_ns"))
  
  # pdf("~/xfer/dds2.pdf") 
  # # ggplot(dds_et, aes(x=log2FoldChange, y=logFC)) +
  # ggplot(dds_et, aes(x=log2FoldChange, y=logFC, color=sig, group=sig)) +
  #   geom_point()
  # ggplot(wt_et, aes(x=Log2FC, y=logFC, color=sig, group=sig)) +
  #   geom_point()
  # dev.off()
  return(list("edger"=et_res,  "wilcox"=wilcox_res))
}

rmVstBatchEffects <- function(vsd, dds, samples_meta, 
                              condition='condition', batchcol='batch'){
  mat <- assay(vsd)
  mm <- model.matrix(as.formula(paste0("~", condition)), colData(dds))
  mat <- limma::removeBatchEffect(mat, design=mm, 
                                  batch=samples_meta[colnames(dds),batchcol])
  assay(vsd) <- mat
  return(vsd)
}


##############
#### Main ####
dir.create(file.path(pdir, "manual", "objs"), showWarnings = F, recursive = T)
dir.create(file.path(pdir, "manual", "pca"), showWarnings = F, recursive = T)
dir.create(file.path(pdir, "manual", "loupe"), showWarnings = F, recursive = T)
dir.create(file.path(pdir, "manual", "deg"), showWarnings = F, recursive = T)
dir.create(file.path(pdir, "manual", "gsea"), showWarnings = F, recursive = T)
dir.create(file.path(pdir, "manual", "ssgsea"), showWarnings = F, recursive = T)

rm.samples <- c('XS_A', 'XS_9', 'XS_11')

deganaldir <- file.path(pdir, "manual", "differential_expression")
file <- 'all.tsv'
min_expr <- 3
expr_frac_samples <- 0.2

# Load bulk RNAseq data
data=read.table(file.path(pdir, "counts", file), header = T, check.names=FALSE, 
                stringsAsFactors = FALSE, na.strings = "")
genes <- gm$ENSEMBL$SYMBOL[data$gene]
genes[is.na(genes)] <- 'NA'

data <- data %>% 
  dplyr::select(-c(gene)) %>% 
  split(., f=genes) %>%
  sapply(., colMeans)  %>% 
  t  %>% as.data.frame %>%
  rename_with(., ~gsub("genesresults$", "", .)) 

if(length(rm.samples)>0){
  rm.samples <- sapply(rm.samples, grep, x=colnames(data), value=T)
  keep.samples <- setdiff(colnames(data), rm.samples)
  data <- data[,keep.samples]
}



# Remove genes where [expr_frac_sample] of the dataset has fewer than [min_expr] 
# reads linked to that gene
# e.g. Remove: LOC101929116 is expressed at higher than 3 counts in only 2/42 samples
low_expr_idx <- which(rowSums(data >= min_expr) >= 2 ) #(ncol(data) * expr_frac_samples))
data <- data[low_expr_idx,]



#### 1. Create metadata ####
# Create metadata data frame
metadata <- colnames(data)
if(dataset == 'bonemarrow'){
  metad <- strsplit(gsub("_S[0-9]*$", "", metadata), split="_") %>% 
    do.call(rbind, .) %>% 
    magrittr::set_colnames(., c('Marker', 'Lys', 'Treatment', 'sample_number')) %>%
    as.data.frame %>% 
    dplyr::mutate('sample'=as.character(metadata),
                  "condition"=paste0(Marker, "_", Lys, "_", Treatment))
} else if(dataset=='peritoneal'){
  metad <- strsplit(gsub("_S[0-9]*$", "", metadata), split="_") %>% 
    do.call(rbind, .) %>% 
    magrittr::set_colnames(., c('XS', 'Sample')) %>%
    as.data.frame %>% 
    dplyr::mutate('sample'=as.character(metadata),
                  'condition'=grepl("[0-9]+", Sample))
}


saveRDS(metad, file.path(pdir, "manual", "objs", "metad.rds"))

data <- as.matrix(data)
storage.mode(data) <- 'integer'
dds_all <- DESeqDataSetFromMatrix(countData=data,
                                  colData=metad,
                                  design=as.formula("~condition"))
saveRDS(dds_all, file=file.path(pdir, "manual", "objs", "deseq_all.rds"))


#### 2. ssGSEA ####
dds_all <- readRDS(file=file.path(pdir, "manual", "objs", "deseq_all.rds"))
metad <- readRDS(file=file.path(pdir, "manual", "objs", "metad.rds"))

# Read in geneset database
gmt <- GSA::GSA.read.gmt(gprofiler_f)
gprof_ds <-setNames(gmt$genesets, 
                    paste0(unlist(gmt$geneset.names), "_", unlist(gmt$geneset.descriptions)))
gprof_ds <- gprof_ds[grep("^CORUM", names(gprof_ds), invert=T)]


msig_ds <- lapply(names(gprof_ds), function(sublvl){
  data.frame("gs_name"=sublvl,
             "entrez_gene"=gprof_ds[[sublvl]])
}) %>% 
  do.call(rbind,.) %>% 
  mutate(classification =  gsub(":.*", "", gs_name))
msig_ds <- msig_ds %>%
  dplyr::mutate(entrez_gene = gm$ENSEMBL$SYMBOL[entrez_gene]) %>% 
  dplyr::filter(!classification %in% c('HP', 'HPA'),
                !is.na(entrez_gene))

# Obtain normalized counts
counts <- vst(dds_all, blind=T)
metad2 <- metad %>%
  tibble::column_to_rownames("sample")
tcnts <- rmVstBatchEffects(counts, dds_all, metad2, 
                           condition='condition') %>%
  assay %>% as.data.frame

# Run GSVA ssGSEA
sig_ens_gs <- split(setNames(msig_ds$entrez_gene, msig_ds$entrez_gene), 
                    f=msig_ds$gs_name)
gsvapar <- GSVA::ssgseaParam(as.matrix(tcnts), sig_ens_gs)
ssgsea_dat <- gsva(gsvapar)
write.table(ssgsea_dat, file=file.path(pdir, "manual", "ssgsea", "ssgsea.csv"),
            sep=",", col.names = T, row.names = T, quote = F)
saveRDS(ssgsea_dat, file=file.path(pdir, "manual", "ssgsea", "ssgsea.rds"))

#### 3. PCA reduction ####
metad <- readRDS(file.path(pdir, "manual", "objs", "metad.rds"))
max_pc <- 10
metad2 <- metad %>%
  tibble::column_to_rownames("sample")


#--- a) PCA based on GENE expression ----
dds_all <- readRDS(file=file.path(pdir, "manual", "objs", "deseq_all.rds"))

# Run the main PCA analysis on vst-counts
# obtain normalized counts
counts <- vst(dds_all, blind=T)
counts <- rmVstBatchEffects(counts, dds_all, metad2, 
                            condition='condition')

tcnts <- as.data.frame(t(assay(counts)))
pca <- prcomp(tcnts, scale=F)
pca_x <- as.data.frame(cbind(pca$x[,paste0("PC", c(1:max_pc))],
                             "condition"=as.character(counts$condition)))
if(max_pc > ncol(pca_x)) max_pc <- ncol(pca_x)
for(id in paste0("PC", c(1:max_pc))){
  pca_x[,id] <- as.numeric(as.character(pca_x[,id]))
}
saveRDS(pca_x, file=file.path(pdir, "manual", "pca", "pca_dat.all.rds"))

#--- b) PCA based on PATHWAY expression ----
ssgsea_dat <- readRDS(file=file.path(pdir, "manual", "ssgsea", "ssgsea.rds"))
idx <- grep("^(GO|REAC|WP)", rownames(ssgsea_dat))
ssgsea_dat <- ssgsea_dat[idx,keep.samples]

# Run the main PCA analysis on vst-counts
# obtain normalized counts
texpr <- as.data.frame(t(ssgsea_dat))
pca <- prcomp(texpr, scale=F)
pca_x <- as.data.frame(cbind(pca$x[,paste0("PC", c(1:max_pc))]))
if(max_pc > ncol(pca_x)) max_pc <- ncol(pca_x)
for(id in paste0("PC", c(1:max_pc))){
  pca_x[,id] <- as.numeric(as.character(pca_x[,id]))
}
saveRDS(pca_x, file=file.path(pdir, "manual", "pca", "pca_dat.ssgsea.rds"))

#### 4. Loupe visualization ####
dds_all <- readRDS(file=file.path(pdir, "manual", "objs", "deseq_all.rds"))
pca_x <- readRDS(file=file.path(pdir, "manual", "pca", "pca_dat.all.rds"))
pca_pathway <- readRDS(file=file.path(pdir, "manual", "pca", "pca_dat.ssgsea.rds"))
pca_pathway
ssgsea_dat <- readRDS(file=file.path(pdir, "manual", "ssgsea", "ssgsea.rds"))
idx <- grep("^(GO|REAC|WP)", rownames(ssgsea_dat))
ssgsea_dat <- ssgsea_dat[idx,]
metad <- readRDS(file=file.path(pdir, "manual", "objs", "metad.rds"))
metad2 <- metad %>% 
  tibble::column_to_rownames("sample")

counts <- vst(dds_all, blind=T)
counts <- rmVstBatchEffects(counts, dds_all, metad2, 
                            condition='condition')
makeLoupe(mat=assay(counts), meta=metad2, 
          projections=list("gene"=dplyr::select(pca_x, -condition),"ssgsea"=pca_pathway),
          output_dir=file.path(pdir, "manual", "loupe"),
          output_name="pca_all",
          barcodes_f=barcodes_f)
makeLoupe(mat=ssgsea_dat, meta=metad2, 
          projections=list("gene"=dplyr::select(pca_x, -condition),"ssgsea"=pca_pathway),
          output_dir=file.path(pdir, "manual", "loupe"),
          output_name="pca_pathway",
          barcodes_f=barcodes_f)
file.copy(file.path(pdir, "manual", "loupe", "pca_all.cloupe"), to="~/xfer", overwrite = T)
file.copy(file.path(pdir, "manual", "loupe", "pca_pathway.cloupe"), to="~/xfer", overwrite = T)

#### 5. Differential Analysis ####
dds_all <- readRDS(file=file.path(pdir, "manual", "objs", "deseq_all.rds"))
deseqobj <- dds_all
deseqobj$condition <- factor(sample(letters[1:4], size=ncol(deseqobj), replace=T))

pairwise_comparisons(deseqobj, save_dir='/path/to/outdir')
pairwise_comparisons <- function(dds, save_dir=NULL){
  conditions <- levels(dds$condition)
  ## Run all pairwise tests
  all_comps <- lapply(conditions, function(condition_i){
    # Iterate through each sample, setting one condition as a reference and all
    # other conditions as the comparative
    condition_not_i <- setdiff(conditions, condition_i)
    coef_names <- paste("condition", condition_not_i, "vs", condition_i, sep = "_")
    
    dds$condition <- relevel(dds$condition, ref = condition_i)
    
    # Calculate differnetial expression using DESeq, then adjust LFC using apeglm
    dds <- DESeq(dds)
    # res <- results(dds)
    resl <- lapply(coef_names, function(coef_name_j){
      ape_res <- lfcShrink(dds, coef=coef_name_j, type="apeglm")
      ape_res$comparison <- coef_name_j
      return(ape_res)
    })
    
    # Identify the comparisons made
    names(resl) <- gsub("condition_", "", coef_names)
    return(resl)
  })
  
  # Return all comparisons as a single list rather than a list of lists
  all_comps <- unlist(all_comps, recursive=F)
  if(!is.null(save_dir)){
    message(paste0("Saving each DEG comparison to: ", save_dir))
    for(resid in names(all_comps)){
      message(paste0("....writing: ", resid))
      write.table(all_comps[[resid]], 
                  file=file.path(save_dir, paste0("deg.", resid, ".csv")),
                  sep=",", col.names = T, row.names = T, quote = F)
    }
  }
  
  return(all_comps)
}
#### 5. Differential analysis ####
# Comparing H1048-CM samples to H841-CM samples in HD143 and HD151
idx <- intersect(grep("HD143|HD151", metad$HDsample),
                 grep("H1048|H841", metad$Cellline))
dds <- DESeqDataSetFromMatrix(countData=data[,idx],
                              colData=metad[idx,],
                              design=as.formula("~Cellline"))
# Setting H841-CM treated as the baseline (i.e. a +ve LFC will indicate that H1048 is larger than H841 samples)
dds$Cellline <- relevel(dds$Cellline, ref = "H841")
levels(dds$Cellline) # ensure proper levelling

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- lfcShrink(dds, coef="Cellline_H1048_vs_H841", type="apeglm") %>% # Cellline_H1048_vs_Donor = resultsNames(dds)[2]
  as.data.frame %>%
  arrange(padj)
saveRDS(res, file=file.path(pdir, "manual", "deg", "Cellline_H1048_vs_H841.rds"))

#### 6. GSEA ####
res <- readRDS(file=file.path(pdir, "manual", "deg", "Cellline_H1048_vs_H841.rds"))
res <- res %>% 
  tibble::rownames_to_column("symbol") %>%
  mutate(entrez = gm$SYMBOL$ENSEMBL[symbol])

# Read in geneset database
gmt <- GSA::GSA.read.gmt(gprofiler_f)
gprof_ds <-setNames(gmt$genesets, 
                    paste0(unlist(gmt$geneset.names), "_", unlist(gmt$geneset.descriptions)))
gprof_ds <- gprof_ds[grep("^CORUM", names(gprof_ds), invert=T)]


msig_ds <- lapply(names(gprof_ds), function(sublvl){
  data.frame("gs_name"=sublvl,
             "entrez_gene"=gprof_ds[[sublvl]])
}) %>% 
  do.call(rbind,.) %>% 
  mutate(classification =  gsub(":.*", "", gs_name))

# run gsea on the geneset
lfc_v <- setNames(res$log2FoldChange,
                  res$entrez)
gsea <- clusterProfiler::GSEA(sort(na.omit(lfc_v), decreasing = T), 
             TERM2GENE = msig_ds, pvalueCutoff = 1)
gsea_short <- as.data.frame(gsea) %>% 
  dplyr::select(-c(ID, Description, core_enrichment)) %>%
  arrange(p.adjust)
saveRDS(gsea, file=file.path(pdir, "manual", "gsea", "Cellline_H1048_vs_H841.rds"))

