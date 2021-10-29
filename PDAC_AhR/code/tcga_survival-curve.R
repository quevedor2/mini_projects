###########################
#### Read in TCGA Data ####
library(reshape2)
library(ggplot2)
library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)
library(gridExtra)
library(cowplot)
library(scales)
library(MASS)
library(taRifx)
source("~/git/mini_projects/PDAC_AhR/code/functions/TCGAanalyze_SurvivalKM2.R")

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/pdac_ahr/tcga'
outdir <- file.path(PDIR, 'results')
# icdo3: downloaded from https://seer.cancer.gov/icd-o-3/
icdo3_codes <- file.path(PDIR, "ref", "sitetype.icdo3.d20210607.csv")
# icd10: downloaded from https://raw.githubusercontent.com/kamillamagna/ICD-10-CSV/master/codes.csv
icd10_codes <- file.path(PDIR, "ref", "icd10_codes.csv")
# Downloaded from TCGAbiolinks
tcgadir <- file.path(PDIR, 'input/data/fpkm/obj')
clindir <- file.path(PDIR, 'input/data/clinical')
visualize <- TRUE

icdo3 <- read.csv(icdo3_codes, sep=",", header=TRUE, fill=FALSE, stringsAsFactors = FALSE, check.names = FALSE)
icd10 <- read.csv(icd10_codes, sep=",", header=FALSE, fill=FALSE, stringsAsFactors = FALSE, check.names = FALSE)
icd10$V3 <- paste0(icd10$V1, ".", icd10$V2)
projects <- c("TCGA-PAAD")
workflow_type <- list("counts"="HTSeq - Counts", "fpkm"="HTSeq - FPKM-UQ")

list.gene <- c('IDO1' = 'ENSG00000131203',
               'IDO2' = "ENSG00000188676",
               'TDO2' = "ENSG00000151790",
               "MRC1" = "ENSG00000260314")
assay_type <- "HTSeq - FPKM-UQ"

## Read in TCGA expression and metadata
setwd(tcgadir)
tcga_files <- list.files(path = tcgadir, pattern="obj.rds$")
tcga_idx <- sapply(projects, grepl, x=tcga_files)
tcga_files <- tcga_files[which(tcga_idx)]
objs <- lapply(tcga_files, readRDS)
names(objs) <- gsub("\\..*", "", tcga_files)
exprs <- lapply(objs, function(i) assay(i, assay_type))

## Read in TCGA clinical data
setwd(clindir)
clin_files <- list.files(path = clindir, pattern="clinical.rds$")
clin_idx <- sapply(projects, grepl, x=clin_files)
clin_files <- clin_files[which(clin_idx)]
clins <- lapply(clin_files, readRDS)
names(clins) <- gsub("\\..*", "", clin_files)

###################
#### Functions ####
annotateIcdCode <- function(p_df, icd_df, icd_code_col, icd_anno_col){
  m_idx <- match(p_df$Var2, icd_df[,icd_code_col])
  p_df$Var2 <- paste0("[", p_df$Var2, "] ",
                      icd_df[m_idx, icd_anno_col])
  
  m_idx <- match(p_df$Var1, icd_df[,icd_code_col])
  p_df$Var1 <- paste0("[", p_df$Var1, "] ",
                      icd_df[m_idx, icd_anno_col])
  p_df
}

###########################
#### Survival Analysis ####
dir.create(file.path(outdir, "survival_curves"), showWarnings = FALSE)

## Obtain the survival curves for each cancer type in TCGA
survivals <- lapply(projects, function(proj){
  # get Log2 FPKM + 1
  print(paste0(proj, " (", grep(proj, projects), "/", length(projects), ")"))
  log2exprs <- log10(exprs[[proj]] + 1)
  log2exprs <- log2exprs[list.gene,]
  rownames(log2exprs) <- names(list.gene)
  
  # Run survival analysis and get log-rank p-values
  if(visualize) pdf(file.path(outdir, "survival_curves", paste0(proj, "_skm.pdf")), height=10, width=6)
  skm <- TCGAanalyze_SurvivalKM2(clins[[proj]], log2exprs, names(list.gene), Survresult=TRUE,
                         threshcuts = c(0.5), caption=proj, p.cut = 1,
                         add.legend=T, add.pval=T)
  if(visualize) dev.off()
  skm <- as.data.frame(do.call(rbind, skm))
  rownames(skm) <- list.gene
  skm$HUGO <- names(list.gene[match(rownames(skm), list.gene)])
  return(skm)
})
# Format the data structures
names(survivals) <- projects
survival_pvals <- as.data.frame(survivals[[1]])
write.table(survival_pvals, file=file.path(outdir, "survival_curves", "survival_pvals.tsv"),
            quote=F, sep="\t", col.names=T, row.names = F)

###########################################
#### Visualize Expression and Survival ####
# Extract the log FPKM of expression for the genes of interest
i <- objs[[projects]]
expr <- assay(i, assay_type)
expr <- expr[match(list.gene, rownames(expr)),]
expr <- cbind(data.frame("id"=colnames(expr)),
              as.data.frame(t(log10(expr + 1))))
colnames(expr)[-1] <- names(list.gene)


# Reorder cancer types based on survival analysis pvals
melt_expr <- melt(as.data.frame(expr))

# Visualize gene expression in FPKM+1 ordered by survival pvals
scaleFUN <- function(x) sprintf("%.1f", x)
expr_dp <- ggplot(melt_expr, aes(x=variable, y=value)) + 
  geom_violin(trim=FALSE, fill="#74a9cf", scale = "width") +
  labs(title=paste0("Gene expression in ", projects), 
       x="Gene", y = "Log10(FPKM)") +
  geom_boxplot(width=0.2, fill="white") +
  theme_minimal() +
  scale_y_continuous(labels=scaleFUN) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

if(visualize){
  pdf(file.path(outdir, "survival_curves", "expr.pdf"), width=5, height=5)
  print(expr_dp)
  dev.off()
}


###################################################
#### Comparison to ssGSEA score of AhR Geneset ####
library(GSVA)
library(org.Hs.eg.db)

genome_gse <- org.Hs.eg.db
txby <- keys(genome_gse, 'SYMBOL')
gene_ids <- mapIds(genome_gse, keys=txby, column='ENSEMBL',
                   keytype='SYMBOL', multiVals="first")
ens_ids <- setNames(names(gene_ids), gene_ids)


ahr_geneset <- c('CYP1A1', 'CYP1B1', 'ABCA12', 'CERS3', 'MAP3K9', 'ITSN2',
                 'PNPT1', 'SDR9C7', 'EPN3', 'C17ORF109', 'ANO8', 'CD226', 
                 'AHRR', 'AC007639.1', 'SLC45A4', 'SECTM1', 'CA4')
ahr_geneset <- setNames(gene_ids[ahr_geneset], ahr_geneset)
ahr_geneset['C17ORF109'] <- 'ENSG00000204323'
ahr_geneset['AC007639.1'] <- 'ENSG00000263680'


# 
# .filterFeatures <- function(expr, method) {
#   
#   ## filter out genes with constant expression values
#   ## DelayedMatrixStats::rowSds() works for both base and 
#   ## DelayedArray matrices
#   sdGenes <- DelayedMatrixStats::rowSds(expr)
#   ## the following fixes this bug, see issues
#   ## https://github.com/rcastelo/GSVA/issues/54
#   ## https://github.com/HenrikBengtsson/matrixStats/issues/204
#   sdGenes[sdGenes < 1e-10] <- 0
#   if (any(sdGenes == 0) || any(is.na(sdGenes))) {
#     warning(sum(sdGenes == 0 | is.na(sdGenes)),
#             " genes with constant expression values throuhgout the samples.")
#     if (method != "ssgsea") {
#       warning("Since argument method!=\"ssgsea\", genes with constant expression values are discarded.")
#       expr <- expr[sdGenes > 0 & !is.na(sdGenes), ]
#     }
#   } 
#   
#   if (nrow(expr) < 2)
#     stop("Less than two genes in the input assay object\n")
#   
#   if(is.null(rownames(expr)))
#     stop("The input assay object doesn't have rownames\n")
#   
#   expr
# }
# X <- .filterFeatures(expr, 'ssgsea')
# geneSets <- GSVA:::.mapGeneSetsToFeatures(ahr_geneset, rownames(X))
# geneSets <- filterGeneSets(geneSets, min.sz=max(1, 1), max.sz=Inf)
# 
# n <- ncol(X)
# X <- matrix(ceiling(runif(50, 1, 1000)), ncol=10)
# R <- apply(X, 2, function(x, p) as.integer(rank(x)))
# Ra <- abs(R)^alpha
# j <- 1
# # order from highest ranked number (highest expr) to lowest
# geneRanking <- order(R[, j], decreasing=TRUE)
# sort.int(c(4,2))
# 
# .fastRndWalk <- function(gSetIdx, geneRanking, j, Ra) {
#   # gSetIdx <- geneSets[[1]]
#   n <- length(geneRanking)
#   k <- length(gSetIdx)
#   idxs <- sort.int(match(gSetIdx, geneRanking)) #Get the ranking of the gois
#   
#   # Ra[geneRanking[idxs], j] : adjusted geneRankings of GOIs
#   # (n-idxs + 1) : distance from the highest expr genes
#   stepCDFinGeneSet2 <- 
#     sum(Ra[geneRanking[idxs], j] * (n - idxs + 1)) /
#     sum((Ra[geneRanking[idxs], j]))    
#   
#   
#   # Random ranking from a geneset of equal size
#   stepCDFoutGeneSet2 <- 
#     (n * (n + 1) / 2 - sum(n - idxs + 1)) /
#     (n - k)
#   
#   walkStat <- stepCDFinGeneSet2 - stepCDFoutGeneSet2
#   
#   walkStat
# }


# Calculate ssGSEA for the AhR geneset
expr <- log10(assay(objs[[projects]], assay_type)+1)
non_symbol_idx <- which(is.na(ens_ids[rownames(expr)]))
expr <- expr[-non_symbol_idx,]
gsva_es <- gsva(expr, list("AhR"=ahr_geneset), method='ssgsea', 
                ssgsea.norm=T, verbose=T)
summary(gsva_es[1,])

# Calculate ssGSEA for the AhR geneset
expr <- expr[match(list.gene, rownames(expr)),]
expr <- cbind(data.frame("id"=colnames(expr)),
              as.data.frame(t(expr)))
colnames(expr)[-1] <- names(list.gene)

# Do the plotty plotties
melt_expr <- melt(expr)
melt_expr <- merge(melt_expr, data.frame("id"=colnames(gsva_es), "ssgsea"=gsva_es[1,]), by='id')
sp <- ggscatter(melt_expr, x = "ssgsea", y = "value",
                add = "reg.line", 
                add.params = list(color = "blue", fill = "lightgray"), 
                conf.int = TRUE) +
  facet_wrap(vars(variable), ncol=2) +
  ylab("Log10(FPKM)") + xlab("ssGSEA")
if(visualize){
  pdf(file.path(outdir, "survival_curves", "ssgsea.pdf"))
  print(sp + stat_cor(method = "spearman", label.x = -1, label.y = 6.5))
  dev.off()
}

