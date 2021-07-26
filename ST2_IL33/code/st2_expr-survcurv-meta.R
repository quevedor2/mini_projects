###########################
#### Read in TCGA Data ####
library(reshape2)
library(ggplot2)
library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)
library(gridExtra)
library(cowplot)

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/tcga'
outdir <- file.path(PDIR, 'results')
tcgadir <- file.path(PDIR, 'input/data/fpkm/obj')
clindir <- file.path(PDIR, 'input/data/clinical')

projects <- c("TCGA-BRCA", "TCGA-GBM", "TCGA-OV", "TCGA-LUAD", "TCGA-UCEC",
              "TCGA-KIRC", "TCGA-HNSC", "TCGA-LGG", "TCGA-THCA", "TCGA-LUSC",
              "TCGA-PRAD", "TCGA-SKCM", "TCGA-COAD", "TCGA-STAD", "TCGA-BLCA",
              "TCGA-LIHC", "TCGA-CESC", "TCGA-KIRP", "TCGA-SARC", "TCGA-LAML", 
              "TCGA-PAAD", "TCGA-ESCA", "TCGA-PCPG", "TCGA-READ", "TCGA-TGCT",
              "TCGA-THYM", "TCGA-KICH", "TCGA-ACC", "TCGA-MESO", "TCGA-UVM",
              "TCGA-DLBC", "TCGA-UCS", "TCGA-CHOL")
workflow_type <- list("counts"="HTSeq - Counts", "fpkm"="HTSeq - FPKM-UQ")

list.gene <- c('IL1RL1' = 'ENSG00000115602', # IL1RL1
               'PDCD1_1' = "ENSG00000188389", # PDCD1
               'PDCD1_2' = "ENSG00000276977", # PDCD1
               'CTLA4' = "ENSG00000163599", # CTLA4
               'TOX' = "ENSG00000198846", # TOX
               'KLRG1' = "ENSG00000139187") # KLRG1
assay_type <- "HTSeq - FPKM-UQ"

## Read in TCGA expression and metadata
setwd(tcgadir)
tcga_files <- list.files(path = tcgadir, pattern="obj.rds$")
objs <- lapply(tcga_files, readRDS)
names(objs) <- gsub("\\..*", "", tcga_files)
exprs <- lapply(objs, function(i) assay(i, assay_type))

## Read in TCGA clinical data
setwd(clindir)
clin_files <- list.files(path = clindir, pattern="clinical.rds$")
clins <- lapply(clin_files, readRDS)
names(clins) <- gsub("\\..*", "", clin_files)

###########################
#### Survival Analysis ####
## Obtain the survival curves for each cancer type in TCGA
survivals <- lapply(projects, function(proj){
  # get Log2 FPKM + 1
  print(paste0(proj, " (", grep(proj, projects), "/", length(projects), ")"))
  log2exprs <- log2(exprs[[proj]] + 1)
  
  # Run survival analysis and get log-rank p-values
  dir.create(file.path(outdir, "survival_curves"), showWarnings = FALSE)
  pdf(file.path(outdir, "survival_curves", paste0(proj, "_skm.pdf")))
  skm <- TCGAanalyze_SurvivalKM(clins[[proj]], log2exprs, list.gene, Survresult=TRUE,
                         ThreshTop = 0.67, ThreshDown = 0.33, p.cut = 1)
  dev.off()
  skm$HUGO <- names(list.gene[match(rownames(skm), list.gene)])
  return(skm)
})
# Format the data structures
names(survivals) <- projects
st2_pvals <- as.data.frame(t(sapply(survivals, function(s) s[list.gene['IL1RL1'],])))
st2_pvals <- st2_pvals[order(unlist(st2_pvals$pvalue)),]

#########################################
#### Isolate Expression and Metadata ####
groups <- c('state', 'tissue_type', 'ajcc_pathologic_stage', 'tumor_stage', 
            'ajcc_pathologic_t', 'ajcc_pathologic_n', 'ajcc_pathologic_m',
            'morphology', 'tumor_grade', 'progression_or_recurrence', 
            'paper_Tumor.Type', 'paper_pathologic_stage', 'paper_Tumor_Grade')
goi <-'IL1RL1'
st2_dat <- lapply(objs, function(i) {
  metadata <- colData(i)
  expr <- assay(i, assay_type)
  expr <- expr[match(list.gene[goi], rownames(expr)),]
  cbind(data.frame("index"=seq_along(expr), 
                   "expr"=log2(expr+1)), 
        metadata[,which(colnames(metadata) %in% groups)])
})
# Reduce expression and metadat a by cancer types
st2_exprs <- Reduce(function(x,y) merge(x,y,by='index', all=TRUE), lapply(st2_dat, function(i) i[,1:2]))
st2_exprs <- st2_exprs[,-1]
colnames(st2_exprs) <- names(st2_dat)
st2_exprs <- st2_exprs[,rownames(st2_pvals)]

###########################################
#### Visualize Expression and Survival ####
# Reorder cancer types based on survival analysis pvals
melt_st2_exprs <- melt(st2_exprs)
melt_st2_exprs <- melt_st2_exprs[-which(is.na(melt_st2_exprs$value)),]
melt_st2_exprs$value <- as.numeric(as.character(melt_st2_exprs$value))

# Visualize gene expression in FPKM+1 ordered by survival pvals
scaleFUN <- function(x) sprintf("%.1f", x)
expr_dp <- ggplot(melt_st2_exprs, aes(x=variable, y=value)) + 
  geom_violin(trim=FALSE, fill="#74a9cf", scale = "width") +
  labs(title="IL1RL1 gene expression",x="TCGA Projects", y = "Log2(FPKM + 1)") +
  geom_boxplot(width=0.2, fill="white") +
  theme_minimal() +
  scale_y_continuous(labels=scaleFUN) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(file.path("~/xfer", "expr.pdf"), width=12, height=5)
expr_dp
dev.off()

# Visualize survival pval ordered by survival pvals
st2_pvals$groups <- factor(rownames(st2_pvals), levels=rownames(st2_pvals))
st2_pvals$sig <- st2_pvals$pvalue < 0.05
surv_dp <- ggplot(st2_pvals, aes(x=groups, y=pvalue, fill=sig)) + 
  geom_bar(stat="identity", color="black") + 
  scale_fill_manual(values=c("#ef8a62", "#67a9cf")) +
  labs(title="IL1RL1 Survival p-values",x="TCGA Projects", y = "Log-rank p-value") +
  theme_minimal() +
  coord_cartesian(ylim=c(0,1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none")
  
pdf(file.path(outdir, "expr_skm.pdf"), width=12, height=9)
grid.arrange(expr_dp + theme(axis.title.x=element_blank(),
                             axis.text.x=element_blank(),
                             axis.ticks.x=element_blank()), 
             surv_dp, 
          ncol = 1, nrow = 2)
dev.off()


##########################################
#### Expression and Metadata analysis ####
aov_pvals <- lapply(st2_dat, function(st2){
  grpids <- colnames(st2)[-1]
  aovs <- lapply(setNames(grpids, grpids), function(g){
    res.aov <- tryCatch({
      aov(as.formula(paste0("expr ~ ", g)), data = st2)
    }, error=function(e){NA})
    res.aov
  })
  pval <- sapply(aovs, function(i)  tryCatch(summary(i)[[1]]$`Pr(>F)`[1], error=function(e){NA}))
  data.frame("category"=grpids, "p"=pval)
})
pvals <- Reduce(function(x,y) as.data.frame(merge(x,y,by='category', all=TRUE)), aov_pvals)
pvals <- matrix(unlist(pvals), ncol=ncol(pvals))
rownames(pvals) <- pvals[,1]
pvals <- pvals[,-1]
colnames(pvals) <-names(st2_dat)
storage.mode(pvals) <- 'numeric'
pvals <- round(pvals, 4)

## Heatmap visualization of ANOVA analysis
meta_pvals <- pvals[,rownames(st2_pvals)]
melt_meta_pvals <- melt(meta_pvals)
meta_dp <- ggplot(melt_meta_pvals, aes(Var2, Var1, fill= value)) + 
  geom_tile() +
  labs(title="ANOVA p-value for IL1RL1 expression",x="TCGA Projects", y = "Metadata") +
  scale_fill_gradientn(colours = c("#cb181d","#2171b5", "#6baed6", "#c6dbef", "#f7fbff"), 
                       values = rescale(c(0, 0.05, 0.1, 0.15, 1)),
                       guide = "colorbar", limits=c(0,1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom")


pdf(file.path("~/xfer", "expr_skm_meta.pdf"), width=12, height=14)
plot_grid(expr_dp + theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()), 
          surv_dp + theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()),
          meta_dp, 
          ncol=1, align="v")
dev.off()


# 
# kws <- lapply(setNames(groups, groups), function(g){
#   res.kw <- tryCatch({
#     kruskal.test(as.formula(paste0("expr ~ ", g)), data = df)
#   }, error=function(e){NA})
#   res.kw
# })
# aovs <- lapply(setNames(groups, groups), function(g){
#   res.aov <- tryCatch({
#     aov(as.formula(paste0("expr ~ ", g)), data = df)
#   }, error=function(e){NA})
#   res.aov
# })
# pvals <- sapply(aovs, function(g){
#   tryCatch({
#     summary(g)[[1]]$`Pr(>F)`[1]
#   }, error=function(e){NA})
# })
# paired <- lapply(aovs, function(g){
#   tryCatch({
#     TukeyHSD(res.aov)
#   }, error=function(e){NA})
# })
# paired_t <- lapply(setNames(groups, groups), function(g){
#   pairwise.t.test(df$expr, df[,g], p.adjust.method = "fdr")
# })
# paired[names(which(pvals < 0.05))]
# paired_t[names(which(pvals < 0.05))]
# 
# summary(res.aov)
# 
# # For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
# BRCARnaseq_CorOutliers <- TCGAanalyze_Preprocessing(BRCARnaseqSE)