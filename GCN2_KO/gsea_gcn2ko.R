library(reshape2)
library(enrichplot)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library("clusterProfiler")
library("enrichplot")
library(DOSE)
library(ggrepel)
library(msigdbr)
library(org.Mm.eg.db)

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/Xin/GCN2'
gcn2 <- read.csv(file.path(PDIR, "data", "2020-10-25_BMDM_Triplicate_Normalised_p-value.csv"))
colnames(gcn2)[c(11,12)] <- c('pval', 'fold_change')
gcn2$qval <- p.adjust(gcn2$pval, method='fdr')
gcn2$sig <- gcn2$qval < 0.1 & abs(gcn2$fold_change) > 3
gcn2$log10q <- -1*log10(gcn2$qval)
gcn2_sig <- gcn2[which(gcn2$sig),]

## Volcano Plot
dir.create(file.path(PDIR, "results", "plots"), recursive = T, showWarnings = F)
gp <- ggplot(data=gcn2, aes(x=fold_change, y=log10q, color=sig)) +
  geom_point(alpha = 0.6) + 
  theme_minimal() +
  xlim(-10,10) +
  geom_text_repel(data=gcn2_sig, aes(label=Gene.names))
pdf(file.path(PDIR, "results", "plots", "volcano_plot.pdf"))
gp
dev.off()


## GSEA
m_df <- msigdbr(species = "Mus musculus")
msig_ds <- list("pathway"=list("C2", "CP:REACTOME"))
msig_gsl <- lapply(msig_ds, function(mds){
  msigdbr(species = "Mus musculus", category = mds[[1]], subcategory = mds[[2]]) %>% 
    dplyr::select(gs_name, entrez_gene)
})


genome_gse <- org.Mm.eg.db
txby <- keys(genome_gse, 'SYMBOL')
gene_ids <- mapIds(genome_gse, keys=txby, column='ENTREZID',
                   keytype='SYMBOL', multiVals="first")

# Make gene list
gene_list <- gcn2$fold_change
names(gene_list) <- gene_ids[gcn2$Gene.names]
gene_list<-sort(na.omit(gene_list), decreasing = TRUE)

gsea <- GSEA(gene_list, TERM2GENE =msig_gsl[[1]], 
             pvalueCutoff = 1)
# gsea@result[,-c(1,2,11)]

pdf(file.path(PDIR, "results", "plots", "gsea_nes.pdf"), width = 12)
gsea %>%
  filter(grepl("MTOR|AKT|PI3K|S6K|GCN2", ID) | p.adjust < 0.01) %>%
  mutate(direction=NES>0) %>%
  group_by(direction) %>%
  slice(1:15) %>%
  enrichplot:::barplot.enrichResult(x='NES', showCategory=40) + 
  theme(text = element_text(size=8),
        axis.text.y = element_text(size=8))  +
  xlim(-4,4) +
  xlab("NES") +
  ggtitle("PI3K/AKT/MTOR gene sets")

gsea %>%
  filter(grepl("MTOR|AKT|PI3K|S6K|GCN2", ID)) %>%
  enrichplot:::barplot.enrichResult(x='NES', showCategory=30) + 
  theme(text = element_text(size=8),
        axis.text.y = element_text(size=8))  +
  xlim(-2,2) +
  xlab("NES") +
  ggtitle("PI3K/AKT/MTOR gene sets")
dev.off()

## Gene by Gene-Set Heatmap of LFC
gseax <- setReadable(gsea, 'org.Mm.eg.db', 'ENTREZID')
colors <- c('navyblue', 'lightgrey', 'darkred')
maxval <- max(abs(gene_list))
b <- c(-ceiling(maxval), 0, ceiling(maxval))
pdf(file.path(PDIR, "results", "plots", "gsea_gene_heatmap.pdf"), width = 12)
gseax %>%
  filter(grepl("MTOR|AKT|PI3K|S6K|GCN2", ID)) %>%
  heatplot(showCategory=20, foldChange=gene_list) +
  theme_minimal() +
  scale_fill_gradientn(limits = c(min(b),max(b)), colors = colors, 
                       breaks = b, labels = format(b)) +
  theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1))
dev.off()


## Gene by Gene-Set Network Map
pdf(file.path(PDIR, "results", "plots", "gsea_emap.pdf"), width = 10, height = 10)
gsea %>%
  filter(grepl("MTOR|AKT|PI3K|S6K|GCN2", ID) | p.adjust < 0.01) %>%
  pairwise_termsim() %>%
  emapplot(layout='nicely', cex_category=1, 
           cex_label_category=0.4, showCategory = 100)

gsea %>%
  filter(grepl("MTOR|AKT|PI3K|S6K|GCN2", ID)) %>%
  pairwise_termsim() %>%
  emapplot(layout='nicely', cex_category=1, 
           cex_label_category=0.5)
dev.off()
