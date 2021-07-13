# This code is meant to look at the bulkRNAseq data from the Withers
# mDC, pDC, and Tfh cohort. It was processed using the rna-seq-star-deseq2
# workflow: [https://github.com/quevedor2/rna-seq-star-deseq2]
#
# It assumes you're in the main directory and that the counts matrix
# can be found in results/deseq2/all.rds
library(DESeq2)
library(ggplot2)
library(umap)
library(ggdendro)

inpath <- 'results/deseq2'
outpath <- 'results/adhoc/pca'
dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

# Read and transform the count data using regularized logarithm
dds <- readRDS(file.path(inpath, "all.rds"))
counts2 <- counts <- rlog(dds, blind=FALSE)
counts2$condition <- factor(gsub("-.*", "", names(counts$sizeFactor)))  # Relabel groups based on mDC/pDC/Tfh

###########################
#### Base PCA Analysis ####
# Principle component analysis using the top 500 most variable genes
rv <- rowVars(assay(counts))
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
pca <- prcomp(t(assay(counts)[select,]))
percent_var <- round(pca$sdev^2 / sum( pca$sdev^2 ),4)  # Calculate variance for each PC
pc_cutoff <- min(which(cumsum(percent_var) > 0.95))     # Identify 95% variance cutoff

# Generate the PCA plots using DESeq2 plotting function
outs <- lapply(list(counts, counts2), function(counts){
  plotPCA(counts, intgroup='condition', returnData = TRUE)
})
pdf(file.path(outpath, "pca.pdf"))
## Scatterplot for top 2 PCs
ps <- lapply(outs, function(o){
  ggplot(o, aes(x= PC1, y= PC2, color=group, label=name))+
    geom_point() +
    geom_text(aes(label=name),hjust=0, vjust=0) +
    xlim(-60, 60) +
    ylim(-50, 100)
})
ps[[1]]
ps[[2]]

## Scree plot
barplot(percent_var[1:50], ylim=c(0,0.5), ylab="% variance", xlab="PCs", 
        names.arg=paste0("PC", c(1:length(percent_var[1:50]))), las=2, cex.names=0.7)
dev.off()


#################################
#### 95% variance clustering ####
# Run a simple hierarchical clustering using euclidean distances
colors <- c("#1b9e77", "#d95f02", "#7570b3")
hc <- hclust(dist(pca$x[,1:pc_cutoff]))
dd.row <- as.dendrogram(hc)
ddata_x <- dendro_data(dd.row)

# Visualize the hclust and colour according to mDC/pDC/Tfh
p2 <- ggplot(segment(ddata_x)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
labs <- label(ddata_x)
labs$group <- gsub("-.*", "", labs$label)
pdf(file.path(outpath, "hclust.pdf"))
p2 + geom_text(data=label(ddata_x),
               aes(label=label, x=x, y=0, colour=labs$group),
               size=3,hjust=1,nudge_y = -5) +
  scale_colour_manual(values=c("#1b9e77", "#d95f02", "#7570b3")) +
  coord_flip() + 
  ylim(-50, 300) +
  theme_bw()
dev.off()


# Run a UMAP dimensional reduction using varying number of neighbours
custom.config = umap.defaults
custom.config$random_state = 123
custom.config$min_dist = 0.3
custom.config$n_epochs = 500

pdf(file.path(outpath, "umap.pdf"))
sapply(c(2:15), function(n){
  custom.config$n_neighbors = n
  immune_umap <- umap(pca$x[,1:pc_cutoff], config=custom.config)
  layout <- immune_umap$layout
  immune_labels <- factor(gsub("-.*", "", rownames(layout)))
  xylim <- range(layout)
  
  plot(layout[,1], layout[,2], xlim=xylim+c(-1, 1), ylim=xylim, type='n', 
       xlab='UMAP1', ylab='UMAP2', main=paste0("n_neigh=", n))
  points(layout, col=alpha(colors[as.integer(immune_labels)],0.70), pch=16)
  text(layout+0.1, labels=rownames(layout), cex=0.7, 
       col = colors[as.integer(immune_labels)], pos = 4)
  legend("bottomleft", legend=as.character(unique(immune_labels)), inset=0.03,
         col=colors[as.integer(unique(immune_labels))],
         bty="n", pch=16, cex=0.7)
})
dev.off()
