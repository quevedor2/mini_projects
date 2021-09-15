library(GenomicRanges)
library(ggplot2)
library(lsa)
library(gplots)
library(ggdendro)

caller <- 'seacr' # 'macs2'

if(caller=='macs2'){
  PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/sabelo_irf2/cutandrun/results/alignment/macs'
  pattern <- "_peaks.narrowPeak$"
  cols <- c('chr', 'start', 'end', 'name', 'score', 'strand',
            'max', 'p', 'q', 'peak')
} else {
  PDIR <- "/cluster/projects/mcgahalab/data/mcgahalab/sabelo_irf2/cutandrun/results/peaks/seacr_single"
  pattern <- ".stringent.bed$"
  cols <- c('chr', 'start', 'end', 'total', 'max', 'max_range')
}


setwd(PDIR)
bed_files <- list.files(pattern=pattern)
merge_bgr <- read.table("merge.bed", header = FALSE)
colnames(merge_bgr) <- c('chr', 'start', 'end')
merge_bgr <- makeGRangesFromDataFrame(merge_bgr)
bgr <- lapply(bed_files, function(bed){
  bed <- read.table(bed, header=FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(bed) <- cols
  if(caller=='macs2') bed$max_range <- with(bed, paste0(chr, ':', start, '-', end))
  maxbed <- data.frame('chr'=gsub(":.*", "", bed$max_range),
                       "start"=gsub("^.*:(.*)-.*$", "\\1", bed$max_range),
                       "end"=gsub("^.*-", "", bed$max_range),
                       "peak"=as.integer(bed$max))
  # gr <- makeGRangesFromDataFrame(bed, keep.extra.columns = FALSE)
  gr <- makeGRangesFromDataFrame(maxbed, keep.extra.columns = TRUE)
  return(gr)
})
names(bgr) <- gsub(pattern, "", bed_files)

peaks <- rep(0, length(merge_bgr))
peak_mat <- sapply(bgr, function(gr0){
  ol <- findOverlaps(gr0, merge_bgr)
  peaks[subjectHits(ol)] <- gr0[queryHits(ol),]$peak
  peaks
})
keep_row_idx <- which(rowSums(peak_mat != 0) > 3)
peak_mat_filt <- peak_mat[keep_row_idx,]
peak_mat_adj <- apply(peak_mat_filt, 2, scale, center=TRUE, scale=TRUE)
pca <- prcomp(t(peak_mat_adj))
percent_var <- round(pca$sdev^2 / sum( pca$sdev^2 ),4)  # Calculate variance for each PC
pc_cutoff <- min(which(cumsum(percent_var) > 0.95))     # Identify 95% variance cutoff

## Fragment plot
frags <- lapply(bed_files, function(bed){
  bed <- read.table(bed, header=FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(bed) <- cols
  if(caller=='macs2') bed$max_range <- with(bed, paste0(chr, ':', start, '-', end))
  maxbed <- data.frame('chr'=gsub(":.*", "", bed$max_range),
                       "start"=gsub("^.*:(.*)-.*$", "\\1", bed$max_range),
                       "end"=gsub("^.*-", "", bed$max_range),
                       "peak"=as.integer(bed$max))
  maxgr <- makeGRangesFromDataFrame(maxbed, keep.extra.columns = FALSE)
  gr <- makeGRangesFromDataFrame(bed, keep.extra.columns = FALSE)
  return(data.frame('peak'=width(gr), 'max'=width(maxgr)))
})
names(frags) <- gsub(pattern, "", bed_files)

frag_df <- do.call(rbind, frags)
frag_df$sample <- rep(names(frags), sapply(frags, nrow))
winsorize <- function(a, q){
  thresh <- quantile(a,q)
  a[a > thresh] <- thresh
  return(a)
}
frag_df$peak <- winsorize(frag_df$peak, 0.99)
frag_df$max <- winsorize(frag_df$max, 0.99)


pfrag <-  ggplot(data=frag_df, aes(x=peak, color=sample)) +
  geom_density(alpha=0.3,size=1)+ 
  labs(x= "Peak size") +
  theme_minimal()
pmax <-  ggplot(data=frag_df, aes(x=max, color=sample)) +
  geom_density(alpha=0.3,size=1)+ 
  labs(x= "Peak size") +
  theme_minimal()
pdf("~/xfer/frags.pdf")
pfrag
pmax
dev.off()

## Cosine similarity plots
cosine_mat <- apply(peak_mat_filt, 2, function(x){
  apply(peak_mat_filt, 2, function(y){
    cosine(x, y)
  })
})
summary(as.numeric(cosine_mat))
df <- as.data.frame(t(apply(cosine_mat,2,function(i) names(head(sort(i, decreasing = TRUE))))))[,-1]

pdf("~/xfer/cosine.pdf")
heatmap.2(cosine_mat, trace='none', col=viridis::viridis,
          cexRow=rel(0.8), cexCol=rel(0.8), margins=c(8,8))
dev.off()

## Principal component analysis
o <- as.data.frame(pca$x)
o$name <- rownames(o)
o$group <- gsub("^.*_(.*)_.*$", "\\1", rownames(o))
o$group <- gsub("^.*_(.*)$", "\\1", rownames(o))
p <-   ggplot(o, aes(x= PC1, y= PC2, color=group, label=name))+
  geom_point() +
  geom_text(aes(label=name),vjust="inward",hjust="inward", size=3) +
  scale_x_continuous(expand = c(.1, .1))
pdf("~/xfer/pca.pdf")
p
barplot(percent_var, ylim=c(0,0.5), ylab="% variance", xlab="PCs", 
        names.arg=paste0("PC", c(1:length(percent_var))), las=2, cex.names=0.7)
dev.off()


hc <- hclust(dist(pca$x[,1:pc_cutoff]))
dd.row <- as.dendrogram(hc)
ddata_x <- dendro_data(dd.row)

# Visualize the hclust and colour according to mDC/pDC/Tfh
p2 <- ggplot(segment(ddata_x)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
labs <- label(ddata_x)
labs$group <- gsub("^.*_(.*)$", "\\1", labs$label)

pdf(file.path("~/xfer", "hclust.pdf"))
n <- length(unique(labs$group))
max_y <- ceiling(max(ddata_x$segments$yend) /10) * 10
p2 + geom_text(data=label(ddata_x),
               aes(label=label, x=x, y=0, colour=labs$group),
               size=3,hjust=1,nudge_y = -5) +
  coord_flip() +
  ylim(-1 * ceiling(max_y/2), max_y) +
  theme_bw()
dev.off()


library(rjson)
json_file <- "~/Desktop/chr17/irf2.json"
json_data <- fromJSON(paste(readLines(json_file), collapse=""))
genes <- sapply(json_data$associations, function(i) i$gene$symbol)

hom <- read.table("~/Desktop/chr17/HOM_MouseHumanSequence.rpt", sep="\t",
                  header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
human_hom <- hom[which(hom$Symbol %in% genes),]
hm_hom <- hom[which(hom$`DB Class Key` %in% human_hom$`DB Class Key`),]
hm_hom <- split(hm_hom, hm_hom$`Common Organism Name`)
mouse_hom <- hm_hom$`mouse, laboratory`
mouse_hom$chr <- gsub(":.*", "", mouse_hom$`Genomic Coordinates (mouse: , human: )`)
mouse_hom_chr <- split(mouse_hom, mouse_hom$chr)
sapply(mouse_hom_chr, nrow)


x <- read.table("x.bed", header=FALSE)
x <- x[order(x$V9, decreasing = TRUE),]
x$coord <- paste0(rep("chr17:", nrow(x)), x$V2, "-", x$V3)
cat(paste(head(split(x, x$V1)[['WT_Tumor_IRF2_peaks.narrowPeak:17']]$coord, n=6), collapse="\n"))

x <- x[order(x$V5, decreasing = TRUE),]
x$coord <- paste0(rep("chr17:", nrow(x)), x$V2-100, "-", x$V3+100)
cat(paste(head(split(x, x$V1)[['WT_LN_IRF2.stringent.bed:17']]$coord, n=6), collapse="\n"))

id <- 'WT_Spleen_IRF2'
fr <- frags[[id]]
q <- quantile(fr$peak, 0.999)
which(fr$peak > q)
bgr[[id]][which(fr$peak > q),]
chr17 <- bgr[[id]]
chr17 <- chr17[which(seqnames(chr17) == '17'),]
chr17[order(chr17$peak, decreasing = TRUE),]