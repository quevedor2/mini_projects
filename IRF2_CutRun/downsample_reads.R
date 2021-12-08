library(GenomicRanges)
library(ggplot2)
library(reshape2)
library(dplyr)

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/downsample_cutandrun'
outdir <- file.path(PDIR, "results/plots")
seacr_dir <- file.path(PDIR, "results/peaks/seacr_single")
setwd(seacr_dir)
dir.create(outdir)


params <- list("heni"=list("pattern"="Heni.*CTCF",
                           "parent"="28M"),
               "hpb"=list("pattern"="HPB",
                          "parent"="2M"))
dataset <- 'hpb'

## Order the files
files <- list.files(pattern = params[[dataset]]$pattern)
parent_file <- grep(params[[dataset]]$parent, files, value = T)

ord <- order(sapply(strsplit(files, "_"), function(i) {
  gsub("K", "000", i[[2]]) %>%
    gsub("M", "000000", .) %>%
    as.integer(.)
}))
files <- files[ord]

## Read in all peak files
grl <- lapply(files, function(f){
  f <- read.table(f, sep="\t", stringsAsFactors = F)
  colnames(f) <- c('chr', 'start', 'end', 'integral', 'height', 'range')
  makeGRangesFromDataFrame(f, keep.extra.columns = T)
})
names(grl) <- files

## Assemble the parent file
parent_gr <- grl[[parent_file]]
parent_gr <- parent_gr[order(parent_gr$height),]
parent_gr <- parent_gr[seqnames(parent_gr) %in% standardChromosomes(parent_gr),]

## Create a cnt matrix of overlap
peaks_mat <- sapply(grl, function(gr){
  peaks_found <- rep(0, length(parent_gr))
  ov <- findOverlaps(parent_gr, gr)
  if(length(ov)>0) peaks_found[unique(queryHits(ov))] <- 1
  return(peaks_found)
})
colnames(peaks_mat) <- sapply(strsplit(files, "_"), function(i) i[[2]])
rownames(peaks_mat) <- make.unique(as.character(ceiling(parent_gr$height)))

## Do the plotties - Raw overlap
peaks_melt <- melt(peaks_mat)
peaks_melt$value <- factor(as.character(peaks_melt$value))
levels(peaks_melt$value) <- c('No-overlap', 'Overlap')
peaks_melt$Var2 <- factor(as.character(peaks_melt$Var2),
                          levels = colnames(peaks_mat))
peaks_melt$Var1 <- factor(as.character(peaks_melt$Var1),
                          levels=rownames(peaks_mat))

pdf(file.path(outdir, paste0(dataset, "_overlap_boolean.pdf")), height = 15)
ggplot(peaks_melt, aes(Var2, Var1, fill=value)) + 
  geom_tile()
dev.off()

## Split into quantiles and plot overlap
peak_height <- as.numeric(rownames(peaks_mat))
decile <- seq(0, 0.9, by=0.1)
frac_mat <- sapply(decile, function(dec){
  peaks_sub <- peaks_mat[peak_height > quantile(peak_height, dec),]
  frac <- colSums(peaks_sub)/nrow(peaks_sub)
  c(nrow(peaks_sub), frac)
})
frac_mat <- as.data.frame(frac_mat)
colnames(frac_mat) <- decile
n_mat <- frac_mat[1,]
frac_mat <- frac_mat[-1,]
frac_melt <- melt(t(frac_mat))

## Plotties
pdf(file.path(outdir, paste0(dataset, "_fractional_overlap.pdf")))
ggplot(frac_melt, aes(x=Var1, y=value, color=Var2)) +
  geom_line(size=1, alpha=0.9, linetype=1) + 
  theme_minimal() +
  xlab("decile-cutoff") + ylab("")
ggplot(melt(n_mat), aes(x=variable, y=value)) +
  geom_bar(stat='identity') +
  theme_minimal() +
  xlab("decile-cutoff") + ylab("# of peaks")
dev.off()