library(dplyr)
library(ggplot2)
PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/wt_cutandrun/snakemake/results/alignment/blast'
reads  <- read.csv(file.path(PDIR, 'data', 'map_unmap', 'cnts.txt'), 
                   header=FALSE)
reads$V1 <- gsub("\\/", "", reads$V1)
reads_sample <- split(reads, reads$V1)
map_unmaps <- lapply(reads_sample, function(rs){
  sample <- unique(rs$V1)
  print(sample)
  blastDIR <- file.path(PDIR, "blast", sample)
  clip <- 10
  
  map_unmap <- lapply(c("mapped", "unmapped"), function(m){
    print(m)
    mapped <- read.table(file.path(blastDIR, paste0(m, ".counts.txt")), sep=" ", 
                         header=FALSE, fill=TRUE)
    hs_idx <- c(grep("Homo_sapiens", mapped$V1),
                grep("Human", mapped$V1))
    if(length(hs_idx)>0){
      mapped <- rbind(c("Homo_sapiens", sum(na.omit(mapped$V2[hs_idx]))),
                      mapped[-hs_idx,])
      mapped$V2 <- as.integer(mapped$V2)
    }
    
    mapped <- rbind(mapped[c(1:9),],
                    c("Other", sum(na.omit(mapped[-c(1:9),'V2']))))
    mapped$V2 <- as.integer(mapped$V2)
    
    mapped <- as.data.frame(mapped)
    mapped$frac <- (mapped$V2) / sum(na.omit(mapped$V2))
    midx <- match(m, rs$V2)
    reads <- rs[midx,]$V3
    mapped$est <- ceiling(reads * mapped$frac)
    mapped$sample <- sample
    mapped$mapping <- m
    return(mapped)
  })
  return(as.data.frame(do.call(rbind, map_unmap)))
})

map_unmaps <- as.data.frame(do.call(rbind, map_unmaps))
map_unmaps$V1 <- gsub("^.*PREDICTED:_", "", map_unmaps$V1)
map_unmaps$V1 <- gsub("^(.*?_.*?)_.*", "\\1", map_unmaps$V1)

x <- split(map_unmaps, map_unmaps$V1)
as.matrix(sort(sapply(x, function(i) sum(i$est)), decreasing=TRUE))
as.matrix(sort(sapply(x, function(i) mean(i$frac)), decreasing=TRUE))

p <- ggplot(map_unmaps, aes(fill=V1, y=est, x=sample)) + 
  geom_bar(position="stack", stat="identity") +
  ylab("Read counts") + xlab("Sample") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(rows = vars(mapping))

pdf(file.path(PDIR, "blast", "map-unmap_blast.pdf"), width=15)
p
dev.off()