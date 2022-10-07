library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(dplyr)

pdir <- '/cluster/projects/pughlab/projects/KOMBAT/test/Codex_segs'
setwd(pdir)
seg_dir <- pdir

#Parameters
bin_size <- 1000000          # bin size (bp) if merging based on uniform bins 
                             # across the genome
merge_mode <- 'intersect'    # intersect segments or create a tiled genome
seg_col <- c('Segment_Mean') # column name of seg files containing value (iterates)
sample_col <- 'Sample'       # column name of seg file containing its name

# Get reference genome
genome <- BSgenome.Hsapiens.UCSC.hg38
genome_gr <- GRanges(genome@seqinfo) %>% 
  keepStandardChromosomes(., pruning.mode = 'coarse')
seqlevelsStyle(genome_gr) <- 'UCSC'

# Read in all seg files
segl <- lapply(list.files(seg_dir, pattern='seg$|.seg.txt$'), function(seg_f){
  read.table(seg_f, header=T, check.names = F, stringsAsFactors = F) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = T)
})

# Create a reference seg object
if(merge_mode == 'tile'){
  #WARNING, NOT FULLY IMPLEMENTED YET
  ref_gr <- tileGenome(seqlengths(genome_gr), tilewidth=bin_size) %>%
    unlist
} else if (merge_mode == 'intersect'){
  # Find all unique breakpoints
  grl <- c(list("ref"=genome_gr), segl) %>% 
    as(., "GRangesList") %>% 
    unlist
  
  # Create a master reference based on the unique breakpoints
  grl_chr <- split(grl, seqnames(grl))
  ref_gr <- lapply(grl_chr, function(gr_chr) {
    bp <- c(start(gr_chr), end(gr_chr)) %>% sort
    data.frame("chr"=rep(seqnames(gr_chr)@values, length(bp)-1) %>%
                 as.character,
               "start"=bp[-length(bp)],
               "end"=bp[-1]-1) %>%
      makeGRangesFromDataFrame
  }) %>%
    as(., "GRangesList") %>%
    unlist
}

# Map the seg files to the reference seg object
merge_seg_gr <- ref_gr
segmean_mat <- lapply(segl, function(seg_i){
  ov_idx <- findOverlaps(ref_gr, seg_i)
  sample <- elementMetadata(seg_i)[,sample_col] %>% 
    unique
  
  # Find max number of segments that overlap one reference bin
  max_n <- max(table(queryHits(ov_idx)))
  
  # Populate segmean matrix, iterating over multiple intervals per ref_bin
  segmean <- data.frame(matrix(nrow=length(ref_gr), ncol=max_n))
  
  ov_idx_i <- ov_idx
  idx <- 1
  dup_status <- TRUE
  while(dup_status){
    # Split based on bin with multiple overlaps 
    dup_status <- any(duplicated(queryHits(ov_idx_i)))
    spl_ov <- split(ov_idx_i, duplicated(queryHits(ov_idx_i)))
    if(dup_status) ov_idx_i <- spl_ov[['TRUE']]
    
    # Fill in segmean_matrix with all non-duplicated entries
    non_dup_ov <- spl_ov[['FALSE']]
    segmean[queryHits(non_dup_ov),idx] <- elementMetadata(seg_i)[subjectHits(non_dup_ov),
                                                                 seg_col]
    
    idx <- idx+1
  }
  
  # Reduce segmean into one value
  segmean <- apply(segmean, 1, median, na.rm=T) %>% 
    as.data.frame %>%
    setNames(., sample)
  
  return(segmean)
}) %>% 
  do.call(cbind, .) %>%
  as.data.table
elementMetadata(merge_seg_gr) <- segmean_mat

saveRDS(merge_seg_gr, file=file.path(seg_dir, "merge.rds"))
write.table(data.frame(merge_seg_gr), file=file.path(seg_dir, "merge.seg"),
            sep="\t", col.names = T, row.names = F, quote = F)
