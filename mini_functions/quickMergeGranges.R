# Map every GRanges object ot a reference GRange object and 
# combine the median value of the metacol to its elementMetadata
quickMergeGranges <- function(segl, ref_gr, metacol){
  require(GenomicRanges)
  require(data.table)
  require(matrixStats)
  merge_seg_gr <- ref_gr
  
  seg_mat <- lapply(names(segl), function(sample_id){
    seg_i <- segl[[sample_id]]
    ov_idx <- findOverlaps(ref_gr, sort(seg_i))
    
    # Find max number of segments that overlap one reference bin
    max_n <- max(table(queryHits(ov_idx)))
    
    # Populate segmean matrix, iterating over multiple intervals per ref_bin
    segmat <- data.frame(matrix(nrow=length(ref_gr), ncol=max_n))
    
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
      segmat[queryHits(non_dup_ov),idx] <- elementMetadata(seg_i)[subjectHits(non_dup_ov),
                                                                   metacol]
      
      idx <- idx+1
    }
    
    # Reduce segmat into one value
    segmat[is.na(segmat)] <- 0
    segmat <- segmat %>%
      as.matrix %>% 
      # apply(., 1, max, na.rm=T) %>%
      # qlcMatrix::rowMax(.) %>%
      # matrixStats::rowMedians(., na.rm=T) %>% 
      matrixStats::rowMeans2(.[i!=0], na.rm = T) %>%
      as.data.frame %>%
      setNames(., sample_id)

    # segmat[is.nan(segmat[,sample_id]),sample_id] <- NA
    return(segmat)
  }) %>% 
    do.call(cbind, .) %>%
    as.data.table
  
  elementMetadata(merge_seg_gr) <- seg_mat
  return(merge_seg_gr)
}

