# An expansion to the ChIPseeker::annotatePeak() function that will annotate
# all genes that overlap a TSS
#
# seqlevelsStyle(cnr_catalogue) <- 'NCBI'
# cnr_catalogue <- annotatePeak(peak=cnr_catalogue, TxDb=txdb, level='gene',
#                               tssRegion=c(-3000, 3000), verbose=FALSE)@anno
# cnr_catalogue <- annotateMultipleTSS(peak=cnr_catalogue, TxDb=txdb, level='gene')
# cnr_catalogue$symbol <- ens2sym_ids[cnr_catalogue$geneId]
# cnr_catalogue$entrez <- ens2entrez_ids[cnr_catalogue$geneId]
# multi_tss_idx <- grep(",", cnr_catalogue$multipleTSS)
# sym_ids <- sapply(cnr_catalogue[multi_tss_idx,]$multipleTSS, function(ens_ids){
#   ens2sym_ids[strsplit(ens_ids, ",")[[1]]] %>% 
#     paste(., collapse=",")
# })
# cnr_catalogue$symbol[multi_tss_idx] <- sym_ids
# seqlevelsStyle(cnr_catalogue) <- 'UCSC'
annotateMultipleTSS <- function(peak, TxDb, level='gene'){
  # Copy annotatePeak() function for select features
  TxDb <- ChIPseeker:::loadTxDb(TxDb)
  if (level == "transcript") {
    features <- ChIPseeker:::getGene(TxDb, by = "transcript")
  }else {
    features <- ChIPseeker:::getGene(TxDb, by = "gene")
  }
  
  # subset for just the TSS start site
  features_tss <- resize(features, width=1) # faster
  
  # Find overlaps between peaks and TSS
  hit <- findOverlaps(peak, unstrand(features_tss)) %>%
    as.data.frame %>% 
    mutate(genes=features[subjectHits]$gene_id)
  
  # Append a comma-separated listed of all genes overlapping each TSS
  hit_qh <- unique(hit$queryHits)
  peak$multipleTSS <- NA
  peak[hit_qh]$multipleTSS <- sapply(split(hit$genes, hit$queryHits), paste, collapse=",")
  
  return(peak)
}


