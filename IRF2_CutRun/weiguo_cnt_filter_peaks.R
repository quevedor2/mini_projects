library(lsa)
library(GenomicRanges)
library(dplyr)
library(reshape2)
library(ggplot2)

library(ChIPseeker)
library(org.Mm.eg.db)
library(GenomicFeatures)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(msigdbr)


pdir <- '/cluster/projects/mcgahalab/data/mcgahalab/sabelo_irf2/weiguo_cutandtag'
gtf_file <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
gbuild <- 'GRCm38'
species <- 'org.Mm.eg.db'
genome <- org.Mm.eg.db

ext_range <- 100  # BP to add to a peak to extend the peak search overlap
quantile_vals <- c(seq(0.1, 0.8, by=0.1), 
                   seq(0.81, 0.9, by=0.01),
                   seq(0.905, 1, by=0.005))

irf2_bed <- file.path(pdir, "ref", "IRF2_MA0051.1.mm10.genome.bed")
peaks_dir <- file.path(pdir, "results", "peaks", "seacr_single")

###################
#### Functions ####
# Reads in the peak bed files and formats it into GR objects
readPeak <- function(pathto, ext_range=0){
  grpeak <- setNames(read.table(pathto, sep="\t"),
                     c("chr", "start_raw", "end_raw", "AUC", "max", "range")) %>%
    mutate(start=start_raw - ext_range,
           end=end_raw + ext_range) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = T)
  seqlevelsStyle(grpeak) <- 'UCSC'
  grpeak
}

# Populate the target-control signal matrix
populateCatalogue <- function(cnr_catalogue, ov_idx, grirf2, gr_target, gr_ctrl){
  targ_ctrl_df <- matrix(nrow=length(cnr_catalogue), ncol=3) %>%
    as.data.frame() %>%
    setNames(., c("target", "control", "irf2"))
  targ_ctrl_df[queryHits(ov_idx),]$irf2 <- grirf2[subjectHits(ov_idx)]$sig
  
  ov <- findOverlaps(cnr_catalogue, gr_target)
  targ_ctrl_df[queryHits(ov),1] <- round(gr_target[subjectHits(ov)]$AUC,2)
  ov <- findOverlaps(cnr_catalogue, gr_ctrl)
  targ_ctrl_df[queryHits(ov),2] <- round(gr_ctrl[subjectHits(ov)]$AUC,2)
  return(targ_ctrl_df)
}

##################
#### 0. Setup ####
# Read in IRF2 motif matrix
grirf2 <- setNames(read.table(irf2_bed, sep="\t", check.names = F),
                   c("chr", "start", "end", "strand", "sig")) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = T)

# Generate Cut&Run peak catalogue
cnr_catalogue <- sapply(list.files(peaks_dir, pattern="bed"), function(pfile){
  readPeak(file.path(peaks_dir, pfile),
           ext_range=ext_range)
})
cnr_catalogue <- as(cnr_catalogue, "GRangesList") %>%
  unlist() %>%
  reduce()

# Create a mapping between catalogue and IRF2 motifs
irf2_cat_ov <- findOverlaps(cnr_catalogue, grirf2)

# Create listing of TARGET - CONTROL pair of files
f_pairs <- list.files(peaks_dir, pattern="bed")
f_grps <- gsub("_[a-zA-Z0-9]*\\..*", "", f_pairs)
f_pairs <- split(f_pairs, f_grps)
f_pairs <- lapply(f_pairs, function(fp){
  data.frame("target"=grep("igg", fp, ignore.case = T, invert = T, value = T)) %>%
    mutate(control=grep("igg", fp, ignore.case = T, value=T))
}) %>%
  do.call(rbind, .)
rownames(f_pairs) <- gsub(".stringent.*", "", f_pairs$target)

# Import GTF file
txdb <- makeTxDbFromGFF(file = gtf_file, format = "gtf")
txby <- keys(genome, 'ENSEMBL')
ens2sym_ids <- mapIds(genome, keys=txby, column='SYMBOL',
                      keytype='ENSEMBL', multiVals="first")

###################################
#### 1. IRF2 Overlapping Peaks ####
dist <- sapply(list.files(peaks_dir, pattern="bed"), function(pfile){
  grpeak <- readPeak(file.path(peaks_dir, pfile),
                     ext_range=ext_range)
  
  quantiles <- quantile(grpeak$AUC, quantile_vals)
  sapply(quantiles, function(q){
    grpeak_q <- grpeak[grpeak$AUC > q,]
    ov_idx <- findOverlaps(grpeak_q, grirf2)
    length(unique(queryHits(ov_idx))) / length(grpeak_q)
  })
})

dist2 <- round(dist, 3)
colnames(dist2) <- gsub(".stringent.*", "", colnames(dist2))

mdist2 <- melt(dist2) %>%
  rename_with(., ~ c("signalAUC_percentile", "sample", "IRF2_proportion")) %>%
  mutate(signalAUC_percentile = gsub("%", "", signalAUC_percentile),
         group = gsub("^.*_", "", sample),
         marker = gsub("(.*_.*)_.*", "\\1", sample))
mdist2$signalAUC_percentile <- as.numeric(as.character(mdist2$signalAUC_percentile))
pdf("~/xfer/irf2_motif_peak.all.pdf", width = 10)
ggplot(mdist2, aes(x=signalAUC_percentile, IRF2_proportion, color=marker, group=sample)) +
  geom_line() + 
  facet_grid(rows=vars(group)) +
  theme_classic() +
  theme(axis.text.x=element_text(size=6, angle=45))
dev.off()

############################################
#### 2. Cosine similarity between peaks ####
# Fill in the Signal matrix for all instances of the Peaks catalogue
signal_mat <- sapply(list.files(peaks_dir, pattern="bed"), function(pfile){
  gr0 <- readPeak(file.path(peaks_dir, pfile), ext_range = 0)
  
  setNames(read.table(file.path(peaks_dir, pfile), sep="\t"),
                  c("chr", "start", "end", "AUC", "max", "range")) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = T)
  signal_vec <- matrix(nrow=length(cnr_catalogue), ncol=1)
  
  ov <- findOverlaps(cnr_catalogue, gr0)
  signal_vec[queryHits(ov),] <- round(gr0[subjectHits(ov)]$AUC,2)
  return(signal_vec)
})

# Generate a cosine similarity matrix between the catalogue
jaccard <- function(a, b) {
  sum(!is.na(a) & !is.na(b)) / sum(!is.na(a) | !is.na(b))
}

metric <- 'cosine'
signal_mat[is.na(signal_mat)] <- 0
colnames(signal_mat) <- gsub(".stringent.*", "", colnames(signal_mat))
sim_mats <- lapply(quantile_vals, function(min_percentile){
  sim_mat <- apply(signal_mat, 2, function(i){
    i[which(i < quantile(i[i!=0], min_percentile))] <- switch(metric, 
                                                              jaccard=NA,
                                                              cosine=0)
    apply(signal_mat, 2, function(j){
      j[which(j < quantile(j[j!=0], min_percentile))] <- switch(metric, 
                                                                jaccard=NA,
                                                                cosine=0)
      if(metric == 'cosine'){
        lsa::cosine(i,j)
      } else {
        jaccard(i,j)
      }
    })
  })
  return(sim_mat)
})
names(sim_mats) <- quantile_vals
  
melt_sims <- melt(sim_mats) %>%
  mutate(L1 = gsub("%", "", L1),
         target_j = gsub("^.*_", "", Var2),
         sample = gsub("(.*_.*)_.*", "\\1", Var2)) %>%
  rename_with(., ~ c("sample_i", "sample_j", "cosine_similarity", 
                     "percentile", "target", "sample")) %>% 
  filter(., as.numeric(cosine_similarity) < 0.99)
melt_sims$percentile <- as.numeric(as.character(melt_sims$percentile)) 
pdf("~/xfer/test.pdf", height = 12)
ggplot(melt_sims, aes(x=percentile, y=cosine_similarity, 
                      group=sample_j, col=target, linetype=sample)) +
  facet_grid(rows=vars(sample_i), scales="free") +
  ylim(0,0.5) +
  geom_line()
dev.off()

#######################################################
#### 3. Target/Control overlap and Motif breakdown ####
# cnr_catalogue   # peak-catalogue across all samples
# grirf2          # interval catalogue for IRF2 motifs
irf2_props <- apply(f_pairs, 1, function(pfiles){
  gr_target <- readPeak(file.path(peaks_dir, pfiles['target']), ext_range = 0)
  gr_ctrl <- readPeak(file.path(peaks_dir, pfiles['control']), ext_range = 0)
  
  # Populate the target-control signal matrix
  targ_ctrl_df <- populateCatalogue(cnr_catalogue, irf2_cat_ov, grirf2, gr_target, gr_ctrl)
  
  ## Create percentile-cutoffs for each sample to look at only peaks exceeding
  # a set signal value
  quantiles <- apply(targ_ctrl_df[,c('target', 'control')], 2, quantile, 
                     probs=quantile_vals[-length(quantile_vals)], na.rm=T)
  irf2_prop <- apply(quantiles, 1, function(qs){
    # Remove peaks that don't meet the threshold
    targ_ctrl_mat <- as.matrix(targ_ctrl_df[,c(1,2)])
    thresh_mat <- matrix(qs, byrow = T, 
                         nrow=nrow(targ_ctrl_mat), ncol=ncol(targ_ctrl_mat))
    thresh_mat_lowCtrl <- thresh_mat
    targ_ctrl_mat[targ_ctrl_mat < thresh_mat] <- NA
    
    # Keep low signal control peaks to compare against high-signal target peaks
    targ_ctrl_mat_lowCtrl <-  as.matrix(targ_ctrl_df[,c(1,2)])
    thresh_mat_lowCtrl[,2] <- 0.1
    targ_ctrl_mat_lowCtrl[targ_ctrl_mat_lowCtrl < thresh_mat_lowCtrl] <- NA
    
    # Subset for peaks that are found only in target, or target+control
    keep_idx <- rowSums(is.na(targ_ctrl_mat))
    keep_idx <- setNames(split(c(1:nrow(targ_ctrl_mat)), f=keep_idx),
                         c('bothHigh', 'one', 'none'))
    
    # Subset for peaks that are found both in thresholded-target and no-threshold control
    keep_idx_lowCtrl <- rowSums(is.na(targ_ctrl_mat_lowCtrl))
    keep_idx_lowCtrl <- setNames(split(c(1:nrow(targ_ctrl_mat)), f=keep_idx_lowCtrl),
                                 c('bothLow', 'one', 'none'))['bothLow']
    
    # Subset for peaks that are found in target-only or control-only
    ss_keep_idx <- is.na(targ_ctrl_mat[keep_idx[['one']],1])
    ss_keep_idx <- setNames(split(keep_idx[['one']], ss_keep_idx),
                            c('trgPk', 'ctrlPk'))
    keep_idx2 <- c(keep_idx['bothHigh'], keep_idx_lowCtrl['bothLow'],
                   ss_keep_idx[c('trgPk', 'ctrlPk')])
    
    sapply(keep_idx2, function(kidx){
      sum(!is.na(targ_ctrl_df$irf2[kidx])) / length(kidx)
    })
  })
  
  melt_irf2_prop <- irf2_prop %>%
    as.data.frame() %>%
    rename_with(., ~ gsub("%", "", colnames(irf2_prop))) %>%
    t() %>% melt() %>%
    setNames(., c('percentile', 'category', 'irf2_proportion'))
  melt_irf2_prop$percentile <- as.numeric(as.character(melt_irf2_prop$percentile))
  gg <- ggplot(melt_irf2_prop, aes(x=percentile, y=irf2_proportion, group=category, col=category)) + 
    geom_line() + 
    ggtitle(paste(gsub(".stringent.*", "", pfiles), collapse=" - ")) + 
    theme_classic() + ylim(0,0.25)
  
  return(list("gg"=gg, "mat"=irf2_prop))
})
pdf("~/xfer/target_control_irf2.pdf", height = 10)
plot_grid(plotlist=lapply(irf2_props[c(1,3,4,2,5)], function(i) i$gg), ncol=1)
dev.off()

#####################################################################
#### 4. Target/Control identifying thresholds for IRF2 detection ####
# cnr_catalogue   # Sec2: peak-catalogue across all samples
# grirf2          # Sec1: interval catalogue for IRF2 motifs
# irf2_cat_ov     # Sec3: Overlap between Peak-catalogue and IRF2 motifs
# f_pairs         # Sec3: list of all target-control sample pairs
irf2_props_mat <- apply(f_pairs, 1, function(pfiles){
  gr_target <- readPeak(file.path(peaks_dir, pfiles['target']), ext_range = 0)
  gr_ctrl <- readPeak(file.path(peaks_dir, pfiles['control']), ext_range = 0)
  
  # Populate the target-control signal matrix
  targ_ctrl_df <- populateCatalogue(cnr_catalogue, irf2_cat_ov, grirf2, gr_target, gr_ctrl)
  
  # Create index of peaks where both Target and Control overlap
  tc_ov <- as.matrix(targ_ctrl_df[,1:2]) %>%
    is.na() %>%
    rowSums()
  tc_ov_idx <- split(c(1:nrow(targ_ctrl_df)), f=tc_ov)[['0']]
  
  ## Create percentile-cutoffs for each sample to look at only peaks exceeding
  # a set signal value
  quantiles <- apply(targ_ctrl_df[,c('target', 'control')], 2, quantile, 
                     probs=quantile_vals[-length(quantile_vals)], na.rm=T) %>%
    as.data.frame()
  irf2_ov <- targ_ctrl_df$irf2[tc_ov_idx]
  prop_mat <- sapply(setNames(quantiles$target, rownames(quantiles)), function(q_t){
    targ <- targ_ctrl_df$target[tc_ov_idx]
    targ[targ < q_t] <- NA
    sapply(setNames(quantiles$control, rownames(quantiles)), function(q_c){
      ctrl <- targ_ctrl_df$control[tc_ov_idx]
      ctrl[ctrl < q_c] <- NA
      both_idx <- which(!is.na(targ) & !is.na(ctrl))
      sum(!is.na(irf2_ov[both_idx])) / length(both_idx)
    })
  })
  return(list("matrix"=prop_mat, "quantiles"=quantiles))
})
names(irf2_props_mat) <- apply(f_pairs, 1, paste, collapse=" - ") %>%
  gsub(".stringent[a-zA-Z0-9.]*", "", .)

pdf("~/xfer/irf2_prop_matrix.pdf")
b <- c(0.1, 0.15, 0.2)
lapply(names(irf2_props_mat), function(id){
  ggplot(melt(irf2_props_mat[[id]]), aes(Var1, Var2)) +
    # geom_raster(aes(fill = value), interpolate = FALSE) +
    geom_tile(aes(fill=value)) +
    ggtitle(id) + ylab("Signal percentile threshold (target)") + 
    xlab("Signal percentile threshold (control)") +
    theme(axis.text.x = element_text(angle=45, size=6),
          axis.text.y = element_text(size=6)) +
    scale_fill_gradientn(colours=c("black", "grey", "yellow", "red"))
})
dev.off()

#####################################
#### 5. Refining confident peaks ####
thresholds <- setNames(c(0.985, 0.985, 0.915), 
                       paste0("Act_unfixed_", c("IRF2", "H3K4me3", "IgG"), ".stringent.bed"))
grs_anno <- apply(f_pairs[c('Act_unfixed_IRF2', 'Act_unfixed_H3K4me3'),], 1, function(pfiles){
  gr_target <- readPeak(file.path(peaks_dir, pfiles['target']), ext_range = 0)
  gr_ctrl <- readPeak(file.path(peaks_dir, pfiles['control']), ext_range = 0)
  targ_thresh <- thresholds[pfiles['target']]
  ctrl_thresh <- thresholds[pfiles['control']]
  
  # Populate the target-control signal matrix
  targ_ctrl_df <- populateCatalogue(cnr_catalogue, irf2_cat_ov, grirf2, gr_target, gr_ctrl)
  
  # Subset for peaks that are found only in target, or target+control
  targ_ctrl_mat <- as.matrix(targ_ctrl_df[,c(1,2)])
  keep_idx <- rowSums(is.na(targ_ctrl_mat))
  keep_idx <- setNames(split(c(1:nrow(targ_ctrl_mat)), f=keep_idx),
                       c('both', 'one', 'none'))
  
  # Remove peaks that don't meet the threshold
  thresh_mat <- matrix(c(quantile(targ_ctrl_df$target, targ_thresh, na.rm=T),
                         quantile(targ_ctrl_df$control, ctrl_thresh, na.rm=T)), 
                       byrow = T, nrow=nrow(targ_ctrl_mat), ncol=ncol(targ_ctrl_mat))
  targ_ctrl_mat[targ_ctrl_mat[,1] < thresh_mat[,1],1] <- NA
  targ_ctrl_mat[targ_ctrl_mat[,2] >= thresh_mat[,2],2] <- NA

  # Identify peaks that are found  in targetOnly
  idx <- list('both'=intersect(which(rowSums(!is.na(targ_ctrl_mat)) == 2), keep_idx[['both']]), 
              'trgPk'=intersect(which(!is.na(targ_ctrl_mat[,1])), keep_idx[['one']]) )
  
  ## Annotate the peaks with their gene and genomic annotation
  grs <- lapply(names(idx), function(idx_id){
    idx_i <- idx[[idx_id]]
    gr0 <- cnr_catalogue[idx_i]
    gr0@elementMetadata <- as(setNames(targ_ctrl_df[idx_i,],
                                       c('targetSignal', 'ctrlSignal', 'IRF2')), 
                              "DataFrame")
    gr0$classification <- idx_id
    seqlevelsStyle(gr0) <- 'NCBI'
    grl <- setNames(split(gr0, is.na(gr0$IRF2)),
                    c('IRF2pos', 'IRF2neg'))
    
    gr0_anno <- lapply(grl, annotatePeak, TxDb=txdb, level='gene',
                       tssRegion=c(-3000, 3000), verbose=FALSE)
    return(gr0_anno)
  })
  names(grs) <- names(idx)
  
  return(grs)
})

# Writing the annotated peaks to a file
lapply(names(grs_anno), function(id){
  grl <- unlist(grs_anno[[id]], recursive=F)
  lapply(names(grl), function(id_i){
    as.data.frame(grl[[id_i]]@anno) %>%
      mutate(annotation=gsub(",", "-", annotation),
             symbol=ens2sym_ids[geneId]) %>%
      write.table(., file=paste0("~/xfer/", id, ".", gsub("\\.", "_", id_i), ".csv"),
                  sep=",", row.names = F, col.names = T, quote = F)
  })
})


