library(reshape2)
library(ggplot2)
library(dplyr)
wideScreen()
setwd("/cluster/projects/mcgahalab/data/mcgahalab/xin_GCN2/results/rsem")
kogene <- 'Tsc1'

#########################################
#### Differential isoform expression ####
isoform_stat <- read.csv("isoform_status.csv", header=F) %>% 
  tibble::column_to_rownames(., "V2")

id <- 'SIGNR1_LysGCN2_UN_3_S49-merged.isoforms.results'
id.df <- read.table(id, sep="\t", header=T, nrows = 2)
x <- read.table(paste0(kogene, ".isoforms.txt"), sep="\t", header=F) %>%
  mutate(sample=gsub("-merged.*", "", V1),
         V1=gsub("^.*\\:", "", V1))
colnames(x) <- c(colnames(id.df), 'sample')

xl <- split(x[,c('sample', 'transcript_id', 'TPM')], x$sample) 
x.df <- lapply(xl, function(x_i){
  colnames(x_i)[3] <- unique(x_i$sample)
  x_i[,-1]
}) %>%
  purrr::reduce(., full_join, by='transcript_id') %>%
  tibble::column_to_rownames(., "transcript_id") %>%
  select(grep("CD169|SIGNR1", colnames(.), value=T))

# Calculate the fractional composition of transcripts per sample
xfrac.df <- apply(x.df, 2, function(i) i/sum(i))


# Store metadata
xmeta <- data.frame(celltype=gsub("_.*$", "\\1", colnames(x.df)),
                    knockout=gsub("^.*?_(Lys[a-zA-Z0-9]*|WT)_.*$", "\\1", colnames(x.df)),
                    activated=gsub("^.*_(AC|UN)_.*", "\\1", colnames(x.df)),
                    id=colnames(x.df)) %>%
  tibble::column_to_rownames(., 'id')

# Calculate the fraction of functional and non-functional isoforms
xfunc.l <- split(as.data.frame((xfrac.df)), f=isoform_stat[rownames(xfrac.df),]$V3)
xfunc_melt <- sapply(xfunc.l, function(xfunc_i){
  i <- apply(xfunc_i, 2, sum)
}) %>% 
  cbind(xmeta, .) %>% 
  tibble::rownames_to_column(., "sample") %>%
  melt

# Test the difference between fractional composition per transcript
# Splikt by celltype, then activation status, then test for difference between 
# knockout condition and wt
xfrac.ct <- split(as.data.frame(t(xfrac.df)), 
                  f=xmeta$celltype)
trx_frac_dist <- lapply(xfrac.ct, function(xfrac.ct_i){
  xmeta.ct <- xmeta[rownames(xfrac.ct_i),]
  xfrac.ct.act <- lapply(split(xfrac.ct_i, f=xmeta.ct$activated), function(x_i){
    xmeta.ct.act <- xmeta[rownames(x_i),]
    x.ko_i <- split(x_i, f=xmeta.ct.act$knockout)
    sapply(x.ko_i, function(i){
      colMeans(i) - colMeans(x.ko_i$WT)
    }) %>% 
      as.data.frame %>% 
      mutate(activation=unique(xmeta.ct.act$activated))
  }) %>%
    do.call(rbind, .) %>%
    as.data.frame %>%
    mutate(celltype=unique(xmeta.ct$celltype))
  xfrac.ct.act
}) %>%
  do.call(rbind, .) %>%   as.data.frame %>%
  tibble::rownames_to_column(., 'transcript') %>%
  mutate(transcript=gsub("^.*\\.", "", transcript))

trxdist.melt <- melt(trx_frac_dist) 

x.melt <- melt(t(xfrac.df)) %>%
  mutate(group=gsub("^.*?_?(Lys[A-Za-z0-9]*|WT)_.*$", "\\1", Var1, perl = T),
         Var2=factor(Var2, levels=names(rowSums(xfrac.df) %>% sort)),
         celltype=gsub("_.*", "", Var1))

pdf(paste0("~/xfer/xin_", kogene, ".pdf"), width = 9)
ggplot(trxdist.melt, aes(x=value, y=transcript, fill=variable)) +
  geom_bar(stat='identity', position='dodge') +
  facet_grid(variable ~ celltype+activation, scales='free', space='free') +
  theme_classic() + 
  ggtitle(kogene) +
  xlim(-0.4, 0.4)

ggplot(xfunc_melt, aes(x=sample, y=value, fill=variable, group=variable)) +
  geom_bar(stat='identity', position='stack') +
  facet_grid(knockout ~ celltype+activated, scales='free', space='free') +
  theme_classic() 

ggplot(x.melt, aes(x=Var1, y=value, fill=Var2)) +
  facet_grid(celltype ~ group, scales='free', space='free') +
  geom_bar(stat='identity', position='stack') +
  theme_classic()
dev.off()

# library(biomaRt)
# GENES = useMart("ENSEMBL_MART_ENSEMBL", 
#                 dataset = "mmusculus_gene_ensembl")   
# 
# x <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol',
#                           'ensembl_transcript_id','refseq_mrna',
#                           'ucsc','chromosome_name',
#                           'transcript_start','transcript_end'), 
#            mart = GENES)
# x[which(x$refseq_mrna %in% c("NM_001306081", "NM_028898", "XM_006534363")),]

##########################################################
#### Differential expression between individual exons ####
# Used bedtools to count the coverage across each exon for the genes we're 
# interested in (RAPTOR, GCN2, TSC1), quantile-normalized across all samples,
# and then use a simple one-sided t.test or wilcox.test to test for significance

library(dplyr)
library(reshape2)
library(ggplot2)
library(preprocessCore)
library(ggrepel)
dir='/cluster/projects/mcgahalab/data/mcgahalab/xin_GCN2/results/star/pe/subset/exon_counts'
setwd(dir)

gene_map <- c('ENSMUSG00000005102'='Eif2ak4',
              'ENSMUSG00000025583'='Rptor', 
              'ENSMUSG00000026812'='Tsc1')

exon_cnt <- lapply(list.files(pattern="(CD169|SIGNR1).*.bed"), function(f){
  id <- gsub("-merged.*", "", f)
  df <- read.table(f, sep="\t", header=F) %>%
    mutate(gene=gsub("^.*?gene_id ([a-zA-Z0-9]*?);.*", "\\1", V4),
           transcript=gsub("^.*?transcript_id ([a-zA-Z0-9]*?);.*", "\\1", V4),
           exon=gsub("^.*?exon_number ([a-zA-Z0-9]*?);.*", "\\1", V4),
           pos=paste0(V1, ":", V2, "-", V3)) %>%
    dplyr::select(-c(V4, V1, V2, V3)) %>% 
    rename_with(., ~c(id, 'gene', 'transcript', 'exon', "pos"))
  return(df)
}) %>% purrr::reduce(., full_join, by=c('gene', 'transcript', 'exon', "pos"))
exon_cnt_meta <- exon_cnt[,c('gene', 'transcript', 'exon', "pos")] %>%
  mutate(symbol=gene_map[gene])
exon_cnt <- exon_cnt %>% dplyr::select(-c('gene', 'transcript', 'exon', "pos"))
exon_norm <- normalize.quantiles(as.matrix(exon_cnt)) %>%
  as.data.frame %>%
  rename_with(., ~colnames(exon_cnt))

celltype <- gsub("^(.*?)_.*$", "\\1", colnames(exon_norm))
exon_norm_celltypes <- split(as.data.frame(t(exon_norm)), f=celltype)
exon_qvals <- lapply(exon_norm_celltypes, function(exon_ct_j){
  
  ids <- gsub("^.*?(Lys[a-zA-Z0-9]*?|WT)_.*$", "\\1", rownames(exon_ct_j))
  exon_norm_spl <- split(exon_ct_j, f=ids)
  
  pmat <- lapply(exon_norm_spl[1:3], function(exon_i){
    pvals <- sapply(seq_along(colnames(exon_i)), function(exon_idx){
      tryCatch({
        res <- t.test(exon_i[,exon_idx], 
                      exon_norm_spl$WT[,exon_idx], 
                      alternative='less')
        c('t'=res$statistic, 'p'=res$p.value)
      }, error=function(e){c('t'=0, 'p'=1)})
    }) %>%
      t %>% as.data.frame %>% 
      rename_with(., ~c('t', 'p')) %>%
      mutate(q=p.adjust(p, method='fdr'))
    return(pvals)
    # return(p.adjust(pvals, method='fdr'))
  }) %>% do.call(cbind, .) %>% 
    cbind(exon_cnt_meta, .) %>%
    as.data.frame %>%
    mutate(gene_exon=paste0(symbol, "_", exon))
  
  pmat
})

threshold <- 0.15
dat_melt <- lapply(names(exon_qvals), function(celltype){
  dat <- exon_qvals[[celltype]]
  dat_melt_i <- lapply(c('LysGCN2', 'LysRAPTOR', 'LysTSC'), function(ko){
    print(ko)
    dat[,c('gene_exon', grep(ko, colnames(dat),value=T))] %>% 
      rename_with(., ~c('gene_exon', 't', 'p', 'q')) %>%
      mutate(knockout=ko,
             sig=ifelse(q < 0.1, paste0('sig (q<', threshold, ')'), 'nonsig'),
             q=(-1*log10(q)),
             label=t <= max(head(sort(t), 10)))
  }) %>% 
    do.call(rbind, .) %>%
    as.data.frame %>%
    mutate(celltype=celltype)
  
  return(dat_melt_i)
}) %>% 
  do.call(rbind, .)


pdf("~/xfer/volcano.pdf", width=10, height = 9)
ggplot(dat_melt, aes(x=t, y=q, color=sig)) +
  geom_point() +
  geom_text_repel(data=subset(dat_melt, label), #q > (-1*log10(threshold))),
                  aes(t,q,label=gene_exon), size=3,
                  box.padding = 0.5, max.overlaps = Inf) +
  facet_grid(celltype ~ knockout, space='free') +
  ylim(0,4) + xlim(-10,10) +
  geom_vline(xintercept=0, lty=2) +
  theme_classic()
dev.off()

celltype_i <- 'SIGNR1'
exon_qvals[[celltype_i]] %>% arrange(LysGCN2) %>% head(.,20) filter(LysGCN2 < threshold) 
exon_qvals[[celltype_i]] %>% arrange(LysRAPTOR.q) %>% head(.,20) filter(LysRAPTOR < threshold)
exon_qvals[[celltype_i]] %>% arrange(LysTSC) %>% head(.,20) filter(LysTSC < threshold)