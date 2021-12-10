library("ChIPseeker")
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("dplyr")
library("reshape2")
library("GenomicFeatures")

library("KEGG.db")
library("GO.db")
library("Category")
library("GOstats")

library(DESeq2)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(VennDiagram)
library("msigdbr")


PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/sabelo_irf2/irf2_meta'
outdir <- file.path(PDIR, 'results/philip')
dir.create(outdir, recursive = T, showWarnings = F)
mouse_goi <- c('Tox', 'Pdcd1', 'Batf', 'Irf2')
human_goi <- toupper(mouse_goi)
gbuild <- 'hg38'
human_txdb <- '/cluster/projects/mcgahalab/ref/genomes/human/hg38/GTF/genome.gtf'
mouse_txdb <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'

###################
#### Functions ####
GRangesToString <- function(grange, sep = c("-", "-")) {
  regions <- paste0(
    as.character(x = seqnames(x = grange)),
    sep[[1]],
    start(x = grange),
    sep[[2]],
    end(x = grange)
  )
  return(regions)
}

getVennCnts <- function(blank_out, irf2_genes, samples){
  venn_cnts <- lapply(blank_out, function(bo){
    irf <- irf2_genes[samples]
    if(!is.null(bo)){
      for(b in bo){
        irf[[b]] <- ''
      }
    }
    venn <- get.venn.partitions(irf)
    venn[,c('..set..', '..count..')]
  })
  venn_cnts <- Reduce(function(x,y) merge(x,y,by='..set..', all=TRUE), venn_cnts)
  rownames(venn_cnts) <- venn_cnts[,1]
  venn_cnts <- venn_cnts[,-1]
  return(venn_cnts)
}

ggVennCnt <- function(vc, title){
  melt_venn <- melt(t(vc))
  melt_venn$Var1 <- factor(melt_venn$Var1, levels=c('N', 'D8', 'D27'))
  melt_venn$Var2 <- factor(melt_venn$Var2, 
                           levels=rev(c('No', 'D8nN', 'All', 'D8o', 'D27nD8', 'D27o')))
  
  ggplot(melt_venn, aes(x=as.integer(Var1), y=value, fill=Var2)) + 
    geom_area() +
    scale_x_continuous(breaks=c(1,2,3), labels=levels(melt_venn$Var1)) +
    ylab("Number of IRF2-linked Genes") + xlab("") + ggtitle(title) +
    theme_minimal()
}

getMotif <- function(p, pfm, bsgenome){
  # Top scoring motif per position 
  motif_pos <- matchMotifs(pwms = pfm,
                           subject = p@anno,
                           out = 'positions',
                           genome = bsgenome)
  names(motif_pos) <- sapply(names(motif_pos), function(i) pfm[[i]]@name)
  motif_pos <- lapply(motif_pos, function(motif){
    ov <- findOverlaps(motif, p@anno)
    motif$symbol <- as.character(p@anno[subjectHits(ov),]$symbol)
    motif$distanceToTSS <-  as.character(p@anno[subjectHits(ov),]$distanceToTSS)
    motif$annotation <-  as.character(p@anno[subjectHits(ov),]$annotation)
    motif$geneId <-  as.character(p@anno[subjectHits(ov),]$geneId)
    motif
  })
  
  # All possible motifs per location - matrix
  motif_ix <- matchMotifs(pwms = pfm,
                          subject = p@anno,
                          out = 'scores',
                          genome = bsgenome)
  motif_scores <- motifScores(object = motif_ix)
  rownames(motif_scores) <- GRangesToString(grange = p@anno, sep=c(":", "-"))
  colnames(motif_scores) <- sapply(colnames(motif_scores), function(i) pfm[[i]]@name)
  
  return(list("scores"=motif_scores, "pos"=motif_pos))
}

######################
#### Human ####
setwd(file.path(PDIR, 'data/philip/human_subset'))
txdb <- makeTxDbFromGFF(file = human_txdb, format = "gtf")

# Create a reference map of ENSEMBL to SYMBOL
if(gbuild %in% c('GRCh38', 'GRCh37', 'hg19', 'hg38')){
  if(any(gbuild %in% c('hg19', 'GRCh37'))) {
    library(BSgenome.Hsapiens.UCSC.hg37)
    bsgenome <- BSgenome.Hsapiens.UCSC.hg37
  } else if(any(gbuild %in% c('hg38', 'GRCh38'))){
    library(BSgenome.Hsapiens.UCSC.hg38)
    bsgenome <- BSgenome.Hsapiens.UCSC.hg38
  }
  
  library("org.Hs.eg.db")
  species <- 'org.Hs.eg.db'
  jaspar_species <- 'Homo sapiens'
  genome <- org.Hs.eg.db
} else if(gbuild %in% c('mm10', 'mm9', 'GRCm38')){
  if(gbuild %in% c('mm9')) {
    library(BSgenome.Mmusculus.UCSC.mm9)
    bsgenome <- BSgenome.Mmusculus.UCSC.mm9
  } else if(any(gbuild %in% c('mm10', 'GRCm38'))){
    library(BSgenome.Mmusculus.UCSC.mm10)
    bsgenome <- BSgenome.Mmusculus.UCSC.mm10
  }
  
  library("org.Mm.eg.db")
  species <- 'org.Mm.eg.db'
  jaspar_species <- 'Mus musculus'
  genome <- org.Mm.eg.db
}
txby <- keys(genome, 'ENSEMBL')
gene_ids <- mapIds(genome, keys=txby, column='SYMBOL',
                   keytype='ENSEMBL', multiVals="first")

## Read a file containing all the files and groupings
files <- list.files(pattern='normalizedCounts.txt.gz$')
file_list <- data.frame('file'=files,
                        'group'=gsub("^.*ATAC_(.+)_[1-9]_normalized.*$", "\\1", files))

#########################
#### Peak Annotation ####
## Read in all the peak files based on their group ID
file_list_split <- split(file_list, file_list$group)
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

## Annotate the peaks with their closest gene
peaks <- lapply(file_list_split, function(files_df){
  file_sizes <- unlist(sapply(files_df$file, file.info)['size',])
  if(any(file_sizes == 0)) files_df <- files_df[-which(file_sizes<1000),]
  files <- as.list(files_df$file)
  names(files) <- gsub("^.*ATAC_(.+)_normalized.*$", "\\1", files_df$file)
  
  ## Categorize the peaks based on their genomic annotation (CEAS style plots)
  peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, level = "gene",
                         tssRegion=c(-3000, 3000), verbose=FALSE)
  for(idx in seq_along(peakAnnoList)){
    peakAnnoList[[idx]]@anno$symbol2 <- gene_ids[peakAnnoList[[idx]]@anno$geneId]
  }
  
  return(peakAnnoList)
})

peaks_goi <- lapply(unlist(peaks, recursive=F), function(p){
  dp <- p@anno
  symbol_idx <- grep("symbol$", colnames(dp@elementMetadata), value = T)
  dp_goi <- dp[unlist(sapply(paste0("^", human_goi, "$"), grep, x=dp@elementMetadata[,symbol_idx])),]
  return(dp_goi)
})

x <- lapply(names(peaks_goi), function(id){
  if(!file.exists(file.path(outdir, paste0(id, "-targetGenes.tsv")))){
    df <- as.data.frame(peaks_goi[[id]])
    df <- do.call(rbind, lapply(split(df, df$symbol), function(x) x[order(abs(x$distanceToTSS)),]))
    
    colids <- c('seqnames', 'start', 'end', 'refseqID', colnames(df)[9], 'peak_annotation',
                'annotation', 'geneId', 'distanceToTSS', 'symbol', 'symbol2')
    write.table(df[,colids], file=file.path(outdir, paste0(id, "-targetGenes.tsv")),
                quote = F, col.names = T, row.names = F, sep="\t")
  }
})


#########################
#### Motif inference ####
# Add motif information from JASPAR2020
irf2_motif <- 'MA0051.1'
pfm <- getMatrixSet(x = JASPAR2020,
                    opts = list(collection = "CORE", species=jaspar_species))
pfm_sel <- getMatrixSet(x = JASPAR2020,
                        opts = list(collection = "CORE", ID=irf2_motif))
pfm <- c(pfm, pfm_sel)
irf2_motif <- setNames(irf2_motif, sapply(pfm_sel, function(i) i@name))

# Since Philip et al. aggregated on conserved peaks, we can just run one motif search
motif <- getMotif(peaks[[1]][[1]], pfm, bsgenome)

###############################################
#### DESeq2 differential abundance testing ####
# Get the mean peak height for each group
isolate_peak <- 'promoter'
peaks_heights <- lapply(peaks, function(peak_grp){
  peak_height_grp <- sapply(peak_grp, function(peak_i){
    peak_subset <- as.data.frame(peak_i@anno)
    peak_idx <- which(colnames(peak_subset) %in% names(peak_grp))
    if(!is.null(isolate_peak)){
      peak_isolate <- split(peak_subset, peak_subset$peak_annotation)[[isolate_peak]]
      peak_height <- peak_isolate[,peak_idx]
      peak_height <- setNames(peak_height, (peak_isolate$symbol))
    } else {
      peak_height <- peak_subset[,peak_idx]
      peak_height <- setNames(peak_height, (peak_subset$symbol))
    }
    return(peak_height)
  })
  
  return(list("raw"=peak_height_grp, "mean"=rowMeans(peak_height_grp)))
})
peaks_heights_raw <- lapply(peaks_heights, function(i) i$raw)
peaks_heights_mean <- sapply(peaks_heights, function(i) i$mean)

## Assemble a matrix to feed into DESeq
## and run DESeq2 on the peak height matrix
peaks_heights_mat <- as.data.frame(as.matrix(do.call(cbind, peaks_heights_raw)))
peaks_heights_mat <- ceiling(peaks_heights_mat)
peaks_heights_mat <- cbind(data.frame("gene"=rownames(peaks_heights_mat)),
                           peaks_heights_mat)
peaks_meta <- data.frame("id"=colnames(peaks_heights_mat)[-1],
                         "grp"=gsub("_[1-9]$", "", colnames(peaks_heights_mat)[-1]))

# DESeq2 for differential abundance testing
dds <- DESeqDataSetFromMatrix(countData=peaks_heights_mat, 
                              colData=peaks_meta, 
                              design=~grp, tidy = TRUE)
dds <- DESeq(dds)
res <- results(dds)
res <- lfcShrink(dds, coef="grp_Hum_PD1hi_vs_HN", type="apeglm")

# output the DESEq2 results and volcano plot
resdf <- as.data.frame(res)
resdf$gene <- rownames(resdf)
resdf$log10_padj <- -1*log10(resdf$padj)
resdf$sig <- resdf$padj < 0.001 & abs(resdf$log2FoldChange) > 2
resdf$sig <- setNames(c('q<0.001 & log2LFC>2', ''), c('TRUE', 'FALSE'))[as.character(resdf$sig)]
resdf_sig <- resdf[which(nchar(resdf$sig)>0),]

write.table(resdf, file=file.path(outdir, "DESeq2_dag.tsv"),
            row.names=T, col.names=T, sep="\t", quote=F)

gp <- ggplot(data=resdf, aes(x=log2FoldChange, y=log10_padj, color=sig)) +
  geom_point(alpha = 0.6) + 
  theme_bw() +
  xlim(-10,10) +
  geom_text_repel(data=resdf_sig, aes(label=gene))
pdf(file.path(outdir, "volcano_plot.pdf"))
gp
dev.off()


########################################################
#### Gene Set enrichment of the DAGs: JASPAR Motifs ####
# Assemble GeneSets based on motifs and run GSEA using lfc
# motif_gs <- lapply(names(motif$pos), function(id){
motif_gs <- lapply(c(names(motif$pos), 'IRF2'), function(id){
  data.frame("gs_name"=id, "gene"=motif$pos[[id]]$symbol)
})
motif_gs <- as.data.frame(do.call(rbind, motif_gs))
motif_gs <- unique(motif_gs)

# Gene set enrichment to identify up/downregulated motifs between PD1Hi and healthy
# lfc_v <- setNames(as.numeric(unlist(lfc[,1])), rownames(lfc))
lfc_v <- setNames(as.numeric(res$log2FoldChange), rownames(res))
dup_idx <- duplicated(names(lfc_v))
if(any(dup_idx)) lfc_v <- lfc_v[-which(dup_idx)]
if(any(is.infinite(lfc_v))) lfc_v[is.infinite(lfc_v)] <- max(lfc_v[-which(is.infinite(lfc_v))])
gse <- GSEA(sort(na.omit(lfc_v), decreasing = T), TERM2GENE = motif_gs, 
     pvalueCutoff = 1, maxGSSize=15000)
gse_df <- as.data.frame(gse)[,1:10]
rownames(gse_df) <- NULL

write.table(gse_df[,-c(2)], file=file.path(outdir, "GSEA_jaspar_motifs.tsv"),
            sep="\t", col.names=T, row.names = F, quote = F)
# gse_df[gse_df$p.adjust < 0.05,]

pdf(file.path(outdir, "GSEA_jaspar_motifs.pdf"), width = 12)
gse %>%
  enrichplot:::barplot.enrichResult(x='NES', showCategory=40) + 
  theme(text = element_text(size=8),
        axis.text.y = element_text(size=8))  +
  xlim(-2,2) + xlab("NES") +
  ggtitle("Top JASPAR motifs")

  maxval <- max(abs(lfc_v))
  colors <- c('navyblue', 'lightgrey', 'darkred')
  b <- c(-ceiling(maxval), 0, ceiling(maxval))
  gse %>%
    filter(grepl("IRF2", ID)) %>%
    heatplot(showCategory=20, foldChange=lfc_v) +
    theme_bw() +
    scale_fill_gradientn(limits = c(min(b),max(b)), colors = colors, 
                         breaks = b, labels = format(b)) +
    theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1))
dev.off()


################################################################
#### Gene Set enrichment of the IRF2 DAGs: mSigDb Gene Sets ####
## Identify the genes that are > X LFC in PD1hi to helathy IRF2 genes
irf2 <- split(motif_gs, motif_gs$gs_name)[['IRF2']]
ov_idx <- which(irf2$gene %in% rownames(res))
irf2_lfc <- res[match(irf2$gene[ov_idx], rownames(res)),]
irf2_lfc <- irf2_lfc[order(irf2_lfc[,1]),]
lfc_v <- setNames(irf2_lfc$log2FoldChange, rownames(irf2_lfc))
padj_v <- setNames(irf2_lfc$padj, rownames(irf2_lfc))

irf2_lfc_sig <- irf2_lfc[order(irf2_lfc[,'padj']),]
write.table(irf2_lfc_sig, file=file.path(outdir, "sig_irf2_genes.tsv"),
            col.names=T, row.names = T, sep="\t", quote = F)
# irf2_lfc_sig <- irf2_lfc_sig[irf2_lfc_sig$padj < 0.15,]

## GSEA on the IRF2 LFC using msigdb
txby <- keys(genome, 'SYMBOL')
gene_ids <- mapIds(genome, keys=txby, column='ENTREZID',
                   keytype='SYMBOL', multiVals="first")
names(lfc_v) <- gene_ids[names(lfc_v)]
names(padj_v) <- gene_ids[names(padj_v)]

msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                  'C2'=list('CP:REACTOME'),             # curated gene sets
                  'C5'=list('GO:BP', 'GO:MF'), # ontology gene sets
                  'C7'=list('IMMUNESIGDB'))             # immunologic signature gene sets
## Over representation analysis
ora_gsea <- lapply(names(msig_lvls), function(mlvl){
  sub_ora_gra <- lapply(msig_lvls[[mlvl]], function(sublvl){
    print(paste0(">", mlvl, ":", sublvl, "..."))
    msig_ds <- msigdbr(species = jaspar_species, category = mlvl, subcategory = sublvl) %>% 
      dplyr::select(gs_name, entrez_gene) %>% 
      as.data.frame()
    
    # overrepresentation analysis
    sig_ora <- tryCatch({
      enricher(gene = na.omit(names(padj_v)[padj_v < 0.05]), 
               TERM2GENE = msig_ds, maxGSSize=2000)@result
    }, error=function(e){NULL})
    sig_ora$ID <- paste0(c(mlvl, sublvl), collapse="-")
    rownames(sig_ora) <- NULL
    
    # GSEA analysis
    msig_gsea <- tryCatch({
      GSEA(sort(na.omit(lfc_v), decreasing = T), TERM2GENE = msig_ds, 
           pvalueCutoff = 1, maxGSSize = 2000)
    }, error=function(e){NULL})
    msig_gsea_df <- as.data.frame(msig_gsea)[,1:10]
    msig_gsea_df$ID <- paste0(c(mlvl, sublvl), collapse="-")
    rownames(msig_gsea_df) <- NULL
    return(list("gsea"=msig_gsea,
                "gsea_df"=msig_gsea_df,
                "ora"=sig_ora,
                "lfc"=sort(na.omit(lfc_v), decreasing = T)))
  })
})
saveRDS(ora_gsea, file=file.path(outdir, "ora_gsea.rds"))

ora_gsea <- unlist(ora_gsea, recursive=F)
names(ora_gsea) <- unlist(sapply(names(msig_lvls), function(i) paste(i, unlist(msig_lvls[[i]]), sep="_")))
gsea_df <- do.call(rbind, lapply(ora_gsea, function(i) i$gsea_df))
ora <- do.call(rbind, lapply(ora_gsea, function(i) i$ora))
write.table(gsea_df, file=file.path(outdir, "irf2-gsea.tsv"),
            col.names=T, row.names = F, quote = F, sep="\t")
write.table(ora, file=file.path(outdir, "irf2-ora.tsv"),
            col.names=T, row.names = F, quote = F, sep="\t")



pdf(file.path(outdir, "GSEA_irf2_msigdb.pdf"), width = 11)
lapply(names(ora_gsea), function(id){
  gse <- ora_gsea[[id]]$gsea
  
  ggp <- gse %>%
    enrichplot:::barplot.enrichResult(x='NES', showCategory=40) + 
    theme(text = element_text(size=8),
          axis.text.y = element_text(size=8))  +
    xlim(-4,4) + xlab("NES") +
    scale_fill_gradient(trans = "log", low = "light green",  high = "black", 
                        breaks = c(0, 0.001, 0.05, 1)) +
    ggtitle(paste0("IRF2: ", id))
  print(ggp)
  NULL
})
dev.off()

pdf(file.path(outdir, "GSEA_irf2_msigdb_heatmap.pdf"), width = 20, height = 5)
lapply(names(ora_gsea), function(id){
  gse <- ora_gsea[[id]]$gsea
  lfc_v <- ora_gsea[[id]]$lfc
  gsex <- setReadable(gse, 'org.Hs.eg.db', 'ENTREZID')
  maxval <- max(abs(lfc_v))
  colors <- c('navyblue', 'lightgrey', 'darkred')
  b <- c(-ceiling(maxval), 0, ceiling(maxval))
  ggp <- tryCatch({
    gsex %>%
    filter(p.adjust < 0.15) %>%
    heatplot(showCategory=10, foldChange=lfc_v) +
    theme_bw() + 
      ggtitle(paste0("Top 10 IRF2: ", id, " (padj < 0.15)")) +
    scale_fill_gradientn(limits = c(min(b),max(b)), colors = colors, 
                         breaks = b, labels = format(b)) +
    theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1))
  }, error=function(e){NULL})
  print(ggp)
  NULL
})
dev.off()




