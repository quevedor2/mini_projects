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

library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(VennDiagram)
library(DESeq2)
library(msigdbr)

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/sabelo_irf2/irf2_meta'
outdir <- file.path(PDIR, 'results/philip/mouse')
dir.create(outdir, recursive = T, showWarnings = F)
mouse_goi <- c('Tox', 'Pdcd1', 'Batf', 'Irf2')
human_goi <- toupper(mouse_goi)
gbuild <- 'mm9'
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

###############
#### Mouse ####
setwd(file.path(PDIR, 'data/philip/mouse_subset'))
txdb <- makeTxDbFromGFF(file = human_txdb, format = "gtf")

# Create a reference map of ENSEMBL to SYMBOL
gbuild <- 'mm10'
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
                        'group'=gsub("^.*ATAC_(.*?)_.*normalized.*$", "\\1", files))
file_list$group <- gsub("N[1-9]*", "N", file_list$group)

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
  dp_goi <- dp[unlist(sapply(paste0("^", mouse_goi, "$"), grep, x=dp@elementMetadata[,symbol_idx])),]
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

# Get the mean peak height for each group
isolate_peak <- 'promoter'
peaks_heights <- lapply(peaks, function(peak_grp){
  peak_height_grp <- sapply(peak_grp, function(peak_i){
    peak_subset <- as.data.frame(peak_i@anno)
    if(any(grepl("L28_L4", colnames(peak_subset)))){
      colnames(peak_subset)[grep("L28_L4", colnames(peak_subset))] <- 'L28_3'
    }
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

# Run DESeq2 on the peak height matrix
peaks_heights_mat <- as.data.frame(as.matrix(do.call(cbind, peaks_heights_raw)))
peaks_heights_mat <- ceiling(peaks_heights_mat)
peaks_heights_mat <- cbind(data.frame("gene"=rownames(peaks_heights_mat)),
                           peaks_heights_mat)
peaks_meta <- data.frame("id"=colnames(peaks_heights_mat)[-1],
                         "grp"=gsub("_?[1-9]$", "", colnames(peaks_heights_mat)[-1]))

peaks_meta$grp <- factor(as.character(peaks_meta$grp), levels=c('N', 'L7', 'L28'))
peaks_meta <- peaks_meta[order(peaks_meta$grp),]
dds <- DESeqDataSetFromMatrix(countData=peaks_heights_mat, 
                              colData=peaks_meta, 
                              design=~grp, tidy = TRUE)
dds <- DESeq(dds)
resl <- lapply(resultsNames(dds)[-1], function(conditions){
  lfcShrink(dds, coef = conditions, type = 'normal', res = res)
})

# Assemble GeneSets based on motifs and run GSEA using lfc
# motif_gs <- lapply(names(motif$pos), function(id){
motif_gs <- lapply(c(names(motif$pos), 'IRF2'), function(id){
  data.frame("gs_name"=id, "gene"=motif$pos[[id]]$symbol)
})
motif_gs <- as.data.frame(do.call(rbind, motif_gs))
motif_gs <- unique(motif_gs)

# Gene set enrichment to identify up/downregulated motifs between PD1Hi and healthy
# lfc_v <- setNames(as.numeric(unlist(lfc[,1])), rownames(lfc))
gsel <- lapply(resl, function(res){
  ids <- gsub("^.*grp ", "", res@elementMetadata$description[2]) %>% gsub(" ", "_", .)
  
  lfc_v <- setNames(as.numeric(res$log2FoldChange), rownames(res))
  dup_idx <- duplicated(names(lfc_v))
  if(any(dup_idx)) lfc_v <- lfc_v[-which(dup_idx)]
  if(any(is.infinite(lfc_v))) lfc_v[is.infinite(lfc_v)] <- max(lfc_v[-which(is.infinite(lfc_v))])
  gse <- GSEA(sort(na.omit(lfc_v), decreasing = T), TERM2GENE = motif_gs, 
              pvalueCutoff = 1, maxGSSize=15000)
  gse_df <- as.data.frame(gse)[,1:10]
  rownames(gse_df) <- NULL
  
  write.table(gse_df, file=file.path(outdir, paste0("GSEA_motifs.", ids, ".tsv")),
              sep="\t", col.names=T, row.names = F, quote = F)
  return(list("gse"=gse, "df"=gse_df))
})


pdf(file.path(outdir, "gsea_motif.pdf"))
lapply(gsel, function(gse){
  gg_gsea <- gseaplot2(gse, geneSetID = c(head(gse_df$ID, 5), "IRF2"))
  print(gg_gsea)
})
dev.off()

## Identify the genes that are > X LFC in PD1hi to helathy IRF2 genes
irf2 <- split(motif_gs, motif_gs$gs_name)[['IRF2']]
ov_idx <- which(irf2$gene %in% rownames(res))
irf2_lfc <- res[match(irf2$gene[ov_idx], rownames(res)),]
irf2_lfc <- irf2_lfc[order(irf2_lfc[,1]),]
lfc_v <- setNames(irf2_lfc$log2FoldChange, rownames(irf2_lfc))
padj_v <- setNames(irf2_lfc$padj, rownames(irf2_lfc))

irf2_lfc_sig <- irf2_lfc[order(irf2_lfc[,'padj']),]
irf2_lfc_sig <- irf2_lfc_sig[irf2_lfc_sig$padj < 0.15,]
write.table(irf2_lfc_sig, file=file.path(outdir, "sig_irf2_genes.tsv"),
            col.names=T, row.names = T, sep="\t", quote = F)

## GSEA on the IRF2 LFC using msigdb
library("org.Hs.eg.db")
species <- 'Homo sapiens'
genome <- org.Hs.eg.db

txby <- keys(genome, 'SYMBOL')
gene_ids <- mapIds(genome, keys=txby, column='ENTREZID',
                   keytype='SYMBOL', multiVals="first")
names(lfc_v) <- gene_ids[names(lfc_v)]
names(padj_v) <- gene_ids[names(padj_v)]

msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                  'C2'=list('CP:REACTOME'),             # curated gene sets
                  'C5'=list('GO:BP', 'GO:MF'), # ontology gene sets
                  'C7'=list('IMMUNESIGDB'),             # immunologic signature gene sets
                  'C8'=list(NULL))                      # cell type signature gene sets
ora_gsea <- lapply(names(msig_lvls), function(mlvl){
  sub_ora_gra <- lapply(msig_lvls[[mlvl]], function(sublvl){
    print(paste0(">", mlvl, ":", sublvl, "..."))
    msig_ds <- msigdbr(species = species, category = mlvl, subcategory = sublvl) %>% 
      dplyr::select(gs_name, entrez_gene) %>% 
      as.data.frame()
    
    # overrepresentation analysis
    sig_ora <- tryCatch({
      enricher(gene = na.omit(names(padj_v)[padj_v < 0.05]), TERM2GENE = msig_ds)@result
    }, error=function(e){NULL})
    sig_ora$ID <- paste0(c(mlvl, sublvl), collapse="-")
    rownames(sig_ora) <- NULL
    
    # GSEA analysis
    msig_gsea <- tryCatch({
      GSEA(sort(na.omit(lfc_v), decreasing = T), TERM2GENE = msig_ds, pvalueCutoff = 1)
    })
    msig_gsea_df <- as.data.frame(msig_gsea)[,1:10]
    msig_gsea_df$ID <- paste0(c(mlvl, sublvl), collapse="-")
    rownames(msig_gsea_df) <- NULL
    return(list("gsea"=msig_gsea_df[msig_gsea_df$p.adjust<0.15,],
                "ora"=sig_ora[sig_ora$p.adjust < 0.15,-c(8)]))
  })
})
saveRDS(ora_gsea, file=file.path(outdir, "ora_gsea.rds"))

ora_gsea <- unlist(ora_gsea, recursive=F)
gsea <- do.call(rbind, lapply(ora_gsea, function(i) i$gsea))
ora <- do.call(rbind, lapply(ora_gsea, function(i) i$ora))
write.table(gsea, file=file.path(outdir, "irf2-gsea.tsv"),
            col.names=T, row.names = F, quote = F, sep="\t")
write.table(ora, file=file.path(outdir, "irf2-ora.tsv"),
            col.names=T, row.names = F, quote = F, sep="\t")







########################################################################

# Calculate enrichment score for each motif
gene_list = sort(gene_list, decreasing = TRUE)

gseGO
GSEA(gene_list, TERM2GENE = msig_gs)
TERM2GENE


motif_diff <- as.data.frame(sapply(motifs, function(i) sapply(i$pos, length)))
naive_idx <- grep("HN", colnames(motif_diff))
motif_diff$N <- rowMeans(motif_diff[,naive_idx]) 
motif_diff$T <- rowMeans(motif_diff[,-naive_idx]) 
motif_diff$delta_frac <- motif_diff$delta / rowMeans(motif_diff[,naive_idx])
motif_diff <- motif_diff[order(motif_diff$delta_frac),]

#######################################################
#### Differential genes associated with IRF2 peaks ####
library(msigdbr)

# Load in msigdb gene sets
msig_hallmark <- msigdbr(species = jaspar_species, category = "H") %>% 
  dplyr::select(gs_name, entrez_gene) %>% 
  as.data.frame()
msig_immune <- msigdbr(species = jaspar_species, category = "C7") %>% 
  dplyr::select(gs_name, entrez_gene) %>% 
  as.data.frame()

txby <- keys(genome, 'ENSEMBL')
ens2ent <- mapIds(genome, keys=txby, column='ENTREZID',
                   keytype='ENSEMBL', multiVals="first")
txby <- keys(genome, 'ENTREZID')
ent2sym <- mapIds(genome, keys=txby, column='SYMBOL',
                   keytype='ENTREZID', multiVals="first")

irf2_genes <- lapply(motifs, function(m){
  m_gr <- m$pos[[names(irf2_motif)]]
  promoter_idx <- grep("Promoter.*1kb", m_gr$annotation)
  ens2ent[unique(m_gr[promoter_idx,]$geneId)]
})


comp_groups <- list("Acute27-8"=c('AcuteD27', 'AcuteD8'),
                    "Acute8-N"=c('AcuteD8', 'Naive'),
                    "Chronic27-8"=c('ChronicD27', 'Chronic8'),
                    "Chronic8-N"=c('ChronicD8', 'Naive'),
                    'Naive'=c('Naive'))

msig_gs <- msig_hallmark
# msig_gs <- msig_immune
topx <- 10
ggp_enrs <- lapply(list('hallmark'=msig_hallmark, 'immune'=msig_immune), function(msig_gs){
  enrs <- lapply(comp_groups, function(grp){
    print(paste(grp, collapse="-"))
    grp1 <- irf2_genes[[grp[1]]]
    ## Calculate genes in gene-set-1 but not gene-set-2 (Gain) and vice versa (lost)
    if(length(grp)>1) {
      grp2 <- irf2_genes[[grp[2]]]
      gene_set_gained <- setdiff(grp1, grp2)
      gene_set_lost <- setdiff(grp2, grp1)
    } else {
      gene_set_gained <- grp1
      gene_set_lost <- c()
    }
    stat <- if(length(gene_set_lost)>0) TRUE else FALSE
    
    ## Do enrichment analysis using msigdb geneset
    gain_enr <- enricher(gene = gene_set_gained, TERM2GENE = msig_gs)@result
    if(stat) lost_enr <- enricher(gene = gene_set_lost, TERM2GENE = msig_gs)@result
    
    ## Convert the entrezID in msigdb geneset to Symbol
    convert_ent2sym <- function(genes){
      sapply(genes, function(i){
        paste(ent2sym[strsplit(i, split = "\\/")[[1]]], collapse=",")
      })
    }
    gain_enr$geneID <- convert_ent2sym(gain_enr$geneID)
    if(stat) lost_enr$geneID <- convert_ent2sym(lost_enr$geneID)
    
    ## Subset the top X or q<0.15
    gain_enr <- if(topx>0) head(gain_enr, 25) else gain_enr[gain_enr$p.adjust < 0.15,]
    if(stat) lost_enr <- if(topx>0) head(lost_enr, 25) else lost_enr[lost_enr$p.adjust < 0.15,]
    
    gain_enr$stat <- "GAIN"
    if(stat) lost_enr$stat <- "LOST"
    if(!stat) lost_enr <- NULL
    return(list("gained"=gain_enr, "lost"=lost_enr, 
                'gene_gain'=gene_set_gained, 'gene_lost'=gene_set_lost))
  })
  
  ## Bind all the separate dataframes into one single dataframe
  enrs_m <- lapply(names(enrs), function(id){ 
    enr <- enrs[[id]]
    enr_m <- do.call(rbind, enr[c(1:2)])
    enr_m$ID <- id
    enr_m
  })
  enrs_m <- do.call(rbind, enrs_m)
  enrs_m$qvalue[which(enrs_m$qvalue < (1*10^-10))] <- (1*10^-10)
  
  ## Do the plotties
  ggp <- ggplot(enrs_m, aes(x=ID, y=Description, color=stat)) +
    geom_point(aes(size = abs(log10(qvalue)))) + 
    theme_minimal() +
    theme(axis.text.y = element_text(size=7),
          axis.text.x = element_text(size=7, angle=90))
  
  return(list("gg"=ggp, "enr"=enrs))
})

## Plot the enrichment analysis of msigdb for IRF2 accessible genes
pdf(file.path('~/xfer', "immune_enr.pdf"), width = 12, height = 13)
lapply(ggp_enrs, function(i) i$gg)
dev.off()

#############################################
#### Counts of IRF2 genes longitudinally ####
## Plot the change in number of genes for the IRF2 accessible genes for N->D8->D27
library(VennDiagram)
## Chronic LCMV
blank_out <- list(NULL,
                  c('ChronicD27'),
                  c('ChronicD27', 'ChronicD8'))
samples <- c('Naive', 'ChronicD8', 'ChronicD27')
venn_cnts <- getVennCnts(blank_out, irf2_genes, samples)
rownames(venn_cnts) <- c('D27o', 'D8o', 'D27nD8', 'No', 'D27nN', 'D8nN', 'All')
colnames(venn_cnts) <- c('D27', 'D8', 'N')
venn_cnts <- venn_cnts[-match('D27nN', rownames(venn_cnts)),]

# Blank out -only- counts
venn_cnts_bo <- venn_cnts
for(id in colnames(venn_cnts)){
  venn_cnts_bo[grep(paste0(id, "|All"), rownames(venn_cnts_bo), invert=T), id] <- 0
}

ggchronic <- lapply(list("all"=venn_cnts, "bo"=venn_cnts_bo), ggVennCnt, title='Chronic')


## Acute LCMV
blank_out <- list(NULL,
                  c('AcuteD27'),
                  c('AcuteD27', 'AcuteD8'))
samples <- c('Naive', 'AcuteD8', 'AcuteD27')
venn_cnts <- getVennCnts(blank_out, irf2_genes, samples)
rownames(venn_cnts) <- c('D27o', 'D8o', 'D27nD8', 'No', 'D27nN', 'D8nN', 'All')
colnames(venn_cnts) <- c('D27', 'D8', 'N')
venn_cnts <- venn_cnts[-match('D27nN', rownames(venn_cnts)),]

# Blank out -only- counts
venn_cnts_bo <- venn_cnts
for(id in colnames(venn_cnts)){
  venn_cnts_bo[grep(paste0(id, "|All"), rownames(venn_cnts_bo), invert=T), id] <- 0
}

ggacute <- lapply(list("all"=venn_cnts, "bo"=venn_cnts_bo), ggVennCnt, title='Acute')

pdf(file.path("~/xfer", "irf2_gene_counts.pdf"))
ggchronic
ggacute
dev.off()


