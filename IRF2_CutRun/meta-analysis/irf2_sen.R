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

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/sabelo_irf2/irf2_meta'
outdir <- file.path(PDIR, 'results/sen')
setwd(file.path(PDIR, 'data/sen'))
dir.create(outdir, recursive = T, showWarnings = F)
goi <- c('Tox', 'Pdcd1', 'Batf', 'Irf2')
gbuild <- 'mm9'
txdb <- '/cluster/projects/mcgahalab/ref/genomes/mouse/mm9/GTF/genome.gtf'

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
##############
#### Main ####
txdb <- makeTxDbFromGFF(file = txdb, format = "gtf")

# Create a reference map of ENSEMBL to SYMBOL
gbuild <- 'mm9'
if(gbuild %in% c('GRCh38', 'GRCh37', 'hg19', 'hg38')){
  library("org.Hs.eg.db")
  if(any(gbuild %in% c('hg19', 'GRCh37'))) {
    library(BSgenome.Hsapiens.UCSC.hg37)
    bsgenome <- BSgenome.Hsapiens.UCSC.hg37
  } else if(any(gbuild %in% c('hg38', 'GRCh38'))){
    library(BSgenome.Hsapiens.UCSC.hg38)
    bsgenome <- BSgenome.Hsapiens.UCSC.hg38
  }
  
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
files <- strsplit(paste(list.files(pattern='Mouse.*bed.gz$'), collapse=","), ",")[[1]]
file_list <- data.frame('file'=files,
                        'group'=gsub("^.*Mouse_([a-zA-Z0-9]+)_.*$", "\\1", files))

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
  names(files) <- gsub("^.*Mouse_([a-zA-Z0-9]+)_.*$", "\\1", files_df$file)
  
  ## Categorize the peaks based on their genomic annotation (CEAS style plots)
  peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                         tssRegion=c(-3000, 3000), verbose=FALSE)
  peakAnnoList[[1]]@anno$symbol <- gene_ids[peakAnnoList[[1]]@anno$geneId]
  
  return(peakAnnoList)
})

peaks_goi <- lapply(peaks, function(p){
  dp <- p[[1]]@anno
  symbol_idx <- grep("symbol$", colnames(dp@elementMetadata), value = T)
  dp_goi <- dp[unlist(sapply(paste0("^", goi, "$"), grep, x=dp@elementMetadata[,symbol_idx])),]
  return(dp_goi)
})

lapply(names(peaks_goi), function(id){
  df <- as.data.frame(peaks_goi[[id]])
  colids <- c('seqnames', 'start', 'end', paste0("V", c(4:9)), 'annotation', 'geneId', 'distanceToTSS', 'symbol')
  write.table(df[,colids], file=file.path(outdir, paste0(id, "-targetGenes.tsv")),
              quote = F, col.names = T, row.names = F, sep="\t")
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

motifs <- lapply(peaks, function(p){
  print(names(p))
  # Top scoring motif per position 
  motif_pos <- matchMotifs(pwms = pfm,
                           subject = p[[1]]@anno,
                           out = 'positions',
                           genome = bsgenome)
  names(motif_pos) <- sapply(names(motif_pos), function(i) pfm[[i]]@name)
  motif_pos <- lapply(motif_pos, function(motif){
    ov <- findOverlaps(motif, p[[1]]@anno)
    motif$symbol <- as.character(p[[1]]@anno[subjectHits(ov),]$symbol)
    motif$distanceToTSS <-  as.character(p[[1]]@anno[subjectHits(ov),]$distanceToTSS)
    motif$annotation <-  as.character(p[[1]]@anno[subjectHits(ov),]$annotation)
    motif$geneId <-  as.character(p[[1]]@anno[subjectHits(ov),]$geneId)
    motif
  })
  
  # All possible motifs per location - matrix
  motif_ix <- matchMotifs(pwms = pfm,
                          subject = p[[1]]@anno,
                          out = 'scores',
                          genome = bsgenome)
  motif_scores <- motifScores(object = motif_ix)
  rownames(motif_scores) <- GRangesToString(grange = p[[1]]@anno, sep=c(":", "-"))
  colnames(motif_scores) <- sapply(colnames(motif_scores), function(i) pfm[[i]]@name)
  
  return(list("scores"=motif_scores, "pos"=motif_pos))
})

motif_diff <- as.data.frame(sapply(motifs, function(i) sapply(i$pos, length)))
naive_idx <- grep("Naive", colnames(motif_diff))
motif_diff$delta <- motif_diff[,naive_idx] - rowMeans(motif_diff[,-naive_idx])
motif_diff$delta_frac <- motif_diff$delta / motif_diff[,naive_idx]
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


