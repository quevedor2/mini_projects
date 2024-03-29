## Build Transcription Factor BED file
library("ChIPseeker")
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("dplyr")
library("reshape2")
library("GenomicFeatures")
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)

## Set motif of interest on line 66 ##

gbuild <- 'mm10'
human_txdb <- '/cluster/projects/mcgahalab/ref/genomes/human/hg38/GTF/genome.gtf'
mouse_txdb <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
motif_dir <- '/cluster/projects/mcgahalab/ref/tf/JASPAR2020'
motifs <- list("IRF2"=c('Mus musculus'='MA0051.1', 'Homo sapiens'='MA0051.1'),
               "IRF1"=c('Mus musculus'='MA0050.2', 'Homo sapiens'='MA0050.2'),
               'IRF4'=c('Mus musculus'='MA1419.1', 'Homo sapiens'='MA1419.1'),
               'IRF8'=c('Mus musculus'='MA0652.1', 'Homo sapiens'='MA0652.1'),
               'STAT1'=c('Mus musculus'='MA0137.1', 'Homo sapiens'='MA0137.1'))


########################
#### Reference Data ####
# Create a reference map of ENSEMBL to SYMBOL
if(gbuild %in% c('GRCh38', 'GRCh37', 'hg19', 'hg38')){
  txdb <- makeTxDbFromGFF(file = human_txdb, format = "gtf")
  
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
  txdb <- makeTxDbFromGFF(file = mouse_txdb, format = "gtf")
  
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

######################################
#### Get all predicted IRF2 sites ####
idx <- 5 # Set Motif of interest
motif_interest <- setNames(motifs[[idx]][jaspar_species],
                           names(motifs)[idx]) 

## get unique genes
genes <- ChIPseeker:::getGene(txdb, 'gene')
tss <- ifelse(strand(genes) == "+", start(genes), end(genes))
pr <- GRanges(seqnames = seqnames(genes), ranges = IRanges(tss, tss), strand = strand(genes))
pr_idx <- which(duplicated(pr))
genes <- genes[-pr_idx,]


## Get Motif data from JASPAR2020
if(is.null(motif_interest)){
  pfm <- getMatrixSet(x = JASPAR2020,
                      opts = list(collection = "CORE", species=jaspar_species))
} else {
  pfm <- getMatrixSet(x = JASPAR2020,
                      opts = list(collection = "CORE", ID=motif_interest))
  if(length(pfm) == 0){
    pfm <- getMatrixSet(x = JASPAR2020,
                        opts = list(collection = "CORE", name=names(motif_interest)))
    motif_interest <- setNames(names(pfm), names(motifs)[idx])
  }
}

# get promoters
dtss <- 3000
chrs <- GRanges(seqinfo(bsgenome))
chrs <- keepStandardChromosomes(chrs, pruning.mode='coarse')
promoter <- promoters(x=genes(txdb), upstream=dtss, downstream=dtss) %>%  
  sort %>%
  keepStandardChromosomes(., pruning.mode='coarse')
seqlevelsStyle(promoter) <- 'UCSC'

# motif inference from JASPAR2020
motif_ix <- matchMotifs(pwms = pfm,
                        subject = chrs,
                        out = 'positions',
                        bg='genome',
                        p.cutoff = 0.00005, #0.0001,
                        genome = bsgenome)
motif <- sort(motif_ix[[1]])

# Intersect motif with promoter (Not Saved)
ov_idx <- findOverlaps(promoter, motif)
uniq_idx <- which(!duplicated(queryHits(ov_idx)))
promoter$motif <- FALSE
promoter$motif[queryHits(ov_idx)] <- TRUE


# Print out all motifs
motif_df <- as.data.frame(motif)[,-4]
motif_df[,4] <- "."
write.table(motif_df, file=file.path(motif_dir, paste0(names(motif_interest), "_", motif_interest, ".", gbuild, ".genome.bed")),
            col.names = F, row.names = F, sep="\t", quote = F)
motif_gr <- sort(makeGRangesFromDataFrame(motif_df, keep.extra.columns = T))
saveRDS(motif_gr, file=file.path(motif_dir, "rds", paste0(names(motif_interest), "_", motif_interest, ".", gbuild, ".genome.rds")))
