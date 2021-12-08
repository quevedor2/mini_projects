library(GenomicFeatures)
library(ChIPseeker)
library(clusterProfiler)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)

human_txdb <- '/cluster/projects/mcgahalab/ref/genomes/human/hg38/GTF/genome.gtf'
mouse_txdb <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
pdir <- '/cluster/projects/mcgahalab/data/mcgahalab/sabelo_irf2/irf2_sites'

######################
#### Human ####
txdb <- makeTxDbFromGFF(file = human_txdb, format = "gtf")

# Create a reference map of ENSEMBL to SYMBOL
gbuild <- 'hg38'
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

##############
#### Main ####

# get unique genes
genes <- ChIPseeker:::getGene(txdb, 'gene')
tss <- ifelse(strand(genes) == "+", start(genes), end(genes))
pr <- GRanges(seqnames = seqnames(genes), ranges = IRanges(tss, tss), strand = strand(genes))
pr_idx <- which(duplicated(pr))
genes <- genes[-pr_idx,]

# get PFM
irf2_motif <- 'MA0051.1'
pfm_sel <- getMatrixSet(x = JASPAR2020,
                        opts = list(collection = "CORE", ID=irf2_motif))

# get promoters
dist_tss <- c(50, seq(500, 3000, by=500))
irf2_scores <- sapply(dist_tss, function(dtss){
  # get promoter +/- distance
  promoter <- getPromoters(TxDb=txdb, upstream=dtss, downstream=dtss)

  # motif inference from JASPAR2020
  motif_ix <- matchMotifs(pwms = pfm_sel,
                          subject = promoter,
                          out = 'scores',
                          genome = bsgenome)
  motif_scores <- motifScores(object = motif_ix)
  rownames(motif_scores) <- genes$gene_id
  return(motif_scores)
})
irf2 <- as.data.frame(as.matrix(do.call(cbind, irf2_scores)))
colnames(irf2) <- paste0("tss_D", dist_tss)
irf2$symbol <- gene_ids[rownames(irf2)]
irf2$ensembl <- rownames(irf2)

write.table(irf2, file=file.path(pdir, "irf2_sites.tsv"),
            col.names=T, row.names = F, quote = F, sep="\t")
