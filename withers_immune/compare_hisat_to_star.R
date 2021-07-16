# This code is meant to look at the bulkRNAseq data from the Withers
# mDC, pDC, and Tfh cohort and compare expression between two aligners. 
#  - The STAR counts matrix was processed using the 
#   rna-seq-star-deseq2 workflow: [https://github.com/quevedor2/rna-seq-star-deseq2]
#  - The HISAT counts matrix was processed previously by RAM, the previous
#   bioinformatician. The samples were split across multiple files and needed
#   to be aggregated

library(org.Hs.eg.db)
library(DESeq2)
library(reshape2)

##################
#### Paramers ####
## Aligners path
star_cnts <- '/cluster/projects/mcgahalab/data/mcgahalab/wither_mcgaha_ram/snakemake_workflow/results/deseq2/all.rds'
hisat_dir <- '/cluster/projects/mcgahalab/data/mcgahalab/wither_mcgaha_ram/Merged_read_counts'

#####################
#### STAR Counts ####
star_cnts <- readRDS(star_cnts)
star_cnts <- counts(star_cnts)

######################
#### HISAT Counts ####
# Counts are spread across multiple files
cnts_df <- lapply(grep("read_counts", list.files(hisat_dir), value = TRUE), function(f){
  read.table(file.path(hisat_dir, f), sep=" ", header = TRUE)
})
# Aggregate all the samples into one count matrix and remove duplicate samples
hisat_cnts  <- Reduce(function(x,y) merge(x,y,by='Geneid'), cnts_df)
sample_ids  <- gsub("^.*(mDC|pDC|Tfh)_([a-zA-Z0-9]*).*$", "\\1_\\2", colnames(hisat_cnts))
dup_idx     <- which(duplicated(sample_ids))

hisat_cnts  <- hisat_cnts[,-dup_idx]
sample_ids  <- sample_ids[-dup_idx]
colnames(hisat_cnts) <- gsub("\\.", "_", sample_ids)


#######################################
#### Standardize Gene Names (ROWS) ####
# HISAT matrices have fewer genes are annotated in SYMBOL format
# STAR matrices have more genes are annotated in ENSEMBL format

# Create a reference map of ENSEMBL to SYMBOL
txby <- keys(org.Hs.eg.db, 'ENSEMBL')
gene_ids <- mapIds(org.Hs.eg.db, keys=txby, column='SYMBOL',
                   keytype='ENSEMBL', multiVals="first")

# Convert all HISAT SYMBOLs to ENSEMBL gene IDs
idx <- match(hisat_cnts$Geneid, gene_ids)
hisat_cnts$ensembl   <- names(gene_ids[idx])
hisat_cnts           <- hisat_cnts[-which(is.na(hisat_cnts$ensembl)),]

# Reduce the STAR matrices to just the genes found in the HISAT ones
idx <- match(hisat_cnts$ensembl, rownames(star_cnts))
star_cnts2 <- as.data.frame(star_cnts[idx,])
colnames(star_cnts2) <- gsub("-", "_", colnames(star_cnts2))


###############################
#### Calculate Concordance ####
# Calculate the correlation matrix between HISAT and STAR counts
sample_ids <- grep("_", unique(sort(c(colnames(star_cnts2), colnames(hisat_cnts)))), value=TRUE)
cormat <- sapply(sample_ids, function(sid){
  tryCatch({
    cor(star_cnts2[,sid], hisat_cnts[,sid], use = 'complete.obs')
  }, error=function(e){NA})
})

# Reorganize the correlations into a celltype by sample matrix (e.g. mDC,pDC,Tfh x 367,360,355)
cormat    <- as.matrix(cormat)
cormat    <- as.data.frame(cbind(cormat, do.call(rbind, strsplit(rownames(cormat), "_"))))
cormat$V1 <- round(as.numeric(cormat$V1), 4)
cormat    <- dcast(data=cormat, V2 ~ V3, value.var = "V1")
rownames(cormat) <- cormat[,1]
cormat    <- cormat[,-1]

# Identify samples that are found only in one run of the data
ids <- list("all"=data.frame("sample"=sample_ids),
            "star"=data.frame("sample"=colnames(star_cnts2), "star"=rep(TRUE, ncol(star_cnts2))),
            "hisat"=data.frame("sample"=colnames(hisat_cnts), "hisat"=rep(TRUE, ncol(hisat_cnts))))
ids <- Reduce(function(x,y) merge(x,y,by='sample',all=TRUE), ids)
ids <- ids[which(is.na(rowSums(ids[,-1]))),]
print(ids)
