library(org.Hs.eg.db)
library(DESeq2)
PDIR <- "/cluster/projects/mcgahalab/data/mcgahalab/wither_mcgaha_ram/snakemake_workflow/results"
setwd(PDIR)
star_cnts <- file.path("deseq2", "all.rds")

#####################
#### STAR Counts ####
star_cnts <- readRDS(star_cnts)
star_cnts <- counts(star_cnts)

# Create a reference map of ENSEMBL to SYMBOL
txby <- keys(org.Hs.eg.db, 'ENSEMBL')
gene_ids <- mapIds(org.Hs.eg.db, keys=txby, column='SYMBOL',
                   keytype='ENSEMBL', multiVals="first")

# Convert Ensembl to Symbol
symbol_ids <- gene_ids[rownames(star_cnts)]
na_idx <- which(is.na(symbol_ids))  # replace NA's with ensembl ids
symbol_ids[na_idx] <- rownames(star_cnts)[na_idx]
symbol_ids <- make.names(symbol_ids, unique=TRUE) # make each duplicate gene unique 

rownames(star_cnts) <- symbol_ids

# Save to a tsv or RDS
write.table(star_cnts, file=file.path("counts", "all.SYMBOL.tsv"), col.names = TRUE,
            row.names = TRUE, quote = FALSE)
