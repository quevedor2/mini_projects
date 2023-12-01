## Sara's re-analysis of Rahul's BMDM vs apBMDM 2016 data
library(cowplot)
library(igraph)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(msigdbr)
library(ggrepel)

# Read in gtf file to map ENS->Biotype
gtf <- '/cluster/projects/mcgahalab/ref/genomes/human/GRCh38/GTF/genome.gtf'
GTF <- rtracklayer::import(gtf)
ens2biotype_ids <- with(GTF, setNames(gene_biotype, gene_id))

# Params
ccle_f <- '/cluster/projects/mcgahalab/ext_data/ccle/CCLE_depmmap_expr.tsv'
ccle_meta_f <- '/cluster/projects/mcgahalab/ext_data/ccle/CCLE_sample_info.csv'

pdir <- "/cluster/projects/mcgahalab/data/mcgahalab/cell_line"
published_dir <- file.path("published", "GSE140026")
data_dedup <- file.path("results", "rsem")
dedir <- file.path(pdir, "results", "diffexp")
outdir <- file.path(pdir, "results", "manual")
dir.create(file.path(outdir, "rds"), showWarnings = F, recursive = F)
setwd(pdir)

# Read in deseq2 data
rsem_dat <- read.table(file = file.path(data_dedup, "MOLM6-merged.genes.results"),
                       sep="\t", header = T, check.names = T, stringsAsFactors = F)


#######################
#### Functions Git ####
source("~/git/mini_projects/mini_functions/wgcnaComplexHeatmap.R")
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")
source("~/git/mini_projects/mini_functions/linearModelEigen.R")
source("~/git/st2_il33/functions/msm_ssm_functions.R")
source("~/git/mini_projects/mini_functions/geneMap.R")

###########################
#### 0. Metadata setup ####
gm <- geneMap()

rsem_dat <- rsem_dat %>% 
  mutate(gene=gm$ENSEMBL$SYMBOL[gene_id],
         entrez=gm$ENSEMBL$ENTREZID[gene_id],
         biotype=ens2biotype_ids[gene_id])
rsem_molm6 <- rsem_dat %>%
  select(TPM, entrez)

ccle_sample_info <-  read.csv(ccle_meta_f, sep=",", header=T, 
                              stringsAsFactors = F, check.names = F)
molm6_id <- ccle_sample_info %>% 
  filter(grepl("MOLM6", CCLE_Name)) %>%
  pull(DepMap_ID)

ccle <- read.table(ccle_f, sep=",", header=T, 
                   stringsAsFactors = F, check.names = F)
colnames(ccle)[1] <- 'DepMap_ID'
ccle_molm6 <- ccle %>% 
  filter(DepMap_ID == molm6_id) %>% 
  tibble::column_to_rownames('DepMap_ID') %>%
  t %>% as.data.frame %>%
  tibble::rownames_to_column("GeneId") %>%
  mutate(#gene=gsub(" \\(.*", "", GeneId),
         entrez=gsub("^.*\\((.*)\\)$", "\\1", GeneId)) %>%
  select(-GeneId)

ccle_rsem <- left_join(ccle_molm6, rsem_molm6, by='entrez') %>% 
  filter(!duplicated(entrez)) %>% 
  tibble::column_to_rownames("entrez") %>%
  mutate(TPM=log2(TPM+1)) %>% 
  rename_with(., ~c('CCLE_MOLM6', 'Local_MOLM6'))

pdf("~/xfer/molm6.pdf")
ggplot(ccle_rsem, aes(x=Local_MOLM6, y=CCLE_MOLM6)) + 
  geom_point()
dev.off()


