library(org.Hs.eg.db)
library(KEGG.db)

library(GO.db)
library(Category)
library(GOstats)

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/wither_mcgaha_ram/snakemake_workflow/results'
setwd(PDIR)
dir.create(file.path("adhoc", "go"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("adhoc", "gsea"), recursive = TRUE, showWarnings = FALSE)

# Create a reference map of ENSEMBL to SYMBOL
txby <- keys(org.Hs.eg.db, 'ENSEMBL')
gene_ids <- mapIds(org.Hs.eg.db, keys=txby, column='ENTREZID',
                   keytype='ENSEMBL', multiVals="first")

## Read gene expr and DEG
infile <- 'ana-vs-control.diffexp.tsv'
gseas <- for(infile in list.files(file.path("diffexp"), pattern=".tsv")){
  res <- read.table(file.path("diffexp", infile), sep="\t", 
                    header=TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(res) <- res$gene
  
  ## Select the significant genes
  min_expr <- 50 # mean(res$baseMean)
  max_p <- 0.05
  resSigind = res[ which(res$padj < max_p & res$log2FoldChange > 0 & res$baseMean > min_expr), ]
  resSigrep = res[ which(res$padj < max_p & res$log2FoldChange < 0 & res$baseMean > min_expr), ]
  resSig = rbind(resSigind, resSigrep)
  resFilt <- res[res$baseMean > min_expr,]
  
  ## GO ontology enrichment analysis
  params=new("GOHyperGParams",
             geneIds=unique(na.omit(gene_ids[rownames(resSig)])),
             universeGeneIds=unique(na.omit(gene_ids[resFilt$gene])),
             annotation="org.Hs.eg.db",
             ontology="BP",
             pvalueCutoff=0.001,
             conditional=TRUE,
             testDirection="over")
  overRepresented=hyperGTest(params)
  go_summ <- summary(overRepresented)[,c(1,2,5,6,7)]
  
  write.table(go_summ, file=file.path("adhoc", "go", gsub("diffexp", "GO", infile)),
              sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  
  ## GSEA
  # Sort the Log2FC genes
  gene_list <- resFilt$log2FoldChange
  names(gene_list) <- resFilt$gene
  gene_list<-sort(na.omit(gene_list), decreasing = TRUE)
  
  gse <- gseGO(geneList=gene_list, 
               ont ="ALL", 
               keyType = "ENSEMBL", 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Hs.eg.db, 
               pAdjustMethod = "fdr")
  
  pdf(file.path("adhoc", "gsea", gsub("diffexp.*", "gsea.pdf", infile)), height = 20, width = 20)
  print(emapplot(pairwise_termsim(gse), showCategory = 50))
  print(ridgeplot(gse) + labs(x = "enrichment distribution"))
  #gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
  dev.off()
  
  write.table(gse@result[,1:11], file=file.path("adhoc", "gsea", gsub("diffexp.*", "gsea.tsv", infile)),
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  
  
  ## KEGG GSEA
  kegg_gene_list <- resFilt$log2FoldChange
  names(kegg_gene_list) <- gene_ids[resFilt$gene]
  kegg_gene_list <- kegg_gene_list[-which(is.na(names(kegg_gene_list)))]
  kegg_gene_list<-sort(na.omit(kegg_gene_list), decreasing = TRUE)
  
  
  
  kk2 <- gseKEGG(geneList     = kegg_gene_list,
                 organism     = 'hsa',
                 nPerm        = 10000,
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "fdr",
                 keyType       = "ncbi-geneid",
                 use_internal_data=TRUE)
  
  pdf(file.path("adhoc", "gsea", gsub("diffexp.*", "gsea-kegg.pdf", infile)), height = 20, width = 20)
  print(dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + 
    facet_grid(.~.sign))
  print(emapplot(pairwise_termsim(kk2), showCategory = 50))
  print(ridgeplot(kk2) + labs(x = "enrichment distribution"))
  dev.off()
  
  write.table(kk2@result[,1:11], file=file.path("adhoc", "gsea", gsub("diffexp.*", "gsea-kegg.tsv", infile)),
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  
  #return(list("go"=gse, "kegg"=kk2))
}
