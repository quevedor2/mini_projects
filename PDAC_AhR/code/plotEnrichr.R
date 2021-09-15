library(dplyr)
library(ggplot2)
setwd("~/Projects/mcgaha/pdac_ahr/enrichr")

hubmap <- read.table("HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression_table.txt", 
                     header=TRUE, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)
hga <- read.table("Human_Gene_Atlas_table.txt", header = TRUE,
                  sep="\t", stringsAsFactors = FALSE, check.names = FALSE)

df <- hubmap

plotBarplot <- function(df){
  colnames(df) <- gsub(" ", "_", colnames(df)) %>% gsub("-", "_", .)
  df <- df[order(df$Combined_Score),]
  df$Term <- gsub(" \\(id:.*", "", as.character(df$Term))
  df$Term <- factor(as.character(df$Term), as.character(df$Term))
  df$Sig <- df$P_value > 0.05
  
  ggp <- ggplot(df, aes(y=Combined_Score, x=Term, fill=Sig)) +
    geom_bar(stat='identity', show.legend = FALSE) + 
    theme_minimal() + 
    scale_fill_manual("legend", 
                      values = c("FALSE" = "#d6604d", 
                                 "TRUE" = "#878787")) +
    coord_flip() +
    theme(axis.title.y=element_blank(),
          axis.text = element_text(size=rel(0.5)),
          axis.title.x = element_text(size=rel(0.7)))
  
  return(ggp)
}

pdf("enrichr_barplots.pdf", height = 2.5)
plotBarplot(hubmap)
plotBarplot(hga)
dev.off()
