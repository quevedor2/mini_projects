# AutomatedClusterMarkerTable returns FindAllMarkers table with extra bits of useful information
# and an educated guess about cluster identity

AutomatedClusterMarkerTable <- function(Seurat_Object, ClusterList, cell_markers=NULL){
  if(is.null(cell_markers)){
    cell_markers <- list("Epithelial"=c("KRT7","KRT8","KRT18","KRT19","EPCAM","CDH1"),
                         "T-cell"=c("CD3E","CD3G","CD3D","CD4","IL7R","CD8A","LEF1"),
                         "Myeloid"=c("CD14","ITGAM","MNDA","MPEG1","ITGAX"),
                         "B-cell"=c("CD79A","MS4A1","CD19"),
                         "Fibroblasts"=c("CDH11","PDGFRA","PDGFRB","ACTA2"),
                         "RBC"=c("HBA1","HBB","HBA2"),
                         "NK"=c("NCR3","FCGR3A","NCAM1","KLRF1","KLRC1","CD38","KLRC1"),
                         "Endothelial"=c("CDH5","PECAM1"),
                         "Acinar"=c("TRY4","SPINK1","AMY2A"))
  }
  
  library(dplyr)
  library(pracma)
  library(tibble)
  library(Seurat)
  # ClusterList <- list()
  # Idents(object = Seurat_Object) <- "seurat_clusters"
  current.cluster.ids <- sort(as.numeric(levels(Seurat_Object@active.ident)))
  new.cluster.ids <- c()
  
  
  for(i in current.cluster.ids){
    List_Position <- i + 1
    # ClusterList[[List_Position]] <- FindMarkers(object = Seurat_Object, ident.1 = i, min.pct = 0.25, only.pos = TRUE)
    Positive_Genes <- rownames(ClusterList[[List_Position]])
    Num_Positive_Genes <- length(Positive_Genes)
    
    RPS_Num <- length(grep(pattern = "^RPS", x = Positive_Genes))
    RPL_Num <- length(grep(pattern = "^RPL", x = Positive_Genes))
    RP_Percent <- sum(RPS_Num, RPL_Num)/length(Positive_Genes)*100
    RP_Label <- paste("RP%:", RP_Percent, sep = " ")
    
    Mito_Num <- length(grep(pattern = "^MT-", x = Positive_Genes))
    Mito_Percent <- Mito_Num/length(Positive_Genes)*100
    Mito_Label <- paste("Mito%:", RP_Percent, sep = " ")
    
    ClusterCells <- WhichCells(object = Seurat_Object, idents = i)
    Cell_Barcodes <- unlist(Seurat_Object@assays$RNA@counts@Dimnames[2])
    Cell_Number <- match(ClusterCells, Cell_Barcodes)
    # 
    # for(k in 1:length(ClusterCells)){
    #   Cell_Position <- grep(pattern = ClusterCells[k], x = Cell_Barcodes, value = FALSE)
    #   Cell_Number <- c(Cell_Number,Cell_Position)
    # }
    
    
    S_Score <- Seurat_Object@meta.data$S.Score
    G2M_Score <- Seurat_Object@meta.data$G2M.Score
    Cluster_S_Score <- S_Score[Cell_Number]
    Cluster_G2M_Score <- G2M_Score[Cell_Number]
    Avg_Cluster_S_Score <- mean(Cluster_S_Score)
    Avg_Cluster_G2M_Score <- mean(Cluster_G2M_Score)
    Cluster_S_Score_Range <- range(Cluster_S_Score)
    Cluster_G2M_Score_Range <- range(Cluster_G2M_Score)
    
    nFeature <- Seurat_Object@meta.data$nFeature_RNA
    nCount <- Seurat_Object@meta.data$nCount_RNA
    Mito <- Seurat_Object@meta.data$percent.mt
    
    Cluster_nFeature <- nFeature[Cell_Number]
    Cluster_nCount <- nCount[Cell_Number]
    Cluster_Mito <- Mito[Cell_Number]
    
    Avg_Cluster_nFeature <- as.integer(mean(Cluster_nFeature))
    Avg_Cluster_nCount <- as.integer(mean(Cluster_nCount))
    Max_Cluster_Mito <- max(Cluster_Mito)
    
    # Find least common multiple based on length of all cell marker types
    lcmbool <- TRUE; lcmval=0
    while(lcmbool){
      lcmval <- lcmval + 1
      lcms <- pracma::Lcm(lcmval, sapply(cell_markers, length))
      if(length(unique(lcms))==1) lcmbool <- FALSE
    }
    Score_Weights <- lcmval / sapply(cell_markers, length)
    clust_thresh <- (lcmval/2) + 0.5
    
    All_Scores <- as.list(setNames(rep(0, length(cell_markers)), names(cell_markers)))
    # Cell_Types <- c("Epi","T Cell","Myeloid","B Cell","Fibroblast","RBC","NK", "Endo","Acinar")
    # Epi_Markers <- c("KRT7","KRT8","KRT18","KRT19","EPCAM","CDH1")
    # T_Cell_Markers <- c("CD3E","CD3G","CD3D","CD4","IL7R","CD8A","LEF1")
    # Myeloid_Markers <- c("CD14","ITGAM","MNDA","MPEG1","ITGAX")
    # B_Cell_Markers <- c("CD79A","MS4A1","CD19")
    # Fibroblast_Markers <- c("CDH11","PDGFRA","PDGFRB","ACTA2")
    # RBC_Markers <- c("HBA1","HBB","HBA2")
    # NK_Markers <- c("NCR3","FCGR3A","NCAM1","KLRF1","KLRC1","CD38","KLRC1")
    # Endo_Markers <- c("CDH5","PECAM1")
    # Acinar_Markers <- c("TRY4","SPINK1","AMY2A")
    # All_Markers <- list(Epi_Markers,T_Cell_Markers,Myeloid_Markers,B_Cell_Markers,Fibroblast_Markers,RBC_Markers,NK_Markers,Endo_Markers,Acinar_Markers)
    # 
    # Epi_Score <- 0
    # T_Cell_Score <- 0
    # Myeloid_Score <- 0
    # B_Cell_Score <- 0
    # Fibroblast_Score <- 0
    # RBC_Score <- 0
    # NK_Score <- 0
    # Endo_Score <- 0
    # Acinar_Score <- 0 
    # All_Scores <- list(Epi_Score,T_Cell_Score,Myeloid_Score,B_Cell_Score,Fibroblast_Score,RBC_Score,NK_Score,Endo_Score,Acinar_Score)
    # Score_Weights <- c(1.85,1.85,2.22,3.7,2.78,3.7,1.85,5.56,3.7) 
    
    
    Weighted_Scores <- c()
    for(h in 1:length(cell_markers)){
      Markers_to_Test<- cell_markers[[h]]
      Marker_Row <- h
      for(j in 1:length(Markers_to_Test)){
        Gene_Found <- 0
        Gene_Found <- length(grep(pattern = Markers_to_Test[j], x = Positive_Genes))
        if(Gene_Found > 0 ){
          All_Scores[[Marker_Row]] <- All_Scores[[Marker_Row]]+1
        }
      }
      Weighted_Scores[Marker_Row] <- All_Scores[[Marker_Row]]*Score_Weights[Marker_Row]
    }
    
    ClusterID <- which(Weighted_Scores >= clust_thresh)
    if(length(ClusterID) > 0){
      if(length(ClusterID) > 1){
        ID <- "Multiple"
      }else{
        ID <- names(cell_markers)[ClusterID]
      }
    }else{
      ID <- i
    } 
    if(RP_Percent > 30){
      ID <- paste("RP_",ID,sep = "")
    }
    if(Avg_Cluster_S_Score > 0.01 | Avg_Cluster_G2M_Score > 0.01){
      CellCycleID <- "Cycling"
      ID <- paste("Cycling_",ID,sep = "")
    }else{
      CellCycleID <- "N/A"
    }
    if(Avg_Cluster_nCount < 700){
      ID <- paste("G_",ID,sep = "")
    }
    new.cluster.ids <- c(new.cluster.ids,ID)
    
    # Label_Row <- length(Positive_Genes) + 1
    # Label_Row2 <- length(Positive_Genes) + 2
    # Label_Row3 <- length(Positive_Genes) + 3
    # Label_Row4 <- length(Positive_Genes) + 4
    # Label_Row5 <- length(Positive_Genes) + 5
    # 
    # Label1 <- c("Summary:",paste("Cluster",i, sep = " "), paste("ID:",ID, sep = " "),paste("Mito%:",Mito_Percent, sep = " "),paste("RP%:",RP_Percent, sep = " "))
    # Label2 <- c("Immune Summary",paste("T Cell Score:",All_Scores[[2]], sep = " "),paste("Myeloid Score:",All_Scores[[3]], sep = " "),paste("B Cell Score:",All_Scores[[4]], sep = " "),
    #             paste("NK Score:",All_Scores[[7]], sep = " "))
    # Label3 <- c(paste("Epi Score:",All_Scores[[1]], sep = " "),paste("Fib Score:",All_Scores[[5]], sep = " "),paste("Acinar Score:",All_Scores[[9]], sep = " "),
    #             paste("Endo Score:",All_Scores[[8]], sep = " "),paste("RBC Score:",All_Scores[[6]], sep = " "))
    # Label4 <- c("Avg S Score:", Avg_Cluster_S_Score, "Avg G2M Score:", Avg_Cluster_G2M_Score, CellCycleID)
    # Label5 <- c("Filter Info",paste("Avg. nGene:", Avg_Cluster_nFeature, sep = " "),paste("Avg. nCounts:", Avg_Cluster_nCount, sep = " ")
    #             , paste("Highest Mito:",Max_Cluster_Mito, sep = " "), paste("# Cells:",length(Cell_Number), sep = " "))
    # 
    # ClusterList[[List_Position]][Label_Row,] <- c(Label1, rep(NA,2))
    # ClusterList[[List_Position]][Label_Row2,] <- c(Label2, rep(NA,2))
    # ClusterList[[List_Position]][Label_Row3,] <- c(Label3, rep(NA,2))
    # ClusterList[[List_Position]][Label_Row4,] <- c(Label4, rep(NA,2))
    # ClusterList[[List_Position]][Label_Row5,] <- c(Label5, rep(NA,2))
    # 
    # ClusterList[[List_Position]] <- rownames_to_column(.data = ClusterList[[List_Position]],var = "Gene")
    # ClusterList[[List_Position]][Label_Row,"Gene"] <- "Summary1"
    # ClusterList[[List_Position]][Label_Row2,"Gene"] <- "Summary2"
    # ClusterList[[List_Position]][Label_Row3,"Gene"] <- "Summary3"
    # ClusterList[[List_Position]][Label_Row4,"Gene"] <- "Summary4"
    # ClusterList[[List_Position]][Label_Row5,"Gene"] <- "Summary5"
    # 
    
  }
  
  ClusterDataFrame <- bind_rows(ClusterList, .id = "column_label")
  ClusterDataFrame <- ClusterDataFrame[,-1]
  ClusterPackage <- list(ClusterDataFrame, new.cluster.ids)
  return(ClusterPackage)
}
