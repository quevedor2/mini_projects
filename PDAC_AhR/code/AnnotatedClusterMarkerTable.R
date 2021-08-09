# AutomatedClusterMarkerTable returns FindAllMarkers table with extra bits of useful information
# and an educated guess about cluster identity

AutomatedClusterMarkerTable <- function(Seurat_Object, ClusterList, cell_markers=NULL,
                                        exact.match=TRUE, cyc=TRUE, rp=TRUE){
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
  
  new.cluster.ids <- sapply(current.cluster.ids, assignClusterId,
                            ClusterList=ClusterList, Seurat_Object=Seurat_Object, 
                            cell_markers=cell_markers, exact.match=exact.match,
                            cyc=cyc, rp=rp)

  ClusterDataFrame <- bind_rows(ClusterList, .id = "column_label")
  ClusterDataFrame <- ClusterDataFrame[,-1]
  ClusterPackage <- list(ClusterDataFrame, new.cluster.ids)
  return(ClusterPackage)
}


assignClusterId <- function(i, ClusterList, Seurat_Object, cell_markers, 
                            exact.match=FALSE, cyc=TRUE, rp=TRUE){
  List_Position <- as.integer(as.character(i)) + 1
  # ClusterList[[List_Position]] <- FindMarkers(object = Seurat_Object, ident.1 = i, min.pct = 0.25, only.pos = TRUE)
  Positive_Genes <- rownames(ClusterList[[List_Position]])
  Num_Positive_Genes <- length(Positive_Genes)
  
  if(rp){
    # Count number of differentially expressed ribosomal genes
    RPS_Num <- length(grep(pattern = "^RPS", x = Positive_Genes))
    RPL_Num <- length(grep(pattern = "^RPL", x = Positive_Genes))
    RP_Percent <- sum(RPS_Num, RPL_Num)/length(Positive_Genes)*100
    RP_Label <- paste("RP%:", RP_Percent, sep = " ")
  } else {
    RP_Percent <- 0
  }
  
  if(cyc){
    # Assign Cell barcodes to cluster IDs
    ClusterCells <- WhichCells(object = Seurat_Object, idents = i)
    Cell_Barcodes <- unlist(Seurat_Object@assays$RNA@counts@Dimnames[2])
    Cell_Number <- match(ClusterCells, Cell_Barcodes)
    
    # Identify cells that are currently cell cycling
    S_Score <- Seurat_Object@meta.data$S.Score
    G2M_Score <- Seurat_Object@meta.data$G2M.Score
    Cluster_S_Score <- S_Score[Cell_Number]
    Cluster_G2M_Score <- G2M_Score[Cell_Number]
    Avg_Cluster_S_Score <- mean(Cluster_S_Score)
    Avg_Cluster_G2M_Score <- mean(Cluster_G2M_Score)
    Cluster_S_Score_Range <- range(Cluster_S_Score)
    Cluster_G2M_Score_Range <- range(Cluster_G2M_Score)
  } else {
    Avg_Cluster_S_Score <- Avg_Cluster_G2M_Score <- 0
  }
  
  # Find least common multiple based on length of all cell marker types
  # in order to set a common threshold across all cell marker types
  lcmbool <- TRUE; lcmval=0
  while(lcmbool){
    lcmval <- lcmval + 1
    lcms <- pracma::Lcm(lcmval, sapply(cell_markers, length))
    if(length(unique(lcms))==1) lcmbool <- FALSE
  }
  Score_Weights <- lcmval / sapply(cell_markers, length)
  clust_thresh <- (lcmval/2) + 0.5
  
  
  # Calculate a score based on whether the marker is found in the DEGs
  Weighted_Scores <- sapply(names(cell_markers), function (celltype){
    Markers_to_Test <- cell_markers[[celltype]]
    if(exact.match){
      score <- sum(Markers_to_Test %in% Positive_Genes)
    } else {
      score <- length(unlist(sapply(Markers_to_Test, grep, x=Positive_Genes)))
    }
    score * Score_Weights[[celltype]]
  })
  names(Weighted_Scores) <- names(cell_markers)
  
  # Assign cellType ID based on the weighted score of marker genes in DEG
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
  
  # Add ribosomal protein tag
  if(RP_Percent > 30 & rp){ ID <- paste("RP_",ID,sep = "") }
  
  # Add cycling tag
  if((Avg_Cluster_S_Score > 0.01 | Avg_Cluster_G2M_Score > 0.01) & cyc){
    CellCycleID <- "Cycling"
    ID <- paste("Cycling_",ID,sep = "")
  }
  
  return(ID)
}
