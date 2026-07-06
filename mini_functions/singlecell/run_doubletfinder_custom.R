# https://biostatsquid.com/doubletfinder-tutorial/
# run_doubletfinder_custom runs Doublet_Finder() and returns a dataframe with the cell IDs and a column with either 'Singlet' or 'Doublet'
run_doubletfinder_custom <- function(seu_sample_subset, multiplet_rate = NULL, sampleid='orig.ident'){
  # for debug
  #seu_sample_subset <- samp_split[[1]]
  # Print sample number
  print(paste0("Sample ", unique(seu_sample_subset[[sampleid]]), '...........')) 
  
  if(is.null(multiplet_rate)){
    print('multiplet_rate not provided....... estimating multiplet rate from cells in dataset')
    
    # 10X multiplet rates table
    #https://rpubs.com/kenneditodd/doublet_finder_example
    multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                      'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
                                      'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))
    
    print(multiplet_rates_10x)
    
    multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(seu_sample_subset@meta.data)) %>% 
      dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
      dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
    
    print(paste('Setting multiplet rate to', multiplet_rate))
  }
  
  # Pre-process seurat object with standard seurat workflow --- 
  sample <- NormalizeData(seu_sample_subset)
  sample <- FindVariableFeatures(sample)
  sample <- ScaleData(sample)
  sample <- RunPCA(sample, nfeatures.print = 10)
  
  # Find significant PCs
  stdv <- sample[["pca"]]@stdev
  percent_stdv <- (stdv/sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1] 
  co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - 
                       percent_stdv[2:length(percent_stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min_pc <- min(co1, co2)
  
  # Finish pre-processing with min_pc
  sample <- RunUMAP(sample, dims = 1:min_pc)
  sample <- FindNeighbors(object = sample, dims = 1:min_pc)              
  sample <- FindClusters(object = sample, resolution = 0.1)
  
  # pK identification (no ground-truth) 
  #introduces artificial doublets in varying props, merges with real data set and 
  # preprocesses the data + calculates the prop of artficial neighrest neighbours, 
  # provides a list of the proportion of artificial nearest neighbours for varying
  # combinations of the pN and pK
  sweep_list <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE)   
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats) # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  optimal.pk <- bcmvn %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
  
    ## Homotypic doublet proportion estimate
    annotations <- sample@meta.data$seurat_clusters # use the clusters as the user-defined cell types
    homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets
    
    nExp.poi <- round(multiplet_rate * nrow(sample@meta.data)) # multiply by number of cells to get the number of expected multiplets
    nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
    
    # run DoubletFinder
    sample <- doubletFinder(seu = sample, 
                            PCs = 1:min_pc, 
                            pK = optimal.pk, # the neighborhood size used to compute the number of artificial nearest neighbours
                            nExp = nExp.poi.adj) # number of expected real doublets
    # change name of metadata column with Singlet/Doublet information
    colnames(sample@meta.data)[grepl('DF.classifications.*', colnames(sample@meta.data))] <- "doublet_finder"
    
    # Subset and save
    # head(sample@meta.data['doublet_finder'])
    # singlets <- subset(sample, doublet_finder == "Singlet") # extract only singlets
    # singlets$ident
    double_finder_res <- sample@meta.data['doublet_finder'] # get the metadata column with singlet, doublet info
    double_finder_res <- rownames_to_column(double_finder_res, "row_names") # add the cell IDs as new column to be able to merge correctly
    return(double_finder_res)
}

# seusub <- append_doubletfinder(seusub, samp_split)
append_doubletfinder <- function(seu, res, col.name = 'doublet_finder'){
  sglt_dblt_metadata <- data.frame(bind_rows(res)) # merge to a single dataframe
  rownames(sglt_dblt_metadata) <- sglt_dblt_metadata$row_names # assign cell IDs to row names to ensure match
  sglt_dblt_metadata$row_names <- NULL
  
  seu <- AddMetaData(seu, sglt_dblt_metadata, col.name = col.name)
  tryCatch({
    doublets_summary <- seu@meta.data %>%
    group_by(batch, doublet_finder) %>%
    summarise(total_count = n(),.groups = 'drop') %>% as.data.frame() %>% ungroup() %>%
    group_by(batch) %>%
    mutate(countT = sum(total_count)) %>%
    group_by(doublet_finder, .add = TRUE) %>%
    mutate(percent = paste0(round(100 * total_count/countT, 2),'%')) %>%
    dplyr::select(-countT)
    print(doublets_summary)
  }, error=function(e){NULL})
  return(seu)
}