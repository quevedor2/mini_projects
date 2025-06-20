
#' immdataToBarcode
#' @author Rene Q
#' @description Takes in the input from immunarch::repLoad()
#' and aims to simplify the data for integration with Seurat object by
#' removing long data (e.g. barcode and contig), and reduce 
#' to a single non-repetitive barcode and information
#' 
#' @param immdata_10x_data Output from immunarch::repLoad()$data$filtered_contig_annotations (v_0.9.0)
#' @param rm_dup Whether to remove duplicate barcode entries and separate them
#' @param colprefix Append a prefix to the column titles to identify the 
#' VDJ data in your seurat object
#'
#' @return dataframe containing a simplified immunarch::repLoad()$data$filtered_contig_annotations
#' @export
immdataToBarcode <- function(immdata_10x_data, rm_dup=TRUE, colprefix=''){
  # Simplify data: remove long data (e.g. barcode and contig), and reduce 
  # to a single non-repetitive barcode and information
  immdata_data <- apply(immdata_10x_data, 1, function(i){
    as.data.frame(cbind(
      data.frame(t(i[-which(names(i) %in% c('Barcode', 'ContigID'))])),  #remove barcode and contigid
      data.frame("Barcode"=strsplit(i['Barcode'], split=";")[[1]])
    ))
  })
  immdata_data <- as.data.frame(do.call(rbind, immdata_data))
  dup_immdata_data <- NULL
  if(rm_dup & any(duplicated(immdata_data$Barcode))){
    ## Some duplicate entries exist because one barcode has a raw_clonotype_id, while
    ## the duplicate entry doesnt.  This function will select the barcode entry id that
    ## has a corresponding raw_clonotype_id. The "removed" duplicate entries are 
    ## reported in a separate dataframe
    dupentries <- unique(immdata_data$Barcode[which(duplicated(immdata_data$Barcode))])
    immdata_data[(immdata_data$Barcode %in% dupentries),] 
    ilsplit <- immdata_data %>% 
      dplyr::filter(Barcode %in% dupentries) %>%
      dplyr::arrange((raw_clonotype_id)) %>%
      split(., f=.$Barcode)
    
    immdata_data <- immdata_data %>% 
      dplyr::filter(!Barcode %in% dupentries) %>%
      rbind(., do.call(rbind, lapply(ilsplit, function(i) i[1,,drop=F])))
    dup_immdata_data <- do.call(rbind, lapply(ilsplit, function(i) i[-1,,drop=F]))
  }
  if(colprefix != '') colnames(immdata_data) <- paste0(colprefix, ".", colnames(immdata_data))
  # rownames(immdata_data) <- gsub("-[1-9]$", "", immdata_data$Barcode)
  return(list("immdat"=immdata_data, "duplicates"=dup_immdata_data))
}

#' appendImmdataMeta
#' @author Rene Q
#' @description Appends the data from immdataToBarcode() function to the 
#' Seurat object passed in. 
#' Currently, this function is designed to run on the output of
#' cellranger multi + cellranger aggr, so an aggregate file mapping the 
#' index to the sample names is required (i.e. -1 => SampleX)
#' 
#' @param seu Seurat object with a column in the metadata that has sample IDs
#' @param immdat output from immdataToBarcode()
#' @param aggfile Aggregate file that can map samples to index numbers
#' @param samplecol Sample column in seu@meta.data
#'
#' @return Seurat object with immune data appended in metadata
#' @export
appendImmdataMeta <- function(seu, immdat, aggfile, samplecol='batch'){
  aggdata <- read.csv(aggfile) %>%
    pull(sample_id) %>%
    setNames(., paste0("-", seq_along(.)))
  
  metatmp <- seu@meta.data %>% 
    tibble::rownames_to_column(., "barcode") %>% 
    dplyr::mutate(barcode = gsub("-.*", "", barcode))
  
  samples <- immdat %>% 
    dplyr::pull(grep("Barcode", colnames(.), value=T)) %>% 
    gsub("^.*-", "-", .)
  immdat <- immdat %>% 
    dplyr::mutate(!!sym(samplecol) := aggdata[samples],
                  barcode=gsub("-.*", "", !!sym(grep("Barcode", colnames(.), value=T))))
  
  exclusive_to_immdat <- dplyr::anti_join(immdat, metatmp, by = c('barcode', samplecol))
  if(nrow(exclusive_to_immdat)>0) warning(paste0(nrow(exclusive_to_immdat), " barcodes had no match in seurat object"))
  metatmp <- dplyr::left_join(metatmp, immdat, by=c('barcode', samplecol)) %>%
    magrittr::set_rownames(Cells(seu)) %>%
    dplyr::select(-barcode)
  seu@meta.data <- metatmp
  
  
  return(list("seu"=seu, "immdat_only"=exclusive_to_immdat))
}


#' .convertMetadataToImmunarch
#' @author Rene Q
#' @description Takes the VDJ appended to the Seurat object and cleans it up
#' so it'll work with the standard immunarch functions
#' 
#' @param metadat seu@meta.data after running appendImmdataMeta() function
#' @param sample Sample identity column.  If present, will split metadata
#' into sample-specific dataframes; else, it'll treat the given metadata as a 
#' single sample
#' @param pattern VDJ pattern prefix in the metadata
#' @param rm.empty.entries Removes sample-specific entries that have no VDJ data
#' @param rm.nontcr.dat Removes samples where there is no VDJ info, needed 
#' to work with some immunarch functions [Default=TRUE]
#'
#' @return A metadata list containing cleaned up data for immunarch
#' @export
.convertMetadataToImmunarch <- function(metadat, sample=NULL, pattern="VDJT", 
                                        rm.empty.entries=TRUE, rm.nontcr.dat=TRUE){
  message(paste0("Operating on: ", pattern))
  colnames(metadat) <- gsub(paste0("^", pattern ,"."), "", colnames(metadat))
  if(rm.nontcr.dat) {
    metadat <- metadat %>% 
      dplyr::filter(!is.na(!!rlang::sym(immunarch:::IMMCOL$count)))
  }
  
  # Splits metadata into sample-specific batches
  if(!is.null(sample)) {
    metadat <- split(metadat, f=metadat[,sample])
  }
  
  # Removes sample-specific entries that have no VDJ data
  if(rm.empty.entries){
    cl <- class(metadat)
    if(cl != 'list') metadat <- list("Sample"=metadat)
    vcheck <- sapply(seq_along(metadat), function(i){
      metadat[[i]] %>% dplyr::select(immunarch:::IMMCOL$count) %>% 
        table %>% dim
    })
    if(any(vcheck == 0)) metadat <- metadat[-which(vcheck==0)]
    if(cl != 'list') metadat <- metadat[[1]]
  }
  return(metadat)
}

#' .makeImmunarchMeta
#' @author Rene Q
#' @description Small function needed to create the metadata for .immunarchViz()
#' and immunarch::vis group splitting
#' 
#' @param seu Seurat object
#' @param column metadata column containing the condition to split by
#' @param samplecol Sample column
#'
#' @return
#' @export
.makeImmunarchMeta <- function(seu, column='condition', samplecol='batch'){
  .meta <- seu@meta.data[,c(samplecol, column)] %>% 
    # dplyr::select(!!rlang::sym(samplecol), !!rlang::sym(column)) %>% 
    unique %>%
    dplyr::rename_with(., ~"Sample", .cols=!!rlang::sym(samplecol))
  return(.meta)
}

#' .runImmunarchAnalysis
#' @author Rene Q
#' @description Runs the gauntlet of immunarch functions from the immunarch
#' tutorial. Count, diversity, etc..
#' 
#' @param seu Seurat object with VDJ data appended via appendImmdataMeta()
#' @param pattern VDJ prefix pattern identifying the specific immunedata to work on
#' @param sample sample column id in the meta.data
#'
#' @return A list containing all the outputs from immunarch
#' @export
.runImmunarchAnalysis <- function(seu, pattern='VDJT', sample='batch'){
  .data <- seu@meta.data %>% 
    .convertMetadataToImmunarch(., pattern=pattern, sample=sample)
  
  res <- list()
  res[['exp_vol']] <- .data %>%
    immunarch::repExplore(., .method = "volume")
  res[['exp_len']] <- .data %>%
    immunarch::repExplore(.,  .method = "len", .col = "aa")
  res[['exp_cnt']] <-  .data %>%
    immunarch::repExplore(., .method = "count")
  
  res[['imm_top']] <- .data %>%
    immunarch::repClonality(., .method = "top", .head = c(10, 100, 1000, 3000, 10000))
  res[['imm_rare']] <- .data %>%
    immunarch::repClonality(., .method = "rare")
  res[['imm_hom']] <- .data %>%
    immunarch::repClonality(.,
                            .method = "homeo",
                            .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1))
  return(res)
}

#' .immunarchViz
#' @author Rene Q
#' @description Visualizes all the results from .runImmunarchAnalysis
#'
#' @param res Output from .runImmunarchAnalysis()
#' @param .meta Outut from .makeImmunarchMeta 
#' @param .by group to split by
#' @param outf out filename
#' @param width width of pdf
#'
#' @return
#' @export
.immunarchViz <- function(res, .meta=NULL, .by=NULL, outf=NULL, width=20){
  p1 <- immunarch::vis(res$exp_len)
  p2 <- immunarch::vis(res$exp_cnt)
  p3 <- immunarch::vis(res$exp_vol)
  p4 <- immunarch::vis(res$exp_len, .by = .by, .meta = .meta)
  p5 <- immunarch::vis(res$exp_cnt, .by = .by, .meta = .meta)
  p6 <- immunarch::vis(res$exp_vol, .by = .by, .meta = .meta)
  
  p7 <- immunarch::vis(res$imm_top)
  p7b <- immunarch::vis(res$imm_top, .by = .by, .meta = .meta)
  p8 <- immunarch::vis(res$imm_rare) 
  p8b <- immunarch::vis(res$imm_rare, .by = .by, .meta = .meta)
  p9 <- immunarch::vis(res$imm_hom) 
  p9b <- immunarch::vis(res$imm_hom, .by = .by, .meta = .meta)
  
  pdf(outf, width=width)
  if(is.null(.meta)){
    plot(p1)
    plot(p2 + p3)
    plot(p7)
    plot(p8)
    plot(p9)
  } else {
    plot(p1)
    plot(p2 + p3)
    plot(p4)
    plot(p5 + p6)
    plot(p7 + p7b)
    plot(p8 + p8b)
    plot(p9 + p9b)
  }
  dev.off()
}
