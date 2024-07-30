## Added the option to use a "taxa_code" such as "f__"
# This will extract all families prefixed with f__ and not include
# anything downstream that also has f__ in its classification (such as
# f__familyX|g__genusY).  if the taxa_Code is left as NULL, it will
# operate mp_import_metaphlan in the exaact same manner as defined
# in the MicrobiotaProcess package
mp_import_metaphlan_custom <- function (profile, mapfilename = NULL,
                                        treefile = NULL, linenum = NULL, 
                                        taxa_code=NULL,...) {
  skipnrow <- MicrobiotaProcess:::guess_skip_nrow(profile)
  if (!is.null(linenum)) {
    sampleda <- utils::read.table(profile, sep = "\t", skip = skipnrow, 
                                  nrow = linenum, comment.char = "", quote = "") %>% 
      dplyr::mutate_at("V1", as.character) %>% dplyr::mutate(V1 = tidyr::replace_na(.data$V1, 
                                                                                    paste0("unknown", seq_len(sum(is.na(.data$V1)))))) %>% 
      tibble::column_to_rownames(var = "V1") %>% t() %>% 
      as.data.frame()
    da <- utils::read.table(profile, sep = "\t", skip = linenum + 
                              skipnrow, comment.char = "", quote = "")
    sampleindx <- da %>% vapply(., function(x) is.numeric(x), 
                                logical(1)) %>% which(.)
    colnames(da)[sampleindx] <- paste0("sample", sampleindx - 
                                         1)
    sampleda %<>% dplyr::slice(sampleindx - 1) %>% magrittr::set_rownames(paste0("sample", 
                                                                                 sampleindx - 1))
  } else {
    da <- utils::read.table(profile, sep = "\t", head = TRUE, 
                            skip = skipnrow, comment.char = "", quote = "", check.names = FALSE)
    sampleda <- NULL
  }
  
  
  if (is.null(linenum) && !is.null(mapfilename)) {
    if (inherits(mapfilename, "character") && file.exists(mapfilename)) {
      sampleda <- MicrobiotaProcess:::read_qiime_sample(mapfilename)
    } else if (!inherits(mapfilename, "data.frame")) {
      rlang::abort("The mapfilename should be a file or data.frame contained sample information!")
    } else {
      sampleda <- mapfilename
    }
  }
  
  if (!is.null(treefile)) {
    otutree <- ape::read.tree(treefile)
    otutree$tip.label %<>% strsplit("\\|") %>% lapply(., 
                                                      function(x) x[length(x)]) %>% unlist()
    otutree %<>% treeio::as.treedata()
  } else {
    otutree <- NULL
  }
  
  clnm <- colnames(da)
  # if(!is.null(taxa_code)){
  #   da %<>% dplyr::filter(!!as.symbol(clnm[1]) %>%
  #                           grepl(paste0("^.*", taxa_code, "[^|]*$"), .))
  # }
  
  if(is.null(taxa_code)){
    max.sep <- da %>% dplyr::pull(1) %>% 
      gregexpr("\\|", .) %>% 
      lapply(length) %>% unlist %>% max
    
    dat <- da %>% 
      dplyr::filter(!!as.symbol(clnm[1]) %>% 
                      gregexpr("\\|", .) %>% 
                      lapply(length) %>% 
                      unlist %>% 
                      magrittr::equals(max.sep))
    dat %<>% magrittr::set_rownames(paste0("row", seq_len(nrow(dat))))
    
    taxatab <- dat %>% 
      dplyr::select(1) %>% 
      split_str_to_list(sep = "\\|") %>% 
      magrittr::set_colnames(c(MicrobiotaProcess:::taxaclass[seq_len(max.sep)], "OTU")) %>% 
      MicrobiotaProcess:::fillNAtax()
    
    assay <- dat %>% 
      dplyr::mutate(`:=`(!!as.symbol(clnm[1]), 
                         strsplit(!!as.symbol(clnm[1]), split = "\\|") %>% 
                           lapply(function(x) x %>% magrittr::extract2(max.sep + 1)))) %>% 
      dplyr::select(-1)
  } else {
    dat <- da %>% dplyr::filter(!!as.symbol(clnm[1]) %>%
                                  grepl(paste0("^.*", taxa_code, "[^|]*$"), .))
    dat %<>% magrittr::set_rownames(paste0("row", seq_len(nrow(dat))))

    max.sep <- dat %>% dplyr::pull(1) %>%
      gregexpr("\\|", .) %>%
      lapply(length) %>% unlist %>% max

    taxatab <- dat %>%
      dplyr::select(1) %>%
      split_str_to_list(sep = "\\|") %>%
      magrittr::set_colnames(c(MicrobiotaProcess:::taxaclass[seq_len(max.sep)], "OTU")) %>%
      MicrobiotaProcess:::fillNAtax()
    assay <- dat %>%
      dplyr::mutate(`:=`(!!as.symbol(clnm[1]),
                         strsplit(!!as.symbol(clnm[1]), split = "\\|") %>%
                           lapply(function(x) x[length(x)]))) %>%
      dplyr::select(-1)
  }
  
  featureIndex <- intersect(rownames(assay), rownames(taxatab))
  taxatab %<>% 
    dplyr::mutate(OTU=make.unique(OTU)) %>% 
    magrittr::extract(featureIndex, ) %>% 
    magrittr::set_rownames(NULL) %>% 
    tibble::column_to_rownames(var = "OTU")
  assay %<>% magrittr::extract(featureIndex, ) %>% 
    magrittr::set_rownames(rownames(taxatab))
  rowdaindx <- assay %>% vapply(., function(x) !is.numeric(x), 
                                logical(1)) %>% which() %>% unname()
  rowdaindx <- union(rowdaindx, which(grepl("NCBI_tax_id", 
                                            colnames(assay), ignore.case = TRUE)))
  if (length(rowdaindx) > 0) {
    rowda <- assay %>% dplyr::select(rowdaindx)
    assay %<>% dplyr::select(-rowdaindx)
  } else {
    rowda <- NULL
  }
  res <- list(otutab = assay, otutree = otutree, sampleda = sampleda, 
              taxatab = taxatab, otu.metada = rowda)
  mpse <- MicrobiotaProcess:::.build_mpse(res)
  return(mpse)
}
