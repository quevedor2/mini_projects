# Maps the labels of different VDJ packages to each other
# For example, given a seurat object of scRNA+scVDJ that was integrated 
# using immunarch, you can output a seurat object labeled using the 
# scRepertoire schema
mapTCRlabels <- function(seu, seu_scr, seu_djv, seu_imma, 
                         from.format='immunarch', to.format='immunarch'){
  # djvdj, scRepertoire, immunarch
  ord <- c('djvdj', 'scRepertoire', 'immunarch')
  map_keys <- list(
    "cdr3_nt" = setNames(c("cdr3_nt", "CTnt", "CDR3.nt"), ord),
    "cdr3_aa" = setNames(c("cdr3", "CTaa", 'CDR3.aa'), ord),
    "v"=setNames(c("v_gene", "CTgene", "V.name"), ord),
    "d"=setNames(c("d_gene", "CTgene", "D.name"), ord),
    "j"=setNames(c("j_gene", "CTgene", "J.name"), ord),
    "c"=setNames(c("c_gene", "CTgene", NA), ord),
    "chain"=setNames(c('chains', NA, 'chain'), ord),
    "clonotype_id"=setNames(c('clonotype_id', NA, 'raw_clonotype_id'), ord),
    "barcode"=setNames(rep("barcode", 3), ord)
  ) %>% do.call(rbind, .)
  
  tcr_vals <- lapply(setNames(rownames(map_keys),rownames(map_keys)), function(id){
    print(id)
    .getVal <- function(seu_i, map){
      tryCatch({
        if(map=='barcode'){
          colnames(seu_i)
        } else {
          seu_i@meta.data[,map]
        }
      }, error=function(e){NULL})
    }
    list("imma"=.getVal(seu_imma, map_keys[id,'immunarch']),
         "scr"=.getVal(seu_scr, map_keys[id,'scRepertoire']),
         "djv"=.getVal(seu_djv, map_keys[id,'djvdj']))
  })
  
  sel_common <- switch(from.format,
                       "scRepertoire"=.screpertoire2common(tcr_vals),
                       "immunarch"=.immunarch2common(tcr_vals),
                       "djvdj"=.djvdj2common(tcr_vals))
  vdj_df <- switch(to.format,
                   "scRepertoire"=.common2screpertoire(sel_common),
                   "immunarch"=.common2immunarch(sel_common),
                   "djvdj"=.common2djvdj(sel_common))
  
  seu_vdj <- AddMetaData(object = seu, metadata = vdj_df)
  return(seu_vdj)
}

##-- Aggregate TCR data ----
.screpertoire2common <- function(tcr_vals){
  ## Helper functions
  .getVDJC <- function(vdjc, gene){
    sapply(vdjc, function(vdjc_i){
      sapply(vdjc_i, function(vdjc_ij) vdjc_ij[gene]) %>% 
        paste(., collapse=";")
    })
  }
  .makeChain <- function(v){
    sapply(v, function(v_i){
      strsplit(v_i, split = ";") %>% 
        sapply(., function(i) gsub("^(TR[AB]).*", "\\1", i)) %>%
        paste(., collapse=";")
    })
  }
  
  ## Load in data
  scr_tcr <- lapply(tcr_vals, function(i) i$scr)
  
  ## Simple cdr3_nt and cdr3_aa
  common_tcr <- list()
  common_tcr[['nt']] <- gsub("_", ";", scr_tcr$cdr3_nt)  # ACGT;ACGT
  common_tcr[['aa']] <- gsub("_", ";", scr_tcr$cdr3_aa)  # CVNM;DTQY
  
  ## Parse VDJC
  genes <- strsplit(scr_tcr$v, split="_|;")
  vdjc <- lapply(genes, function(gene_i){
    lapply(gene_i, function(gene_ij){
      vdjc_ij <- strsplit(gene_ij, split="\\.")[[1]]
      .grepvdj <- function(pattern, vdjc_ij){
        res <- grep(pattern, vdjc_ij, value=T)
        ifelse(length(res)==0, NA, res)
      }
      vdjc_v <- c('v'=.grepvdj("TR[AB]V", vdjc_ij),
                  'd'=.grepvdj("TR[AB]D", vdjc_ij),
                  'j'=.grepvdj("TR[AB]J", vdjc_ij),
                  'c'=.grepvdj("TR[AB]C", vdjc_ij))
    })
  }) 
  common_tcr[['v']] <- .getVDJC(vdjc, 'v')               # TRAV;TRBV
  common_tcr[['d']] <- .getVDJC(vdjc, 'd')               # TRAD;TRBD
  common_tcr[['j']] <- .getVDJC(vdjc, 'j')               # TRAJ;TRBJ
  common_tcr[['c']] <- .getVDJC(vdjc, 'c')               # TRAC;TRBC
  
  ## Get clonotype and chains
  common_tcr[['chain']] <- .makeChain(common_tcr[['v']]) # TRA;TRB
  common_tcr[['clonotype']] <- rep(NA, length(scr_tcr$cdr3_nt)) # clonotype50
  common_tcr[['barcode']] <- scr_tcr$barcode
  common_tcr
}
.immunarch2common <- function(tcr_vals){
  ## Helper functions
  .makeId <- function(id){
    not_na <- which(!is.na(id))
    id[not_na] <- paste0("clonotype", id[not_na])
    return(id)
  }
  
  ## Load in data
  imma_tcr <- lapply(tcr_vals, function(i) i$imma)
  
  ## Simple cdr3_nt and cdr3_aa
  common_tcr <- list()
  common_tcr[['nt']] <- imma_tcr$cdr3_nt                  # ACGT;ACGT
  common_tcr[['aa']] <- imma_tcr$cdr3_aa                  # CVNM;DTQY
  
  ## Parse VDJC
  common_tcr[['v']] <- imma_tcr$v                         # TRAV;TRBV
  common_tcr[['d']] <- imma_tcr$d                         # TRAD;TRBD
  common_tcr[['j']] <- imma_tcr$j                         # TRAJ;TRBJ
  common_tcr[['c']] <- imma_tcr$c                         # TRAC;TRBC
  
  ## Get clonotype and chains
  common_tcr[['chain']] <- imma_tcr$chain                 # TRA;TRB
  common_tcr[['clonotype']] <- .makeId(imma_tcr$clonotype_id) # clonotype50
  common_tcr[['barcode']] <- imma_tcr$barcode
  common_tcr
}
.djvdj2common <- function(tcr_vals){
  ## Helper functions
  
  ## Load in data
  djv_tcr <- lapply(tcr_vals, function(i) i$djv)
  
  ## Simple cdr3_nt and cdr3_aa
  common_tcr <- list()
  common_tcr[['nt']] <- djv_tcr$cdr3_nt                  # ACGT;ACGT
  common_tcr[['aa']] <- djv_tcr$cdr3_aa                  # CVNM;DTQY
  
  ## Parse VDJC
  common_tcr[['v']] <- djv_tcr$v                         # TRAV;TRBV
  common_tcr[['d']] <- djv_tcr$d                         # TRAD;TRBD
  common_tcr[['j']] <- djv_tcr$j                         # TRAJ;TRBJ
  common_tcr[['c']] <- djv_tcr$c                         # TRAC;TRBC
  
  ## Get clonotype and chains
  common_tcr[['chain']] <- djv_tcr$chain                 # TRA;TRB
  common_tcr[['clonotype']] <- djv_tcr$clonotype_id      # clonotype50
  common_tcr[['barcode']] <- djv_tcr$barcode
  
  common_tcr
}

# used to make scRepertoire Frequency column
.getFreq <- function(mvdjc){
  na_idx <- is.na(mvdjc)
  freq <- rep(NA, length(mvdjc))
  freq[!na_idx] <- sapply(mvdjc[!na_idx], function(i) sum(i == mvdjc, na.rm = T))
  return(freq)
}
.common2screpertoire <- function(x_common, 
                                 clrange=c(Single=1, Small=5, Medium=20, 
                                           Large=100, Hyperexpanded=500)){
  ## Helper Functions
  .makeVDJC <- function(v,d,j,c){
    if(is.null(c)){
      c <- gsub("TR[AB][A-Z][0-9|-]*", "NA", v)
    }
    vdjc_l <- setNames(list(v,d,j,c), c("v", "d", "j", "c"))
    vdjc_spl <- lapply(vdjc_l, function(i){
      strsplit(i, split=";")
    })
    
    tra_trb <- lapply(seq_along(v), function(idx){
      is_na <- all(c(vdjc_spl$v[[idx]], vdjc_spl$d[[idx]],
                     vdjc_spl$j[[idx]], vdjc_spl$c[[idx]]) == 'NA')
      if(is_na) return(NA)
      vdj_grp_tra <- paste(vdjc_spl$v[[idx]], vdjc_spl$j[[idx]], 
                           vdjc_spl$c[[idx]], sep=".")
      vdj_grp_trb <- paste(vdjc_spl$v[[idx]], vdjc_spl$d[[idx]], 
                           vdjc_spl$j[[idx]], vdjc_spl$c[[idx]],
                           sep=".")
      
      # TRA is made of VJC;  TRB is made of VDJC
      # Remove the 'D' from TRAs
      vdj_grp <- setNames(vdj_grp_trb, gsub(".*(TR[A|B]).*", "\\1", vdj_grp_trb))
      tra_idx <- which(names(vdj_grp) == 'TRA')
      vdj_grp[tra_idx] <- vdj_grp_tra[tra_idx]
      
      return(vdj_grp)
    })
    return(tra_trb)
  }
  # Used to make scRepertoire CTgene column
  .mergeVDJC <- function(vdjc){
    sapply(vdjc, function(vdjc_i){
      if(is.na(vdjc_i)) return(NA)
      
      # merge multiple TRA's and TRB's with ;
      vdjc_spl <- split(vdjc_i, names(vdjc_i))
      tra_trb <- lapply(vdjc_spl, function(i) paste(i, collapse=";"))
      
      # merge TRA with TRB using _
      paste(tra_trb$TRA, tra_trb$TRB, sep="_")
    })
  }
  # used to make scRepertoire CTscript column
  .mergeVDJCnt <- function(vdjc, nt){
    ntl <- strsplit(nt, split="_")
    
    sapply(seq_along(vdjc), function(idx){
      vdjc_i <- vdjc[[idx]]
      nt_i <- ntl[[idx]]
      
      if(is.na(vdjc_i)) return(NA)
      
      # merge multiple TRA's and TRB's with ;
      vdjc_spl <- split(vdjc_i, names(vdjc_i))
      nt_spl <- split(nt_i, names(vdjc_i))
      
      vdjc_tra_trb <- lapply(vdjc_spl, function(i) paste(i, collapse=";"))
      nt_tra_trb <- lapply(nt_spl, function(i) paste(i, collapse=";"))
      
      # merge TRA with TRB using _
      paste(vdjc_tra_trb$TRA, nt_tra_trb$TRA, tra_trb$TRB, nt_tra_trb$TRB, sep="_")
    })
  }
  # used to make scRepertoire cloneType column
  .getClonetype <- function(freq){
    na_idx <- is.na(freq)
    ctype <- rep(NA, length(freq))
    freq_cut <- cut(freq, c(0, clrange, Inf))
    cl_ids <- paste0(c(names(clrange), "Infinite"), 
                     " (", c(0, clrange), 
                     " < X <= ", c(clrange, Inf), ")")
    levels(freq_cut) <- cl_ids
    return(as.character(freq_cut))
  }
  
  ctnt <- gsub(";", "_", x_common$nt)
  ctaa <- gsub(";", "_", x_common$aa)
  
  vdjc <- .makeVDJC(x_common$v, x_common$d, x_common$j, x_common$c)
  ctgene <- .mergeVDJC(vdjc)
  ctstrict <- .mergeVDJCnt(vdjc, ctnt)
  freq <- .getFreq(ctgene)
  clonetype <- .getClonetype(freq)
  
  df <- data.frame("CTgene"=ctgene,
                   "CTnt"=ctnt,
                   "CTaa"=ctaa,
                   "CTstrict"=ctstrict,
                   "Frequency"=freq,
                   "cloneType"=clonetype)
  rownames(df) <- x_common$barcode
  return(df)
}
.common2immunarch <- function(x_common){
  ## Helper Functions
  
  ## Main
  freq <- .getFreq(x_common$aa)
  prop <- freq / sum(!is.na(freq))
  
  df <- data.frame("Clones"=freq,
                   "Proportion"=prop,
                   "CDR3.nt"=x_common$nt,
                   "CDR3.aa"=x_common$aa,
                   "V.name"=x_common$v,
                   "D.name"=x_common$d,
                   "J.name"=x_common$j,
                   "V.end"=NA,
                   "D.start"=NA,
                   "D.end"=NA,
                   "J.start"=NA,
                   "VJ.ins"=NA,
                   "VD.ins"=NA,
                   "DJ.ins"=NA,
                   "Sequence"=x_common$nt,
                   "chain"=x_common$chain,
                   "raw_clonotype_id"=gsub("clonotype", "", x_common$clonotype))
  rownames(df) <- x_common$barcode
  return(df)
}
.common2djvdj <- function(x_common){
  ## Helper Functions
  .getNChains <- function(chains){
    na_idx <- (chains=='NA')
    nchains <- rep(NA, length(chains))
    nchains[(!na_idx)] <- sapply(strsplit(chains, split=";")[(!na_idx)], length)
    return(nchains)
  }
  .getLen <- function(seq_i){
    na_idx <- is.na(seq_i) | (seq_i=='NA')
    slen <- rep(NA, length(seq_i))
    slen[!na_idx] <- sapply(strsplit(seq_i[!na_idx], split=";"), function(i){
      paste(nchar(i), collapse=";")
    })
    return(slen)
  }
  .getPaired <- function(chains){
    na_idx <- is.na(chains) | (chains=='NA')
    paired <- rep(NA, length(chains))
    paired[!na_idx] <- (grepl("TRA", chains) & grepl("TRB", chains))[!na_idx]
    return(paired)
  }
  .listTRUE <- function(chains){
    na_idx <- is.na(chains) | (chains=='NA')
    truev <- rep(NA, length(chains))
    truev[!na_idx] <- gsub("TR[AB]", "TRUE", chains[!na_idx])
    return(truev)
  }
  
  ## Main
  nchains <- .getNChains(x_common$chain)
  aa_len <- .getLen(x_common$aa)
  nt_len <- .getLen(x_common$nt)
  chains <- x_common$chain
  paired <- .getPaired(chains)
  truev <- .listTRUE(chains)
  
  df <- data.frame("clonotype_id"=x_common$clonotype,
                   "chains"=chains,
                   "n_chains"=nchains,
                   "cdr3"=x_common$aa,
                   "cdr3_nt"=x_common$nt,
                   "cdr3_length"=aa_len,
                   "cdr3_nt_length"=nt_len,
                   "v_gene"=x_common$v,
                   "d_gene"=x_common$d,
                   "j_gene"=x_common$j,
                   "c_gene"=ifelse(is.null(x_common$c), NA, x$common),
                   "reads"=NA,
                   "umis"=NA,
                   "productive"=truev,
                   "full_length"=truev,
                   "paired"=paired)
  rownames(df) <- x_common$barcode
  return(df)
}


