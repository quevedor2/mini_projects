#' orthogene
#' @description converts a vector of gene names from mouse to human
#' or vice versa
#' 
#' @param x character [n]: vector of gene names
#' @param from character [1]: either 'mouse' or 'human'
#' @param to character [1]: either 'mouse' or 'human'
#' @param verbose boolean: identify which genes could not be converted
#'
#' @return character [n]: a vector of converted gene names
#'
#' @examples orthogene(c('Mtor', 'Pten', 'Tp53'), from='mouse', to='human')
orthogene <- function(x, from='mouse', to='human', verbose=FALSE){
  stopifnot((from == 'mouse' | from == 'human'))
  stopifnot((to == 'mouse' | to == 'human'))
  to <- ifelse(to=='human', 'Gene.HS', 'Gene.MM')
  from <- ifelse(from=='human', 'Gene.HS', 'Gene.MM')
  
  require(ProjecTILs)
  data(Hs2Mm.convert.table)
  ortholog_table <- Hs2Mm.convert.table
  
  ortho_g <- ortholog_table[[to]][match(x, ortholog_table[[from]])]
  nas <- is.na(ortho_g)
  if(any(nas)){
    alt.from <- ifelse(from=='Gene.HS', 'Alt.symbol.HS', 'Alt.symbol')
    ortho_g2 <- ortholog_table[[to]][match(x[which(nas)], ortholog_table[[alt.from]])]
    ortho_g[which(nas)] <- ortho_g2
  }
  
  nas <- is.na(ortho_g)
  if(any(nas) & verbose){
    cat(paste0("The following genes were not found and were dropped: \n",
               paste(paste0("\t- ", x[nas]), collapse="\n"), "\n" ))
  }
  ortho_g <- ortho_g[!nas]
  return(ortho_g)
}

#' orthogene.scgate
#' @description
#' @param obj 
#' @param from 
#' @param to 
#'
#' @return
#' @export
#'
#' @examples
orthogene.scgate <- function(obj, from='mouse', to='human'){
  df <- obj@misc$scGate[[from]]
  ortho_gs <- sapply(df$signature, function(genelist){
    x <- strsplit(genelist, split=";")[[1]]
    ortho_g <- orthogene(x, from=from, to=to)
    paste(ortho_g, collapse=";")
  })
  df$signature <- as.character(ortho_gs)
  obj@misc$scGate[[to]] <- df
  return(obj)
}