

list2df <- function(lst, id.col=NULL, value.col=NULL){
  if(any(sapply(lst, class) == 'NULL')){
    lst <- lst[-which(sapply(lst, class) == 'NULL')]
  }
  if(all(sapply(lst, class) == 'list')){
    lst <- lapply(lst, list2df)
  }
  if(all(sapply(lst, class) == 'data.frame')){
    newcol <- make.unique(c(colnames(lst[[1]]), 'Var'))[ncol(lst[[1]]) + 1]
    do.call(rbind, lst) %>% 
      mutate(!!rlang::sym(newcol) := rep(names(lst), sapply(lst, nrow))) %>% 
      tibble::remove_rownames() 
  } else if(all(sapply(lst, class) == 'character')) {
    data.frame("Var1"=rep(names(lst), sapply(lst, length)),
               "Var2"=as.character(unlist(lst)))
  }
}