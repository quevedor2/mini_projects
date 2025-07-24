gptcelltype <- function(input, tissuename=NULL, model='gpt-4', topgenenumber = 10) {
  if (class(input)=='list') {
    input <- sapply(input,paste,collapse=',')
  } else {
    input <- input[input$avg_log2FC > 0,,drop=FALSE]
    input <- tapply(input$gene,list(input$cluster),function(i) paste0(i[1:topgenenumber],collapse=','))
  }
  
  message = paste0('Identify cell types of ',tissuename,' cells using the following markers separately for each\n row. Only provide the cell type name. Do not show numbers before the name.\n Some can be a mixture of multiple cell types. ',  "\n", paste0(names(input), ':',unlist(input),collapse = "\n"))
  
  return(message)
}