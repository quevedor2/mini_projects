mp_plot_diff_barplot_custom <- function(.data, .group, .size = 2, .taxas=NULL,
                                        errorbar.xmin = NULL, 
                                        errorbar.xmax = NULL,
                                        point.x = NULL, 
                                        taxa.class = 'all', group.abun = FALSE, 
                                        plotgrid=TRUE,
                                        removeUnknown=FALSE, ...){
  taxa.class <- rlang::enquo(taxa.class)
  .group <- rlang::enquo(.group)
  .size <- rlang::enquo(.size)
  errorbar.xmin <- rlang::enquo(errorbar.xmin)
  errorbar.xmax <- rlang::enquo(errorbar.xmax)
  point.x <- rlang::enquo(point.x)
  
  params <- list(...)
  taxa.class <- MicrobiotaProcess:::quo.vect_to_str.vect(taxa.class)
  if (!is.null(suppressMessages(taxatree(.data)))){
    tbl <- .data %>% mp_extract_abundance(taxa.class = 'all')
    xx <- .data %>% mp_extract_feature()
    xx %<>% dplyr::select(setdiff(colnames(xx), colnames(tbl)))
    tbl %<>% dplyr::left_join(xx, by = c('label'='OTU'))
  }else if (is.null(suppressMessages(taxatree(.data))) || taxa.class == 'OTU'){
    tbl <- .data %>% mp_extract_feature()
    if (!any(grepl("AbundanceBySample$", colnames(tbl)))){
      trash <- try(silent = TRUE,  
                   expr = {
                     .data <- suppressMessages(mp_cal_abundance(.data, .abundance = "Abundance"))  
                   }
      )
      if (inherits(trash, "try-error")) {
        .data <- mp_cal_abundance(.data, .abundance = "Abundance", force = TRUE)
      }
      tbl <- .data %>% mp_extract_feature()
    }
    tbl %<>% dplyr::rename(label = 'OTU') 
  }
  if(!is.null(.taxas)){
    .taxas_int <- intersect(.taxas, tbl$label)
    tbl <- data.frame("label"=.taxas_int) %>%
      dplyr::left_join(tbl, by='label') %>%
      as_tibble
    
  } 
  if (!any(taxa.class %in% c('all', 'ALL', 'All')) && !is.null(suppressMessages(taxatree(.data)))){
    tbl %<>% dplyr::filter(.data$nodeClass %in% taxa.class)
  }
  if (removeUnknown){
    tbl %<>% dplyr::filter(!grepl('__un_', .data$label))
  }
  nmda <- colnames(tbl)
  nm.abun <- nmda[grepl('BySample', nmda)][1]
  tbl %<>% tidyr::unnest(nm.abun)
  nmda <- colnames(tbl)
  if (rlang::quo_is_missing(.group)){
    if (any(grepl('^Sign_', nmda))){
      sign.group <- nmda[grepl('^Sign_', nmda)][1]
      .group <- gsub('Sign_', "", sign.group)
    }else{
      stop_wrap('The .group name should be specified manually.')
    }
  }else{
    .group <- rlang::as_name(.group)
    if (!grepl('^Sign_', .group)){
      if (paste0('Sign_', .group) %in% nmda){
        sign.group <- paste0('Sign_', .group)
      }else{
        #stop('Please check the mp_diff_analysis(..., action="add") has been run.')
        sign.group <- .group
      }
    }else{
      sign.group <- .group
      .group <- gsub('Sign_', "", .group)
    }
  }
  tbl %<>% dplyr::filter(!is.na(!!rlang::sym(sign.group)))
  plot.p2 <- TRUE
  # tbl %>% dplyr::filter(!is.na(Sign_Classification))
  if (any(grepl('LDA', nmda)) && rlang::quo_is_null(errorbar.xmin) && 
      rlang::quo_is_null(errorbar.xmax) && rlang::quo_is_null(point.x)){
    xlabtext <- bquote(paste(Log[10],"(",.("LDA"), ")"))
    xtext <- "LDAmean"
    xmintext <- "LDAlower"
    xmaxtext <- "LDAupper"
  }else if ('MDA' %in% nmda && rlang::quo_is_null(errorbar.xmin) && 
            rlang::quo_is_null(errorbar.xmax) && rlang::quo_is_null(point.x)){
    xlabtext <- "MDA"
    xtext <- "MDAmean"
    xmintext <- "MDAlower"
    xmaxtext <- "MDAupper"
  }else if (!rlang::quo_is_null(errorbar.xmin) && !rlang::quo_is_null(errorbar.xmax) && !rlang::quo_is_null(point.x)){
    xlabtext <- gsub("^\"|\"$", "", rlang::as_label(point.x))
    xtext <- xlabtext
    xmintext <- gsub("^\"|\"$", "", rlang::as_label(errorbar.xmin))
    xmaxtext <- gsub("^\"|\"$", "", rlang::as_label(errorbar.xmax))
    
  }else if (!rlang::quo_is_null(point.x) && any(rlang::quo_is_null(errorbar.xmin) || rlang::quo_is_null(errorbar.xmax))){
    xlabtext <- gsub("^\"|\"$", "", rlang::as_label(point.x))
    xtext <- xlabtext
    xmintext <- NULL
    xmaxtext <- NULL
  } else {
    warning("LDA or MDA has not been run; will be foregoing p2")
    plot.p2 <- FALSE
    xlabtext <- gsub("^\"|\"$", "", rlang::as_label(point.x))
    xtext <- xlabtext
    xmintext <- NULL
    xmaxtext <- NULL
  }
  if (any(grepl('Rel.*BySample', nmda))){
    abunda <- nmda[grepl('Rel.*BySample', nmda)]
  }else{
    abunda <- gsub('BySample$', '', nm.abun) 
  }
  if (group.abun){
    mapping1 <- aes(x = !!rlang::sym(abunda), 
                    y = !!rlang::sym("label"),
                    group = !!rlang::sym(.group),
                    fill = !!rlang::sym(.group)
    )
    panel1.geom <- 'geom_bar'
    panel1.args <- MicrobiotaProcess:::.extract_args(geom = 'geom_bar')
    panel1.args <- params[names(params) %in% panel1.args]
    panel1.args <- panel1.args[!names(panel1.args) %in% c('fill', 'group')]
    panel1.args$fun <- mean
    panel1.args$stat <- "summary"
    panel1.args$position <- "dodge"
    panel1.args$orientation <- 'y'
  }else{
    mapping1 <- aes(
      x = !!rlang::sym(abunda),
      y = !!rlang::sym("label"),
      fill = !!rlang::sym(.group)
    )
    panel1.geom <- 'geom_boxplot'
    panel1.args <- MicrobiotaProcess:::.extract_args(geom="geom_boxplot")
    panel1.args <- params[names(params) %in% panel1.args]
    panel1.args <- panel1.args <- panel1.args[!names(panel1.args) %in% c('fill')]
  }
  
  if (!is.null(xmintext)){
    mapping2 <- aes(
      xmin = !!rlang::sym(xmintext),
      xmax = !!rlang::sym(xmaxtext),
      y = !!rlang::sym('label')
    )
    panel2.geom <- 'geom_errorbarh'
  }else{
    mapping2 <- aes(x = 0, 
                    xend = !!rlang::sym(xtext), 
                    y = !!rlang::sym('label'),
                    yend = !!rlang::sym('label')
    )
    panel2.geom <- 'geom_segment'
  }
  
  panel2.args <- MicrobiotaProcess:::.extract_args(geom = panel2.geom)
  panel2.args <- params[names(params) %in% panel2.args]
  if (panel2.geom == 'geom_errorbarh'){
    panel2.args$height <- .3
    if ('height' %in% names(params)){
      panel2.args$height <- params$height        
    }
  }
  
  mapping3 <- aes_string(
    x = xtext, 
    y = 'label', 
    color = sign.group#, 
    #size = paste0("-log10(", rlang::as_name(.size), ")")
  )
  
  point.args <- MicrobiotaProcess:::.extract_args(geom='geom_point')
  point.args <- params[names(params) %in% point.args]
  
  point.args <- point.args[!names(point.args) %in% c("color", 'colour', "shape")]
  
  panel1.args$mapping <- mapping1
  panel2.args$mapping <- mapping2
  point.args$mapping <- mapping3
  point.args$show.legend <- c(colour=FALSE)
  
  if (MicrobiotaProcess:::is.quo.numeric(.size)){
    point.args <- point.args[!names(point.args) %in% c("size")]
    point.args$size <- as.numeric(rlang::as_label(.size))
  }else{
    if (inherits(rlang::quo_get_expr(.size), 'call')){
      size.mapping <- aes(size=!!.size)
    }else{
      if(any(grepl("\\(*.\\)$", rlang::quo_get_expr(.size)))){
        warning_wrap('The .y argument might be a call, you should remove the \' or \\" quote symbols.')
        size.mapping <- aes_string(size=gsub("^\"|\"$", "", rlang::as_label(.size)))
      }else{
        size.mapping <- aes_string(size=paste0("-log10(", gsub("^\"|\"$", "", rlang::as_label(.size)), ")"))
      }
    }
    point.args$mapping <- modifyList(point.args$mapping, size.mapping)
  }
  
  p <- ggplot(tbl)
  
  panel1.geom <- do.call(panel1.geom, panel1.args)
  p1 <- p + 
    panel1.geom + 
    ylab(NULL) + 
    xlab('Abundance') +
    scale_x_continuous(expand = c(0, 0, 0, .2))
  p1 <- p1 + ggfun::theme_blinds(colour = c("grey90", "white"), axis.ticks.y=element_blank(), 
                          panel.grid.major.x = element_blank(), 
                          panel.grid.minor.x = element_blank()
  )
  
  panel2.geom <- do.call(panel2.geom, panel2.args)
  point.geom <- do.call(geom_point, point.args)
  p2 <- p + 
    panel2.geom +
    point.geom +
    ylab(NULL) +
    xlab(xlabtext) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  p2 <- p2 + ggfun::theme_blinds(colour = c("grey90", "white"))
  fix.group <- levels(tbl[[.group]])
  if (is.null(fix.group)){
    fix.group <- unique(sort(tbl[[.group]]))
  }
  #fix.group.color <- p1$scales$get_scales('fill')$palette.cache
  fix.group.color <- scales::hue_pal()(length(fix.group))
  names(fix.group.color) <- fix.group
  fix.group.color <- fix.group.color[match(unique(as.character(tbl[[sign.group]])), names(fix.group.color))]
  p2 <- p2 + scale_color_manual(values=fix.group.color)
  
  # p <- p1 %>% aplot::insert_right(p2, width = .8)
  if(plotgrid){
    if(plot.p2){
      p <- cowplot::plot_grid(p1, p2, nrow=1)
    } else {
      p <- p1
    }
  } else {
    p <- list(p1, p2)
  }
  
  return(p)
}