# layout = "radial"
# tree.type = "taxatree"
# .taxa.class = NULL
# barplot.x = NULL
# point.size = NULL
# sample.num = 50
# tiplab.size = 2
# offset.abun = 0.04
# pwidth.abun = 0.8
# offset.effsize = 0.3
# pwidth.effsize = 0.5
# group.abun = FALSE
# tiplab.linetype = 3
# 
# #@@
# .data = mpse3
# .group = !!rlang::sym(effect)
# .group = (effect)
# pwidth.abun=0.1
# 
# # mpse3.2 <- mpse3.bkup %>% 
# #   mp_diff_analysis(.abundance=RareAbundance, 
# #                    .group=!!rlang::sym(effect), 
# #                    strict=F,
# #                    cl.min=3,
# #                    ldascore=0) 
# # p <- mpse3.2 %>%
# #   mp_plot_diff_res(
# #     group.abun = TRUE,
# #     pwidth.abun=0.1
# #   ) 
# 

mp_plot_diff_res_custom <- function(.data,
                                    .group,
                                    layout = "radial",
                                    tree.type = "taxatree",
                                    .taxa.class = NULL,
                                    barplot.x = NULL,
                                    point.size = NULL,
                                    .taxons=NULL,
                                    sample.num = 50,
                                    tiplab.size = 2,
                                    offset.abun = 0.04,
                                    pwidth.abun = 0.8,
                                    offset.effsize = 0.3,
                                    pwidth.effsize = 0.5,
                                    group.abun = FALSE,
                                    tiplab.linetype = 3,
                                    ...
){
  .taxa.class <- rlang::enquo(.taxa.class)
  .group <- rlang::enquo(.group)
  barplot.x <- rlang::enquo(barplot.x)
  point.size <- rlang::enquo(point.size)
  
  layout %<>% match.arg(c("rectangular", "roundrect", "ellipse", "circular", 
                          "slanted", "radial", "inward_circular"))
  
  tree.type %<>% match.arg(c("taxatree", "otutree"))
  
  if (tree.type == 'otutree'){
    anno.tree <- .data %>% mp_extract_tree(type="otutree") %>% suppressMessages()
    if (is.null(anno.tree)){
      stop_wrap("The otutree slot of the MPSE class is empty, you can try to 
                       select taxatree by setting tree.type to 'taxatree'.")
    }else{
      taxada <- .data %>% mp_extract_taxonomy() %>% suppressMessages()
      if (!is.null(taxada)){
        anno.tree %<>% dplyr::left_join(taxada, by=c("label"="OTU"))
        if (rlang::quo_is_null(.taxa.class)){
          .taxa.class <- rlang::sym(colnames(taxada)[3])
        }
      }
    }
  }
  
  if (tree.type == "taxatree"){
    anno.tree <- .data %>% mp_extract_tree() %>% suppressMessages()
    if (is.null(anno.tree)){
      stop_wrap("The taxatree slot of the MPSE class is empty, you can try to 
                       select otutree by setting tree.type to 'otutree'.")
    }else{
      if (rlang::quo_is_null(.taxa.class)){
        .taxa.class <- anno.tree %>% 
          tidytree::filter(!!rlang::sym("nodeDepth")==2, keep.td=FALSE) %>% 
          pull(!!rlang::sym("nodeClass")) %>% 
          unique()
        .taxa.class <- rlang::sym(.taxa.class)
      }
    }
  }
  
  tbl.f <- .data %>% mp_extract_feature()
  tbl.f %<>% dplyr::select(c(colnames(tbl.f)[1], 
                             setdiff(colnames(tbl.f)[-1], tidytree::get.fields(anno.tree)))
  )
  
  anno.tree %<>% dplyr::left_join(tbl.f, by=c('label' = 'OTU'))
  
  nsample <- .data %>% mp_extract_sample() %>% nrow()
  field.da.nm <- tidytree::get.fields(anno.tree)
  
  if (any(grepl("LDAmean", field.da.nm)) && rlang::quo_is_null(barplot.x)){
    x.bar <- 'LDAmean'
    x.bar.title <- "log10(LDA)"
  }else if (any(grepl("MDAmean", field.da.nm)) && rlang::quo_is_null(barplot.x)){
    x.bar <- "MDAmean"
    x.bar.title <- "MDA"
  }else if (!rlang::quo_is_null(barplot.x)){
    x.bar <- x.bar.title <- gsub("^\"|\"$", "", rlang::as_label(barplot.x))
  }else{
    stop_wrap("Please provide barplot.x or verify the 'mp_diff_analysis' is done.")
  }
  
  if (rlang::quo_is_missing(.group)){
    if (any(grepl('^Sign_', field.da.nm))){
      sign.field <- field.da.nm[grepl('^Sign_', field.da.nm)][1]
      group.nm <- gsub('Sign_', "", sign.field)
    }else{
      stop_wrap('The .group name should be specified manually.')
    }
  }else{
    group.nm <- rlang::as_name(.group)
    if (!grepl('^Sign_', group.nm)){
      if (paste0('Sign_', group.nm) %in% field.da.nm){
        sign.field <- paste0('Sign_', group.nm)
      }else{
        #stop('Please check the mp_diff_analysis(..., action="add") has been run.')
        sign.field <- group.nm
      }
    }else{
      sign.field <- group.nm
      group.nm <- gsub('Sign_', "", group.nm)
    }
  }
  
  if (rlang::quo_is_null(point.size) && 'fdr' %in% field.da.nm){
    size.mapping <- '-log10(fdr)'
  }else if (!rlang::quo_is_null(point.size)){
    if (inherits(rlang::quo_get_expr(point.size), 'call')){
      size.mapping <- rlang::as_label(point.size)
    }else{
      size.mapping <- paste0("-log10(", gsub("^\"|\"$", "", rlang::as_label(point.size)), ")")
    }
  }else{
    stop_wrap("Please specify the 'point.size' argument manually !")
  }
  
  flag <- grepl(paste0("By", group.nm), field.da.nm)
  
  if (nsample > sample.num || group.abun){
    if (!any(flag)){
      stop_wrap("The relative abundance of each group will be displayed, but the
                       relative abundance of each group is not calculated, please run 
                       the mp_cal_abundance specified group argument before !")
    }
    abun.col <- field.da.nm[flag]
  }else{
    abun.col <- field.da.nm[grepl("BySample", field.da.nm)]
  }
  abun.col <- abun.col[1]
  x.abun.col <- anno.tree %>% 
    dplyr::select(!!rlang::sym(abun.col)) %>%
    tidyr::unnest(!!rlang::sym(abun.col)) %>%
    colnames()
  if (any(grepl('BySample$',  abun.col))){
    if (any(grepl('^Rel', x.abun.col))){
      x.abun.col <- paste0('Rel', abun.col)
    }else{
      x.abun.col <- gsub('BySample$', '', abun.col)
    }
  }else{
    if (any(grepl("^Rel", x.abun.col))){
      x.abun.col <- paste0('Rel', abun.col)
    }else{
      x.abun.col <- abun.col
    }
  }
  #x.abun.col <- x.abun.col[grepl("^Rel", x.abun.col)]
  gplot.pck <- "ggplot2"
  require(gplot.pck, character.only=TRUE) %>% suppressMessages()
  if (tree.type == "otutree"){
    p1 <- ggtree(
      anno.tree,
      layout = layout,
      size = 0.3
    )
    if (!is.null(taxada)){
      p1 <- p1 +
        geom_tiplab(
          align = TRUE,
          size = 0,
          linetype = tiplab.linetype
        ) +
        geom_fruit(
          data = td_filter(!!rlang::sym("isTip")),
          geom = geom_point,
          mapping = aes(colour = !!.taxa.class),
          size = 1.5,
          offset = 0
        )    
    }
  }
  
  if (tree.type == "taxatree"){
    p1 <- suppressWarnings(
      ggtree::ggtree(
        anno.tree,
        layout = layout,
        size = 0.3
      ) +
        geom_point(
          data = td_filter(!.data$isTip),
          fill = "white",
          size = 1,
          shape = 21
        )
    )
    
    p1 <- suppressWarnings(
      p1 +
        ggtree::geom_hilight(
          data = td_filter(!!rlang::sym("nodeClass")==rlang::as_name(.taxa.class)),
          mapping = aes(
            node = !!rlang::sym("node"), 
            fill = !!rlang::sym("label")
          )
        )
    )
    # pdf("~/xfer/f.pdf"); p1; dev.off()
  }
  
  if (nsample > sample.num || group.abun){
    mapping <- aes(x=!!rlang::sym(group.nm), size = !!rlang::sym(x.abun.col), 
                   fill = !!rlang::sym(group.nm), subset = !!rlang::sym(x.abun.col) > 0) 
    n.pwidth <- .data %>%
      mp_extract_sample() %>%
      dplyr::pull(!!rlang::sym(group.nm)) %>%
      unique() %>% length()    
  }else{
    mapping <- aes(x=forcats::fct_reorder(!!rlang::sym("Sample"), !!rlang::sym(group.nm), .fun=min),
                   size = !!rlang::sym(x.abun.col), 
                   fill = !!rlang::sym(group.nm), 
                   subset = !!rlang::sym(x.abun.col) > 0
    )
    n.pwidth <- ncol(.data)
  }
  ggstar <- "ggstar"
  require(ggstar, character.only=TRUE) %>% suppressMessages()
  p2 <- suppressWarnings(
    p1 + 
      ggnewscale::new_scale_fill() +
      geom_fruit(
        data = td_unnest(!!rlang::sym(abun.col)),
        geom = geom_star,
        mapping = mapping,
        starshape = 13,
        starstroke = 0.05,
        offset = offset.abun,
        pwidth = pwidth.abun,
        grid.params = list(linetype=2)
      ) +  
      scale_size_continuous(
        name= ifelse(grepl('^Rel', abun.col), "Relative Abundance (%)", gsub("By.*", "", abun.col)),
        range = c(.5, 3),
        guide = guide_legend(override.aes = list(fill="black"))
      )
  )
  
  if (nsample > 50 || group.abun){
    p3 <- suppressWarnings(
      p2 + 
        geom_tiplab(size=tiplab.size, offset = max(p2$data$xmaxtmp, na.rm=TRUE) - 0.98*max(p2$data$x, na.rm=TRUE), 
                    align = TRUE, 
                    linetype=NA)
    )
  }else{
    p3 <- suppressWarnings(
      p2 + 
        geom_tiplab(size=tiplab.size, offset = max(p2$data$xmaxtmp, na.rm=TRUE) - 0.95*max(p2$data$x, na.rm=TRUE))
    )
  }
  
  ## Custom section, select "Significant OTU" based on input taxons
  if(!is.null(.taxons)){
    .taxons=taxons
    dat <- p3$data
    dat[,sign.field]
    rm_sig_idx <- which(!dat$label %in% .taxons)
    dat[rm_sig_idx,sign.field] <- NA
    p3$data <- dat
  }
  
  # display the LDA of significant OTU.
  title.height <- 4.4e-06 * sum(p3$data$isTip) 
  p4 <- suppressWarnings(
    p3 +
      ggnewscale::new_scale_fill() +
      geom_fruit(
        #data = td_filter(!is.na(!!rlang::sym(x.bar))),
        data = td_filter(!is.na(!!rlang::sym(sign.field))),
        geom = "geom_col",
        mapping = aes(
          x = !!rlang::sym(x.bar),
          fill = !!rlang::sym(sign.field)
        ),
        orientation = "y",
        offset = offset.effsize,
        pwidth = pwidth.effsize,
        axis.params = list(axis = "x",
                           title = x.bar.title,
                           title.height = title.height,
                           title.size = 2,
                           text.size = 1.8,
                           vjust = 1),
        grid.params = list(linetype = 2)
      )
  )
  
  # display the significant (FDR) taxonomy after kruskal.test (default)
  p5 <- suppressWarnings(
    p4 +
      ggnewscale::new_scale("size") +
      geom_point(
        data=td_filter(!is.na(!!rlang::sym(sign.field))),
        mapping = aes_string(
          size = size.mapping,
          fill = sign.field,
        ),
        shape = 21
      ) +
      scale_size_continuous(range=c(1, 3)) #+
    #scale_fill_manual(values=c("#1B9E77", "#D95F02"))
  )
  
  p6 <- suppressWarnings(
    p5 + theme(
      legend.key.height = unit(0.3, "cm"),
      legend.key.width = unit(0.3, "cm"),
      legend.spacing.y = unit(0.02, "cm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 9),
    )
  )
  
  return (suppressWarnings(p6))
}