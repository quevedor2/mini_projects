library(dplyr)
pc_mat <- as.data.frame(matrix(rnorm(n=1000), nrow=10)) %>%
  magrittr::set_colnames(., paste0("PC", c(1:100))) %>%
  magrittr::set_rownames(., paste0("Sample_", c(1:10)))
meta_mat <- data.frame("Group"=c(rep("A",5), rep("B", 5)),
                       "Group2"=sample(c("X", "Y"), size=10, replace=T)) %>%
  magrittr::set_rownames(., paste0("Sample_", c(1:10)))
groups <- c("Group", "Group2")

pcs_meta <- lapply(groups, function(group){
  print(group)
  if(length(unique(meta_mat[,group]))<=1) return(NULL)
  
  lm_stat <- linearModelEigen(
    eigen=t(pc_mat), 
    model.matrix(as.formula(paste0("~ meta_mat$", group)))
  )
  
  gg_df <- pc_mat[,lm_stat$module[1:2]] %>%
    tibble::rownames_to_column('sample') %>%
    left_join(meta_mat[,group,drop=F] %>% 
                tibble::rownames_to_column('sample'),
              by=c("sample"="sample"))
  
  ggp <- ggplot(gg_df, aes_string(x=lm_stat$module[1], y=lm_stat$module[2], 
                                  fill=group, color=group)) +
    geom_point() + 
    xlab(paste0(lm_stat$module[1], " (", round(lm_stat$adj.P.Val[1],4), ")")) +
    ylab(paste0(lm_stat$module[2], " (", round(lm_stat$adj.P.Val[2],4), ")")) +
    theme_classic() +
    theme(legend.position='none')
  
  ggp_pcx <- ggplot(gg_df, aes_string(x=lm_stat$module[1], fill=group)) +
    geom_density(alpha=0.5) + 
    theme_classic() +
    ggtitle(group) +
    xlab("") + ylab(lm_stat$module[1]) +
    theme(legend.position = "none",
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  ggp_pcy <- ggplot(gg_df, aes_string(x=lm_stat$module[2], fill=group)) +
    geom_density(alpha=0.5) + 
    theme_classic() +
    coord_flip() +
    xlab(lm_stat$module[2]) + ylab("") +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  ggp_pg <- plot_grid(ggp_pcx, NULL,ggp, ggp_pcy, labels = "AUTO", 
                      rel_widths = c(4, 2),
                      rel_heights = c(1,4))
  
  return(list("stat"=lm_stat, "gg"=ggp_pg))
})
