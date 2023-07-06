pathways <- list("Pathway_x"=paste0("Gene", c(1:5)),
                 "Pahtway_y"=paste0("Gene", c(4:10)),
                 "Pathway_z"=paste0("Gene", c(1, 5, 10:60)))
mutated_genes <- paste0("Gene", c(2,2,2,3,3,5,50:70))

pathway_df <- lapply(names(pathways), function(p){
  pdf <- data.frame("Pathway"=p,
                    "Genes"=pathways[[p]])
  pdf %>% 
    mutate("Mutated"=ifelse(Genes %in% mutated_genes, 1, 0))
}) %>% 
  do.call(rbind, .)

cont_tbl <- with(pathway_df, table(Pathway, Mutated))
chisq_val <- chisq.test(cont_tbl)
chisq_posthoc <- fifer::chisq.post.hoc(cont_tbl) %>%
  arrange(adj.p)
