DEAnalysisMAST <- function(scdata,id,path){
  pseudo.count = 0.1
  data.used.log2   <- log2(scdata+pseudo.count)
  colnames(data.used.log2)<-make.unique(colnames(data.used.log2))
  diff.cutoff=0.5
  for(i in unique(id)){
    ## 2024-Jan-14:  Added a simple check to see if the files already exist and skip 
    # regenerating it if it's already there
    if(file.exists(paste(path,"/",i,"_MIST.RData", sep=""))){
      print(paste0(i, " exists.. skipping."))
      next
    }
    print(paste0("Analyzing: ", i))
    cells.symbol.list2     = colnames(data.used.log2)[which(id==i)]
    cells.coord.list2      = match(cells.symbol.list2, colnames(data.used.log2))                          
    cells.symbol.list1     = colnames(data.used.log2)[which(id != i)]
    cells.coord.list1      = match(cells.symbol.list1, colnames(data.used.log2))   
    data.used.log2.ordered  = cbind(data.used.log2[,cells.coord.list1], data.used.log2[,cells.coord.list2])
    group.v <- c(rep(0,length(cells.coord.list1)), rep(1, length(cells.coord.list2)))
    #ouput
    log2.stat.result <- stat.log2(data.used.log2.ordered, group.v, pseudo.count)
    Auc <- m.auc(data.used.log2.ordered, group.v)
    bigtable <- data.frame(cbind(log2.stat.result, Auc))
    
    DE <- bigtable[bigtable$log2_fc >diff.cutoff,] 
    dim(DE)
    if(dim(DE)[1]>1){
      data.1                 = data.used.log2[,cells.coord.list1]
      data.2                 = data.used.log2[,cells.coord.list2]
      genes.list = rownames(DE)
      log2fold_change        = cbind(genes.list, DE$log2_fc)
      colnames(log2fold_change) = c("gene.name", "log2fold_change")
      counts  = as.data.frame(cbind( data.1[genes.list,], data.2[genes.list,] ))
      groups  = c(rep("Cluster_Other", length(cells.coord.list1) ), rep(i, length(cells.coord.list2) ) )
      groups  = as.character(groups)
      data_for_MIST <- as.data.frame(cbind(rep(rownames(counts), dim(counts)[2]), melt(counts),rep(groups, each = dim(counts)[1]), rep(1, dim(counts)[1] * dim(counts)[2]) ))
      colnames(data_for_MIST) = c("Gene", "Subject.ID", "Et", "Population", "Number.of.Cells")
      vbeta = data_for_MIST
      vbeta.fa <-FromFlatDF(vbeta, idvars=c("Subject.ID"),
                            primerid='Gene', measurement='Et', ncells='Number.of.Cells',
                            geneid="Gene",  cellvars=c('Number.of.Cells', 'Population'),
                            phenovars=c('Population'), id='vbeta all')
      vbeta.1 <- subset(vbeta.fa,Number.of.Cells==1)
      # .3 MAST 
      head(colData(vbeta.1))
      zlm.output <- zlm(~ Population, vbeta.1, method='bayesglm', ebayes=TRUE)
      show(zlm.output)
      coefAndCI <- summary(zlm.output, logFC=TRUE)
      zlm.lr <- lrTest(zlm.output, 'Population')
      zlm.lr_pvalue <- melt(zlm.lr[,,'Pr(>Chisq)'])
      zlm.lr_pvalue <- zlm.lr_pvalue[which(zlm.lr_pvalue$test.type == 'hurdle'),]
      
      
      
      lrTest.table <-  merge(zlm.lr_pvalue, DE, by.x = "primerid", by.y = "row.names")
      colnames(lrTest.table) <- c("Gene", "test.type", "p_value", paste("log2.mean.", "Cluster_Other", sep=""), paste("log2.mean.",i,sep=""), "log2fold_change", "Auc")
      cluster_lrTest.table <- lrTest.table[rev(order(lrTest.table$Auc)),]
      
      #. 4 save results
      write.csv(cluster_lrTest.table, file=paste(path,"/",i,"_lrTest.csv", sep=""))
      save(cluster_lrTest.table, file=paste(path,"/",i,"_MIST.RData", sep=""))
    }
  }
}