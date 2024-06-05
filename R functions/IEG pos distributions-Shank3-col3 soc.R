distribution_IEG.col3 <- function(Sub, Socmd, numb.iters= 1000, IEG.pctl = 0.95, IEG.mult = 3){
  ##Create distribution for 10 minute
  distrib.Soc <- matrix(nrow = length(levels(Sub)), ncol = numb.iters)
  for(cellt in 1:length(levels(Sub))){
    print(cellt)
    rm(md.CT, CT.subset)
    ##Subset one celltype
    CT.subset <- subset(Sub, idents =levels(Sub)[cellt])
    CT.subset <- subset(CT.subset, subset = status == "Soc")
    md.CT <- Socmd[(Socmd$celltype == levels(Sub)[cellt]),]
    
    for(d in 1:numb.iters){
      rm(Soc, IEG.pos.KOsoc, IEG.pos.WTsoc)
      ##Sample the number of cells for WT social condition and find the proportion IEG+
      md.sub <- md.CT[(md.CT$condition == "Soc WT"),]
      set.seed(as.integer(Sys.time()) %% 10000)
      sampled.cells <- sample(colnames(CT.subset), size = md.sub$N, replace = F)
      obj.sub <- subset(CT.subset, cells = sampled.cells)
      
      mat <- matrix(ncol=length(genes), nrow= ncol(obj.sub))
      for(i in 1:length(genes)){
        tmp <- GetAssayData(object=obj.sub, slot = "data")[genes[i], ]
        q <- (tmp > quantile(tmp, IEG.pctl))
        q1 <- as.integer(q)
        mat[,i] <- q
      }
      name <- names(q)
      multIEG <- rowSums(mat)
      aa <- as.numeric(as.integer(multIEG>IEG.mult))
      if(sum(aa) == 0){
        IEG.pos.WTsoc <- NA
      } else{
        IEG.pos.WTsoc <- sum(aa)/ncol(obj.sub)
      }
      rm(aa, multIEG, name, mat, obj.sub, sampled.cells, md.sub,i)
      
      ##Sample the number of cells for WT social condition and find the proportion IEG+
      md.sub <- md.CT[(md.CT$condition == "Soc KO"),]
      set.seed(as.integer(Sys.time()) %% 10000)
      sampled.cells <- sample(colnames(CT.subset), size = md.sub$N, replace = F)
      obj.sub <- subset(CT.subset, cells = sampled.cells)
      
      mat <- matrix(ncol=length(genes), nrow= ncol(obj.sub))
      for(i in 1:length(genes)){
        tmp <- GetAssayData(object=obj.sub, slot = "data")[genes[i], ]
        q <- (tmp > quantile(tmp, IEG.pctl))
        q1 <- as.integer(q)
        mat[,i] <- q
      }
      name <- names(q)
      multIEG <- rowSums(mat)
      aa <- as.numeric(as.integer(multIEG>IEG.mult))
      if(sum(aa) == 0){
        IEG.pos.KOsoc <- NA
      } else{
        IEG.pos.KOsoc <- sum(aa)/ncol(obj.sub)
      }
      rm(aa, multIEG, name, mat, obj.sub, sampled.cells, md.sub,i)
      
      ##Calculate a proportion for each of the social conditions
      Soc <- IEG.pos.KOsoc/IEG.pos.WTsoc
      
      #Save to use for a distribution where row is celltype and column is the iteration
      distrib.Soc[cellt,d] <- Soc
    }
  }
  rownames(distrib.Soc) <- levels(Sub)
  return(distrib.Soc)
}