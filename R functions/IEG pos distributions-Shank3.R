distribution_IEG.col2 <- function(Sub, md, numb.iters= 1000, IEG.pctl = 0.95, IEG.mult = 3){
  ##Create distribution for 10 minute
  distrib.KO <- matrix(nrow = length(levels(Sub)), ncol = numb.iters)
  distrib.WT <- matrix(nrow = length(levels(Sub)), ncol = numb.iters)
  for(cellt in 1:length(levels(Sub))){
    print(cellt)
    rm(md.CT, CT.subset)
    ##Subset one celltype
    CT.subset <- subset(Sub, idents =levels(Sub)[cellt])
    md.CT <- md[(md$celltype == levels(Sub)[cellt]),]
    WT.subset <- subset(CT.subset, subset = genotype == "WT")
    KO.subset <- subset(CT.subset, subset = genotype == "KO")
    
    for(d in 1:numb.iters){
      ##Sample the number of cells for not social and find the proportion IEG+ for WT
      md.sub <- md.CT[(md.CT$condition == "Not WT"),]
      set.seed(as.integer(Sys.time()) %% 10000)
      sampled.cells <- sample(colnames(WT.subset), size = md.sub$N , replace = F)
      obj.sub <- subset(WT.subset, cells = sampled.cells)
      
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
        #Save to use for a distribution where row is cell type and column is the iteration
        IEG.pos.WTnot <- NA
      } else{
        IEG.pos.WTnot <- sum(aa)/ncol(obj.sub)
      }
      rm(aa, multIEG, name, mat, obj.sub, sampled.cells, md.sub,i)
      ##Sample the number of cells for WT social condition and find the proportion IEG+
      md.sub <- md.CT[(md.CT$condition == "Soc WT"),]
      set.seed(as.integer(Sys.time()) %% 10000)
      sampled.cells <- sample(colnames(WT.subset), size = md.sub$N, replace = F)
      obj.sub <- subset(WT.subset, cells = sampled.cells)
      
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
      
      ##Sample the number of cells for not social and find the proportion IEG+ for KO
      md.sub <- md.CT[(md.CT$condition == "Not KO"),]
      set.seed(as.integer(Sys.time()) %% 10000)
      a <- md.sub$N
      if(a == 0){
        IEG.pos.KOnot <- 0
        }
      else{sampled.cells <- sample(colnames(KO.subset), size = md.sub$N , replace = F)
      obj.sub <- subset(KO.subset, cells = sampled.cells)
      
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
        #Save to use for a distribution where row is cell type and column is the iteration
        IEG.pos.KOnot <- NA
      } else{
        IEG.pos.KOnot <- sum(aa)/ncol(obj.sub)
      }
      rm(aa, multIEG, name, mat, obj.sub, sampled.cells, md.sub,i)
      }
      
      ##Sample the number of cells for WT social condition and find the proportion IEG+
      md.sub <- md.CT[(md.CT$condition == "Soc KO"),]
      set.seed(as.integer(Sys.time()) %% 10000)
      sampled.cells <- sample(colnames(KO.subset), size = md.sub$N, replace = F)
      obj.sub <- subset(KO.subset, cells = sampled.cells)
      
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
      WT <- IEG.pos.WTsoc/IEG.pos.WTnot
      KO <- IEG.pos.KOsoc/IEG.pos.KOnot
      
      #Save to use for a distribution where row is celltype and column is the iteration
      distrib.WT[cellt,d] <- WT
      distrib.KO[cellt,d] <- KO
    }
  }
  rownames(distrib.WT) <- levels(Sub)
  rownames(distrib.KO) <- levels(Sub)
  distrib <- list(distrib.WT, distrib.KO)
  names(distrib) <- c("distrib.WT", "distrib.KO")
  return(distrib)
}