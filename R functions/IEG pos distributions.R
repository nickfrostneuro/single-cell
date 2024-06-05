distribution_IEG <- function(Sub, md, numb.iters= 500, IEG.pctl = 0.95, IEG.mult = 3){
  ##Create distribution for 10 minute
  distrib.10min.over.not <- matrix(nrow = length(levels(Sub)), ncol = numb.iters)
  distrib.35min.over.not <- matrix(nrow = length(levels(Sub)), ncol = numb.iters)
  for(cellt in 1:length(levels(Sub))){
    print(cellt)
    rm(md.CT, CT.subset)
    ##Subset one celltype
    CT.subset <- subset(Sub, idents =levels(Sub)[cellt])
    md.CT <- md[(md$celltype == levels(Sub)[cellt]),]
    
    for(d in 1:numb.iters){
      ##Sample the number of cells for not social and find the proportion IEG+
      md.sub <- md.CT[(md.CT$status == "Not"),]
      set.seed(as.integer(Sys.time()) %% 10000)
      sampled.cells <- sample(colnames(CT.subset), size = md.sub$N , replace = F)
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
        #Save to use for a distribution where row is cell type and column is the iteration
        IEG.pos.not <- NA
      } else{
        IEG.pos.not <- sum(aa)/ncol(obj.sub)
      }
      rm(aa, multIEG, name, mat, obj.sub, sampled.cells, md.sub,i)
      ##Sample the number of cells for 10 min social and find the proportion IEG+
      md.sub <- md.CT[(md.CT$status == "ShS"),]
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
        IEG.pos.10min <- NA
      } else{
        IEG.pos.10min <- sum(aa)/ncol(obj.sub)
      }
      rm(aa, multIEG, name, mat, obj.sub, sampled.cells, md.sub,i)
      ##Sample the number of cells for 35 min social and find the proportion IEG+
      md.sub <- md.CT[(md.CT$status == "Soc"),]
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
        IEG.pos.35min <- NA
      } else{
        IEG.pos.35min <- sum(aa)/ncol(obj.sub)
      }
      rm(aa, multIEG, name, mat, obj.sub, sampled.cells, md.sub,i)
      
      ##Calculate a proportion for each of the social conditions
      Early <- IEG.pos.10min/IEG.pos.not
      Late <- IEG.pos.35min/IEG.pos.not
      
      #Save to use for a distribution where row is celltype and column is the iteration
      distrib.10min.over.not[cellt,d] <- Early
      distrib.35min.over.not[cellt,d] <- Late
    }
  }
  rownames(distrib.10min.over.not) <- levels(Sub)
  rownames(distrib.35min.over.not) <- levels(Sub)
  distrib <- list(distrib.10min.over.not, distrib.35min.over.not)
  names(distrib) <- c("distrib.10min", "distrib.35min")
  return(distrib)
}