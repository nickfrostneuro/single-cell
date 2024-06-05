distribution_sem <- function(Sub, md, numb.iters= 1000, IEG.pctl = NA, IEG.mult = NA){
  ##Create distribution for 10 minute
  distrib.KO <- matrix(nrow = length(levels(Sub)), ncol = numb.iters)
  distrib.WT <- matrix(nrow = length(levels(Sub)), ncol = numb.iters)
  IEG.KO.soc <- matrix(nrow = length(levels(Sub)), ncol = numb.iters)
  IEG.WT.soc <- matrix(nrow = length(levels(Sub)), ncol = numb.iters)
  for(cellt in 1:length(levels(Sub))){
    print(cellt)
    rm(md.CT, CT.subset)
    ##Subset one celltype
    CT.subset <- subset(Sub, idents =levels(Sub)[cellt])
    md.CT <- md[(md$celltype == levels(Sub)[cellt]),]
    WT.Not <- subset(CT.subset, subset = condition == "Not WT")
    WT.Not[["CellName"]] <- colnames(WT.Not)
    WT.Soc <- subset(CT.subset, subset = condition == "Soc WT")
    WT.Soc[["CellName"]] <- colnames(WT.Soc)
    KO.Not <- subset(CT.subset, subset = condition == "Not KO")
    KO.Not[["CellName"]] <- colnames(KO.Not)
    KO.Soc <- subset(CT.subset, subset = condition == "Soc KO")
    KO.Soc[["CellName"]] <- colnames(KO.Soc)
    
    for(d in 1:numb.iters){
      ##Sample the number of cells for not social and find the proportion IEG+ for WT
      md.sub <- md.CT[(md.CT$condition == "Not WT"),]
      set.seed(as.integer(Sys.time()) %% 10000)
      sampled.cells <- sample(WT.Not$CellName, size = md.sub$N, replace = TRUE)
      obj.sub <- subset(WT.Not, subset = CellName %in% sampled.cells)
      
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
      IEG.pos.WTnot <- sum(aa)/ncol(obj.sub)
      
      rm(aa, multIEG, name, mat, obj.sub, sampled.cells, md.sub,i)
      ##Sample the number of cells for WT social condition and find the proportion IEG+
      md.sub <- md.CT[(md.CT$condition == "Soc WT"),]
      set.seed(as.integer(Sys.time()) %% 10000)
      sampled.cells <- sample(WT.Soc$CellName, size = md.sub$N, replace = TRUE)
      obj.sub <- subset(WT.Soc, subset = CellName %in% sampled.cells)
      
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
      IEG.pos.WTsoc <- sum(aa)/ncol(obj.sub)

      rm(aa, multIEG, name, mat, obj.sub, sampled.cells, md.sub,i)
      
      ##Sample the number of cells for not social and find the proportion IEG+ for KO
      md.sub <- md.CT[(md.CT$condition == "Not KO"),]
      set.seed(as.integer(Sys.time()) %% 10000)
      a <- md.sub$N
      sampled.cells <- sample(KO.Not$CellName, size = md.sub$N, replace = TRUE)
      obj.sub <- subset(KO.Not, subset = CellName %in% sampled.cells)
      
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
      IEG.pos.KOnot <- sum(aa)/ncol(obj.sub)
      
      rm(aa, multIEG, name, mat, obj.sub, sampled.cells, md.sub,i)
      
      
      ##Sample the number of cells for WT social condition and find the proportion IEG+
      md.sub <- md.CT[(md.CT$condition == "Soc KO"),]
      set.seed(as.integer(Sys.time()) %% 10000)
      sampled.cells <- sample(KO.Soc$CellName, size = md.sub$N, replace = TRUE)
      obj.sub <- subset(KO.Soc, subset = CellName %in% sampled.cells)
      
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
      IEG.pos.KOsoc <- sum(aa)/ncol(obj.sub)
      
      rm(aa, multIEG, name, mat, obj.sub, sampled.cells, md.sub,i)
      
      ##Calculate a proportion for each of the social conditions
      WT <- IEG.pos.WTsoc/IEG.pos.WTnot
      KO <- IEG.pos.KOsoc/IEG.pos.KOnot
      
      #Save to use for a distribution where row is celltype and column is the iteration
      IEG.WT.soc[cellt,d] <- IEG.pos.WTsoc
      IEG.KO.soc[cellt,d] <- IEG.pos.KOsoc
      distrib.WT[cellt,d] <- WT
      distrib.KO[cellt,d] <- KO
    }
  }
  rownames(distrib.WT) <- levels(Sub)
  rownames(distrib.KO) <- levels(Sub)
  
  distrib.KOWT <- distrib.KO/distrib.WT
  distrib.WTKO <- distrib.WT/distrib.KO
  
  distrib.KOWT[sapply(distrib.KOWT,is.infinite)] <- NA
  distrib.WTKO[sapply(distrib.WTKO,is.infinite)] <- NA
  
  stats.KOWT <- matrix(ncol = 3, nrow = nrow(distrib.KOWT))
  rm(n)
  for(n in 1:nrow(distrib.KOWT)){
    m <- mean(as.numeric(distrib.KOWT[n,]), na.rm = TRUE)
    sdev <- sd(as.numeric(distrib.KOWT[n,]), na.rm = TRUE)
    md.CT <- md[(md$celltype == rownames(distrib.KOWT)[n]),]
    sem <- sdev/sqrt(mean(md.CT$N))
    stats.KOWT[n,1] <- m
    stats.KOWT[n,2] <- sdev
    stats.KOWT[n,3] <- sem
  }
  rownames(stats.KOWT) <- rownames(distrib.KOWT)
  colnames(stats.KOWT) <- c("mean", "sd", "sem")
  
  stats.WTKO <- matrix(ncol = 3, nrow = nrow(distrib.WTKO))
  rm(n)
  for(n in 1:nrow(distrib.WTKO)){
    m <- mean(as.numeric(distrib.WTKO[n,]), na.rm = TRUE)
    sdev <- sd(as.numeric(distrib.WTKO[n,]), na.rm = TRUE)
    md.CT <- md[(md$celltype == rownames(distrib.WTKO)[n]),]
    sem <- sdev/sqrt(mean(md.CT$N))
    stats.WTKO[n,1] <- m
    stats.WTKO[n,2] <- sdev
    stats.WTKO[n,3] <- sem
  }
  rownames(stats.WTKO) <- rownames(distrib.WTKO)
  colnames(stats.WTKO) <- c("mean", "sd", "sem")
 
  stats <- list(stats.KOWT, stats.WTKO, distrib.KO, distrib.WT)
  names(stats) <- c("stats.KOWT", "stats.WTKO","distrib.KO", "distrib.WT")
  return(stats)
}