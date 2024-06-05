##Make Shank3 proportion Barplot
#mIEG is table of condition, cell type, N IEG+
#md is table of condition, cell type, N cells
#Sub is the neuron subset of the original "CellFile" Seurat object
#distrib.KO is the distribution of shuffles soc/not values for KO calculated by the "distribution_IEG.col2" function
#distrib.WT is the distribution of shuffles soc/not values for WT calculated by the "distribution_IEG.col2" function
#object is a data frame with row names = cell types in same order as rows of distrib.KO and 3 columns (mean, sd, and sem of the shuffled data bootstrapping)

Shank3_barplot <- function(mIEG, md, Sub, distrib.KO, distrib.WT, object){
  val <- as.numeric(mIEG[,"N"])/as.numeric(md[,"N"])
  mIEG <- cbind(mIEG, val)

  mat <- matrix(nrow=nrow(mIEG), ncol = 1)
  for(i in 1:nrow(mIEG)){
    if((i %% 2) == 0)next
    l <- mIEG$val[(i+1)]/mIEG$val[i]
    mat[i,] = l
  }
  mIEG1 <- cbind(as.character(mIEG$condition), as.character(mIEG$celltype), as.numeric(mat))
  row_odd <- seq_len(nrow(mIEG1)) %% 2 
  mIEG1 <- mIEG1[row_odd == 1, ]
  colnames(mIEG1) <- c("status", "celltype", "val")
  mIEG1 <- as.data.frame(mIEG1)
  
  mIEG1$celltype <- factor(mIEG1$celltype, levels = rev(levels(Sub)))
  mIEG1$genotype <- rep(c("WT", "KO"), nrow(mIEG1)/2)
  mIEG1$genotype <- factor(mIEG1$genotype, levels = c("WT", "KO"))
  mIEG1 <- mIEG1[order(mIEG1$genotype), ]
  
  mIEG1_KO <- mIEG1[(length(levels(Sub))+1):nrow(mIEG1),]
  mIEG1_KO$celltype <- factor(mIEG1_KO$celltype, levels = levels(Sub))
  mIEG1_KO <- mIEG1_KO[order(mIEG1_KO$celltype), ]
  mIEG1_WT <- mIEG1[1:length(levels(Sub)),]
  mIEG1_WT$celltype <- factor(mIEG1_WT$celltype, levels = levels(Sub))
  mIEG1_WT <- mIEG1_WT[order(mIEG1_WT$celltype), ]
  
  #Plot distributions by cell type with red line to show the real proportion increase value
  pval.KO <- matrix(ncol = 1, nrow = nrow(distrib.KO))
  pval.WT <- matrix(ncol = 1, nrow = nrow(distrib.WT))
  for(i in 1:nrow(distrib.KO)){
    rm(pval, pval1)
    #rm(tfmin, pval, tenmin, pval1)
    #tfmin <- density(na.omit(as.numeric(distrib.KO[i,])))
    #plot(tfmin, xlim = c(0,max(as.numeric(mIEG1_KO$val)[is.finite(as.numeric(mIEG1_KO$val))])), main = paste0(rownames(distrib.KO)[i], " KO"))
    #abline(v = as.numeric(mIEG1_KO$val[i]), col = "red")
    pval <- length(distrib.KO[i,][(distrib.KO[i,] > as.numeric(mIEG1_KO$val[i]))])/length(distrib.KO[i,])
    pval[pval == 0] =0.001
    #tfmin <- density(na.omit(as.numeric(distrib.WT[i,])))
    #plot(tfmin, xlim = c(0,max(as.numeric(mIEG1_WT$val)[is.finite(as.numeric(mIEG1_WT$val))])), main = paste0(rownames(distrib.WT)[i], " WT"))
    #abline(v = as.numeric(mIEG1_WT$val[i]), col = "red")
    pval.KO[i,] = pval
    pval1 <- length(distrib.WT[i,][(distrib.WT[i,] > as.numeric(mIEG1_WT$val[i]))])/length(distrib.WT[i,])
    pval1[pval1 == 0] =0.001
    pval.WT[i,] = pval1
  }
  rownames(pval.KO) <- rownames(distrib.KO)
  rownames(pval.WT) <- rownames(distrib.WT)
  pval.KO <- rev(pval.KO)
  pval.WT <- rev(pval.WT)
  pvals <- c(pval.WT,pval.KO)
  mIEG1 <- cbind(mIEG1, pvals)
  
  mIEG1 <- mIEG1[order(mIEG1$celltype), ]
  
  #Calculate what the value of KO/WT from col2 would be for shuffled
  #distrib <- distrib.KO/distrib.WT
  rm(d)
  distrib <- list()
  for(d in 1:ncol(distrib.KO)){
    distrib.WT <- distrib.WT[,c(2:ncol(distrib.WT),1)]
    dist <- distrib.KO/distrib.WT
    distrib[[d]] <- dist
  }
  distrib <- do.call(cbind, distrib)
  rownames(distrib) <- rownames(distrib.KO)
  
  #Calculate real values of KO col2/WT col2
  mat <- matrix(nrow=nrow(mIEG1), ncol = 1)
  for(i in 1:nrow(mIEG1)){
    if((i %% 2) == 0)next
    l <- as.numeric(mIEG1$val[(i+1)])/as.numeric(mIEG1$val[i])
    mat[i,] = l
  }
  KO.WT <- cbind(as.character(mIEG1$celltype), as.numeric(mat), rep(rev(object[,2]), each =2))
  row_odd <- seq_len(nrow(KO.WT)) %% 2 
  KO.WT <- KO.WT[row_odd == 1, ]
  colnames(KO.WT) <- c("celltype", "val", "sd")
  KO.WT <- as.data.frame(KO.WT)
  
  KO.WT <- KO.WT[nrow(KO.WT):1,]
  
  #Plot distributions by cell type with red line to show the real proportion increase value
  pval <- matrix(ncol = 1, nrow = nrow(distrib))
  for(i in 1:nrow(distrib)){
    #tfmin <- density(na.omit(as.numeric(distrib[i,])))
    #plot(tfmin, xlim = c(0,max(as.numeric(KO.WT$val)[is.finite(as.numeric(KO.WT$val))])), main = paste0(rownames(distrib.WT)[i], " KO/WT"))
    #abline(v = as.numeric(KO.WT$val[i]), col = "red")
    pval.tmp <- length(distrib[i,][(distrib[i,] > as.numeric(KO.WT$val[i]))])/length(distrib[i,])
    pval.tmp[pval.tmp == 0] =0.001
    pval[i,] = pval.tmp
  }
  rownames(pval) <- rownames(levels(Sub))
  KO.WT <- cbind(KO.WT, pval)
  KO.WT$condition <- rep("KO/WT", nrow(KO.WT))
  KO.WT <- KO.WT[nrow(KO.WT):1,]
  
  mIEG1.save <- mIEG1
  KO.WT.save <- KO.WT
  KO.WT$celltype <- factor(KO.WT$celltype, levels = rev(levels(Sub)))
  
  new_list <- list(mIEG1.save, KO.WT.save)
  names(new_list) <- c("mIEG1", "KO.WT")
  return(new_list)
}