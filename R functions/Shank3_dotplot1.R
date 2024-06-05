##Make Shank3 proportion DotPlot
Shank3_dotplot1 <- function(mIEG, md, Sub, distrib.KO, distrib.WT){
  val <- as.numeric(mIEG[,"N"])/as.numeric(md[,"N"])
  mIEG <- cbind(mIEG, val)
  
  dot_plot <- ggplot(mIEG, aes(x=condition, y=celltype)) +
    geom_point(aes(fill = val), color="black", shape=22, size  = 8) +
    #scale_size("Proportion", range = c(0,6)) +
    scale_fill_gradientn(colours = viridisLite::inferno(100, direction = -1),
                         guide = guide_colorbar(ticks.colour = "black",
                                                frame.colour = "black"),
                         name = "Proportion") +
    ylab("") + xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
          axis.text.y = element_text(size=12, color="black"),
          axis.title = element_text(size=14))
  
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
  KO.WT <- cbind(as.character(mIEG1$celltype), as.numeric(mat))
  row_odd <- seq_len(nrow(KO.WT)) %% 2 
  KO.WT <- KO.WT[row_odd == 1, ]
  colnames(KO.WT) <- c("celltype", "val")
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
  mIEG1$pvals[mIEG1$pvals > 0.05] = NA
  dot_plot2 <- ggplot(mIEG1, aes(x=genotype, y=celltype)) +
    geom_point(aes(size = -1*log10(pvals), fill = as.numeric(val)), color="black", shape=21) +
    scale_size("-logP", range = c(1,10)) +
    scale_fill_gradientn(colours = viridisLite::inferno(100, direction =-1),
                         guide = guide_colorbar(ticks.colour = "black",
                                                frame.colour = "black"),
                         name = "Soc/Not") +
    ylab("") + xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(size=10, angle=0, color="black"),
          axis.text.y=element_blank(),
          axis.title = element_text(size=14))
  
  KO.WT.save <- KO.WT
  KO.WT$celltype <- factor(KO.WT$celltype, levels = rev(levels(Sub)))
  KO.WT$pval[KO.WT$pval > 0.05] = NA
  dot_plot3 <- ggplot(KO.WT, aes(x=condition, y=celltype)) +
    geom_point(aes(size = -1*log10(pval), fill = as.numeric(val)), color="black", shape=21) +
    scale_size("-logP", range = c(2,12)) +
    scale_fill_gradientn(colours = viridisLite::inferno(100, direction =-1),
                         guide = guide_colorbar(ticks.colour = "black",
                                                frame.colour = "black"),
                         name = "KO/WT") +
    ylab("") + xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(size=10, angle=0, color="black"),
          axis.text.y=element_blank(),
          axis.title = element_text(size=14))
    
  library(ggpubr)
  dot_plot1 <- ggarrange(dot_plot, dot_plot2, dot_plot3, legend = "bottom",ncol = 3) 
  dot_plot1
  
  new_list <- list(dot_plot1, mIEG1.save, KO.WT.save)
  names(new_list) <- c("dot_plot1", "mIEG1", "KO.WT")
  return(new_list)
}