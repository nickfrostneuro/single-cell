##Make Shank3 proportion DotPlot
Shank3_dotplot1 <- function(mIEG, md, Sub, distrib.KO, distrib.WT, distrib.Soc){
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
    pval <- length(distrib.KO[i,][(distrib.KO[i,] > as.numeric(mIEG1_KO$val[i]))])/length(distrib.KO[i,])
    pval[pval == 0] =0.001
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
  
  mat <- matrix(nrow=nrow(Soc), ncol = 1)
  for(i in 1:nrow(Soc)){
    if((i %% 2) == 0)next
    l <- Soc$val[(i+1)]/Soc$val[i]
    mat[i,] = l
  }
  Soc$genotype <- substr(Soc$condition, 5, 6)
  Soc1 <- cbind(Soc$genotype, as.character(Soc$celltype), as.numeric(mat))
  row_odd <- seq_len(nrow(Soc1)) %% 2 
  Soc1 <- Soc1[row_odd == 1, ]
  colnames(Soc1) <- c("genotype", "celltype", "val")
  Soc1 <- as.data.frame(Soc1)
  Soc1$celltype <- factor(Soc1$celltype, levels = levels(Sub))
  Soc1$genotype <- rep("Soc (KO/WT)", nrow(Soc1))
  Soc1 = Soc1[nrow(Soc1):1,]
  
  pval.Soc <- matrix(ncol = 1, nrow = nrow(distrib.Soc))
  for(i in 1:nrow(distrib.Soc)){
    rm(pval)
    dens <- density(na.omit(distrib.Soc[i,]))
    plot(dens,xlim = c(0,6), main = paste0(rownames(distrib.Soc)[i], " Soc"))
    abline(v = as.numeric(Soc1$val[i]), col = "red")
    pval <- length(distrib.Soc[i,][(distrib.Soc[i,] > as.numeric(Soc1$val[i]))])/length(distrib.Soc[i,])
    pval[pval == 0] =0.001
    pval.Soc[i,] = pval
  }
  rownames(pval.Soc) <- rownames(distrib.Soc)
  pvals <- pval.Soc
  Soc1 <- cbind(Soc1, pvals)
  Soc <- Soc1
  
  Soc1 <- Soc1[nrow(Soc1):1,]
  Soc1$pvals[Soc1$pvals > 0.05] = NA
  Soc1$celltype <- factor(Soc1$celltype, levels = rev(levels(Sub)))
  Soc1 <- Soc1[order(Soc1$celltype), ]
  dot_plot3 <- ggplot(Soc1, aes(x=genotype, y=celltype)) +
    geom_point(aes(size = -1*log10(pvals), fill = as.numeric(val)), color="black", shape=21) +
    scale_size("-logP", range = c(2,12)) +
    scale_fill_gradientn(colours = viridisLite::inferno(100, direction =-1),
                         guide = guide_colorbar(ticks.colour = "black",
                                                frame.colour = "black"),
                         name = "Ratio") +
    ylab("") + xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(size=10, angle=0, color="black"),
          axis.text.y=element_blank(),
          axis.title = element_text(size=14))
  
  library(ggpubr)
  dot_plot1 <- ggarrange(dot_plot, dot_plot2, dot_plot3, legend = "bottom",ncol = 3) 
  dot_plot1
  
  new_list <- list(dot_plot1, mIEG1, Soc)
  names(new_list) <- c("dot_plot1", "mIEG1", "Soc1")
  return(new_list)
}