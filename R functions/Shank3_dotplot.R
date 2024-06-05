##Make Shank3 proportion DotPlot
Shank3_dotplot <- function(mIEG){
  source("/home/hailee/10X/Final Soc Analysis/R functions/IEG pos distributions-Shank3.R")
  source("/home/hailee/10X/Final Soc Analysis/R functions/IEG pos distributions2-Shank3.R")
  row_odd <- seq_len(nrow(mIEG)) %% 2 
  Soc <- mIEG[row_odd == 0, ]
  Socmd <- md[row_odd ==0,]
  val <- as.numeric(Soc[,"N"])/as.numeric(Socmd[,"N"])
  Soc <- cbind(Soc, val)
  
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
  
  source("/home/hailee/10X/Final Soc Analysis/R functions/Chi2_multiIEG.R")
  chitab <- chi2test(mIEG, md)
  row_odd <- seq_len(nrow(chitab)) %% 2 
  chitab <- chitab[row_odd == 1, ]
  
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
  
  chitab <- replace(chitab, chitab==0, (1*10^-20))
  chitab = -1* log10(chitab)
  mIEG1 <- cbind(mIEG1, chitab)
  
  mIEG1$genotype <- factor(mIEG1$genotype, levels = c("WT", "KO"))
  dot_plot2 <- ggplot(mIEG1, aes(x=genotype, y=celltype)) +
    geom_point(aes(size = chitab, fill = as.numeric(val)), color="black", shape=21) +
    scale_size("-logP", range = c(1,10)) +
    scale_fill_gradientn(colours = viridisLite::inferno(100, direction =-1),
                         guide = guide_colorbar(ticks.colour = "black",
                                                frame.colour = "black"),
                         name = "Ratio") +
    ylab("") + xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(size=10, angle=0, color="black"),
          axis.text.y=element_blank(),
          axis.title = element_text(size=14))
  
  chitab <- chi2test(Soc, Socmd)
  row_odd <- seq_len(nrow(chitab)) %% 2 
  chitab <- chitab[row_odd == 1, ]
  
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
  Soc1$celltype <- factor(Soc1$celltype, levels = rev(levels(Sub)))
  Soc1$genotype <- rep("Soc (KO/WT)", nrow(Soc1))
  Soc1 <- Soc1[order(Soc1$celltype, Soc1$genotype), ]
  chitab <- replace(chitab, chitab==0, (1*10^-20))
  chitab = -1* log10(chitab)
  Soc1 <- cbind(Soc1, chitab)
  
  dot_plot3 <- ggplot(Soc1, aes(x=genotype, y=celltype)) +
    geom_point(aes(size = chitab, fill = as.numeric(val)), color="black", shape=21) +
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
  
  par(mfrow = c(1,3))
  library(ggpubr)
  dot_plot1 <- ggarrange(dot_plot, dot_plot2, dot_plot3, legend = "bottom",ncol = 3) 
  dot_plot1
  new_list <- list(dot_plot1, mIEG1, Soc1)
  names(new_list) <- c("dot_plot1", "mIEG1", "Soc1")
  return(new_list)
}