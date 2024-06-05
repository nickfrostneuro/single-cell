mIEG_dotplot <- function(mIEG, mIEG1){
  ##Make Dotplots to view the data
  mIEG$celltype <- factor(mIEG$celltype, levels = rev(levels(Sub)))
  dot_plot <- ggplot(mIEG, aes(x=status, y=celltype)) +
    geom_point(aes(fill = val), color="black", shape=22, size  = 10) +
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
  
  mIEG1$pvals[mIEG1$pvals > 0.05] = NA
  dot_plot2 <- ggplot(mIEG1, aes(x=status, y=celltype)) +
    geom_point(aes(size = -1 * log10(pvals), fill = log10(as.numeric(val))), color="black", shape=21) +
    scale_size("-LogP", range = c(2,12)) +
    scale_fill_gradientn(colours = viridisLite::inferno(100, direction = -1),
                         guide = guide_colorbar(ticks.colour = "black",
                                                frame.colour = "black"),
                         name = "log10(Soc/Not)") +
    ylab("") + xlab("") +
    theme_bw() +
    theme(axis.text.x = element_text(size=10, angle=0, color="black"),
          axis.text.y=element_blank(),
          axis.title = element_text(size=14))
  library(ggpubr)
  dot_plot1 <- ggarrange(dot_plot, dot_plot2, legend = "bottom",ncol = 3) 
  new_list <- list(dot_plot1, mIEG1)
  names(new_list) <- c("dot_plot1", "mIEG1")
  return(new_list)
}
