###Figure2 - Heatmaps
## Hailee Walker - 07142023
#Makes heatmaps and outputs csvs of data used to make figures in Fig. 2
#Adapted from Seurat vignette protocols
rm(list = ls())
graphics.off()

library(Seurat)
library(readr)
library(ComplexHeatmap)
library(dplyr)
library(magrittr)
library(viridis)
library(circlize)
library(rlist)

## Figure 2A PFC ## _________________________________________________________________________________________________________
##Rerun but split by mouse
##Create heatmaps for PFC
load("~/10X/Final Soc Analysis/PFC_All_analysis/Comb_PFC.Rdata")

Sub <- subset(CellFile, subset= genotype == "WT")
Sub <- subset(Sub, idents = c("Astro", "Oligo", "OPC","Microglia","Endothelial", "VLMC","Pericyte"), invert = TRUE)
#Full IEG panel
#genes <- (as.list(na.omit(read_csv("/home/hailee/10X/IEG list (Hong,W. 2017).csv", col_names = FALSE)[,1])))[[1]]
#2 fold or greater IEG panel
genes <- c("Fos", "Arc","Junb", "Egr4", "Btg2", "Dusp1", "Ccn1", "Fosb", "Egr3", "Ier5", "Egr2", "Fosl2", "Nr4a3", "Egr1", "Dusp6", "Per1", "Jun", "Nr4a2", "Ier2", "Gadd45g", "Cebpb", 
           "Bhlhe40", "Hes1", "Sik1", "Ccn2", "Ier3", "Gadd45b", "Nrn1", "Npas4", "F3", "Trib1", "Nr4a1", "Klf10", "Plk2", "Nfkbid", "Bdnf", "Zfp36l2", "Nfkbia", "Gem", "Ptgs2", "Klf2", 
           "Zfp36l1", "Apold1", "Dusp2", "Pim1", "Maff", "Ifit3", "Inhba", "Arhgef3", "Nfib", "Ifit2") #PFC + Cereb genes

#Pull out the proportion of the cluster that expresses each gene from the calculations made for making dotplots (Can be done many ways, but this was easiest)
plot <- DotPlot(object = Sub, features = genes, split.by = "orig.ident", cols = viridis(300))
plot_data <- plot$data %>% 
  select(pct.exp, id)
plot_data$status <- substr(plot_data$id, nchar(as.character(plot_data$id))-7, nchar(as.character(plot_data$id))-5)
plot_data$celltype <- factor(substr(plot_data$id, 1, nchar(as.character(plot_data$id))-9), levels = (levels(Sub)))
plot_data$orig.ident <- substr(plot_data$id, nchar(as.character(plot_data$id))-7, nchar(as.character(plot_data$id)))
start <- seq(from =1, to =length(genes)*length(levels(Sub)), by = length(genes))
end <- seq(from =length(genes), to =length(genes)*length(levels(Sub)), by = length(genes))

#Reorder into IEG x Cell type table by condition for PFC
DFSoc <- plot_data[plot_data$status == "Soc",]
DFSoc <- DFSoc[order(DFSoc$orig.ident, DFSoc$celltype),]
start_m <- seq(from = 1, to = length(unique(Sub$orig.ident[Sub$status == "Soc"]))*length(genes)*length(levels(Sub)), by = nrow(DFSoc)/length(unique(Sub$orig.ident[Sub$status == "Soc"])))
end_m <- seq(from = length(genes)*length(levels(Sub)), to = length(unique(Sub$orig.ident[Sub$status == "Soc"]))*length(genes)*length(levels(Sub)), by = nrow(DFSoc)/length(unique(Sub$orig.ident[Sub$status == "Soc"])))
Soc_list_table <- list()
for(m in 1:length(unique(Sub$orig.ident[Sub$status == "Soc"]))){
  rm(tmp, DF_table_Soc)
  tmp <- DFSoc$pct.exp[start_m[m]:end_m[m]]
  DF_table_Soc <- matrix(ncol = length(genes), nrow= length(levels(Sub)))
  for(i in 1:length(levels(Sub))){
    row <- tmp[start[i]:end[i]]
    DF_table_Soc[i,] <- row
  }
  colnames(DF_table_Soc) <- genes
  rownames(DF_table_Soc) <- levels(Sub)
  Soc_list_table[[m]] <- DF_table_Soc
}
names(Soc_list_table) <- names(unique(Sub$orig.ident[Sub$status == "Soc"]))

DFNot <- plot_data[plot_data$status == "Not",]
DFNot <- DFNot[order(DFNot$orig.ident, DFNot$celltype),]
start_m <- seq(from = 1, to = length(unique(Sub$orig.ident[Sub$status == "Not"]))*length(genes)*length(levels(Sub)), by = nrow(DFNot)/length(unique(Sub$orig.ident[Sub$status == "Not"])))
end_m <- seq(from = length(genes)*length(levels(Sub)), to = length(unique(Sub$orig.ident[Sub$status == "Not"]))*length(genes)*length(levels(Sub)), by = nrow(DFNot)/length(unique(Sub$orig.ident[Sub$status == "Not"])))
Not_list_table <- list()
for(m in 1:length(unique(Sub$orig.ident[Sub$status == "Not"]))){
  rm(tmp, DF_table_Not)
  tmp <- DFNot$pct.exp[start_m[m]:end_m[m]]
  DF_table_Not <- matrix(ncol = length(genes), nrow= length(levels(Sub)))
  for(i in 1:length(levels(Sub))){
    row <- tmp[start[i]:end[i]]
    DF_table_Not[i,] <- row
  }
  colnames(DF_table_Not) <- genes
  rownames(DF_table_Not) <- levels(Sub)
  Not_list_table[[m]] <- DF_table_Not
}
names(Not_list_table) <- names(unique(Sub$orig.ident[Sub$status == "Not"]))

DFShS <- plot_data[plot_data$status == "ShS",]
DFShS <- DFShS[order(DFShS$orig.ident, DFShS$celltype),]
start_m <- seq(from = 1, to = length(unique(Sub$orig.ident[Sub$status == "ShS"]))*length(genes)*length(levels(Sub)), by = nrow(DFShS)/length(unique(Sub$orig.ident[Sub$status == "ShS"])))
end_m <- seq(from = length(genes)*length(levels(Sub)), to = length(unique(Sub$orig.ident[Sub$status == "ShS"]))*length(genes)*length(levels(Sub)), by = nrow(DFShS)/length(unique(Sub$orig.ident[Sub$status == "ShS"])))
ShS_list_table <- list()
for(m in 1:length(unique(Sub$orig.ident[Sub$status == "ShS"]))){
  rm(tmp, DF_table_ShS)
  tmp <- DFShS$pct.exp[start_m[m]:end_m[m]]
  DF_table_ShS <- matrix(ncol = length(genes), nrow= length(levels(Sub)))
  for(i in 1:length(levels(Sub))){
    row <- tmp[start[i]:end[i]]
    DF_table_ShS[i,] <- row
  }
  colnames(DF_table_ShS) <- genes
  rownames(DF_table_ShS) <- levels(Sub)
  ShS_list_table[[m]] <- DF_table_ShS
}
names(ShS_list_table) <- names(unique(Sub$orig.ident[Sub$status == "ShS"]))

rm(i)
ShS.Not.list.pfc <- list()
for(i in 1:length(ShS_list_table)){
  ShS_mouse <- ShS_list_table[[i]]
  rm(j)
  int.list <- list()
  for(j in 1:length(Not_list_table)){
    rm(tab, Not_mouse)
    Not_mouse <- Not_list_table[[j]]
    tab <- (ShS_mouse[1:15,]+1)/(Not_mouse[1:15,]+1)
    int.list[[j]] <- tab
  }
  tmp <- (Reduce("+", int.list))/(length(int.list))
  ShS.Not.list.pfc[[i]] <- tmp
}
ShS.PFC.avg.mouse <- (Reduce("+", ShS.Not.list.pfc))/(length(ShS.Not.list.pfc))
ShS.PFC.avg.mouse <- ShS.PFC.avg.mouse[,order(colnames(ShS.PFC.avg.mouse))]

rm(i)
Soc.Not.list.pfc <- list()
for(i in 1:length(Soc_list_table)){
  Soc_mouse <- Soc_list_table[[i]]
  rm(j)
  int.list <- list()
  for(j in 1:length(Not_list_table)){
    rm(tab, Not_mouse)
    Not_mouse <- Not_list_table[[j]]
    tab <- (Soc_mouse[1:15,]+1)/(Not_mouse[1:15,]+1)
    int.list[[j]] <- tab
  }
  tmp <- (Reduce("+", int.list))/(length(int.list))
  Soc.Not.list.pfc[[i]] <- tmp
}
Soc.PFC.avg.mouse <- (Reduce("+", Soc.Not.list.pfc))/(length(Soc.Not.list.pfc))
Soc.PFC.avg.mouse <- Soc.PFC.avg.mouse[,order(colnames(Soc.PFC.avg.mouse))]

##Create heatmaps for Cereb
load("/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/Cereb-All.Rdata")

Sub <- subset(CellFile, subset= genotype == "WT")
Sub <- subset(Sub, idents = c("Astro", "Oligo", "OPC", "Microglia", "Endo Mural", "Endo Stalk", "Bergmann", "Choroid","Fibroblast", "Ependymal"), invert = TRUE)
#Full IEG panel
#genes <- (as.list(na.omit(read_csv("/home/hailee/10X/Wu, 2017 IEGs.csv", col_names = FALSE)[,1])))[[1]]
#2 fold or greater IEG panel
genes <- c("Fos", "Arc","Junb", "Egr4", "Btg2", "Dusp1", "Ccn1", "Fosb", "Egr3", "Ier5", "Egr2", "Fosl2", "Nr4a3", "Egr1", "Dusp6", "Per1", "Jun", "Nr4a2", "Ier2", "Gadd45g", "Cebpb", 
           "Bhlhe40", "Hes1", "Sik1", "Ccn2", "Ier3", "Gadd45b", "Nrn1", "Npas4", "F3", "Trib1", "Nr4a1", "Klf10", "Plk2", "Nfkbid", "Bdnf", "Zfp36l2", "Nfkbia", "Gem", "Ptgs2", "Klf2", 
           "Zfp36l1", "Apold1", "Dusp2", "Pim1", "Maff", "Ifit3", "Inhba", "Arhgef3", "Nfib", "Ifit2")#PFC + Cereb genes

#Pull out the proportion of the cluster that expresses each gene from the calculations made for making dotplots (Can be done many ways, but this was easiest)
plot <- DotPlot(object = Sub, features = genes, split.by = "orig.ident", cols = viridis(300))
plot_data <- plot$data %>% 
  select(pct.exp, id)
plot_data$status <- substr(plot_data$id, nchar(as.character(plot_data$id))-7, nchar(as.character(plot_data$id))-5)
plot_data$celltype <- factor(substr(plot_data$id, 1, nchar(as.character(plot_data$id))-9), levels = (levels(Sub)))
plot_data$orig.ident <- substr(plot_data$id, nchar(as.character(plot_data$id))-7, nchar(as.character(plot_data$id)))
start <- seq(from =1, to =length(genes)*length(levels(Sub)), by = length(genes))
end <- seq(from =length(genes), to =length(genes)*length(levels(Sub)), by = length(genes))

#Reorder into IEG x Cell type table by condition for Cerebellum
DFSoc <- plot_data[plot_data$status == "Soc",]
DFSoc <- DFSoc[order(DFSoc$orig.ident, DFSoc$celltype),]
start_m <- seq(from = 1, to = length(unique(Sub$orig.ident[Sub$status == "Soc"]))*length(genes)*length(levels(Sub)), by = nrow(DFSoc)/length(unique(Sub$orig.ident[Sub$status == "Soc"])))
end_m <- seq(from = length(genes)*length(levels(Sub)), to = length(unique(Sub$orig.ident[Sub$status == "Soc"]))*length(genes)*length(levels(Sub)), by = nrow(DFSoc)/length(unique(Sub$orig.ident[Sub$status == "Soc"])))
Soc_list_table <- list()
for(m in 1:length(unique(Sub$orig.ident[Sub$status == "Soc"]))){
  rm(tmp, DF_table_Soc)
  tmp <- DFSoc$pct.exp[start_m[m]:end_m[m]]
  DF_table_Soc <- matrix(ncol = length(genes), nrow= length(levels(Sub)))
  for(i in 1:length(levels(Sub))){
    row <- tmp[start[i]:end[i]]
    DF_table_Soc[i,] <- row
  }
  colnames(DF_table_Soc) <- genes
  rownames(DF_table_Soc) <- levels(Sub)
  Soc_list_table[[m]] <- DF_table_Soc
}
names(Soc_list_table) <- names(unique(Sub$orig.ident[Sub$status == "Soc"]))

DFNot <- plot_data[plot_data$status == "Not",]
DFNot <- DFNot[order(DFNot$orig.ident, DFNot$celltype),]
start_m <- seq(from = 1, to = length(unique(Sub$orig.ident[Sub$status == "Not"]))*length(genes)*length(levels(Sub)), by = nrow(DFNot)/length(unique(Sub$orig.ident[Sub$status == "Not"])))
end_m <- seq(from = length(genes)*length(levels(Sub)), to = length(unique(Sub$orig.ident[Sub$status == "Not"]))*length(genes)*length(levels(Sub)), by = nrow(DFNot)/length(unique(Sub$orig.ident[Sub$status == "Not"])))
Not_list_table <- list()
for(m in 1:length(unique(Sub$orig.ident[Sub$status == "Not"]))){
  rm(tmp, DF_table_Not)
  tmp <- DFNot$pct.exp[start_m[m]:end_m[m]]
  DF_table_Not <- matrix(ncol = length(genes), nrow= length(levels(Sub)))
  for(i in 1:length(levels(Sub))){
    row <- tmp[start[i]:end[i]]
    DF_table_Not[i,] <- row
  }
  colnames(DF_table_Not) <- genes
  rownames(DF_table_Not) <- levels(Sub)
  Not_list_table[[m]] <- DF_table_Not
}
names(Not_list_table) <- names(unique(Sub$orig.ident[Sub$status == "Not"]))

DFShS <- plot_data[plot_data$status == "ShS",]
DFShS <- DFShS[order(DFShS$orig.ident, DFShS$celltype),]
start_m <- seq(from = 1, to = length(unique(Sub$orig.ident[Sub$status == "ShS"]))*length(genes)*length(levels(Sub)), by = nrow(DFShS)/length(unique(Sub$orig.ident[Sub$status == "ShS"])))
end_m <- seq(from = length(genes)*length(levels(Sub)), to = length(unique(Sub$orig.ident[Sub$status == "ShS"]))*length(genes)*length(levels(Sub)), by = nrow(DFShS)/length(unique(Sub$orig.ident[Sub$status == "ShS"])))
ShS_list_table <- list()
for(m in 1:length(unique(Sub$orig.ident[Sub$status == "ShS"]))){
  rm(tmp, DF_table_ShS)
  tmp <- DFShS$pct.exp[start_m[m]:end_m[m]]
  DF_table_ShS <- matrix(ncol = length(genes), nrow= length(levels(Sub)))
  for(i in 1:length(levels(Sub))){
    row <- tmp[start[i]:end[i]]
    DF_table_ShS[i,] <- row
  }
  colnames(DF_table_ShS) <- genes
  rownames(DF_table_ShS) <- levels(Sub)
  ShS_list_table[[m]] <- DF_table_ShS
}
names(ShS_list_table) <- names(unique(Sub$orig.ident[Sub$status == "ShS"]))

rm(i)
ShS.Not.list.cereb <- list()
for(i in 1:length(ShS_list_table)){
  ShS_mouse <- ShS_list_table[[i]]
  rm(j)
  int.list <- list()
  for(j in 1:length(Not_list_table)){
    rm(tab, Not_mouse)
    Not_mouse <- Not_list_table[[j]]
    tab <- (ShS_mouse[1:7,]+1)/(Not_mouse[1:7,]+1)
    int.list[[j]] <- tab
  }
  tmp <- (Reduce("+", int.list))/(length(int.list)) #Average by mouse
  ShS.Not.list.cereb[[i]] <- tmp
}
ShS.Cereb.avg.mouse <- (Reduce("+", ShS.Not.list.cereb))/(length(ShS.Not.list.cereb)) #Average across mice
ShS.Cereb.avg.mouse <- ShS.Cereb.avg.mouse[,order(colnames(ShS.Cereb.avg.mouse))]

rm(i)
Soc.Not.list.cereb <- list()
for(i in 1:length(Soc_list_table)){
  Soc_mouse <- Soc_list_table[[i]]
  rm(j)
  int.list <- list()
  for(j in 1:length(Not_list_table)){
    rm(tab, Not_mouse)
    Not_mouse <- Not_list_table[[j]]
    tab <- (Soc_mouse[1:7,]+1)/(Not_mouse[1:7,]+1)
    int.list[[j]] <- tab
  }
  tmp <- (Reduce("+", int.list))/(length(int.list))
  Soc.Not.list.cereb[[i]] <- tmp
}
Soc.Cereb.avg.mouse <- (Reduce("+", Soc.Not.list.cereb))/(length(Soc.Not.list.cereb))
Soc.Cereb.avg.mouse <- Soc.Cereb.avg.mouse[,order(colnames(Soc.Cereb.avg.mouse))]

##Plot heatmaps of real data
#Replace Inf with a high value as it is caused by #/0 and NaN with 1 as it's caused by 0/0 (no change)
#Columns ordered by 10min pfc
col_fun = colorRamp2(c(0, 1.1, 2.5), c("blue", "white", "red"))
col.order.pfc <-  c("Fos", "Arc","Junb", "Egr4", "Btg2", "Dusp1", "Ccn1", "Fosb", "Egr3", "Ier5", "Egr2", "Fosl2", "Nr4a3", "Egr1", "Dusp6", "Per1", "Jun", "Nr4a2", "Ier2", "Gadd45g", "Cebpb", 
                    "Bhlhe40", "Hes1", "Sik1", "Ccn2", "Ier3", "Gadd45b", "Nrn1", "Npas4", "F3", "Trib1", "Nr4a1", "Klf10", "Plk2", "Nfkbid", "Bdnf", "Zfp36l2", "Nfkbia", "Gem", "Ptgs2", "Klf2", 
                    "Zfp36l1", "Apold1", "Dusp2", "Pim1", "Maff", "Ifit3", "Inhba", "Arhgef3", "Nfib", "Ifit2")
ShS.PFC.avg.mouse <- ShS.PFC.avg.mouse[,col.order.pfc]
Soc.PFC.avg.mouse <- Soc.PFC.avg.mouse[,col.order.pfc]
Heatmap((ShS.PFC.avg.mouse), cluster_rows = F,name = "ShS/Not prop", cluster_columns = F, 
        column_title = "10 min+1/ Baseline+1 - PFC", col = col_fun)
Heatmap((Soc.PFC.avg.mouse), cluster_rows = F, name = "Soc/Not prop",cluster_columns = F, 
        column_title = "35 min+1/Baseline+1 - PFC", col = col_fun)
ShS.Cereb.avg.mouse <- ShS.Cereb.avg.mouse[,col.order.pfc]
Soc.Cereb.avg.mouse <- Soc.Cereb.avg.mouse[,col.order.pfc]
Heatmap((ShS.Cereb.avg.mouse), cluster_rows = F,  
        name = "ShS/Not prop", cluster_columns = F, column_title = "10 min+1/Baseline+1 - PFC", col = col_fun)
Heatmap((Soc.Cereb.avg.mouse), cluster_rows = F, 
        name = "Soc/Not prop",cluster_columns = F, column_title = "35 min+1/Baseline+1 - PFC", col = col_fun)
#Columns ordered by 10min cereb
col.order.cereb <-  c("Nr4a3", "Fos", "Nr4a1", "Nr4a2", "Per1", "Sik1", "Dusp1", "Fosl2", "Ier2", "Inhba", "Nfkbid", "Pim1", "Fosb", "Ccn1", "Dusp6", "Arc", "Npas4", "Cebpb", "Plk2", 
                      "Gem", "Gadd45g", "Junb", "Klf2", "Zfp36l2", "Trib1", "Ccn2", "Bdnf", "Gadd45b", "Klf10", "Hes1", "Ifit2", "Btg2", "Jun", "Nrn1", "Zfp36l1", "F3", "Ier5", "Nfib", "Maff", 
                      "Egr2", "Dusp2", "Ptgs2", "Ifit3", "Nfkbia", "Bhlhe40", "Apold1", "Egr3", "Egr4", "Ier3", "Egr1", "Arhgef3")
ShS.PFC.avg.mouse <- ShS.PFC.avg.mouse[,col.order.cereb]
Soc.PFC.avg.mouse <- Soc.PFC.avg.mouse[,col.order.cereb]
Heatmap((ShS.PFC.avg.mouse), cluster_rows = F,name = "ShS/Not prop", cluster_columns = F, 
        column_title = "10 min/ Baseline - Cereb", col = col_fun)
Heatmap((Soc.PFC.avg.mouse), cluster_rows = F, name = "Soc/Not prop",cluster_columns = F, 
        column_title = "35 min/Baseline - Cereb", col = col_fun)
ShS.Cereb.avg.mouse <- ShS.Cereb.avg.mouse[,col.order.cereb]
Soc.Cereb.avg.mouse <- Soc.Cereb.avg.mouse[,col.order.cereb]
Heatmap((ShS.Cereb.avg.mouse), cluster_rows = F,  
        name = "ShS/Not prop", cluster_columns = F, column_title = "10 min+1/Baseline+1 - Cereb", col = col_fun)
Heatmap((Soc.Cereb.avg.mouse), cluster_rows = F, 
        name = "Soc/Not prop",cluster_columns = F, column_title = "35 min+1/Baseline+1 - Cereb", col = col_fun)


##Figure 2D ##_________________________________________________________________________________________
#Calculate correlation matrix between each mouse
#rerun all of above but add one to numerator and denominator before doing the division to make ShS.Not.lists and Soc.Not.lists
corr.list.10 <- list()
corr.list.10.na <- list()
rm(i)
for(i in 1:length(ShS.Not.list.pfc)){
  tmp1 <- ShS.Not.list.pfc[[i]]
  tmp2 <- ShS.Not.list.cereb[[i]]
  rm(j)
  corr.list.by.mouse <- list()
  for(j in 1:length(Soc.Not.list.pfc)){
    tmp3 <- Soc.Not.list.pfc[[j]]
    tmp4 <- Soc.Not.list.cereb[[j]]
    comb <- rbind(tmp1, tmp2, tmp3, tmp4)
    corr<- cor(t(comb))
    corr.list.by.mouse[[j]] <- corr
  }
  corr.list.10.na[[i]] <- corr.list.by.mouse
  tmp <- (Reduce("+", corr.list.by.mouse))/(length(corr.list.by.mouse))
  corr.list.10[[i]] <- tmp
}

corr.list.35 <- list()
corr.list.35.na <- list()
rm(i)
for(i in 1:length(Soc.Not.list.pfc)){
  tmp1 <- Soc.Not.list.pfc[[i]]
  tmp2 <- Soc.Not.list.cereb[[i]]
  rm(j)
  corr.list.by.mouse <- list()
  for(j in 1:length(ShS.Not.list.pfc)){
    tmp3 <- ShS.Not.list.pfc[[j]]
    tmp4 <- ShS.Not.list.cereb[[j]]
    comb <- rbind(tmp3, tmp4, tmp1, tmp2)
    corr<- cor(t(comb))
    corr.list.by.mouse[[j]] <- corr
  }
  corr.list.35.na[[i]] <- corr.list.by.mouse
  tmp <- (Reduce("+", corr.list.by.mouse))/(length(corr.list.by.mouse))
  corr.list.35[[i]] <- tmp
}
corr.avg <- (Reduce("+", corr.list.35))/(length(corr.list.35))
col_fun = colorRamp2(c(-0.25, 0.1, 0.55,1), c("blue", "white","red", "red"))
pheatmap(corr.avg, cluster_rows = F, cluster_cols = F, col = col_fun) ## Plot in Fig.3C

##Figure 3 E ##____________________________________________________________________________________________________________________________________
##calculate correlation with 10 min PFC and other conditions
pfc.10.corr.exc <- list()
pfc.10.corr.inh <- list()
rm(i)
for(i in 1:length(ShS.Not.list.pfc)){
  pfc.early <- ShS.Not.list.pfc[[i]]
  cereb10.pfc10.exc <- matrix(ncol=1, nrow = length(ShS.Not.list.cereb))
  cereb10.pfc10.inh <- matrix(ncol=1, nrow = length(ShS.Not.list.cereb))
  pfc10.pfc10.exc <- matrix(ncol=1, nrow = length(ShS.Not.list.pfc))
  pfc10.pfc10.inh <- matrix(ncol=1, nrow = length(ShS.Not.list.pfc))
  rm(j)
  for(j in 1:length(ShS.Not.list.cereb)){
    cereb.early <- ShS.Not.list.cereb[[j]]
    pfc.early.tmp <- ShS.Not.list.pfc[[j]]
    comb <- rbind(pfc.early, cereb.early, pfc.early.tmp)
    corr<- cor(t(comb))
    cereb10.pfc10.exc[j,] <- mean(corr[1:8,16:17])
    cereb10.pfc10.inh[j,] <- mean(corr[9:15,18:22])
    pfc10.pfc10.exc[j,] <- mean(corr[1:8,23:30])
    pfc10.pfc10.inh[j,] <- mean(corr[9:15,31:37])
  }
  pfc10.pfc35.exc <-matrix(ncol=1, nrow = length(Soc.Not.list.pfc))
  pfc10.pfc35.inh <-matrix(ncol=1, nrow = length(Soc.Not.list.pfc))
  rm(c)
  for(c in 1:length(Soc.Not.list.pfc)){
    pfc.late <- Soc.Not.list.pfc[[c]]
    comb <- rbind(pfc.early, pfc.late)
    corr<- cor(t(comb))
    pfc10.pfc35.exc[c,] <- mean(corr[1:8,16:23])
    pfc10.pfc35.inh[c,] <- mean(corr[9:15,24:30])
  }
  pfc10.cereb35.exc <-matrix(ncol=1, nrow = length(Soc.Not.list.cereb))
  pfc10.cereb35.inh <-matrix(ncol=1, nrow = length(Soc.Not.list.cereb))
  rm(n)
  for(n in 1:length(Soc.Not.list.cereb)){
    cereb.late <- Soc.Not.list.cereb[[n]]
    comb <- rbind(pfc.early, cereb.late)
    corr<- cor(t(comb))
    pfc10.cereb35.exc[n,] <- mean(corr[1:8,16:17])
    pfc10.cereb35.inh[n,] <- mean(corr[9:15,18:22])
  }
  pfc.10.corr.exc[[i]] <- list(pfc10.pfc10.exc, cereb10.pfc10.exc, pfc10.pfc35.exc, pfc10.cereb35.exc)
  pfc.10.corr.inh[[i]] <- list(pfc10.pfc10.inh, cereb10.pfc10.inh, pfc10.pfc35.inh, pfc10.cereb35.inh)
}
#Get all of the excitatory values
pfc10.pfc10.vals.exc <- rbind(pfc.10.corr.exc[[1]][[1]], pfc.10.corr.exc[[2]][[1]],pfc.10.corr.exc[[3]][[1]])
View(pfc10.pfc10.vals.exc)
pfc10.cereb10.vals.exc <- rbind(pfc.10.corr.exc[[1]][[2]], pfc.10.corr.exc[[2]][[2]],pfc.10.corr.exc[[3]][[2]])
View(pfc10.cereb10.vals.exc)
pfc10.pfc35.vals.exc <- rbind(pfc.10.corr.exc[[1]][[3]], pfc.10.corr.exc[[2]][[3]],pfc.10.corr.exc[[3]][[3]])
View(pfc10.pfc35.vals.exc)
pfc10.cereb35.vals.exc <- rbind(pfc.10.corr.exc[[1]][[4]], pfc.10.corr.exc[[2]][[4]],pfc.10.corr.exc[[3]][[4]])
View(pfc10.cereb35.vals.exc)
##I then put the values into a CSV and used prism to plot because the differing lengths make it difficult to graph here

#Get all of the inhibitory values
pfc10.pfc10.vals.inh <- rbind(pfc.10.corr.inh[[1]][[1]], pfc.10.corr.inh[[2]][[1]],pfc.10.corr.inh[[3]][[1]])
View(pfc10.pfc10.vals.inh)
pfc10.cereb10.vals.inh <- rbind(pfc.10.corr.inh[[1]][[2]], pfc.10.corr.inh[[2]][[2]],pfc.10.corr.inh[[3]][[2]])
View(pfc10.cereb10.vals.inh)
pfc10.pfc35.vals.inh <- rbind(pfc.10.corr.inh[[1]][[3]], pfc.10.corr.inh[[2]][[3]],pfc.10.corr.inh[[3]][[3]])
View(pfc10.pfc35.vals.inh)
pfc10.cereb35.vals.inh <- rbind(pfc.10.corr.inh[[1]][[4]], pfc.10.corr.inh[[2]][[4]],pfc.10.corr.inh[[3]][[4]])
View(pfc10.cereb35.vals.inh)
##I then put the values into a CSV and used prism to plot because the differing lengths make it difficult to graph here

##calculate correlation with 10 min Cereb and other conditions
cereb.10.corr.exc <- list()
cereb.10.corr.inh <- list()
rm(i)
for(i in 1:length(ShS.Not.list.cereb)){
  cereb.early <- ShS.Not.list.cereb[[i]]
  cereb10.pfc10.exc <- matrix(ncol=1, nrow = length(ShS.Not.list.cereb))
  cereb10.pfc10.inh <- matrix(ncol=1, nrow = length(ShS.Not.list.cereb))
  cereb10.cereb10.exc <- matrix(ncol=1, nrow = length(ShS.Not.list.cereb))
  cereb10.cereb10.inh <- matrix(ncol=1, nrow = length(ShS.Not.list.cereb))
  rm(j)
  for(j in 1:length(ShS.Not.list.pfc)){
    pfc.early <- ShS.Not.list.pfc[[j]]
    cereb.early.tmp <- ShS.Not.list.cereb[[j]]
    comb <- rbind(cereb.early, pfc.early, cereb.early.tmp)
    corr<- cor(t(comb))
    cereb10.pfc10.exc[j,] <- mean(corr[1:2,8:15])
    cereb10.pfc10.inh[j,] <- mean(corr[3:7,16:22])
    cereb10.cereb10.exc[j,] <- mean(corr[1:2,23:24])
    cereb10.cereb10.inh[j,] <- mean(corr[3:7,25:29])
  }
  cereb10.pfc35.exc <-matrix(ncol=1, nrow = length(Soc.Not.list.pfc))
  cereb10.pfc35.inh <-matrix(ncol=1, nrow = length(Soc.Not.list.pfc))
  rm(c)
  for(c in 1:length(Soc.Not.list.pfc)){
    pfc.late <- Soc.Not.list.pfc[[c]]
    comb <- rbind(cereb.early, pfc.late)
    corr<- cor(t(comb))
    cereb10.pfc35.exc[c,] <- mean(corr[1:2,8:15])
    cereb10.pfc35.inh[c,] <- mean(corr[3:7,16:22])
  }
  cereb10.cereb35.exc <-matrix(ncol=1, nrow = length(Soc.Not.list.cereb))
  cereb10.cereb35.inh <-matrix(ncol=1, nrow = length(Soc.Not.list.cereb))
  rm(n)
  for(n in 1:length(Soc.Not.list.cereb)){
    cereb.late <- Soc.Not.list.cereb[[n]]
    comb <- rbind(cereb.early, cereb.late)
    corr<- cor(t(comb))
    cereb10.cereb35.exc[n,] <- mean(corr[1:2,8:9])
    cereb10.cereb35.inh[n,] <- mean(corr[3:7,10:nrow(corr)])
  }
  cereb.10.corr.exc[[i]] <- list(cereb10.cereb10.exc, cereb10.pfc10.exc, cereb10.pfc35.exc, cereb10.cereb35.exc)
  cereb.10.corr.inh[[i]] <- list(cereb10.cereb10.inh, cereb10.pfc10.inh, cereb10.pfc35.inh, cereb10.cereb35.inh)
}
cereb10.cereb10.vals.exc <- rbind(cereb.10.corr.exc[[1]][[1]], cereb.10.corr.exc[[2]][[1]],cereb.10.corr.exc[[3]][[1]])
View(cereb10.cereb10.vals.exc)
cereb10.pfc10.vals.exc <- rbind(cereb.10.corr.exc[[1]][[2]], cereb.10.corr.exc[[2]][[2]],cereb.10.corr.exc[[3]][[2]])
View(cereb10.pfc10.vals.exc)
cereb10.pfc35.vals.exc <- rbind(cereb.10.corr.exc[[1]][[3]], cereb.10.corr.exc[[2]][[3]],cereb.10.corr.exc[[3]][[3]])
View(cereb10.pfc35.vals.exc)
cereb10.cereb35.vals.exc <- rbind(cereb.10.corr.exc[[1]][[4]], cereb.10.corr.exc[[2]][[4]],cereb.10.corr.exc[[3]][[4]])
View(cereb10.cereb35.vals.exc)
##I then put the values into a CSV and used prism to plot because the differing lengths make it difficult to graph here

cereb10.cereb10.vals.inh <- rbind(cereb.10.corr.inh[[1]][[1]], cereb.10.corr.inh[[2]][[1]],cereb.10.corr.inh[[3]][[1]])
View(cereb10.cereb10.vals.inh)
cereb10.pfc10.vals.inh <- rbind(cereb.10.corr.inh[[1]][[2]], cereb.10.corr.inh[[2]][[2]],cereb.10.corr.inh[[3]][[2]])
View(cereb10.pfc10.vals.inh)
cereb10.pfc35.vals.inh <- rbind(cereb.10.corr.inh[[1]][[3]], cereb.10.corr.inh[[2]][[3]],cereb.10.corr.inh[[3]][[3]])
View(cereb10.pfc35.vals.inh)
cereb10.cereb35.vals.inh <- rbind(cereb.10.corr.inh[[1]][[4]], cereb.10.corr.inh[[2]][[4]],cereb.10.corr.inh[[3]][[4]])
View(cereb10.cereb35.vals.inh)
##I then put the values into a CSV and used prism to plot because the differing lengths make it difficult to graph here


##calculate correlation with 35 min Cereb and other conditions
cereb.35.corr.exc <- list()
cereb.35.corr.inh <- list()
rm(i)
for(i in 1:length(Soc.Not.list.cereb)){
  cereb.late <- Soc.Not.list.cereb[[i]]
  cereb35.pfc10.exc <- matrix(ncol=1, nrow = length(ShS.Not.list.cereb))
  cereb35.pfc10.inh <- matrix(ncol=1, nrow = length(ShS.Not.list.cereb))
  cereb35.cereb10.exc <- matrix(ncol=1, nrow = length(ShS.Not.list.cereb))
  cereb35.cereb10.inh <- matrix(ncol=1, nrow = length(ShS.Not.list.cereb))
  rm(j)
  for(j in 1:length(ShS.Not.list.pfc)){
    pfc.early <- ShS.Not.list.pfc[[j]]
    cereb.early <- ShS.Not.list.cereb[[j]]
    comb <- rbind(cereb.late, pfc.early, cereb.early)
    corr<- cor(t(comb))
    cereb35.pfc10.exc[j,] <- mean(corr[1:2,8:15])
    cereb35.pfc10.inh[j,] <- mean(corr[3:7,16:22])
    cereb35.cereb10.exc[j,] <- mean(corr[1:2,23:24])
    cereb35.cereb10.inh[j,] <- mean(corr[3:7,25:29])
  }
  cereb35.pfc35.exc <-matrix(ncol=1, nrow = length(Soc.Not.list.pfc))
  cereb35.pfc35.inh <-matrix(ncol=1, nrow = length(Soc.Not.list.pfc))
  rm(c)
  for(c in 1:length(Soc.Not.list.pfc)){
    pfc.late <- Soc.Not.list.pfc[[c]]
    comb <- rbind(cereb.late, pfc.late)
    corr<- cor(t(comb))
    cereb35.pfc35.exc[c,] <- mean(corr[1:2,8:15])
    cereb35.pfc35.inh[c,] <- mean(corr[3:7,16:22])
  }
  cereb35.cereb35.exc <-matrix(ncol=1, nrow = length(Soc.Not.list.cereb))
  cereb35.cereb35.inh <-matrix(ncol=1, nrow = length(Soc.Not.list.cereb))
  rm(n)
  for(n in 1:length(Soc.Not.list.cereb)){
    cereb.late.tmp <- Soc.Not.list.cereb[[n]]
    comb <- rbind(cereb.late, cereb.late.tmp)
    corr<- cor(t(comb))
    cereb35.cereb35.exc[n,] <- mean(corr[1:2,8:9])
    cereb35.cereb35.inh[n,] <- mean(corr[3:7,10:nrow(corr)])
  }
  cereb.35.corr.exc[[i]] <- list(cereb35.cereb10.exc, cereb35.pfc10.exc, cereb35.pfc35.exc, cereb35.cereb35.exc)
  cereb.35.corr.inh[[i]] <- list(cereb35.cereb10.inh, cereb35.pfc10.inh, cereb35.pfc35.inh, cereb35.cereb35.inh)
}
cereb35.cereb10.vals.exc <- rbind(cereb.35.corr.exc[[1]][[1]], cereb.35.corr.exc[[2]][[1]],cereb.35.corr.exc[[3]][[1]], cereb.35.corr.exc[[4]][[1]], cereb.35.corr.exc[[5]][[1]],cereb.35.corr.exc[[6]][[1]])
View(cereb35.cereb10.vals.exc)
cereb35.pfc10.vals.exc <- rbind(cereb.35.corr.exc[[1]][[2]], cereb.35.corr.exc[[2]][[2]],cereb.35.corr.exc[[3]][[2]], cereb.35.corr.exc[[4]][[2]], cereb.35.corr.exc[[5]][[2]],cereb.35.corr.exc[[6]][[2]])
View(cereb35.pfc10.vals.exc)
cereb35.pfc35.vals.exc <- rbind(cereb.35.corr.exc[[1]][[3]], cereb.35.corr.exc[[2]][[3]],cereb.35.corr.exc[[3]][[3]], cereb.35.corr.exc[[4]][[3]], cereb.35.corr.exc[[5]][[3]],cereb.35.corr.exc[[6]][[3]])
View(cereb35.pfc35.vals.exc)
cereb35.cereb35.vals.exc <- rbind(cereb.35.corr.exc[[1]][[4]], cereb.35.corr.exc[[2]][[4]],cereb.35.corr.exc[[3]][[4]], cereb.35.corr.exc[[4]][[4]], cereb.35.corr.exc[[5]][[4]],cereb.35.corr.exc[[6]][[4]])
View(cereb35.cereb35.vals.exc)
##I then put the values into a CSV and used prism to plot because the differing lengths make it difficult to graph here

cereb35.cereb10.vals.inh <- rbind(cereb.35.corr.inh[[1]][[1]], cereb.35.corr.inh[[2]][[1]],cereb.35.corr.inh[[3]][[1]], cereb.35.corr.exc[[4]][[1]], cereb.35.corr.exc[[5]][[1]],cereb.35.corr.exc[[6]][[1]])
View(cereb35.cereb10.vals.inh)
cereb35.pfc10.vals.inh <- rbind(cereb.35.corr.inh[[1]][[2]], cereb.35.corr.inh[[2]][[2]],cereb.35.corr.inh[[3]][[2]], cereb.35.corr.exc[[4]][[2]], cereb.35.corr.exc[[5]][[2]],cereb.35.corr.exc[[6]][[2]])
View(cereb35.pfc10.vals.inh)
cereb35.pfc35.vals.inh <- rbind(cereb.35.corr.inh[[1]][[3]], cereb.35.corr.inh[[2]][[3]],cereb.35.corr.inh[[3]][[3]], cereb.35.corr.exc[[4]][[3]], cereb.35.corr.exc[[5]][[3]],cereb.35.corr.exc[[6]][[3]])
View(cereb35.pfc35.vals.inh)
cereb35.cereb35.vals.inh <- rbind(cereb.35.corr.inh[[1]][[4]], cereb.35.corr.inh[[2]][[4]],cereb.35.corr.inh[[3]][[4]], cereb.35.corr.exc[[4]][[4]], cereb.35.corr.exc[[5]][[4]],cereb.35.corr.exc[[6]][[4]])
View(cereb35.cereb35.vals.inh)
##I then put the values into a CSV and used prism to plot because the differing lengths make it difficult to graph here

##calculate correlation with 35 min PFC and other conditions
pfc.35.corr.exc <- list()
pfc.35.corr.inh <- list()
rm(i)
for(i in 1:length(Soc.Not.list.pfc)){
  pfc.late <- Soc.Not.list.pfc[[i]]
  pfc35.pfc10.exc <- matrix(ncol=1, nrow = length(ShS.Not.list.pfc))
  pfc35.pfc10.inh <- matrix(ncol=1, nrow = length(ShS.Not.list.pfc))
  pfc35.cereb10.exc <- matrix(ncol=1, nrow = length(ShS.Not.list.cereb))
  pfc35.cereb10.inh <- matrix(ncol=1, nrow = length(ShS.Not.list.cereb))
  rm(j)
  for(j in 1:length(ShS.Not.list.pfc)){
    pfc.early <- ShS.Not.list.pfc[[j]]
    cereb.early <- ShS.Not.list.cereb[[j]]
    comb <- rbind(pfc.late, pfc.early, cereb.early)
    corr<- cor(t(comb))
    pfc35.pfc10.exc[j,] <- mean(corr[1:8,16:23])
    pfc35.pfc10.inh[j,] <- mean(corr[9:15,24:30])
    pfc35.cereb10.exc[j,] <- mean(corr[1:8,31:32])
    pfc35.cereb10.inh[j,] <- mean(corr[9:15,33:37])
  }
  pfc35.cereb35.exc <-matrix(ncol=1, nrow = length(Soc.Not.list.cereb))
  pfc35.cereb35.inh <-matrix(ncol=1, nrow = length(Soc.Not.list.cereb))
  rm(c)
  for(c in 1:length(Soc.Not.list.cereb)){
    cereb.late <- Soc.Not.list.cereb[[c]]
    comb <- rbind(cereb.late, pfc.late)
    corr<- cor(t(comb))
    pfc35.cereb35.exc[c,] <- mean(corr[1:2,8:15])
    pfc35.cereb35.inh[c,] <- mean(corr[3:7,16:22])
  }
  pfc35.pfc35.exc <-matrix(ncol=1, nrow = length(Soc.Not.list.pfc))
  pfc35.pfc35.inh <-matrix(ncol=1, nrow = length(Soc.Not.list.pfc))
  rm(n)
  for(n in 1:length(Soc.Not.list.pfc)){
    pfc.late.tmp <- Soc.Not.list.pfc[[n]]
    comb <- rbind(pfc.late, pfc.late.tmp)
    corr<- cor(t(comb))
    pfc35.pfc35.exc[n,] <- mean(corr[1:8,16:23])
    pfc35.pfc35.inh[n,] <- mean(corr[9:15,24:30])
  }
  pfc.35.corr.exc[[i]] <- list(pfc35.pfc35.exc, pfc35.pfc10.exc, pfc35.cereb35.exc, pfc35.cereb10.exc)
  pfc.35.corr.inh[[i]] <- list(pfc35.pfc35.inh, pfc35.pfc10.inh, pfc35.cereb35.inh, pfc35.cereb10.inh)
}
pfc35.pfc35.vals.exc <- rbind(pfc.35.corr.exc[[1]][[1]], pfc.35.corr.exc[[2]][[1]],pfc.35.corr.exc[[3]][[1]], pfc.35.corr.exc[[4]][[1]], pfc.35.corr.exc[[5]][[1]])
View(pfc35.pfc35.vals.exc)
pfc35.pfc10.vals.exc <- rbind(pfc.35.corr.exc[[1]][[2]], pfc.35.corr.exc[[2]][[2]],pfc.35.corr.exc[[3]][[2]], pfc.35.corr.exc[[4]][[2]], pfc.35.corr.exc[[5]][[2]])
View(pfc35.pfc10.vals.exc)
pfc35.cereb35.vals.exc <- rbind(pfc.35.corr.exc[[1]][[3]], pfc.35.corr.exc[[2]][[3]],pfc.35.corr.exc[[3]][[3]], pfc.35.corr.exc[[4]][[3]], pfc.35.corr.exc[[5]][[3]])
View(pfc35.cereb35.vals.exc)
pfc35.cereb10.vals.exc <- rbind(pfc.35.corr.exc[[1]][[4]], pfc.35.corr.exc[[2]][[4]],pfc.35.corr.exc[[3]][[4]], pfc.35.corr.exc[[4]][[4]], pfc.35.corr.exc[[5]][[4]])
View(pfc35.cereb10.vals.exc)
##I then put the values into a CSV and used prism to plot because the differing lengths make it difficult to graph here

pfc35.pfc35.vals.inh <- rbind(pfc.35.corr.inh[[1]][[1]], pfc.35.corr.inh[[2]][[1]],pfc.35.corr.inh[[3]][[1]], pfc.35.corr.exc[[4]][[1]], pfc.35.corr.exc[[5]][[1]])
View(pfc35.pfc35.vals.inh)
pfc35.pfc10.vals.inh <- rbind(pfc.35.corr.inh[[1]][[2]], pfc.35.corr.inh[[2]][[2]],pfc.35.corr.inh[[3]][[2]], pfc.35.corr.exc[[4]][[2]], pfc.35.corr.exc[[5]][[2]])
View(pfc35.pfc10.vals.inh)
pfc35.cereb35.vals.inh <- rbind(pfc.35.corr.inh[[1]][[3]], pfc.35.corr.inh[[2]][[3]],pfc.35.corr.inh[[3]][[3]], pfc.35.corr.exc[[4]][[3]], pfc.35.corr.exc[[5]][[3]])
View(pfc35.cereb35.vals.inh)
pfc35.cereb10.vals.inh <- rbind(pfc.35.corr.inh[[1]][[4]], pfc.35.corr.inh[[2]][[4]],pfc.35.corr.inh[[3]][[4]], pfc.35.corr.exc[[4]][[4]], pfc.35.corr.exc[[5]][[4]])
View(pfc35.cereb10.vals.inh)
##I then put the values into a CSV and used prism to plot because the differing lengths make it difficult to graph here
