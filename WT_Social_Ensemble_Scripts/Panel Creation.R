###Figure3 - Heatmaps
## Hailee Walker - 07142023
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

## Figure 3A and 3B PFC ## _________________________________________________________________________________________________________
##Rerun but split by mouse
##Create heatmaps for PFC
load("~/10X/Final Soc Analysis/PFC_All_analysis/Comb_PFC.Rdata")

Sub <- subset(CellFile, subset= genotype == "WT")
Sub <- subset(Sub, idents = c("Astro", "Oligo", "OPC","Microglia","Endothelial", "VLMC","Pericyte"), invert = TRUE)
#Full IEG panel
genes <- (as.list(na.omit(read_csv("/home/hailee/10X/IEG list (Hong,W. 2017).csv", col_names = FALSE)[,1])))[[1]]

#Pull out the proportion of the cluster that expresses each gene from the calculations made for making dotplots (Can be done many ways, but this was easiest)
plot <- DotPlot(object = Sub, features = genes, split.by = "status", cols = viridis(300))
plot_data <- plot$data %>% 
  select(pct.exp, id)
plot_data$status <- substr(plot_data$id, nchar(as.character(plot_data$id))-2, nchar(as.character(plot_data$id)))
plot_data$celltype <- factor(substr(plot_data$id, 1, nchar(as.character(plot_data$id))-4), levels = levels(Sub))
start <- seq(from =1, to =length(genes)*length(levels(Sub)), by = length(genes))
end <- seq(from =length(genes), to =length(genes)*length(levels(Sub)), by = length(genes))

#Reorder into IEG x Cell type table by condition for PFC
DFSoc <- plot_data[plot_data$status == "Soc",]
DFSoc <- DFSoc[order(DFSoc$celltype),]
DF_table_Soc <- matrix(ncol = length(genes), nrow= length(levels(Sub)))
  for(i in 1:length(levels(Sub))){
    row <- DFSoc$pct.exp[start[i]:end[i]]
    DF_table_Soc[i,] <- row
  }
colnames(DF_table_Soc) <- genes
rownames(DF_table_Soc) <- levels(Sub)

DFNot <- plot_data[plot_data$status == "Not",]
DFNot <- DFNot[order(DFNot$celltype),]
DF_table_Not <- matrix(ncol = length(genes), nrow= length(levels(Sub)))
  for(i in 1:length(levels(Sub))){
    row <- DFNot$pct.exp[start[i]:end[i]]
    DF_table_Not[i,] <- row
  }
colnames(DF_table_Not) <- genes
rownames(DF_table_Not) <- levels(Sub)

DFShS <- plot_data[plot_data$status == "ShS",]
DFShS <- DFShS[order(DFShS$celltype),]
DF_table_ShS <- matrix(ncol = length(genes), nrow= length(levels(Sub)))
  for(i in 1:length(levels(Sub))){
    row <- DFShS$pct.exp[start[i]:end[i]]
    DF_table_ShS[i,] <- row
  }
colnames(DF_table_ShS) <- genes
rownames(DF_table_ShS) <- levels(Sub)

rm(i)

ShSPFC <- (DF_table_ShS)/(DF_table_Not+1)
ShSPFC[is.infinite(ShSPFC)] <-NaN 
shsmt <- matrix(nrow=1, ncol = ncol(ShSPFC))
for(i in 1:ncol(ShSPFC)){
  mx <- max(na.omit(ShSPFC[,i]))
  shsmt[i] <- mx
}
colnames(shsmt) = colnames(ShSPFC)
colnames(shsmt)[shsmt >=2]

SocPFC <- (DF_table_Soc)/(DF_table_Not+1)
SocPFC[is.infinite(SocPFC)] <-NaN 
Socmt <- matrix(nrow=1, ncol = ncol(SocPFC))
for(i in 1:ncol(SocPFC)){
  mx <- max(na.omit(SocPFC[,i]))
  Socmt[i] <- mx
}
colnames(Socmt) = colnames(SocPFC)
colnames(Socmt)[Socmt >=2]

##Create heatmaps for Cereb
load("/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/Cereb-All.Rdata")

Sub <- subset(CellFile, subset= genotype == "WT")
Sub <- subset(Sub, idents = c("Astro", "Oligo", "OPC", "Microglia", "Endo Mural", "Endo Stalk", "Bergmann", "Choroid","Fibroblast", "Ependymal"), invert = TRUE)
#Full IEG panel
genes <- (as.list(na.omit(read_csv("/home/hailee/10X/Wu, 2017 IEGs.csv", col_names = FALSE)[,1])))[[1]]
#2 fold or greater IEG panel

#Pull out the proportion of the cluster that expresses each gene from the calculations made for making dotplots (Can be done many ways, but this was easiest)
plot <- DotPlot(object = Sub, features = genes, split.by = "status", cols = viridis(300))
plot_data <- plot$data %>% 
  select(pct.exp, id)
plot_data$status <- substr(plot_data$id, nchar(as.character(plot_data$id))-2, nchar(as.character(plot_data$id)))
plot_data$celltype <- factor(substr(plot_data$id, 1, nchar(as.character(plot_data$id))-4), levels = levels(Sub))
start <- seq(from =1, to =length(genes)*length(levels(Sub)), by = length(genes))
end <- seq(from =length(genes), to =length(genes)*length(levels(Sub)), by = length(genes))

#Reorder into IEG x Cell type table by condition for Cereb
DFSoc <- plot_data[plot_data$status == "Soc",]
DFSoc <- DFSoc[order(DFSoc$celltype),]
DF_table_Soc <- matrix(ncol = length(genes), nrow= length(levels(Sub)))
for(i in 1:length(levels(Sub))){
  row <- DFSoc$pct.exp[start[i]:end[i]]
  DF_table_Soc[i,] <- row
}
colnames(DF_table_Soc) <- genes
rownames(DF_table_Soc) <- levels(Sub)

DFNot <- plot_data[plot_data$status == "Not",]
DFNot <- DFNot[order(DFNot$celltype),]
DF_table_Not <- matrix(ncol = length(genes), nrow= length(levels(Sub)))
for(i in 1:length(levels(Sub))){
  row <- DFNot$pct.exp[start[i]:end[i]]
  DF_table_Not[i,] <- row
}
colnames(DF_table_Not) <- genes
rownames(DF_table_Not) <- levels(Sub)

DFShS <- plot_data[plot_data$status == "ShS",]
DFShS <- DFShS[order(DFShS$celltype),]
DF_table_ShS <- matrix(ncol = length(genes), nrow= length(levels(Sub)))
for(i in 1:length(levels(Sub))){
  row <- DFShS$pct.exp[start[i]:end[i]]
  DF_table_ShS[i,] <- row
}
colnames(DF_table_ShS) <- genes
rownames(DF_table_ShS) <- levels(Sub)

rm(i)

ShSCereb <- (DF_table_ShS)/(DF_table_Not+1)
ShSCereb[is.infinite(ShSCereb)] <-NaN 
shsmt <- matrix(nrow=1, ncol = ncol(ShSCereb))
for(i in 1:ncol(ShSCereb)){
  mx <- max(na.omit(ShSCereb[,i]))
  shsmt[i] <- mx
}
colnames(shsmt) = colnames(ShSCereb)
colnames(shsmt)[shsmt >=2]

SocCereb <- (DF_table_Soc)/(DF_table_Not+1)
SocCereb[is.infinite(SocCereb)] <-NaN 
Socmt <- matrix(nrow=1, ncol = ncol(SocCereb))
for(i in 1:ncol(SocCereb)){
  mx <- max(na.omit(SocCereb[,i]))
  Socmt[i] <- mx
}
colnames(Socmt) = colnames(SocCereb)
colnames(Socmt)[Socmt >=2]