##Supp. Figure 2&3 - Social Cerebellum
#Hailee Walker 07182023
#Finds the proportional increase in cells expressing specific panels of IEGs and plots them as a colored dotplot
#You must also have the directory to the functions sourced below entered correctly for this script to run.
#additionally the shuffling portion of this script used to determine p-values is very time consuming
#Also allows you to make downsampled UMAPs to visualize more accurately the difference in proportion between condi

rm(list=ls())
graphics.off()

library(Seurat)
library(data.table)
library(magrittr)
library(ggplot2)
source("/home/hailee/10X/Final Soc Analysis/R functions/Early_v_Late.R")
source("/home/hailee/10X/Final Soc Analysis/R functions/mIEG_dotplot.R")
source("/home/hailee/10X/Final Soc Analysis/R functions/IEG pos distributions.R")

#Import Seurat object and subset to just WT and neurons
load("~/10X/Final Soc Analysis/Cereb_All_analysis/Cereb-All.Rdata")
Sub <- subset(CellFile, subset = genotype == "WT")
Sub <- subset(Sub, idents = c("Astro", "Oligo", "OPC", "Microglia", "Endo Mural", "Endo Stalk", "Bergmann", "Choroid","Fibroblast", "Ependymal"), invert = TRUE)
Not <- subset(Sub, subset = status == "Not")### Find cells with multiple IEGs

#Find number of cells by cell type and status in the Seurat metadata
md <- Sub@meta.data %>% as.data.table
md$celltype <- factor(md$celltype, levels = levels(Sub))
md <- md[, .N, by = c("status", "celltype")]

md <- md[order(md$celltype, md$status), ]
md <- do.call(cbind.data.frame, md)

color <- c("tomato1", "darkorange", "cyan4", "aquamarine3", "darkturquoise", "turquoise", "lightblue2")

## Supp. Figure 2 ##________________________________________________________________________________________________________________

genes <- "Fos"
#Find cells that express the gene above the 95% of the control condition expression of that gene 
mat <- matrix(ncol=length(genes), nrow= ncol(Sub))
p <- GetAssayData(object=Not, slot = "data")[genes, ]
tmp <- GetAssayData(object=Sub, slot = "data")[genes, ]
q <- (tmp > quantile(p, 0.95))
q1 <- as.integer(q)
mat[,1] <- q
name <- names(q)
multIEG <- rowSums(mat)
aa <- as.numeric(as.integer(multIEG>0))
tab <- cbind(name, aa)
pos_ids = tab[tab[,2]==1]
v1 <- subset(Sub, cells = pos_ids)

Sub1 <- Sub
nz <- colSums(table(Sub1$celltype, Sub1@meta.data$condition))
nz <-min(nz[nz>0])
Idents(Sub1) <- "status"
ds <- subset(Sub1, downsample = nz)

Idents(ds) <- "celltype"
levels(ds) <- c("Granule", "UBC", "Purkinje", "MLI1", "MLI2", "PLI", "Golgi")
Not <- subset(subset(ds, subset = genotype == "WT"), subset = status == "Not")
mat <- matrix(ncol=length(genes), nrow= ncol(ds))
for(i in 1:length(genes)){
  p <- GetAssayData(object=Not, slot = "data")[genes[i], ]
  tmp <- GetAssayData(object=ds, slot = "data")[genes[i], ]
  q <- (tmp > quantile(p, 0.95))
  q1 <- as.integer(q)
  mat[,i] <- q
}

name <- names(q)
multIEG <- rowSums(mat)
aa <- as.numeric(as.integer(multIEG>0))
tab <- cbind(name, aa)
pos_ids = tab[tab[,2]==1]
v2 <- subset(ds, cells = pos_ids)

#Make UMAP of cells that meet IEG+ requirements to overlay Grey UMAP of all cells
DimPlot(v2, split.by = "status", cols = color) + NoLegend() + xlim(-10,17) + ylim(-17,14)
DimPlot(ds, split.by = "status", cols = rep("grey", length(levels(Sub)))) + NoLegend() + xlim(-10,17) + ylim(-17,14)


#Calculate number of cells of each type that meet the 95%ile IEG requirement
mIEG <- v1@meta.data %>% as.data.table
mIEG$celltype <- factor(mIEG$celltype, levels = levels(Sub))
mIEG <- mIEG[, .N, by = c("status", "celltype")]

mIEG <- mIEG[order(mIEG$celltype, mIEG$status), ]
mIEG <- do.call(cbind.data.frame, mIEG)

#Calculate proportion of cells that meet IEG requirement by cluster and condition
val <- as.numeric(mIEG[,"N"])/as.numeric(md[,"N"])
mIEG <- cbind(mIEG, val)

##Compare the Early (10 min) and Late (35 min) conditions each to baseline (i.e. 35- Baseline)
mIEG1 <- Early_v_Late(mIEG, md, Sub)
##Split into 35 and 10 minutes to make it easier to graph distributions below
mIEG1_35min <- mIEG1[8:nrow(mIEG1),]
mIEG1_35min$celltype <- factor(mIEG1_35min$celltype, levels = levels(Sub))
mIEG1_35min <- mIEG1_35min[order(mIEG1_35min$celltype), ]
mIEG1_10min <- mIEG1[1:7,]
mIEG1_10min$celltype <- factor(mIEG1_10min$celltype, levels = levels(Sub))
mIEG1_10min <- mIEG1_10min[order(mIEG1_10min$celltype), ]
##Create distributions for each cell type and for each condition
#This is pretty slow as it involves quite a few subset steps. Biggest clusters will take longest
distrib <- distribution_IEG(Sub, md, numb.iters = 1000, IEG.pctl = 0.95, IEG.mult = 0)
distrib.35min <- distrib[["distrib.35min"]]
distrib.10min <- distrib[["distrib.10min"]]
write.csv(distrib.35min, file = paste0("/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/distrib_35min -", genes, ".csv"))
write.csv(distrib.10min, file = paste0("/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/distrib_10min", genes, ".csv"))
rm(i)
#Plot distributions by cell type with red line to show the real proportion increase value as well as caculate the p-value
pval.35min <- matrix(ncol = 1, nrow = nrow(distrib.35min))
pval.10min <- matrix(ncol = 1, nrow = nrow(distrib.10min))
for(i in 1:nrow(distrib.35min)){
  rm(tfmin, pval, tenmin, pval1)
  #tfmin <- density(na.omit(distrib.35min[i,]))
  #plot(tfmin, xlim = c(0,max(as.numeric(mIEG1_35min$val)[is.finite(as.numeric(mIEG1_35min$val))])), main = paste0(rownames(distrib.35min)[i], " 35 min"))
  #abline(v = as.numeric(mIEG1_35min$val[i]), col = "red")
  pval <- length(distrib.35min[i,][(distrib.35min[i,] > as.numeric(mIEG1_35min$val[i]))])/length(distrib.35min[i,])
  pval[pval == 0] =0.001
  pval.35min[i,] = pval
  #tenmin <- density(na.omit(distrib.10min[i,]))
  #plot(tenmin, xlim = c(0,max(as.numeric(mIEG1_10min$val)[is.finite(as.numeric(mIEG1_10min$val))])), main = paste(rownames(distrib.10min)[i], " 10 min"))
  #abline(v = as.numeric(mIEG1_10min$val[i]), col = "red")
  pval1 <- length(distrib.10min[i,][(distrib.10min[i,] > as.numeric(mIEG1_10min$val[i]))])/length(distrib.10min[i,])
  pval1[pval1 == 0] =0.001
  pval.10min[i,] = pval1
}
rownames(pval.35min) <- rownames(distrib.35min)
rownames(pval.10min) <- rownames(distrib.10min)
pval.35min <- rev(pval.35min)
pval.10min <- rev(pval.10min)
pvals <- c(pval.10min,pval.35min)
mIEG1 <- cbind(mIEG1, pvals)
write.csv(mIEG1, file = paste0("/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/", genes, "dp vals.csv"))

##Write Dotplot of proportions to show difference between social groups and baseline expression
list1 <- mIEG_dotplot(mIEG, mIEG1)
dot_plot1 <- list1[["dot_plot1"]]
pdf(paste0("/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/",genes,"-DotPlot.pdf"))
print(dot_plot1)
dev.off()

genes <- "Arc"
#Find cells that express the gene above the 95% of the control condition expression of that gene 
mat <- matrix(ncol=length(genes), nrow= ncol(Sub))
p <- GetAssayData(object=Not, slot = "data")[genes, ]
tmp <- GetAssayData(object=Sub, slot = "data")[genes, ]
q <- (tmp > quantile(p, 0.95))
q1 <- as.integer(q)
mat[,1] <- q
name <- names(q)
multIEG <- rowSums(mat)
aa <- as.numeric(as.integer(multIEG>0))
tab <- cbind(name, aa)
pos_ids = tab[tab[,2]==1]
v1 <- subset(Sub, cells = pos_ids)

#Calculate number of cells of each type that meet the 95%ile IEG requirement
mIEG <- v1@meta.data %>% as.data.table
mIEG$celltype <- factor(mIEG$celltype, levels = levels(Sub))
mIEG <- mIEG[, .N, by = c("status", "celltype")]
zero_row <- cbind("ShS", "Purkinje", "0")
zero_row1 <- cbind("Not", "Purkinje", "0")
zero_row2 <- cbind("ShS", "MLI2", "0")
mIEG <- rbind(mIEG, zero_row, zero_row1, zero_row2, use.names = FALSE)

mIEG <- mIEG[order(mIEG$celltype, mIEG$status), ]
mIEG <- do.call(cbind.data.frame, mIEG)

#Calculate proportion of cells that meet IEG requirement by cluster and condition
val <- as.numeric(mIEG[,"N"])/as.numeric(md[,"N"])
mIEG <- cbind(mIEG, val)

##Compare the Early (10 min) and Late (35 min) conditions each to baseline (i.e. 35- Baseline)
mIEG1 <- Early_v_Late(mIEG, md, Sub)
##Split into 35 and 10 minutes to make it easier to graph distributions below
mIEG1_35min <- mIEG1[8:nrow(mIEG1),]
mIEG1_35min$celltype <- factor(mIEG1_35min$celltype, levels = levels(Sub))
mIEG1_35min <- mIEG1_35min[order(mIEG1_35min$celltype), ]
mIEG1_10min <- mIEG1[1:7,]
mIEG1_10min$celltype <- factor(mIEG1_10min$celltype, levels = levels(Sub))
mIEG1_10min <- mIEG1_10min[order(mIEG1_10min$celltype), ]
##Create distributions for each cell type and for each condition
#This is pretty slow as it involves quite a few subset steps. Biggest clusters will take longest
distrib <- distribution_IEG(Sub, md, numb.iters = 1000, IEG.pctl = 0.95, IEG.mult = 0)
distrib.35min <- distrib[["distrib.35min"]]
distrib.10min <- distrib[["distrib.10min"]]
write.csv(distrib.35min, file = paste0("/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/distrib_35min -", genes, ".csv"))
write.csv(distrib.10min, file = paste0("/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/distrib_10min", genes, ".csv"))
rm(i)
#Plot distributions by cell type with red line to show the real proportion increase value as well as caculate the p-value
pval.35min <- matrix(ncol = 1, nrow = nrow(distrib.35min))
pval.10min <- matrix(ncol = 1, nrow = nrow(distrib.10min))
for(i in 1:nrow(distrib.35min)){
  rm(tfmin, pval, tenmin, pval1)
  #tfmin <- density(na.omit(distrib.35min[i,]))
  #plot(tfmin, xlim = c(0,max(as.numeric(mIEG1_35min$val)[is.finite(as.numeric(mIEG1_35min$val))])), main = paste0(rownames(distrib.35min)[i], " 35 min"))
  #abline(v = as.numeric(mIEG1_35min$val[i]), col = "red")
  pval <- length(distrib.35min[i,][(distrib.35min[i,] > as.numeric(mIEG1_35min$val[i]))])/length(distrib.35min[i,])
  pval[pval == 0] =0.001
  pval.35min[i,] = pval
  #tenmin <- density(na.omit(distrib.10min[i,]))
  #plot(tenmin, xlim = c(0,max(as.numeric(mIEG1_10min$val)[is.finite(as.numeric(mIEG1_10min$val))])), main = paste(rownames(distrib.10min)[i], " 10 min"))
  #abline(v = as.numeric(mIEG1_10min$val[i]), col = "red")
  pval1 <- length(distrib.10min[i,][(distrib.10min[i,] > as.numeric(mIEG1_10min$val[i]))])/length(distrib.10min[i,])
  pval1[pval1 == 0] =0.001
  pval.10min[i,] = pval1
}
rownames(pval.35min) <- rownames(distrib.35min)
rownames(pval.10min) <- rownames(distrib.10min)
pval.35min <- rev(pval.35min)
pval.10min <- rev(pval.10min)
pvals <- c(pval.10min,pval.35min)
mIEG1 <- cbind(mIEG1, pvals)
write.csv(mIEG1, file = paste0("/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/", genes, "dp vals.csv"))

##Write Dotplot of proportions to show difference between social groups and baseline expression
list1 <- mIEG_dotplot(mIEG, mIEG1)
dot_plot1 <- list1[["dot_plot1"]]
pdf(paste0("/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/",genes,"-DotPlot.pdf"))
print(dot_plot1)
dev.off()

genes <- "Nr4a1"
#Find cells that express the gene above the 95% of the control condition expression of that gene 
mat <- matrix(ncol=length(genes), nrow= ncol(Sub))
p <- GetAssayData(object=Not, slot = "data")[genes, ]
tmp <- GetAssayData(object=Sub, slot = "data")[genes, ]
q <- (tmp > quantile(p, 0.95))
q1 <- as.integer(q)
mat[,1] <- q
name <- names(q)
multIEG <- rowSums(mat)
aa <- as.numeric(as.integer(multIEG>0))
tab <- cbind(name, aa)
pos_ids = tab[tab[,2]==1]
v1 <- subset(Sub, cells = pos_ids)

#Calculate number of cells of each type that meet the 95%ile IEG requirement
mIEG <- v1@meta.data %>% as.data.table
mIEG$celltype <- factor(mIEG$celltype, levels = levels(Sub))
mIEG <- mIEG[, .N, by = c("status", "celltype")]
zero_row <- cbind("ShS", "Purkinje", "0")
zero_row1 <- cbind("Not", "Purkinje", "0")
mIEG <- rbind(mIEG, zero_row, zero_row1, use.names = FALSE)

mIEG <- mIEG[order(mIEG$celltype, mIEG$status), ]
mIEG <- do.call(cbind.data.frame, mIEG)

#Calculate proportion of cells that meet IEG requirement by cluster and condition
val <- as.numeric(mIEG[,"N"])/as.numeric(md[,"N"])
mIEG <- cbind(mIEG, val)

##Compare the Early (10 min) and Late (35 min) conditions each to baseline (i.e. 35- Baseline)
mIEG1 <- Early_v_Late(mIEG, md, Sub)
##Split into 35 and 10 minutes to make it easier to graph distributions below
mIEG1_35min <- mIEG1[8:nrow(mIEG1),]
mIEG1_35min$celltype <- factor(mIEG1_35min$celltype, levels = levels(Sub))
mIEG1_35min <- mIEG1_35min[order(mIEG1_35min$celltype), ]
mIEG1_10min <- mIEG1[1:7,]
mIEG1_10min$celltype <- factor(mIEG1_10min$celltype, levels = levels(Sub))
mIEG1_10min <- mIEG1_10min[order(mIEG1_10min$celltype), ]
##Create distributions for each cell type and for each condition
#This is pretty slow as it involves quite a few subset steps. Biggest clusters will take longest
distrib <- distribution_IEG(Sub, md, numb.iters = 1000, IEG.pctl = 0.95, IEG.mult = 0)
distrib.35min <- distrib[["distrib.35min"]]
distrib.10min <- distrib[["distrib.10min"]]
write.csv(distrib.35min, file = paste0("/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/distrib_35min -", genes, ".csv"))
write.csv(distrib.10min, file = paste0("/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/distrib_10min", genes, ".csv"))
rm(i)
#Plot distributions by cell type with red line to show the real proportion increase value as well as caculate the p-value
pval.35min <- matrix(ncol = 1, nrow = nrow(distrib.35min))
pval.10min <- matrix(ncol = 1, nrow = nrow(distrib.10min))
for(i in 1:nrow(distrib.35min)){
  rm(tfmin, pval, tenmin, pval1)
  #tfmin <- density(na.omit(distrib.35min[i,]))
  #plot(tfmin, xlim = c(0,max(as.numeric(mIEG1_35min$val)[is.finite(as.numeric(mIEG1_35min$val))])), main = paste0(rownames(distrib.35min)[i], " 35 min"))
  #abline(v = as.numeric(mIEG1_35min$val[i]), col = "red")
  pval <- length(distrib.35min[i,][(distrib.35min[i,] > as.numeric(mIEG1_35min$val[i]))])/length(distrib.35min[i,])
  pval[pval == 0] =0.001
  pval.35min[i,] = pval
  #tenmin <- density(na.omit(distrib.10min[i,]))
  #plot(tenmin, xlim = c(0,max(as.numeric(mIEG1_10min$val)[is.finite(as.numeric(mIEG1_10min$val))])), main = paste(rownames(distrib.10min)[i], " 10 min"))
  #abline(v = as.numeric(mIEG1_10min$val[i]), col = "red")
  pval1 <- length(distrib.10min[i,][(distrib.10min[i,] > as.numeric(mIEG1_10min$val[i]))])/length(distrib.10min[i,])
  pval1[pval1 == 0] =0.001
  pval.10min[i,] = pval1
}
rownames(pval.35min) <- rownames(distrib.35min)
rownames(pval.10min) <- rownames(distrib.10min)
pval.35min <- rev(pval.35min)
pval.10min <- rev(pval.10min)
pvals <- c(pval.10min,pval.35min)
mIEG1 <- cbind(mIEG1, pvals)
write.csv(mIEG1, file = paste0("/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/", genes, "dp vals.csv"))

##Write Dotplot of proportions to show difference between social groups and baseline expression
list1 <- mIEG_dotplot(mIEG, mIEG1)
dot_plot1 <- list1[["dot_plot1"]]
pdf(paste0("/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/",genes,"-DotPlot.pdf"))
print(dot_plot1)
dev.off()

genes <- "Nr4a3"
#Find cells that express the gene above the 95% of the control condition expression of that gene 
mat <- matrix(ncol=length(genes), nrow= ncol(Sub))
p <- GetAssayData(object=Not, slot = "data")[genes, ]
tmp <- GetAssayData(object=Sub, slot = "data")[genes, ]
q <- (tmp > quantile(p, 0.95))
q1 <- as.integer(q)
mat[,1] <- q
name <- names(q)
multIEG <- rowSums(mat)
aa <- as.numeric(as.integer(multIEG>0))
tab <- cbind(name, aa)
pos_ids = tab[tab[,2]==1]
v1 <- subset(Sub, cells = pos_ids)

#Calculate number of cells of each type that meet the 95%ile IEG requirement
mIEG <- v1@meta.data %>% as.data.table
mIEG$celltype <- factor(mIEG$celltype, levels = levels(Sub))
mIEG <- mIEG[, .N, by = c("status", "celltype")]
mIEG <- mIEG[order(mIEG$celltype, mIEG$status), ]
mIEG <- do.call(cbind.data.frame, mIEG)

#Calculate proportion of cells that meet IEG requirement by cluster and condition
val <- as.numeric(mIEG[,"N"])/as.numeric(md[,"N"])
mIEG <- cbind(mIEG, val)

##Compare the Early (10 min) and Late (35 min) conditions each to baseline (i.e. 35- Baseline)
mIEG1 <- Early_v_Late(mIEG, md, Sub)
##Split into 35 and 10 minutes to make it easier to graph distributions below
mIEG1_35min <- mIEG1[8:nrow(mIEG1),]
mIEG1_35min$celltype <- factor(mIEG1_35min$celltype, levels = levels(Sub))
mIEG1_35min <- mIEG1_35min[order(mIEG1_35min$celltype), ]
mIEG1_10min <- mIEG1[1:7,]
mIEG1_10min$celltype <- factor(mIEG1_10min$celltype, levels = levels(Sub))
mIEG1_10min <- mIEG1_10min[order(mIEG1_10min$celltype), ]
##Create distributions for each cell type and for each condition
#This is pretty slow as it involves quite a few subset steps. Biggest clusters will take longest
distrib <- distribution_IEG(Sub, md, numb.iters = 1000, IEG.pctl = 0.95, IEG.mult = 0)
distrib.35min <- distrib[["distrib.35min"]]
distrib.10min <- distrib[["distrib.10min"]]
write.csv(distrib.35min, file = paste0("/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/distrib_35min -", genes, ".csv"))
write.csv(distrib.10min, file = paste0("/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/distrib_10min", genes, ".csv"))
rm(i)
#Plot distributions by cell type with red line to show the real proportion increase value as well as caculate the p-value
pval.35min <- matrix(ncol = 1, nrow = nrow(distrib.35min))
pval.10min <- matrix(ncol = 1, nrow = nrow(distrib.10min))
for(i in 1:nrow(distrib.35min)){
  rm(tfmin, pval, tenmin, pval1)
  #tfmin <- density(na.omit(distrib.35min[i,]))
  #plot(tfmin, xlim = c(0,max(as.numeric(mIEG1_35min$val)[is.finite(as.numeric(mIEG1_35min$val))])), main = paste0(rownames(distrib.35min)[i], " 35 min"))
  #abline(v = as.numeric(mIEG1_35min$val[i]), col = "red")
  pval <- length(distrib.35min[i,][(distrib.35min[i,] > as.numeric(mIEG1_35min$val[i]))])/length(distrib.35min[i,])
  pval[pval == 0] =0.001
  pval.35min[i,] = pval
  #tenmin <- density(na.omit(distrib.10min[i,]))
  #plot(tenmin, xlim = c(0,max(as.numeric(mIEG1_10min$val)[is.finite(as.numeric(mIEG1_10min$val))])), main = paste(rownames(distrib.10min)[i], " 10 min"))
  #abline(v = as.numeric(mIEG1_10min$val[i]), col = "red")
  pval1 <- length(distrib.10min[i,][(distrib.10min[i,] > as.numeric(mIEG1_10min$val[i]))])/length(distrib.10min[i,])
  pval1[pval1 == 0] =0.001
  pval.10min[i,] = pval1
}
rownames(pval.35min) <- rownames(distrib.35min)
rownames(pval.10min) <- rownames(distrib.10min)
pval.35min <- rev(pval.35min)
pval.10min <- rev(pval.10min)
pvals <- c(pval.10min,pval.35min)
mIEG1 <- cbind(mIEG1, pvals)
write.csv(mIEG1, file = paste0("/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/", genes, "dp vals.csv"))

##Write Dotplot of proportions to show difference between social groups and baseline expression
list1 <- mIEG_dotplot(mIEG, mIEG1)
dot_plot1 <- list1[["dot_plot1"]]
pdf(paste0("/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/",genes,"-DotPlot.pdf"))
print(dot_plot1)
dev.off()

## Supp Figure 3 ##_______________________________________________________________________________________________________________

genes <- c("Arc", "Dusp1", "Fos","Ier2", "Junb","Klf2", "Nfkbid", "Nr4a1", "Nr4a2", "Nr4a3", "Per1", "Pim1",
                    "Plk2", "Sik1", "Arhgef3", "Bdnf", "Inhba", "Nfib")
mat <- matrix(ncol=length(genes), nrow= ncol(Sub))
for(i in 1:length(genes)){
  p <- GetAssayData(object=Not, slot = "data")[genes[i], ]
  tmp <- GetAssayData(object=Sub, slot = "data")[genes[i], ]
  q <- (tmp > quantile(p, 0.95))
  q1 <- as.integer(q)
  mat[,i] <- q
  
}
name <- names(q)
multIEG <- rowSums(mat)
aa <- as.numeric(as.integer(multIEG>1))
tab <- cbind(name, aa)
pos_ids = tab[tab[,2]==1]
v1 <- subset(Sub, cells = pos_ids)

##Count the number of cells that belong to each cluster by condition
md <- Sub@meta.data %>% as.data.table
md$celltype <- factor(md$celltype, levels = levels(Sub))
md <- md[, .N, by = c("status", "celltype")]

md <- md[order(md$celltype, md$status), ]
md <- do.call(cbind.data.frame, md)

write.csv(md, file = "/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/md.csv")

##Count number of cells that meets the IEG requirements by cell type and condition
mIEG <- v1@meta.data %>% as.data.table
mIEG$celltype <- factor(mIEG$celltype, levels = levels(Sub))
mIEG <- mIEG[, .N, by = c("status", "celltype")]
zero_row1 <- cbind ("Not", "UBC", "0")
mIEG <- rbind(mIEG, zero_row1, use.names = FALSE)

mIEG <- mIEG[order(mIEG$celltype, mIEG$status), ]
mIEG <- do.call(cbind.data.frame, mIEG)

write.csv(mIEG, file = "/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/mIEG.csv")

##Calculate proportion of IEG+ cells/total cells by cluster and condition
val <- as.numeric(mIEG[,"N"])/as.numeric(md[,"N"])
mIEG <- cbind(mIEG, val)

#Supp. Figure 3D
##Compare the Early (10 min) and Late (35 min) conditions each to baseline (i.e. 35- Baseline)
mIEG1 <- Early_v_Late(mIEG, md, Sub)
##Split into 35 and 10 minutes to make it easier to graph distributions below
mIEG1_35min <- mIEG1[8:nrow(mIEG1),]
mIEG1_35min$celltype <- factor(mIEG1_35min$celltype, levels = levels(Sub))
mIEG1_35min <- mIEG1_35min[order(mIEG1_35min$celltype), ]
mIEG1_10min <- mIEG1[1:7,]
mIEG1_10min$celltype <- factor(mIEG1_10min$celltype, levels = levels(Sub))
mIEG1_10min <- mIEG1_10min[order(mIEG1_10min$celltype), ]
##Create distributions for each cell type and for each condition
#This is pretty slow as it involves quite a few subset steps. Biggest clusters will take longest
distrib <- distribution_IEG(Sub, md, numb.iters = 1000, IEG.pctl = 0.95, IEG.mult = 1)
distrib.35min <- distrib[["distrib.35min"]]
distrib.10min <- distrib[["distrib.10min"]]
write.csv(distrib.35min, file = "/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/distrib_35min.csv")
write.csv(distrib.10min, file = "/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/distrib_10min.csv")
rm(i)
#Plot distributions by cell type with red line to show the real proportion increase value (Subset of plots in supplemental 5)
pval.35min <- matrix(ncol = 1, nrow = nrow(distrib.35min))
pval.10min <- matrix(ncol = 1, nrow = nrow(distrib.10min))
for(i in 1:nrow(distrib.35min)){
  rm(tfmin, pval, tenmin, pval1)
  tfmin <- density(na.omit(distrib.35min[i,]))
  plot(tfmin, xlim = c(0,max(as.numeric(mIEG1_35min$val)[is.finite(as.numeric(mIEG1_35min$val))])), main = paste0(rownames(distrib.35min)[i], " 35 min"))
  abline(v = as.numeric(mIEG1_35min$val[i]), col = "red")
  pval <- length(distrib.35min[i,][(distrib.35min[i,] > as.numeric(mIEG1_35min$val[i]))])/length(distrib.35min[i,])
  pval[pval == 0] =0.001
  pval.35min[i,] = pval
  tenmin <- density(na.omit(distrib.10min[i,]))
  plot(tenmin, xlim = c(0,max(as.numeric(mIEG1_10min$val)[is.finite(as.numeric(mIEG1_10min$val))])), main = paste(rownames(distrib.10min)[i], " 10 min"))
  abline(v = as.numeric(mIEG1_10min$val[i]), col = "red")
  pval1 <- length(distrib.10min[i,][(distrib.10min[i,] > as.numeric(mIEG1_10min$val[i]))])/length(distrib.10min[i,])
  pval1[pval1 == 0] =0.001
  pval.10min[i,] = pval1
}
rownames(pval.35min) <- rownames(distrib.35min)
rownames(pval.10min) <- rownames(distrib.10min)
pval.35min <- rev(pval.35min)
pval.10min <- rev(pval.10min)
pvals <- c(pval.10min,pval.35min)
mIEG1 <- cbind(mIEG1, pvals)
##Supplementary values for dotplot
write.csv(mIEG1, file = "/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/panel_dp_vals.csv")
write.csv(mIEG, file = "/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/panel_dp_props.csv")
write.csv(md, file = "/home/hailee/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/md.csv")

##Write Dotplot of proportions to show difference between social groups and baseline expression
list1 <- mIEG_dotplot(mIEG, mIEG1)
dot_plot1 <- list1[["dot_plot1"]]
dot_plot1

## Figure 5B&C ## _______________________________________________________________________________________________________________
color <- c("tomato1", "darkorange", "cyan4", "aquamarine3", "darkturquoise", "turquoise", "lightblue2")

Sub1 <- Sub
nz <- colSums(table(Sub1$celltype, Sub1@meta.data$condition))
nz <-min(nz[nz>0])
Idents(Sub1) <- "status"
ds <- subset(Sub1, downsample = nz)

Idents(ds) <- "celltype"
levels(ds) <- c("Granule", "UBC", "Purkinje", "MLI1", "MLI2", "PLI", "Golgi")
Not <- subset(subset(ds, subset = genotype == "WT"), subset = status == "Not")
mat <- matrix(ncol=length(genes), nrow= ncol(ds))
for(i in 1:length(genes)){
  p <- GetAssayData(object=Not, slot = "data")[genes[i], ]
  tmp <- GetAssayData(object=ds, slot = "data")[genes[i], ]
  q <- (tmp > quantile(p, 0.95))
  q1 <- as.integer(q)
  mat[,i] <- q
}

name <- names(q)
multIEG <- rowSums(mat)
aa <- as.numeric(as.integer(multIEG>1))#change this value to zero to get Supp fig. 3B
tab <- cbind(name, aa)
pos_ids = tab[tab[,2]==1]
v2 <- subset(ds, cells = pos_ids)

#Make UMAP of cells that meet IEG+ requirements to overlay Grey UMAP of all cells
DimPlot(v2, split.by = "status", cols = color) + NoLegend() + xlim(-10,17) + ylim(-17,14)
DimPlot(ds, split.by = "status", cols = rep("grey", length(levels(Sub)))) + NoLegend() + xlim(-10,17) + ylim(-17,14)
