##Figure 5&6- Social ensemble PFC
#Hailee Walker - 07182023
#Finds the proportional increase in cells expressing specific panels of IEGs and plots them as a colored dotplot
  #You must also have the directory to the functions sourced below entered correctly for this script to run.
  #additionally the shuffling portion of this script used to determine p-values is very time consuming
#Also allows you to make downsampled UMAPs to visualize more accurately the difference in proportion between conditions.

rm(list =ls())
graphics.off()

library(Seurat)
library(data.table)
library(magrittr)
library(ggplot2)
source("/home/hailee/10X/Final Soc Analysis/R functions/Early_v_Late.R")
source("/home/hailee/10X/Final Soc Analysis/R functions/mIEG_dotplot.R")
source("/home/hailee/10X/Final Soc Analysis/R functions/IEG pos distributions.R")

#Import Seurat object and subset to just WT and neurons
load("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/Comb_PFC.Rdata")
Sub <- subset(CellFile, subset = genotype == "WT")
Sub <- subset(Sub, idents = c("Astro", "Oligo", "OPC", "Microglia", "Endothelial", "VLMC", "Pericyte"), invert = TRUE)
Not <- subset(Sub, subset = status == "Not")

#Find number of cells by cell type and status in the Seurat metadata
md <- Sub@meta.data %>% as.data.table
md$celltype <- factor(md$celltype, levels = levels(Sub))
md <- md[, .N, by = c("status", "celltype")]

md <- md[order(md$celltype, md$status), ]
md <- do.call(cbind.data.frame, md)

color <- c("brown2", "tomato", "darkorange1", "darkorange2", "darkgoldenrod2", "darkgoldenrod3", "yellow3", "tan2", "deepskyblue", 
           "aquamarine4", "cyan4", "aquamarine3", "darkturquoise", "cornflowerblue", "lightskyblue")

## Figure 5 ## _____________________________________________________________________________

genes <- "Fos"

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
levels(ds) <- c("L2/3_Exc", "L2/3/5_Exc", "L5_Exc", "L5_NP_CTX" ,"L6_Exc_Syt6", "L6b_CTX", "L6_Exc_Oprk1", "L6 Car3", "PV", "Som", "Som_Chodl",
                "VIP", "Lamp5", "Meis2,Pbx3-1", "Meis2,Pbx3-2")
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

DimPlot(v2, split.by = "status", cols = color) +NoLegend() + xlim(-18,18) + ylim(-18, 18)
DimPlot(ds, split.by = "status", cols = rep("grey", 15))+NoLegend() + xlim(-18,18) + ylim(-18, 18)

mIEG <- v1@meta.data %>% as.data.table
mIEG$celltype <- factor(mIEG$celltype, levels = levels(Sub))
mIEG <- mIEG[, .N, by = c("status", "celltype")]

mIEG <- mIEG[order(mIEG$celltype, mIEG$status), ]
mIEG <- do.call(cbind.data.frame, mIEG)

val <- as.numeric(mIEG[,"N"])/as.numeric(md[,"N"])
mIEG <- cbind(mIEG, val)

##Compare the Early (10 min) and Late (35 min) conditions each to baseline (i.e. 35- Baseline)
mIEG1 <- Early_v_Late(mIEG, md, Sub)
##Split into 35 and 10 minutes to make it easier to graph distributions below
mIEG1_35min <- mIEG1[16:nrow(mIEG1),]
mIEG1_35min$celltype <- factor(mIEG1_35min$celltype, levels = levels(Sub))
mIEG1_35min <- mIEG1_35min[order(mIEG1_35min$celltype), ]
mIEG1_10min <- mIEG1[1:15,]
mIEG1_10min$celltype <- factor(mIEG1_10min$celltype, levels = levels(Sub))
mIEG1_10min <- mIEG1_10min[order(mIEG1_10min$celltype), ]
##Create distributions for each cell type and for each condition
#This is pretty slow as it involves quite a few subset steps. Biggest clusters will take longest
distrib <- distribution_IEG(Sub, md, numb.iters = 1000, IEG.pctl = 0.95, IEG.mult = 0)
distrib.35min <- distrib[["distrib.35min"]]
distrib.10min <- distrib[["distrib.10min"]]
write.csv(distrib.35min, file = paste0("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/distrib_35min -", genes, "-2.csv"))
write.csv(distrib.10min, file = paste0("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/distrib_10min", genes, "-2.csv"))
rm(i)
#Plot distributions by cell type with red line to show the real proportion increase value
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
write.csv(mIEG1, file = paste0("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/", genes, "dp vals-2.csv"))

##Write Dotplot of proportions to show difference between social groups and baseline expression
list1 <- mIEG_dotplot(mIEG, mIEG1)
dot_plot1 <- list1[["dot_plot1"]]
pdf(paste0("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/",genes,"-DotPlot-2.pdf"))
print(dot_plot1)
dev.off()

genes <- "Arc"
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

DimPlot(v1, split.by = "status", cols = color) +NoLegend() + xlim(-18,18) + ylim(-18, 18)

mIEG <- v1@meta.data %>% as.data.table
mIEG$celltype <- factor(mIEG$celltype, levels = levels(Sub))
mIEG <- mIEG[, .N, by = c("status", "celltype")]
zero_row <- cbind("Not", "Lamp5", "0")
mIEG <- rbind(mIEG, zero_row, use.names = FALSE)

mIEG <- mIEG[order(mIEG$celltype, mIEG$status), ]
mIEG <- do.call(cbind.data.frame, mIEG)

val <- as.numeric(mIEG[,"N"])/as.numeric(md[,"N"])
mIEG <- cbind(mIEG, val)

##Compare the Early (10 min) and Late (35 min) conditions each to baseline (i.e. 35- Baseline)
mIEG1 <- Early_v_Late(mIEG, md, Sub)
##Split into 35 and 10 minutes to make it easier to graph distributions below
mIEG1_35min <- mIEG1[16:nrow(mIEG1),]
mIEG1_35min$celltype <- factor(mIEG1_35min$celltype, levels = levels(Sub))
mIEG1_35min <- mIEG1_35min[order(mIEG1_35min$celltype), ]
mIEG1_10min <- mIEG1[1:15,]
mIEG1_10min$celltype <- factor(mIEG1_10min$celltype, levels = levels(Sub))
mIEG1_10min <- mIEG1_10min[order(mIEG1_10min$celltype), ]
##Create distributions for each cell type and for each condition
#This is pretty slow as it involves quite a few subset steps. Biggest clusters will take longest
distrib <- distribution_IEG(Sub, md, numb.iters = 1000, IEG.pctl = 0.95, IEG.mult = 0)
distrib.35min <- distrib[["distrib.35min"]]
distrib.10min <- distrib[["distrib.10min"]]
write.csv(distrib.35min, file = paste0("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/distrib_35min -", genes, ".csv"))
write.csv(distrib.10min, file = paste0("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/distrib_10min", genes, ".csv"))
rm(i)
#Plot distributions by cell type with red line to show the real proportion increase value
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
write.csv(mIEG1, file = paste0("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/", genes, "dp vals.csv"))

##Write Dotplot of proportions to show difference between social groups and baseline expression
list1 <- mIEG_dotplot(mIEG, mIEG1)
dot_plot1 <- list1[["dot_plot1"]]
pdf(paste0("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/",genes,"-DotPlot.pdf"))
print(dot_plot1)
dev.off()

genes <- "Nr4a1"
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

DimPlot(v1, split.by = "status", cols = color) +NoLegend() + xlim(-18,18) + ylim(-18, 18)

mIEG <- v1@meta.data %>% as.data.table
mIEG$celltype <- factor(mIEG$celltype, levels = levels(Sub))
mIEG <- mIEG[, .N, by = c("status", "celltype")]
zero_row <- cbind("Not", "L6 Car3", "0")
zero_row1 <- cbind("Soc", "L6 Car3", "0")
zero_row2 <- cbind("Soc", "Som", "0")
mIEG <- rbind(mIEG, zero_row, zero_row1, zero_row2, use.names = FALSE)

mIEG <- mIEG[order(mIEG$celltype, mIEG$status), ]
mIEG <- do.call(cbind.data.frame, mIEG)

val <- as.numeric(mIEG[,"N"])/as.numeric(md[,"N"])
mIEG <- cbind(mIEG, val)

##Compare the Early (10 min) and Late (35 min) conditions each to baseline (i.e. 35- Baseline)
mIEG1 <- Early_v_Late(mIEG, md, Sub)
##Split into 35 and 10 minutes to make it easier to graph distributions below
mIEG1_35min <- mIEG1[16:nrow(mIEG1),]
mIEG1_35min$celltype <- factor(mIEG1_35min$celltype, levels = levels(Sub))
mIEG1_35min <- mIEG1_35min[order(mIEG1_35min$celltype), ]
mIEG1_10min <- mIEG1[1:15,]
mIEG1_10min$celltype <- factor(mIEG1_10min$celltype, levels = levels(Sub))
mIEG1_10min <- mIEG1_10min[order(mIEG1_10min$celltype), ]
##Create distributions for each cell type and for each condition
#This is pretty slow as it involves quite a few subset steps. Biggest clusters will take longest
distrib <- distribution_IEG(Sub, md, numb.iters = 1000, IEG.pctl = 0.95, IEG.mult = 0)
distrib.35min <- distrib[["distrib.35min"]]
distrib.10min <- distrib[["distrib.10min"]]
write.csv(distrib.35min, file = paste0("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/distrib_35min -", genes, ".csv"))
write.csv(distrib.10min, file = paste0("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/distrib_10min", genes, ".csv"))
rm(i)
#Plot distributions by cell type with red line to show the real proportion increase value
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
write.csv(mIEG1, file = paste0("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/", genes, "dp vals.csv"))

##Write Dotplot of proportions to show difference between social groups and baseline expression
list1 <- mIEG_dotplot(mIEG, mIEG1)
dot_plot1 <- list1[["dot_plot1"]]
pdf(paste0("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/",genes,"-DotPlot.pdf"))
print(dot_plot1)
dev.off()

genes <- "Nr4a3"
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

DimPlot(v1, split.by = "status", cols = color) +NoLegend() + xlim(-18,18) + ylim(-18, 18)

mIEG <- v1@meta.data %>% as.data.table
mIEG$celltype <- factor(mIEG$celltype, levels = levels(Sub))
mIEG <- mIEG[, .N, by = c("status", "celltype")]
zero_row2 <- cbind("Not", "Som", "0")
mIEG <- rbind(mIEG, zero_row2, use.names = FALSE)

mIEG <- mIEG[order(mIEG$celltype, mIEG$status), ]
mIEG <- do.call(cbind.data.frame, mIEG)

val <- as.numeric(mIEG[,"N"])/as.numeric(md[,"N"])
mIEG <- cbind(mIEG, val)

##Compare the Early (10 min) and Late (35 min) conditions each to baseline (i.e. 35- Baseline)
mIEG1 <- Early_v_Late(mIEG, md, Sub)
##Split into 35 and 10 minutes to make it easier to graph distributions below
mIEG1_35min <- mIEG1[16:nrow(mIEG1),]
mIEG1_35min$celltype <- factor(mIEG1_35min$celltype, levels = levels(Sub))
mIEG1_35min <- mIEG1_35min[order(mIEG1_35min$celltype), ]
mIEG1_10min <- mIEG1[1:15,]
mIEG1_10min$celltype <- factor(mIEG1_10min$celltype, levels = levels(Sub))
mIEG1_10min <- mIEG1_10min[order(mIEG1_10min$celltype), ]
##Create distributions for each cell type and for each condition
#This is pretty slow as it involves quite a few subset steps. Biggest clusters will take longest
distrib <- distribution_IEG(Sub, md, numb.iters = 1000, IEG.pctl = 0.95, IEG.mult = 0)
distrib.35min <- distrib[["distrib.35min"]]
distrib.10min <- distrib[["distrib.10min"]]
write.csv(distrib.35min, file = paste0("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/distrib_35min -", genes, ".csv"))
write.csv(distrib.10min, file = paste0("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/distrib_10min", genes, ".csv"))
rm(i)
#Plot distributions by cell type with red line to show the real proportion increase value
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
write.csv(mIEG1, file = paste0("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/", genes, "dp vals.csv"))

##Write Dotplot of proportions to show difference between social groups and baseline expression
list1 <- mIEG_dotplot(mIEG, mIEG1)
dot_plot1 <- list1[["dot_plot1"]]
pdf(paste0("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/",genes,"-DotPlot.pdf"))
print(dot_plot1)
dev.off()

##Figure 6 ## ___________________________________________________________________________________________________
genes <- c("Apold1", "Arc", "Bhlhe40", "Btg2", "Ccn1", "Ccn2", "Cebpb", "Dusp1", "Dusp2", "Dusp6", "Egr1", "Egr2", "Egr3", "Egr4", "F3", "Fos", "Fosb", 
           "Fosl2", "Gadd45b", "Gadd45g", "Gem", "Hes1", "Ier2", "Ier3", "Ier5", "Ifit2", "Ifit3", "Jun", "Junb", "Klf10", "Maff", "Nfkbia", "Npas4", 
           "Nr4a1", "Nr4a2", "Nr4a3", "Nrn1", "Per1", "Ptgs2", "Sik1", "Trib1", "Zfp36l1", "Zfp36l2") #PFC genes

mat <- matrix(ncol=length(genes), nrow= ncol(Sub))
##Find cells that express 4+ IEGs above the 95% of expression in not social for that IEG
for(i in 1:length(genes)){
  p <- GetAssayData(object=Not, slot = "data")[genes[i], ]
  tmp <- GetAssayData(object=Sub, slot = "data")[genes[i], ]
  q <- (tmp > quantile(p, 0.95))
  q1 <- as.integer(q)
  mat[,i] <- q
}
rm(Not)
name <- names(q)
multIEG <- rowSums(mat)
aa <- as.numeric(as.integer(multIEG>3))
tab <- cbind(name, aa)
pos_ids = tab[tab[,2]==1]
v1 <- subset(Sub, cells = pos_ids)

##Pull out the metadata of the total number of cells which belong to each cluster
md <- Sub@meta.data %>% as.data.table
md$celltype <- factor(md$celltype, levels = levels(Sub))
md <- md[, .N, by = c("status", "celltype")]

md <- md[order(md$celltype, md$status), ]
md <- do.call(cbind.data.frame, md)

write.csv(md, file = "/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/md.csv")

##Pull out the metadata for the total number of cells that meet the IEG requirements as specified above
mIEG <- v1@meta.data %>% as.data.table
mIEG$celltype <- factor(mIEG$celltype, levels = levels(Sub))
mIEG <- mIEG[, .N, by = c("status", "celltype")]

mIEG <- mIEG[order(mIEG$celltype, mIEG$status), ]
mIEG <- do.call(cbind.data.frame, mIEG)

write.csv(mIEG, file = "/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/mIEG.csv")

##Calculate proportion IEG+
val <- as.numeric(mIEG[,"N"])/as.numeric(md[,"N"])
mIEG <- cbind(mIEG, val)

#Figure 6D
##Compare the Early (10 min) and Late (35 min) conditions each to baseline (i.e. 35- Baseline)
mIEG1 <- Early_v_Late(mIEG, md, Sub)
##Split into 35 and 10 minutes to make it easier to graph distributions below
mIEG1_35min <- mIEG1[16:nrow(mIEG1),]
mIEG1_35min$celltype <- factor(mIEG1_35min$celltype, levels = levels(Sub))
mIEG1_35min <- mIEG1_35min[order(mIEG1_35min$celltype), ]
mIEG1_10min <- mIEG1[1:15,]
mIEG1_10min$celltype <- factor(mIEG1_10min$celltype, levels = levels(Sub))
mIEG1_10min <- mIEG1_10min[order(mIEG1_10min$celltype), ]
##Create distributions for each cell type and for each condition
#This is pretty slow as it involves quite a few subset steps. Biggest clusters will take longest
distrib <- distribution_IEG(Sub, md, numb.iters = 1000)
distrib.35min <- distrib[["distrib.35min"]]
distrib.10min <- distrib[["distrib.10min"]]
write.csv(distrib.35min, file = paste0("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/distrib_35min -panel.csv"))
write.csv(distrib.10min, file = paste0("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/WT only/distrib_10min - panel.csv"))
rm(i)
#Plot distributions by cell type with red line to show the real proportion increase value (Plots for Supplemental)
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
write.csv(mIEG1, file = "DP vals col2.csv")
write.csv(mIEG, file = "DP vals col1.csv")

##Write Dotplot of proportions to show difference between social groups and baseline expression
list1 <- mIEG_dotplot(mIEG, mIEG1)
dot_plot1 <- list1[["dot_plot1"]]
dot_plot1

## Figure 4C,B ##_________________________________________________________________________________________________________________
#UMAP of cells that meet requirements
color <- c("brown2", "tomato", "darkorange1", "darkorange2", "darkgoldenrod2", "darkgoldenrod3", "yellow3", "tan2", "deepskyblue", 
           "aquamarine4", "cyan4", "aquamarine3", "darkturquoise", "cornflowerblue", "lightskyblue")
Sub1 <- Sub
nz <- colSums(table(Sub1$celltype, Sub1@meta.data$condition))
nz <-min(nz[nz>0])
Idents(Sub1) <- "status"
ds <- subset(Sub1, downsample = nz)

Idents(ds) <- "celltype"
levels(ds) <- c("L2/3_Exc", "L2/3/5_Exc", "L5_Exc", "L5_NP_CTX" ,"L6_Exc_Syt6", "L6b_CTX", "L6_Exc_Oprk1", "L6 Car3", "PV", "Som", "Som_Chodl",
                "VIP", "Lamp5", "Meis2,Pbx3-1", "Meis2,Pbx3-2")
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
aa <- as.numeric(as.integer(multIEG>3)) #Adjust this parameter to make Figure 6B
tab <- cbind(name, aa)
pos_ids = tab[tab[,2]==1]
v2 <- subset(ds, cells = pos_ids)

DimPlot(v2, split.by = "status", cols = color) +NoLegend() + xlim(-18,18) + ylim(-18, 18)
#Grey UMAP of all Cells to overlay colored on
DimPlot(ds, split.by = "status", cols = rep("grey", length(levels(Sub)))) + NoLegend() + xlim(-18,18) + ylim(-18, 18)
