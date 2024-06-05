##Supp. Figure 4 - Make pie chart of the proportion excess for each significantly socially modulated cell type
#Hailee Walker - 07252023

rm(list = ls())
graphics.off()
library(viridis)

## Calculate the proportional excess for cerebellum
#import mIEG1 and md saved from calculations made for Figure 5C
md <- read.csv("~/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/md.csv", row.names=1)
mIEG1 <- read.csv("~/10X/Final Soc Analysis/Cereb_All_analysis/WT_only/panel_dp_vals.csv", row.names=1)

## For 35 min
mIEG1.35 <- mIEG1[8:14,]
mIEG1.35<- mIEG1.35[nrow(mIEG1.35):1,]
md.soc <- md[md$status == "Soc",]
comb <- cbind(mIEG1.35, md.soc$N)
combC.35 <- comb[(comb$pvals < 0.05),]
numPieC.35 <- combC.35$val * combC.35$`md.soc$N`

## For 10 min
mIEG1.10 <- mIEG1[1:7,]
mIEG1.10<- mIEG1.10[nrow(mIEG1.10):1,]
md.soc <- md[md$status == "ShS",]
comb <- cbind(mIEG1.10, md.soc$N)
combC.10 <- comb[(comb$pvals < 0.05),]
numPieC.10 <- combC.10$val * combC.10$`md.soc$N`

##Calculate the proportional excess for PFC
#import mIEG1 and md saved from calculations made for Figure 5C
md <- read.csv("~/10X/Final Soc Analysis/PFC_All_analysis/WT only/md.csv", row.names=1)
mIEG1 <- read.csv("~/10X/Final Soc Analysis/PFC_All_analysis/WT only/WT dotplot values- C2.csv", row.names=1)

## For 35 min
mIEG1.35 <- mIEG1[16:30,]
mIEG1.35<- mIEG1.35[nrow(mIEG1.35):1,]
md.soc <- md[md$status == "Soc",]
comb <- cbind(mIEG1.35, md.soc$N)
comb.35 <- comb[(comb$pvals < 0.05),]
numPie.35 <- comb.35$val * comb.35$`md.soc$N`

## For 10 min
mIEG1.10 <- mIEG1[1:15,]
mIEG1.10<- mIEG1.10[nrow(mIEG1.10):1,]
md.soc <- md[md$status == "ShS",]
comb <- cbind(mIEG1.10, md.soc$N)
comb.10 <- comb[(comb$pvals < 0.05),]
numPie.10 <- comb.10$val * comb.10$`md.soc$N`

##Make matrix of combined PFC and Cerebellum excess for 10 min condition
mat <- as.matrix(c(numPie.10, numPieC.10))
mat<- as.matrix(c(mat[1:16], mat[18:21])) #Had to remove the infinite value for UBC or can't make pie chart

labels <-as.matrix((c(comb.10$celltype, combC.10$celltype)))
labels <- as.matrix(c(labels[1:16], labels[18:21]))

pieval <- (mat/sum(mat))* 100

x <- substr(pieval, 1, 4) 
pie(as.numeric(pieval), labels = c(x,"%"),  cex = 0.8, radius = 1, col = inferno(length(pieval)), main = "10 minute Social Ensemble")
legend("topleft", labels, cex = 0.8, col = inferno(length(pieval)), pch =16)

##Make matrix of combined PFC and Cerebellum excess for 35 min condition
mat <- as.matrix(c(numPie.35, numPieC.35))
mat<- as.matrix(c(mat[1:12], mat[14:15])) #Had to remove the infinite value for UBC or can't make pie chart

labels <-as.matrix((c(comb.35$celltype, combC.35$celltype)))
labels <- as.matrix(c(labels[1:12], labels[14:15]))

pieval <- (mat/sum(mat))* 100

x <- substr(pieval, 1, 4) 
pie(as.numeric(pieval), labels = c(x,"%"),  cex = 0.8, radius = 1, col = inferno(length(pieval)), main = "35 minute Social Ensemble")
legend("topleft", labels, cex = 0.8, col = inferno(length(pieval)), pch =16)

