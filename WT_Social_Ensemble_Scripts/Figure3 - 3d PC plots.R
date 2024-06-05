##Figure 3 - IEG PCA
#Hailee Walker
#Uses Deseq2 to run PCA on pseudobulked reads grouped by region (PFC or Cerebellum), broad neuronal class
#(excitatory or inhibitory) and experimental group (Control, 10 min post-social or 35 min post-social)
#PCs are calculated based only on reads of the PFC+Cerebellum IEG panel

#Plots in 3D PC space
#Then PCs are normalized to the mean of the control group for each broad region-neuron type they belong to
#the baseline normalized PCs are then plotted to examine the shift away at 10 and 35 min post social.
#cosine similarity of all groups to every other group is then calculated and exported as csvs which can be used to make plots in prism

rm(list = ls())
graphics.off()

library(readr)
library(Seurat)
library(doParallel)
library(foreach)
library(DESeq2)
library(scales)
library(ggplot2)
library(dplyr)
library(rstatix)
library(lsa)

load("~/10X/Final Soc Analysis/PFC_All_analysis/Comb_PFC.Rdata")
CellFile_PFC <- CellFile
load("~/10X/Final Soc Analysis/Cereb_All_analysis/Cereb-All.Rdata")
CellFile_Cereb <- CellFile
rm(CellFile)

genes <- c("Dusp1", "Nr4a3", "Nr4a1", "Per1", "Nr4a2", "Fos", "Pim1", "Sik1", "Zfp36l2", "Cebpb", "Btg2", "Jun","Dusp6",
           "Ier2", "Gem", "Fosl2", "Apold1", "Plk2", "Gadd45g", "Inhba", "Nfib", "Arhgef3", "Bdnf", "Trib1", "Nfkbia", "Bhlhe40",
           "Hes1", "Zfp36l1", "Ifit2", "Egr1", "Klf2", "Ier3", "Ier5", "Ccn1", "Nrn1", "F3", "Ccn2", "Gadd45b", "Ifit3", "Npas4",
           "Arc", "Nfkbid", "Junb", "Klf10", "Maff", "Egr3", "Ptgs2", "Fosb", "Dusp2", "Egr2", "Egr4")

CellFile_Cereb <- subset(CellFile_Cereb, features = genes)
CellFile_Cereb <- subset(CellFile_Cereb, subset = genotype == "WT")
CellFile_PFC <- subset(CellFile_PFC, features = genes)
CellFile_PFC <- subset(CellFile_PFC, subset = genotype == "WT")
exc.pfc <- subset(CellFile_PFC, idents = c("L2/3_Exc", "L2/3/5_Exc", "L5_Exc", "L5_NP_CTX", "L6_Exc_Syt6", "L6b_CTX", "L6_Exc_Oprk1", "L6 Car3"))
exc.pfc$CT <- rep("exc.pfc", length(exc.pfc$status))
exc.pfc$actCTs <- factor(paste0(exc.pfc$CT, exc.pfc$status, exc.pfc$orig.ident))
inh.pfc <- subset(CellFile_PFC, idents = c("VIP", "PV", "Lamp5", "Som", "Som_Chodl", "Meis2,Pbx3-1", "Meis2,Pbx3-2"))
inh.pfc$CT <- rep("inh.pfc", length(inh.pfc$status))
inh.pfc$actCTs <- factor(paste0(inh.pfc$CT, inh.pfc$status, inh.pfc$orig.ident))
exc.cereb <- subset(CellFile_Cereb, idents = c("Granule", "UBC"))
exc.cereb$CT <- rep("exc.cereb", length(exc.cereb$status))
exc.cereb$actCTs <- factor(paste0(exc.cereb$CT, exc.cereb$status, exc.cereb$orig.ident))
inh.cereb <- subset(CellFile_Cereb, idents = c("Purkinje", "MLI1", "MLI2", "PLI", "Golgi"))
inh.cereb$CT <- rep("inh.cereb", length(inh.cereb$status))
inh.cereb$actCTs <- factor(paste0(inh.cereb$CT, inh.cereb$status, inh.cereb$orig.ident))
rm(CellFile_PFC, CellFile_Cereb)

Exc <- merge(exc.pfc, exc.cereb)
Inh <- merge(inh.pfc, inh.cereb)
PFC <- merge(exc.pfc, inh.pfc)
Cereb <- merge(exc.cereb, inh.cereb)
Cells <- merge(Exc, Inh)

totalCores = detectCores()
cluster <- makeCluster(totalCores[1]-25) #If nothing else is running (-1) if some is reduce cores by more (-15 or more)
registerDoParallel(cluster)

#Build count matrix by orig.ident from Seurat File
stim.list <- SplitObject(Cells, "actCTs")
y1 <- stim.list[[1]]
y2 <- y1@assays$RNA@counts
mat1 <- matrix(ncol = length(stim.list), nrow = nrow(y2) )
for(i in 1:length(stim.list)){
  print(i)
  y1 <- stim.list[[i]]
  y <- as.matrix(y1@assays$RNA@counts)
  mat2 <- foreach(c = 1:nrow(y), .combine = rbind) %dopar% {
    gene = sum(y[c,])
  }
  mat1[,i] = mat2
}
stopCluster(cluster)
nmc <- names(stim.list)
nmr <- row.names(y2)
colnames(mat1) <- nmc
rownames(mat1) <- nmr

rm(y, y1, mat2, y2)

#Run DESeq to get differential expression
Counts <- mat1[which(rowSums(mat1) > 1),] #Removing any genes that have fewer than 50 (arbitrary) reads for all samples
status <- factor(substr(names(stim.list), nchar(names(stim.list))-7,nchar(names(stim.list))-5))
CT <- factor(substr(names(stim.list), 1, 7))
coldata <- data.frame(row.names = colnames(Counts), paste0(CT, status))
colnames(coldata) <- "group"
dds <- DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~group)
dds <- DESeq(dds)
vsdata <- varianceStabilizingTransformation(dds, blind = FALSE)
my_color_palette <- hue_pal()(length(unique(coldata$group)))
plotPCA(vsdata, intgroup = "group")+ 
  scale_color_manual(values = my_color_palette)

#get table of PC data
project.pca <- prcomp(t(assay(vsdata)), scale=FALSE)
summary(project.pca)
scores = as.data.frame(project.pca$x)
scores$name <- c(paste0(CT,status))
scores <- scores[order(row.names(scores)),]

library(scatterplot3d)
#Plot PCs 1-3 in 3D PC plot - Figure 3A
colors <- cbind(unique(scores$name), my_color_palette)
colnames(colors) <- c("name", "color")
colors <- as.data.frame(colors)
scores<-scores %>%
  mutate(color = case_when(scores$name == colors$name[1] ~ colors$color[1], scores$name == colors$name[2] ~ colors$color[2], scores$name == colors$name[3] ~ colors$color[3],
                           scores$name == colors$name[4] ~ colors$color[4], scores$name == colors$name[5] ~ colors$color[5], scores$name == colors$name[6] ~ colors$color[6],
                           scores$name == colors$name[7] ~ colors$color[7], scores$name == colors$name[8] ~ colors$color[8], scores$name == colors$name[9] ~ colors$color[9],
                           scores$name == colors$name[10] ~ colors$color[10], scores$name == colors$name[11] ~ colors$color[11], scores$name == colors$name[12] ~ colors$color[12]))
scatterplot3d(scores[,1], scores[,3], scores[,2],pch = 19, color = scores$color, cex.symbols = 2)
legend("topright", legend = colors$name, col = colors$color, pch=19, cex = .7)



##Split by CT and region to identify the mean of the baseline values - Figure 3 B,C
scores_split<-split(scores, scores$name)
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
exc.cereb.bl.means <- colMeans(as.matrix(scores_split$exc.cerNot[,1:51]))
exc.pfc.bl.means <- colMeans(as.matrix(scores_split$exc.pfcNot[,1:51]))
inh.cereb.bl.means <- colMeans(as.matrix(scores_split$inh.cerNot[,1:51]))
inh.pfc.bl.means <- colMeans(as.matrix(scores_split$inh.pfcNot[,1:51]))

##Subtract the baseline from the social values and then plot them to show their variation from the mean
scores$CTR <- substr(scores$name, 1, 7)
scores_split_CTR <-split(scores[,1:51], scores$CTR)
exc.cer <- scores_split_CTR$exc.cer - rep.row(exc.cereb.bl.means, nrow(scores_split_CTR$exc.cer))
exc.pfc <- scores_split_CTR$exc.pfc - rep.row(exc.pfc.bl.means, nrow(scores_split_CTR$exc.pfc))
inh.cer <- scores_split_CTR$inh.cer - rep.row(inh.cereb.bl.means, nrow(scores_split_CTR$inh.cer))
inh.pfc <- scores_split_CTR$inh.pfc - rep.row(inh.pfc.bl.means, nrow(scores_split_CTR$inh.pfc))
bl.normed <- rbind(exc.cer, exc.pfc, inh.cer, inh.pfc)

#bl.normed$name <- c(paste0(CT,status))
#bl.normed <- bl.normed[order(row.names(bl.normed)),]
bl_split<-split(bl.normed, scores$name)

##Make table of means of each condition by PC
require(purrr)
scores.small <- lapply(bl_split, function(x) x%>% select(PC1,PC2,PC3))
means <- map_df((lapply(scores.small, colMeans)), ~as.data.frame(t(.)))
rownames(means) <- names(scores.small)

library(scatterplot3d)
#Plot PCs 1-3 in 3D PC plot
colors <- cbind(unique(scores$name), my_color_palette)
colnames(colors) <- c("name", "color")
colors <- as.data.frame(colors)
bl.normed<-bl.normed %>%
  mutate(color = case_when(scores$name == colors$name[1] ~ colors$color[1], scores$name == colors$name[2] ~ colors$color[2], scores$name == colors$name[3] ~ colors$color[3],
                           scores$name == colors$name[4] ~ colors$color[4], scores$name == colors$name[5] ~ colors$color[5], scores$name == colors$name[6] ~ colors$color[6],
                           scores$name == colors$name[7] ~ colors$color[7], scores$name == colors$name[8] ~ colors$color[8], scores$name == colors$name[9] ~ colors$color[9],
                           scores$name == colors$name[10] ~ colors$color[10], scores$name == colors$name[11] ~ colors$color[11], scores$name == colors$name[12] ~ colors$color[12]))
s3d <- scatterplot3d(bl.normed[,1], bl.normed[,3], bl.normed[,2],pch = 19, color = bl.normed$color, cex.symbols = 2)
legend("topright", legend = colors$name, col = colors$color, pch=19, cex = .7)
# Add lines connecting means
s3d <- plot(bl.normed[1:14,1], bl.normed[1:14,2], pch = 19, col = bl.normed$color[1:14], main = "Cerebellum Excitatory Neurons", xlim = c(-3,2), ylim = c(-1.5,2))
points(means[c(3,1,2),1], means[c(3,1,2),2],
             col = "black", type = "b", pch = c(4,8,6))
legend("topright", legend = c(colors$name[1:3],"CTRL-mean", "10 min.", "35 min."), col = c(colors$color[1:3], "black", "black", "black"), pch=c(19, 19, 19, 8, 6, 4))
s3d <- plot(bl.normed[15:28,1], bl.normed[15:28,2],pch = 19, col = bl.normed$color[15:28], main = "PFC Excitatory Neurons", xlim = c(-3,2), ylim = c(-1.5,2))
points(means[c(6,4,5),1], means[c(6,4,5),2],
             col = "black", type = "b", pch = c(4,8,6))
legend("topright", legend = c(colors$name[4:6],"CTRL-mean", "10 min.", "35 min."), col = c(colors$color[4:6], "black", "black", "black"), pch=c(19, 19, 19, 8, 6, 4))
s3d <- plot(bl.normed[29:42,1], bl.normed[29:42,2],pch = 19, col = bl.normed$color[29:42], main = "Cerebellum Inhibitory Neurons", xlim = c(-3,2), ylim = c(-1.5,2))
points(means[c(9,7,8),1], means[c(9,7,8),2],
             col = "black", type = "b", pch = c(4,8,6))
legend("topright", legend = c(colors$name[7:9],"CTRL-mean", "10 min.", "35 min."), col = c(colors$color[7:9], "black", "black", "black"), pch=c(19, 19, 19, 8, 6, 4))
s3d <- plot(bl.normed[43:56,1], bl.normed[43:56,2],pch = 19, col = bl.normed$color[43:56], main = "PFC Inhibitory Neurons", xlim = c(-3,2), ylim = c(-1.5,2))
points(means[c(12,10,11),1], means[c(12,10,11),2],
             col = "black", type = "b", pch = c(4,8,6))
legend("topright", legend = c(colors$name[10:12],"CTRL-mean", "10 min.", "35 min."), col = c(colors$color[10:12], "black", "black", "black"), pch=c(19, 19, 19, 8, 6, 4))

#Plot PC 1 & 2 for all groups
plot(bl.normed[,1], bl.normed[,2], col = bl.normed$color, pch = 19)

##Find Cosine similarity between samples and export values to csvs to be plotted in prism - Figure 4 D-G
status <- substr(rownames(bl.normed), nchar(rownames(bl.normed))-7, nchar(rownames(bl.normed))-5)
bl.normed.shs <- bl.normed[status == "ShS",]
bl.normed.soc <- bl.normed[status == "Soc",]
bl.normed.soc <- rbind(bl.normed.shs, bl.normed.soc)
cosine.values <- cosine(as.matrix(t(bl.normed.soc[,1:50])))

pfc10.pfc10.exc <- cosine.values[4:6,4:6]
pfc10.pfc10.exc <- c(as.matrix(pull_lower_triangle(pfc10.pfc10.exc, diagonal = TRUE))[,2:4])
# Replace lower triangular part as they are duplicate values
pfc10.pfc10.exc.nms <- rep("pfc10.pfc10.exc", length(pfc10.pfc10.exc))
pfc10.pfc35.exc <- c(cosine.values[4:6,19:23])
pfc10.pfc35.exc.nms <- rep("pfc10.pfc35.exc", length(pfc10.pfc35.exc))
pfc10.cereb10.exc <- c(cosine.values[4:6,1:3])
pfc10.cereb10.exc.nms <- rep("pfc10.cereb10.exc", length(pfc10.cereb10.exc))
pfc10.cereb35.exc <- c(cosine.values[4:6,13:18])
pfc10.cereb35.exc.nms <- rep("pfc10.cereb35.exc", length(pfc10.cereb35.exc))
pfc35.pfc35.exc <- cosine.values[19:23,19:23]
pfc35.pfc35.exc <- c(as.matrix(pull_lower_triangle(pfc35.pfc35.exc, diagonal = TRUE))[,2:6])
pfc35.pfc35.exc.nms <- rep("pfc35.pfc35.exc", length(pfc35.pfc35.exc))
pfc35.cereb10.exc <- c(cosine.values[19:23,1:3])
pfc35.cereb10.exc.nms <- rep("pfc35.cereb10.exc", length(pfc35.cereb10.exc))
pfc35.cereb35.exc <- c(cosine.values[19:23,13:18])
pfc35.cereb35.exc.nms <- rep("pfc35.cereb35.exc", length(pfc35.cereb35.exc))
cereb10.cereb10.exc<- cosine.values[1:3,1:3]
cereb10.cereb10.exc <- c(as.matrix(pull_lower_triangle(cereb10.cereb10.exc, diagonal = TRUE))[,2:4])
cereb10.cereb10.exc.nms <- rep("cereb10.cereb10.exc", length(cereb10.cereb10.exc))
cereb10.cereb35.exc <- c(cosine.values[1:3,13:18])
cereb10.cereb35.exc.nms <- rep("cereb10.cereb35.exc", length(cereb10.cereb35.exc))
cereb35.cereb35.exc <- cosine.values[13:18,13:18]
cereb35.cereb35.exc <- c(as.matrix(pull_lower_triangle(cereb35.cereb35.exc, diagonal = TRUE))[,2:7])
cereb35.cereb35.exc.nms <- rep("cereb35.cereb35.exc", length(cereb35.cereb35.exc))
pfc10.pfc10.inh <- cosine.values[10:12,10:12]
pfc10.pfc10.inh <- c(as.matrix(pull_lower_triangle(pfc10.pfc10.inh, diagonal = TRUE))[,2:4])
pfc10.pfc10.inh.nms <- rep("pfc10.pfc10.inh", length(pfc10.pfc10.inh))
pfc10.pfc35.inh <- c(cosine.values[10:12,30:34])
pfc10.pfc35.inh.nms <- rep("pfc10.pfc35.inh", length(pfc10.pfc35.inh))
pfc10.cereb10.inh <- c(cosine.values[10:12,7:9])
pfc10.cereb10.inh.nms <- rep("pfc10.cereb10.inh", length(pfc10.cereb10.inh))
pfc10.cereb35.inh <- c(cosine.values[10:12,24:29])
pfc10.cereb35.inh.nms <- rep("pfc10.cereb35.inh", length(pfc10.cereb35.inh))
pfc35.pfc35.inh <-cosine.values[30:34,30:34]
pfc35.pfc35.inh <- c(as.matrix(pull_lower_triangle(pfc35.pfc35.inh, diagonal = TRUE))[,2:6])
pfc35.pfc35.inh.nms <- rep("pfc35.pfc35.inh", length(pfc35.pfc35.inh))
pfc35.cereb10.inh <- c(cosine.values[30:34,7:9])
pfc35.cereb10.inh.nms <- rep("pfc35.cereb10.inh", length(pfc35.cereb10.inh))
pfc35.cereb35.inh <- c(cosine.values[30:34, 24:29])
pfc35.cereb35.inh.nms <- rep("pfc35.cereb35.inh", length(pfc35.cereb35.inh))
cereb10.cereb10.inh <- cosine.values[7:9,7:9]
cereb10.cereb10.inh <- c(as.matrix(pull_lower_triangle(cereb10.cereb10.inh, diagonal = TRUE))[,2:4])
cereb10.cereb10.inh.nms <- rep("cereb10.cereb10.inh", length(cereb10.cereb10.inh))
cereb10.cereb35.inh <- c(cosine.values[7:9,24:29])
cereb10.cereb35.inh.nms <- rep("cereb10.cereb35.inh", length(cereb10.cereb35.inh))
cereb35.cereb35.inh <- cosine.values[24:29,24:29]
cereb35.cereb35.inh <- c(as.matrix(pull_lower_triangle(cereb35.cereb35.inh, diagonal = TRUE))[,2:7])
cereb35.cereb35.inh.nms <- rep("cereb35.cereb35.inh", length(cereb35.cereb35.inh))

pfc.10.exc <- as.numeric(c(pfc10.pfc10.exc, pfc10.pfc35.exc, pfc10.cereb10.exc, pfc10.cereb35.exc))
pfc.35.exc <- as.numeric(c(pfc35.pfc35.exc, pfc10.pfc35.exc, pfc35.cereb10.exc, pfc35.cereb35.exc))
cereb.10.exc <- as.numeric(c(cereb10.cereb10.exc, cereb10.cereb35.exc, pfc10.cereb10.exc, pfc35.cereb10.exc))
cereb.35.exc <- as.numeric(c(cereb35.cereb35.exc, cereb10.cereb35.exc, pfc10.cereb35.exc, pfc35.cereb35.exc))
pfc.10.exc.nms <- c(pfc10.pfc10.exc.nms, pfc10.pfc35.exc.nms, pfc10.cereb10.exc.nms, pfc10.cereb35.exc.nms)
pfc.35.exc.nms <- c(pfc35.pfc35.exc.nms, pfc10.pfc35.exc.nms, pfc35.cereb10.exc.nms, pfc35.cereb35.exc.nms)
cereb.10.exc.nms <- c(cereb10.cereb10.exc.nms, cereb10.cereb35.exc.nms, pfc10.cereb10.exc.nms, pfc35.cereb10.exc.nms)
cereb.35.exc.nms <- c(cereb35.cereb35.exc.nms, cereb10.cereb35.exc.nms, pfc10.cereb35.exc.nms, pfc35.cereb35.exc.nms)

pfc.10.exc <- data.frame(cosine.similarity = pfc.10.exc, comparisons = pfc.10.exc.nms)
pfc.35.exc <- data.frame(cosine.similarity = pfc.35.exc, comparisons = pfc.35.exc.nms)
cereb.10.exc <- data.frame(cosine.similarity = cereb.10.exc, comparisons = cereb.10.exc.nms)
cereb.35.exc <- data.frame(cosine.similarity = cereb.35.exc, comparisons = cereb.35.exc.nms)

library(ggpubr)
my_comparisons <- rev(list(c("pfc10.pfc10.exc","pfc10.pfc35.exc"),c("pfc10.pfc10.exc","pfc10.cereb10.exc"),c("pfc10.pfc10.exc","pfc10.cereb35.exc")))
ggbarplot(pfc.10.exc, x = "comparisons", y = "cosine.similarity",
          add = c("mean_se", "point"),
          color = "comparisons", fill = "comparisons", alpha = 0.5, main = "PFC 10 exc - Cosine Similarity of PCA") + stat_compare_means(comparison = my_comparisons, method = "t.test")
write.csv(pfc.10.exc, file = "/home/hailee/10X/Final Soc Analysis/PCA comparisons/bl normalized/Cosine Similarity plots/pfc10exc_cos_similarity.csv")

my_comparisons <- rev(list(c("pfc35.pfc35.exc","pfc10.pfc35.exc"),c("pfc35.pfc35.exc","pfc35.cereb10.exc"),c("pfc35.pfc35.exc","pfc35.cereb35.exc")))
ggbarplot(pfc.35.exc, x = "comparisons", y = "cosine.similarity",
          add = c("mean_se", "point"),
          color = "comparisons", fill = "comparisons", alpha = 0.5, main = "PFC 35 exc - Cosine Similarity of PCA") + stat_compare_means(comparison = my_comparisons, method = "t.test")
write.csv(pfc.35.exc, file = "/home/hailee/10X/Final Soc Analysis/PCA comparisons/bl normalized/Cosine Similarity plots/pfc35exc_cos_similarity.csv")

my_comparisons <- rev(list(c("cereb35.cereb35.exc","cereb10.cereb35.exc"),c("cereb35.cereb35.exc","pfc10.cereb35.exc"),c("cereb35.cereb35.exc","pfc35.cereb35.exc")))
ggbarplot(cereb.35.exc, x = "comparisons", y = "cosine.similarity",
          add = c("mean_se", "point"),
          color = "comparisons", fill = "comparisons", alpha = 0.5, main = "Cereb 35 exc - Cosine Similarity of PCA") + stat_compare_means(comparison = my_comparisons, method = "t.test")
write.csv(cereb.35.exc, file = "/home/hailee/10X/Final Soc Analysis/PCA comparisons/bl normalized/Cosine Similarity plots/cereb35exc_cos_similarity.csv")

my_comparisons <- rev(list(c("cereb10.cereb10.exc","cereb10.cereb35.exc"),c("cereb10.cereb10.exc","pfc10.cereb10.exc"),c("cereb10.cereb10.exc","pfc35.cereb10.exc")))
ggbarplot(cereb.10.exc, x = "comparisons", y = "cosine.similarity",
          add = c("mean_se", "point"),
          color = "comparisons", fill = "comparisons", alpha = 0.5, main = "Cereb 10 exc - Cosine Similarity of PCA") + stat_compare_means(comparison = my_comparisons, method = "t.test")
write.csv(cereb.10.exc, file = "/home/hailee/10X/Final Soc Analysis/PCA comparisons/bl normalized/Cosine Similarity plots/cereb10exc_cos_similarity.csv")

pfc.10.inh <- as.numeric(c(pfc10.pfc10.inh, pfc10.pfc35.inh, pfc10.cereb10.inh, pfc10.cereb35.inh))
pfc.35.inh <- as.numeric(c(pfc35.pfc35.inh, pfc10.pfc35.inh, pfc35.cereb10.inh, pfc35.cereb35.inh))
cereb.10.inh <- as.numeric(c(cereb10.cereb10.inh, cereb10.cereb35.inh, pfc10.cereb10.inh, pfc35.cereb10.inh))
cereb.35.inh <- as.numeric(c(cereb35.cereb35.inh, cereb10.cereb35.inh, pfc10.cereb35.inh, pfc35.cereb35.inh))
pfc.10.inh.nms <- c(pfc10.pfc10.inh.nms, pfc10.pfc35.inh.nms, pfc10.cereb10.inh.nms, pfc10.cereb35.inh.nms)
pfc.35.inh.nms <- c(pfc35.pfc35.inh.nms, pfc10.pfc35.inh.nms, pfc35.cereb10.inh.nms, pfc35.cereb35.inh.nms)
cereb.10.inh.nms <- c(cereb10.cereb10.inh.nms, cereb10.cereb35.inh.nms, pfc10.cereb10.inh.nms, pfc35.cereb10.inh.nms)
cereb.35.inh.nms <- c(cereb35.cereb35.inh.nms, cereb10.cereb35.inh.nms, pfc10.cereb35.inh.nms, pfc35.cereb35.inh.nms)

pfc.10.inh <- data.frame(cosine.similarity = pfc.10.inh, comparisons = pfc.10.inh.nms)
pfc.35.inh <- data.frame(cosine.similarity = pfc.35.inh, comparisons = pfc.35.inh.nms)
cereb.10.inh <- data.frame(cosine.similarity = cereb.10.inh, comparisons = cereb.10.inh.nms)
cereb.35.inh <- data.frame(cosine.similarity = cereb.35.inh, comparisons = cereb.35.inh.nms)

library(ggpubr)
my_comparisons <- rev(list(c("pfc10.pfc10.inh","pfc10.pfc35.inh"),c("pfc10.pfc10.inh","pfc10.cereb10.inh"),c("pfc10.pfc10.inh","pfc10.cereb35.inh")))
ggbarplot(pfc.10.inh, x = "comparisons", y = "cosine.similarity",
          add = c("mean_se", "point"),
          color = "comparisons", fill = "comparisons", alpha = 0.5, main = "PFC 10 inh - Cosine Similarity of PCA") + stat_compare_means(comparison = my_comparisons, method = "t.test")
write.csv(pfc.10.inh, file = "/home/hailee/10X/Final Soc Analysis/PCA comparisons/bl normalized/Cosine Similarity plots/pfc10inh_cos_similarity.csv")

my_comparisons <- rev(list(c("pfc35.pfc35.inh","pfc10.pfc35.inh"),c("pfc35.pfc35.inh","pfc35.cereb10.inh"),c("pfc35.pfc35.inh","pfc35.cereb35.inh")))
ggbarplot(pfc.35.inh, x = "comparisons", y = "cosine.similarity",
          add = c("mean_se", "point"),
          color = "comparisons", fill = "comparisons", alpha = 0.5, main = "PFC 35 inh - Cosine Similarity of PCA") + stat_compare_means(comparison = my_comparisons, method = "t.test")
write.csv(pfc.35.inh, file = "/home/hailee/10X/Final Soc Analysis/PCA comparisons/bl normalized/Cosine Similarity plots/pfc35inh_cos_similarity.csv")

my_comparisons <- rev(list(c("cereb35.cereb35.inh","cereb10.cereb35.inh"),c("cereb35.cereb35.inh","pfc10.cereb35.inh"),c("cereb35.cereb35.inh","pfc35.cereb35.inh")))
ggbarplot(cereb.35.inh, x = "comparisons", y = "cosine.similarity",
          add = c("mean_se", "point"),
          color = "comparisons", fill = "comparisons", alpha = 0.5, main = "Cereb 35 inh - Cosine Similarity of PCA")+ stat_compare_means(comparison = my_comparisons, method = "t.test")
write.csv(cereb.35.inh, file = "/home/hailee/10X/Final Soc Analysis/PCA comparisons/bl normalized/Cosine Similarity plots/cereb35inh_cos_similarity.csv")

my_comparisons <- rev(list(c("cereb10.cereb10.inh","cereb10.cereb35.inh"),c("cereb10.cereb10.inh","pfc10.cereb10.inh"),c("cereb10.cereb10.inh","pfc35.cereb10.inh")))
ggbarplot(cereb.10.inh, x = "comparisons", y = "cosine.similarity",
          add = c("mean_se", "point"),
          color = "comparisons", fill = "comparisons", alpha = 0.5, main = "Cereb 10 inh - Cosine Similarity of PCA") + stat_compare_means(comparison = my_comparisons, method = "t.test")
write.csv(cereb.10.inh, file = "/home/hailee/10X/Final Soc Analysis/PCA comparisons/bl normalized/Cosine Similarity plots/cereb10inh_cos_similarity.csv")