##Figure 4 - IEGs by cell type
#Hailee Walker - 04082024
#Reruns UMAP to cluster PFC neurons based only on their IEG expression
#Creates heatmap of cell types most likely to belong to each IEG program cluster
#allows you to determine the marker genes which allow for clustering of cells into these clusters

rm(list =ls())
graphics.off()

library(Seurat)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(foreach)
library(doParallel)
library(pheatmap)
source("/home/hailee/10X/Final Soc Analysis/R functions/Early_v_Late.R")
source("/home/hailee/10X/Final Soc Analysis/R functions/mIEG_dotplot.R")
source("/home/hailee/10X/Final Soc Analysis/R functions/IEG pos distributions.R")

#Import Seurat object and subset to just WT and neurons
load("/home/hailee/10X/Final Soc Analysis/PFC_All_analysis/Comb_PFC.Rdata") #Import Seurat Object
Sub <- subset(CellFile, subset = genotype == "WT") #Subset WT
Neurons <- subset(Sub, idents = c("L2/3_Exc", "L2/3/5_Exc", "L5_Exc", "L5_NP_CTX", "L6_Exc_Syt6", "L6b_CTX", "L6_Exc_Oprk1", "L6 Car3", "PV", "Som", "Som_Chodl", 
                                  "VIP", "Lamp5", "Meis2,Pbx3-1", "Meis2,Pbx3-2"))

genes <- c("Apold1", "Arc", "Bhlhe40", "Btg2", "Ccn1", "Ccn2", "Cebpb", "Dusp1", "Dusp2", "Dusp6", "Egr1", "Egr2", "Egr3", "Egr4", "F3", "Fos", "Fosb", 
           "Fosl2", "Gadd45b", "Gadd45g", "Gem", "Hes1", "Ier2", "Ier3", "Ier5", "Ifit2", "Ifit3", "Jun", "Junb", "Klf10", "Maff", "Nfkbia", "Npas4", 
           "Nr4a1", "Nr4a2", "Nr4a3", "Nrn1", "Per1", "Ptgs2", "Sik1", "Trib1", "Zfp36l1", "Zfp36l2") #PFC genes

IEG_Sub_s <- subset(Neurons, features = genes) #Subset seurat object of neurons to only include IEGs
Not <- subset(IEG_Sub_s, subset = status == "Not")

#set up to run parallel loop
totalCores = detectCores()
cluster <- makeCluster(totalCores[1]-10)
registerDoParallel(cluster)

#Get counts matrix of IEGs by celltype for not social
stim.list <- SplitObject(Not, "celltype")
y1 <- stim.list[[1]]
y <- y1@assays$RNA@counts
mat_not <- matrix(ncol = length(stim.list), nrow = nrow(y) )
for(i in 1:length(stim.list)){
  print(i)
  y1 <- stim.list[[i]]
  y <- y1@assays$RNA@counts
  mat2 <- foreach(c = 1:nrow(y), .combine = rbind) %dopar% {
    gene = sum(y[c,])
  }
  mat_not[,i] = mat2
}
nmc <- names(stim.list)
nmr <- rownames(y)
colnames(mat_not) <- nmc
rownames(mat_not) <- nmr
rm(y, y1)

#Get counts matrix of IEGs by celltype for 10 min
IEG_Sub_soc <- subset(IEG_Sub_s, subset = status == "Soc")
stim.list <- SplitObject(IEG_Sub_soc, "celltype")
y1 <- stim.list[[1]]
y <- y1@assays$RNA@counts
mat_soc <- matrix(ncol = length(stim.list), nrow = nrow(y) )
for(i in 1:length(stim.list)){
  print(i)
  y1 <- stim.list[[i]]
  y <- y1@assays$RNA@counts
  mat2 <- foreach(c = 1:nrow(y), .combine = rbind) %dopar% {
    gene = sum(y[c,])
  }
  mat_soc[,i] = mat2
}
nmc <- names(stim.list)
nmr <- rownames(y)
colnames(mat_soc) <- nmc
rownames(mat_soc) <- nmr
rm(y, y1)

#Get counts matrix of IEGs by celltype for 10 min
IEG_Sub_shs <- subset(IEG_Sub_s, subset = status == "ShS")
stim.list <- SplitObject(IEG_Sub_shs, "celltype")
y1 <- stim.list[[1]]
y <- y1@assays$RNA@counts
mat_shs <- matrix(ncol = length(stim.list), nrow = nrow(y) )
for(i in 1:length(stim.list)){
  print(i)
  y1 <- stim.list[[i]]
  y <- y1@assays$RNA@counts
  mat2 <- foreach(c = 1:nrow(y), .combine = rbind) %dopar% {
    gene = sum(y[c,])
  }
  mat_shs[,i] = mat2
}
stopCluster(cluster)
nmc <- names(stim.list)
nmr <- rownames(y)
colnames(mat_shs) <- nmc
rownames(mat_shs) <- nmr
rm(y, y1)

#divide social groups by not social
mat_shs_not <- mat_shs/mat_not
mat_soc_not <- mat_soc/mat_not
mat_soc_not[mat_soc_not == Inf] <- NaN
mat_shs_not[mat_shs_not == Inf] <- NaN

#Supp. Table 7
#Get table of top 5 upregulated IEGs per cell type for 10 min condition
IEG_tab <- matrix(ncol = ncol(mat_shs_not), nrow = 5)
for(i in 1:ncol(mat_shs_not)){
  top <- rownames(mat_shs_not[order(mat_shs_not[, i], decreasing=TRUE)[1:5], i, drop = FALSE])
  IEG_tab[,i] <- top
}
colnames(IEG_tab) <- colnames(mat_shs_not)
IEG_tab <- as.data.frame(IEG_tab)
IEG_tab<- IEG_tab[,c("L2/3_Exc", "L2/3/5_Exc", "L5_Exc", "L5_NP_CTX", "L6_Exc_Syt6", "L6b_CTX", "L6_Exc_Oprk1", "L6 Car3", "PV", "Som", "Som_Chodl", 
           "VIP", "Lamp5", "Meis2,Pbx3-1", "Meis2,Pbx3-2")]
write.csv(IEG_tab, file = "~/10X/Final Soc Analysis/Originals by Figure/IEG_top_shs.csv")

#Get table of top 5 upregulated IEGs per cell type for 35 min condition
IEG_tab <- matrix(ncol = ncol(mat_soc_not), nrow = 5)
for(i in 1:ncol(mat_soc_not)){
  top <- rownames(mat_soc_not[order(mat_soc_not[, i], decreasing=TRUE)[1:5], i, drop = FALSE])
  IEG_tab[,i] <- top
}
colnames(IEG_tab) <- colnames(mat_soc_not)
IEG_tab <- as.data.frame(IEG_tab)
IEG_tab<- IEG_tab[,c("L2/3_Exc", "L2/3/5_Exc", "L5_Exc", "L5_NP_CTX", "L6_Exc_Syt6", "L6b_CTX", "L6_Exc_Oprk1", "L6 Car3", "PV", "Som", "Som_Chodl", 
                     "VIP", "Lamp5", "Meis2,Pbx3-1", "Meis2,Pbx3-2")]
write.csv(IEG_tab, file = "~/10X/Final Soc Analysis/Originals by Figure/IEG_top_soc.csv")

genes <- (as.list(na.omit(read_csv("/home/hailee/10X/IEG list (Hong,W. 2017).csv", col_names = FALSE)[,1])))[[1]]
IEG_Sub_s <- subset(Neurons, features = genes) #Subset seurat object of neurons to only include IEGs
Not <- subset(IEG_Sub_s, subset = status == "Not")
mat <- matrix(ncol=length(genes), nrow= ncol(IEG_Sub_s))
##Find cells that express 4+ IEGs above the 95% of expression in not social for that IEG
for(i in 1:length(genes)){
  p <- GetAssayData(object=Not, slot = "data")[genes[i], ]
  tmp <- GetAssayData(object=IEG_Sub_s, slot = "data")[genes[i], ]
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
v1 <- subset(IEG_Sub_s, cells = pos_ids)

color <- c("brown2", "tomato", "darkorange1", "darkorange2", "darkgoldenrod2", "darkgoldenrod3", "yellow3", "tan2", "deepskyblue", 
           "aquamarine4", "cyan4", "aquamarine3", "darkturquoise", "cornflowerblue", "lightskyblue")
#Rerun PCA and nearest neighbor to create new UMAP based only on IEGs
IEG_Sub <- RunPCA(v1, verbose = FALSE, npcs = 10)
ElbowPlot(IEG_Sub)
IEG_Sub <- FindNeighbors(IEG_Sub, dims = 1:10, verbose = FALSE)
IEG_Sub <- FindClusters(IEG_Sub, verbose = TRUE, resolution = 0.3)
IEG_Sub <- RunUMAP(IEG_Sub, dims = 1:10, verbose = FALSE)

#UMAP orig clusters, Figure 4 A
IEG_Sub <- subset(IEG_Sub, subset = status == "Not", invert = TRUE)
#Figure 4A
DimPlot(IEG_Sub, label = FALSE)
#Figure 4 C
FeaturePlot(IEG_Sub, features = c("Fos", "Nr4a3","Arc", "Npas4", "Sgk1", "Junb"), cols = c("lightgrey", "#660000"))

#Supplemental Table 2
##Find the percentage of each cluster made up by each cell type
tb <- as.data.frame(table(IEG_Sub$seurat_clusters, IEG_Sub$celltype))
tb <- dcast(tb, Var1~Var2)
rownames(tb) <- tb[,1]
tb <- tb[,2:ncol(tb)]
tb_prop_cluster<-t(apply(tb,1, function(x) x/sum(x)))*100
tb_prop_CT<-t(t(tb)/colSums(tb))*100
write.csv(tb_prop_CT, file = "/home/hailee/10X/Final Soc Analysis/IEG only UMAPs/prop_of_CT.csv")
write.csv(tb_prop_cluster, file = "/home/hailee/10X/Final Soc Analysis/IEG only UMAPs/prop_of_cluster.csv")

#Figure 4 E
tb_prop_CT[is.nan(tb_prop_CT)] <- 0
tb_prop_CT <- tb_prop_CT[,colSums(tb_prop_CT) >0]
tb_hm <- tb_prop_CT[,c("L2/3_Exc", "L2/3/5_Exc", "L5_Exc", "L5_NP_CTX", "L6_Exc_Syt6", "L6b_CTX", "L6_Exc_Oprk1", "L6 Car3", "PV", "Som", "Som_Chodl", 
                            "VIP", "Lamp5", "Meis2,Pbx3-1", "Meis2,Pbx3-2")]
col_fun = colorRamp2(c(0, 10,20, 30, 40), c("white", "#FF9999", "#CC6666", "#993333", "#660000"))
Heatmap(tb_hm, cluster_rows = F,name = "Prop of CT", cluster_columns = F, col = col_fun)
#log scale heatmap
tb_hm[tb_hm == 0] <- 0.00001
log_tb_hm <- log10(tb_hm)
col_fun = colorRamp2(c(-0.25, 0,0.699, 1,1.176, 1.301, 1.398), c( "white", brewer.pal(n = 6, name = "OrRd")))
Heatmap(log_tb_hm, cluster_rows = F,name = "log(CT %age)", cluster_columns = F, col = col_fun)

#Suplemental Table 1
# Make Table of top IEGs by cluster
lev <- as.character(levels(IEG_Sub))

top.genes <- list()
for(i in 1:length(lev)){
  clust <- FindMarkers(IEG_Sub, ident.1 = lev[i], min.pct = 0.25)
  clust <- clust[order(clust$avg_log2FC, decreasing = TRUE), ]
  clust<- clust[,c("avg_log2FC", "p_val_adj")]
  top.15 <- head(clust, n=15)
  top.genes[[i]] <- clust
}
names(top.genes) <- lev

#UMAP by cell type - Figure 4B
Idents(IEG_Sub) <- "celltype"
levels(IEG_Sub) <- c("L2/3_Exc","L2/3/5_Exc","L5_Exc","L5_NP_CTX","L6_Exc_Syt6","L6b_CTX","L6_Exc_Oprk1","L6 Car3",
             "PV", "Som", "Som_Chodl", "VIP", "Lamp5", "Meis2,Pbx3-1", "Meis2,Pbx3-2")
DimPlot(IEG_Sub, label = FALSE, cols = color)

#Figure 4D. Change out correct spot to match cell type you want to look at and make it the correct color for that cell type
#Plot individual cell types
color <- c("lightgrey", "lightgrey", "lightgrey", "darkorange2", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", 
           "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey") #L5 NP CTX
DimPlot(IEG_Sub, label = FALSE, cols = color) + NoLegend()

color <- c("brown2", "tomato", "darkorange1", "darkorange2", "darkgoldenrod2", "darkgoldenrod3", "yellow3", "tan2", "deepskyblue", 
           "aquamarine4", "cyan4", "aquamarine3", "darkturquoise", "cornflowerblue", "lightskyblue")
