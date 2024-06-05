## Figure 1 - PFC cluster Metrics
## Hailee Walker - 07142023
# Runs clustering using Seurat
#uses marker genes and Allen brain through Enrichr to allow for identification and renaming of PFC cell types
#Code to create all panels in figure 1 except A
#Adapted from Seurat vignette protocols

##Import Files
rm(list=ls())
graphics.off()

##Open libraries and source any functions you'll need
library(Seurat)
library(SeuratDisk)
library(sctransform)
library(ggplot2)
library(ggsignif)
library(patchwork)
library(tidyverse)
library(devtools)
library(data.table)
library(magrittr)

#Edit to make a csv that matches your files
samples <- read_csv("id,name
20042X1,Soc WT 1
20042X2,Not WT 1
20042X3,Soc WT 2
20042X4,Not WT 2
20042X5,Soc WT 3
20042X6,Not WT 3
19957X3,Not WT 4
20778X1,ShS WT 1
20778X3,ShS WT 2
20778X5,ShS WT 3
20251X1,Soc WT 4
20251X4,Soc KO 1
20375X1,Not KO 1
20375X6,Not WT 5
20521X1,Soc WT 5
20521X5,Soc KO 2
20523X1,Not WT 6
20523X5,Not KO 2")
n <- nrow(samples)
s1 <- vector("list", n)
for(i in 1:6){
  message("Loading ", samples$id[i])
  x <- Read10X(paste0(data.dir ="/home/hailee/10X/20042/", "MEX/", samples$id[i]))
  s1[[i]] <- CreateSeuratObject(counts = x, project = samples$name[i])
}
for(i in 7){
  message("Loading ", samples$id[i])
  x <- Read10X(paste0(data.dir ="/home/hailee/10X/19957/", samples$id[i],"/outs/filtered_feature_bc_matrix/"))
  s1[[i]] <- CreateSeuratObject(counts = x, project = samples$name[i])
}
for(i in 8:10){
  message("Loading ", samples$id[i])
  x <- Read10X(paste0(data.dir ="/home/hailee/10X/20778R/", samples$id[i],"/outs/filtered_feature_bc_matrix/"))
  s1[[i]] <- CreateSeuratObject(counts = x, project = samples$name[i])
}
for(i in 11:12){
  message("Loading ", samples$id[i])
  x <- Read10X(paste0(data.dir ="/home/hailee/10X/Shank 3/20251/", samples$id[i],"/outs/filtered_feature_bc_matrix/"))
  s1[[i]] <- CreateSeuratObject(counts = x, project = samples$name[i])
}
for(i in 13:14){
  message("Loading ", samples$id[i])
  x <- Read10X(paste0(data.dir ="/home/hailee/10X/Shank 3/20375R/", samples$id[i],"/outs/filtered_feature_bc_matrix/"))
  s1[[i]] <- CreateSeuratObject(counts = x, project = samples$name[i])
}
for(i in 15:16){
  message("Loading ", samples$id[i])
  x <- Read10X(paste0(data.dir ="/home/hailee/10X/Shank 3/20521R/", samples$id[i],"/outs/filtered_feature_bc_matrix/"))
  s1[[i]] <- CreateSeuratObject(counts = x, project = samples$name[i])
}
for(i in 17:n){
  message("Loading ", samples$id[i])
  x <- Read10X(paste0(data.dir ="/home/hailee/10X/Shank 3/20523R/", samples$id[i],"/outs/filtered_feature_bc_matrix/"))
  s1[[i]] <- CreateSeuratObject(counts = x, project = samples$name[i])
}

# Create a merged Seurat object.
CellFile <- merge(x = s1[[1]], y = s1[2:n], add.cell.ids = samples$name, project="All PFC")

# Format metadata
slotNames(CellFile)
CellFile_meta <- function(CellFile){
  rownames_to_column( CellFile@meta.data, "barcode") %>%
    mutate(barcode=gsub(".*_", "", barcode)) %>%
    dplyr::rename(sample = "orig.ident") %>%
    as_tibble()
}
m1 <- CellFile_meta(CellFile)

##Normalize data
CellFile <- PercentageFeatureSet(CellFile, pattern = "^mt-", col.name = "MT")
CellFile <-  subset(CellFile, MT < 5 & nFeature_RNA > 800 & nFeature_RNA < 6000)
CellFile <- SCTransform(CellFile, vars.to.regress = "MT", verbose = FALSE)

# run PCA and UMAP, plot clusters
CellFile <- RunPCA(CellFile, verbose = FALSE)
ElbowPlot(CellFile)
CellFile <- FindNeighbors(CellFile, dims = 1:30, verbose = FALSE)
CellFile <- FindClusters(CellFile, verbose = TRUE)
CellFile <- RunUMAP(CellFile, dims = 1:30, verbose = FALSE)
#CellFile <- RunTSNE(CellFile, dims =1:30, verbose = FALSE)
DimPlot(CellFile, label = TRUE) + NoLegend()

F##Find and remove any clusters that doen't exist in at least % of samples
table(CellFile@meta.data$SCT_snn_res.0.8, CellFile@meta.data$orig.ident)
rowSums(table(CellFile@meta.data$SCT_snn_res.0.8, CellFile@meta.data$orig.ident))
CellFile <- subset(CellFile, idents = c(29,30,37,42,43), invert = TRUE)

#create orig.id, genotype and status columns to use in split.by or group.by
CellFile$status <- substr(CellFile$orig.ident, 1,3)
CellFile$orig.ident <- factor(CellFile$orig.ident, levels=unique(m1$sample))
CellFile$genotype <- substr(CellFile$orig.ident, 5,6)
CellFile$condition <- substr(CellFile$orig.ident, 1,6)

#Identify and Rename Clusters (Allen brain, Bhatacharjee, A. 2019) (genes listed in S1)
cell.type.markers <- c("Nrn1","Gad2", "Slc17a7","Cux2", "Calb1", "Etv1", "Syt6", "Oprk1", "Car3", "Pvalb", "Sst", "Chodl", "Vip", "Chat", "Lamp5", "Meis2", "Pbx3", "Gfap", "Aqp4", "Aspa", "Mbp", "Bmp4", "Pdgfra", "C1qa", "Tmem119", "Ogn", "Flt1", "Foxj1","Cspg4")
DotPlot(CellFile, features = cell.type.markers, dot.scale = 8) +
  RotatedAxis()

#Find Markers
res <- FindAllMarkers(object = CellFile, only.pos = TRUE, logfc.threshold = 0.5)
library(hciRdata)
clust1 <- dplyr::select(res, 6,7,2,3,4,1,5) %>%
  left_join( unique(mouse98[, c(2,8)]), by=c(gene="gene_name")) %>%
  as_tibble() %>% filter(p_val_adj < 0.05) %>% arrange(cluster, p_val_adj, desc(abs(pct.1-pct.2)))
clust1
x <- AverageExpression(object = CellFile, assay="SCT", slot="data")$SCT
c1 <- gather(rownames_to_column(as.data.frame(x), "gene"), "cluster", "mean.1", -gene) %>%
  as_tibble()  %>% mutate(mean.1 = log2(mean.1 + 1))
c2 <- inner_join(clust1, c1, by=c("cluster", "gene"))
clust1 <- mutate(c2, mean.2 = mean.1 - avg_log2FC) %>% dplyr::select(1:2, 9:10, 3:8)

source("/home/hailee/10X/Enrichr/enrichr_clusters.R")
source("/home/hailee/10X/Enrichr/plot_top_enrichr.R")
#enriched pathways and cell types
library(enrichR)
setEnrichrSite("Enrichr")
# dbs <- as_tibble(listEnrichrDbs())
celldbs <- c("Allen_Brain_Atlas_10x_scRNA_2021")
e1 <- enrichr_clusters(clust1, celldbs)
data.frame(n =sapply(e1, nrow))
#Plot top from any dbs (change the number in spot 2)
plot_top_enrichr(e1, 1, top=3, min_p=7)

### Change to correct identities - Make sure each cluster has the correct cell type named
CellFile <- RenameIdents(CellFile, `0` = "L2/3 Exc", `1` = "L6 Exc Syt6", `2` = "Oligo", `3` = "L2/3 Exc", `4`= "Astro", `5` = "L6_Exc_Oprk1", 
                         `6` = "L5 Exc", `7`="L5 Exc",`8` = "L5 Exc", `9` = "PV",`10` = "L2/3/5 Exc", `11` = "Astro", `12` = "Mesi2,Pbx3-1", 
                         `13` = "Microglia", `14` = "L5 NP CTX", `15` = "Som", `16` = "OPC", `17` = "VIP", `18` = "Meis2,Pbx3-2", `19` = "Lamp5", 
                         `20` = "VLMC", `21` = "Astro", `22` = "Oligo", `23`= "L6b CTX", `24` = "L2/3 Exc", `25` = "Som_Chodl", 
                         `26` = "VLMC", `27` = "L2/3/5 Exc", `28` = "Microglia", `31` = "Oligo", `32` = "L6 Car3", 
                         `33` = "Endothelial", `34` = "Pericyte", `35` = "Oligo", `36` = "L6 Exc Syt6", `38` = "L6 Exc Syt6", 
                         `39` = "VLMC", `40` = "PV", `41` = "Oligo")
CellFile$celltype <- Idents(CellFile)
Idents(CellFile) <- "celltype"

#Reorganize cell types
levels(CellFile) <-c("L2/3 Exc", "L2/3/5 Exc", "L5 Exc", "L5 NP CTX", "L6 Exc Syt6", "L6b CTX", "L6 Exc Oprk1", "L6 Car3", "PV", "Som", "Som Chodl", 
                     "VIP", "Lamp5", "Meis2,Pbx3-1", "Meis2,Pbx3-2", "Astro", "Oligo", "OPC", "Microglia", "Endothelial", "VLMC", "Pericyte")
save(CellFile, file = "/home/hailee/10X/Final Soc Analysis/Comb_PFC.Rdata")

### WT Analysis
Sub <- subset(CellFile, subset = genotype == "WT")

## Figure 1E ## ___________________________________________________________________________
#Create Feature plot of some marker genes
marks <- c("Slc17a7", "Cux2", "Etv1", "Gad2", "Vip", "Pvalb", "Sst", "Mbp", "Gfap", "Oprk1")
FeaturePlot(Sub, marks, cols = c("grey", "darkred"), keep.scale = NULL)
p <- GetAssayData(object=Sub, slot = "data")[marks[9],]
pos_ids = names(which(p>1.2)) #Change this value for each gene to about the mean expression
pos_cells = subset(Sub,cells=pos_ids)
DimPlot(pos_cells, cols = rep("blue", length(levels(pos_cells))))+NoLegend() + xlim(-18,18) + ylim(-18, 18)
DimPlot(Sub, cols = rep("grey", length(levels(Sub))))+NoLegend() + xlim(-18,18) + ylim(-18, 18)

## Figure 1B ## ___________________________________________________________________________
##Replot UMAP with new names and colored by cell type
color <- c("brown2", "tomato", "darkorange1", "darkorange2", "darkgoldenrod2", "darkgoldenrod3", "yellow3", "tan2", "deepskyblue", "aquamarine4", 
           "cyan4", "aquamarine3", "darkturquoise", "cornflowerblue", "lightskyblue", "blueviolet", "mediumpurple3", "mediumpurple1", "mediumorchid1","magenta1", 
           "deeppink", "hotpink2")
DimPlot(Sub, label = TRUE, cols = color)

## Figure 1C ## ____________________________________________________________________________
##Make Violin plot of marker genes
cell.type.markers <- c("Slc17a7","Gad2", "Cux2", "Calb1", "Etv1", "Syt6", "Oprk1", "Car3", "Pvalb", "Sst", "Chodl", "Vip", "Lamp5", 
                       "Pbx3", "Aqp4", "Mbp", "Pdgfra", "C1qa", "Flt1", "Ogn", "Cspg4")
VlnPlot(Sub, features = cell.type.markers, stack = TRUE, flip = TRUE, fill.by = "ident", cols = color) + NoLegend()

## Figure 1D ## ____________________________________________________________________________
##Make Bargraph of proportion of cells by type
md <- Sub@meta.data %>% as.data.table
md$celltype <- factor(md$celltype, levels = rev(levels(Sub)))

mc <- md[, mean(nFeature_RNA), by = c("celltype")]
md <- md[, .N, by = c("celltype")]
md <- md[order(md$celltype)]
mc <- mc[order(mc$celltype)]
WTsum <- sum(md$N)
md$totals <- rep(WTsum, nrow(md))
md$val <- md$N/md$totals

barplot <- barplot(as.numeric(md$val), names=md$celltype, horiz = TRUE, cex.names = 0.7, las = 1, col = rev(color), xlim = c(0,0.23)) 
for(b in 1:nrow(md)){
  l <- paste(substr(mc$V1[b],1,6))
  text(0.20,barplot[b], l, cex = 0.7)
  
}

