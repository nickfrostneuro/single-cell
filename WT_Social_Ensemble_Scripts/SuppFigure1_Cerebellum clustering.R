## Supplemental Figure  1 - Cerebellum Clustering
## Hailee Walker 07142023
# Runs clustering using Seurat
#uses marker genes to allow for identification and renaming of Cerebellum cell types
#Code to create all panels in supp. figure 1

#Clear Environment
rm(list=ls())
graphics.off()

#Open necessary libraries and functions
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
20053X1,Soc WT 1
20053X2,Not WT 1
20053X3,Soc WT 2
20053X4,Not WT 2
20053X5,Soc WT 3
20053X6,Not WT 3
19957X2,Soc WT 4
20778X2,ShS WT 1
20778X4,ShS WT 2
20778X6,ShS WT 3
20251X3,Soc WT 5
20251X6,Soc KO 1
20375X5,Not KO 1
20375X10,Not WT 4
20521X4,Soc WT 6
20521X8,Soc KO 2
20523X4,Not WT 5
20523X8,Not KO 2")
n <- nrow(samples)
s1 <- vector("list", n)
for(i in 1:6){
  message("Loading ", samples$id[i])
  x <- Read10X(paste0(data.dir ="/home/hailee/10X/20053/MEX/", samples$id[i]))
  s1[[i]] <- CreateSeuratObject(counts = x, project = samples$name[i])
}
for(i in 7){
  message("Loading ", samples$id[i])
  x <- Read10X(paste0(data.dir ="/home/hailee/10X/19957/19957X2/outs/filtered_feature_bc_matrix"))
  s1[[i]] <- CreateSeuratObject(counts = x, project = samples$name[i])
}
for(i in 8:10){
  message("Loading ", samples$id[i])
  x <- Read10X(paste0(data.dir ="/home/hailee/10X/20778R/",samples$id[i],"/outs/filtered_feature_bc_matrix"))
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
CellFile <- merge(x = s1[[1]], y = s1[2:n], add.cell.ids = samples$name, project="All Cerebellum")

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

#Rename Clusters based on gene expression (genes listed in S2)
cell.type.markers <- c( "Gabra6", "Eomes","Ppp1r17", "Lypd6", "Prkcd", "Klhl1", "Lgi2", "Aqp4", "Mobp", "Ppfibp1", "C1qa", "Gdf10", "Ttr", "Dcn", "Kcnj8", "Flt1", "Foxj1")
DotPlot(CellFile, features = cell.type.markers, dot.scale = 8) +
  RotatedAxis()

##Find and remove any clusters that don't exist in at least 75% of samples or disproportianately belong to one sample (>85% come from one sample)
tb <- table(CellFile@meta.data$SCT_snn_res.0.8, CellFile@meta.data$orig.ident)
rs <- rowSums(table(CellFile@meta.data$SCT_snn_res.0.8, CellFile@meta.data$orig.ident))
(tb/rs)*100
CellFile <- subset(CellFile, idents = c(8), invert = TRUE)

### Change to correct identities
CellFile <- RenameIdents(CellFile, `0` = "Granule", `1` = "Granule", `2` = "Granule",`3` = "Granule", `4` = "Granule", 
                         `5` = "Granule", `6` = "Bergmann", `7` = "MLI1", `9` = "Astro",
                         `10` = "Granule", `11` = "MLI1",`12` = "Oligo", `13` = "PLI", `14` = "MLI2", `15` = "Oligo", 
                         `16` = "Purkinje2",`17` = "Oligo", `18` = "Fibroblast", `19` = "PLI", `20` = "Granule", `21` = "Golgi",
                         `22` = "Granule",`23` = "Granule", `24` = "Granule", `25` = "Choroid", `26` = "OPC",
                         `27` = "UBC",`28` = "Microglia" ,`29` = "Astro", `30` = "MLI2", `31` = "Ependymal", `32` = "Endo Stalk", 
                         `33` = "Granule", `34`= "Endo Mural", `35` = "Fibroblast")
CellFile$celltype <- Idents(CellFile)
Idents(CellFile) <- "celltype"

levels(CellFile) <- c("Granule", "UBC", "Purkinje", "MLI1", "MLI2", "PLI", "Golgi", "Astro", "Oligo", "OPC", "Microglia","Bergmann", "Choroid", 
                             "Fibroblast", "Endo Mural", "Endo Stalk", "Ependymal")

#create orig.id, genotype and status columns to use in split.by or group.by
CellFile$status <- substr(CellFile$orig.ident, 1,3)
CellFile$orig.ident <- factor(CellFile$orig.ident, levels=unique(m1$sample))
CellFile$genotype <- substr(CellFile$orig.ident, 5,6)
CellFile$condition <- substr(CellFile$orig.ident, 1,6)

save(CellFile, file = "/home/hailee/10X/Final Soc Analysis/Cereb-All.Rdata")

### WT Analysis
Sub <- subset(CellFile, subset = genotype == "WT")

## Figure 2D ## __________________________________________________________________________
#Create Feature plot of some marker genes
marks <- c("Slc17a7", "Gad2", "Gabra6", "Ppp1r17", "Lypd6", "Prkcd", "Klhl1", "Lgi2", "Aqp4", "Mobp", "C1qa", "Flt1")
FeaturePlot(Sub, marks, cols = c("grey", "darkred"))
  p <- GetAssayData(object=Sub, slot = "data")[marks[1],]
  pos_ids = names(which(p>1.4))
  pos_cells = subset(Sub,cells=pos_ids)
  DimPlot(pos_cells, cols = rep("blue", length(levels(pos_cells))))+NoLegend() + xlim(-10,16) + ylim(-17,15)
DimPlot(Sub, cols = rep("grey", length(levels(Sub))))+NoLegend() + xlim(-10,16) + ylim(-17,15)


color <- c("tomato1", "darkorange", "cyan4", "aquamarine3", "darkturquoise", "turquoise", "lightblue2", "darkmagenta", "blueviolet", "mediumpurple3", "mediumpurple1", "mediumorchid1","magenta3", 
           "magenta", "deeppink", "hotpink", "pink")

## Figure 2B ## _____________________________________________________________________________
#Make Violin Plot of gene expression
VlnPlot(Sub, cols = color, stack = TRUE, flip = TRUE, fill.by = "ident", features = cell.type.markers) + NoLegend()

## Figure 2A ## ______________________________________________________________________________
#Make UMAP with cell types labeled and colored
DimPlot(Sub, cols = color) + xlim(-10,16) + ylim(-17,15)

## Figure 2C ## ______________________________________________________________________________
##Make Bargraph of proportion of cells by type
md <- Sub@meta.data %>% as.data.table
md$celltype <- factor(md$celltype, levels = rev(levels(Sub)))

mc <- md[, mean(nFeature_RNA), by = c("celltype")]
md <- md[, .N, by = c("celltype")]
md <- md[order(md$celltype)]
mc <- mc[order(mc$celltype)]
WTsum <- sum(md$N)
md$totals <- rep(WTsum, nrow(md))
md$val <- (md$N/md$totals)

barplot <- barplot(as.numeric(md$val), logpropnames=md$celltype, horiz = TRUE, cex.names = 0.7, las = 1, col = rev(color), xlim = c(0,0.8)) 
for(b in 1:nrow(md)){
  l <- paste(substr(mc$V1[b],1,6))
  text(0.750,barplot[b], l, cex = 0.7)
  
}
