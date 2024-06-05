## Function to build counts matrix out of list of variables extracted from Seurat Object
##stim.list = the seurat object split by the original sample identities

Seurat_Counts_Matrix <- function(stim.list){
  library(foreach)
  library(doParallel)
  totalCores = detectCores()
  cluster <- makeCluster(totalCores[1]-15) #If nothing else is running (-1) if some is reduce cores by more (-15 or more)
  registerDoParallel(cluster)
  
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
  return(Counts)
}