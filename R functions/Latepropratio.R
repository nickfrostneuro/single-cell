##md.CT = subset of metadata for N by cell type
##CT.subset = subset of Seurat for cell type
Latepropratio <- function(md.CT, CT.subset){
  ##Sample the number of cells for not social and find the proportion IEG+
  md.sub <- md.CT[(md.CT$status == "Not"),]
  set.seed(as.integer(Sys.time()) %% 10000)
  sampled.cells <- sample(colnames(CT.subset), size = md.sub$N , replace = F)
  obj.sub <- subset(CT.subset, cells = sampled.cells)
  
  mat <- matrix(ncol=length(genes), nrow= ncol(obj.sub))
  for(i in 1:length(genes)){
    tmp <- GetAssayData(object=obj.sub, slot = "data")[genes[i], ]
    q <- (tmp > quantile(tmp, 0.95))
    q1 <- as.integer(q)
    mat[,i] <- q
  }
  name <- names(q)
  multIEG <- rowSums(mat)
  aa <- as.numeric(as.integer(multIEG>3))
  tab <- cbind(name, aa)
  if(sum(aa) == 0){
    #Save to use for a distribution where row is cell type and column is the iteration
    IEG.pos.not <- NA
  } else{
    pos_ids = tab[tab[,2]==1]
    v1 <- subset(Sub, cells = pos_ids)
    
    IEG.pos.not <- ncol(v1)/ncol(obj.sub)}
  rm(v1, pos_ids, tab, aa, multIEG, name, mat, obj.sub, sampled.cells, md.sub,i)
  ##Sample the number of cells for 35 min social and find the proportion IEG+
  md.sub <- md.CT[(md.CT$status == "Soc"),]
  set.seed(as.integer(Sys.time()) %% 10000)
  sampled.cells <- sample(colnames(CT.subset), size = md.sub$N, replace = F)
  obj.sub <- subset(CT.subset, cells = sampled.cells)
  
  mat <- matrix(ncol=length(genes), nrow= ncol(obj.sub))
  for(i in 1:length(genes)){
    tmp <- GetAssayData(object=obj.sub, slot = "data")[genes[i], ]
    q <- (tmp > quantile(tmp, 0.95))
    q1 <- as.integer(q)
    mat[,i] <- q
  }
  name <- names(q)
  multIEG <- rowSums(mat)
  aa <- as.numeric(as.integer(multIEG>3))
  tab <- cbind(name, aa)
  if(sum(aa) == 0){
    #Save to use for a distribution where row is cell type and column is the iteration
    IEG.pos.35min <- NA
  } else{
    pos_ids = tab[tab[,2]==1]
    v1 <- subset(Sub, cells = pos_ids)
    
    IEG.pos.35min <- ncol(v1)/ncol(obj.sub)}
  rm(v1, pos_ids, tab, aa, multIEG, name, mat, obj.sub, sampled.cells, md.sub,i)
  
  ##Calculate a proportion for each of the social conditions
  Late <- IEG.pos.35min/IEG.pos.not
  return(Late)
}