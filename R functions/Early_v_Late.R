## Compares both the Early and Late time points to baseline
##Output will be a list of multiple outputs. You can then pull each out using list[["name"]]

Early_v_Late <- function(mIEG, md, Sub){
  Not <- mIEG[(mIEG$status == "Not"),]
  ShS <- mIEG[(mIEG$status == "ShS"),]
  Soc <- mIEG[(mIEG$status == "Soc"),]
  Early <- rbind(Not, ShS)
  Early <- Early[order(Early$celltype, Early$status), ]
  Late <- rbind(Not, Soc)
  Late <- Late[order(Late$celltype, Late$status), ]
  
  Not <- md[(md$status == "Not"),]
  ShS <- md[(md$status == "ShS"),]
  Soc <- md[(md$status == "Soc"),]
  Earlymd <- rbind(Not, ShS)
  Earlymd <- Earlymd[order(Earlymd$celltype, Earlymd$status), ]
  Latemd <- rbind(Not, Soc)
  Latemd <- Latemd[order(Latemd$celltype, Latemd$status), ]
  
  mat <- matrix(nrow=nrow(Early), ncol = 1)
  for(i in 1:nrow(Early)){
    if((i %% 2) == 0)next
    l <- Early$val[(i+1)]/Early$val[i]
    mat[i,] = l
  }
  Early1 <- cbind(Early$status, as.character(Early$celltype), as.numeric(mat))
  row_odd <- seq_len(nrow(Early1)) %% 2 
  Early1 <- Early1[row_odd == 1, ]
  colnames(Early1) <- c("status", "celltype", "val")
  Early1 <- as.data.frame(Early1)
  
  Early1$celltype <- factor(Early1$celltype, levels = rev(levels(Sub)))
  Early1$status <- rep("10 min. / Not", nrow(Early1))
  Early1 <- Early1[order(Early1$celltype, Early1$status), ]
  
  mat <- matrix(nrow=nrow(Late), ncol = 1)
  for(i in 1:nrow(Late)){
    if((i %% 2) == 0)next
    l <- Late$val[(i+1)]/Late$val[i]
    mat[i,] = l
  }
  Late1 <- cbind(Late$status, as.character(Late$celltype), as.numeric(mat))
  row_odd <- seq_len(nrow(Late1)) %% 2 
  Late1 <- Late1[row_odd == 1, ]
  colnames(Late1) <- c("status", "celltype", "val")
  Late1 <- as.data.frame(Late1)
  
  Late1$celltype <- factor(Late1$celltype, levels = rev(levels(Sub)))
  Late1$status <- rep("35 min. / Not", nrow(Late1))
  colnames(Late1) <- c("status", "celltype", "val")
  Late1 <- Late1[order(Late1$celltype, Late1$status), ]
  
  mIEG1 <- rbind(Early1, Late1)
  return(mIEG1)
}
