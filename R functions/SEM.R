SEM <- function(object, object2){
  data <- data.frame(NA_col = rep(NA, nrow(object)))
  
  for(i in 1:nrow(object)){  
    p1 <- object[i,]
    n1 <- object2[i,]
    p <- p1$val ##proportion of social
    n <- n1$N
    if (is.nan(p)) next
  sem = ((as.numeric(p)*(1-as.numeric(p)))/as.numeric(n))^0.5
  data[ i,] <- sem
  }
  return(data)
}