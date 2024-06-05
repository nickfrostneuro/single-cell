chi2test <- function(object, object2){
  
  data <- data.frame(NA_col = rep(NA, nrow(object)))
  
  for(i in 1:nrow(object)){
    if((i %% 2) == 0)next  
    c <- object[i,]
    n <- object[(i+1),]
    n1 <- c$val ##proportion of social
    n2 <- n$val ##proportion of NS
    
    c <- object2[i,]
    n <- object2[(i+1),]
    N1 <- as.numeric(c$N) ##number of social
    N2 <- as.numeric(n$N) ##number of NS
    
    #Observed data
    n1 <- as.numeric(n1) * as.numeric(N1)
    n2 <- as.numeric(n2)*as.numeric(N2)
    # Pooled estimate of proportion
    p0 <- (n1+n2) / (N1+N2)
    # Expected counts under H0 (null hypothesis)
    n10 <- N1 * p0
    n20 <- N2 * p0
    # Chi-square test, by hand
    observed <- c(n1,N1-n1, n2, N2-n2)
    
    expected <- c(n10, N1-n10, n20, N2-n20)
    
    chi2stat <- sum((observed-expected)^2 / expected)
    
    p <- 1 - pchisq(chi2stat,1)
    data[ i,] <- p                     # Adding new variable to data
    
  }
  return(data)
}
