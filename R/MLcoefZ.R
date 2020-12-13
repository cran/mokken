## Update 27-11-2020: 
# Argument type.z is added to distinguish between range-preserving and wald-based confidence intervals (Koopman et al 2020 a two-step, test-guided MSA for nonclustered and clustered data (Quality of Life Research))


"MLcoefZ" <- function(X, lowerbound = 0, type.z = "WB"){
  #X <- check.ml.data(X)
  if (type.z != "WB" & type.z != "RP"){
    warning("type.z needs to be 'WB' (Wald-based) or 'RP' (range-preserving), the default 'WB' is used.")
    type.z <- "WB"
  }  
  Hs <- MLcoefH(X, nice.output = FALSE)
  if(type.z == "WB") {
    Zij <- (Hs[[1]][, c(1, 3, 5)] - lowerbound) / 
      Hs[[1]][, c(2, 4, 6)]
    Zi <- (Hs[[2]][, c(1, 3, 5)] - lowerbound) / 
      Hs[[2]][, c(2, 4, 6)]
    Z <- (Hs[[3]][, c(1, 3, 5)] - lowerbound) / 
      Hs[[3]][, c(2, 4, 6)]
  } else {
    Zij <- -(log(1 - Hs[[1]][, c(1, 3)]) - log(1 - lowerbound)) / 
      (Hs[[1]][, c(2, 4)] / (1 - Hs[[1]][, c(1, 3)]))
    diag(Zij) <- 0
    Zi <- -(log(1 - Hs[[2]][, c(1, 3)]) - log(1 - lowerbound)) / 
      (Hs[[2]][, c(2, 4)] / (1 - Hs[[2]][, c(1, 3)]))
    Z <- -(log(1 - Hs[[3]][, c(1, 3)]) - log(1 - lowerbound)) / 
      (Hs[[3]][, c(2, 4)] / (1 - Hs[[3]][, c(1, 3)]))
  }
  return(list(Zij=Zij,Zi=Zi,Z=Z))
}
