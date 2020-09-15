# Aangepast op 1 september 2020.

"MLcoefZ" <- function(X, lowerbound = 0){
  #X <- check.ml.data(X)
  Hs <- MLcoefH(X, nice.output = FALSE)
  Zij <- -(log(1 - Hs[[1]][, c(1, 3)]) - log(1 - lowerbound)) / 
    (Hs[[1]][, c(2, 4)] / (1 - Hs[[1]][, c(1, 3)]))
  diag(Zij) <- 0
  Zi <- -(log(1 - Hs[[2]][, c(1, 3)]) - log(1 - lowerbound)) / 
    (Hs[[2]][, c(2, 4)] / (1 - Hs[[2]][, c(1, 3)]))
  Z <- -(log(1 - Hs[[3]][, c(1, 3)]) - log(1 - lowerbound)) / 
    (Hs[[3]][, c(2, 4)] / (1 - Hs[[3]][, c(1, 3)]))
  return(list(Zij=Zij,Zi=Zi,Z=Z))
}
