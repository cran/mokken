# Aangepast op 6 maart 2020

"coefZ" <- function(X, lowerbound = 0, type.se = "delta", level.two.var = NULL){
  if(type.se == "Z") {
    if(!is.null(level.two.var)) warning("level.two.var is ignored for type.se = 'Z'.")
    X <- check.data(X)
    N <- nrow(X)
    S <- var(X)
    Smax <- var(apply(X,2,sort))
    Sij <- outer(apply(X,2,var),apply(X,2,var),"*")
    Zij <- sqrt(N-1) * (S/Smax - lowerbound)/(sqrt(Sij)/Smax)
    diag(S) <- diag(Sij) <- diag(Zij) <- diag(Smax) <- 0
    Zi <- sqrt(N-1) * (apply(S,1,sum) / apply(Smax,1,sum) - lowerbound)/ (sqrt(apply(Sij,1,sum)) / (apply(Smax,1,sum)))
    Z  <- sqrt(N-1) * (sum(S) / (sum(Smax)) - lowerbound)/ (sqrt(sum(Sij)/2) / (sum(Smax) / 2))
  } else if(type.se == "delta") {
    Hs <- coefH(X, nice.output = FALSE, level.two.var = level.two.var)
    Zij <- (Hs[[1]] - lowerbound) / Hs[[2]]
    diag(Zij) <- 0
    Zi <- matrix((Hs[[3]] - lowerbound) / Hs[[4]], nrow = 1)
    Z <- (Hs[[5]] - lowerbound) / Hs[[6]]
  } else {
    stop("please specify type.se as 'Z' or 'delta'")  
  }
  return(list(Zij=Zij,Zi=Zi,Z=Z))
}

"coefZ.old" <- function(X){
  X <- check.data(X)
  N <- nrow(X)
  S <- var(X)
  Sij <- outer(apply(X,2,var),apply(X,2,var),"*")
  Zij <- (S * sqrt(N-1))/sqrt(Sij)
  diag(S) <- diag(Sij) <- diag(Zij) <- 0
  Zi <- (apply(S,1,sum) * sqrt(N-1))/ sqrt(apply(Sij,1,sum))       
  Z  <- (sum(S)/2 * sqrt(N-1))/ sqrt(sum(Sij)/2)
  return(list(Zij=Zij,Zi=Zi,Z=Z))
}
