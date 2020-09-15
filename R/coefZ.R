# Aangepast op 12 september 2020

"coefZ" <- function(X, lowerbound = 0, type.se = "Z", level.two.var = NULL){
  X <- check.data(X)
  if (type.se != "Z" & type.se != "delta"){
    warning(paste("type.se = '", type.se, "' unknown. Type of standard error has been changed to typse.se = 'Z'", sep = ""))
    type.se <- "Z"
  }  
  if(!is.null(level.two.var)) {
    if (nrow(as.matrix(level.two.var)) != nrow(X)) {
      level.two.var <- NULL
      warning("level.two.var not the same length/nrow as X: level.two.var is ignored.")
    } 
    if(any(is.na(level.two.var))) {
      level.two.var <- NULL
      warning("level.two.var contains missing value(s): level.two.var is ignored.")
    } 
    if(type.se == "Z") {
      type.se <- "delta"
      warning("type.se has been changed to 'delta' to enable testing in multilevel data.")
    }
    X <- X[order(level.two.var), ]
    level.two.var <- sort(level.two.var)
    Rs <- as.numeric(table(level.two.var))
    level.two.var <- rep(1:length(Rs), Rs)
    # Ensure each subject has > 1 rater
    if(any(Rs == 1)){ 
      warning('For at least one group there is only 1 respondent. The Z-scores are computed without this (these) group(s).') 
      cases <- !(level.two.var %in% which(Rs == 1))
      X <- X[cases, ]
      level.two.var <- level.two.var[cases]
    }
  }
  if(type.se == "delta") {
    Hs <- coefH(X, nice.output = FALSE, level.two.var = level.two.var, print.to.screen = FALSE)
    Zij <- -(log(1 - Hs[[1]]) - log(1 - lowerbound)) / (Hs[[2]] / (1 - Hs[[1]]))
    diag(Zij) <- 0
    Zi <- -matrix((log(1 - Hs[[3]]) - log(1 - lowerbound)) / (Hs[[4]] / (1 - Hs[[3]])), nrow = 1)
    Z <- -(log(1 - Hs[[5]]) - log(1 - lowerbound)) / (Hs[[6]] / (1 - Hs[[5]]))
  } else {
    if(lowerbound > 0) {
      type.se <- "delta"
      warning("type.se has been changed to 'delta' to enable testing for lowerbound > 0.")
    } else {
      N <- nrow(X)
      S <- var(X)
      Smax <- var(apply(X,2,sort))
      Sij <- outer(apply(X,2,var),apply(X,2,var),"*")
      Zij <- (S * sqrt(N-1))/sqrt(Sij)
      diag(S) <- diag(Sij) <- diag(Zij) <- 0
      Zi <- (apply(S,1,sum) * sqrt(N-1))/ sqrt(apply(Sij,1,sum))       
      Z  <- (sum(S)/2 * sqrt(N-1))/ sqrt(sum(Sij)/2)
    }
  } 
  return(list(Zij = Zij, Zi = Zi, Z = Z))
}


#"coefZ.wald" <- function(X, lowerbound = 0, type.se = "delta", level.two.var = NULL){
#  X <- check.data(X)
#  if(!is.null(level.two.var)) {
#    if (nrow(as.matrix(level.two.var)) != nrow(X)) {
#      level.two.var <- NULL
#      warning("level.two.var not the same length/nrow as X: level.two.var is ignored.")
#    } 
#    if(any(is.na(level.two.var))) {
#      level.two.var <- NULL
#      warning("level.two.var contains missing value(s): level.two.var is ignored.")
#    } 
#    if(type.se == "Z") {
#      type.se <- "delta"
#      warning("type.se has been changed to 'delta' to enable testing in multilevel data.")
#    }
#    X <- X[order(level.two.var), ]
#    level.two.var <- sort(level.two.var)
#    Rs <- as.numeric(table(level.two.var))
#    level.two.var <- rep(1:length(Rs), Rs)
#    # Ensure each subject has > 1 rater
#    if(any(Rs == 1)){ 
#      warning('For at least one group there is only 1 respondent. The Z-scores are computed without this (these) group(s).') 
#      cases <- !(level.two.var %in% which(Rs == 1))
#      X <- X[cases, ]
#      level.two.var <- level.two.var[cases]
#    }
#  }
#  if(type.se == "Z") {
#    if(lowerbound > 0) {
#      type.se <- "delta"
#      warning("type.se has been changed to 'delta' to enable testing for lowerbound > 0.")
#    } else {
#      X <- check.data(X)
#      N <- nrow(X)
#      S <- var(X)
#      Smax <- var(apply(X,2,sort))
#      Sij <- outer(apply(X,2,var),apply(X,2,var),"*")
#      Zij <- (S * sqrt(N-1))/sqrt(Sij)
#      diag(S) <- diag(Sij) <- diag(Zij) <- 0
#      Zi <- (apply(S,1,sum) * sqrt(N-1))/ sqrt(apply(Sij,1,sum))       
#      Z  <- (sum(S)/2 * sqrt(N-1))/ sqrt(sum(Sij)/2)
#    }
#  } 
#  if(type.se == "delta") {
#    Hs <- coefH(X, nice.output = FALSE, level.two.var = level.two.var, print.to.screen = FALSE)
#    Zij <- (Hs[[1]] - lowerbound) / Hs[[2]]
#    diag(Zij) <- 0
#    Zi <- matrix((Hs[[3]] - lowerbound) / Hs[[4]], nrow = 1)
#    Z <- (Hs[[5]] - lowerbound) / Hs[[6]]
#  } else if(type.se != "Z") {
#    stop("please specify type.se as 'Z' or 'delta'")  
#  }
#  return(list(Zij=Zij,Zi=Zi,Z=Z))
#}

#"coefZ.old" <- function(X){
#  X <- check.data(X)
#  N <- nrow(X)
#  S <- var(X)
#  Sij <- outer(apply(X,2,var),apply(X,2,var),"*")
#  Zij <- (S * sqrt(N-1))/sqrt(Sij)
#  diag(S) <- diag(Sij) <- diag(Zij) <- 0
#  Zi <- (apply(S,1,sum) * sqrt(N-1))/ sqrt(apply(Sij,1,sum))       
#  Z  <- (sum(S)/2 * sqrt(N-1))/ sqrt(sum(Sij)/2)
#  return(list(Zij=Zij,Zi=Zi,Z=Z))
#}
