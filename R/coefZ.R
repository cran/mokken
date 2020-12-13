## Update 27-11-2020: 
# Argument type.z is added to distinguish between range-preserving and wald-based confidence intervals (Koopman et al 2020 a two-step, test-guided MSA for nonclustered and clustered data (Quality of Life Research))


"coefZ" <- function(X, lowerbound = 0, type.z = "Z", level.two.var = NULL){
  X <- check.data(X, FALSE)
  if (type.z != "Z" & type.z != "WB" & type.z != "RP"){
    type.z.t <- ifelse(lowerbound > 0, "WB", "Z")
    warning(paste("type.z = '", type.z, "' is unknown. ",
    "Type of z-score was changed to type.z = 'Z' if lowerbound = 0 and to type.z = 'WB' if lowerbound > 0. ",
    "In this case, type.z was changed to '", type.z.t, "'", sep = ""))
    type.z <- type.z.t
  }  
  if(lowerbound > 0 & type.z == "Z") {
    type.z <- "WB"
    warning("type.z has been changed to 'WB' to enable testing for lowerbound > 0.")
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
    if(type.z == "Z") {
      type.z <- "WB"
      warning("type.z has been changed to 'WB' to enable testing in multilevel data.")
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
  if(type.z == "WB") {
    Hs <- suppressWarnings(coefH(X, nice.output = FALSE, level.two.var = level.two.var, results = FALSE))
    Zij <- (Hs[[1]] - lowerbound) / Hs[[2]]
    diag(Zij) <- 0
    Zi <- matrix((Hs[[3]] - lowerbound) / Hs[[4]], nrow = 1)
    Z <- (Hs[[5]] - lowerbound) / Hs[[6]]
  } else if(type.z == "RP") {
    Hs <- suppressWarnings(coefH(X, nice.output = FALSE, level.two.var = level.two.var, results = FALSE))
    Zij <- -(log(1 - Hs[[1]]) - log(1 - lowerbound)) / (Hs[[2]] / (1 - Hs[[1]]))
    diag(Zij) <- 0
    Zi <- -matrix((log(1 - Hs[[3]]) - log(1 - lowerbound)) / (Hs[[4]] / (1 - Hs[[3]])), nrow = 1)
    Z <- -(log(1 - Hs[[5]]) - log(1 - lowerbound)) / (Hs[[6]] / (1 - Hs[[5]]))
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
