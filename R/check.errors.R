"check.errors" <- function(X, returnGplus = TRUE, returnOplus = FALSE){

medCouple <- function(x){
  x <- sort(x)
  median.x <- median(x)
  xp <- x[x >=  median.x]
  xn <- x[x <=  median.x]
  p <- length(xp)
  n <- length(xn)
  rn <- 1 : n
  rp <- 1 : p
  h <- outer(xn, xp, function(xn, xp, median.x){((xp - median.x) - (median.x - xn))/(xp - xn)}, median.x)
  S <- outer(xn, xp, "==") 
  h[S] <- sign(outer(rn, rp, function(rn, rp, p){p - 1 - rp - rn}, p))[S]
  return(median(h))
}

  X <- check.data(X)
  maxx <- max(X)
  minx <- min(X)
  N <- nrow(X)
  J <- ncol(X)

  if(returnGplus) {
     # Y = vec(X') 
     Y <- matrix(t(X), 1, N * J)
     # Z = 1 %*% Y
     Z <- matrix(rep(Y, maxx), maxx, N * J, TRUE)
     # Z = Z uitgeschreven als X = 4 = 1 1 0 0 
     Z <- ifelse(Z < row(Z), 0, 1)
     Z <- matrix(as.vector(Z), N, maxx * J, TRUE)
     if (maxx == 1) tmp.1 <- matrix(apply(X, 2, tabulate, maxx), nrow = 1) else tmp.1 <- apply(X, 2, tabulate, maxx)
     tmp.2 <- apply(tmp.1, 2, function(x) rev(cumsum(rev(x))))
     tmp.3 <- as.numeric(matrix(rank(-tmp.2), 1, maxx * J))
     # tmp.3 is a vector with the order of the ISRFs
     if(any(duplicated(tmp.3))) {
        Gplusx <- matrix(NA, N, 1000)
        set.seed(1)
        for (i in 1 : 1000){ 
          tmp.2x <- tmp.2 + runif(length(tmp.2), -.001, .001)
          tmp.3x <- as.numeric(matrix(rank(-tmp.2x), 1, maxx * J))
          Zx <- Z[, order(tmp.3x)]
          Gplusx[, i] <- apply(Zx, 1, function(x){sum(x * cumsum(abs(x - 1)))})
        }
        Gplus <- round(apply(Gplusx, 1, mean))  
        Gplus
     }  else { 
        Zx <- Z[, order(tmp.3)]
        Gplus <- apply(Zx, 1, function(x){sum(x * cumsum(abs(x - 1)))})
     }   
     Q3 <- summary(Gplus)[[5]]
     IQR <- Q3 - summary(Gplus)[[2]]
     U1Gplus <- Q3 + 1.5 * IQR
     U2Gplus <- U1Gplus * exp(3.87 * medCouple(Gplus))
  }
  if(returnOplus) {
     O <- apply(apply(apply(X + 1, 2, tabulate), 2, rank), 2, rev) - 1
     OScores <- X
     for (j in 1 : J) OScores[, j] <- O[X[, j] + 1, j]
     Oplus <- apply(OScores, 1, sum)
     Q3 <- summary(Oplus)[[5]]
     IQR <- Q3 - summary(Oplus)[[2]]
     U1Oplus <- Q3 + 1.5 * IQR
     U2Oplus <- U1Oplus * exp(3.87 * medCouple(Oplus))
  }
  if(returnGplus && returnOplus) return(list(Gplus = Gplus, UGplus = list(U1Gplus = U1Gplus, U2Gplus = U2Gplus), Oplus = Oplus, UOplus = list(U1Oplus = U1Oplus, U2Oplus = U2Oplus))) else {
    if (returnGplus) return(list(Gplus = Gplus, UGplus = list(U1Gplus = U1Gplus, U2Gplus = U2Gplus))) else {
      if (returnOplus) return(list(Oplus = Oplus, UOplus = list(U1Oplus = U1Oplus, U2Oplus = U2Oplus))) else return(list()) }}
}


     ##  items <- rep(1:J, each = maxx)
     ##  steps <- rep(1 : maxx, J)
     ##  ranks <- matrix(c(items, steps, tmp.3), nrow = 3, byrow = T)
     ##  which.ties <- outer(tmp.3, tmp.3, "==")
     ##  which.ties[row(which.ties) >= col(which.ties)] <- FALSE
     ##  sameItem <- matrix(1, maxx, maxx) 
     ##  for (i in 1 : (J - 1)) sameItem <- direct.sum(sameItem, matrix(1, maxx, maxx))
     ##  which.ties <- which.ties  & (sameItem == 0)
     ##  n.ties <-     
