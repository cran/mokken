
coefHi <- function(X){
  S <- var(X)
  Smax <- var(apply(X, 2, sort))
  diag(S) <- 0
  diag(Smax) <- 0
  return(apply(S, 1, sum)/apply(Smax, 1, sum))
}

compute.defaultminsize <- function(Y){
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  N <- 0
  if(class(Y) == "matrix") N <- nrow(Y) 
  if(class(Y) == "numeric") if(is.wholenumber(Y)) N <- Y 
  default.minsize <- ifelse(N > 500, floor(N / 10), floor(N / 5))
  default.minsize <- ifelse(N <= 250, floor(N / 3), default.minsize)
  default.minsize <- ifelse(N <  150, 50, default.minsize)
  default.minsize <- ifelse(N <  50, NA, default.minsize)
  return(default.minsize)
}  
  

joinRSG <- function(R, minsize = default.minsize){
  N <- length(R)
  default.minsize <- compute.defaultminsize(N)
  if (is.na(minsize)) stop("N < 50. Too few respondents")
  sorted.R <- sort(R)
  group <- max(which(sorted.R == sorted.R[minsize]))
  repeat{
    if(N - max(group) < minsize) break
    group <- c(group, max(which(sorted.R == sorted.R[minsize + max(group)])))
  }
  group <- group[-length(group)]
  group <- c(sorted.R[group], max(sorted.R))
  L <- length(group)
  return(apply(1 - outer(R, group, "<="), 1, sum) + 1)
}

covweight <- function(Y, N, option = "pnorm"){
  if(any(apply(Y, 2, sd) == 0)) weight <- 0 else {
    rab <- cor(Y[, 1], Y[, 2])
    m <- 0.5 * log((1 + rab) / (1 - rab + 1e-40))
    s <- 1 / sqrt(N - 3)
    weight <- switch(option,
     pnorm    = pnorm(-m / s),
     noweight = as.numeric(rab < 0))
  }   
  return(weight)
}

nweight <- function(N, option = "noweight") switch(option, sqrt = sqrt(N), 1)

# Collect the rest scores groups for each item pair in R
compute.restscores <- function(X, minsize = default.minsize){
  J <- ncol(X)
  default.minsize <- compute.defaultminsize(nrow(X))
  if (is.na(minsize)) stop("N < 50. Too few respondents")
  R <- list()
  for (i in 1 : J) R[[i]] <- list()
  for (i in 1 : (J - 1)) for (j in (i + 1) : J)
  R[[i]][[j]] <- R[[j]][[i]] <- joinRSG(apply(X[, -c(i, j)], 1, sum), minsize) 
  return(R)
}  

compute.W1 <- function(X, mingroup = 4, covweightoption = "pnorm", nweightoption = "noweight"){
  N <- nrow(X)
  J <- ncol(X)
  scores <- sort(unique(as.numeric(X)))
  M <- length(scores)
  W1 <- matrix(NA, J, J)
  for (a in 1 : J) for (c in 1 : J) if (a != c) {
    RES <- matrix(0, J, M)
    for (j in 1 : J) if (j != a & j != c) for (x in 1 : M){
      Group <- X[, c] == scores[x]
      GroupN <- sum(Group) 
      RES[j, x] <- ifelse(GroupN >= mingroup, covweight(X[Group, c(a, j)], GroupN, covweightoption), 0) * nweight(GroupN, nweightoption)
    }
    W1[a, c] <- sum(RES)
  }
  return(W1)
}

compute.W2 <- function(X, R, mingroup = 4, covweightoption = "pnorm", nweightoption = "noweight"){
  N <- nrow(X)
  J <- ncol(X)
  W2 <- matrix(NA, 1, J)
  for (a in 1 : J) {
    RES <- matrix(0, J, 1)
    for (j in 1 : J) if (j != a){
      RM <- max(R[[a]][[j]]) 
      RES2 <- matrix(0, RM, 1)
      for (r in 1 : RM){
         Group <- R[[a]][[j]] == r
         GroupN <- sum(Group) 
         RES2[r, 1] <- ifelse(GroupN >= mingroup, covweight(X[Group, c(a, j)], GroupN, covweightoption), 0) * nweight(GroupN, nweightoption)
      }
      RES[j, 1] <- sum(RES2) 
    }  
    W2[1, a] <- sum(RES)
  }
  return(W2)
}

compute.W3 <- function(X, R, mingroup = 4, covweightoption = "pnorm", nweightoption = "noweight"){
  N <- nrow(X)
  J <- ncol(X)
  W3 <- matrix(NA, J, J)
  for (a in 1 : J) for (b in 1 : J) if (a < b) {
    RM <- max(R[[a]][[b]]) 
    RES <- matrix(0, RM, 1)
    for (r in 1 : RM){
      Group <- R[[a]][[b]] == r
      GroupN <- sum(Group) 
      RES[r, 1] <- ifelse(GroupN >= mingroup, covweight(X[Group, c(a, b)], GroupN, covweightoption), 0) * nweight(GroupN, nweightoption)
    }
    W3[a, b] <- W3[b, a] <- sum(RES)
  }
  return(W3)
}

flag <- function(Y){
  flag.criterion <- quantile(Y, na.rm = TRUE)[4] + 3 * (quantile(Y, na.rm = TRUE)[4] - quantile(Y, na.rm = TRUE)[3]) 
  Y[is.na(Y)] <- 0
  return(Y > flag.criterion)
}  

"check.ca" <- function(X, Windex = FALSE, MINSIZE = 4, NWEIGHTOPTION = "noweight", COVWEIGHTOPTION = "pnorm", MINGROUP = 4){
 X <- check.data(X)  
 item.labels <- dimnames(X)[[2]]
 J <- ncol(X)
 inset <- rep(TRUE, J)
 RES <- list()
 RES$Flagged <- RES$Index <- RES$InScale <- list()
 k <- 0
 NAmatrix  <- matrix(NA, J, J)
 dimnames(NAmatrix) <- list(item.labels, item.labels)
 NAvector  <- matrix(NA, 1, J)
 dimnames(NAmatrix)[[2]] <- item.labels
 
 repeat{
    k <- k + 1
    RES[[1]][[k]] <- inset
    RES[[2]][[k]]<- list()
    RES[[3]][[k]]<- list()
    Xin <- X[, inset]
    R  <- compute.restscores(Xin, minsize = MINSIZE)

    W1 <- compute.W1(Xin, mingroup = MINGROUP, covweightoption = COVWEIGHTOPTION, nweightoption = NWEIGHTOPTION) 
    RES[[2]][[k]]$W1 <- NAmatrix
    RES[[2]][[k]]$W1[inset, inset] <- W1

    W2 <- compute.W2(Xin, R, mingroup = MINGROUP, covweightoption = COVWEIGHTOPTION, nweightoption = NWEIGHTOPTION) 
    RES[[2]][[k]]$W2 <- NAvector
    RES[[2]][[k]]$W2[1, inset] <- W2

    W3 <- compute.W3(Xin, R, mingroup = MINGROUP, covweightoption = COVWEIGHTOPTION, nweightoption = NWEIGHTOPTION) 
    RES[[2]][[k]]$W3 <- NAmatrix
    RES[[2]][[k]]$W3[inset, inset] <- W3

    F1 <- flag(W1)
    RES[[3]][[k]]$F1 <- NAmatrix
    RES[[3]][[k]]$F1[inset, inset] <- as.numeric(F1)

    F2 <- flag(W2)
    RES[[3]][[k]]$F2 <- NAvector
    RES[[3]][[k]]$F2[1, inset] <- as.numeric(F2)

    F3 <- flag(W3)
    RES[[3]][[k]]$F3 <- NAmatrix
    RES[[3]][[k]]$F3[inset, inset] <- as.numeric(F3)
    
    flags <- apply(F1, 1, sum) + apply(F1, 2, sum) + F2 + apply(F3, 1, sum)
    if (sum(flags) == 0) break else{ candidates <-  which(flags == max(flags))}
    if (length(candidates) > 1){
       H <- coefHi(Xin) 
       candidates <- candidates[which(H[candidates] == min(H[candidates]))]
       if (length(candidates) > 1) candidates <- sample(candidates, 1)
    }
    inset[inset][candidates] <- FALSE
    if(sum(inset) < 4) {inset <- rep(FALSE, J); break}
  }  
  if (Windex) return(RES) else return(RES[[1]])
}  
