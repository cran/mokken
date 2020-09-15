"check.data" <- function(X){
   if (data.class(X) != "matrix" && data.class(X) != "data.frame")
     stop("Data are not matrix or data.frame")
    matrix.X <- as.matrix(X)
    if (any(is.na(matrix.X))) stop("Missing values are not allowed")
    if (mode(matrix.X)!="numeric") stop("Data must be numeric")
    if (any(matrix.X < 0)) stop("All scores should be nonnegative")
    if (any(matrix.X %% 1 !=0)) stop("All scores must be integers")
    matrix.X <- matrix.X - min(matrix.X)
    if (max(matrix.X > 9)) stop("Some items have more than 10 categories. mokken cannot only handle up to 10 categories")
    score.check <- apply(apply(matrix.X, 2, function(x) seq(0, max(matrix.X), 1) %in% x), 2, all)
    if (!all(score.check)) 
    warning("Varying numbers of item scores were observed across the items.\n  Either the items have the same number of response categories but some item categories were not endorsed;\n  or the items have different numbers of categories by design. \n  In the latter case, the sum score cannot be used for (ordinal) measurement.")
    return(matrix.X)
}

"check.ml.data" <- function(X){
  if (data.class(X) != "matrix" && data.class(X) != "data.frame")
    stop("Data are not matrix or data.frame")
  matrix.X <- as.matrix(X)
  if (any(is.na(matrix.X))) stop("Missing values are not allowed")
  if (mode(matrix.X)!="numeric") stop("Data must be numeric")
  if (any(matrix.X < 0)) stop("All scores should be nonnegative")
  if (any(matrix.X %% 1 !=0)) stop("All scores must be integers")
  matrix.X[, -1] <- matrix.X[, -1] - min(matrix.X[, -1])
    score.check <- apply(apply(matrix.X[, -1], 2, function(x) seq(0, max(matrix.X[, -1]), 1) %in% x), 2, all)
    if (!all(score.check)) 
    warning("Varying numbers of item scores were observed across the items.\n  Either the items have the same number of response categories but some item categories were not endorsed;\n  or the items have different numbers of categories by design. \n  In the latter case, the sum score cannot be used for (ordinal) measurement.")
  return(matrix.X)
}

"coefHTiny" <- function(X){
    S <- var(X)
    Smax <- var(apply(X, 2, sort))
    Hij <- S/Smax
    diag(S) <- 0
    diag(Smax) <- 0
    Hi <- apply(S, 1, sum)/apply(Smax, 1, sum)
    H <- sum(S)/sum(Smax)
    return(list(Hij = Hij, Hi = Hi, H = H))
}

"all.patterns" <- function(J,m){
  grid <- list()
  j <- 0;
  p <- m^J
  for (j in 1:J){
    grid <- c(grid, j)
    grid[[j]] <- 0:(m-1)
  }
  X <- t(expand.grid(grid))
  dimnames(X) <- NULL
  return(X[J:1,])
}

"phi" <- function(A,f, action){
  # Numerical values are translations h(A %*% f) = A %*% f -
  eps = 1E-80;
  switch(action,
         "identity" = A %*% f,
         "exp"      = A %*% exp(f),
         "log"      = A %*% log(abs(f)+eps),
         "sqrt"     = A %*% sqrt(f),
         "xlogx"    = A %*% (-f*log(f+eps)),
         "xbarx"    = A %*% (f*(1-f))  # x(1-x)
  )
}

"dphi" <- function(A,f,df, action){
  eps=1E-80;
  switch(action,
         "identity" = A %*% df,
         "exp"      = A %*% (as.numeric(exp(f)) * df),
         "log"      = A %*% (as.numeric(1/(f+eps)) * df),
         "sqrt"     = A %*% (as.numeric(1/(2*sqrt(f))) * df),
         "xlogx"    = A %*% (as.numeric(-1-log(f+eps)) * df),
         "xbarx"    = A %*% (as.numeric(1-2*f) * df)  #, x(1-x)
  )
}

"string2integer" <- function(s) as.numeric(unlist(strsplit(s,NULL)))

"oldweights" <- function(X, maxx = max.x, minx = 0, itemstep.order = NULL){
 max.x <- max(X)
 g <- maxx + 1
 N <- nrow(X)
 if (ncol(X) != 2){
   warning('X contains more than two columns. Only first two columns will be used')
   X <- X[,1:2]
 }
# Compute order of the ISRFs
 if (maxx == 1) tmp.1 <- matrix(apply(X,2,tabulate, maxx), nrow=1) else tmp.1 <- apply(X,2,tabulate, maxx)
 tmp.2 <- apply(tmp.1,2,function(x) rev(cumsum(rev(x))))+runif(2*maxx,0,1e-3)

 # runif is added to avoid equal ranks
 if (is.null(itemstep.order)) order.of.ISRFs <- matrix(rank(-tmp.2), 1, maxx * 2) else order.of.ISRFs <- matrix(rank(itemstep.order), 1, maxx * 2) 
# Compute
 Y <- matrix(all.patterns(2,g),nrow=1)
 Z <- matrix(rep(Y, maxx), nrow = maxx, byrow = TRUE)
 Z <- ifelse(Z < row(Z),0,1)
 Z <- matrix(as.vector(Z), ncol = maxx*2, byrow = T)
# COMPUTE WEIGHTS
 Z <- Z[, order(order.of.ISRFs)]
 w <- matrix(apply(Z, 1, function(x){sum(x * cumsum(abs(x - 1)))}), nrow = 1)
 return(w)
}
# Decrepit as of 25-11-2019

## weights function 21-11-2019 by Letty Koopman, adjusted original supporting function by Renske Kuijpers in coefH(). 
"weights" <- function(X, maxx = max.x, minx = 0, itemstep.order = NULL){
  # Computes the Guttman weights in Mokken Scale Analysis.
  #
  # Args:
  #   X: Data matrix N x 2 of integer scores [0,1, ..., maxx]
  #   maxx: The highest possible answer category. If not specified it is determined by using the highest item score.
  #   minx: The lowest possible answer category, default 0.
  #   w: Guttman weights 1 x g^2
  #   Depends on "all.patterns".
  #
  # Returns:
  #   Guttman weights
  
  # Error handling:
  max.x <- max(X)
  g <- maxx + 1
  N <- nrow(X)
  if (ncol(X) != 2){
    warning('X contains more than two columns. Only first two columns will be used')
    X <- X[,1:2]
  }
  
  # Compute frequencies of the item scores
  Rel1 <- table(factor(X[, 1], levels = minx:maxx))
  names(Rel1) <- paste0(1, minx:maxx)
  Rel2 <- table(factor(X[, 2], levels = minx:maxx))
  names(Rel2) <- paste0(2, minx:maxx)
  
  # Cumulative frequencies, dealing with equal ranks
  CumRel <- c(rev(cumsum(rev(Rel1[-1]))), rev(cumsum(rev(Rel2[-1]))))
  names(CumRel) <- 1:length(CumRel)
  y <- sort(CumRel, decreasing = T)
  
  if (any(duplicated(y)) & is.null(itemstep.order)) {
    perm <- function(n, r, v = 1:n, set = TRUE) {
      if (r == 1) 
        matrix(v, n, 1)
      else if (n == 1) 
        matrix(v, 1, r)
      else {
        X <- NULL
        for (i in 1:n) X <- rbind(X, cbind(v[i], Recall(n - 
                                                          1, r - 1, v[-i])))
        X
      }
    }
    o <- lapply(unique(y), function(val) {
      m <- as.numeric(names(y[y == val]))
      if (length(m) <= 1) {
        return(m)
      }
      as.data.frame(perm(length(m), length(m), m))
    })
    # Ensure fixed ordering within items
    for (i in 1:length(o)) {
      g <- o[[i]]
      if (length(unique(g)) > 1) {
        select <- matrix(0, nrow(g))
        for (j in 1:nrow(g)) {
          h <- g[j, ]
          if (any(h >= 1 & h <= maxx) & any(h >= maxx + 
                                            1 & h <= maxx * 2)) {
            select[j] <- (all(h[which(h >= 1 & h <= maxx)] == 
                                sort(h[which(h >= 1 & h <= maxx)])) & all(h[which(h >= 
                                                                                    maxx + 1 & h <= maxx * 2)] == sort(h[which(h >= 
                                                                                                                                 maxx + 1 & h <= maxx * 2)]))) * j
          }
          else {
            i1 <- ifelse(any(h >= 1 & h <= maxx), all(h[which(h >= 
                                                                1 & h <= maxx)] == sort(h[which(h >= 1 & 
                                                                                                  h <= maxx)])), 0)
            i2 <- ifelse(any(h >= maxx + 1 & h <= maxx * 
                               2), all(h[which(h >= maxx + 1 & h <= maxx * 
                                                 2)] == sort(h[which(h >= maxx + 1 & h <= 
                                                                       maxx * 2)])), 0)
            select[j] <- (i1 + i2) * j
          }
        }
        o[[i]] <- matrix(apply(o[[i]][select, ], 1, paste0, 
                               collapse = "."))
      }
    }
    out <- matrix(apply(expand.grid(o), 1, paste0, collapse = "."))
    out <- matrix(unlist(strsplit(out, "[.]")), nrow = nrow(out), 
                  byrow = T)
    w <- NULL
    Z <- matrix(rep(matrix(all.patterns(2, maxx + 1), nrow = 1), 
                    maxx), nrow = maxx, byrow = TRUE)
    Z <- matrix(ifelse(Z < row(Z), 0, 1), ncol = (maxx) * 
                  2, byrow = TRUE)
    for (i in 1:nrow(out)) {
      ords <- as.numeric(out[i, ])
      # Compute Z matrix for each possible item-response pattern
      Z1 <- Z[, (ords)]
      # Compute weights
      w <- rbind(w, apply(Z1, 1, function(x) {
        sum(x * cumsum(abs(x - 1)))
      }))
    }
    # Compute average weight per pattern
    wr <- matrix(colMeans(w), nrow = 1)
  }
  else {
    if (is.null(itemstep.order)) ords <- as.numeric(names(y)) else ords <- matrix(rank(itemstep.order), 1, maxx * 2)
    Z <- matrix(rep(matrix(all.patterns(2, maxx + 1), nrow = 1), 
                    maxx), nrow = maxx, byrow = TRUE)
    Z <- matrix(ifelse(Z < row(Z), 0, 1), ncol = (maxx) * 
                  2, byrow = TRUE)
    Z <- Z[, (ords)]
    wr <- apply(Z, 1, function(x) {
      sum(x * cumsum(abs(x - 1)))
    })
  }
  return(wr)
  
}

"complete.observed.frequencies" <- function(data,J,m, order.items=FALSE){
  if(order.items) order <- rev(order(apply(data,2,mean))) else order <- 1:J
  data <- as.matrix(data[,order])
  t.R <- cbind(t(all.patterns(J,m)),0)
  p <- m^J
  N <- nrow(data)
  for (i in 1:p){
    size <- abs(data - matrix(1,N,1) %*% t.R[i,1:J]) %*% matrix(1,J,1) == 0
    t.R[i,J+1] <- length(size[size==TRUE])
  }
  return(matrix(t.R[,J+1]))
}

"direct.sum" <- function (...){
     p.tr = 0;p.ll = 0;
     matlist = list(...);
     nmat = length(matlist);
     m1 = matlist[[1]];
     matlist = if(nmat==1 && is.list(m1)) m1 else matlist # check if list of matrices is given and amend accordingly
     nmat = length(matlist);                              # ,,
     m1 = matlist[[1]];                                   # ,,
     if(nmat==1) return(m1);
     for(i in 2:nmat){
        m2 = matlist[[i]];
        topleft <- m1
        topright <- matrix(p.tr, nrow(m1), ncol(m2))
        colnames(topright) <- colnames(m2)
        lowleft <- matrix(p.ll, nrow(m2), ncol(m1))
        lowright <- m2
        m1 = rbind(cbind(topleft, topright), cbind(lowleft, lowright))
     }
     return(m1)
}

"print.scalCoefsSE" <- function(x){
  maxl <- length(x)
  print(x[1:(maxl - 3)])
  invisible(x[(maxl - 2):maxl])
}  
