"check.data" <-
function(X){
   if (data.class(X) != "matrix" && data.class(X) != "data.frame")
     stop("Data are not matrix or data.frame")
    matrix.X <- as.matrix(X)
    if (any(is.na(matrix.X))) stop("Missing values are not allowed")
    if (mode(matrix.X)!="numeric") stop("Data must be numeric")
    if (any(matrix.X < 0)) stop("All scores should be nonnegative")
    if (any(matrix.X %% 1 !=0)) stop("All scores must be integers")
    matrix.X <- matrix.X - min(matrix.X)
    return(matrix.X)
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

"weights" <-
# X: Data matrix N x 2 of integer scores [0,1, ..., maxx]
# w: Guttman weights 1 x g^2
# depends on "all.patterns"

function(X, maxx=max.x, minx=0){
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
 order.of.ISRFs <- matrix(rank(-tmp.2),1,maxx*2)
# Compute
 Y <- matrix(all.patterns(2,g),nrow=1)
 Z <- matrix(rep(Y, maxx), nrow = maxx, byrow = TRUE)
 Z <- ifelse(Z < row(Z),0,1)
 Z <- matrix(as.vector(Z), ncol = maxx*2, byrow = T)
# COMPUTE WEIGHTS
 Z <- Z[,order(order.of.ISRFs)]
 w <- matrix(apply(Z,1,function(x){sum(x*cumsum(abs(x-1)))}),nrow=1)
 return(w)
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
