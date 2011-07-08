"coefH" <-
function(X, se = TRUE, nice.output = TRUE){
    X <- check.data(X)
    eps <- 1e-40
    labels <- dimnames(X)[[2]]
    if (se == FALSE){
      S <- var(X)
      Smax <- var(apply(X, 2, sort))
      Hij <- S/Smax
      diag(S) <- 0
      diag(Smax) <- 0
      Hi <- apply(S, 1, sum)/apply(Smax, 1, sum)
      H <- sum(S)/sum(Smax)
      return(list(Hij=Hij,Hi=Hi,H=H))  
    } else {
      g <- max(X) - min(X) + 1
      J <- ncol(X)
      P <- choose(J,2)
      N <- nrow(X)
      if(any(apply(X,2,var) < eps))stop("One or more variables have zero variance")
      n <- as.matrix(table(apply(X,1,paste, collapse="")))
      lab.n <- matrix(names(table(apply(X,1,paste, collapse=""))))
      lab.b <- apply(all.patterns(2,g),2,paste,collapse = "")
      lab.u <- as.character(0:(g-1))
      Bi <- substr(lab.b,1,1)
      Bj <- substr(lab.b,2,2)

      R <- t(apply(lab.n,1,string2integer))
      r <- length(lab.n)

      U <- list()     # Univariate frequencies per item (1 x g)
      for (j in 1:J) U[[j]] <- count(X[,j],0:(g-1))

      W <- list()     # Item weights per item-pair (1 x g^2)
      WA <- list()
      WY <- list()
      WE <- list()
      WF <- list()
      for (i in 1:J){
        W[[i]] <- list()
        WA[[i]] <- list()
        WY[[i]] <- list()
        WE[[i]] <- list()
        WF[[i]] <- list()
        for(j in i:J) if (j > i){
           W[[i]][[j]] <- weights(X[,c(i,j)])                                                            # [1]
           A1a <- NULL
           for (a in 0:(g-1)) for (b in 0:(g-1))  A1a <- rbind(A1a,as.numeric(R[,i]==a & R[,j]==b))      # [2]
           WA[[i]][[j]] <- W[[i]][[j]] %*% A1a                                                           # [3]
           # Y22 (B x U)
           Eij <- matrix(U[[i]][as.numeric(Bi)+1],nrow=g^2,ncol=1) *  matrix(U[[j]][as.numeric(Bj)+1],nrow=g^2,ncol=1)/N
           Y22 <- cbind(outer(Bi,lab.u,"=="),outer(Bj,lab.u,"==")) * matrix(Eij, nrow=g^2, ncol=2*g)
 
           Ri <- substr(lab.n,i,i)                                                                       # [4]
           Rj <- substr(lab.n,j,j)
           Z2 <- rbind(outer(lab.u, Ri,"==")[,,1], outer(lab.u, Rj,"==")[,,1])* c(1/U[[i]],1/U[[j]])
           Z2[is.nan(Z2)] <- 1/eps
           YZ2 <- Y22 %*% Z2 - Eij %*% matrix(1/N,1,r)

          WY[[i]][[j]] <- W[[i]][[j]] %*% YZ2                                                           # [5]
           Fij <- complete.observed.frequencies(X[,c(i,j)],2,g)
           WF[[i]][[j]] <- W[[i]][[j]] %*% Fij
           WE[[i]][[j]] <- W[[i]][[j]] %*% Eij
        }
      }
      # Hij
      g3 <- matrix(c(unlist(WF[[1]][[2]]),unlist(WF),unlist(WE)),nrow = 2*P + 1,byrow=TRUE)
      A4 <- rbind(matrix(c(1,-1,rep(0,(J*(J-1))-1)),1,(J*(J-1))+1),cbind(matrix(0,P,1),diag(P),-1*diag(P)))
      A5 <- cbind(matrix(1,P,1), -1*diag(P))
      g4 <- phi(A4,g3,"log")
      g5 <- phi(A5,g4,"exp")
      G3 <- matrix(c(unlist(WA[[1]][[2]]),unlist(WA),unlist(WY)),nrow = 2*P + 1,byrow=TRUE)
      G4 <- dphi(A4,g3,G3,"log")
      G5 <- dphi(A5,g4,G4,"exp")
      Hij <- se.Hij <- matrix(0,J,J)
      Hij[lower.tri(Hij)] <- g5
      Hij <- Hij + t(Hij)
      se.Hij[lower.tri(se.Hij)] <- sqrt(diag(G5 %*% (as.numeric(n) * t(G5))))
      se.Hij <- se.Hij + t(se.Hij)
      dimnames(se.Hij)  <- dimnames(Hij) <- list(labels,labels)
      # H
      if (P > 1) G3 <- rbind(apply(G3[2:(P+1),],2,sum),apply(G3[2:(P+1),],2,sum),apply(G3[(P+2):(2*P+1),],2,sum))
      g3 <- matrix(c(sum(g3[2:(P+1),]),sum(g3[2:(P+1),]),sum(g3[(P+2):(2*P+1),])),ncol=1)
      A4 <-cbind(matrix(1,2,1),-1*diag(2))
      A5 <- matrix(c(1,-1),1,2)
      g4 <- phi(A4,g3,"log")
      g5 <- phi(A5,g4,"exp")
      G4 <- dphi(A4,g3,G3,"log")
      G5 <- dphi(A5,g4,G4,"exp")
      H <- g5
      se.H <- sqrt(diag(G5 %*% (as.numeric(n) * t(G5))))
      
      # Hi
      G3 <- matrix(0,2*J + 1,r)
      g3 <- matrix(0,2*J + 1,1)
      for (j in 1:J) for (k in 1:J) if (k > j) {
         g3[j+1,] <- g3[j+1,] + WF[[j]][[k]]
         g3[J+j+1,] <- g3[J+j+1,] + WE[[j]][[k]]
         G3[j+1,] <- G3[j+1,] + WA[[j]][[k]]
         G3[J+j+1,] <- G3[J+j+1,] + WY[[j]][[k]]
      } else { if (k < j) {
         g3[j+1,] <- g3[j+1,] + WF[[k]][[j]]
         g3[J+j+1,] <- g3[J+j+1,] + WE[[k]][[j]]
         G3[j+1,] <- G3[j+1,] + WA[[k]][[j]]
         G3[J+j+1,] <- G3[J+j+1,] + WY[[k]][[j]]
      }}
      g3[1,] <- g3[2,]
      G3[1,] <- G3[2,]
      A4 <- rbind(matrix(c(1,-1,rep(0,(J*2)-1)),1,(J*2)+1),cbind(matrix(0,J,1),diag(J),-1*diag(J)))
      A5 <- cbind(matrix(1,J,1),-1*diag(J))
      g4 <- phi(A4,g3,"log")
      g5 <- phi(A5,g4,"exp")
      G4 <- dphi(A4,g3,G3,"log")
      G5 <- dphi(A5,g4,G4,"exp")
      Hi <- matrix(g5)
      se.Hi <- matrix(sqrt(diag(G5 %*% (as.numeric(n) * t(G5)))))
      dimnames(se.Hi)[[1]] <- dimnames(Hi)[[1]] <- labels
          
      if (nice.output){
        output.matrix.Hij <- matrix(NA,J,J*2)
        for (j in 2*(1:J)){
          output.matrix.Hij[,j-1] <- format(paste(" ",formatC(round(Hij[,j/2]   ,3),digits=3,format='f')," ",sep=""), width = 7, justify="right")
          output.matrix.Hij[,j] <-   format(paste("(",formatC(round(se.Hij[,j/2],3),digits=3,format='f'),")",sep=""), width = 7, justify="right")
        }  
        new.labels <- rep(labels,each=2)
        new.labels[2*(1:J)] <- "se"
        dimnames(output.matrix.Hij)[[1]] <- labels      
        dimnames(output.matrix.Hij)[[2]] <- new.labels      
        output.matrix.Hij[row(output.matrix.Hij) ==  .5 * col(output.matrix.Hij)] <- format("",width = 7, justify="right")
        output.matrix.Hij[(row(output.matrix.Hij) ) ==  (.5 * col(output.matrix.Hij) + .5)] <- format("",width = 7, justify="right") 
        output.matrix.Hij <- noquote(output.matrix.Hij)
        output.matrix.Hi <- matrix(NA,J,2)
        output.matrix.Hi[,1] <- format(formatC(round(Hi,3),digits=3,format='f'),width=7, justify="right")
        output.matrix.Hi[,2] <- format(paste("(",formatC(round(se.Hi,3),digits=3,format='f'),")",sep=""), width = 7, justify="right")
        dimnames(output.matrix.Hi) <- list(labels, c("Item H","se"))
        output.matrix.Hi <- noquote(output.matrix.Hi)
        output.matrix.H  <- matrix(NA,1,2)
        output.matrix.H[,1] <- format(formatC(round(H,3),digits=3,format='f'),width=7, justify="right")
        output.matrix.H[,2] <- format(paste("(",formatC(round(se.H,3),digits=3,format='f'),")",sep=""), width = 7, justify="right")
        dimnames(output.matrix.H) <- list("", c("Scale H","se"))
        output.matrix.H <- noquote(output.matrix.H)
        return(list(Hij = output.matrix.Hij,Hi = output.matrix.Hi,H = output.matrix.H))
      }  else return(list(Hij=Hij,se.Hij=se.Hij,Hi=Hi,se.Hi=se.Hi,H=H,se.H=se.H))
   }
}
          

string2integer <- function(s) as.numeric(unlist(strsplit(s,NULL)))

all.patterns <- function(J,m){
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

weights <-
# X: Data matrix N x 2 of integer scores [0,1, ..., maxx]
# w: Guttman weights 1 x g^2
# depends on "count", "all.patterns"
function(X, maxx=max.x, minx=0){
 max.x <- max(X)
 g <- maxx + 1
 N <- nrow(X)
 if (ncol(X) != 2){
   warning('X contains more than two columns. Only first two columns will be used')
   X <- X[,1:2]
 }
# Compute order of the ISRFs
 tmp.1 <- apply(X,2,count, seq(minx,maxx))
 tmp.1 <- matrix(tmp.1[-1,],maxx,2)
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

count <-
# counts the number of elements in x equal to elm; e.g.
# > x <- c(1,1,2,2,2,3,3)
# > y <- count(x,c(1,2))
# > y
# [1] 2  3
function(x, elm=sort(unique(x))){
  j <- 0;  res <- 0
  if(any(is.na(x)))stop("data contain missing values")
  unique.elm <- unique(elm)
  for (j in 1:length(unique.elm))res[j] <- length(x[x==unique.elm[j]])
  return(res)
}

phi <- function(A,f, action){
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
dphi <- function(A,f,df, action){
  eps=1E-80;
  switch(action,
    "identity" = A %*% df,
    "exp"      = A %*% (as.numeric(exp(f)) * df),
    "log"      = A %*% (as.numeric(1/(f+eps)) * df),
    "sqrt"     = A %*% (as.numeric(1/(2*sqrt(f))) * df),
    "xlogx"    = A %*% (as.numeric(-1-log(f+eps)) * df),
    "xbarx"    = A %*% (as.numeric(1-2*f) * df),  # x(1-x)
  )
}

complete.observed.frequencies <- function(data,J,m, order.items=FALSE){
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
