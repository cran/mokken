"check.pmatrix" <-
function(X, minvi = .03){

   compute.pmatrix <- function(X,P1,N,J,m){
      label <- as.vector(t(outer(paste("P(X",1:J,">=",sep=""),paste(1:(m-1),")",sep=""), paste, sep="")))
      Pmm <- Ppp <- matrix(0,J*(m-1),J*(m-1))
      i <- 0
      j <- 0
      for(i in 1:(J-1)){
         for(j in (i+1):J){
            Ppp[((i-1)*(m-1)+1):((i-1)*(m-1)+(m-1)),((j-1)*(m-1)+1):((j-1)*(m-1)+(m-1))] <-
            t(outer(X[,i],0:(m-2),">")) %*% outer(X[,j],0:(m-2),">")/N
            Pmm[((i-1)*(m-1)+1):((i-1)*(m-1)+(m-1)),((j-1)*(m-1)+1):((j-1)*(m-1)+(m-1))] <-
            t(outer(X[,i],0:(m-2),"<=")) %*% outer(X[,j],0:(m-2),"<=")/N
         }
      }
      Ppp <- Ppp + t(Ppp) + kronecker(diag(J),matrix(-1,m-1,m-1))
      Pmm <- Pmm + t(Pmm) + kronecker(diag(J),matrix(-1,m-1,m-1))
      Ppp[Ppp < -.5] <- NA
      Pmm[Pmm < -.5] <- NA
      dimnames(Ppp) <- dimnames(Pmm) <- list(label,label)
      Ppp <- Ppp[order(P1),order(P1)]
      Pmm <- Pmm[order(P1),order(P1)]
      return(list(Ppp=Ppp,Pmm=Pmm))
   }

  X <- check.data(X)
  J <- ncol(X)
  N <- nrow(X)
  m <- max(X) + 1
  P1 <- matrix(t(apply(outer(as.matrix(X), 1:(m-1), ">=")*1,c(2,3),mean)),nrow=(m-1)*J)
  I.item <- rep(1:J,each=m-1)[order(P1)]
  I.step <- as.vector(t(outer(paste("X",1:J,">=",sep=""),paste(1:(m-1),sep=""), paste, sep="")))[order(P1)]
  I.labels <- dimnames(X)[[2]]
  if(length(I.labels)==0) I.labels <- paste("C",1:ncol(X))

  P2 <- compute.pmatrix(X,P1,N,J,m)
  Pmm <- P2$Pmm
  Ppp <- P2$Ppp
  Hi <- coefH(X)$Hi
  pmatrix.list <- list(Ppp=Ppp,Pmm=Pmm, I.item=I.item, I.step=I.step, I.labels=I.labels, Hi=Hi, minvi=minvi)
  class(pmatrix.list) <- "pmatrix.class"  
  return(pmatrix.list)
  }

