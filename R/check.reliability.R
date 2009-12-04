"check.reliability" <-
function(X){
   X <- check.data(X)
   compute.PP <- function(X,P,N,J,m){
     label <- as.vector(t(outer(paste("P(X",1:J,">=",sep=""),paste(1:(m-1),")",sep=""), paste, sep="")))
     PP <- matrix(0,J*(m-1),J*(m-1))
     i <- 0
     j <- 0
     for(i in 1:(J-1)) for(j in (i+1):J)
     PP[((i-1)*(m-1)+1):((i-1)*(m-1)+(m-1)),((j-1)*(m-1)+1):((j-1)*(m-1)+(m-1))] <-
      t(outer(X[,i],0:(m-2),">")) %*% outer(X[,j],0:(m-2),">")/N
     PP <- PP + t(PP) + kronecker(diag(J),matrix(-1,m-1,m-1))
     PP[PP < -.5] <- NA
     dimnames(PP) <- list(label,label)
     PP <- PP[order(P),order(P)]
     return(PP)
   }

   J <- ncol(X)
   N <- nrow(X)
   m <- max(X) + 1
   P <- matrix(t(apply(outer(as.matrix(X), 1:(m-1), ">=")*1,c(2,3),mean)),nrow=(m-1)*J)
   PP <- compute.PP(X,P,N,J,m)
   P <- matrix(sort(P))
   off.boundary <- P > 0 & P < 1
   P <- matrix(P[off.boundary,1])
   PP <- PP[off.boundary,off.boundary]
   km <- length(P)
   is.na.PP <- is.na(PP)
   lower.bound.PP <- P %*% t(P)
   upper.bound.PP <- outer(as.numeric(P),as.numeric(P),FUN = "pmin")
   dimnames(P) <- list(dimnames(PP)[[1]],"")

   OO <- is.na(PP)
   set.matrix <- outer(as.numeric(P),as.numeric(P),"==")*1
   unique.cells <- which(apply(set.matrix,1,sum)==1)
   Type <- set.matrix
   Type[unique.cells,unique.cells] <- 0
   set.vector <- sign(apply(Type,1,sum))
   Type <- outer(set.vector,set.vector) + outer(rep(1,km),set.vector)*2 + outer(set.vector,rep(1,km))

   for (i in 1:km) for (j in 1:km) if(is.na(PP[i,j])){
      if (Type[i,j]==4) PP[i,j] <- mean(PP[set.matrix[i,]==1,set.matrix[,j]==1],na.rm=T)
      else{
        if(Type[i,j]==1){
           RightN <- ifelse(any(!is.na(PP[i,j:km])), j + min(which(!is.na(PP[i,j:km]))) - 1, NA)
           RightPP <- ifelse(is.na(RightN),NA,mean(PP[set.matrix[i,]==1,set.matrix[,RightN]==1],na.rm=T))
           LeftN <- ifelse(any(!is.na(PP[i,1:j])), max(which(!is.na(PP[i,1:j]))),NA)
           LeftPP <- ifelse(is.na(LeftN),NA,mean(PP[set.matrix[i,]==1,set.matrix[,LeftN]==1],na.rm=T))
           LowerN <- UpperN <- i
           LowerPP <- UpperPP <- mean(PP[set.matrix[i,]==1,set.matrix[,j]==1],na.rm=T)
        }
        if(Type[i,j]==2){
           LowerN <- ifelse(any(!is.na(PP[i:km,j])), i + min(which(!is.na(PP[i:km,j]))) - 1, NA)
           LowerPP <- ifelse(is.na(LowerN),NA,mean(PP[set.matrix[LowerN,]==1,set.matrix[,j]==1],na.rm=T))
           UpperN <- ifelse(any(!is.na(PP[1:i,j])), max(which(!is.na(PP[1:i,j]))),NA)
           UpperPP <- ifelse(is.na(UpperN),NA,mean(PP[set.matrix[UpperN,]==1,set.matrix[,j]==1],na.rm=T))
           RightN <- LeftN <- j
           RightPP <- LeftPP <- mean(PP[set.matrix[i,]==1,set.matrix[,j]==1],na.rm=T)
        }
        if(Type[i,j]==0){
           RightN <- ifelse(any(!is.na(PP[i,j:km])), j + min(which(!is.na(PP[i,j:km]))) - 1, NA)
           RightPP <- ifelse(is.na(RightN),NA,mean(PP[set.matrix[i,]==1,set.matrix[,RightN]==1],na.rm=T))
           LeftN <- ifelse(any(!is.na(PP[i,1:j])), max(which(!is.na(PP[i,1:j]))),NA)
           LeftPP <- ifelse(is.na(LeftN),NA,mean(PP[set.matrix[i,]==1,set.matrix[,LeftN]==1],na.rm=T))
           LowerN <- ifelse(any(!is.na(PP[i:km,j])), i + min(which(!is.na(PP[i:km,j]))) - 1, NA)
           LowerPP <- ifelse(is.na(LowerN),NA,mean(PP[set.matrix[LowerN,]==1,set.matrix[,j]==1],na.rm=T))
           UpperN <- ifelse(any(!is.na(PP[1:i,j])), max(which(!is.na(PP[1:i,j]))),NA)
           UpperPP <- ifelse(is.na(UpperN),NA,mean(PP[set.matrix[UpperN,]==1,set.matrix[,j]==1],na.rm=T))
        }

        E17a <- ifelse(is.na(LowerN),NA,LowerPP * P[i]/P[LowerN])
        E17b <- ifelse(is.na(RightN),NA,RightPP * P[j]/P[RightN])
        E17c <- ifelse(is.na(UpperN),NA,UpperPP * P[i]/P[UpperN])
        E17d <- ifelse(is.na(LeftN),NA, LeftPP  * P[j]/P[LeftN] )
        E21a <- ifelse(is.na(LowerN),NA,LowerPP  * (1 - P[i])/(1 - P[LowerN]) -
                P[j] * (P[LowerN] - P[i])/(1 - P[LowerN]))
        E21b <- ifelse(is.na(RightN),NA,RightPP  * (1 - P[j])/(1 - P[RightN]) -
                P[i] * (P[RightN] - P[j])/(1 - P[RightN]))
        E21c <- ifelse(is.na(UpperN),NA,UpperPP  * (1 - P[i])/(1 - P[UpperN]) +
                P[j] * (P[i] - P[UpperN])/(1 - P[UpperN]))
        E21d <- ifelse(is.na(LeftN) ,NA,LeftPP   * (1 - P[j])/(1 - P[LeftN] ) +
                P[i] * (P[j] -  P[LeftN])/(1 - P[LeftN] ))

        PP[i,j] <- mean(c(E17a,E17b,E17c,E17d,E21a,E21b,E21c,E21d),na.rm=T)
     }
   }
   PP[is.nan(PP)] <- 0
   PP[PP > upper.bound.PP & OO] <- upper.bound.PP[PP > upper.bound.PP & OO]
   PP[PP < lower.bound.PP & OO] <- lower.bound.PP[PP < lower.bound.PP & OO]
   MS <-  sum(PP-lower.bound.PP)/(var(apply(X,1,sum))*((N-1)/N))
   varX <- var(X)
   alpha <- J/(J-1) * (sum(varX) - sum(diag(varX)))/sum(varX)
   lambda.2 <-  ((sum(varX)-sum(diag(varX))) + (sqrt((J/(J-1))*(sum(varX^2)-sum(diag(varX^2))))))/sum(varX)
   return(list(MS = MS, alpha = alpha, lambda.2 = lambda.2))
}
