"summary.pmatrix.class" <-
function(object, ...){

  ppp.vi. <- function(X){
    i <- 0
    nr.vi <- 0
    for (i in 1: nrow(X)){
      x <- X[i,]
      z <- outer(x,x, function(x,y,criterion){x-y > criterion},criterion=minvi)
      z[is.na(z)] <- F
      z[row(z) >= col(z)] <- F 
      nr.vi <- nr.vi + sum(z)
    }
   return(nr.vi) 
  }  

  pmm.vi. <- function(X){
    i <- 0
    nr.vi <- 0
    for (i in 1: nrow(X)){
      x <- X[i,]
      z <- outer(x,x, function(x,y,criterion){x-y < -criterion},criterion=minvi)
      z[is.na(z)] <- F
      z[row(z) >= col(z)] <- F 
      nr.vi <- nr.vi + sum(z)
    }
   return(nr.vi) 
  }  

  ppp.vi.. <- function(X,minvi){
    i <- 0
    max.vi <- 0
    sum.vi <- 0
    for (i in 1: nrow(X)){
      x <- outer(X[i,],X[i,],"-")
      x <- x[row(x) < col(x)]
      max.vi <- max(max.vi,max(x))
      sum.vi <- sum(sum.vi,sum(x[x > minvi]))
    }
  return(list(max.vi=max.vi,sum.vi=sum.vi)) 
  }    

  pmm.vi.. <- function(X,minvi){
    i <- 0
    max.vi <- 0
    sum.vi <- 0
    for (i in 1: nrow(X)){
      x <- outer(X[i,],X[i,],"-")
      x <- x[row(x) > col(x)]
      max.vi <- max(max.vi,max(x))
      sum.vi <- sum(sum.vi,sum(x[x > minvi]))
    }
  return(list(max.vi=max.vi,sum.vi=sum.vi)) 
  }    

  minvi <- object$minvi 
  J <- max(object$I.item)

  ppp.summary.matrix <- matrix(nrow=J,ncol=7)
  dimnames(ppp.summary.matrix) <- list(object$I.labels,c("ItemH","#ac","#vi","#vi/#ac","maxvi","sum","sum/#ac"))
  ppp.summary.matrix[,1] <- round(object$Hi,2) 
  Ppp <- object$Ppp
  I.item <- object$I.item
  items <- 1:max(I.item)
  j <- 1
  for (j in items){
     Ppp.j <- Ppp[I.item==j,I.item!=j] 
     if(!is.matrix(Ppp.j)) Ppp.j <- t(as.matrix(Ppp.j))
     ppp.summary.matrix[j,2] <- length(Ppp.j) # ac
     ppp.summary.matrix[j,3] <- ppp.vi.(Ppp.j)
     tmp <- ppp.vi..(Ppp.j,minvi)
     ppp.summary.matrix[j,5] <- ifelse(tmp$max.vi > minvi,round(tmp$max.vi,2),0)
     ppp.summary.matrix[j,6] <- round(tmp$sum.vi,2)
  }     
  ppp.summary.matrix[,4] <- round(ppp.summary.matrix[,5]/ppp.summary.matrix[,2],3)    
  ppp.summary.matrix[,7] <- round(ppp.summary.matrix[,6]/ppp.summary.matrix[,2],3)    
  
  pmm.summary.matrix <- matrix(nrow=J,ncol=7)
  dimnames(pmm.summary.matrix) <- list(object$I.labels,c("ItemH","#ac","#vi","#vi/#ac","maxvi","sum","sum/#ac"))
  pmm.summary.matrix[,1] <- round(object$Hi,2) 
  Pmm <- object$Pmm
  I.item <- object$I.item
  items <- 1:max(I.item)
  j <- 1
  for (j in items){
     Pmm.j <- Pmm[I.item==j,I.item!=j] 
     if(!is.matrix(Pmm.j)) Pmm.j <- t(as.matrix(Pmm.j))
     pmm.summary.matrix[j,2] <- length(Pmm.j) # ac
     pmm.summary.matrix[j,3] <- pmm.vi.(Pmm.j)
     tmp <- pmm.vi..(Pmm.j,minvi)
     pmm.summary.matrix[j,5] <- ifelse(tmp$max.vi > minvi,round(tmp$max.vi,2),0)
     pmm.summary.matrix[j,6] <- round(tmp$sum.vi,2)
  }     
  pmm.summary.matrix[,4] <- round(pmm.summary.matrix[,5]/pmm.summary.matrix[,2],3)    
  pmm.summary.matrix[,7] <- round(pmm.summary.matrix[,6]/pmm.summary.matrix[,2],3)    
  
  return(list(ppp.summary.matrix=ppp.summary.matrix,pmm.summary.matrix=pmm.summary.matrix))  
}

