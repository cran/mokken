"check.iio" <-
function(X, minvi = m * .03, minsize = default.minsize, alpha = .05){

  X <- check.data(X)
  J <- ncol(X)
  N <- nrow(X)
  m <- max(X) + 1
  default.minsize <- ifelse(N > 500, floor(N/10), floor(N/5))
  default.minsize <- ifelse(N <= 250, floor(N/3), default.minsize)
  default.minsize <- ifelse(N <  150, 50, default.minsize)

  # Is the sample size large enough?
  if (N < minsize) stop("Sample size less than minsize")

  # Are there enough items?
  if (J < 3) stop("Less than 3 items. Restscore cannot be computed")

  # Put the items in descending order
  item.order <- rev(order(colMeans(X)))
  X <- X[,item.order]
  I.labels <- dimnames(X)[[2]]
  if(length(I.labels)==0) I.labels <- paste("C",1:ncol(X))
  vi.matrix <- matrix(0,J,J)
  dimnames(vi.matrix) <- list(I.labels,I.labels)

  g <- 0; i <- 0; j <- 0; k <- 0
  # Check violations of Method Restscore for each item pair
  for (i in 1:(J-1)){
    R <- as.matrix(X) %*% (matrix(1,J,J) - diag(J)) - X[,i]
    R[,i] <- 0
    for (j in (i+1):J){
      k <- k + 1
      sorted.R <- sort(R[,j])
      group <- max(which(sorted.R==sorted.R[minsize]))
      repeat{
        if(N - max(group) < minsize)break
        group <- c(group,max(which(sorted.R==sorted.R[minsize+max(group)])))
      }
      group <- group[-length(group)]  
      group <- c(sorted.R[group],max(sorted.R))
      L <- length(group)
      member <- apply(1 - outer(R[,j], group, "<="),1,sum) + 1

      for (g in 1:L){
         EXi.R <- mean(X[member==g,i])
         EXj.R <- mean(X[member==g,j])
         greater.than.minvi <- ifelse(EXj.R - EXi.R <= minvi,FALSE,TRUE) 
         if(greater.than.minvi){ 
            tt <- t.test(X[member==g,i],X[member==g,j],alternative="less")
            vi.matrix[i,j] <- vi.matrix[i,j] + 1*(tt$p.value < alpha)
#           cat(I.labels[i],I.labels[j],"g=",g,"d=",format(round(EXj.R - EXi.R,2),nsmall=2), "T=",format(round(tt$statistic,2),nsmall=2),"p=",format(tt$p.value,nsmall=3),fill=T)
         }  
      }
   }
 }
 vi.matrix <- vi.matrix + t(vi.matrix)
 vi.matrix.a <- vi.matrix
 VI <- matrix(apply(sign(vi.matrix),1,sum)) 
 items.removed <- NULL
 repeat{
  nvi <- apply(sign(vi.matrix.a),2,sum)
  if (sum(nvi)==0) break
  maxvi <- which(nvi ==max(nvi))
  if(length(maxvi) > 1){ 
     H <- rep(0,length(maxvi))
     for (i in 1:length(maxvi)) H[i] <- coefH(X[,-c(items.removed,maxvi[i])])$H
     maxvi <- maxvi[which(H==min(H))[1]]
  } 
  vi.matrix.a[maxvi,] <- vi.matrix.a[,maxvi] <- 0
  items.removed <- c(items.removed,maxvi)
  VI. <- matrix(apply(sign(vi.matrix.a),1,sum)) 
  VI.[items.removed,] <- NA
  VI <- cbind(VI,VI.)
 }
dimnames(VI) <- list(I.labels,paste("step",1:ncol(VI)))
HT <- ifelse(length(items.removed)==0,coefH(t(X))$H,coefH(t(X[,-items.removed]))$H) 
return(list(violations=VI,items.removed=items.removed,HT = HT))
}
