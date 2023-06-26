# Adapted Feb 19, 2020
"check.iio" <- 
function (X, method="MIIO", minvi = default.minvi, minsize = default.minsize, 
          alpha = .05, item.selection=TRUE, verbose=FALSE, level.two.var = NULL)
{
   # CHECK DATA
   X <- check.data(X)

   # READ 
   J <- ncol(X)
   N <- nrow(X)
   ncat <- max(X) + 1 # CHECK VOOR IT EN MS-CPM (WAARSCHIJNLIJK GOED) VOOR MIIO WAARSCHIJNLIJK m <- max(X)
  
   if(!is.null(level.two.var)) {
     X <- X[order(level.two.var), ]
     level.two.var <- sort(level.two.var)
     if(method != "MIIO"){
       warning("Only method MIIO is available for level.two.var, method is changed to MIIO.")
       method <- "MIIO"
     }
   } 
   
   # DETERMINE METHOD
   if (substr(method,1,1)=="I" || substr(method,1,1)=="i") method="IT" else
   if (substr(method,1,2)=="MS"|| substr(method,1,2)=="ms" ||substr(method,1,1)=="C" ||substr(method,1,1)=="c") method="MSCPM" else method <- "MIIO"
   
   # DETERMINE MINVI and MINSIZE
   default.minvi <- ifelse(method=="MIIO",(ncat-1)*.03,.03)
   default.minsize <- ifelse(N >= 500, floor(N/10), floor(N/5))
   default.minsize <- ifelse(N <= 250, floor(N/3), default.minsize)
   default.minsize <- ifelse(N < 150, 50, default.minsize)
   # minvi <- default.minvi
   # minsize <- default.minsize
    
   # STOP IF THERE ARE NOT ENOUGH ITEMS OR RESPONDENTS
   if (N < minsize) stop("Sample size less than Minsize")
   if (J < 3) stop("Less than 3 items. Restscore cannot be computed")

   # INITIAL VALUES
   results <- list()
   vi.matrix <- matrix(0,J,J)
   item.order <- rev(order(colMeans(X)))
   if(is.null(dimnames(X)[[2]])) dimnames(X)[[2]] <- paste("V", 1:ncol(X),sep="")
   X <- X[, item.order]
   I.labels <- dimnames(X)[[2]]
   dimnames(vi.matrix) <- list(I.labels,I.labels)
   Hi <- coefHTiny(X)$Hi
   g <- 0
   gg <- 0
   h <- 0
   i <- 1
   j <- 2
   k <- 0
   Ll1 <- NULL
   item.mean <- colMeans(X)
   
   # METHOD IT
   if (method=="IT"){
     ii <- jj <- kk <- 0
     uv <- matrix(0,nrow=choose(ncat,2),ncol=2)
     for (ii in 0:(ncat-2)) for (jj in (ii+1):(ncat-1)) {kk <- kk + 1; uv[kk,] <- c(ii,jj)}
     for (i in 1:(J - 1)) {
        R <- as.matrix(X) %*% (matrix(1, J, J) - diag(J)) - X[,i] # R_(ij) for j != i
        R[, i] <- 0
        for (j in (i + 1):J) {
           k <- k + 1
           rvm <- nrow(uv) + 1 #D#
           violation.matrix <- matrix(0, nrow = rvm, ncol = 8)
           dimnames(violation.matrix) <- list(c(t(paste("P(X",i,"=",uv[,2],",X",j,"=",uv[,1],") P(X",i,"=",uv[,1],",X",j,"=",uv[,2],")" ,sep="")),"Total"), #D#
                                             c("#ac", "#vi", "#vi/#ac", "maxvi", "sum", "sum/#ac","X^2 max", "#X^2 sig"))
           results[[k]] <- list()
           results[[k]][[1]] <- list()
           results[[k]][[1]][1] <- I.labels[i]
           results[[k]][[1]][2] <- I.labels[j]
           sorted.R <- sort(R[, j])
           group <- max(which(sorted.R == sorted.R[minsize]))
           repeat {
              if (N - max(group) < minsize) break
              group <- c(group, max(which(sorted.R == sorted.R[minsize + max(group)])))
           } 
           group <- group[-length(group)]
           summary.matrix <- matrix(nrow = length(group) + 1, ncol = 5 + 2*rvm)
           dimnames(summary.matrix)[[2]] <- c("Group", "Lo","Hi", "N", "n", paste("E(X", i, ")", sep = ""), paste("E(X",j, ")", sep = ""), paste("P(X",i,"=",uv[,1],",X",j,"=",uv[,2],")",sep=""), paste("P(X",i,"=",uv[,2],",X",j,"=",uv[,1],")" ,sep=""))
           summary.matrix[, 1] <- 1:nrow(summary.matrix)
           summary.matrix[, 4] <- c(group, N) - c(0, group)
           group <- c(sorted.R[group], max(sorted.R))
           L <- length(group)
           summary.matrix[, 3] <- group
           summary.matrix[, 2] <- c(min(sorted.R), group[-L] + 1)
           member <- apply(1 - outer(R[, j], group, "<="), 1, sum) + 1
           for (g in 1:L) {
              summary.matrix[g, 6] <- mean(X[member == g, i])
              summary.matrix[g, 7] <- mean(X[member == g, j])
              u <- 0
              for (u in 1:nrow(uv)){ 
                summary.matrix[g,7+u]           <- length(X[member==g & X[,i]==uv[u,1] & X[,j]==uv[u,2],])/J
                summary.matrix[g,7+nrow(uv)+ u] <- length(X[member==g & X[,i]==uv[u,2] & X[,j]==uv[u,1],])/J 
              }  
           }# end g-loop
           summary.matrix[,5] <- apply(summary.matrix[,8:(7+2*nrow(uv))],1,sum)
           results[[k]][[2]] <- summary.matrix
           ac <- violation.matrix[1:(rvm - 1), 1] <- L - 1
           for (g in 1:nrow(uv)) {
              n.uv <- summary.matrix[, 7 + g]
              n.vu <- summary.matrix[, 7 + nrow(uv) + g]
                Ng <- summary.matrix[,4]
                d <- (n.uv - n.vu)/Ng
                d[d <= minvi] <- 0
                vi <- length(d[d > minvi/2])
                sum.vi <- sum(d)
                chi.statistic <- rep(0, L)
                chi.pvalue <- rep(1, L)
                if (any(d > 0)) {
                    for (gg in 1:L) {
                    if (d[gg] > 0) {
                        e <- (n.uv[gg] + n.vu[gg])/2  
                        chi.statistic[gg] <- (n.uv[gg] - e)^2/e + (n.vu[gg] - e)^2/e
                        chi.pvalue[gg] <- 1- pchisq(chi.statistic[gg],1)
                    }#endif
                    }#end gg-loop
                }#endif
                violation.matrix[g, 2:8] <- c(vi,vi/ac, max(d), sum.vi, sum.vi/ac, max(chi.statistic),sum(chi.pvalue < alpha))
              vi.matrix[i,j] <- vi.matrix[i,j] + sum(chi.pvalue < alpha)
           }#end g-loop
           violation.matrix[rvm,c(1,2,5,8)] <- apply(violation.matrix[1:(rvm-1),c(1,2,5,8)],2,sum)
           violation.matrix[rvm,3] <- violation.matrix[rvm,2]/violation.matrix[rvm,1]
           violation.matrix[rvm,6] <- violation.matrix[rvm,5]/violation.matrix[rvm,1]
           violation.matrix[rvm,c(4,7)] <- apply(violation.matrix[1:(rvm-1),c(4,7)],2,max)
           results[[k]][[3]] <- violation.matrix
        }#end j-loop
     }#end i-loop
   }#end if (method IT)


   # METHOD MS-CPM
   if (method=="MSCPM"){
      for (i in 1:(J - 1)) {
         R <- as.matrix(X) %*% (matrix(1, J, J) - diag(J)) - X[,i]
         R[, i] <- 0
         for (j in (i + 1):J) {
            k <- k + 1
            rvm <- (ncat - 1) + 1
            violation.matrix <- matrix(0, nrow = rvm, ncol = 8)
            dimnames(violation.matrix) <- list(c(t(paste("P(X",i, ">=", 1:(ncat - 1), ")  P(X",j, ">=", 1:(ncat - 1), ")", sep = "")), "Total"), 
                                               c("#ac", "#vi", "#vi/#ac", "maxvi", "sum", "sum/#ac","zmax", "#zsig"))
            results[[k]] <- list()
            results[[k]][[1]] <- list()
            results[[k]][[1]][1] <- I.labels[i]
            results[[k]][[1]][2] <- I.labels[j]
            sorted.R <- sort(R[, j])
            group <- max(which(sorted.R == sorted.R[minsize]))
            repeat {
               if (N - max(group) < minsize) break
               group <- c(group, max(which(sorted.R == sorted.R[minsize + max(group)])))
            } 
            group <- group[-length(group)]
            summary.matrix <- matrix(nrow = length(group) + 1, ncol = 6 + 2 * (ncat - 1))
            dimnames(summary.matrix)[[2]] <- c("Group", "Lo","Hi", "N", paste("E(X", i, ")", sep = ""), paste("E(X",j, ")", sep = ""), paste("P(X", i, ">=", 1:(ncat -1), ")", sep = ""), paste("P(X", j, ">=", 1:(ncat -1), ")", sep = ""))
            summary.matrix[, 1] <- 1:nrow(summary.matrix)
            Ng <- summary.matrix[, 4] <- c(group, N) - c(0, group)
            group <- c(sorted.R[group], max(sorted.R))
            L <- length(group)
            summary.matrix[, 3] <- group
            summary.matrix[, 2] <- c(min(sorted.R), group[-L] + 1)
            member <- apply(1 - outer(R[, j], group, "<="), 1, sum) + 1
             
            for (g in 1:L) {
               summary.matrix[g, 5] <- mean(X[member == g, i])
               summary.matrix[g, 6] <- mean(X[member == g, j])
               freqi <- tabulate(X[member == g, i] + 1, ncat)
               freqj <- tabulate(X[member == g, j] + 1, ncat)
               cum.freqi <- rev(cumsum(rev(freqi))/Ng[g])
               cum.freqj <- rev(cumsum(rev(freqj))/Ng[g])
               summary.matrix[g, 7:(5 + ncat)] <- cum.freqi[2:ncat]
               summary.matrix[g, (6 + ncat):(4 + 2 * ncat)] <- cum.freqj[2:ncat]
            }# end g-loop
            results[[k]][[2]] <- summary.matrix

            ac <- violation.matrix[1:(rvm - 1), 1] <- L - 1
            for (g in 1:(ncat - 1)) {
               p.1 <- summary.matrix[, 6 + g]
               p.2 <- summary.matrix[, 5 + ncat + g]
               d <- p.2 - p.1
               d[d <= minvi] <- 0
               vi <- length(d[d > minvi/2])
               sum.vi <- sum(d)
               z.statistic <- rep(0, L)
               if (any(d > 0)) {
                 for (gg in 1:L) {
                   if (d[gg] > 0) {
                     Xgg <- X[member == gg, ]
                     f.01 <- length(which(Xgg[, i]>=g& Xgg[, j] < (g+ncat)))
                     f.10 <- length(which(Xgg[, i]< g& Xgg[, j]>= (g+ncat)))
                     f.k <- min(f.01, f.10)
                     f.n <- f.01 + f.10
                     f.b <- ((2 * f.k + 1 - f.n)^2 - 10 * f.n)/(12 * f.n)
                     z.statistic[gg] <- abs(sqrt(2 * f.k + 2 + f.b) - sqrt(2 * f.n - 2 * f.k + f.b))
                   }#endif
                 }#end gg-loop
               }#endif
               violation.matrix[g, 2:8] <- c(vi,vi/ac, max(d), sum.vi, sum.vi/ac, max(z.statistic),length(z.statistic[abs(z.statistic) > qnorm(1-alpha)]))
               vi.matrix[i,j] <- vi.matrix[i,j] + length(z.statistic[abs(z.statistic) > qnorm(1-alpha)])
            }#end g-loop
            violation.matrix[rvm,c(1,2,5,8)] <- apply(violation.matrix[1:(rvm-1),c(1,2,5,8)],2,sum)
            violation.matrix[rvm,3] <- violation.matrix[rvm,2]/violation.matrix[rvm,1]
            violation.matrix[rvm,6] <- violation.matrix[rvm,5]/violation.matrix[rvm,1]
            violation.matrix[rvm,c(4,7)] <- apply(violation.matrix[1:(rvm-1),c(4,7)],2,max)
            results[[k]][[3]] <- violation.matrix
         }#end j-loop
      }#end i-loop
   }#end if (method MS-CPM)

   # METHOD MIIO
   dich=2
   if (method=="MIIO"){
      for (i in 1:(J-1)){
         R <- as.matrix(X) %*% (matrix(1,J,J) - diag(J)) - X[,i]
         R[,i] <- 0
         for (j in (i+1):J){
            k <- k + 1
            rvm <- 2
            violation.matrix <- matrix(0, nrow = 1, ncol = 8)
            dimnames(violation.matrix) <- list(c(t(paste("E(X",i,")  E(X",j,")", sep = ""))), 
                                               c("#ac", "#vi", "#vi/#ac", "maxvi", "sum", "sum/#ac",ifelse(ncat == dich,"zmax", "tmax"), ifelse(ncat == dich, "#zsig", "#tsig")))
            results[[k]] <- list()
            results[[k]][[1]] <- list()
            results[[k]][[1]][1] <- I.labels[i]
            results[[k]][[1]][2] <- I.labels[j]

            sorted.R <- sort(R[,j])
            group <- max(which(sorted.R==sorted.R[minsize]))
            repeat{
              if(N - max(group) < minsize)break
              group <- c(group,max(which(sorted.R==sorted.R[minsize+max(group)])))
            }
            group <- group[-length(group)]  
            summary.matrix <- matrix(nrow = length(group) + 1, ncol = 8)
            dimnames(summary.matrix)[[2]] <- c("Group", "Lo","Hi", "N", paste("E(X", i, ")", sep = ""), paste("E(X",j, ")", sep = ""),paste("SD(X", i, ")", sep = ""), paste("SD(X",j, ")", sep = ""))
            summary.matrix[, 1] <- 1:nrow(summary.matrix)
            summary.matrix[, 4] <- c(group, N) - c(0, group)
            group <- c(sorted.R[group],max(sorted.R))
            L <- Ll1[j] <- length(group)
            summary.matrix[, 3] <- group
            summary.matrix[, 2] <- c(min(sorted.R), group[-L] + 1)
            member <- apply(1 - outer(R[,j], group, "<="),1,sum) + 1

            for (g in 1:L){
               summary.matrix[g, 5] <- mean(X[member == g, i])
               summary.matrix[g, 6] <- mean(X[member == g, j])
               summary.matrix[g, 7] <- sd(X[member == g, i])
               summary.matrix[g, 8] <- sd(X[member == g, j])
            }# end g-loop
            results[[k]][[2]] <- summary.matrix

            ac <- violation.matrix[1:(rvm - 1), 1] <- L - 1
            E.1 <- summary.matrix[, 5]
            E.2 <- summary.matrix[, 6]
            d <- E.2 - E.1
            d[d <= minvi] <- 0
            vi <- length(d[d > minvi/2])
            sum.vi <- sum(d)
            t.statistic <- rep(0, L)
            t.pvalue <- rep(1, L)
            if (any(d > 0)) {
              for (gg in 1:L) {
                if (d[gg] > 0) {
                   if (ncat > dich){
                      tt <- t.test(X[member==gg,i],X[member==gg,j],alternative="less", paired = TRUE)
                      t.statistic[gg] <- abs(tt$statistic)
                      t.pvalue[gg] <- tt$p.value
                   } else {
                      Xgg <- X[member == gg, ]
                      f.01 <- length(which(Xgg[, i]>=1& Xgg[, j] < 1)) 
                      f.10 <- length(which(Xgg[, i]< 1& Xgg[, j] >= 1))
                      f.k <- min(f.01, f.10)
                      f.n <- f.01 + f.10
                      f.b <- ((2 * f.k + 1 - f.n)^2 - 10 * f.n)/(12 * f.n)
                      t.statistic[gg] <- abs(sqrt(2 * f.k + 2 + f.b) - sqrt(2 * f.n - 2 * f.k + f.b)) # even though it is a z-statistic
                      t.pvalue[gg] <- 1-pnorm(t.statistic[gg])
                   } 
                }#endif
              }#end gg-loop
            }#endif
            violation.matrix[1, 2:8] <- c(vi,vi/ac, max(d), sum.vi, sum.vi/ac, max(t.statistic),sum(t.pvalue < alpha))
            vi.matrix[i,j] <- vi.matrix[i,j] + sum(t.pvalue < alpha)
           results[[k]][[3]] <- violation.matrix
         }# end j-loop
      }# end i-loop
     
     if(!is.null(level.two.var)){
       vi.matrixA <- matrix(0,J,J)
       dimnames(vi.matrixA) <- list(I.labels,I.labels)
       
       Rs <- table(level.two.var)
       Xa <- as.matrix(aggregate(X, by = list(level.two.var), FUN = mean)[, -1])
       Xas <- Xa # / Rs
       Na <- nrow(Xa)
       resultsa <- list()
       #Xa <- Xa[, item.order]
       #Xas <- Xas[, item.order]
       His <- MLcoefH(cbind(level.two.var, X), se = F, nice.output = F)$Hi #coefHTiny(X)$Hi
       g <- 0
       gg <- 0
       h <- 0
       i <- 1
       j <- 2
       k <- 0
       
       for (i in 1:(J-1)){
         
         Ras <- as.matrix(Xa) %*% (matrix(1,J,J) - diag(J)) - Xa[,i]
         Ras <- round(Ras, 4)
         Ras[, i] <- 0
         
         Ra <- Ras[rep(1:nrow(Ras), Rs), ]
         Ra[,i] <- 0
         
         for (j in (i+1):J){
           k <- k + 1
           rvm <- 2
           violation.matrix <- matrix(0, nrow = 1, ncol = 8)
           dimnames(violation.matrix) <- list(c(t(paste("E(X",i,")  E(X",j,")", sep = ""))), 
                                              c("#ac", "#vi", "#vi/#ac", "maxvi", "sum", "sum/#ac",ifelse(ncat == dich,"zmax", "tmax"), ifelse(ncat == dich, "#zsig", "#tsig")))
           resultsa[[k]] <- list()
           resultsa[[k]][[1]] <- list()
           resultsa[[k]][[1]][1] <- I.labels[i]
           resultsa[[k]][[1]][2] <- I.labels[j]
           
           ## Aggregated
           sorted.Ra <- sort(round(Ra[,j], 4))
           minsizeA <- floor(nrow(X) / (Ll1[j] + 1.75))
           group <- max(which(sorted.Ra == sorted.Ra[minsizeA]))
           repeat{
             if(N - max(group) < minsizeA)break
             group <- c(group,max(which(sorted.Ra==sorted.Ra[minsizeA+max(group)])))
           }
           group <- group[-length(group)]  
           summary.matrix <- matrix(nrow = length(group) + 1, ncol = 8)
           dimnames(summary.matrix)[[2]] <- c("Group", "Lo","Hi", "N", 
                                              paste("E(X", i, ")", sep = ""), 
                                              paste("E(X",j, ")", sep = ""),
                                              paste("SD(X", i, ")", sep = ""), 
                                              paste("SD(X",j, ")", sep = ""))
           summary.matrix[, 1] <- 1:nrow(summary.matrix)
           group <- c(sorted.Ra[group],max(sorted.Ra))
           L <- length(group)
           summary.matrix[, 3] <- group #/ Rs
           summary.matrix[, 2] <- c(min(sorted.Ra), group[-L])
           member <- apply(1 - outer(Ras[,j], group, "<="),1,sum) + 1
           summary.matrix[, 4] <- table(member)#(c(group, N) - c(0, group)) / Rs
           #member <- apply(1 - outer(Ra[,j], group, "<="),1,sum) + 1
           
           for (g in 1:L){
             summary.matrix[g, 5] <- mean(Xas[member == g, i]) #s
             summary.matrix[g, 6] <- mean(Xas[member == g, j]) #s
             summary.matrix[g, 7] <- sd(Xas[member == g, i]) #s
             summary.matrix[g, 8] <- sd(Xas[member == g, j]) #s
           }# end g-loop
           resultsa[[k]][[2]] <- summary.matrix
           
           ac <- violation.matrix[1:(rvm - 1), 1] <- L - 1
           E.1 <- summary.matrix[, 5]
           E.2 <- summary.matrix[, 6]
           d <- E.2 - E.1
           d[d <= minvi] <- 0
           vi <- length(d[d > minvi/2])
           sum.vi <- sum(d)
           t.statistic <- rep(0, L)
           t.pvalue <- rep(1, L)
           if (any(d > 0)) {
             for (gg in 1:L) {
               if (d[gg] > 0) {
                 tt <- t.test(Xas[member==gg,i],Xas[member==gg,j],alternative="less", paired = TRUE)
                 t.statistic[gg] <- abs(tt$statistic)
                 t.pvalue[gg] <- tt$p.value 
               }#endif
             }#end gg-loop
           }#endif
           violation.matrix[1, 2:8] <- c(vi,vi/ac, max(d), sum.vi, sum.vi/ac, max(t.statistic),sum(t.pvalue < alpha))
           vi.matrixA[i,j] <- vi.matrixA[i,j] + sum(t.pvalue < alpha)
           resultsa[[k]][[3]] <- violation.matrix
         }# end j-loop
       }# end i-loop
       #vi.matrixA <- vi.matrix
       vi.matrixA <- vi.matrixA + t(vi.matrixA)
       VIA <- matrix(apply(sign(vi.matrixA),1,sum)) 
       dimnames(VIA) <- list(I.labels,paste("step",1:ncol(VIA)))
     } # end if(!is.null(level.two.var))
     
   }# end if(method="MIIO") 

   vi.matrix <- vi.matrix + t(vi.matrix)
   vi.matrix.a <- vi.matrix
   VI <- matrix(apply(sign(vi.matrix),1,sum)) 
   items.removed <- NULL
   if (item.selection){
     repeat{
       nvi <- apply(sign(vi.matrix.a),2,sum)
       if (sum(nvi)==0) break
       maxvi <- which(nvi ==max(nvi))
       if(length(maxvi) > 1){ 
          H <- rep(0,length(maxvi))
          for (i in 1:length(maxvi)) H[i] <- coefHTiny(X[,-c(items.removed,maxvi[i])])$H + rnorm(1,0,1e-5)
          maxvi <- maxvi[which(H==max(H))[1]]
       }# end if 
       vi.matrix.a[maxvi,] <- vi.matrix.a[,maxvi] <- 0
       items.removed <- c(items.removed,maxvi)
       VI. <- matrix(apply(sign(vi.matrix.a),1,sum)) 
       VI.[items.removed,] <- NA
       VI <- cbind(VI,VI.)
     }
   }# end if(item.selection)  
   dimnames(VI) <- list(I.labels,paste("step",1:ncol(VI)))
   if(verbose) print(VI)
   
   # Computation of HT

   if (length(items.removed) > 0) X <- X[,-items.removed]
   HT <- coefHT(X)
   
   iio.list <- list(results = results, violations = VI, items.removed = items.removed, Hi = Hi, HT = HT, method = method, item.mean = item.mean, m = ncat-1)
   class(iio.list) <- "iio.class"
   
   if(!is.null(level.two.var)) {
     HTA <- coefHT(Xas)
     iio.lista <- list(results = resultsa, violations = VIA, items.removed = items.removed, Hi = His[, 2], HT = HTA, method = method, item.mean = item.mean, m = ncat-1)
     class(iio.lista) <- "iio.class"
     
     iio.list <- list("Level.1" = iio.list, "Level.2" = iio.lista)
     class(iio.list) <- "iio.class"
   }
   
   return(iio.list)
}

coefHT <- function(Y, largeN = 1000){
   eq.var <- apply(Y, 1, sd)
   Y <- Y[eq.var > 0, ]
   N <- nrow(Y)
   tY <- t(Y)
   tYm <- apply(tY, 2, sort)

   G <- ceiling(N / largeN)
   sum.S <- 0
   sum.Smax <- 0
   repeat{
     if (N %% largeN != 1) break
     largeN <- largeN -1
   }
   for (i in 1 : G) for (j in 1 : G){
     ni <- ifelse(i < G, largeN, N %% largeN)
     nj <- ifelse(j < G, largeN, N %% largeN)
     S <- var(tY[, (i-1) * largeN + 1 : ni], tY[, (j-1) * largeN + 1 : nj])
     Smax <- var(tYm[, (i-1) * largeN + 1 : ni], tYm[, (j-1) * largeN + 1 : nj])
     if (i == j) diag(S) <- diag(Smax) <- 0
     sum.S <- sum.S + sum(S)
     sum.Smax <- sum.Smax + sum(Smax)
   }
   return(sum.S / sum.Smax)
}

# DECREPIT New function above 
#coefHT <- function(Y, largeN = 1000){
#    eij <- function(x) {
#        tab <- tabulate(x[ ,1] + 2 * x[ ,2] + 1L, 4)
#        e <- min(tab[3], tab[2])
#        e0 <- ((tab[1] + e) * (e + tab[4])) / nrow(x)
#        return(c(e, e0))
#    }
#
#   eq.var <- apply(Y, 1, sd)
#   Y <- Y[eq.var > 0, ]
#   N <- nrow(Y)
#   tY <- t(Y)
#   if (N < largeN){ 
#      S    <- var(tY)
#      Smax <- var(apply(tY, 2, sort))
#      diag(S) <- diag(Smax) <- 0
#      HT <- sum(S)/sum(Smax)
#   } else {  
#      e1 <- 0L
#      e0 <- 0L
#      for (i in 1 : (ncol(tY) - 1)){
#         for (j in (i + 1) : ncol(tY)){
#            e <- eij(tY[, c(i,j)])
#            e1 <- e1 + e[1]
#            e0 <- e0 + e[2]
#         }
#         cat("\r", "Large sample size, computation may take several minutes. Progress: ", round((i * 100) / (ncol(tY) - 1), 1), "%      ", sep = "") 
#         flush.console()  # Make visible on RGui
#      }
#      cat("\r", "Large sample size, computation may take several minutes. Progress: DONE    ", sep = "") 
#      cat("\n")
#      HT <- 1 - e1/e0
#   }
#   return(HT)
#}
