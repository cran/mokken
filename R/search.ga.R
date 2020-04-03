"search.ga" <- function(X, popsize, maxgens, alpha, critval, pxover, pmutation){
             
     nitem <- ncol(X)
        
     if(nitem <= 10)      iter <- maxgens
     else if(nitem <= 20) iter <- 5000
     else if(nitem > 20)  iter <- round(4000/(nitem/20),0)
     
     if(maxgens < iter) iter <- maxgens
     
     npers <- nrow(X)
     variance <- var(X)
     max.variance <- var(apply(X, 2, sort))
     SijMatrix <- outer(apply(X, 2, var), apply(X, 2, var), "*")

     ans <- runGeneticAlgorithm(
             POPSIZE = as.integer(popsize), 
             NPERS = as.integer(npers),
             MAXGENS = as.integer(iter),
             PXOVER = as.double(pxover),
             PMUTATION = as.double(pmutation), 
             critval = as.double(critval),
             alpha = as.double(alpha),
             NITEM = as.integer(nitem),
             ITER = as.integer(iter),
             VAR = variance,
             MAXVAR = max.variance, 
             SijMatrix = SijMatrix
     )
     
     InSet <- as.matrix(c(matrix(c(ans),popsize+2,nitem,byrow=T)[popsize+1,]))
     return(InSet)
}
