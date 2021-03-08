"twoway" <- function(X, nCompletedDataSets = 1, minX = defaultMinX, maxX = defaultMaxX, seed = FALSE){

   if (is.numeric(seed)) set.seed(seed)
   if (!is.matrix(X) && !is.data.frame(X)) stop("X must be a matrix or a data.frame")
   Xm <- as.matrix(X)
   if (mode(Xm)!="numeric") stop("X must be numeric")
   if (any(Xm < 0, na.rm = TRUE)) stop("All scores must be nonnegative")
   if (any(Xm %% 1 !=0, na.rm = TRUE)) stop("All scores must be integers")
   defaultMaxX <- max(Xm, na.rm = TRUE)
   defaultMinX <- min(Xm, na.rm = TRUE)
   if (any(Xm < minX, na.rm = TRUE)) stop("All scores must be greater than or equal to minX")
   if (any(Xm > maxX, na.rm = TRUE)) stop("All scores must be less than or equal to maxX")
   N <- nrow(Xm)
   J <- ncol(Xm)

   # Missings
   M <- matrix(FALSE, N, J)
   M[is.na(Xm)] <- TRUE

   # No missings
   noMissing <- !any(M)
   if(noMissing){
      warning("X does not contain missing values")
      return(X)
   }

   # All missings
   if (any(apply(M, 1, sum) == J)) stop("At least one row has no observed scores.") 
   if (any(apply(M, 2, sum) == N)) stop("At least one column has no observed scores.") 
   
   # TWO WAY
   IM <- matrix(apply(Xm, 2, mean, na.rm = TRUE), nrow = N, ncol = J, byrow = TRUE)
   PM <- matrix(apply(Xm, 1, mean, na.rm = TRUE), nrow = N, ncol = J, byrow = FALSE)
   OM <- matrix(mean(Xm, na.rm = TRUE), nrow = N, ncol = J)
   TW <- PM + IM - OM
   Xc <- list()
   for (i in 1 : nCompletedDataSets){ 
      E <- matrix(rnorm(N * J, 0, var((Xm - TW)[!M])), N, J)
      Xtmp <- round(TW + E)
      Xtmp[Xtmp < minX] <- minX
      Xtmp[Xtmp > maxX] <- maxX
      Xtmp[!is.na(Xm)] <- Xm[!is.na(Xm)]
      dimnames(Xtmp) <- dimnames(Xm)
      Xc[[i]] <- Xtmp      
   }
   if (nCompletedDataSets == 1) Xc <- Xc[[1]]
   return(Xc)
}
