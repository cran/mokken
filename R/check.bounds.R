"check.bounds" <- function(X, quant = .90, lower = TRUE, upper = FALSE){
   X <- check.data(X)
   J <- ncol(X)
   if(J < 3) stop("At least three items are required to compute Ellis' upper and lower bounds.")
   R <- cor(X)
   R[abs(R) < 1e-10] <- 1e-10

   if (lower){
     RLB1 <- RLB2 <- matrix(0, J, J)
     for (i in 1 : (J - 1)) for (j in (i + 1) : J){
        arg <- (R[i, ] * R[, j])[-c(i, j)]
        RLB1[j, i] <- RLB1[i, j] <- max(arg)
        RLB2[j, i] <- RLB2[i, j] <- quantile(arg, quant)
      }
      L <- list(L1rij = RLB1, L2rij = RLB2)
   }   
   if (upper){
      RUB1 <- RUB2 <- matrix(1, J, J)
      for (i in 1 : (J - 1)) for (j in (i + 1) : J){
         arg <- pmin(c((R[i, ] / R[, j])[-c(i, j)]), c((R[j, ] / R[, i])[-c(i, j)]))
         RUB1[j, i] <- RUB1[i, j] <- min(arg) 
         RUB2[j, i] <- RUB2[i, j] <- quantile(arg, 1 - quant)
      }
      U <- list(U1rij = RUB1, U2rij = RUB2)
   }
   if (upper  &  lower) output <- list(UpperBounds = U, LowerBounds=L) 
   if (upper  & !lower) output <- list(UpperBounds = U) 
   if (!upper &  lower) output <- list(LowerBounds = L) 
   if (!upper & !lower) output <- NA 
   return(output)
}
