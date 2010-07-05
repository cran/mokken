"summary.monotonicity.class" <-
function(object, ...){
   results <- object$results
   I.labels <- object$I.labels
   Hi <- object$Hi
   m <- object$m
   J <- length(results)
   j <- 0
   summary.matrix <- matrix(nrow=J,ncol=9)
   dimnames(summary.matrix) <- list(I.labels,c("ItemH","#ac","#vi","#vi/#ac","maxvi","sum","sum/#ac","zmax","#zsig"))
   summary.matrix[,1] <- Hi
   for (j in 1:J) summary.matrix[j,2:9] <- results[[j]][[3]][m,c(1:7,10)]
   return(round(summary.matrix,2))
}

