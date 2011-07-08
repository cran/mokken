"summary.iio.class" <-
function(object, ...){
   results <- object$results
   violations <- object$violations
   items.removed <- object$items.removed
   HT <- object$HT
   method <- object$method
   item.mean <- object$item.mean
   J <- nrow(violations)
   I.labels <- dimnames(violations)[[1]]
   summary.matrix <- matrix(0,nrow=J,ncol=9)
   summary.matrix[,1] <- round(item.mean,2)
   tmax <- rep(0,J)
   if(method=="MIIO") {
     dimnames(summary.matrix) <- list(I.labels,c("mean","#ac","#vi","#vi/#ac","maxvi","sum","sum/#ac","tmax","#tsig"))
   }
   if(method=="MSCPM") {
     dimnames(summary.matrix) <- list(I.labels,c("mean","#ac","#vi","#vi/#ac","maxvi","sum","sum/#ac","zmax","#zsig"))
   }
   if(method=="IT") {
     dimnames(summary.matrix) <- list(I.labels,c("mean","#ac","#vi","#vi/#ac","maxvi","sum","sum/#ac","xmax","#xsig"))
   }
   k <- 0; i <- 1; j <- 2
   for(i in 1:(J-1)){for(j in (i+1):J){
     k <- k+1
     input <- results[[k]][[3]]
     rvm <- nrow(input)
     tmax[i] <- max(tmax[i],input[rvm,7])     
     tmax[j] <- max(tmax[j],input[rvm,7])     
     summary.matrix[i,c(2,3,6,9)] <- summary.matrix[i,c(2,3,6,9)] + input[rvm,c(1,2,5,8)]
     summary.matrix[j,c(2,3,6,9)] <- summary.matrix[j,c(2,3,6,9)] + input[rvm,c(1,2,5,8)]
     summary.matrix[i,c(5)] <- max(summary.matrix[i,c(5)], input[rvm,c(4)])
     summary.matrix[j,c(5)] <- max(summary.matrix[j,c(5)], input[rvm,c(4)])
     summary.matrix[i,c(8)] <- max(summary.matrix[i,c(8)], input[rvm,c(7)])
     summary.matrix[j,c(8)] <- max(summary.matrix[j,c(8)], input[rvm,c(7)])
   }}
   summary.matrix[,4] <- summary.matrix[,3]/summary.matrix[,2]
   summary.matrix[,7] <- summary.matrix[,6]/summary.matrix[,2]
   summary.matrix[,8] <- tmax
   
   return(list(method = method, item.summary = round(summary.matrix,2), backward.selection= violations, HT=HT))
}
