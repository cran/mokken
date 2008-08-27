"plot.restscore.class" <-
function(x, item.pairs = all, ...){
  J <- length(x$Hi)
  max.item.pairs <- J*(J-1)/2
  all <- 1:max.item.pairs
  j <- 0; i <- 0
  results <- x$results
  m <- x$m
  for (j in item.pairs){
    plot.matrix <- results[[j]][[2]]
    x.labels <- paste(plot.matrix[,2],"-",plot.matrix[,3],sep="")
    par("ask"=TRUE)
    plot(plot.matrix[,1],plot.matrix[,5]/m,
      ylim=c(0,1),
      xaxt = 'n',
      xlab = "Rest score group",
      ylab = "Item rest function",
      type = "n", lwd=3)
    lines(plot.matrix[,1],plot.matrix[,5]/m, lwd=5, lty=1)
    lines(plot.matrix[,1],plot.matrix[,6]/m, lwd=5, lty=3)
    title(paste(results[[j]][[1]][1],"(solid)",results[[j]][[1]][2],"(dashed)"))
    axis(1, at=1:nrow(plot.matrix),labels=x.labels)
    for(i in 1:(m-1)){
     lines(plot.matrix[,1],plot.matrix[,(6+i)], col=4, lwd=2)
     lines(plot.matrix[,1],plot.matrix[,(6+(m-1)+i)], col=4, lwd=2, lty=3)
    }
  }
 invisible()
}

