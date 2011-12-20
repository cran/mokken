"plot.iio.class" <-
function(x, item.pairs = all.pairs, ask = TRUE, ...){
  J <- length(x$item.mean)
  m <- x$m
  max.item.pairs <- J*(J-1)/2
  all.pairs <- 1:max.item.pairs
  j <- 0; i <- 0
  results <- x$results
  if (ask==TRUE) par("ask"=TRUE) else par("ask"=FALSE)
  for (j in item.pairs){
    plot.matrix <- results[[j]][[2]]
    x.labels <- paste(plot.matrix[,2],"-",plot.matrix[,3],sep="")
    plot(plot.matrix[,1],plot.matrix[,5],
      ylim=c(0,m),
      xaxt = 'n',
      xlab = "Rest score group",
      ylab = "Item rest function",
      type = "n", lwd=3)
    lines(plot.matrix[,1],plot.matrix[,5], lwd=2, lty=1)
    lines(plot.matrix[,1],plot.matrix[,6], lwd=2, lty=3)
    title(paste(results[[j]][[1]][1],"(solid)",results[[j]][[1]][2],"(dashed)"))
    axis(1, at=1:nrow(plot.matrix),labels=x.labels)
  }
 invisible()
}
