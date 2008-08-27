"plot.monotonicity.class" <-
function(x, items = all, ...){
  results <- x$results
  m <- x$m
  all <- 1:length(x$I.labels)
  i <- 0; j <- 0
  for (j in items){
    plot.matrix <- results[[j]][[2]]
    x.labels <- paste(plot.matrix[,2],"-",plot.matrix[,3],sep="")
    par("ask"=TRUE)
    plot(plot.matrix[,1],plot.matrix[,m+5]/(m-1),
      ylim=c(0,1),
      xaxt = 'n',
      xlab = "Rest score group",
      ylab = "Item rest function",
      type = "l", 
      lwd=3)
    title(results[[j]][[1]])
    axis(1, at=1:nrow(plot.matrix),labels=x.labels)
    for(i in 2:m){
     lines(plot.matrix[,1],plot.matrix[,(m+4+i)], col=4, lwd=2,lty=3)
    }
  }
 invisible()
}


