"plot.pmatrix.class" <-
function(x, items = all, pmatrix = "both", ...){
  all <- 1:max(x$I.item)
  m <- length(x$I.item)/max(x$I.item)
  j <- 1; i <- 1
  #
  if (pmatrix == "both" || pmatrix == "ppp"){
    I.item <- x$I.item
    I.step <- x$I.step
    plot.matrix <- x$Ppp
    for (j in items){
       plot.matrix.j <- plot.matrix[I.item==j,I.item!=j]
       if(!is.matrix(plot.matrix.j)) plot.matrix.j <- t(as.matrix(plot.matrix.j))
       I.step.j <- I.step[I.item!=j]
       x.axis <- length(I.step.j)
       par("ask"=TRUE)
       plot(1:x.axis,plot.matrix.j[1,],
         ylim=c(0,1),
         xlim=c(1,x.axis),
         xaxt = 'n',
         xlab = "ordered item steps",
         ylab = paste("P(X",j," >= x, item step)",sep=""),
         type = "n", lwd=3)
       title(paste("P(++) matrix: ", x$I.labels[[j]]))
       if (x.axis < 10) axis(1, at=1:x.axis,labels=I.step.j) else axis(1, at=1:x.axis,labels=rep("",x.axis))
       for(i in 1:m) lines(1:x.axis,plot.matrix.j[i,], col=4, lwd=2)
     }
  }
  if (pmatrix == "both" || pmatrix == "pmm"){
    I.item <- x$I.item
    I.step <- x$I.step
    plot.matrix <- x$Pmm
    for (j in items){
       plot.matrix.j <- plot.matrix[I.item==j,I.item!=j]
       I.step.j <- I.step[I.item!=j]
       x.axis <- length(I.step.j)
       plot(1:x.axis,plot.matrix.j[1,],
         ylim=c(0,1),
         xlim=c(1,x.axis),
         xaxt = 'n',
         xlab = "ordered item steps",
         ylab = paste("P(X",j," < x| item step)",sep=""),
         type = "n", lwd=3)
       title(paste("P(--) matrix: ", x$I.labels[[j]]))
       if (x.axis < 10) axis(1, at=1:x.axis,labels=I.step.j) else axis(1, at=1:x.axis,labels=rep("",x.axis))
       for(i in 1:m) lines(1:x.axis,plot.matrix.j[i,], col=4, lwd=2)
     }
  }
 invisible()
}

