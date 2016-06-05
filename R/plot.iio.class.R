"plot.iio.class" <-
function(x, item.pairs = all.pairs, ci = TRUE, alpha = .05, color = c("black", "blue"), transparancy = 20, ask = TRUE, ...){

up.lo.bound.mean <- function(n,alpha=.05){
   n[n < 1e-10] <- 1e-10
   n <- matrix(n)
   p <- length(n)
   scores <- rep(c(0:(p-1)),n)
   m <- mean(scores)
   ase <- sd(scores)/sqrt(sum(n))
   z <- qnorm(1 - alpha/2)
   matrix(c(m - z * ase, m + z * ase),1,2)
}
  
  def.par <- par(no.readonly = TRUE)
  J <- length(x$item.mean)
  max.item.pairs <- J*(J-1)/2
  all.pairs <- 1:max.item.pairs
  results <- x$results
  m <- x$m
  if (ask==TRUE) par("ask"=TRUE) else par("ask"=FALSE)
  i <- 0; j <- 0
  c1 <- as.numeric(col2rgb(color[1]))
  c2 <- as.numeric(col2rgb(color[2]))
  colorCi   = c(rgb(c1[1], c1[2], c1[3], transparancy, maxColorValue = 255), rgb(c2[1], c2[2], c2[3], transparancy, maxColorValue = 255))

  for (j in item.pairs){
    plot.matrix <- results[[j]][[2]]
    x.labels <- paste(plot.matrix[,2],"-",plot.matrix[,3],sep="")
    z <- qnorm(1 - alpha/2)
    mi1 <- plot.matrix[,5]
    mi2 <- plot.matrix[,6]
    if(ci){ 
        lo1 <-  plot.matrix[,5] - z * plot.matrix[,7]/sqrt(plot.matrix[,4])
        lo2 <-  plot.matrix[,6] - z * plot.matrix[,8]/sqrt(plot.matrix[,4])
        up1 <-  plot.matrix[,5] + z * plot.matrix[,7]/sqrt(plot.matrix[,4])
        up2 <-  plot.matrix[,6] + z * plot.matrix[,8]/sqrt(plot.matrix[,4])
    }
    plot(plot.matrix[,1],mi1,
     ylim=c(0,m),
     xaxt = 'n',
     xlab = "Rest score group",
     ylab = "Item response functions",
     type = "n")
    title(paste(results[[j]][[1]][1],"(solid)",results[[j]][[1]][2],"(dashed)"))
    axis(1, at=1:nrow(plot.matrix),labels=x.labels)
    if(ci){
       polygon(c((1:length(up1))[!is.na(up1)],rev((1:length(lo1))[!is.na(lo1)])),c(up1[!is.na(up1)],rev(lo1[!is.na(lo1)])),col = colorCi[1], border=NA)
       polygon(c((1:length(up2))[!is.na(up2)],rev((1:length(lo2))[!is.na(lo2)])),c(up2[!is.na(up2)],rev(lo2[!is.na(lo2)])),col = colorCi[2], border=NA)
    }   
    lines(plot.matrix[,1],mi1, lwd=4, lty=1, col = color[1])
    lines(plot.matrix[,1],mi2, lwd=4, lty=2, col = color[2])
  }
 invisible()
 par(def.par)
}
