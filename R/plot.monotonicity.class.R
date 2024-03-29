"plot.monotonicity.class" <-
function(x, items = all.items, curves = "both", ci = TRUE, alpha = .05, color = "black", transparancy = 20, ask = TRUE, ...){
up.lo.bound.mean <- function(n, alpha = .05){
   n[n < 1e-10] <- 1e-10
   n <- matrix(n)
   p <- length(n)
   scores <- rep(c(0:(p-1)),n)
   m <- mean(scores)
   ase <- sd(scores)/sqrt(sum(n))
   z <- qnorm(1 - alpha/2)
   matrix(c(m - z * ase, m + z * ase),1,2)
}

if(length(x) == 2) {
  x2l <- x
  for(lvl in 1:2){
    x <- x2l[[lvl]]
    def.par <- par(no.readonly = TRUE)
    results <- x$results
    m <- x$m
    all.items <- 1:length(x$I.labels)
    if (ask==TRUE) par("ask"=TRUE) else par("ask"=FALSE)
    if (curves == "both") layout(matrix(c(1,2),1,2))
    i <- 0; j <- 0
    c1 <- as.numeric(col2rgb(color))
    colorCi = rgb(c1[1], c1[2], c1[3], alpha = transparancy, maxColorValue = 255)
    
    
    for (j in items){
      plot.matrix <- results[[j]][[2]]
      x.labels <- paste(plot.matrix[,2],"-",plot.matrix[,3],sep="")
      if(curves=="both" | curves=="ISRF"){
        est <- t(plot.matrix[,(m+4+2):(m+4+m)])
        if(ci){ 
          n = plot.matrix[,4]
          n[n < 1e-10] = 1e-10
          if (m > 2) se  = t(apply(est, 1, function(x) sqrt((x - x^2)/n)))
          if (m==2) se = sqrt((est - est^2)/n)
          lo = (est-qnorm(1-alpha/2) * se)
          lo[lo < 0] <- 0
          up = (est+qnorm(1-alpha/2) * se)
          up[up > 1] <- 1
        }   
        plot(plot.matrix[,1],est[1,],
             ylim=c(0,1),
             xaxt = 'n',
             xlab = "Rest score group",
             ylab = "Item step response function",
             type = "n")
        title(paste0("Level ", lvl, ": ", results[[j]][[1]]))
        axis(1, at=1:nrow(plot.matrix),labels=x.labels)
        if(m==2){
          if(ci) polygon(c((1:length(up))[!is.na(up)],rev((1:length(lo))[!is.na(lo)])),c(up[!is.na(up)],rev(lo[!is.na(lo)])), col=colorCi, border=NA)
          lines(plot.matrix[,1],est, lwd=4, col = color)
        } 
        if(m>2){
          if(ci) for(i in 1:(m-1)) polygon(c((1:length(up[i,]))[!is.na(up[i,])],rev((1:length(lo[i,]))[!is.na(lo[i,])])),c(up[i,!is.na(up[i,])],rev(lo[i,!is.na(lo[i,])])),col=colorCi, border=NA)
          for(i in 1:(m-1)) lines(plot.matrix[,1],est[i,], lwd=3, col = color)
        }
      }  
      if(curves=="both" | curves=="IRF"){
        est <- t(plot.matrix[,m+5])
        if(ci){
          up.lo <- apply(plot.matrix[,1:m+4],1,up.lo.bound.mean, alpha)
          lo <- up.lo[1,]
          up <- up.lo[2,]
        }   
        plot(plot.matrix[,1],est,
             ylim = c(0, m - 1),
             xaxt = 'n',
             xlab = "Rest score group",
             ylab = "Item response function",
             type = "n", 
             lwd=3)
        title(paste0("Level ", lvl, ": ", results[[j]][[1]]))
        axis(1, at=1:nrow(plot.matrix),labels=x.labels)
        if(ci) polygon(c((1:length(up))[!is.na(up)],rev((1:length(lo))[!is.na(lo)])),c(up[!is.na(up)],rev(lo[!is.na(lo)])),col = colorCi, border=NA)
        lines(plot.matrix[,1],est, lwd=4, col = color)
      }
    }
    invisible()
    par(def.par)
  }
} else {
  def.par <- par(no.readonly = TRUE)
  results <- x$results
  m <- x$m
  all.items <- 1:length(x$I.labels)
  if (ask==TRUE) par("ask"=TRUE) else par("ask"=FALSE)
  if (curves == "both") layout(matrix(c(1,2),1,2))
  i <- 0; j <- 0
  c1 <- as.numeric(col2rgb(color))
  colorCi = rgb(c1[1], c1[2], c1[3], alpha = transparancy, maxColorValue = 255)
  
  
  for (j in items){
    plot.matrix <- results[[j]][[2]]
    x.labels <- paste(plot.matrix[,2],"-",plot.matrix[,3],sep="")
    if(curves=="both" | curves=="ISRF"){
      est <- t(plot.matrix[,(m+4+2):(m+4+m)])
      if(ci){ 
        n = plot.matrix[,4]
        n[n < 1e-10] = 1e-10
        if (m > 2) se  = t(apply(est, 1, function(x) sqrt((x - x^2)/n)))
        if (m==2) se = sqrt((est - est^2)/n)
        lo = (est-qnorm(1-alpha/2) * se)
        lo[lo < 0] <- 0
        up = (est+qnorm(1-alpha/2) * se)
        up[up > 1] <- 1
      }   
      plot(plot.matrix[,1],est[1,],
           ylim=c(0,1),
           xaxt = 'n',
           xlab = "Rest score group",
           ylab = "Item step response function",
           type = "n")
      title(results[[j]][[1]])
      axis(1, at=1:nrow(plot.matrix),labels=x.labels)
      if(m==2){
        if(ci) polygon(c((1:length(up))[!is.na(up)],rev((1:length(lo))[!is.na(lo)])),c(up[!is.na(up)],rev(lo[!is.na(lo)])), col=colorCi, border=NA)
        lines(plot.matrix[,1],est, lwd=4, col = color)
      } 
      if(m>2){
        if(ci) for(i in 1:(m-1)) polygon(c((1:length(up[i,]))[!is.na(up[i,])],rev((1:length(lo[i,]))[!is.na(lo[i,])])),c(up[i,!is.na(up[i,])],rev(lo[i,!is.na(lo[i,])])),col=colorCi, border=NA)
        for(i in 1:(m-1)) lines(plot.matrix[,1],est[i,], lwd=3, col = color)
      }
    }  
    if(curves=="both" | curves=="IRF"){
      est <- t(plot.matrix[,m+5])
      if(ci){
        up.lo <- apply(plot.matrix[,1:m+4],1,up.lo.bound.mean, alpha)
        lo <- up.lo[1,]
        up <- up.lo[2,]
      }   
      plot(plot.matrix[,1],est,
           ylim = c(0, m - 1),
           xaxt = 'n',
           xlab = "Rest score group",
           ylab = "Item response function",
           type = "n", 
           lwd=3)
      title(results[[j]][[1]])
      axis(1, at=1:nrow(plot.matrix),labels=x.labels)
      if(ci) polygon(c((1:length(up))[!is.na(up)],rev((1:length(lo))[!is.na(lo)])),c(up[!is.na(up)],rev(lo[!is.na(lo)])),col = colorCi, border=NA)
      lines(plot.matrix[,1],est, lwd=4, col = color)
    }
  }
  invisible()
  par(def.par)
}

  
}
