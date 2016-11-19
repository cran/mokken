check.norms <- function(y, nice.output = TRUE){

direct.sum <- function (...){
     p.tr = 0;p.ll = 0;
     matlist = list(...);
     nmat = length(matlist);
     m1 = matlist[[1]];
     matlist = if(nmat==1 && is.list(m1)) m1 else matlist # check if list of matrices is given and amend accordingly
     nmat = length(matlist);                              # ,,
     m1 = matlist[[1]];                                   # ,,
     if(nmat==1) return(m1);
     for(i in 2:nmat){
        m2 = matlist[[i]];
        topleft <- m1
        topright <- matrix(p.tr, nrow(m1), ncol(m2))
        colnames(topright) <- colnames(m2)
        lowleft <- matrix(p.ll, nrow(m2), ncol(m1))
        lowright <- m2
        m1 = rbind(cbind(topleft, topright), cbind(lowleft, lowright))
     }
     return(m1)
}
R <- sort(unique(y))
K <- length(R)

n <-matrix(0,K,1)
for (i in 1:K){
n[i] <- sum(y==R[i])
}

N <- length(y)
if (min(y)<=0) {
    adjust <- abs(min(y) - 1)
    y <- y + adjust
    R <- R + adjust
} else {
    adjust <- 0
}

A1 <- rbind(R,1)
A2 <- matrix(c(1,-1),1,2)

meany <- exp(A2%*%log(A1%*%n)) -adjust


g0 <- n
g1 <- log(A1%*%n)

G0 <- diag(K)
G1 <- solve(diag(c(A1%*%g0),length(A1%*%g0),length(A1%*%g0)))%*%A1%*%G0
G2 <- diag(c(exp(A2%*%g1)),length(exp(A2%*%g1)),length(exp(A2%*%g1)))%*%A2%*%G1
Gmean <-G2
Vmean <- Gmean%*%diag(c(n),length(n),length(n))%*%t(Gmean)-Gmean%*%(n%*%t(n)/N)%*%t(Gmean)
Se.mean <- sqrt(Vmean)

lomeany <- meany - 1.96*Se.mean
upmeany <- meany + 1.96*Se.mean

A1 <- rbind(R,R^2,1,1)
A2 <- diag(c(2,1,1,-1))
A2[4,3] <- 1
A2[1,3] <- -1
A3 <- matrix(c(-1,0,1,0,0,1,0,-1),2,4)
A4 <- matrix(c(.5,-0.5),1,2)

sdy <- exp(A4 %*% log (A3 %*% exp(A2 %*% log(A1 %*% n))))

g0 <- n
g1 <- log(A1%*%n)
g2 <- exp(A2%*%log(A1%*%n))
g3 <- log (A3 %*% exp(A2 %*% log(A1 %*% n)))

G0 <- diag(K)
G1 <- solve(diag(c(A1%*%g0),length(A1%*%g0),length(A1%*%g0)))%*%A1%*%G0
G2 <- diag(c(exp(A2%*%g1)),length(exp(A2%*%g1)),length(exp(A2%*%g1)))%*%A2%*%G1
G3 <- solve(diag(c(A3%*%g2),length(A3%*%g2),length(A3%*%g2)))%*%A3%*%G2
G4 <- diag(c(exp(A4%*%g3)),length(exp(A4%*%g3)),length(exp(A4%*%g3)))%*%A4%*%G3

Gsd <-G4

Vsd <- Gsd%*%diag(c(n),length(n),length(n))%*%t(Gsd)-Gsd%*%(n%*%t(n)/N)%*%t(Gsd)
Se.sd <- sqrt(Vsd)

losdy <- sdy - 1.96*Se.sd
upsdy <- sdy + 1.96*Se.sd

A1 <- t(cbind(diag(K),rep(1,K),rep(1,K)))
A2 <- direct.sum(diag(K), t(c(1,-1)))
A3 <- direct.sum(rbind(R, R^2, rep(1, K), rep(1, K)), matrix(R))
A4 <- direct.sum(matrix(c(1,2,0,0,0,0,0,1,0,0,0,0,0,0,1,-1,-1,0,1,-1),5,4),diag(K))
A5 <- direct.sum(matrix(1), cbind(-1,1), cbind(1,-1), diag(K))
A6 <- matrix(c(1,rep(0,K),-1/2, rep(1,K)-3/2, 1/2, rep(1,K)-1/2,rbind(rep(0,K),diag(K))),K+1, K+3)
A7 <- cbind(rep(-1,K), diag(K))

Zy <- A7%*%exp(A6%*%log(A5%*%exp(A4%*%log(A3%*%exp(A2%*%log(A1%*%n))))))

g0 <- n
g1 <- log(A1%*%n)
g2 <- exp(A2%*%log(A1%*%n))
g3 <- log(A3%*%exp(A2%*%log(A1%*%n)))
g4 <- exp(A4%*%log(A3%*%exp(A2%*%log(A1%*%n))))
g5 <- log(A5%*%exp(A4%*%log(A3%*%exp(A2%*%log(A1%*%n)))))
g6 <- exp(A6%*%log(A5%*%exp(A4%*%log(A3%*%exp(A2%*%log(A1%*%n))))))

G0 <- diag(K)
G1 <- solve(diag(c(A1%*%g0),length(A1%*%g0),length(A1%*%g0)))%*%A1%*%G0
G2 <- diag(c(exp(A2%*%g1)),length(exp(A2%*%g1)),length(exp(A2%*%g1)))%*%A2%*%G1
G3 <- solve(diag(c(A3%*%g2),length(A3%*%g2),length(A3%*%g2)))%*%A3%*%G2
G4 <- diag(c(exp(A4%*%g3)),length(exp(A4%*%g3)),length(exp(A4%*%g3)))%*%A4%*%G3
G5 <- solve(diag(c(A5%*%g4),length(A5%*%g4),length(A5%*%g4)))%*%A5%*%G4
G6 <- diag(c(exp(A6%*%g5)),length(exp(A6%*%g5)),length(exp(A6%*%g5)))%*%A6%*%G5
G7 <- A7 %*% G6
GZy <-G7

VZy <- GZy%*%diag(c(n),length(n),length(n))%*%t(GZy)-GZy%*%(n%*%t(n)/N)%*%t(GZy)-GZy%*%(n%*%t(n)/N)%*%t(GZy)
Se.Zy <- sqrt(diag(VZy))

loZy <- Zy - 1.96*Se.Zy
upZy <- Zy + 1.96*Se.Zy

A1 <- rbind(R,R^2,rep(1,K),rep(1,K))
A2 <- matrix(c(2,0,0,1,0,0,1,0,0,0,0,0,0,0,1,-1,0,1,-1,-1),5,4)
A3 <- matrix(c(-1,0,0,1,0,0,0,1,0,0,0,1,0,-1,0),3,5)
A4 <- matrix(c(1/2,0,-1/2,0,0,1),2,3)
A5 <- cbind(seq(-1.75,1.75,.50),rep(1,8))

Sty <- A5 %*% exp (A4 %*% log ( A3 %*% exp ( A2 %*% log ( A1 %*% n )))) -adjust

g0 <- n
g1 <- log(A1%*%n)
g2 <- exp(A2%*%log(A1%*%n))
g3 <- log(A3%*%exp(A2%*%log(A1%*%n)))
g4 <- exp(A4%*%log(A3%*%exp(A2%*%log(A1%*%n))))

G0 <- diag(K)
G1 <- solve(diag(c(A1%*%g0),length(A1%*%g0),length(A1%*%g0)))%*%A1%*%G0
G2 <- diag(c(exp(A2%*%g1)),length(exp(A2%*%g1)),length(exp(A2%*%g1)))%*%A2%*%G1
G3 <- solve(diag(c(A3%*%g2),length(A3%*%g2),length(A3%*%g2)))%*%A3%*%G2
G4 <- diag(c(exp(A4%*%g3)),length(exp(A4%*%g3)),length(exp(A4%*%g3)))%*%A4%*%G3
G5 <- A5 %*% G4
GSty <-G5

VSty <- GSty%*%diag(c(n),length(n),length(n))%*%t(GSty)-GSty%*%(n%*%t(n)/N)%*%t(GSty)
Se.Sty <- sqrt(diag(VSty))

loSty <- Sty - 1.96*Se.Sty
upSty <- Sty + 1.96*Se.Sty

A1 <- matrix(1,K+1,K)
A1[upper.tri(A1,diag=FALSE)] <- 0
A2 <- cbind(diag(K),rep(-1,K))
A3 <- matrix(0,K,K)
for (i in 0:K) {
i + 1
A3[i,i-1] <- 50
A3[i,i] <- 50
}

Pry <- A3%*%exp(A2%*%log(A1%*%n))

g0 <- n
g1 <- log(A1%*%n)
g2 <- exp(A2%*%log(A1%*%n))

G0 <- diag(K)
G1 <- solve(diag(c(A1%*%g0),length(A1%*%g0),length(A1%*%g0)))%*%A1%*%G0
G2 <- diag(c(exp(A2%*%g1)),length(exp(A2%*%g1)),length(exp(A2%*%g1)))%*%A2%*%G1
G3 <- A3 %*% G2
GPry <-G3

VPry <- GPry%*%diag(c(n),length(n),length(n))%*%t(GPry)-GPry%*%(n%*%t(n)/N)%*%t(GPry)
Se.Pry <- sqrt(diag(VPry))

loPry <- Pry - 1.96*Se.Pry
upPry <- Pry + 1.96*Se.Pry

if (nice.output){
 output.matrix.mean. <- matrix(NA, 1, 4)
 output.matrix.mean.[, 1] <- format(formatC(round(meany,3), digits = 3, format = "f"), width = 7, justify = "right")
 output.matrix.mean.[, 2] <- format(paste("(",formatC(round(Se.mean, 3), digits = 3, format = "f"),")", sep = ""), width = 7, justify = "right")
 output.matrix.mean.[, 3] <- format(formatC(round(lomeany,3), digits = 3, format = "f"), width = 7, justify = "right")
 output.matrix.mean.[, 4] <- format(formatC(round(upmeany,3), digits = 3, format = "f"), width = 7, justify = "right")
 dimnames(output.matrix.mean.) <- list("", c("Mean","SE", "lo", "up"))
 output.matrix.mean. <- noquote(output.matrix.mean.)

output.matrix.sd. <- matrix(NA, 1, 4)
 output.matrix.sd.[, 1] <- format(formatC(round(sdy,3), digits = 3, format = "f"), width = 7, justify = "right")
 output.matrix.sd.[, 2] <- format(paste("(",formatC(round(Se.sd, 3), digits = 3, format = "f"),")", sep = ""), width = 7, justify = "right")
 output.matrix.sd.[, 3] <- format(formatC(round(losdy,3), digits = 3, format = "f"), width = 7, justify = "right")
 output.matrix.sd.[, 4] <- format(formatC(round(upsdy,3), digits = 3, format = "f"), width = 7, justify = "right")
 dimnames(output.matrix.sd.) <- list("", c("SD","SE", "lo", "up"))
 output.matrix.sd. <- noquote(output.matrix.sd.)

output.matrix.Zy. <- matrix(NA, K, 4)
 output.matrix.Zy.[, 1] <- format(formatC(round(Zy,3), digits = 3, format = "f"), width = 7,justify = "right")
 output.matrix.Zy.[, 2] <- format(paste("(",formatC(round(Se.Zy, 3), digits = 3, format = "f"),")", sep = ""), width = 7, justify = "right")
 output.matrix.Zy.[, 3] <- format(formatC(round(loZy,3), digits = 3, format = "f"), width = 7,justify = "right")
 output.matrix.Zy.[, 4] <- format(formatC(round(upZy,3), digits = 3, format = "f"), width = 7,justify = "right")
 dimnames(output.matrix.Zy.) <- list(R-adjust, c("Zscores","SE", "lo", "up"))
 output.matrix.Zy. <- noquote(output.matrix.Zy.)

output.matrix.Sty. <- matrix(NA, 8, 4)
 output.matrix.Sty.[, 1] <- format(formatC(round(Sty,3), digits = 3, format = "f"), width = 7,justify = "right")
 output.matrix.Sty.[, 2] <- format(paste("(",formatC(round(Se.Sty, 3), digits = 3, format = "f"),")", sep = ""), width = 7, justify = "right")
 output.matrix.Sty.[, 3] <- format(formatC(round(loSty,3), digits = 3, format = "f"), width = 7,justify = "right")
 output.matrix.Sty.[, 4] <- format(formatC(round(upSty,3), digits = 3, format = "f"), width = 7,justify = "right")
 dimnames(output.matrix.Sty.) <- list(c("1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9"), c("Stanines","SE", "lo", "up"))
 output.matrix.Sty. <- noquote(output.matrix.Sty.)

output.matrix.Pry. <- matrix(NA, K, 4)
 output.matrix.Pry.[, 1] <- format(formatC(round(Pry,3), digits = 3, format = "f"), width = 7,justify = "right")
 output.matrix.Pry.[, 2] <- format(paste("(",formatC(round(Se.Pry, 3), digits = 3, format = "f"),")", sep = ""), width = 7, justify = "right")
 output.matrix.Pry.[, 3] <- format(formatC(round(loPry,3), digits = 3, format = "f"), width = 7,justify = "right")
 output.matrix.Pry.[, 4] <- format(formatC(round(upPry,3), digits = 3, format = "f"), width = 7,justify = "right")
 dimnames(output.matrix.Pry.) <- list(R-adjust, c("Percentiles","SE", "lo", "up"))
 output.matrix.Pry. <- noquote(output.matrix.Pry.)
} else {
 output.matrix.mean. <- matrix(NA, 1, 4)
 output.matrix.mean.[, 1] <- meany
 output.matrix.mean.[, 2] <- Se.mean
 output.matrix.mean.[, 3] <- lomeany
 output.matrix.mean.[, 4] <- upmeany
 dimnames(output.matrix.mean.) <- list("", c("Mean","SE", "lo", "up"))

output.matrix.sd. <- matrix(NA, 1, 4)
 output.matrix.sd.[, 1] <- sdy
 output.matrix.sd.[, 2] <- Se.sd
 output.matrix.sd.[, 3] <- losdy
 output.matrix.sd.[, 4] <- upsdy
 dimnames(output.matrix.sd.) <- list("", c("SD","SE", "lo", "up"))

output.matrix.Zy. <- matrix(NA, K, 4)
 output.matrix.Zy.[, 1] <- R-adjust
 output.matrix.Zy.[, 2] <- n
 output.matrix.Zy.[, 3] <- loZy
 output.matrix.Zy.[, 4] <- upZy
 dimnames(output.matrix.Zy.) <- list(R-adjust, c("Zscores","SE", "lo", "up"))

output.matrix.Sty. <- matrix(NA, 8, 4)
 output.matrix.Sty.[, 1] <- Sty
 output.matrix.Sty.[, 2] <- Se.Sty
 output.matrix.Sty.[, 3] <- loSty
 output.matrix.Sty.[, 4] <- upSty
 dimnames(output.matrix.Sty.) <- list(c("1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9"), c("Stanines","SE", "lo", "up"))

output.matrix.Pry. <- matrix(NA, K, 4)
 output.matrix.Pry.[, 1] <- R-adjust
 output.matrix.Pry.[, 2] <- n
 output.matrix.Pry.[, 3] <- loPry
 output.matrix.Pry.[, 4] <- upPry
 dimnames(output.matrix.Pry.) <- list(R-adjust, c("Percentiles","SE", "lo", "up"))
}

return(list(mean = output.matrix.mean., sd = output.matrix.sd., z = output.matrix.Zy., sta9 = output.matrix.Sty., perc = output.matrix.Pry.))

}
