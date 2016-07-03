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

eps <- 1e-8
y <- sort(y)
y[y > -eps & y < eps] <- eps
R <- unique(y)
K <- length(R)
N <- length(y)          
n <- apply(matrix(R), 1, function(x) sum(y %in% x))

A1 <- rbind(R,1)
A2 <- matrix(c(1,-1),1,2)

meany <- exp(A2%*%log(A1%*%n))


g0 <- n
g1 <- log(A1%*%n)

G0 <- diag(K)
G1 <- solve(diag(c(A1%*%g0),length(A1%*%g0),length(A1%*%g0)))%*%A1%*%G0
G2 <- diag(c(exp(A2%*%g1)),length(exp(A2%*%g1)),length(exp(A2%*%g1)))%*%A2%*%G1
Gmean <-G2
Vmean <- Gmean%*%diag(c(n),length(n),length(n))%*%t(Gmean)-Gmean%*%(n%*%t(n)/N)%*%t(Gmean)
Se.mean <- sqrt(Vmean)


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


A1 <- rbind(R,R^2,rep(1,K),rep(1,K))
A2 <- matrix(c(2,0,0,1,0,0,1,0,0,0,0,0,0,0,1,-1,0,1,-1,-1),5,4)
A3 <- matrix(c(-1,0,0,1,0,0,0,1,0,0,0,1,0,-1,0),3,5)
A4 <- matrix(c(1/2,0,-1/2,0,0,1),2,3)
A5 <- cbind(seq(-1.75,1.75,.50),rep(1,8))

Sty <- A5 %*% exp (A4 %*% log ( A3 %*% exp ( A2 %*% log ( A1 %*% n ))))

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

if (nice.output){
   output.matrix.mean. <- matrix(NA, 1, 2)
   output.matrix.mean.[, 1] <- format(formatC(round(meany,3), digits = 3, format = "f"), width = 7, justify = "right")
   output.matrix.mean.[, 2] <- format(paste("(",formatC(round(Se.mean, 3), digits = 3, format = "f"),")", sep = ""), width = 7, justify = "right")
   dimnames(output.matrix.mean.) <- list("", c("Mean","SE"))
   output.matrix.mean. <- noquote(output.matrix.mean.)
   
   output.matrix.sd. <- matrix(NA, 1, 2)
   output.matrix.sd.[, 1] <- format(formatC(round(sdy,3), digits = 3, format = "f"), width = 7, justify = "right")
   output.matrix.sd.[, 2] <- format(paste("(",formatC(round(Se.sd, 3), digits = 3, format = "f"),")", sep = ""), width = 7, justify = "right")
   dimnames(output.matrix.sd.) <- list("", c("SD","SE"))
   output.matrix.sd. <- noquote(output.matrix.sd.)
   
   output.matrix.Zy. <- matrix(NA, K, 4)
   output.matrix.Zy.[, 1] <- format(formatC(round(R,0), digits = 3, format = "f"), width = 7,justify = "right")
   output.matrix.Zy.[, 2] <- format(formatC(round(n,0), digits = 0, format = "f"), width = 7,justify = "right")
   output.matrix.Zy.[, 3] <- format(formatC(round(Zy,3), digits = 3, format = "f"), width = 7,justify = "right")
   output.matrix.Zy.[, 4] <- format(paste("(",formatC(round(Se.Zy, 3), digits = 3, format = "f"),")", sep = ""), width = 7, justify = "right")
   dimnames(output.matrix.Zy.) <- list(R, c("Scores", "Freq", "Zscores","SE"))
   output.matrix.Zy. <- noquote(output.matrix.Zy.)

   output.matrix.Sty. <- matrix(NA, 8, 2)
   output.matrix.Sty.[, 1] <- format(formatC(round(Sty,3), digits = 3, format = "f"), width = 7,justify = "right")
   output.matrix.Sty.[, 2] <- format(paste("(",formatC(round(Se.Sty, 3), digits = 3, format = "f"),")", sep = ""), width = 7, justify = "right")
   dimnames(output.matrix.Sty.) <- list(c("1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9"), c("Stanines","SE"))
   output.matrix.Sty. <- noquote(output.matrix.Sty.)
   
   output.matrix.Pry. <- matrix(NA, K, 4)
   output.matrix.Pry.[, 1] <- format(formatC(round(R,0), digits = 3, format = "f"), width = 7,justify = "right")
   output.matrix.Pry.[, 2] <- format(formatC(round(n,0), digits = 0, format = "f"), width = 7,justify = "right")
   output.matrix.Pry.[, 3] <- format(formatC(round(Pry,3), digits = 3, format = "f"), width = 7,justify = "right")
   output.matrix.Pry.[, 4] <- format(paste("(",formatC(round(Se.Pry, 3), digits = 3, format = "f"),")", sep = ""), width = 7, justify = "right")
   dimnames(output.matrix.Pry.) <- list(R, c("Scores", "Freq","Percentiles","SE"))
   output.matrix.Pry. <- noquote(output.matrix.Pry.)
} else {
   output.matrix.mean. <- matrix(NA, 1, 2)
   output.matrix.mean.[, 1] <- meany
   output.matrix.mean.[, 2] <- Se.mean
   dimnames(output.matrix.mean.) <- list("", c("Mean","SE"))
   
   output.matrix.sd. <- matrix(NA, 1, 2)
   output.matrix.sd.[, 1] <- sdy
   output.matrix.sd.[, 2] <- Se.sd
   dimnames(output.matrix.sd.) <- list("", c("SD","SE"))
   
   output.matrix.Zy. <- matrix(NA, K, 4)
   output.matrix.Zy.[, 1] <- round(R, 0)
   output.matrix.Zy.[, 2] <- round(n, 0)
   output.matrix.Zy.[, 3] <- Zy
   output.matrix.Zy.[, 4] <- Se.Zy
   dimnames(output.matrix.Zy.) <- list(R, c("Scores", "Freq", "Zscores","SE"))
   
   output.matrix.Sty. <- matrix(NA, 8, 2)
   output.matrix.Sty.[, 1] <- Sty
   output.matrix.Sty.[, 2] <- Se.Sty
   dimnames(output.matrix.Sty.) <- list(c("1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9"), c("Stanines","SE"))

   output.matrix.Pry. <- matrix(NA, K, 4)
   output.matrix.Pry.[, 1] <- round(R, 0)
   output.matrix.Pry.[, 2] <- round(n, 0)
   output.matrix.Pry.[, 3] <- Pry
   output.matrix.Pry.[, 4] <- Se.Pry
   dimnames(output.matrix.Pry.) <- list(R, c("Scores", "Freq","Percentiles","SE"))
}
return(list(mean = output.matrix.mean., sd = output.matrix.sd., z = output.matrix.Zy., sta9 = output.matrix.Sty., perc = output.matrix.Pry.))

}
