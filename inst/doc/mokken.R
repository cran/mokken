###################################################
### chunk number 1: intro1
###################################################
#line 122 "mokken.Rnw"
6 - 3


###################################################
### chunk number 2: intro2
###################################################
#line 129 "mokken.Rnw"
options(width=60)
x <- sqrt(2)
x


###################################################
### chunk number 3: intro3
###################################################
#line 140 "mokken.Rnw"
library(mokken)


###################################################
### chunk number 4: intro5 eval=FALSE
###################################################
## #line 148 "mokken.Rnw"
## q()


###################################################
### chunk number 5: spss1 eval=FALSE
###################################################
## #line 195 "mokken.Rnw"
## library(foreign)
## ExampleR <- data.frame(read.spss("C:/ExampleSPSS.sav"))
## fix(ExampleR)


###################################################
### chunk number 6: spss2 eval=FALSE
###################################################
## #line 213 "mokken.Rnw"
## save(ExampleR, file="C:/ExampleR.Rdata")


###################################################
### chunk number 7: spss2 eval=FALSE
###################################################
## #line 221 "mokken.Rnw"
## load("C:/ExampleR.Rdata")


###################################################
### chunk number 8: spss3 eval=FALSE
###################################################
## #line 259 "mokken.Rnw"
## library(foreign)
## write.foreign(ExampleR, datafile="C:/ExampleSPSS.txt", codefile="C:/ExampleSPSS.SPS", package="SPSS")


###################################################
### chunk number 9: sas1 eval=FALSE
###################################################
## #line 272 "mokken.Rnw"
## library(foreign)
## ExampleR <- data.frame(read.xport("C:/ExampleSAS.xpt"))
## fix(ExampleR)


###################################################
### chunk number 10: sas2 eval=FALSE
###################################################
## #line 287 "mokken.Rnw"
## save(ExampleR, file="C:/ExampleR.Rdata")


###################################################
### chunk number 11: sas3 eval=FALSE
###################################################
## #line 295 "mokken.Rnw"
## load("C:/ExampleR.Rdata")


###################################################
### chunk number 12: sas4 eval=FALSE
###################################################
## #line 311 "mokken.Rnw"
## library(foreign)
## write.foreign(ExampleR, datafile="C:/ExampleSAS.txt",codefile="C:/ExampleSAS.XXX", package="SAS")


###################################################
### chunk number 13: stat1 eval=FALSE
###################################################
## #line 324 "mokken.Rnw"
## library(foreign)
## ExampleR <- data.frame(read.dta("C:/ExampleSTATA.dta"))
## fix(ExampleR)


###################################################
### chunk number 14: stat2 eval=FALSE
###################################################
## #line 339 "mokken.Rnw"
## save(ExampleR, file="C:/ExampleR.Rdata")


###################################################
### chunk number 15: stat3 eval=FALSE
###################################################
## #line 347 "mokken.Rnw"
## load("C:/ExampleR.Rdata")


###################################################
### chunk number 16: stat1 eval=FALSE
###################################################
## #line 363 "mokken.Rnw"
## library(foreign)
## write.foreign(ExampleR, datafile="C:/ExampleSTATA.dat", codefile="C:/ExampleSTATA.do", package="Stata")


###################################################
### chunk number 17: splus1 eval=FALSE
###################################################
## #line 378 "mokken.Rnw"
## library(foreign)
## ExampleR <- data.frame(read.s("C:/ExampleSplus.ssc"))
## fix(ExampleR)


###################################################
### chunk number 18: splus2 eval=FALSE
###################################################
## #line 393 "mokken.Rnw"
## save(ExampleR, file="C:/ExampleR.Rdata")


###################################################
### chunk number 19: splus3 eval=FALSE
###################################################
## #line 401 "mokken.Rnw"
## load("C:/ExampleR.Rdata")


###################################################
### chunk number 20: splus4 eval=FALSE
###################################################
## #line 413 "mokken.Rnw"
## dump(ExampleSplus,"C:/Example.dmp")


###################################################
### chunk number 21: splus5 eval=FALSE
###################################################
## #line 421 "mokken.Rnw"
## ExampleR <- dget("C:/Example.dmp")


###################################################
### chunk number 22: mokken1
###################################################
#line 435 "mokken.Rnw"
library(mokken)


###################################################
### chunk number 23: mokken2 eval=FALSE
###################################################
## #line 443 "mokken.Rnw"
## help(mokken)


###################################################
### chunk number 24: mokken3
###################################################
#line 452 "mokken.Rnw"
# help(mokken)


###################################################
### chunk number 25: mokken4
###################################################
#line 462 "mokken.Rnw"
data(acl)
data(cavalini)
data(transreas)


###################################################
### chunk number 26: mokken5 eval=FALSE
###################################################
## #line 472 "mokken.Rnw"
## help(acl)


###################################################
### chunk number 27: mokken5 eval=FALSE
###################################################
## #line 480 "mokken.Rnw"
## fix(cavalini)


###################################################
### chunk number 28: mokken6
###################################################
#line 490 "mokken.Rnw"
X <- acl
Y <- 3
Z <- c(1,2,3,8:11)


###################################################
### chunk number 29: mokken7
###################################################
#line 503 "mokken.Rnw"
Y
Z


###################################################
### chunk number 30: mokken8
###################################################
#line 512 "mokken.Rnw"
X1 <- acl[,1]


###################################################
### chunk number 31: mokken9
###################################################
#line 520 "mokken.Rnw"
X2 <- acl[,11:20]


###################################################
### chunk number 32: mokken10
###################################################
#line 528 "mokken.Rnw"
X3 <- acl[1:10,]


###################################################
### chunk number 33: mokken11
###################################################
#line 536 "mokken.Rnw"
X4 <- acl[232,133]


###################################################
### chunk number 34: mokken12
###################################################
#line 544 "mokken.Rnw"
scale.1 <- c(1,2,4)
X5 <- acl[c(1:100,201:300),scale.1]


###################################################
### chunk number 35: mokken13
###################################################
#line 553 "mokken.Rnw"
X6 <- acl[acl[,1]==2,]


###################################################
### chunk number 36: mokken14
###################################################
#line 567 "mokken.Rnw"
dimnames(acl)[[1]] <- 1:nrow(acl)


###################################################
### chunk number 37: E1
###################################################
#line 592 "mokken.Rnw"
data(acl)
Communality <- acl[,1:10]
scale <- aisp(Communality, verbose=FALSE)
scale


###################################################
### chunk number 38: E1A eval=FALSE
###################################################
## #line 605 "mokken.Rnw"
## scale <- aisp(Communality, search="ga")


###################################################
### chunk number 39: E1B eval=FALSE
###################################################
## #line 613 "mokken.Rnw"
## scale  <- aisp(Communality, lowerbound = .2, alpha =.10)


###################################################
### chunk number 40: E1C eval=FALSE
###################################################
## #line 621 "mokken.Rnw"
## aisp(Communality, verbose=TRUE)


###################################################
### chunk number 41: E1D eval=FALSE
###################################################
## #line 629 "mokken.Rnw"
## help(aisp)


###################################################
### chunk number 42: E2
###################################################
#line 648 "mokken.Rnw"
data(acl)
Communality <- acl[,1:10]
coefficients <- coefH(Communality)
coefficients$Hij
coefficients$Hi
coefficients$H


###################################################
### chunk number 43: E3
###################################################
#line 673 "mokken.Rnw"
data(acl)
Communality <- acl[,1:10]
iio.results <- check.iio(Communality)
summary(iio.results)


###################################################
### chunk number 44: E3A eval=FALSE
###################################################
## #line 688 "mokken.Rnw"
## check.iio(Communality, minvi=0.00, minsize=50)


###################################################
### chunk number 45: E3B eval=FALSE
###################################################
## #line 696 "mokken.Rnw"
## summary(check.iio(Communality, method="MS-CPM"))
## summary(check.iio(Communality, method="IT"))


###################################################
### chunk number 46: E3C eval=FALSE
###################################################
## #line 705 "mokken.Rnw"
## summary(check.iio(Communality, alpha=.01))


###################################################
### chunk number 47: E3D eval=FALSE
###################################################
## #line 713 "mokken.Rnw"
## summary(check.iio(Communality, item.selection=FALSE))
## summary(check.iio(Communality, verbose=TRUE))


###################################################
### chunk number 48: E3D eval=FALSE
###################################################
## #line 722 "mokken.Rnw"
## help(check.iio)


###################################################
### chunk number 49: E4 eval=FALSE
###################################################
## #line 743 "mokken.Rnw"
## data(acl)
## Communality <- acl[,1:10]
## monotonicity.results <- check.monotonicity(Communality)
## summary(monotonicity.results)
## plot(monotonicity.results, items = c(1,2))


###################################################
### chunk number 50: E4
###################################################
#line 753 "mokken.Rnw"
Communality <- acl[,1:10]
monotonicity.results <- check.monotonicity(Communality)
summary(monotonicity.results)


###################################################
### chunk number 51: E4Afig
###################################################
#line 761 "mokken.Rnw"
layout(matrix(1:1,nrow=1,byrow=TRUE))
plot(monotonicity.results, ask=FALSE, items = 1)


###################################################
### chunk number 52: E4Bfig
###################################################
#line 768 "mokken.Rnw"
plot(monotonicity.results, ask=FALSE, items = 2)


###################################################
### chunk number 53: E4Afigplot
###################################################
#line 777 "mokken.Rnw"
#line 761 "mokken.Rnw"
layout(matrix(1:1,nrow=1,byrow=TRUE))
plot(monotonicity.results, ask=FALSE, items = 1)
#line 778 "mokken.Rnw"


###################################################
### chunk number 54: E4Bfigplot
###################################################
#line 784 "mokken.Rnw"
#line 768 "mokken.Rnw"
plot(monotonicity.results, ask=FALSE, items = 2)
#line 785 "mokken.Rnw"


###################################################
### chunk number 55: E4A eval=FALSE
###################################################
## #line 800 "mokken.Rnw"
## check.monotonicity(Communality, minvi=0.00, minsize=50)


###################################################
### chunk number 56: E4B eval=FALSE
###################################################
## #line 808 "mokken.Rnw"
## plot(check.monotonicity(Communality), item=c(1,2))


###################################################
### chunk number 57: E4C eval=FALSE
###################################################
## #line 817 "mokken.Rnw"
## pdf("monotonicity.pdf")
## plot(monotonicity.results, ask=FALSE)
## dev.off()


###################################################
### chunk number 58: E4D eval=FALSE
###################################################
## #line 827 "mokken.Rnw"
## help(check.monotonicity)


###################################################
### chunk number 59: E5 eval=FALSE
###################################################
## #line 850 "mokken.Rnw"
## data(acl)
## Communality <- acl[,1:10]
## pmatrix.results <- check.pmatrix(Communality)
## summary(pmatrix.results)
## plot(pmatrix.results)


###################################################
### chunk number 60: E5
###################################################
#line 860 "mokken.Rnw"
pmatrix.results <- check.pmatrix(Communality)
summary(pmatrix.results)


###################################################
### chunk number 61: E5A eval=FALSE
###################################################
## #line 873 "mokken.Rnw"
## check.pmatrix(Communality, minvi=0.00)


###################################################
### chunk number 62: E5B eval=FALSE
###################################################
## #line 881 "mokken.Rnw"
## plot(check.pmatrix(Communality), pmatrix="ppp", item=c(1,2))
## plot(check.pmatrix(Communality), pmatrix="pmm", item=5)


###################################################
### chunk number 63: E5C eval=FALSE
###################################################
## #line 891 "mokken.Rnw"
## pdf("pmatrix.pdf")
## plot(pmatrix.results, ask=FALSE)
## dev.off()


###################################################
### chunk number 64: E5D eval=FALSE
###################################################
## #line 901 "mokken.Rnw"
## help(check.pmatrix)


###################################################
### chunk number 65: E6
###################################################
#line 917 "mokken.Rnw"
data(acl)
Communality <- acl[,1:10]
check.reliability(Communality)


###################################################
### chunk number 66: E7 eval=FALSE
###################################################
## #line 941 "mokken.Rnw"
## data(acl)
## Communality <- acl[,1:10]
## restscore.results <- check.restscore(Communality)
## summary(restscore.results)
## plot(restscore.results, item.pairs = c(1,2))


###################################################
### chunk number 67: E7
###################################################
#line 951 "mokken.Rnw"
Communality <- acl[,1:10]
restscore.results <- check.restscore(Communality)
summary(restscore.results)


###################################################
### chunk number 68: E7Afig
###################################################
#line 959 "mokken.Rnw"
plot(restscore.results, ask=FALSE, item.pairs = 1)


###################################################
### chunk number 69: E7Bfig
###################################################
#line 965 "mokken.Rnw"
plot(restscore.results, ask=FALSE, item.pairs = 1)


###################################################
### chunk number 70: E7Afigplot
###################################################
#line 974 "mokken.Rnw"
#line 959 "mokken.Rnw"
plot(restscore.results, ask=FALSE, item.pairs = 1)
#line 975 "mokken.Rnw"


###################################################
### chunk number 71: E7Bfigplot
###################################################
#line 981 "mokken.Rnw"
#line 965 "mokken.Rnw"
plot(restscore.results, ask=FALSE, item.pairs = 1)
#line 982 "mokken.Rnw"


###################################################
### chunk number 72: E7A eval=FALSE
###################################################
## #line 998 "mokken.Rnw"
## check.restscore(Communality, minvi=0.00, minsize=50)


###################################################
### chunk number 73: E7B eval=FALSE
###################################################
## #line 1006 "mokken.Rnw"
## plot(check.restscore(Communality))


###################################################
### chunk number 74: E7C eval=FALSE
###################################################
## #line 1015 "mokken.Rnw"
## pdf("restscore.pdf")
## plot(restscore.results, ask=FALSE)
## dev.off()


###################################################
### chunk number 75: E7D eval=FALSE
###################################################
## #line 1025 "mokken.Rnw"
## help(check.restscore)


###################################################
### chunk number 76: E8
###################################################
#line 1041 "mokken.Rnw"
data(acl)
Communality <- acl[,1:10]
Group <- acl[,11]
coefH(Communality[Group==0|Group==1,])$H
coefH(Communality[Group==2,])$H
coefH(Communality[Group==3,])$H
coefH(Communality[Group==4,])$H


###################################################
### chunk number 77: Table.3.1A
###################################################
#line 1060 "mokken.Rnw"
library(mokken)
data(transreas)
grades <- transreas[,1]
item.scores <- transreas[,-1]


###################################################
### chunk number 78: Table.3.1B
###################################################
#line 1070 "mokken.Rnw"
apply(item.scores,2,mean)
apply(item.scores[grades==2,],2,mean)
apply(item.scores[grades==3,],2,mean)
apply(item.scores[grades==4,],2,mean)
apply(item.scores[grades==5,],2,mean)
apply(item.scores[grades==6,],2,mean)


###################################################
### chunk number 79: Table.3.1C
###################################################
#line 1083 "mokken.Rnw"
Total.group <- round(apply(item.scores,2,mean),2)
for (i in 2:6) assign(paste("Grade.",i,sep=""),
 round(apply(item.scores[grades==i,],2,mean),2))
Task <- c(9,12,10,11,4,5,2,7,3,1,8,6)
Property <- attributes(transreas)$property
Format <- attributes(transreas)$format
Objects <- attributes(transreas)$objects
Measures <- attributes(transreas)$measures
Table.3.1 <- data.frame(Task,Property,Format,Objects,Measures, Total.group,Grade.2,Grade.3,Grade.4,Grade.5,Grade.6)
Table.3.1


###################################################
### chunk number 80: Table.3.2B
###################################################
#line 1107 "mokken.Rnw"
coefH(item.scores,FALSE)$Hi
coefH(item.scores,FALSE)$H
coefZ(item.scores)$Zi
coefZ(item.scores)$Z


###################################################
### chunk number 81: Table.3.2C
###################################################
#line 1119 "mokken.Rnw"
coefH(item.scores[,-c(2,4)],FALSE)$Hi
coefH(item.scores[,-c(2,4)],FALSE)$H
coefZ(item.scores[,-c(2,4)])$Zi
coefZ(item.scores[,-c(2,4)])$Z


###################################################
### chunk number 82: Table.3.2D
###################################################
#line 1129 "mokken.Rnw"
Task <- c("9","12","10","11","4","5","2","7","3","1","8","6", "Total item set")
Property <- c(attributes(transreas)$property,"")
Format <- c(attributes(transreas)$format,"")
Table.3.2 <- data.frame(Task,Property,Format,matrix(NA,13,8))
analysis <- list(c(1:12),c(1,3,5:12),c(1,3,6,8:12),c(1,3,8:12))
k <- 3
for (i in 1:4) for (j in 1:2){
 k <- k + 1
 Table.3.2[c(analysis[[i]],13),k] <-
 c(round(coefH(item.scores[,analysis[[i]]],FALSE)$Hi,2),
 round(coefH(item.scores[,analysis[[i]]],FALSE)$H,2))
}
dimnames(Table.3.2)[[2]][4:11] <- paste(c("k=12","k=12",
 "k=10","k=10","k=8","k=8","k=7","k=7"),c("Hi","Zi"))
Table.3.2


###################################################
### chunk number 83: Table.5.1A
###################################################
#line 1158 "mokken.Rnw"
options(width=60)
scale <- aisp(item.scores, verbose=FALSE)


###################################################
### chunk number 84: Table.5.1B
###################################################
#line 1167 "mokken.Rnw"
options(width=60)
scale.1 <- c(12,8,1,11,9,3,10)
scale.2 <- c(7,5)
Hi.top <- matrix(NA,8,6)
for (i in 1:6)
Hi.top[1:(i+1),i] <- round(coefH(item.scores[,scale.1[1:(i+1)]],FALSE)$Hi,2)
for (i in 1:6)
Hi.top[8,i] <- round(coefH(item.scores[,scale.1[1:(i+1)]],FALSE)$H,2)
dimnames(Hi.top)[[2]] <- paste("Step",1:6)
Table.5.1.top <- data.frame(
 Task = c(Task[scale.1],"Total H"),
 Property = c(Property[scale.1],""),
 Format = c(Format[scale.1],""),
 Pi = c( round(apply(item.scores[,scale.1],2,mean),2),NA)
)
 Table.5.1.top <- cbind(Table.5.1.top,Hi.top)
 Table.5.1.top


###################################################
### chunk number 85: Table.5.2A
###################################################
#line 1196 "mokken.Rnw"
data(cavalini)
X <- cavalini
X[cavalini < 2] <- 0
X[cavalini > 1] <- 1
apply(X,2,mean)


###################################################
### chunk number 86: Table.5.2B
###################################################
#line 1208 "mokken.Rnw"
Table.5.2 <- data.frame(1:17, attributes(X)$labels,
        round(apply(X,2,mean),2))
dimnames(Table.5.2)[[2]] <- c("Item.number","Item.text","Pi")
rownames(Table.5.2) <- NULL
Table.5.2


###################################################
### chunk number 87: Table.5.3A eval=FALSE
###################################################
## #line 1228 "mokken.Rnw"
## aisp(X,lowerbound=0.00, verbose = FALSE)
## aisp(X,lowerbound=0.05, verbose = FALSE)
## aisp(X,lowerbound=0.10, verbose = FALSE)
## # etc.


###################################################
### chunk number 88: Table.5.3B
###################################################
#line 1239 "mokken.Rnw"
lower.bound <- seq(0,.6,by=.05)
scaling.results <- matrix(NA,length(lower.bound),ncol(X))
for (i in 1:length(lower.bound)) scaling.results[i,] <- aisp(X, lowerbound=lower.bound[i],verbose=FALSE)
equal <- function(x,n) which(x==n)
scale.1 <- sapply(apply(scaling.results,1,"equal", 1), paste,collapse=" ")
scale.2 <- sapply(apply(scaling.results,1,"equal", 2), paste,collapse=" ")
scale.3 <- sapply(apply(scaling.results,1,"equal", 3), paste,collapse=" ")
scale.4 <- sapply(apply(scaling.results,1,"equal", 4), paste,collapse=" ")
scale.5 <- sapply(apply(scaling.results,1,"equal", 5), paste,collapse=" ")
Table.5.3 <- data.frame(lower.bound, scale.1,scale.2, scale.3,scale.4,scale.5)
Table.5.3


###################################################
### chunk number 89: Table.5.4A
###################################################
#line 1264 "mokken.Rnw"
scale.3 <- aisp(X,lowerbound=0.30)
scale.35 <- aisp(X,lowerbound=0.35)


###################################################
### chunk number 90: Table.5.4B
###################################################
#line 1273 "mokken.Rnw"
scale.30 <- aisp(X,lowerbound=0.30,verbose=F)
max.scale <- max(scale.30)
Table.5.4.left <- data.frame()
for (i in 1:max.scale){
  max.item <- max(length(scale.30[scale.30==i]))
  Scale <- c(i,rep("",max.item-1))
  Item.30 <- which(scale.30==i)
  Hi.30 <- round(coefH(X[,scale.30==i],FALSE)$Hi,2)
  H.30 <- c(rep("",max.item-1),round(coefH(X[,scale.30==i],FALSE)$H,2))
  Table.5.4.left <- rbind(Table.5.4.left,data.frame(Scale=Scale,
     Item=Item.30,Hi=Hi.30,H=H.30),c("","","",""))
}
rownames(Table.5.4.left) <- NULL
Table.5.4.left


###################################################
### chunk number 91: Table.5.4C
###################################################
#line 1292 "mokken.Rnw"
scale.35 <- aisp(X,lowerbound=0.35,verbose=F)
max.scale <- max(scale.35)
Table.5.4.right <- data.frame()
for (i in 1:max.scale){
  max.item <- max(length(scale.35[scale.35==i]))
  Scale <- c(i,rep("",max.item-1))
  Item.35 <- which(scale.35==i)
  Hi.35 <- round(coefH(X[,scale.35==i],FALSE)$Hi,2)
  H.35 <- c(rep("",max.item-1),round(coefH(X[,scale.35==i],FALSE)$H,2))
  Table.5.4.right <- rbind(Table.5.4.right,data.frame(Scale=Scale,
     Item=Item.35,Hi=Hi.35,H=H.35),c("","","",""))
}
rownames(Table.5.4.right) <- NULL
Table.5.4.right


###################################################
### chunk number 92: Table.6.1
###################################################
#line 1322 "mokken.Rnw"
library(mokken)
data(transreas)
X <- transreas[,-c(1,3,5)]
check.restscore(X,minsize=2)$results[[21]]
check.restscore(X,minsize=40)$results[[21]]
plot(check.restscore(X,minsize=2),item.pairs=21)
plot(check.restscore(X,minsize=40),item.pairs=21)
R <- apply(X[,-c(3,7)],1,sum)
table(X[,3],X[,7],R)
as.numeric(table(X[,3][R < 5],X[,7][R < 5]))


###################################################
### chunk number 93: Table.6.2A
###################################################
#line 1346 "mokken.Rnw"
library(mokken)
data(transreas)
X <- transreas[,-c(1,3,5)]
Task <- c(9,10,4,5,2,7,3,1,8,6)
ppp <- check.pmatrix(X)$Ppp
dimnames(ppp) <- list(Task,Task)
round(ppp,2)


###################################################
### chunk number 94: Table.6.2B
###################################################
#line 1359 "mokken.Rnw"
pmm <- check.pmatrix(X)$Pmm
dimnames(pmm) <- list(Task,Task)
round(pmm,2)


