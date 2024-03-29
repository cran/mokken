# Adapted 23 June 2023, Letty added level.two.var
"check.monotonicity" <-
function(X, minvi = .03, minsize = default.minsize, level.two.var = NULL){
  X <- check.data(X)
  N <- nrow(X)
  J <- ncol(X)
  m <- max(X) + 1
  default.minsize <- ifelse(N >= 500, floor(N/10), floor(N/5))
  default.minsize <- ifelse(N <= 250, floor(N/3), default.minsize)
  default.minsize <- ifelse(N <  150, 50, default.minsize)

  if (N < minsize) stop("Sample size less than Minsize")
  if (minsize > N/2) stop("Minsize value is too high")
  
  if(!is.null(level.two.var)) {
    X <- X[order(level.two.var), ]
    level.two.var <- sort(level.two.var)
  } 
  

  # Initial computation
  R <- as.matrix(X) %*% (matrix(1,J,J) - diag(J))
  results <- list()
  Ll1 <- NULL
  
  # results checks per item
  I.labels <- dimnames(X)[[2]]
  if(length(I.labels)==0) I.labels <- paste("C",1:ncol(X))
  for (j in 1:J){
    violation.matrix <- matrix(0,nrow=m,ncol=10)
    dimnames(violation.matrix) <- list(c(paste("P(X >=",1:(m-1),")",sep=""),"Total"), dimnames(violation.matrix)[[2]] <- c("#ac","#vi","#vi/#ac","maxvi","sum","sum/#ac","zmax","group","group","#zsig"))
    results[[j]] <- list()
    results[[j]][1] <- I.labels[j]
    sorted.R <- sort(R[,j])
    group <- max(which(sorted.R==sorted.R[minsize]))
    repeat{
      if(N - max(group) < minsize)break
      group <- c(group,max(which(sorted.R==sorted.R[minsize+max(group)])))
    }
    group <- group[-length(group)]
    summary.matrix <- matrix(nrow = length(group)+1,ncol = 4 + 2* m)
    dimnames(summary.matrix)[[2]] <- c("Group", "Lo Score", "Hi Score", "N", paste("F",0:(m-1)), "Mean", paste("P(X >=",1:(m-1),")",sep=""))
    summary.matrix[,1] <- 1:nrow(summary.matrix)
    summary.matrix[,4] <- c(group,N) - c(0,group)
    group <- c(sorted.R[group],max(sorted.R))
    L <- length(group)
    Ll1[j] <- L
    summary.matrix[,3] <- group
    summary.matrix[,2] <- c(min(sorted.R),group[-L]+1)

    member <- apply(1 - outer(R[,j], group, "<="),1,sum) + 1
    for (i in 1:L){
      Ni <- summary.matrix[i,4]
      freq <- tabulate(X[member==i,j]+1,m)
      summary.matrix[i,5:(m+4)] <- freq
      summary.matrix[i,m+5] <- sum(freq * min(X):max(X)) / Ni
      cum.freq <- rev(cumsum(rev(freq))/Ni)
      summary.matrix[i,(m+6):(2*m+4)] <- cum.freq[2:m]
    }
    results[[j]][[2]] <- summary.matrix
    nac <- rep(0,m-1)
    for (i in 1:(m-1)) nac[i] <- sum(matrix(rep(summary.matrix[,5 + m + i] > 1e-10,L),L,L,byrow=FALSE) * matrix(rep(summary.matrix[,5 + m + i] < .999999999999,L),L,L,TRUE) * upper.tri(matrix(,L,L)))
    violation.matrix[1:(m-1),1] <- nac
    violation.matrix[m,1] <- sum(nac)

    freq <- summary.matrix[,5:(m+4)]
    for (i in 1:(m-1)){
      V <- outer(summary.matrix[,(m+5+i)],summary.matrix[,(m+5+i)],"-")
      V[row(V) <= col(V)] <- 0
      V[V >= -minvi] <- 0
      violation.matrix[i,2] <- sum(ceiling(abs(V)))
      violation.matrix[i,4] <- max(abs(V))
      if(violation.matrix[i,4] > minvi){
        violation.matrix[i,5] <- sum(abs(V))
        freqd <- cbind(apply(as.matrix(freq[,1:i]),1,sum), apply(as.matrix(freq[,(i+1):m]),1,sum))
        Z <- abs(sign(-V) * 2 * (sqrt(outer(freqd[,2]+1,freqd[,1]+1)) - sqrt(outer(freqd[,1],freqd[,2]))) /
              sqrt(outer(freqd[,2],freqd[,1],"+") + outer(freqd[,1],freqd[,2],"+") - 1))
        violation.matrix[i,7] <- max(Z)
        violation.matrix[i,8] <- min(col(Z)[Z==max(Z)])
        violation.matrix[i,9] <- min(row(Z)[Z==max(Z)])
        violation.matrix[i,10] <- sum(sign(Z[Z > 1.6449]))
      }
    }
    violation.matrix[m,2] <- sum(violation.matrix[1:(m-1),2])
    violation.matrix[1:m,3] <- violation.matrix[1:m,2]/violation.matrix[1:m,1]
    violation.matrix[m,4] <- max(violation.matrix[1:(m-1),4])
    violation.matrix[m,5] <- sum(violation.matrix[1:(m-1),5])
    violation.matrix[1:m,6] <- violation.matrix[1:m,5]/violation.matrix[1:m,1]
    violation.matrix[m,7] <- max(violation.matrix[1:(m-1),7])
    violation.matrix[m,10] <- sum(violation.matrix[1:(m-1),10])
    results[[j]][[3]] <- violation.matrix
    results[[j]][[4]] <- paste("Minsize = ",minsize," Minvi = ",minvi,sep="")

  }
  
  if(!is.null(level.two.var)) {
    Xa <- aggregate(X, by = list(level.two.var), FUN = mean)[, -1] # Aggregated dataset
    Ra <- as.matrix(Xa) %*% (matrix(1, J, J) - diag(J)) # Aggregated rest score matrix
    Rs <- (table(level.two.var))
    repRs <- rep(Rs, Rs)
    Ra <- round(Ra[rep(1:nrow(Ra), Rs), ], 4) #round(Ra[,j], 4)
    resultsa <- list()
    
    for (j in 1:J) {
      violation.matrix <- matrix(0, nrow = m, ncol = 10)
      dimnames(violation.matrix) <- list(c(paste("P(X >=", 
                                                 1:(m - 1), ")", sep = ""), "Total"), 
                                         dimnames(violation.matrix)[[2]] <- c("#ac", 
                                                                              "#vi", "#vi/#ac", "maxvi", "sum", "sum/#ac", "zmax", 
                                                                              "group", "group", "#zsig"))
      resultsa[[j]] <- list()
      resultsa[[j]][1] <- I.labels[j]
      # Aggregated
      sorted.Ra <- sort(Ra[, j])
      minsizeA <- floor(nrow(X) / (Ll1[j] + 1.75))
      group <- max(which(sorted.Ra == sorted.Ra[minsizeA]))
      repeat {
        if (N - max(group) < minsizeA) 
          break
        group <- c(group, max(which(sorted.Ra == sorted.Ra[minsizeA + 
                                                             max(group)])))
      }
      group <- group[-length(group)]
      
      summary.matrix <- matrix(nrow = length(group) + 1, ncol = 4 + 2 * m)
      dimnames(summary.matrix)[[2]] <- c("Group", "Lo Score", 
                                         "Hi Score", "N", paste("F", 0:(m - 1)), 
                                         "Mean", paste("P(X >=", 
                                                       1:(m - 1), ")", sep = ""))
      summary.matrix[, 1] <- 1:nrow(summary.matrix)
      #summary.matrix[, 4] <- (c(group, N) - c(0, group)) #/ Rs
      
      group <- c(sorted.Ra[group], max(sorted.Ra))
      L <- length(group)
      summary.matrix[, 3] <- group
      summary.matrix[, 2] <- c(min(sorted.Ra), group[-L])
      member <- apply(1 - outer(Ra[, j], group, "<="), 1, sum) + 1
      #Xas <- Xa[rep(1:nrow(Xa), each = Rs), ]
      
      for (i in 1:L) {
        #Ni <- summary.matrix[i, 4]
        freq <- sapply(1:m, function(x) sum((X[member == i, j] + 1 == x) / repRs[member == i]))#tabulate(X[member == i, j] + 1, m)
        summary.matrix[i, 5:(m + 4)] <- freq 
        Ni <- summary.matrix[i, 4] <- sum(freq)
        summary.matrix[i, m + 5] <- sum(freq * min(X):max(X))/Ni
        cum.freq <- rev(cumsum(rev(freq))/Ni)
        summary.matrix[i, (m + 6):(2 * m + 4)] <- cum.freq[2:m]
      }
      resultsa[[j]][[2]] <- summary.matrix
      nac <- rep(0, m - 1)
      for (i in 1:(m - 1)) nac[i] <- 
        sum(matrix(rep(summary.matrix[,  5 + m + i] > 1e-10, L), L, L, byrow = FALSE) * 
              matrix(rep(summary.matrix[, 5 + m + i] < 0.999999999999, L), L, L, TRUE) * 
              upper.tri(matrix(, L, L)))
      violation.matrix[1:(m - 1), 1] <- nac
      violation.matrix[m, 1] <- sum(nac)
      freq <- summary.matrix[, 5:(m + 4)] 
      for (i in 1:(m - 1)) {
        V <- outer(summary.matrix[, (m + 5 + i)], summary.matrix[, 
                                                                 (m + 5 + i)], "-")
        V[row(V) <= col(V)] <- 0
        V[V >= -minvi] <- 0
        violation.matrix[i, 2] <- sum(ceiling(abs(V)))
        violation.matrix[i, 4] <- max(abs(V))
        if (violation.matrix[i, 4] > minvi) {
          violation.matrix[i, 5] <- sum(abs(V))
          freqd <- cbind(apply(as.matrix(freq[, 1:i]), 
                               1, sum), apply(as.matrix(freq[, (i + 1):m]), 
                                              1, sum))
          #freqd <- freqd * Rs 
          Z <- abs(sign(-V) * 2 * (sqrt(outer(freqd[, 2] + 
                                                1, freqd[, 1] + 1)) - sqrt(outer(freqd[, 1], 
                                                                                 freqd[, 2])))/sqrt(outer(freqd[, 2], freqd[, 
                                                                                                                            1], "+") + outer(freqd[, 1], freqd[, 2], "+") - 
                                                                                                      1))
          violation.matrix[i, 7] <- max(Z)
          violation.matrix[i, 8] <- min(col(Z)[Z == max(Z)])
          violation.matrix[i, 9] <- min(row(Z)[Z == max(Z)])
          violation.matrix[i, 10] <- sum(sign(Z[Z > 1.6449]))
        }
      }
      violation.matrix[m, 2] <- sum(violation.matrix[1:(m - 
                                                          1), 2])
      violation.matrix[1:m, 3] <- violation.matrix[1:m, 2]/violation.matrix[1:m, 
                                                                            1]
      violation.matrix[m, 4] <- max(violation.matrix[1:(m - 
                                                          1), 4])
      violation.matrix[m, 5] <- sum(violation.matrix[1:(m - 
                                                          1), 5])
      violation.matrix[1:m, 6] <- violation.matrix[1:m, 5]/violation.matrix[1:m, 
                                                                            1]
      violation.matrix[m, 7] <- max(violation.matrix[1:(m - 
                                                          1), 7])
      violation.matrix[m, 10] <- sum(violation.matrix[1:(m - 
                                                           1), 10])
      resultsa[[j]][[3]] <- violation.matrix
      resultsa[[j]][[4]] <- paste("Minsize = ", minsize, " Minvi = ", 
                                  minvi, sep = "")
      
      
    }
    
    Hi <- MLcoefH(cbind(level.two.var, X), se = F, nice.output = F)$Hi
    
    monotonicity.list <- list(results = results, I.labels = I.labels, 
                              Hi = Hi[, 1], m = m, X = X)
    
    monotonicity.lista <- list(results = resultsa, I.labels = I.labels, 
                               Hi = Hi[, 2], m = m, X = Xa)
    
    monotonicity.list <- list(monotonicity.list,
                              monotonicity.lista)
  } else {
    Hi <- coefHTiny(X)$Hi
    monotonicity.list <- list(results = results, I.labels = I.labels, Hi = Hi, m = m, X = X)
  }
  class(monotonicity.list) <- "monotonicity.class"
  
 return(monotonicity.list)
}
