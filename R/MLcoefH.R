######## Multilevel Mokken scale analysis - MLweight and MLcoefH
######## Letty Koopman 
######## University of Amsterdam
######## March 23, 2017

# Used functions from the mokken package:
"all.patterns" <- function(J,m){
  grid <- list()
  j <- 0;
  p <- m^J
  for (j in 1:J){
    grid <- c(grid, j)
    grid[[j]] <- 0:(m-1)
  }
  X <- t(expand.grid(grid))
  dimnames(X) <- NULL
  return(X[J:1,])
}

"phi" <- function(A,f, action){
  # Numerical values are translations h(A %*% f) = A %*% f -
  eps = 1E-80;
  switch(action,
         "identity" = A %*% f,
         "exp"      = A %*% exp(f),
         "log"      = A %*% log(abs(f)+eps),
         "sqrt"     = A %*% sqrt(f),
         "xlogx"    = A %*% (-f*log(f+eps)),
         "xbarx"    = A %*% (f*(1-f))  # x(1-x)
  )
}

"dphi" <- function(A,f,df, action){
  eps=1E-80;
  switch(action,
         "identity" = A %*% df,
         "exp"      = A %*% (as.numeric(exp(f)) * df),
         "log"      = A %*% (as.numeric(1/(f+eps)) * df),
         "sqrt"     = A %*% (as.numeric(1/(2*sqrt(f))) * df),
         "xlogx"    = A %*% (as.numeric(-1-log(f+eps)) * df),
         "xbarx"    = A %*% (as.numeric(1-2*f) * df)  #, x(1-x)
  )
}

"check.data" <- function(X){
  if (data.class(X) != "matrix" && data.class(X) != "data.frame")
    stop("Data are not matrix or data.frame")
  matrix.X <- as.matrix(X)
  if (any(is.na(matrix.X))) stop("Missing values are not allowed")
  if (mode(matrix.X)!="numeric") stop("Data must be numeric")
  if (any(matrix.X < 0)) stop("All scores should be nonnegative")
  if (any(matrix.X %% 1 !=0)) stop("All scores must be integers")
  matrix.X <- matrix.X - min(matrix.X)
  return(matrix.X)
}


## New functions to be implemented in mokken package
"MLweight" <- function(X, maxx = NULL, minx = NULL){
  # Computes the two-level Guttman weights for two-level Mokken Scale Analysis.
  #
  # Args:
  #   X: Data matrix with a subject column and two item columns. The subject column is assumed to be the first.
  #   maxx: The highest possible answer category. If not specified it is determined by using the lowest item score.
  #   minx: The lowest possible answer category. If not specified it is determined by using the highest item score.
  #   Depends on "all.patterns".
  #
  # Returns:
  #   Guttman weights for two-level data.
  
  # Error handling:
  if (ncol(X) != 3){
    warning('X contains more than two items. Only first two items will be used')
    X <- X[, 1:3]
  }
  if(is.null(minx)) minx <- min(X[, -1])
  if(is.null(maxx)) maxx <- max(X[, -1])
  X[, -1] <- X[, -1] - minx
  maxx <- maxx - minx
  minx <- 0
  s <- X[, 1] # Subject column
  f1 <- factor(X[, 2], levels = minx:maxx) # First item 
  f2 <- factor(X[, 3], levels = minx:maxx) # Second item
  
  # Compute relative frequencies
  Rel1 <- colSums(table(s, f1) / rowSums(table(s, f1))) 
  names(Rel1) <- paste0(1, minx:maxx)
  Rel2 <- colSums(table(s, f2) / rowSums(table(s, f2)))
  names(Rel2) <- paste0(2, minx:maxx)
  
  # Cumulative relative frequencies, dealing with equal ranks
  CumRel <- c(rev(cumsum(rev(Rel1[-1]))), rev(cumsum(rev(Rel2[-1]))))
  names(CumRel) <- 1:length(CumRel)
  y <- sort(CumRel, decreasing=T)
  
  if (any(duplicated(y))) {
    perm <- function (n, r, v = 1:n, set = TRUE){
      if (r == 1)  
        matrix(v, n, 1)
      else if (n == 1) 
        matrix(v, 1, r)
      else {
        X <- NULL
        for (i in 1:n) X <- rbind(X, cbind(v[i], Recall(n - 1, r - 1, v[-i])))
        X
      }
    }
    o <- lapply(unique(y), function(val) { 
      m <- as.numeric(names(y[y==val]))
      if (length(m) <= 1) {
        return(m)
      }
      # do.call(paste0, as.data.frame(perm(length(m), length(m), m)))
      as.data.frame(perm(length(m), length(m), m))
    })
    
    # Ensure fixed ordering within items
    for (i in 1:length(o)){
      g <- o[[i]]
      if(length(unique(g)) > 1){
        select <- matrix(0, nrow(g))
        for(j in 1:nrow(g)){
          h <- g[j, ]
          if(any(h >= 1 & h <= maxx) & any(h >= maxx + 1 & h <= maxx * 2)){
            select[j] <- (all(h[which(h >= 1 & h <= maxx)] == sort(h[which(h >= 1 & h <= maxx)])) & all(h[which(h >= maxx + 1 & h <= maxx * 2)] == sort(h[which(h >= maxx + 1 & h <= maxx * 2)]))) * j
          } else {
            i1 <- ifelse(any(h >= 1 & h <= maxx), all(h[which(h >= 1 & h <= maxx)] == sort(h[which(h >= 1 & h <= maxx)])), 0)
            i2 <- ifelse(any(h >= maxx + 1 & h <= maxx * 2), 
                         all(h[which(h >= maxx + 1 & h <= maxx * 2)] == sort(h[which(h >= maxx + 1 & h <= maxx * 2)])), 0)
            select[j] <- (i1 + i2) * j
          }
        }
        o[[i]] <- matrix(apply(o[[i]][select, ], 1, paste0, collapse = ""))
      }
    }
    
    out <- matrix(do.call(paste0, expand.grid(o)))
    w <- Z <- NULL
    for(i in 1:nrow(out)){
      ords <- as.numeric(unlist(strsplit(out[i, ],NULL)))
      # Compute Z matrix for each item-response pattern
      Z <- matrix(rep(matrix(all.patterns(2, maxx + 1), nrow = 1), maxx), nrow = maxx, byrow = TRUE)
      Z <- matrix(ifelse(Z < row(Z), 0, 1), ncol = (maxx) * 2, byrow = TRUE)
      Z <- Z[, ords]
      
      # Compute weights
      w <- rbind(w, apply(Z, 1, function(x){sum(x * cumsum(abs(x - 1)))}))
    }
    wr <- matrix(colMeans(w), nrow = 1)
  } else {
    out <- paste(names(y), collapse = "")
    
    ords <- as.numeric(unlist(strsplit(out,NULL)))
    # Compute Z matrix for each item-response pattern
    Z <- matrix(rep(matrix(all.patterns(2, maxx + 1), nrow = 1), maxx), nrow = maxx, byrow = TRUE)
    Z <- matrix(ifelse(Z < row(Z), 0, 1), ncol = (maxx) * 2, byrow = TRUE)
    Z <- Z[, ords]
    
    # Compute weights
    wr <- apply(Z, 1, function(x){sum(x * cumsum(abs(x - 1)))})
    
  }
  return(wr)
}

"MLcoefH" <- function(X, se = TRUE, nice.output = TRUE, subject = 1){
  # Computes the two-level scalability coefficients in Mokken scale analysis
  #
  # Args:
  #   X: Data matrix with a subject column and one column per item. Preferably the subject column consists of integers.
  #   se: If TRUE, computes the standard errors for the coefficients, 
  #       if FALSE, only the coefficients are computed. Default is TRUE.
  #   nice.output: If TRUE, prints the coefficients and standard errors in a matrix with nice lay-out,
  #                if FALSE, they are printed in a regular type matrix which can be used for further computations. Default is TRUE.
  #   Subject: Represents the subject column. Default is column 1. 
  # 
  # Depends on package "data.table" and mokken functions "MLweight", "check.data", "all.patterns", "phi", and "dphi".
  #
  # Returns: 
  #   Two-level scalability coefficients and optionally their standard errors.
  
  # Error handling:
  if(subject != 1){
    X <- cbind(X[, subject], X[, -subject])
  }
  X <- X[order(X[, 1]), ] # Order the data according to S.
  Rs <- as.numeric(table(X[, 1]))
  LS <- length(Rs)
  S <- 1:LS 
  X[, 1] <- rep(S, Rs)
  X <- check.data(X) 
  X[, 1] <- rep(S, Rs) # make sure subject column runs from 1 to S. 
  if(is.null(colnames(X))) colnames(X) <- c("Subs", paste("Item", 1:(ncol(X) - 1)))
  
  if(any(Rs == 1)){ 
    warning('For at least one subject there is only 1 rater. The scalability coefficients are computed without this (these) subject(s).') 
    X <- X[!(X[, 1] %in% which(Rs == 1)), ]
    Rs <- as.numeric(table(X[, 1]))
    LS <- length(Rs)
    S <- 1:LS
  }
  
  X <- X[do.call(order, lapply(1:NCOL(X), function(i) X[, i])), ]
  
  labels <- dimnames(X[, -1])[[2]]
  m <- max(X[, -1]) 
  J <- ncol(X[, -1])
  K <- choose(J, 2) 
  g <- m + 1
  
  nams <- apply(combn(colnames(X)[-1],2), 2, function(z) paste(z, collapse = ' ')) 
  cols <- combn(J, 2)
  Patterns <- cbind("Xa" = rep(0:m, each = m + 1), "Xb" = rep(0:m, m + 1)) 
  
  
  
  if(se == TRUE){
    R <- unique(X)
    n <- as.numeric(table(factor(apply(X, 1, paste, collapse=","), levels=unique(apply(X, 1, paste, collapse=",")))))
    Rss <- table(R[, 1])
    Rd <- rep(Rs, Rss)
    
    G5 <- matrix(0, 3 * K, nrow(R))
    g5 <- matrix(0, K, 3)
    ABs <- AWs <- AEs <- matrix(0, g^2, nrow(R))
    Fw <- Fb <- Fe <- NULL
      
    for(k in 1:K){
      z <- cols[, k]
      Subs <- R[, 1]
      Xa <- R[, z[1] + 1]
      Xb <- R[, z[2] + 1]
      Weights <- MLweight(X[, c(1, z + 1)], minx = 0, maxx = m)# Weights <- MLweight(Rsub, minx = 0, maxx = m)
      for(x in 1:g^2){
        if(Weights[x] > 0){
          i <- Patterns[x, 1]
          j <- Patterns[x, 2]
          nw <- tapply((Xa == i & Xb == j) * n, Subs, sum)
          ni <- tapply((Xa == i) * n, Subs, sum)
          Ni <- sum(ni / Rs)
          nj <- tapply((Xb == j) * n, Subs, sum)
          Nj <- sum(nj / Rs)
          at <- (Xa == i) * (rep(nj, Rss) - (Xb == j))
          nb <- tapply(at * n, Subs, sum)
          
          Fw[x] <- Weights[x] * sum(nw / Rs)
          Fb[x] <- Weights[x] * sum(nb / (Rs * (Rs - 1)))
          Fe[x] <- Weights[x] * Ni * Nj / LS
          
          cat <- as.numeric(Xa == i) + as.numeric(Xb == j)
          ABs[x, ] <- Weights[x] * (at * Rd - rep(nb, Rss)) / (Rd^2 * (Rd - 1))
          AWs[x, ] <- Weights[x] * ifelse(cat == 2, (Rd - rep(nw, Rss)) / Rd^2, -rep(nw, Rss) / Rd^2)
          AEs[x, ] <- Weights[x] * ifelse(cat == 2, (Rd * (Ni + Nj) - Ni * rep(nj, Rss) - Nj * rep(ni, Rss)) / (Rd^2 * LS), 
                                           ifelse(cat == 1 & Xa == i, -(Ni * rep(nj, Rss) + Nj * (rep(ni, Rss) - Rd)) / (Rd^2 * LS), 
                                                  ifelse(cat == 1 & Xb == j, -(Nj * rep(ni, Rss) + Ni * (rep(nj, Rss) - Rd)) / (Rd^2 * LS),
                                                         (-Ni * rep(nj, Rss) - Nj * rep(ni, Rss)) / (Rd^2 * LS))))
         
        } else {
          Fw[x] <- Fb[x] <- Fe[x] <- 0
          ABs[x, ] <- AWs[x, ] <- AEs[x, ] <- 0
        }
      }
      g5[k, 1] <- sum(Fb)
      g5[k, 2] <- sum(Fw)
      g5[k, 3] <- sum(Fe)
      
      G5[k, ] <- colSums(ABs)
      G5[k + K, ] <- colSums(AWs)
      G5[k + 2 * K, ] <- colSums(AEs)
    }

    g5i <- matrix(0, J, 3)
    for(i in 1:J) {
      g5i[i, ] <- colSums(g5[apply(cols, 2, function(x) any(x == i)), ])
    }
    g5i <- rbind(g5i[1, 1], matrix(g5i))
    
    g5ii <- colSums(g5)
    g5ii <- rbind(g5ii[1], matrix(g5ii))
    
    g5 <- rbind(g5[1, 1], matrix(g5))
   
    # for item coefficients
    G5i <- matrix(0, J * 3, nrow(R))
    for(i in 1:J) {
      items <- which(apply(cols, 2, function(x) any(x == i)))
      G5i[i, ] <- colSums(G5[items, ])
      G5i[i + J, ] <- colSums(G5[items + K, ])
      G5i[i + J * 2, ] <- colSums(G5[items + K * 2, ])
    }
    G5i <- rbind(G5i[1, ], G5i)
    
    # for total scale coefficients
    G5ii <- rbind(colSums(G5[1:K, ]), colSums(G5[1:K, ]), 
                  colSums(G5[(K + 1):(2 * K), ]), colSums(G5[(2 * K + 1):(3 * K), ]))

    G5 <- rbind(G5[1, ], G5)
    
    # Create A6 --> To compute the ratio of observed to expected errors
    A6 <- rbind(c(1, -1, rep(0, 3 * K - 1)), 
                cbind(rep(0, 2 * K), diag(2 * K), 
                      rbind(-diag(K), -diag(K))))
    
    A6i <- rbind(c(1, -1, rep(0, 3 * J - 1)), 
                 cbind(rep(0, 2 * J), diag(2 * J), 
                       rbind(-diag(J), -diag(J))))
    
    A6ii <- matrix(c(1, 0, 0, -1, 1, 0, 0, 0, 1, 0, -1, -1), 3)
    
    # Create A7 --> To compute the Hij/Hi/H values
    A7 <- cbind(rep(1, 2 * K), -diag(2 * K))
    
    A7i <- cbind(rep(1, 2 * J), -diag(2 * J))
    
    A7ii <- matrix(c(1, 1, -1, 0, 0, -1), 2)
    
    # Create A8 and A9 --> To compute ratio HB/HW
    A8 <- cbind(diag(K), -diag(K))
    A8i <- cbind(diag(J), -diag(J))
    A8ii <- matrix(c(1, -1), 1)
    
    A9 <- diag(K)
    A9i <- diag(J)
    A9ii <- matrix(1)
    
    #Hij
    
    g6 <- phi(A6, g5, "log")
    G6 <- dphi(A6, g5, G5, "log")
    
    g7 <- phi(A7, g6, "exp") 
    G7 <- dphi(A7, g6, G6, "exp")
    
    g8 <- phi(A8, g7, "log") 
    G8 <- dphi(A8, g7, G7, "log")
    
    G9 <- dphi(A9, g8, G8, "exp")
    
    HBij <- g7[1:K, ]
    HWij <- g7[-c(1:K), ]
    HBWij <- HBij / HWij
    
    se.Hij <- sqrt(diag(G7 %*% (as.numeric(n) * t(G7))))
    se.HBij <- se.Hij[1:K]
    se.HWij <- se.Hij[-c(1:K)]
    se.HBWij <- sqrt(diag(G9 %*% (as.numeric(n) * t(G9))))
    
    # Hi
    
    g6 <- phi(A6i, g5i, "log") 
    G6 <- dphi(A6i, g5i, G5i, "log")
    
    g7 <- phi(A7i, g6, "exp")
    G7 <- dphi(A7i, g6, G6, "exp")
    
    g8 <- phi(A8i, g7, "log") 
    G8 <- dphi(A8i, g7, G7, "log")
    
    g9 <- phi(A9i, g8, "exp") 
    G9 <- dphi(A9i, g8, G8, "exp")
    
    HBi <- g7[1:J, ]
    HWi <- g7[-c(1:J), ]
    HBWi <- HBi / HWi
    
    se.Hi <- sqrt(diag(G7 %*% (as.numeric(n) * t(G7))))
    se.HBi <- se.Hi[1:J]
    se.HWi <- se.Hi[-c(1:J)]
    se.HBWi <- sqrt(diag(G9 %*% (as.numeric(n) * t(G9))))
    
    # H
    
    g6 <- phi(A6ii, g5ii, "log")
    G6 <- dphi(A6ii, g5ii, G5ii, "log")
    
    g7 <- phi(A7ii, g6, "exp")
    G7 <- dphi(A7ii, g6, G6, "exp")
    
    g8 <- phi(A8ii, g7, "log") 
    G8 <- dphi(A8ii, g7, G7, "log")
    
    g9 <- phi(A9ii, g8, "exp") 
    G9 <- dphi(A9ii, g8, G8, "exp")
    
    HB <- g7[1, ]
    HW <- g7[2, ]
    HBW <- HB / HW
    se.H <- sqrt(diag(G7 %*% (as.numeric(n) * t(G7))))
    se.HB <- se.H[1]
    se.HW <- se.H[2]
    se.HBW <- sqrt(diag(G9 %*% (as.numeric(n) * t(G9))))
    
    if(nice.output == TRUE){
      Hij <- HBijt <- matrix(0, J, J, dimnames = list(labels, labels))
      Hij[lower.tri(Hij)] <- HWij
      HBijt[lower.tri(HBijt)] <- HBij
      Hij <- Hij + t(HBijt)
      
      se.Hij <- se.Hijt <- matrix(0, J, J, dimnames = list(labels, labels))
      se.Hij[lower.tri(se.Hij)] <- se.HWij
      se.Hijt[lower.tri(se.Hijt)] <- se.HBij
      se.Hij <- se.Hij + t(se.Hijt)
      
      new.labels <- rep(labels, each = 2)
      new.labels[2 * 1:J] <- "(se)"
      OM.Hij <- matrix(NA, J + 3, J * 2 + 1)
      for (j in 2 * (1:J)) {
        OM.Hij[, j ] <- c("", "", "", format(paste(" ", formatC(round(Hij[, j/2], 3), digits = 3, format = "f"), " ", sep = ""), width = 7, justify = "right"))
        OM.Hij[, j + 1] <- c("", "", "", format(paste("(", formatC(round(se.Hij[, j/2], 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right"))
      }
      OM.Hij[, 1] <- OM.Hij[-c(1:3), -1][row(OM.Hij[-c(1:3), -1]) == (0.5 * col(OM.Hij[-c(1:3), -1]) + 0.5)] <- format("", width = 7, justify = "right")
      OM.Hij[-c(1:3), -1][row(OM.Hij[-c(1:3), -1]) == (0.5 * col(OM.Hij[-c(1:3), -1]))] <- format("", width = 7, justify = "right")
      OM.Hij[round(J / 2) + 3, 1] <- format("(HWij)", width = 7, justify = "centre")
      OM.Hij[2, round(J / 2) * 2] <- format("(HBij)", width = 7, justify = "centre")
      rownames(OM.Hij) <- c("", "", "", labels)
      colnames(OM.Hij) <- c("", new.labels)
      OM.Hij <- noquote(OM.Hij)
      
      # HWi & HBi
      OM.Hi <- matrix(NA, J, 7)
      OM.Hi[, 1] <- format("", width = 7, justify = "right")
      OM.Hi[, 2] <- format(formatC(round(HWi, 3), digits = 3, format = "f"), width = 7, justify = "right")
      OM.Hi[, 3] <- format(paste("(", formatC(round(se.HWi, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
      OM.Hi[, 4] <- format(formatC(round(HBi, 3), digits = 3, format = "f"), width = 7, justify = "right")
      OM.Hi[, 5] <- format(paste("(", formatC(round(se.HBi, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
      OM.Hi[, 6] <- format(formatC(round(HBWi, 3), digits = 3, format = "f"), width = 7, justify = "right")
      OM.Hi[, 7] <- format(paste("(", formatC(round(se.HBWi, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
      dimnames(OM.Hi) <- list(labels, c("", "   HWi", "  (se)  ", "   HBi", "  (se)  ", "   BWi", "  (se)  "))
      OM.Hi <- noquote(OM.Hi)
      
      # HW & HB
      OM.H <- matrix(NA, 1, 7)
      OM.H[, 1] <- format("", width = 7, justify = "right")
      OM.H[, 2] <- format(formatC(round(HW, 3), digits = 3, format = "f"), width = 7, justify = "right")
      OM.H[, 3] <- format(paste("(", formatC(round(se.HW, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
      OM.H[, 4] <- format(formatC(round(HB, 3), digits = 3, format = "f"), width = 7, justify = "right")
      OM.H[, 5] <- format(paste("(", formatC(round(se.HB, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
      OM.H[, 6] <- format(formatC(round(HBW, 3), digits = 3, format = "f"), width = 7, justify = "right")
      OM.H[, 7] <- format(paste("(", formatC(round(se.HBW, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
      dimnames(OM.H) <- list("Scale", c(" ", "   HW", "  (se)  ", "   HB", "  (se)  ", "   BW", "  (se)  "))
      OM.H <- noquote(OM.H)
      
      # Output:
      OL <- list(Hij = OM.Hij, Hi = OM.Hi, H = OM.H)
      
    } else {
      Hij <- data.frame(HWij, se.HWij, HBij, se.HBij, HBWij,  se.HBWij) 
      rownames(Hij) <- nams
      Hi <- data.frame(HWi, se.HWi, HBi = HBi, se.HBi, HBWi, se.HBWi) 
      rownames(Hi) <- labels
      H <- data.frame(HW, se.HW, HB, se.HB, HBW, se.HBW)
      OL <- list(Hij = Hij, Hi = Hi, H = H)
    }
    
  } else {
    Fwt <- Fbt <- Fet <- Fw <- Fb <- Fe <- NULL
    for(k in 1:K){
      z <- cols[, k]
      Subs <- X[, 1]
      Xa <- X[, z[1] + 1]
      Xb <- X[, z[2] + 1]
      Rss <- rep(Rs, Rs)
      Weights <- MLweight(X[, c(1, z + 1)], minx = 0, maxx = m)
      
      for(x in 1:g^2){
        if(Weights[x] > 0){
          i <- Patterns[x, 1]
          j <- Patterns[x, 2]
          nw <- sum((Xa == i & Xb == j) / Rss)
          ni <- sum((Xa == i) / Rss)
          nj <- tapply((Xb == j), Subs, sum)
          at <- (Xa == i) * (rep(nj, Rs) - (Xb == j))
          nb <- tapply(at, Subs, sum)
          nj <- sum(nj / Rs)
          
          Fwt[x] <- Weights[x] * nw
          Fbt[x] <- Weights[x] * sum(nb / (Rs * (Rs - 1)))
          Fet[x] <- Weights[x] * ni * nj / LS
        } else {
          Fwt[x] <- Fbt[x] <- Fet[x] <- 0
        }
      } 
      Fw[k] <- sum(Fwt)
      Fb[k] <- sum(Fbt)
      Fe[k] <- sum(Fet)
    }
    
    Fwi <- Fbi <- Fei <- NULL
    for(i in 1:J) {
      items <- apply(cols, 2, function(x) any(x == i))
      Fwi[i] <- sum(Fw[items])
      Fbi[i] <- sum(Fb[items])
      Fei[i] <- sum(Fe[items])
    }
    
    HBij <- 1 - Fb / Fe
    HWij <- 1 - Fw / Fe
    HBWij <- HBij / HWij
    
    HBi <- 1 - Fbi / Fei
    HWi <- 1 - Fwi / Fei
    HBWi <- HBi / HWi
    
    HB <- 1 - sum(Fb) / sum(Fe)
    HW <- 1 - sum(Fw) / sum(Fe)
    HBW <- HB / HW
    
    if(nice.output == TRUE){
      Hij <- HBijt <- matrix(0, J, J, dimnames = list(labels, labels))
      Hij[lower.tri(Hij)] <- HWij
      HBijt[lower.tri(Hij)] <- HBij
      Hij <- Hij + t(HBijt)
      
      OM.Hij <- matrix(NA, J + 3, J + 1)
      for (j in (1:J)) {
        OM.Hij[, j + 1] <- c("", "", "", format(paste(" ", formatC(round(Hij[, j], 3), digits = 3, format = "f"), " ", sep = ""), width = 7, justify = "right"))
      }
      OM.Hij[, 1] <- OM.Hij[-c(1:2), ][row(OM.Hij[-c(1:2), ]) == col(OM.Hij[-c(1:2), ])] <- format("", width = 7, justify = "right")
      OM.Hij[round(J / 2) + 3, 1] <- format("(HWij)", width = 7, justify = "centre")
      OM.Hij[2, round(J / 2) + 1] <- format("(HBij)", width = 7, justify = "centre")
      rownames(OM.Hij) <- c("", "", "", labels)
      colnames(OM.Hij) <- c("", labels)
      OM.Hij <- noquote(OM.Hij)
      
      # HWi & HBi
      OM.Hi <- matrix(NA, J, 4)
      OM.Hi[, 1] <- format("", width = 7, justify = "right")
      OM.Hi[, 2] <- format(formatC(round(HWi, 3), digits = 3, format = "f"), width = 7, justify = "right")
      OM.Hi[, 3] <- format(formatC(round(HBi, 3), digits = 3, format = "f"), width = 7, justify = "right")
      OM.Hi[, 4] <- format(formatC(round(HBWi, 3), digits = 3, format = "f"), width = 7, justify = "right")
      dimnames(OM.Hi) <- list(labels, c("", "   HWi", "   HBi", "   BWi"))
      OM.Hi <- noquote(OM.Hi)
      
      # HW & HB
      OM.H <- matrix(NA, 1, 4)
      OM.H[, 1] <- format("", width = 7, justify = "right")
      OM.H[, 2] <- format(formatC(round(HW, 3), digits = 3, format = "f"), width = 7, justify = "right")
      OM.H[, 3] <- format(formatC(round(HB, 3), digits = 3, format = "f"), width = 7, justify = "right")
      OM.H[, 4] <- format(formatC(round(HBW, 3), digits = 3, format = "f"), width = 7, justify = "right")
      dimnames(OM.H) <- list("Scale", c(" ", "   HW", "   HB", "   BW"))
      OM.H <- noquote(OM.H)
      
      # Output:
      OL <- list(Hij = OM.Hij, Hi = OM.Hi, H = OM.H)
      
    } else {
      Hij <- data.frame(HWij, HBij, HBWij) 
      rownames(Hij) <- nams
      Hi <- data.frame(HWi, HBi = HBi, HBWi) 
      rownames(Hi) <- labels
      H <- data.frame(HW, HB, HBW)
      OL <- list(Hij = Hij, Hi = Hi, H = H)
    } 
  }
  

  return(OL)
  
  }
