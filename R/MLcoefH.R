######## Letty Koopman 
######## University of Amsterdam
# MLcoefH() version 18-9-2017
#
# Updates: 
## Computes proportions instead of relative frequencies.
## Can handle only two items
## Additional error handling: handling duplicate column names. 

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
  eps <- 1e-40
  X <- X[order(X[, 1]), ] # Order the data according to S.
  Rs <- as.numeric(table(X[, 1]))
  LS <- length(Rs)
  S <- 1:LS 
  X[, 1] <- rep(S, Rs)
  X <- check.data(X) 
  X[, 1] <- rep(S, Rs) # make sure subject column runs from 1 to S. --> deze hoeft wellicht niet meer? 
  if(is.null(colnames(X))) colnames(X) <- c("Subs", paste("Item", 1:(ncol(X) - 1)))
  if(any(table(colnames(X)) > 1)){
    warning("At least two items have the same name.")
    colnames(X) <- make.unique(colnames(X))
  }
  if (any(apply(X, 2, var) < eps)) 
    stop("One or more variables have zero variance. The duplicates have received an index.")
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
  B <- K * g^2 
  U <- J * g 
  
  nams <- apply(combn(colnames(X)[-1],2), 2, paste, collapse = ' ') 
  cols <- combn(J, 2)
  Patterns <- cbind("Xa" = rep(0:m, each = m + 1), "Xb" = rep(0:m, m + 1)) 
  
  
  if(se == TRUE){
    R <- unique(X)
    Rss <- rep(Rs, Rs)
    RRs <- table(R[, 1])
    Rd <- rep(Rs, RRs)
    n <- as.numeric(table(factor(apply(X, 1, paste, collapse=","), levels=unique(apply(X, 1, paste, collapse=","))))) 
    SubsR <- R[, 1]
    Rt <- aggregate(n ~ ., data = R[, -1], sum) 
    np <- Rt[, J + 1] 
    Rt <- Rt[, -(J + 1)]
    Lst <- length(np)
    SubsX <- X[, 1]
    nprel <- aggregate(n / (Rd * LS) ~ ., data = R[, -1], sum)[, J + 1]
    nS <- aggregate(SubsX ~ ., X, unique)[, J + 1]
    nN <- cbind(R, rownr = 1:nrow(R))
    nN <- aggregate(rownr ~ ., nN[, -1], unique)[, J + 1]
    
    # Covariance matrix for nested data
    patsX <- apply(X[, -1], 1, paste, collapse = "")
    patsR <- apply(Rt, 1, paste, collapse = "")
    pSubs <- matrix(0, LS, Lst)
    for(s in S) {
      pSubs[s, ] <- rowSums(outer(patsR, patsX[SubsX == s], "==")) / Rs[s]
    }
    covps <- cov(pSubs)
    p <- colMeans(pSubs)
    covp <- (diag(length(p)) * p - (p %*% t(p))) 
    nu <- mean(Rs)
    
    # Variance covariance matrix needed for delta method 
    covtot <- LS * nu * (nu - 1) * covps + LS * nu * covp
    
    # Creating g3 and G3
    ns <- list()
    nuni <- matrix(0, J, g)
    for(i in 1:J) {
      Xa <- R[, i + 1]
      ns[[i]] <- matrix(0, g, LS)
      for(a in 0:m){
        nst <- tapply((Xa == a) * n, SubsR, sum)
        ns[[i]][a + 1, ] <- nst
        nuni[i, a + 1] <- sum(nst / (Rs * LS))
      }
    }
    
    
    Fwt <- Fbt <- Fet <- Fw <- Fb <- Fe <- eij <- NULL
    G3W <- G3B <- G3E <- matrix(0, K, length(np))
    for(k in 1:K){
      z <- cols[, k]
      Ra <- R[, z[1] + 1]
      Rb <- R[, z[2] + 1]
      Ta <- Rt[, z[1]]
      Tb <- Rt[, z[2]]
      Weights <- MLweight(X[, c(1, z + 1)], minx = 0, maxx = m)
      Wmat <- Bmat <- Emat <- matrix(0, g^2, length(np))
      
      for(x in 1:g^2){
        if(Weights[x] > 0){
          i <- Patterns[x, 1]
          j <- Patterns[x, 2]
          
          temp <- (Ra == i) * (rep(ns[[z[2]]][j + 1, ], RRs) - (Rb == j)) * n
          nw <- sum(n * (Ra == i & Rb == j) / (rep(Rs, RRs) * LS)) 
          at <- temp / (rep(Rs, RRs) * (rep(Rs, RRs) - 1) * LS) 
          at <- sapply(1:length(np), function(x) sum(at[nN[[x]]]))
          
          Wmat[x, ] <- ((Ta == i & Tb == j) * nprel / np) * Weights[x]
          Bmat[x, ] <- (at / np) * Weights[x] 
          Emat[x, ] <- Weights[x] / np * ((Tb == j) * nuni[z[[1]], i + 1] + (Ta == i) * nuni[z[[2]], j + 1]) * nprel
          Fwt[x] <- Weights[x] * nw
          Fbt[x] <- Weights[x] * sum(tapply(temp, SubsR, sum) / (Rs * (Rs - 1) * LS)) 
          Fet[x] <- Weights[x] * nuni[z[[1]], i + 1] * nuni[z[[2]], j + 1] 
        } else {
          Fwt[x] <- Fbt[x] <- Fet[x] <- 0
        }
      } 
      
      G3W[k, ] <- colSums(Wmat)
      G3B[k, ] <- colSums(Bmat)
      G3E[k, ] <- colSums(Emat)
      Fw[k] <- sum(Fwt)
      Fb[k] <- sum(Fbt)
      Fe[k] <- sum(Fet)
    }
    
    Fwi <- Fbi <- Fei <- NULL
    G3Wi <- G3Bi <- G3Ei <- matrix(0, J, length(np))
    for(i in 1:J) {
      items <- apply(cols, 2, function(x) any(x == i))
      Fwi[i] <- sum(Fw[items])
      Fbi[i] <- sum(Fb[items])
      Fei[i] <- sum(Fe[items])
      if(J > 2){
        G3Wi[i, ] <- colSums(G3W[items, ])
        G3Bi[i, ] <- colSums(G3B[items, ])
        G3Ei[i, ] <- colSums(G3E[items, ])
      } else {
        G3Wi[i, ] <- (G3W[items, ])
        G3Bi[i, ] <- (G3B[items, ])
        G3Ei[i, ] <- (G3E[items, ])
      }
    }
    
    g3 <- matrix(c(Fb[1], Fb, Fw, Fe))
    
    g3i <- matrix(c(Fbi[1], Fbi, Fwi, Fei))
    
    g3ii <- matrix(c(sum(Fb), sum(Fb), sum(Fw), sum(Fe)))
    
    # Jacobian Hij, Hi and H coefficients
    G3 <- rbind(G3B[1, ], G3B, G3W, G3E)
    
    G3i <- rbind(G3Bi[1, ], G3Bi, G3Wi, G3Ei)
    
    G3ii <- rbind(colSums(G3B), colSums(G3B), colSums(G3W), colSums(G3E))
    
    # Create A6 --> To compute the ratio of observed to expected errors
    A6 <- rbind(matrix(c(1, -1, rep(0, 3 * K - 1)), 1), 
                cbind(matrix(0, K * 2, 1), diag(K * 2), rbind(-1 * diag(K), -1 * diag(K))))
    A6i <- rbind(matrix(c(1, -1, rep(0, 3 * J - 1)), 1), 
                 cbind(matrix(0, J * 2, 1), diag(J * 2), rbind(-1 * diag(J), -1 * diag(J))))
    A6ii <- rbind(matrix(c(1, -1, rep(0, 3 - 1)), 1), 
                  cbind(matrix(0, 2, 1), diag(2), rbind(-1, -1)))
    
    # Create A7 --> To compute the Hij/Hi/H values
    A7 <- rbind(cbind(rep(1, 2 * K), -diag(2 * K)),
                cbind(rep(1, 2 * K), -diag(2 * K)))
    A7i <- rbind(cbind(rep(1, 2 * J), -diag(2 * J)),
                 cbind(rep(1, 2 * J), -diag(2 * J)))
    A7ii <- matrix(c(1, 1, 1, 1, -1, 0, -1, 0, 0, -1, 0, -1), 4)
    
    # Create A8 and A9 --> To compute ratio HB/HW
    A8 <- cbind(diag(K * 3), rbind(matrix(0, 2 * K, K), -diag(K)))
    A8i <- cbind(diag(J * 3), rbind(matrix(0, 2 * J, J), -diag(J)))
    A8ii <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, -1), 3)
    
    A9 <- diag(K * 3)
    A9i <- diag(J * 3)
    A9ii <- diag(3)
    
    # Hij
    
    g4 <- phi(A6, g3, "log")
    G4 <- dphi(A6, g3, G3, "log")
    
    g5 <- phi(A7, g4, "exp")
    G5 <- dphi(A7, g4, G4, "exp")
    
    g6 <- phi(A8, g5, "log")
    G6 <- dphi(A8, g5, G5, "log")
    
    g7 <- phi(A9, g6, "exp")
    G7 <- dphi(A9, g6, G6, "exp")
    
    se.Hij <- sqrt(diag(G7 %*% (covtot %*% t(G7))))
    
    HBij <- g5[1:K, ]
    HWij <- g5[(K + 1):(K * 2), ]
    HBWij <- HBij / HWij
    
    se.HBij <- se.Hij[1:K]
    se.HWij <- se.Hij[(K + 1):(K * 2)]
    se.HBWij <- se.Hij[-c(1:(K * 2))]
    
    
    # Hi
    
    g4 <- phi(A6i, g3i, "log")
    G4 <- dphi(A6i, g3i, G3i, "log")
    
    g5 <- phi(A7i, g4, "exp")
    G5 <- dphi(A7i, g4, G4, "exp")
    
    g6 <- phi(A8i, g5, "log")
    G6 <- dphi(A8i, g5, G5, "log")
    
    g7 <- phi(A9i, g6, "exp")
    G7 <- dphi(A9i, g6, G6, "exp")
    
    se.Hi <- sqrt(diag(G7 %*% (covtot %*% t(G7))))
    
    HBi <- g5[1:J, ]
    HWi <- g5[(J + 1):(J * 2), ]
    HBWi <- HBi / HWi
    
    se.HBi <- se.Hi[1:J]
    se.HWi <- se.Hi[(J + 1):(J * 2)]
    se.HBWi <- se.Hi[-c(1:(J * 2))]
    
    # H
    
    g4 <- phi(A6ii, g3ii, "log")
    G4 <- dphi(A6ii, g3ii, G3ii, "log")
    
    g5 <- phi(A7ii, g4, "exp")
    G5 <- dphi(A7ii, g4, G4, "exp")
    
    g6 <- phi(A8ii, g5, "log")
    G6 <- dphi(A8ii, g5, G5, "log")
    
    g7 <- phi(A9ii, g6, "exp")
    G7 <- dphi(A9ii, g6, G6, "exp")
    
    se.H <- sqrt(diag(G7 %*% (covtot %*% t(G7))))
    
    HB <- g5[1, ]
    HW <- g5[2, ]
    HBW <- HB / HW
    
    se.HB <- se.H[1]
    se.HW <- se.H[2]
    se.HBW <- se.H[3]
    
    
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
    
    Subs <- X[, 1]
    # Creating g3 and G3
    ns <- list()
    nuni <- matrix(0, J, g)
    for(i in 1:J) {
      Xa <- X[, i + 1]
      ns[[i]] <- matrix(0, g, LS)
      for(a in 0:m){
        nst <- tapply((Xa == a), Subs, sum)
        ns[[i]][a + 1, ] <- nst
        nuni[i, a + 1] <- sum(nst / (Rs * LS))
      }
    }
    
    Rss <- rep(Rs, Rs)
    Fwt <- Fbt <- Fet <- matrix(0, 1, g^2)
    Fw <- Fb <- Fe <- matrix(0, 1, K)
    for(k in 1:K){
      z <- cols[, k]
      Xa <- X[, z[1] + 1]
      Xb <- X[, z[2] + 1]
      Weights <- MLweight(X[, c(1, z + 1)], minx = 0, maxx = m)
      
      for(x in 1:g^2){
        if(Weights[x] > 0){
          i <- Patterns[x, 1]
          j <- Patterns[x, 2]
          
          nw <- sum((Xa == i & Xb == j) / (Rss * LS))
          at <- sum((Xa == i) * (rep(ns[[z[2]]][j + 1, ], Rs) - (Xb == j)) / (Rss * (Rss - 1) * LS))
          
          Fwt[x] <- Weights[x] * nw
          Fbt[x] <- Weights[x] * at
          Fet[x] <- Weights[x] * nuni[z[[1]], i + 1] * nuni[z[[2]], j + 1] 
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
      Hij <- data.frame(t(rbind(HWij, HBij, HBWij)), row.names = nams) 
      Hi <- data.frame(HWi, HBi = HBi, HBWi, row.names = labels) 
      H <- data.frame(HW, HB, HBW)
      OL <- list(Hij = Hij, Hi = Hi, H = H)
    } 
  }
  
  
  return(OL)
  
}

# example
# X <-  data.frame(Subs = c(1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3),
#       Xa   = c(0, 0, 1, 0, 1, 1, 1, 2, 1, 0, 1, 2, 0, 0, 0), 
#       Xb   = c(0, 0, 1, 0, 2, 2, 2, 1, 2, 1, 2, 2, 1, 1, 0))
