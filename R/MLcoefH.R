######## Letty Koopman 
######## University of Amsterdam
# MLcoefH() version 27-11-2020
#
## Update 27-11-2020: 
# Argument type.ci is added to distinguish between range-preserving and wald-based confidence intervals (Koopman et al 2020 a two-step, test-guided MSA for nonclustered and clustered data (Quality of Life Research))
#
# MLcoefH() version 11-09-2020
#
## Updates: 
#  Argument ci is added to compute range-preserving confidence intervals (Koopman et al 2020 range-preserving confidence intervals (IMPS proceedings))
#  The variance-covariance matrices are invisibly printed but can be extracted
## Previous updates:
#  Weighted proportions are used rather than average proportions (see Koopman et al 2020 a two-step, test-guided MSA for nonclustered and clustered data (Quality of Life Research))
#  but averaged proportions can be used if weigh.props = FALSE
#  Uses weights() rather than MLweight() (such that weighted proportions are used)
#  check.data has been changed to check.ml.data (previously the first minx subjects were ignored with minx the minimum score value).
#  A fixed item ordering now can be used as argument
#  More efficient by using vectors with pasted score patterns, rather than matrices
#  Uses the harmonic mean when samples sizes are unequal


"MLcoefH" <- function(X, se = TRUE, ci = FALSE, nice.output = TRUE, 
                      subject = 1, fixed.itemstep.order = NULL, 
                      weigh.props = TRUE, type.ci = "WB", cov.mat = FALSE){
  # Computes the two-level scalability coefficients in Mokken scale analysis
  #
  # Args:
  #   X: Data matrix with a subject column and one column per item. Preferably the subject column consists of integers.
  #   se: If TRUE, computes the standard errors for the coefficients, 
  #       if FALSE, only the coefficients are computed. Default is TRUE.
  #   ci: The confidence level of the 95% Wald-based confidence intervals. 
  #       if FALSE, the intervals are not printed.
  #   nice.output: If TRUE, prints the coefficients and standard errors in a matrix with nice lay-out,
  #                if FALSE, they are printed in a regular type matrix which can be used for further computations. Default is TRUE.
  #   Subject: Represents the subject column. Default is column 1. 
  #   fixed.itemstep.order: 
  #   weight.props: If TRUE: Use weighted proportions across groups to estimate coefficients and standard errors,
  #                 if FALSE: Use averaged proportions across groups to estimate coefficients and standard errors. Default is TRUE
  #   type.ci: The type of confidence interval, can be "WB" for Wald-based or "RP" for range-preserving.
  # 
  # Supporting mokken functions "MLweight", "weights", "check.ml.data", "all.patterns", "phi", and "dphi", "print.scalCoefsSE"
  #
  # Returns: 
  #   Two-level scalability coefficients and optionally their standard errors and confidence intervals.
  
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
  X <- check.ml.data(X) 
  if(is.null(colnames(X))) colnames(X) <- c("Subs", paste("Item", 1:(ncol(X) - 1)))
  
  # Ensure item names are unique
  if(any(table(colnames(X)) > 1)){
    warning("At least two items have the same name. The duplicates have received an index.")
    colnames(X) <- make.unique(colnames(X))
  }
  
  # Ensure variables are not constants
  if (any(apply(X, 2, var) < eps)) 
    stop("One or more variables have zero variance.")
  
  # Ensure given item step order has correct format
  if (!is.null(fixed.itemstep.order) && !is.matrix(fixed.itemstep.order)) {
    fixed.itemstep.order <- NULL
    warning("fixed.itemstep.order is not a matrix: fixed.itemstep.order ignored")
  }
  if (!is.null(fixed.itemstep.order) && is.matrix(fixed.itemstep.order)) 
    if (ncol(fixed.itemstep.order) != ncol(X[, -1]) && nrow(fixed.itemstep.order) != 
        max(X[, -1]) && sort(as.numeric(fixed.itemstep.order)) != 
        1:(max(X[, -1]) * ncol(X[, -1]))) {
      fixed.itemstep.order <- NULL
      warning("fixed.itemstep.order as incorrect dimensions and/or incorrect values: fixed.itemstep.order ignored")
    }
  if (is.null(ci)) {
    warning("ci needs to be the confidence level between zero and one, the default of .95 if used. Use ci = FALSE if the ci should not be printed.")
    ci <- .95
  } else {
    if(type.ci != "WB" & type.ci != "RP" & ci != FALSE) {
      warning("type.ci needs to be 'WB' (Wald-based) or 'RP' (range-preserving), the default 'WB' is used.")
      type.ci <- "WB"
    }
    if (ci < 0L | ci >= 1L | is.na(ci) | (!is.numeric(ci) & ci != FALSE)) {
      warning("ci needs to be the confidence level between zero and one, the default of .95 if used. Use ci = FALSE if the ci should not be printed.")
      ci <- .95
    } 
  }
  # Ensure each subject has > 1 rater
  if(any(Rs == 1)){ 
    warning('For at least one subject there is only 1 rater. The scalability coefficients are computed without this (these) subject(s).') 
    X <- X[!(X[, 1] %in% which(Rs == 1)), ]
    Rs <- as.numeric(table(X[, 1]))
    LS <- length(Rs)
    S <- 1:LS
    X[, 1] <- rep(S, Rs)
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
  
  
  if(se | ci){
    Xred <- apply(X, 1, paste, collapse=",")
    uniqueRows <- which(!duplicated(Xred))
    R <- X[uniqueRows, ]
    Rred <- Xred[uniqueRows]
    Vred <- apply(X[, -1], 1, paste, collapse = ",")
    Tred <- Vred[uniqueRows]
    
    SubsX <- X[, 1]
    SubsR <- SubsX[uniqueRows]
    
    Rss <- rep(Rs, Rs)
    RRs <- table(SubsR)
    Rd <- rep(Rs, RRs)
    
    n <- as.numeric(table(factor(Xred, levels=Rred))) 
    
    npred <- tapply(n, Tred, sum)
    Rtred <- names(npred)
    Lst <- length(npred)
    if(weigh.props == TRUE) {
      nprelred <- tapply(n / (sum(Rs)), Tred, sum)#/ (Rd * LS), Tred, sum) 
    } else {
      nprelred <- tapply(n / (Rd * LS), Tred, sum)
    }
    
    nNred <- tapply(1:length(Rred), Tred, unique)
    
    covps <- matrix(0, Lst, Lst)
    p <- matrix(0, Lst)
    
    if(weigh.props == TRUE) {
      for(s in S) {
        pt <- rowSums(outer(Rtred, Vred[SubsX == s], "==")) 
        prows <- which(pt > 0)
        pt <- pt[prows]
        p[prows] <- p[prows] + pt
        covps[prows, prows] <- covps[prows, prows] + (pt %*% t(pt)) / Rs[s] 
      }
      p <- as.numeric(p / sum(Rs)) # E(p) Eq 52 Koopman et al 2019 Standard Errors
      covps <- covps / sum(Rs) - (p %*% t(p)) # E(p p') 
      covp <- (diag(length(p)) * p - (p %*% t(p))) 
      nu <- LS / sum(1 / Rs)
      
      # Variance covariance matrix needed for delta method 
      covtot <- LS * nu * (nu - 1) * covps + LS * nu * covp 
      
      
      # Creating g3 and G3
      ns <- list()
      nuni <- matrix(0, J, g) # univariate props
      for(i in 1:J) {
        Xa <- R[, i + 1]
        ns[[i]] <- matrix(0, g, LS)
        for(a in 0:m){
          nst <- tapply((Xa == a) * n, SubsR, sum)
          ns[[i]][a + 1, ] <- nst
          nuni[i, a + 1] <- sum(nst / sum(Rs))
        }
      }
      
    } else {
      for(s in S) {
        pt <- rowSums(outer(Rtred, Vred[SubsX == s], "==")) / Rs[s]
        prows <- which(pt > 0)
        pt <- pt[prows]
        p[prows] <- p[prows] + pt
        covps[prows, prows] <- covps[prows, prows] + pt %*% t(pt)
      }
      p <- as.numeric(p / LS)
      covps <- covps / LS - (p %*% t(p))
      covp <- (diag(length(p)) * p - (p %*% t(p))) 
      nu <- LS / sum(1 / Rs)
      
      # Variance covariance matrix needed for delta method 
      covtot <- LS * nu * (nu - 1) * covps + LS * nu * covp
      
      
      # Creating g3 and G3
      ns <- list()
      nuni <- matrix(0, J, g) # univariate props
      for(i in 1:J) {
        Xa <- R[, i + 1]
        ns[[i]] <- matrix(0, g, LS)
        for(a in 0:m){
          nst <- tapply((Xa == a) * n, SubsR, sum)
          ns[[i]][a + 1, ] <- nst
          nuni[i, a + 1] <- sum(nst / (Rs * LS))
        }
      }
      
    }
    
    
    
    
    Rt <- matrix(as.numeric(unlist(strsplit(Rtred, "[,]"))), nrow = Lst, byrow = T)
    
    Fwt <- Fbt <- Fet <- Fw <- Fb <- Fe <- eij <- NULL
    G3W <- G3B <- G3E <- matrix(0, K, length(npred))
    for(k in 1:K){
      z <- cols[, k]
      Ra <- R[, z[1] + 1]
      Rb <- R[, z[2] + 1]
      Ta <- Rt[, z[1]]
      Tb <- Rt[, z[2]]
      if (is.null(fixed.itemstep.order)) {
        if(weigh.props == TRUE) {
          Weights <- weights(X[, z + 1], minx = 0, maxx = m)
        } else {
          Weights <- MLweight(X[, c(1, z + 1)], minx = 0, maxx = m)
        }
      } else {
        if(weigh.props == TRUE) {
          Weights <- weights(X[, z + 1], minx = 0, maxx = m, itemstep.order = fixed.itemstep.order[, z])
        } else {
          Weights <- MLweight(X[, c(1, z + 1)], minx = 0, maxx = m, itemstep.order = fixed.itemstep.order[, z])
        }
      }
      
      Wmat <- Bmat <- Emat <- matrix(0, g^2, length(npred))
      
      for(x in 1:g^2){
        if(Weights[x] > 0){
          i <- Patterns[x, 1]
          j <- Patterns[x, 2]
          
          if(weigh.props == TRUE) {
            temp <- (Ra == i) * (rep(ns[[z[2]]][j + 1, ], RRs) - (Rb == j)) * n
            nw <- sum(n * (Ra == i & Rb == j) / sum(Rs))
            at <- temp / (sum(Rs * (Rs - 1)))
            at <- sapply(1:length(npred), function(x) sum(at[nNred[[x]]]))
            
            Wmat[x, ] <- ((Ta == i & Tb == j) * nprelred / npred) * Weights[x]
            Bmat[x, ] <- (at / npred) * Weights[x] 
            Emat[x, ] <- Weights[x] / npred * ((Tb == j) * nuni[z[[1]], i + 1] + (Ta == i) * nuni[z[[2]], j + 1]) * nprelred
            Fwt[x] <- Weights[x] * nw
            Fbt[x] <- Weights[x] * sum(tapply(temp, SubsR, sum) / (sum(Rs * (Rs - 1))))
            Fet[x] <- Weights[x] * nuni[z[[1]], i + 1] * nuni[z[[2]], j + 1] 
          } else {
            temp <- (Ra == i) * (rep(ns[[z[2]]][j + 1, ], RRs) - (Rb == j)) * n
            nw <- sum(n * (Ra == i & Rb == j) / (rep(Rs, RRs) * LS)) 
            at <- temp / (rep(Rs, RRs) * (rep(Rs, RRs) - 1) * LS) 
            at <- sapply(1:length(npred), function(x) sum(at[nNred[[x]]]))
            
            Wmat[x, ] <- ((Ta == i & Tb == j) * nprelred / npred) * Weights[x]
            Bmat[x, ] <- (at / npred) * Weights[x] 
            Emat[x, ] <- Weights[x] / npred * ((Tb == j) * nuni[z[[1]], i + 1] + (Ta == i) * nuni[z[[2]], j + 1]) * nprelred
            Fwt[x] <- Weights[x] * nw
            Fbt[x] <- Weights[x] * sum(tapply(temp, SubsR, sum) / (Rs * (Rs - 1) * LS)) 
            Fet[x] <- Weights[x] * nuni[z[[1]], i + 1] * nuni[z[[2]], j + 1] 
          }
          
          
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
    G3Wi <- G3Bi <- G3Ei <- matrix(0, J, length(npred))
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
    
    ACM.Hij <- G7 %*% (covtot %*% t(G7))
    se.Hij <- sqrt(diag(ACM.Hij))
    
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
    
    ACM.Hi <- G7 %*% (covtot %*% t(G7))
    se.Hi <- sqrt(diag(ACM.Hi))
    
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
    
    ACM.H <- G7 %*% (covtot %*% t(G7))
    se.H <- sqrt(diag(ACM.H))
    
    HB <- g5[1, ]
    HW <- g5[2, ]
    HBW <- HB / HW
    
    se.HB <- se.H[1]
    se.HW <- se.H[2]
    se.HBW <- se.H[3]
    
    
    if(nice.output){
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
      OM.Hij <- H.ci.matrix.Hij <- matrix(NA, J + 3, J * 2 + 1)
      ci.matrix.Hij <- matrix(NA, J + 3, J + 1)
      if(type.ci == "WB") { # Wald-based
        for (j in 2 * (1:J)) {
          H.ci.matrix.Hij[, j ] <- OM.Hij[, j ] <- c("", "", "", format(paste(" ", formatC(round(Hij[, j/2], 3), digits = 3, format = "f"), " ", sep = ""), width = 7, justify = "right"))
          OM.Hij[, j + 1] <- c("", "", "", format(paste("(", formatC(round(se.Hij[, j/2], 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right"))
          ci.matrix.Hij[, j / 2 + 1] <- 
            H.ci.matrix.Hij[, j + 1] <- 
            c("", "", "", format(paste("[", 
                                       formatC(round(Hij[, j/2] + qnorm((1 - ci) / 2) * se.Hij[, j/2], 3), digits = 3, format = "f"), 
                                       ", ",
                                       formatC(round(Hij[, j/2] - qnorm((1 - ci) / 2) * se.Hij[, j/2], 3), digits = 3, format = "f"), 
                                       "] ", 
                                       sep = ""), width = 7, justify = "right"))
        }
      } else { # Range-preserving
        gHij <- log(1 - Hij)
        VgHij <- se.Hij^2 / (1 - Hij)^2
        
        for (j in 2 * (1:J)) {
          H.ci.matrix.Hij[, j ] <- OM.Hij[, j ] <- c("", "", "", format(paste(" ", formatC(round(Hij[, j/2], 3), digits = 3, format = "f"), " ", sep = ""), width = 7, justify = "right"))
          OM.Hij[, j + 1] <- c("", "", "", format(paste("(", formatC(round(se.Hij[, j/2], 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right"))
          ci.matrix.Hij[, j / 2 + 1] <- 
            H.ci.matrix.Hij[, j + 1] <- 
            c("", "", "", format(paste("[", 
                                       formatC(round(1 - exp(gHij[, j/2] - qnorm((1 - ci) / 2) * sqrt(VgHij[, j/2])), 3), digits = 3, format = "f"), 
                                       ", ",
                                       formatC(round(1 - exp(gHij[, j/2] + qnorm((1 - ci) / 2) * sqrt(VgHij[, j/2])), 3), digits = 3, format = "f"), 
                                       "] ", 
                                       sep = ""), width = 7, justify = "right"))
        }
      }
      OM.Hij[, 1] <- 
        OM.Hij[-c(1:3), -1][row(OM.Hij[-c(1:3), -1]) == (0.5 * col(OM.Hij[-c(1:3), -1]) + 0.5)] <- 
        OM.Hij[-c(1:3), -1][row(OM.Hij[-c(1:3), -1]) == (0.5 * col(OM.Hij[-c(1:3), -1]))] <- 
        H.ci.matrix.Hij[, 1] <- 
        H.ci.matrix.Hij[-c(1:3), -1][row(OM.Hij[-c(1:3), -1]) == (0.5 * col(OM.Hij[-c(1:3), -1]) + 0.5)] <- 
        H.ci.matrix.Hij[-c(1:3), -1][row(OM.Hij[-c(1:3), -1]) == (0.5 * col(OM.Hij[-c(1:3), -1]))] <- 
        ci.matrix.Hij[, 1] <- 
        ci.matrix.Hij[-c(1:3), -1][row(ci.matrix.Hij[-c(1:3), -1]) == (col(ci.matrix.Hij[-c(1:3), -1]))] <- 
        format("", width = 7, justify = "right")
      OM.Hij[round(J / 2) + 3, 1] <- H.ci.matrix.Hij[round(J / 2) + 3, 1] <- ci.matrix.Hij[round(J / 2) + 3, 1] <- format("(HWij)", width = 7, justify = "centre")
      OM.Hij[2, round(J / 2) * 2] <- H.ci.matrix.Hij[2, round(J / 2) * 2] <- ci.matrix.Hij[2, round(J / 2)] <- format("(HBij)", width = 7, justify = "centre")
      rownames(OM.Hij) <- rownames(H.ci.matrix.Hij) <- rownames(ci.matrix.Hij) <- c("", "", "", labels)
      colnames(OM.Hij) <- c("", new.labels)
      new.labels[2 * 1:J] <- paste(ci * 100, "% ci", sep = "")
      colnames(H.ci.matrix.Hij) <- c("", new.labels)
      colnames(ci.matrix.Hij) <- c("", labels)
      OM.Hij <- noquote(OM.Hij)
      H.ci.matrix.Hij <- noquote(H.ci.matrix.Hij)
      ci.matrix.Hij <- noquote(ci.matrix.Hij)
      
      # HWi & HBi
      
      
      OM.Hi <- matrix(NA, J, 10)
      if(type.ci == "WB") { # Wald-based
        OM.Hi[, 1] <- format("", width = 7, justify = "right")
        OM.Hi[, 2] <- format(formatC(round(HWi, 3), digits = 3, format = "f"), width = 7, justify = "right")
        OM.Hi[, 3] <- format(paste("(", formatC(round(se.HWi, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
        OM.Hi[, 4] <- format(paste("[",
                                   formatC(round(HWi + qnorm((1 - ci) / 2) * se.HWi, 
                                                 3), digits = 3, format = "f"),
                                   ", ",
                                   formatC(round(HWi - qnorm((1 - ci) / 2) * se.HWi, 
                                                 3), digits = 3, format = "f"),
                                   "] ", sep = ""), width = 7, 
                             justify = "right")
        OM.Hi[, 5] <- format(formatC(round(HBi, 3), digits = 3, format = "f"), width = 7, justify = "right")
        OM.Hi[, 6] <- format(paste("(", formatC(round(se.HBi, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
        OM.Hi[, 7] <- format(paste("[",
                                   formatC(round(HBi + qnorm((1 - ci) / 2) * se.HBi, 
                                                 3), digits = 3, format = "f"),
                                   ", ",
                                   formatC(round(HBi - qnorm((1 - ci) / 2) * se.HBi, 
                                                 3), digits = 3, format = "f"),
                                   "] ", sep = ""), width = 7, 
                             justify = "right")
        OM.Hi[, 8] <- format(formatC(round(HBWi, 3), digits = 3, format = "f"), width = 7, justify = "right")
        OM.Hi[, 9] <- format(paste("(", formatC(round(se.HBWi, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
        OM.Hi[, 10] <- format(paste("[",
                                    formatC(round(HBWi + qnorm((1 - ci) / 2) * se.HBWi, 3), digits = 3, format = "f"),
                                    ", ",
                                    formatC(round(HBWi - qnorm((1 - ci) / 2) * se.HBWi, 3), digits = 3, format = "f"),
                                    "] ", sep = ""), width = 7, 
                              justify = "right")
      } else { # Range-preserving
        gHBi <- log(1 - HBi)
        gHWi <- log(1 - HWi)
        VgHBi <- se.HBi^2 / (1 - HBi)^2
        VgHWi <- se.HWi^2 / (1 - HWi)^2
        
        OM.Hi[, 1] <- format("", width = 7, justify = "right")
        OM.Hi[, 2] <- format(formatC(round(HWi, 3), digits = 3, format = "f"), width = 7, justify = "right")
        OM.Hi[, 3] <- format(paste("(", formatC(round(se.HWi, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
        OM.Hi[, 4] <- format(paste("[",
                                   formatC(round(1 - exp(gHWi - qnorm((1 - ci) / 2) * sqrt(VgHWi)), 
                                                 3), digits = 3, format = "f"),
                                   ", ",
                                   formatC(round(1 - exp(gHWi + qnorm((1 - ci) / 2) * sqrt(VgHWi)), 
                                                 3), digits = 3, format = "f"),
                                   "] ", sep = ""), width = 7, 
                             justify = "right")
        OM.Hi[, 5] <- format(formatC(round(HBi, 3), digits = 3, format = "f"), width = 7, justify = "right")
        OM.Hi[, 6] <- format(paste("(", formatC(round(se.HBi, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
        OM.Hi[, 7] <- format(paste("[",
                                   formatC(round(1 - exp(gHBi - qnorm((1 - ci) / 2) * sqrt(VgHBi)), 
                                                 3), digits = 3, format = "f"),
                                   ", ",
                                   formatC(round(1 - exp(gHBi + qnorm((1 - ci) / 2) * sqrt(VgHBi)), 
                                                 3), digits = 3, format = "f"),
                                   "] ", sep = ""), width = 7, 
                             justify = "right")
        OM.Hi[, 8] <- format(formatC(round(HBWi, 3), digits = 3, format = "f"), width = 7, justify = "right")
        OM.Hi[, 9] <- format(paste("(", formatC(round(se.HBWi, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
        OM.Hi[, 10] <- format(paste("[",
                                    formatC(round(HBWi + qnorm((1 - ci) / 2) * se.HBWi, 3), digits = 3, format = "f"),
                                    ", ",
                                    formatC(round(HBWi - qnorm((1 - ci) / 2) * se.HBWi, 3), digits = 3, format = "f"),
                                    "] ", sep = ""), width = 7, 
                              justify = "right")
      }
      dimnames(OM.Hi) <- list(labels, c("", 
                                        "   HWi", " (se)  ", paste("    ", ci * 100, "% ci", sep = ""), 
                                        "   HBi", " (se)  ", paste("    ", ci * 100, "% ci", sep = ""), 
                                        "   BWi", " (se)  ", paste("    ", ci * 100, "% ci", sep = "")))
      OM.Hi <- noquote(OM.Hi)
      
      
      
      
      
      # HW & HB
      
      OM.H <- matrix(NA, 1, 10)
      if(type.ci == "WB") { # Wald-based
        OM.H[, 1] <- format("", width = 7, justify = "right")
        OM.H[, 2] <- format(formatC(round(HW, 3), digits = 3, format = "f"), width = 7, justify = "right")
        OM.H[, 3] <- format(paste("(", formatC(round(se.HW, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
        OM.H[, 4] <- format(paste("[",
                                  formatC(round(HW + qnorm((1 - ci) / 2) * se.HW, 
                                                3), digits = 3, format = "f"),
                                  ", ",
                                  formatC(round(HW - qnorm((1 - ci) / 2) * se.HW, 
                                                3), digits = 3, format = "f"),
                                  "] ", sep = ""), width = 7, 
                            justify = "right")
        OM.H[, 5] <- format(formatC(round(HB, 3), digits = 3, format = "f"), width = 7, justify = "right")
        OM.H[, 6] <- format(paste("(", formatC(round(se.HB, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
        OM.H[, 7] <- format(paste("[",
                                  formatC(round(HB + qnorm((1 - ci) / 2) * se.HB, 
                                                3), digits = 3, format = "f"),
                                  ", ",
                                  formatC(round(HB - qnorm((1 - ci) / 2) * se.HB, 
                                                3), digits = 3, format = "f"),
                                  "] ", sep = ""), width = 7, 
                            justify = "right")
        OM.H[, 8] <- format(formatC(round(HBW, 3), digits = 3, format = "f"), width = 7, justify = "right")
        OM.H[, 9] <- format(paste("(", formatC(round(se.HBW, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
        OM.H[, 10] <- format(paste("[",
                                   formatC(round(HBW + qnorm((1 - ci) / 2) * se.HBW, 3), digits = 3, format = "f"),
                                   ", ",
                                   formatC(round(HBW - qnorm((1 - ci) / 2) * se.HBW, 3), digits = 3, format = "f"),
                                   "] ", sep = ""), width = 7, 
                             justify = "right")
      } else { # Range-preserving
        gHB <- log(1 - HB)
        gHW <- log(1 - HW)
        VgHB <- se.HB^2 / (1 - HB)^2
        VgHW <- se.HW^2 / (1 - HW)^2
        
        OM.H[, 1] <- format("", width = 7, justify = "right")
        OM.H[, 2] <- format(formatC(round(HW, 3), digits = 3, format = "f"), width = 7, justify = "right")
        OM.H[, 3] <- format(paste("(", formatC(round(se.HW, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
        OM.H[, 4] <- format(paste("[",
                                  formatC(round(1 - exp(gHW - qnorm((1 - ci) / 2) * sqrt(VgHW)), 
                                                3), digits = 3, format = "f"),
                                  ", ",
                                  formatC(round(1 - exp(gHW + qnorm((1 - ci) / 2) * sqrt(VgHW)), 
                                                3), digits = 3, format = "f"),
                                  "] ", sep = ""), width = 7, 
                            justify = "right")
        OM.H[, 5] <- format(formatC(round(HB, 3), digits = 3, format = "f"), width = 7, justify = "right")
        OM.H[, 6] <- format(paste("(", formatC(round(se.HB, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
        OM.H[, 7] <- format(paste("[",
                                  formatC(round(1 - exp(gHB - qnorm((1 - ci) / 2) * sqrt(VgHB)), 
                                                3), digits = 3, format = "f"),
                                  ", ",
                                  formatC(round(1 - exp(gHB + qnorm((1 - ci) / 2) * sqrt(VgHB)), 
                                                3), digits = 3, format = "f"),
                                  "] ", sep = ""), width = 7, 
                            justify = "right")
        OM.H[, 8] <- format(formatC(round(HBW, 3), digits = 3, format = "f"), width = 7, justify = "right")
        OM.H[, 9] <- format(paste("(", formatC(round(se.HBW, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
        OM.H[, 10] <- format(paste("[",
                                   formatC(round(HBW + qnorm((1 - ci) / 2) * se.HBW, 3), digits = 3, format = "f"),
                                   ", ",
                                   formatC(round(HBW - qnorm((1 - ci) / 2) * se.HBW, 3), digits = 3, format = "f"),
                                   "] ", sep = ""), width = 7, 
                             justify = "right")
      }
      dimnames(OM.H) <- list("Scale", c("", 
                                        "   HW", " (se)  ", paste("    ", ci * 100, "% ci", sep = ""), 
                                        "   HB", " (se)  ", paste("    ", ci * 100, "% ci", sep = ""), 
                                        "   BW", " (se)  ", paste("    ", ci * 100, "% ci", sep = "")))
      OM.H <- noquote(OM.H)
      
      # Output:
      if(se & ci) {
        OLcoefs <- list(Hij = OM.Hij, ci.Hij = ci.matrix.Hij, Hi = OM.Hi, H = OM.H)
        OLcoefs <- setNames(OLcoefs, 
                            c("Hij", paste(ci*100, "% ci", sep = ""), "Hi", "H"))
        
      } else if (!se) {
        OLcoefs <- list(Hij = H.ci.matrix.Hij, Hi = OM.Hi[, -c(3, 6, 9)], H = OM.H[, -c(3, 6, 9)])
      } else {
        OLcoefs <- list(Hij = OM.Hij, Hi = OM.Hi[, -c(4, 7, 10)], H = OM.H[, -c(4, 7, 10)])
      }
      OL <- OLcoefs
      
    } else { # nice output == F
      if(ci) {
        if(type.ci == "WB") { # Wald-based
          ci.HWij <- data.frame("ci.HWij.l" = HWij + qnorm((1 - ci) / 2) * se.HWij, 
                                "ci.HWij.u" = HWij - qnorm((1 - ci) / 2) * se.HWij)
          ci.HBij <- data.frame("ci.HBij.l" = HBij + qnorm((1 - ci) / 2) * se.HBij, 
                                "ci.HBij.u" = HBij - qnorm((1 - ci) / 2) * se.HBij)
          ci.HBWij <- data.frame("ci.HBWij.l" = HBWij + qnorm((1 - ci) / 2) * se.HBWij, 
                                 "ci.HBWij.u" = HBWij - qnorm((1 - ci) / 2) * se.HBWij)
          
          ci.HWi <- data.frame("ci.HWi.l" = HWi + qnorm((1 - ci) / 2) * se.HWi, 
                               "ci.HWi.u" = HWi - qnorm((1 - ci) / 2) * se.HWi)
          ci.HBi <- data.frame("ci.HBi.l" = HBi + qnorm((1 - ci) / 2) * se.HBi, 
                               "ci.HBi.u" = HBi - qnorm((1 - ci) / 2) * se.HBi)
          ci.HBWi <- data.frame("ci.HBWi.l" = HBWi + qnorm((1 - ci) / 2) * se.HBWi, 
                                "ci.HBWi.u" = HBWi - qnorm((1 - ci) / 2) * se.HBWi)
          
         ci.HW <- data.frame("ci.HW.l" = HW + qnorm((1 - ci) / 2) * se.HW, 
                              "ci.HW.u" = HW - qnorm((1 - ci) / 2) * se.HW)
          ci.HB <- data.frame("ci.HB.l" = HB + qnorm((1 - ci) / 2) * se.HB, 
                              "ci.HB.u" = HB - qnorm((1 - ci) / 2) * se.HB)
          ci.HBW <- data.frame("ci.HBW.l" = HBW + qnorm((1 - ci) / 2) * se.HBW, 
                               "ci.HBW.u" = HBW - qnorm((1 - ci) / 2) * se.HBW)
          
        } else { # Range-preserving
          gHWij <- log(1 - HWij)
          gHBij <- log(1 - HBij)
          VgHWij <- se.HWij^2 / (1 - HWij)^2
          VgHBij <- se.HBij^2 / (1 - HBij)^2
          
          ci.HWij <- data.frame("ci.HWij.l" = 1 - exp(gHWij - qnorm((1 - ci) / 2) * sqrt(VgHWij)), 
                                "ci.HWij.u" = 1 - exp(gHWij + qnorm((1 - ci) / 2) * sqrt(VgHWij)))
          ci.HBij <- data.frame("ci.HBij.l" = 1 - exp(gHBij - qnorm((1 - ci) / 2) * sqrt(VgHBij)), 
                                "ci.HBij.u" = 1 - exp(gHBij + qnorm((1 - ci) / 2) * sqrt(VgHBij)))
          ci.HBWij <- data.frame("ci.HBWij.l" = HBWij + qnorm((1 - ci) / 2) * se.HBWij, 
                                 "ci.HBWij.u" = HBWij - qnorm((1 - ci) / 2) * se.HBWij)
          
          gHWi <- log(1 - HWi)
          gHBi <- log(1 - HBi)
          VgHWi <- se.HWi^2 / (1 - HWi)^2
          VgHBi <- se.HBi^2 / (1 - HBi)^2
          
          ci.HWi <- data.frame("ci.HWi.l" = 1 - exp(gHWi - qnorm((1 - ci) / 2) * sqrt(VgHWi)), 
                               "ci.HWi.u" = 1 - exp(gHWi + qnorm((1 - ci) / 2) * sqrt(VgHWi)))
          ci.HBi <- data.frame("ci.HBi.l" = 1 - exp(gHBi - qnorm((1 - ci) / 2) * sqrt(VgHBi)), 
                               "ci.HBi.u" = 1 - exp(gHBi + qnorm((1 - ci) / 2) * sqrt(VgHBi)))
          ci.HBWi <- data.frame("ci.HBWi.l" = HBWi + qnorm((1 - ci) / 2) * se.HBWi, 
                                "ci.HBWi.u" = HBWi - qnorm((1 - ci) / 2) * se.HBWi)
          
          gHW <- log(1 - HW)
          gHB <- log(1 - HB)
          VgHW <- se.HW^2 / (1 - HW)^2
          VgHB <- se.HB^2 / (1 - HB)^2
          
          ci.HW <- data.frame("ci.HW.l" = 1 - exp(gHW - qnorm((1 - ci) / 2) * sqrt(VgHW)), 
                              "ci.HW.u" = 1 - exp(gHW + qnorm((1 - ci) / 2) * sqrt(VgHW)))
          ci.HB <- data.frame("ci.HB.l" = 1 - exp(gHB - qnorm((1 - ci) / 2) * sqrt(VgHB)), 
                              "ci.HB.u" = 1 - exp(gHB + qnorm((1 - ci) / 2) * sqrt(VgHB)))
          ci.HBW <- data.frame("ci.HBW.l" = HBW + qnorm((1 - ci) / 2) * se.HBW, 
                               "ci.HBW.u" = HBW - qnorm((1 - ci) / 2) * se.HBW)
        }
      }
      if(se & ci) {
        Hij <- data.frame(HWij, se.HWij, ci.HWij, 
                          HBij, se.HBij, ci.HBij, 
                          HBWij, se.HBWij, ci.HBWij) 
        rownames(Hij) <- nams
        Hi <- data.frame(HWi, se.HWi, ci.HWi,
                         HBi, se.HBi, ci.HBi,
                         HBWi, se.HBWi, ci.HBWi) 
        rownames(Hi) <- labels
        H <- data.frame(HW, se.HW, ci.HW,
                        HB, se.HB, ci.HB,
                        HBW, se.HBW, ci.HBW)
      } else if(!se) {
        Hij <- data.frame(HWij, ci.HWij, 
                          HBij, ci.HBij, 
                          HBWij, ci.HBWij) 
        rownames(Hij) <- nams
        Hi <- data.frame(HWi, ci.HWi,
                         HBi, ci.HBi,
                         HBWi, ci.HBWi) 
        rownames(Hi) <- labels
        H <- data.frame(HW, ci.HW,
                        HB, ci.HB,
                        HBW, ci.HBW)
      } else {
        Hij <- data.frame(HWij, se.HWij,  
                          HBij, se.HBij,  
                          HBWij,  se.HBWij) 
        rownames(Hij) <- nams
        Hi <- data.frame(HWi, se.HWi,
                         HBi, se.HBi, 
                         HBWi, se.HBWi) 
        rownames(Hi) <- labels
        H <- data.frame(HW, se.HW,
                        HB, se.HB,
                        HBW, se.HBW)
      }
      OL <- list(Hij = Hij, 
                 Hi = Hi, 
                 H = H)
    }
    
  } else { # alleen coefs, geen se/ci
    
    Subs <- X[, 1]
    # Creating g3 and G3
    ns <- list()
    nuni <- matrix(0, J, g)
    if(weigh.props == TRUE) {
      for(i in 1:J) {
        Xa <- X[, i + 1]
        ns[[i]] <- matrix(0, g, LS)
        for(a in 0:m){
          nst <- tapply((Xa == a), Subs, sum)
          ns[[i]][a + 1, ] <- nst
          nuni[i, a + 1] <- sum(nst / sum(Rs))#(Rs * LS))
        }
      }
      Rss <- rep(Rs, Rs)
      Fwt <- Fbt <- Fet <- matrix(0, 1, g^2)
      Fw <- Fb <- Fe <- matrix(0, 1, K)
      for(k in 1:K){
        z <- cols[, k]
        Xa <- X[, z[1] + 1]
        Xb <- X[, z[2] + 1]
        if (is.null(fixed.itemstep.order)) {
          Weights <- weights(X[, z + 1], minx = 0, maxx = m)#MLweight(X[, c(1, z + 1)], minx = 0, maxx = m)
        } else {
          Weights <- weights(X[, z + 1], minx = 0, maxx = m, itemstep.order = fixed.itemstep.order[, z])#MLweight(X[, c(1, z + 1)], minx = 0, maxx = m, itemstep.order = fixed.itemstep.order[, z])
        }
        
        for(x in 1:g^2){
          if(Weights[x] > 0){
            i <- Patterns[x, 1]
            j <- Patterns[x, 2]
            
            nw <- sum((Xa == i & Xb == j) / sum(Rs))#(Rss * LS))
            at <- sum((Xa == i) * (rep(ns[[z[2]]][j + 1, ], Rs) - (Xb == j)) / (sum(Rs * (Rs - 1))))#(Rss * (Rss - 1) * LS))
            
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
    } else {
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
        if (is.null(fixed.itemstep.order)) {
          Weights <- MLweight(X[, c(1, z + 1)], minx = 0, maxx = m)
        } else {
          Weights <- MLweight(X[, c(1, z + 1)], minx = 0, maxx = m, itemstep.order = fixed.itemstep.order[, z])
        }
        
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
    
    if(nice.output){
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
      Hij <- cbind(matrix(HWij), matrix(HBij), matrix(HBWij) )
      rownames(Hij) <- nams
      colnames(Hij) <- c("HWij", "HBij", "HBWij")
      Hi <- data.frame(HWi, HBi = HBi, HBWi, row.names = labels) 
      H <- data.frame(HW, HB, HBW)
      OL <- list(Hij = Hij, Hi = Hi, H = H)
    } 
  }
  if(cov.mat == TRUE) {
    if(ci == FALSE & se == FALSE) {
      warning("cov.mat can only be printed if ci or se = TRUE")
    } else {
      colnames(ACM.Hij) <- rownames(ACM.Hij) <- paste(rep(c("HBij", "HWij", "HBWij"), each = K), nams)
      colnames(ACM.Hi) <- rownames(ACM.Hi) <- paste(rep(c("HBi", "HWi", "HBWi"), each = J), labels)
      colnames(ACM.H) <- rownames(ACM.H) <- c("HB", "HW", "HBW")
      OL <- c(OL, list(covHij = ACM.Hij, 
                       covHi = ACM.Hi, 
                       covH = ACM.H))
    }
  }
  return(OL)
  
}

