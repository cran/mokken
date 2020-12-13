# Aangepast door Letty Koopman op 11 juni 2020
# Aangepast door Andries van der Ark op 2 september 2020

## Update 27-11-2020 by Letty: 
# Argument type.ci is added to distinguish between range-preserving and wald-based confidence intervals 
# (Koopman et al 2020 a two-step, test-guided MSA for nonclustered and clustered data (Quality of Life Research))


coefH <- function (X, se = TRUE, ci = FALSE, nice.output = TRUE, level.two.var = NULL, 
                   group.var = NULL, fixed.itemstep.order = NULL, type.ci = "WB", results = TRUE) {

  X <- check.data(X)
  eps <- 1e-40
  labels <- dimnames(X)[[2]]
  return.list <- print.list <- list()
  
  if (!is.null(group.var) && nrow(as.matrix(group.var)) != nrow(X)) {
    group.var <- NULL
    warning("group.var not the same length/nrow as X: group.var ignored")
  }
  
  if (!is.null(fixed.itemstep.order) && !is.matrix(fixed.itemstep.order)) {
    fixed.itemstep.order <- NULL
    warning("fixed.itemstep.order is not a matrix: fixed.itemstep.order ignored")
  }
  
  computeCI <- ifelse(is.logical(ci) & ci == FALSE, FALSE, TRUE) 
  if (computeCI & !(is.numeric(ci) & !is.na(ci) & ci > 0 & ci < 1)){
    warning("ci requires values between 0 and 1, ci is set to default value ci = .95.")
    ci <- .95
  }   
  if(computeCI) {
    if(type.ci != "WB" & type.ci != "RP") {
      warning("type.ci needs to be 'WB' (Wald-based) or 'RP' (range-preserving), the default 'WB' is used.")
      type.ci <- "WB"
    }
  }
  
  if (!is.null(fixed.itemstep.order) && is.matrix(fixed.itemstep.order)) 
    if (ncol(fixed.itemstep.order) != ncol(X) && nrow(fixed.itemstep.order) != max(X) && sort(as.numeric(fixed.itemstep.order)) != 1:(max(X) * ncol(X))) {
      fixed.itemstep.order <- NULL
      warning("fixed.itemstep.order as incorrect dimensions and/or incorrect values: fixed.itemstep.order ignored")
    }
  
  if(!is.null(level.two.var)) {
    if (nrow(as.matrix(level.two.var)) != nrow(X)) {
      level.two.var <- NULL
      warning("level.two.var not the same length/nrow as X: level.two.var ignored and one-level standard errors are returned")
    } else if(any(is.na(level.two.var))) {
      level.two.var <- NULL
      warning("level.two.var contains missing value(s): level.two.var is ignored.")
    } else {
      X <- X[order(level.two.var), ]
      level.two.var <- sort(level.two.var)
      Rs <- as.numeric(table(level.two.var))
      LS <- length(Rs)
      S <- 1:LS 
      level.two.var <- rep(S, Rs)
      # Ensure each subject has > 1 rater
      if(any(Rs == 1)){ 
        warning('For at least one group there is only 1 respondent. The scalability coefficients are computed without this (these) group(s).') 
        cases <- !(level.two.var %in% which(Rs == 1))
        X <- X[cases, ]
        level.two.var <- level.two.var[cases]
        group.var <- group.var[cases]
        Rs <- as.numeric(table(level.two.var))
        LS <- length(Rs)
        S <- 1:LS
        level.two.var <- rep(S, Rs)
      }
    }
  }
  
 
  if (!se && is.null(group.var) && is.null(fixed.itemstep.order) && !computeCI) {
    S <- var(X)
    Smax <- var(apply(X, 2, sort))
    Hij <- S/Smax
    diag(S) <- 0
    diag(Smax) <- 0
    Hi <- apply(S, 1, sum)/apply(Smax, 1, sum)
    H <- sum(S)/sum(Smax)
    return.list <- list(Hij = Hij, Hi = Hi, H = H)
    print.list <- list(Hij = Hij, Hi = Hi, H = H)
   
  } else {

    g <- max(X) - min(X) + 1
    J <- ncol(X)
    P <- choose(J, 2)
    N <- nrow(X)
    if (any(apply(X, 2, var) < eps)) 
      stop("One or more variables have zero variance")
    n <- as.matrix(table(apply(X, 1, paste, collapse = "")))
    lab.n <- matrix(names(table(apply(X, 1, paste, collapse = ""))))
    lab.b <- apply(all.patterns(2, g), 2, paste, collapse = "")
    lab.u <- as.character(0:(g - 1))
    Bi <- substr(lab.b, 1, 1)
    Bj <- substr(lab.b, 2, 2)
    R <- t(apply(lab.n, 1, string2integer))
    r <- length(lab.n)
    U <- list()
    for (j in 1:J) U[[j]] <- tabulate(X[, j] + 1, g)
    W <- list()
    WA <- list()
    WY <- list()
    WE <- list()
    WF <- list()
    for (i in 1:J) {
      W[[i]] <- list()
      WA[[i]] <- list()
      WY[[i]] <- list()
      WE[[i]] <- list()
      WF[[i]] <- list()
      for (j in i:J) if (j > i) {
        W[[i]][[j]] <- weights(X[, c(i, j)], g - 1)
        if (!is.null(fixed.itemstep.order)) 
          W[[i]][[j]] <- weights(X[, c(i, j)], g - 1, itemstep.order = fixed.itemstep.order[, c(i, j)])
        A1a <- NULL
        for (a in 0:(g - 1)) for (b in 0:(g - 1)) A1a <- rbind(A1a, as.numeric(R[, i] == a & R[, j] == b))
        WA[[i]][[j]] <- W[[i]][[j]] %*% A1a
        Eij <- matrix(U[[i]][as.numeric(Bi) + 1], nrow = g^2, ncol = 1) * matrix(U[[j]][as.numeric(Bj) + 1], nrow = g^2, ncol = 1)/N
        Y22 <- cbind(outer(Bi, lab.u, "=="), outer(Bj, lab.u, "==")) * matrix(Eij, nrow = g^2, ncol = 2 * g)
        Ri <- substr(lab.n, i, i)
        Rj <- substr(lab.n, j, j)
        Z2 <- rbind(outer(lab.u, Ri, "==")[, , 1], outer(lab.u, Rj, "==")[, , 1]) * c(1/U[[i]], 1/U[[j]])
        Z2[is.nan(Z2)] <- 1/eps
        YZ2 <- Y22 %*% Z2 - Eij %*% matrix(1/N, 1, r)
        WY[[i]][[j]] <- W[[i]][[j]] %*% YZ2 
        Fij <- complete.observed.frequencies(X[, c(i, j)], 2, g)
        WF[[i]][[j]] <- W[[i]][[j]] %*% Fij
        WE[[i]][[j]] <- W[[i]][[j]] %*% Eij
      }
    }
    g3 <- matrix(c(unlist(WF[[1]][[2]]), unlist(WF), unlist(WE)), nrow = 2 * P + 1, byrow = TRUE)
    A4 <- rbind(matrix(c(1, -1, rep(0, (J * (J - 1)) - 1)), 1, (J * (J - 1)) + 1), cbind(matrix(0, P, 1), diag(P), -1 * diag(P)))
    A5 <- cbind(matrix(1, P, 1), -1 * diag(P))
    g4 <- phi(A4, g3, "log")
    g5 <- phi(A5, g4, "exp")
    G3 <- matrix(c(unlist(WA[[1]][[2]]), unlist(WA), unlist(WY)), nrow = 2 * P + 1, byrow = TRUE)
    G4 <- dphi(A4, g3, G3, "log")
    G5ij <- dphi(A5, g4, G4, "exp")
    Hij <- se.Hij <- matrix(0, J, J)
    Hij[lower.tri(Hij)] <- g5
    Hij <- Hij + t(Hij)
    if (P > 1) G3 <- rbind(apply(G3[2:(P + 1), ], 2, sum), apply(G3[2:(P + 1), ], 2, sum), apply(G3[(P + 2):(2 * P + 1), ], 2, sum))
    g3 <- matrix(c(sum(g3[2:(P + 1), ]), sum(g3[2:(P + 1), ]), sum(g3[(P + 2):(2 * P + 1), ])), ncol = 1)
    A4 <- cbind(matrix(1, 2, 1), -1 * diag(2))
    A5 <- matrix(c(1, -1), 1, 2)
    g4 <- phi(A4, g3, "log")
    g5 <- phi(A5, g4, "exp")
    G4 <- dphi(A4, g3, G3, "log")
    G5 <- dphi(A5, g4, G4, "exp")
    H <- g5
    G3 <- matrix(0, 2 * J + 1, r)
    g3 <- matrix(0, 2 * J + 1, 1)
    for (j in 1:J) for (k in 1:J) if (k > j) {
      g3[j + 1, ] <- g3[j + 1, ] + WF[[j]][[k]]
      g3[J + j + 1, ] <- g3[J + j + 1, ] + WE[[j]][[k]]
      G3[j + 1, ] <- G3[j + 1, ] + WA[[j]][[k]]
      G3[J + j + 1, ] <- G3[J + j + 1, ] + WY[[j]][[k]]
    } else {
      if (k < j) {
        g3[j + 1, ] <- g3[j + 1, ] + WF[[k]][[j]]
        g3[J + j + 1, ] <- g3[J + j + 1, ] + WE[[k]][[j]]
        G3[j + 1, ] <- G3[j + 1, ] + WA[[k]][[j]]
        G3[J + j + 1, ] <- G3[J + j + 1, ] + WY[[k]][[j]]
      }
    }
    g3[1, ] <- g3[2, ]
    G3[1, ] <- G3[2, ]
    A4 <- rbind(matrix(c(1, -1, rep(0, (J * 2) - 1)), 1, (J * 2) + 1), cbind(matrix(0, J, 1), diag(J), -1 * diag(J)))
    A5 <- cbind(matrix(1, J, 1), -1 * diag(J))
    g4 <- phi(A4, g3, "log")
    g5 <- phi(A5, g4, "exp")
    G4 <- dphi(A4, g3, G3, "log")
    G5i <- dphi(A5, g4, G4, "exp")
    Hi <- matrix(g5)
   
    if (se | computeCI) {

       if(!is.null(level.two.var)) {
          pats.n <- matrix(apply(X, 1, paste, collapse = ""))
          n.p <- length(n)
          covps <- matrix(0, n.p, n.p)
          p <- matrix(0, n.p)
          for(s in S) {
             pt <- rowSums(outer(lab.n, pats.n[level.two.var == s], "==")) 
             prows <- which(pt > 0)
             pt <- pt[prows]
             p[prows] <- p[prows] + pt
             covps[prows, prows] <- covps[prows, prows] + (pt %*% t(pt)) / Rs[s]
          }
          p <- as.numeric(p / sum(Rs)) # E(p) Eq 52 Koopman et al 2019 Standard Errors
          covps <- covps / sum(Rs) - (p %*% t(p)) # E(p p') - E(p) E(p)'
          covp <- (diag(length(p)) * p - (p %*% t(p))) 
          nu <- LS / sum(1 / Rs)
        
          # Variance covariance matrix needed for delta method 
          covtot <- LS * nu * covp + LS * nu * (nu - 1) * covps 
        
          ACM.Hij = G5ij %*% (covtot %*% t(G5ij))
          se.Hij[lower.tri(se.Hij)] <- sqrt(diag(ACM.Hij))
          se.Hij <- se.Hij + t(se.Hij)
          dimnames(se.Hij) <- dimnames(Hij) <- list(labels, labels)
          ACM.Hi = G5i %*% (covtot %*% t(G5i))
          se.Hi <- matrix(sqrt(diag(ACM.Hi)))
          dimnames(se.Hi)[[1]] <- dimnames(Hi)[[1]] <- labels
          ACM.H <- G5 %*% (covtot %*% t(G5))
          se.H <- sqrt(diag(ACM.H))

         
       } else {
          ACM.Hij = G5ij %*% (as.numeric(n) * t(G5ij))
          se.Hij[lower.tri(se.Hij)] <- sqrt(diag(ACM.Hij))
          se.Hij <- se.Hij + t(se.Hij)
          dimnames(se.Hij) <- dimnames(Hij) <- list(labels, labels)
          ACM.Hi = G5i %*% (as.numeric(n) * t(G5i))
          se.Hi <- matrix(sqrt(diag(ACM.Hi)))
          dimnames(se.Hi)[[1]] <- dimnames(Hi)[[1]] <- labels
          ACM.H <- G5 %*% (as.numeric(n) * t(G5))
          se.H <- sqrt(diag(ACM.H))
      }

      
      
      
    }
    
    
    if (nice.output) {
      if (se | computeCI) {
        output.matrix.Hij <- H.ci.matrix.Hij <- matrix(NA, J, J * 2)
        ci.matrix.Hij <- matrix(NA, J, J)
        for (j in 2 * (1:J)) {
          H.ci.matrix.Hij[, j - 1] <- output.matrix.Hij[, j - 1] <- format(paste(" ", formatC(round(Hij[, j/2], 3), digits = 3, format = "f"), " ", sep = ""), width = 7, justify = "right")
          output.matrix.Hij[, j] <- format(paste("(", formatC(round(se.Hij[, j/2], 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
          if(type.ci == "WB") {
            ci.matrix.Hij[, j / 2] <- H.ci.matrix.Hij[, j] <- format(paste("[", 
                                                                           formatC(round(Hij[, j/2] + qnorm((1 - ci) / 2) * se.Hij[, j/2], 3), digits = 3, format = "f"), 
                                                                           ", ", 
                                                                           formatC(round(Hij[, j/2] - qnorm((1 - ci) / 2) * se.Hij[, j/2], 3), digits = 3, format = "f"), 
                                                                           "] ", sep = ""), width = 7, justify = "right")
          } else {
            gHij <- log(1 - Hij)
            VgHij <- se.Hij^2 / (1 - Hij)^2
            
            gHi <- log(1 - Hi)
            VgHi <- se.Hi^2 / (1 - Hi)^2
            
            gH <- log(1 - H)
            VgH <- se.H^2 / (1 - H)^2
            
            ci.matrix.Hij[, j / 2] <- H.ci.matrix.Hij[, j] <- format(paste("[", 
                                                                           formatC(round(1 - exp(gHij[, j/2] - qnorm((1 - ci) / 2) * sqrt(VgHij[, j/2])), 3), digits = 3, format = "f"), 
                                                                           ", ", 
                                                                           formatC(round(1 - exp(gHij[, j/2] + qnorm((1 - ci) / 2) * sqrt(VgHij[, j/2])), 3), digits = 3, format = "f"), 
                                                                           "] ", sep = ""), width = 7, justify = "right")
          }
        }
        new.labels <- rep(labels, each = 2)
        new.labels[2 * (1:J)] <- "se"
        dimnames(output.matrix.Hij)[[1]] <- 
          dimnames(ci.matrix.Hij)[[1]] <- 
          dimnames(ci.matrix.Hij)[[2]] <- 
          dimnames(H.ci.matrix.Hij)[[1]] <- labels
        dimnames(output.matrix.Hij)[[2]] <- new.labels
        new.labels[2 * (1:J)] <- paste(ci * 100, "% ci", sep = "")
        dimnames(H.ci.matrix.Hij)[[2]] <- new.labels
        output.matrix.Hij[row(output.matrix.Hij) == 0.5 * col(output.matrix.Hij)] <- format("", width = 7, justify = "right")
        output.matrix.Hij[(row(output.matrix.Hij)) == (0.5 * col(output.matrix.Hij) + 0.5)] <- format("", width = 7, justify = "right")
        output.matrix.Hij <- noquote(output.matrix.Hij)
        
        H.ci.matrix.Hij[row(H.ci.matrix.Hij) == 0.5 * col(H.ci.matrix.Hij)] <- format("", width = 7, justify = "right")
        H.ci.matrix.Hij[(row(H.ci.matrix.Hij)) == (0.5 * col(H.ci.matrix.Hij) + 0.5)] <- format("", width = 7, justify = "right")
        H.ci.matrix.Hij <- noquote(H.ci.matrix.Hij)
        
        ci.matrix.Hij[row(ci.matrix.Hij) == col(ci.matrix.Hij)] <- format("", width = 7, justify = "right")
        ci.matrix.Hij <- noquote(ci.matrix.Hij)
        
        output.matrix.Hi <- matrix(NA, J, 3)
        output.matrix.Hi[, 1] <- format(formatC(round(Hi, 3), digits = 3, format = "f"), width = 7, justify = "right")
        output.matrix.Hi[, 2] <- format(paste("(", formatC(round(se.Hi, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
        if(type.ci == "WB") {
          output.matrix.Hi[, 3] <- format(paste("[", formatC(round(Hi + qnorm((1 - ci) / 2) * se.Hi, 3), digits = 3, format = "f"), 
                                                ", ", formatC(round(Hi - qnorm((1 - ci) / 2) * se.Hi, 3), digits = 3, format = "f"), "] ", sep = ""), width = 7, justify = "right")
        } else {
          output.matrix.Hi[, 3] <- format(paste("[", formatC(round(1 - exp(gHi - qnorm((1 - ci) / 2) * sqrt(VgHi)), 3), digits = 3, format = "f"), 
                                                ", ", formatC(round(1 - exp(gHi + qnorm((1 - ci) / 2) * sqrt(VgHi)), 3), digits = 3, format = "f"), "] ", sep = ""), width = 7, justify = "right")
        }
        dimnames(output.matrix.Hi) <- list(labels, c("Item H", "se", paste(ci*100, "% ci", sep = "")))
        output.matrix.Hi <- noquote(output.matrix.Hi)
        
        output.matrix.H <- matrix(NA, 1, 3)
        output.matrix.H[, 1] <- format(formatC(round(H, 3), digits = 3, format = "f"), width = 7, justify = "right")
        output.matrix.H[, 2] <- format(paste("(", formatC(round(se.H, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
        if(type.ci == "WB") {
          output.matrix.H[, 3] <- format(paste("[", formatC(round(H + qnorm((1 - ci) / 2) * se.H, 3), digits = 3, format = "f"), 
                                               ", ", formatC(round(H - qnorm((1 - ci) / 2) * se.H, 3), digits = 3, format = "f"), "]", sep = ""), width = 7, justify = "right")
        } else {
          output.matrix.H[, 3] <- format(paste("[", formatC(round(1 - exp(gH - qnorm((1 - ci) / 2) * sqrt(VgH)), 3), digits = 3, format = "f"), 
                                               ", ", formatC(round(1 - exp(gH + qnorm((1 - ci) / 2) * sqrt(VgH)), 3), digits = 3, format = "f"), "]", sep = ""), width = 7, justify = "right")
        }
        dimnames(output.matrix.H) <- list("", c("Scale H", "se", paste(ci*100, "% ci", sep = "")))
        output.matrix.H <- noquote(output.matrix.H)
        
        if(se & computeCI) {
          print.list  <- list(Hij = output.matrix.Hij, ci.Hij = ci.matrix.Hij, Hi = output.matrix.Hi, H = output.matrix.H)
          return.list <- list(Hij = output.matrix.Hij, ci.Hij = ci.matrix.Hij, Hi = output.matrix.Hi, H = output.matrix.H, covHij = ACM.Hij, covHi = ACM.Hi, covH = ACM.H)
        } else if(!se) {
          print.list  <- list(Hij = H.ci.matrix.Hij, Hi = output.matrix.Hi[, -2], H = output.matrix.H[, -2])
          return.list <- list(Hij = H.ci.matrix.Hij, Hi = output.matrix.Hi[, -2], H = output.matrix.H[, -2], covHij = ACM.Hij, covHi = ACM.Hi, covH = ACM.H)
        } else {
          print.list  <- list(Hij = output.matrix.Hij, Hi = output.matrix.Hi[, 1:2], H = output.matrix.H[, 1:2]) 
          return.list <- list(Hij = output.matrix.Hij, Hi = output.matrix.Hi[, 1:2], H = output.matrix.H[, 1:2],  covHij = ACM.Hij, covHi = ACM.Hi, covH = ACM.H)
        }
      } 
    }  else {
     if(computeCI) {
       ci.Hij <- list("Lower" = 1 - exp(gHij - qnorm((1 - ci) / 2) * sqrt(VgHij)), 
                      "Upper" = 1 - exp(gHij + qnorm((1 - ci) / 2) * sqrt(VgHij)))
       
       ci.Hi <- data.frame("Lower" = 1 - exp(gHi - qnorm((1 - ci) / 2) * sqrt(VgHi)), 
                           "Upper" = 1 - exp(gHi + qnorm((1 - ci) / 2) * sqrt(VgHi)))
       
       ci.H <- data.frame("Lower" = 1 - exp(gH - qnorm((1 - ci) / 2) * sqrt(VgH)), 
                          "Upper" = 1 - exp(gH + qnorm((1 - ci) / 2) * sqrt(VgH)))
     }
     if (se & computeCI) {
       print.list  <- list(Hij = Hij, se.Hij = se.Hij, ci.Hij = ci.Hij, Hi = Hi, se.Hi = se.Hi, ci.Hi = ci.Hi, H = H, se.H = se.H, ci.H = ci.H)
       return.list <- list(Hij = Hij, se.Hij = se.Hij, ci.Hij = ci.Hij, Hi = Hi, se.Hi = se.Hi, ci.Hi = ci.Hi, H = H, se.H = se.H, ci.H = ci.H, covHij = ACM.Hij, covHi = ACM.Hi, covH = ACM.H)
     } else if(!se & computeCI) {
       print.list  <- list(Hij = Hij, ci.Hij = ci.Hij, Hi = Hi, ci.Hi = ci.Hi, H = H, ci.H = ci.H)
       return.list <- list(Hij = Hij, ci.Hij = ci.Hij, Hi = Hi, ci.Hi = ci.Hi, H = H, ci.H = ci.H, covHij = ACM.Hij, covHi = ACM.Hi, covH = ACM.H)
     } else if(se & !computeCI){
       print.list  <- list(Hij = Hij, se.Hij = se.Hij, Hi = Hi, se.Hi = se.Hi, H = H, se.H = se.H)
       return.list <- list(Hij = Hij, se.Hij = se.Hij, Hi = Hi, se.Hi = se.Hi, H = H, se.H = se.H, covHij = ACM.Hij, covHi = ACM.Hi, covH = ACM.H)
     }
     if (!se & !computeCI) {
      print.list  <- list(Hij = Hij, Hi = Hi, H = H)
      return.list <- list(Hij = Hij, Hi = Hi, H = H)
     } 
   }
    if (!is.null(group.var)) {
      group.item <- length(return.list) + 1
      return.list[[group.item]] <- list()
      names(return.list)[[group.item]] <- "Groups"
      group.var <- apply(as.matrix(group.var), 1, paste, sep = "", collapse = "/")
      group.names <- sort(unique(group.var))
      K <- length(group.names)
      for (group in 1:K) {
        X. <- X[group.var == group.names[group], ]
        if(!is.null(level.two.var)) {
          level.two.var. <- level.two.var[group.var == group.names[group]]
          S. <- unique(level.two.var.)
          LS. <- length(S.)
          Rs. <- table(level.two.var.)
          if(any(Rs. == 1)){ 
            warning('For at least one subgroup (level.two.var) in group ', group.names[group], ' there is only 1 respondent. The scalability coefficients are computed without this (these) group(s).') 
            cases. <- !(level.two.var. %in% which(Rs. == 1))
            X. <- X.[cases., ]
            level.two.var. <- level.two.var.[cases.]
            #group.var. <- group.var.[cases.]
            Rs. <- as.numeric(table(level.two.var.))
            LS. <- length(Rs.)
            S. <- 1:LS.
            level.two.var. <- rep(S., Rs.)
          }
        }
        
        if (length(X.) == ncol(X)) {
          warning(paste("No scalability coefficients computed for group", group.names[group], ". Group contains less than two cases."))
        }
        else {
          X. <- check.data(X.)
          return.list[[group.item]][[group]] <- list()
          names(return.list[[group.item]])[[group]] <- group.names[group]
          N. <- nrow(X.)
          if (any(apply(X., 2, var) < eps)) warning(paste("In group", group.names[group], ", some variables have zero variance"))
          n. <- as.matrix(table(apply(X., 1, paste, collapse = "")))
          lab.n. <- matrix(names(table(apply(X., 1, paste, collapse = ""))))
          lab.b. <- apply(all.patterns(2, g), 2, paste, collapse = "")
          Bi. <- substr(lab.b., 1, 1)
          Bj. <- substr(lab.b., 2, 2)
          R. <- t(apply(lab.n., 1, string2integer))
          r. <- length(lab.n.)
          U. <- list()
          for (j in 1:J) U.[[j]] <- tabulate(X.[, j] + 1, g)
          WA. <- list()
          WY. <- list()
          WE. <- list()
          WF. <- list()
          for (i in 1:J) {
            WA.[[i]] <- list()
            WY.[[i]] <- list()
            WE.[[i]] <- list()
            WF.[[i]] <- list()
            for (j in i:J) if (j > i) {
              A1a. <- NULL
              for (a in 0:(g - 1)) for (b in 0:(g - 1)) A1a. <- rbind(A1a., as.numeric(R.[, i] == a & R.[, j] == b))
              WA.[[i]][[j]] <- W[[i]][[j]] %*% A1a.
              Eij. <- matrix(U.[[i]][as.numeric(Bi.) + 1], nrow = g^2, ncol = 1) * matrix(U.[[j]][as.numeric(Bj.) + 1], nrow = g^2, ncol = 1)/N.
              Y22. <- cbind(outer(Bi., lab.u, "=="), outer(Bj., lab.u, "==")) * matrix(Eij., nrow = g^2, ncol = 2 * g)
              Ri. <- substr(lab.n., i, i)
              Rj. <- substr(lab.n., j, j)
              Z2. <- rbind(outer(lab.u, Ri., "==")[, , 1], outer(lab.u, Rj., "==")[, , 1]) * c(1/U.[[i]], 1/U.[[j]])
              Z2.[is.nan(Z2.)] <- 1/eps
              YZ2. <- Y22. %*% Z2. - Eij. %*% matrix(1/N., 1, r.)
              WY.[[i]][[j]] <- W[[i]][[j]] %*% YZ2.
              Fij. <- complete.observed.frequencies(X.[, c(i, j)], 2, g)
              WF.[[i]][[j]] <- W[[i]][[j]] %*% Fij.
              WE.[[i]][[j]] <- W[[i]][[j]] %*% Eij.
            }
          }
          g3. <- matrix(c(unlist(WF.[[1]][[2]]), unlist(WF.), unlist(WE.)), nrow = 2 * P + 1, byrow = TRUE)
          A4 <- rbind(matrix(c(1, -1, rep(0, (J * (J - 1)) - 1)), 1, (J * (J - 1)) + 1), cbind(matrix(0, P, 1), diag(P), -1 * diag(P)))
          A5 <- cbind(matrix(1, P, 1), -1 * diag(P))
          g4. <- phi(A4, g3., "log")
          g5. <- phi(A5, g4., "exp")
          G3. <- matrix(c(unlist(WA.[[1]][[2]]), unlist(WA.), unlist(WY.)), nrow = 2 * P + 1, byrow = TRUE)
          G4. <- dphi(A4, g3., G3., "log")
          G5ij. <- dphi(A5, g4., G4., "exp")
          Hij. <- matrix(0, J, J)
          Hij.[lower.tri(Hij.)] <- g5.
          Hij. <- Hij. + t(Hij.)
          if (P > 1) G3. <- rbind(apply(G3.[2:(P + 1), ], 2, sum), apply(G3.[2:(P + 1), ], 2, sum), apply(G3.[(P + 2):(2 * P + 1), ], 2, sum))
          g3. <- matrix(c(sum(g3.[2:(P + 1), ]), sum(g3.[2:(P + 1), ]), sum(g3.[(P + 2):(2 * P + 1), ])), ncol = 1)
          A4 <- cbind(matrix(1, 2, 1), -1 * diag(2))
          A5 <- matrix(c(1, -1), 1, 2)
          g4. <- phi(A4, g3., "log")
          g5. <- phi(A5, g4., "exp")
          G4. <- dphi(A4, g3., G3., "log")
          G5. <- dphi(A5, g4., G4., "exp")
          H. <- g5.
          G3. <- matrix(0, 2 * J + 1, r.)
          g3. <- matrix(0, 2 * J + 1, 1)
          for (j in 1:J) for (k in 1:J) if (k > j) {
            g3.[j + 1, ] <- g3.[j + 1, ] + WF.[[j]][[k]]
            g3.[J + j + 1, ] <- g3.[J + j + 1, ] + WE.[[j]][[k]]
            G3.[j + 1, ] <- G3.[j + 1, ] + WA.[[j]][[k]]
            G3.[J + j + 1, ] <- G3.[J + j + 1, ] + WY.[[j]][[k]]
          }
          else {
            if (k < j) {
              g3.[j + 1, ] <- g3.[j + 1, ] + WF.[[k]][[j]]
              g3.[J + j + 1, ] <- g3.[J + j + 1, ] + WE.[[k]][[j]]
              G3.[j + 1, ] <- G3.[j + 1, ] + WA.[[k]][[j]]
              G3.[J + j + 1, ] <- G3.[J + j + 1, ] + WY.[[k]][[j]]
            }
          }
          g3.[1, ] <- g3.[2, ]
          G3.[1, ] <- G3.[2, ]
          A4 <- rbind(matrix(c(1, -1, rep(0, (J * 2) - 1)), 1, (J * 2) + 1), cbind(matrix(0, J, 1), diag(J), -1 * diag(J)))
          A5 <- cbind(matrix(1, J, 1), -1 * diag(J))
          g4. <- phi(A4, g3., "log")
          g5. <- phi(A5, g4., "exp")
          G4. <- dphi(A4, g3., G3., "log")
          G5i. <- dphi(A5, g4., G4., "exp")
          Hi. <- matrix(g5.)
          if (se) {
            if(!is.null(level.two.var)) {
              pats.n. <- matrix(apply(X., 1, paste, collapse = ""))
              n.p. <- length(n.)
              covps. <- matrix(0, n.p., n.p.)
              p. <- matrix(0, n.p.)
              for(s in S.) {
                pt. <- rowSums(outer(lab.n., pats.n.[level.two.var. == s], "==")) 
                prows. <- which(pt. > 0)
                pt. <- pt.[prows.]
                p.[prows.] <- p.[prows.] + pt.
                covps.[prows., prows.] <- covps.[prows., prows.] + (pt. %*% t(pt.)) / Rs.[s]
              }
              p. <- as.numeric(p. / sum(Rs.)) 
              covps. <- covps. / sum(Rs.) - (p. %*% t(p.)) 
              covp. <- (diag(length(p.)) * p. - (p. %*% t(p.))) 
              nu. <- LS. / sum(1 / Rs.)
              
              # Variance covariance matrix needed for delta method 
              covtot. <- LS. * nu. * covp. + LS. * nu. * (nu. - 1) * covps. 
              
              se.Hij. <- matrix(0, J, J)
              ACM.Hij. = G5ij. %*% (covtot. %*% t(G5ij.))
              se.Hij.[lower.tri(se.Hij.)] <- sqrt(diag(ACM.Hij.))
              se.Hij. <- se.Hij. + t(se.Hij.)
              dimnames(se.Hij.) <- dimnames(Hij.) <- list(labels, labels)
              ACM.Hi. = G5i. %*% (covtot. %*% t(G5i.))
              se.Hi. <- matrix(sqrt(diag(ACM.Hi.)))
              dimnames(se.Hi.)[[1]] <- dimnames(Hi.)[[1]] <- labels
              ACM.H. <- G5. %*% (covtot. %*% t(G5.))
              se.H. <- sqrt(diag(ACM.H.))
            } else {
              se.Hij. <- matrix(0, J, J)
              ACM.Hij. = G5ij. %*% (as.numeric(n.) * t(G5ij.))
              se.Hij.[lower.tri(se.Hij.)] <- sqrt(diag(ACM.Hij.))
              se.Hij. <- se.Hij. + t(se.Hij.)
              dimnames(se.Hij.) <- dimnames(Hij.) <- list(labels, labels)
              ACM.Hi. = G5i. %*% (as.numeric(n.) * t(G5i.))
              se.Hi. <- matrix(sqrt(diag(ACM.Hi.)))
              dimnames(se.Hi.)[[1]] <- dimnames(Hi.)[[1]] <- labels
              ACM.H. <- G5. %*% (as.numeric(n.) * t(G5.))
              se.H. <- sqrt(diag(ACM.H.))
            }
            
            #se.Hij. <- matrix(0, J, J)
            #ACM.Hij = G5ij. %*% (as.numeric(n.) * t(G5ij.))
            #se.Hij.[lower.tri(se.Hij.)] <- sqrt(diag(ACM.Hij))
            #se.Hij. <- se.Hij. + t(se.Hij.)
            #dimnames(se.Hij.) <- dimnames(Hij.) <- list(labels, 
            #                                            labels)
            #ACM.Hi = G5i. %*% (as.numeric(n.) * t(G5i.))
            #se.Hi. <- matrix(sqrt(diag(ACM.Hi)))
            #dimnames(se.Hi)[[1]] <- dimnames(Hi.)[[1]] <- labels
            #ACM.H = G5. %*% (as.numeric(n.) * t(G5.))
            #se.H. <- sqrt(diag(ACM.H))
          }
          if (nice.output && se) {
            output.matrix.Hij. <- matrix(NA, J, J * 2)
            for (j in 2 * (1:J)) {
              output.matrix.Hij.[, j - 1] <- format(paste(" ", formatC(round(Hij.[, j/2], 3), digits = 3, format = "f"), " ", sep = ""), width = 7, justify = "right")
              output.matrix.Hij.[, j] <- format(paste("(", formatC(round(se.Hij.[, j/2], 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
            }
            dimnames(output.matrix.Hij.)[[1]] <- labels
            dimnames(output.matrix.Hij.)[[2]] <- new.labels
            output.matrix.Hij.[row(output.matrix.Hij.) == 0.5 * col(output.matrix.Hij.)] <- format("", width = 7, justify = "right")
            output.matrix.Hij.[(row(output.matrix.Hij.)) == (0.5 * col(output.matrix.Hij.) + 0.5)] <- format("", width = 7, justify = "right")
            output.matrix.Hij. <- noquote(output.matrix.Hij.)
            output.matrix.Hi. <- matrix(NA, J, 2)
            output.matrix.Hi.[, 1] <- format(formatC(round(Hi., 3), digits = 3, format = "f"), width = 7, justify = "right")
            output.matrix.Hi.[, 2] <- format(paste("(", formatC(round(se.Hi., 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
            dimnames(output.matrix.Hi.) <- list(labels, c("Item H", "se"))
            output.matrix.Hi. <- noquote(output.matrix.Hi.)
            output.matrix.H. <- matrix(NA, 1, 2)
            output.matrix.H.[, 1] <- format(formatC(round(H., 3), digits = 3, format = "f"), width = 7, justify = "right")
            output.matrix.H.[, 2] <- format(paste("(", formatC(round(se.H., 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
            dimnames(output.matrix.H.) <- list("", c("Scale H", "se"))
            output.matrix.H. <- noquote(output.matrix.H.)
            return.list[[group.item]][[group]] <- list(Hij = output.matrix.Hij., Hi = output.matrix.Hi., H = output.matrix.H.)
          }
          else {
            if (se) 
              return.list[[group.item]][[group]] <- list(Hij = Hij., se.Hij = se.Hij., Hi = Hi., se.Hi = se.Hi., H = H., se.H = se.H.)
            if (!se) 
              return.list[[group.item]][[group]] <- list(Hij = Hij., Hi = Hi., H = H.)
          }
        }
      }
    print.list <- return.list  
    }
  }
  if(results) print(print.list)
  invisible(return.list)    
}

