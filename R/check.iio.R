# Adapted May 28, 2024
"check.iio" <- 
  function (X, method = "MIIO", minvi = default.minvi, minsize = default.minsize, 
            alpha = .05, item.selection = TRUE, verbose = FALSE, fixed.item.order = NULL,
            level.two.var = NULL) 
  {
    # CHECK DATA
    X <- check.data(X)
    
    # READ
    J <- ncol(X)
    N <- nrow(X)
    ncat <- max(X) + 1
    if (!is.null(level.two.var)) {
      X <- X[order(level.two.var), ]
      level.two.var <- sort(level.two.var)
      if (method != "MIIO") {
        warning("Only method MIIO is available for level.two.var, method is changed to MIIO.")
        method <- "MIIO"
      }
      if(!is.null(fixed.item.order)) {
        if(any(table(fixed.item.order) > 1)) {
          warning("level.two.var can currently not be combined with clustered items. fixed.item.order is ingored. Please let the package authors know if you need to combine them.")
          fixed.item.order <- NULL
        }
      }
    }
  
    # DETERMINE METHOD
    if (substr(method, 1, 1) == "I" || substr(method, 1, 1) == 
        "i") {
      method = "IT"
    } else if (substr(method, 1, 2) == "MS" || substr(method, 1, 
              2) == "ms" || substr(method, 1, 1) == "C" || substr(method, 
                        1, 1) == "c")  {
      method = "MSCPM"
    } else {
      method <- "MIIO"
    }
      
    # DETERMINE MINVI and MINSIZE
    default.minvi <- ifelse(method == "MIIO", (ncat - 1) * 0.03, 
                            0.03)
    default.minsize <- ifelse(N >= 500, floor(N/10), floor(N/5))
    default.minsize <- ifelse(N <= 250, floor(N/3), default.minsize)
    default.minsize <- ifelse(N < 150, 50, default.minsize)
    
    # STOP IF THERE ARE NOT ENOUGH ITEMS OR RESPONDENTS
    if (N < minsize) 
      stop("Sample size less than Minsize")
    if (J < 3) 
      stop("Less than 3 items. Restscore cannot be computed")
    
    # INITIAL VALUES
    results <- list()
    clusters <- FALSE
    vi.matrix <- matrix(0, J, J)
    #item.order <- rev(order(colMeans(X)))
    
    if (is.null(dimnames(X)[[2]])) 
      dimnames(X)[[2]] <- paste("V", 1:ncol(X), sep = "")
    if(is.null(fixed.item.order)) {
      item.order <- rev(order(colMeans(X)))
    } else if  (!is.numeric(fixed.item.order)) {
     warning("fixed.item.order is not numeric and is ignored. Please rerun analyses with a numeric input for fixed.item.order")
      item.order <- rev(order(colMeans(X)))
    } else {
      #fixed.item.order <- as.factor(fixed.item.order)
      #orig.item.order <- fixed.item.order
      #fixed.item.order <- as.numeric(fixed.item.order)
      if(any(table(fixed.item.order) > 1)) {
        if(length(unique(table(fixed.item.order))) > 1) {
          warning("Clusters have an unequal number of items, estimation of HTB may be biased. Use function MLcoefH(cbind(fixed.item.order, t(X)), se = F)$H as alternative. Note that computation may be slow, hence it is not incorporated in check.iio.")
        }
        if (method != "MIIO") {
          warning("Only method MIIO is available for clustered items in fixed.item.order, method is changed to MIIO.")
          method <- "MIIO"
        }
        if(item.selection == TRUE) {
          warning("For clustered items and fixed item order item selection is not available. Argument item.selection is set to FALSE.")
          item.selection = FALSE
        }
        
        clusters <- TRUE
        #orig.cluster.order <- orig.item.order
        item.cluster.order <- fixed.item.order
        cmns <- colMeans(X)
        item.order <- order(fixed.item.order, -cmns)
        item.cluster.order <- item.cluster.order[order(fixed.item.order)]
        item.cluster.order.o <- item.cluster.order
        #orig.cluster.order <- orig.cluster.order[order(fixed.item.order)]
        
      } else {
        item.order <- order(fixed.item.order)
      }
      
    }
    
    X <- X[, item.order]
    
    if(clusters) {
      Xaggclus <- (aggregate(t(X), list(item.cluster.order), mean)) 
      C.labels <- unique(item.cluster.order)#Xaggclus[, 1]
      Xaggclus <- t(Xaggclus[, -1])
      colnames(Xaggclus) <- C.labels
      nC <- ncol(Xaggclus)
      #HiC <- coefHTiny(Xaggclus)$Hi
      HiC <- 0
    } 
    
    I.labels <- dimnames(X)[[2]]
    dimnames(vi.matrix) <- list(I.labels, I.labels)
    Hi <- coefHTiny(X)$Hi
    g <- 0
    gg <- 0
    h <- 0
    i <- 1
    j <- 2
    k <- 0
    Ll1 <- NULL
    item.mean <- colMeans(X)
    # METHOD IT
    if (method == "IT") {
      ii <- jj <- kk <- 0
      uv <- matrix(0, nrow = choose(ncat, 2), ncol = 2)
      for (ii in 0:(ncat - 2)) for (jj in (ii + 1):(ncat - 
                                                    1)) {
        kk <- kk + 1
        uv[kk, ] <- c(ii, jj)
      }
      for (i in 1:(J - 1)) {
        R <- as.matrix(X) %*% (matrix(1, J, J) - diag(J)) - 
          X[, i]
        R[, i] <- 0
        for (j in (i + 1):J) {
          k <- k + 1
          rvm <- nrow(uv) + 1
          violation.matrix <- matrix(0, nrow = rvm, ncol = 8)
          dimnames(violation.matrix) <- list(c(t(paste("P(X", 
                                                       i, "=", uv[, 2], ",X", j, "=", uv[, 1], ") P(X", 
                                                       i, "=", uv[, 1], ",X", j, "=", uv[, 2], ")", 
                                                       sep = "")), "Total"), c("#ac", "#vi", "#vi/#ac", 
                                                                               "maxvi", "sum", "sum/#ac", "X^2 max", "#X^2 sig"))
          results[[k]] <- list()
          results[[k]][[1]] <- list()
          results[[k]][[1]][1] <- I.labels[i]
          results[[k]][[1]][2] <- I.labels[j]
          sorted.R <- sort(R[, j])
          group <- max(which(sorted.R == sorted.R[minsize]))
          repeat {
            if (N - max(group) < minsize) 
              break
            group <- c(group, max(which(sorted.R == sorted.R[minsize + 
                                                               max(group)])))
          }
          group <- group[-length(group)]
          summary.matrix <- matrix(nrow = length(group) + 
                                     1, ncol = 5 + 2 * rvm)
          dimnames(summary.matrix)[[2]] <- c("Group", "Lo", 
                                             "Hi", "N", "n", paste("E(X", i, ")", sep = ""), 
                                             paste("E(X", j, ")", sep = ""), paste("P(X", 
                                                                                   i, "=", uv[, 1], ",X", j, "=", uv[, 2], ")", 
                                                                                   sep = ""), paste("P(X", i, "=", uv[, 2], 
                                                                                                    ",X", j, "=", uv[, 1], ")", sep = ""))
          summary.matrix[, 1] <- 1:nrow(summary.matrix)
          summary.matrix[, 4] <- c(group, N) - c(0, group)
          group <- c(sorted.R[group], max(sorted.R))
          L <- length(group)
          summary.matrix[, 3] <- group
          summary.matrix[, 2] <- c(min(sorted.R), group[-L] + 
                                     1)
          member <- apply(1 - outer(R[, j], group, "<="), 
                          1, sum) + 1
          for (g in 1:L) {
            summary.matrix[g, 6] <- mean(X[member == g, 
                                           i])
            summary.matrix[g, 7] <- mean(X[member == g, 
                                           j])
            u <- 0
            for (u in 1:nrow(uv)) {
              summary.matrix[g, 7 + u] <- length(X[member == 
                                                     g & X[, i] == uv[u, 1] & X[, j] == uv[u, 
                                                                                           2], ])/J
              summary.matrix[g, 7 + nrow(uv) + u] <- length(X[member == 
                                                                g & X[, i] == uv[u, 2] & X[, j] == uv[u, 
                                                                                                      1], ])/J
            }
          }
          summary.matrix[, 5] <- apply(summary.matrix[, 
                                                      8:(7 + 2 * nrow(uv))], 1, sum)
          results[[k]][[2]] <- summary.matrix
          ac <- violation.matrix[1:(rvm - 1), 1] <- L - 
            1
          for (g in 1:nrow(uv)) {
            n.uv <- summary.matrix[, 7 + g]
            n.vu <- summary.matrix[, 7 + nrow(uv) + g]
            Ng <- summary.matrix[, 4]
            d <- (n.uv - n.vu)/Ng
            d[d <= minvi] <- 0
            vi <- length(d[d > minvi/2])
            sum.vi <- sum(d)
            chi.statistic <- rep(0, L)
            chi.pvalue <- rep(1, L)
            if (any(d > 0)) {
              for (gg in 1:L) {
                if (d[gg] > 0) {
                  e <- (n.uv[gg] + n.vu[gg])/2
                  chi.statistic[gg] <- (n.uv[gg] - e)^2/e + 
                    (n.vu[gg] - e)^2/e
                  chi.pvalue[gg] <- 1 - pchisq(chi.statistic[gg], 
                                               1)
                }
              }
            }
            violation.matrix[g, 2:8] <- c(vi, vi/ac, max(d), 
                                          sum.vi, sum.vi/ac, max(chi.statistic), sum(chi.pvalue < 
                                                                                       alpha))
            vi.matrix[i, j] <- vi.matrix[i, j] + sum(chi.pvalue < 
                                                       alpha)
          }
          violation.matrix[rvm, c(1, 2, 5, 8)] <- apply(violation.matrix[1:(rvm - 
                                                                              1), c(1, 2, 5, 8)], 2, sum)
          violation.matrix[rvm, 3] <- violation.matrix[rvm, 
                                                       2]/violation.matrix[rvm, 1]
          violation.matrix[rvm, 6] <- violation.matrix[rvm, 
                                                       5]/violation.matrix[rvm, 1]
          violation.matrix[rvm, c(4, 7)] <- apply(violation.matrix[1:(rvm - 
                                                                        1), c(4, 7)], 2, max)
          results[[k]][[3]] <- violation.matrix
        }
      }
    }
    # METHOD MSCPM
    if (method == "MSCPM") {
      for (i in 1:(J - 1)) {
        R <- as.matrix(X) %*% (matrix(1, J, J) - diag(J)) - 
          X[, i]
        R[, i] <- 0
        for (j in (i + 1):J) {
          k <- k + 1
          rvm <- (ncat - 1) + 1
          violation.matrix <- matrix(0, nrow = rvm, ncol = 8)
          dimnames(violation.matrix) <- list(c(t(paste("P(X", 
                                                       i, ">=", 1:(ncat - 1), ")  P(X", j, ">=", 1:(ncat - 
                                                                                                      1), ")", sep = "")), "Total"), c("#ac", "#vi", 
                                                                                                                                       "#vi/#ac", "maxvi", "sum", "sum/#ac", "zmax", 
                                                                                                                                       "#zsig"))
          results[[k]] <- list()
          results[[k]][[1]] <- list()
          results[[k]][[1]][1] <- I.labels[i]
          results[[k]][[1]][2] <- I.labels[j]
          sorted.R <- sort(R[, j])
          group <- max(which(sorted.R == sorted.R[minsize]))
          repeat {
            if (N - max(group) < minsize) 
              break
            group <- c(group, max(which(sorted.R == sorted.R[minsize + 
                                                               max(group)])))
          }
          group <- group[-length(group)]
          summary.matrix <- matrix(nrow = length(group) + 
                                     1, ncol = 6 + 2 * (ncat - 1))
          dimnames(summary.matrix)[[2]] <- c("Group", "Lo", 
                                             "Hi", "N", paste("E(X", i, ")", sep = ""), 
                                             paste("E(X", j, ")", sep = ""), paste("P(X", 
                                                                                   i, ">=", 1:(ncat - 1), ")", sep = ""), paste("P(X", 
                                                                                                                                j, ">=", 1:(ncat - 1), ")", sep = ""))
          summary.matrix[, 1] <- 1:nrow(summary.matrix)
          Ng <- summary.matrix[, 4] <- c(group, N) - c(0, 
                                                       group)
          group <- c(sorted.R[group], max(sorted.R))
          L <- length(group)
          summary.matrix[, 3] <- group
          summary.matrix[, 2] <- c(min(sorted.R), group[-L] + 
                                     1)
          member <- apply(1 - outer(R[, j], group, "<="), 
                          1, sum) + 1
          for (g in 1:L) {
            summary.matrix[g, 5] <- mean(X[member == g, 
                                           i])
            summary.matrix[g, 6] <- mean(X[member == g, 
                                           j])
            freqi <- tabulate(X[member == g, i] + 1, ncat)
            freqj <- tabulate(X[member == g, j] + 1, ncat)
            cum.freqi <- rev(cumsum(rev(freqi))/Ng[g])
            cum.freqj <- rev(cumsum(rev(freqj))/Ng[g])
            summary.matrix[g, 7:(5 + ncat)] <- cum.freqi[2:ncat]
            summary.matrix[g, (6 + ncat):(4 + 2 * ncat)] <- cum.freqj[2:ncat]
          }
          results[[k]][[2]] <- summary.matrix
          ac <- violation.matrix[1:(rvm - 1), 1] <- L - 
            1
          for (g in 1:(ncat - 1)) {
            p.1 <- summary.matrix[, 6 + g]
            p.2 <- summary.matrix[, 5 + ncat + g]
            d <- p.2 - p.1
            d[d <= minvi] <- 0
            vi <- length(d[d > minvi/2])
            sum.vi <- sum(d)
            z.statistic <- rep(0, L)
            if (any(d > 0)) {
              for (gg in 1:L) {
                if (d[gg] > 0) {
                  Xgg <- X[member == gg, ]
                  f.01 <- length(which(Xgg[, i] >= g & 
                                         Xgg[, j] < (g + ncat)))
                  f.10 <- length(which(Xgg[, i] < g & Xgg[, 
                                                          j] >= (g + ncat)))
                  f.k <- min(f.01, f.10)
                  f.n <- f.01 + f.10
                  f.b <- ((2 * f.k + 1 - f.n)^2 - 10 * 
                            f.n)/(12 * f.n)
                  z.statistic[gg] <- abs(sqrt(2 * f.k + 
                                                2 + f.b) - sqrt(2 * f.n - 2 * f.k + 
                                                                  f.b))
                }
              }
            }
            violation.matrix[g, 2:8] <- c(vi, vi/ac, max(d), 
                                          sum.vi, sum.vi/ac, max(z.statistic), length(z.statistic[abs(z.statistic) > 
                                                                                                    qnorm(1 - alpha)]))
            vi.matrix[i, j] <- vi.matrix[i, j] + length(z.statistic[abs(z.statistic) > 
                                                                      qnorm(1 - alpha)])
          }
          violation.matrix[rvm, c(1, 2, 5, 8)] <- apply(violation.matrix[1:(rvm - 
                                                                              1), c(1, 2, 5, 8)], 2, sum)
          violation.matrix[rvm, 3] <- violation.matrix[rvm, 
                                                       2]/violation.matrix[rvm, 1]
          violation.matrix[rvm, 6] <- violation.matrix[rvm, 
                                                       5]/violation.matrix[rvm, 1]
          violation.matrix[rvm, c(4, 7)] <- apply(violation.matrix[1:(rvm - 
                                                                        1), c(4, 7)], 2, max)
          results[[k]][[3]] <- violation.matrix
        }
      }
    }
    dich = 2
    # METHOD MIIO
    if (method == "MIIO") {
      for (i in 1:(J - 1)) {
        if(clusters) {
          iclus <- (item.cluster.order != item.cluster.order[i])
          Ragg <- as.matrix(Xaggclus) %*% (matrix(1, nC, nC) - diag(nC))# - Xaggclus[, i] 
          #Ragg[, i * jclus] <- Ragg[, i * jclus] - Xaggclus[, i * jclus] 
          R <- Ragg[, item.cluster.order]
          R[, iclus] <- R[, iclus] - Xaggclus[, item.cluster.order][, i]
          #Ragg <- Ragg - Xaggclus[, i] 
        } else {
        R <- as.matrix(X) %*% (matrix(1, J, J) - diag(J)) - X[, i] #R_(ij) for j != i
        }
        R[, i] <- 0
        for (j in (i + 1):J) {
          k <- k + 1
          rvm <- 2
          violation.matrix <- matrix(0, nrow = 1, ncol = 8)
          dimnames(violation.matrix) <- list(c(t(paste("E(X", 
                  i, ")  E(X", j, ")", sep = ""))), c("#ac", 
                   "#vi", "#vi/#ac", "maxvi", "sum", "sum/#ac", 
                    ifelse(ncat == dich, "zmax", "tmax"), ifelse(ncat == 
                                 dich, "#zsig", "#tsig")))
          results[[k]] <- list()
          results[[k]][[1]] <- list()
          results[[k]][[1]][1] <- I.labels[i]
          results[[k]][[1]][2] <- I.labels[j]
          sorted.R <- sort(R[, j])
          group <- max(which(sorted.R == sorted.R[minsize]))
          repeat {
            if (N - max(group) < minsize) 
              break
            group <- c(group, max(which(sorted.R == sorted.R[minsize + 
                                                               max(group)])))
          }
          group <- group[-length(group)]
          summary.matrix <- matrix(nrow = length(group) + 
                                     1, ncol = 8)
          dimnames(summary.matrix)[[2]] <- c("Group", "Lo", 
                     "Hi", "N", paste("E(X", i, ")", sep = ""), 
                     paste("E(X", j, ")", sep = ""), paste("SD(X", 
                                  i, ")", sep = ""), paste("SD(X", j, ")", 
                                          sep = ""))
          summary.matrix[, 1] <- 1:nrow(summary.matrix)
          summary.matrix[, 4] <- c(group, N) - c(0, group)
          group <- c(sorted.R[group], max(sorted.R))
          L <- Ll1[j] <- length(group)
          summary.matrix[, 3] <- group
          summary.matrix[, 2] <- c(min(sorted.R), group[-L] + 
                                     1)
          member <- apply(1 - outer(R[, j], group, "<="), 
                          1, sum) + 1
          for (g in 1:L) {
            summary.matrix[g, 5] <- mean(X[member == g, 
                                           i])
            summary.matrix[g, 6] <- mean(X[member == g, 
                                           j])
            summary.matrix[g, 7] <- sd(X[member == g, i])
            summary.matrix[g, 8] <- sd(X[member == g, j])
          } # end g-loop
          results[[k]][[2]] <- summary.matrix
          ac <- violation.matrix[1:(rvm - 1), 1] <- L - 
            1
          E.1 <- summary.matrix[, 5]
          E.2 <- summary.matrix[, 6]
          d <- E.2 - E.1
          d[d <= minvi] <- 0
          vi <- length(d[d > minvi/2])
          sum.vi <- sum(d)
          t.statistic <- rep(0, L)
          t.pvalue <- rep(1, L)
          if (any(d > 0)) {
            for (gg in 1:L) {
              if (d[gg] > 0) {
                if (ncat > dich) {
                  tt <- t.test(X[member == gg, i], X[member == 
                                                       gg, j], alternative = "less", paired = TRUE)
                  t.statistic[gg] <- abs(tt$statistic)
                  t.pvalue[gg] <- tt$p.value
                }
                else {
                  Xgg <- X[member == gg, ]
                  f.01 <- length(which(Xgg[, i] >= 1 & 
                                         Xgg[, j] < 1))
                  f.10 <- length(which(Xgg[, i] < 1 & Xgg[, 
                                                          j] >= 1))
                  f.k <- min(f.01, f.10)
                  f.n <- f.01 + f.10
                  f.b <- ((2 * f.k + 1 - f.n)^2 - 10 * 
                            f.n)/(12 * f.n)
                  t.statistic[gg] <- abs(sqrt(2 * f.k + 
                       2 + f.b) - sqrt(2 * f.n - 2 * f.k + 
                       f.b)) # even though it is a z-statistic
                  t.pvalue[gg] <- 1 - pnorm(t.statistic[gg])
                }
              } # endif
            } # end gg-loop
          } #endif
          violation.matrix[1, 2:8] <- c(vi, vi/ac, max(d), 
                                        sum.vi, sum.vi/ac, max(t.statistic), sum(t.pvalue < 
                                                                                   alpha))
          vi.matrix[i, j] <- vi.matrix[i, j] + sum(t.pvalue < 
                                                     alpha)
          results[[k]][[3]] <- violation.matrix
        }
        }
      if (!is.null(level.two.var)) {
        vi.matrixA <- matrix(0, J, J)
        dimnames(vi.matrixA) <- list(I.labels, I.labels)
        Rs <- table(level.two.var)
        Xa <- as.matrix(aggregate(X, by = list(level.two.var), 
                                  FUN = mean)[, -1])
        Xas <- Xa
        Na <- nrow(Xa)
        resultsa <- list()
        His <- MLcoefH(cbind(level.two.var, X), se = F, nice.output = F)$Hi
        g <- 0
        gg <- 0
        h <- 0
        i <- 1
        j <- 2
        k <- 0
        for (i in 1:(J - 1)) {
          Ras <- as.matrix(Xa) %*% (matrix(1, J, J) - diag(J)) - 
            Xa[, i]
          Ras <- round(Ras, 4)
          Ras[, i] <- 0
          Ra <- Ras[rep(1:nrow(Ras), Rs), ]
          Ra[, i] <- 0
          for (j in (i + 1):J) {
            k <- k + 1
            rvm <- 2
            violation.matrix <- matrix(0, nrow = 1, ncol = 8)
            dimnames(violation.matrix) <- list(c(t(paste("E(X", 
                                                         i, ")  E(X", j, ")", sep = ""))), c("#ac", 
                                                                                             "#vi", "#vi/#ac", "maxvi", "sum", "sum/#ac", 
                                                                                             ifelse(ncat == dich, "zmax", "tmax"), ifelse(ncat == 
                                                                                                                                            dich, "#zsig", "#tsig")))
            resultsa[[k]] <- list()
            resultsa[[k]][[1]] <- list()
            resultsa[[k]][[1]][1] <- I.labels[i]
            resultsa[[k]][[1]][2] <- I.labels[j]
            sorted.Ra <- sort(round(Ra[, j], 4))
            minsizeA <- floor(nrow(X)/(Ll1[j] + 1.75))
            group <- max(which(sorted.Ra == sorted.Ra[minsizeA]))
            repeat {
              if (N - max(group) < minsizeA) 
                break
              group <- c(group, max(which(sorted.Ra == 
                                            sorted.Ra[minsizeA + max(group)])))
            }
            group <- group[-length(group)]
            summary.matrix <- matrix(nrow = length(group) + 
                                       1, ncol = 8)
            dimnames(summary.matrix)[[2]] <- c("Group", 
                                               "Lo", "Hi", "N", paste("E(X", i, ")", sep = ""), 
                                               paste("E(X", j, ")", sep = ""), paste("SD(X", 
                                                                                     i, ")", sep = ""), paste("SD(X", j, ")", 
                                                                                                              sep = ""))
            summary.matrix[, 1] <- 1:nrow(summary.matrix)
            group <- c(sorted.Ra[group], max(sorted.Ra))
            L <- length(group)
            summary.matrix[, 3] <- group
            summary.matrix[, 2] <- c(min(sorted.Ra), group[-L])
            member <- apply(1 - outer(Ras[, j], group, 
                                      "<="), 1, sum) + 1
            summary.matrix[, 4] <- table(member)
            for (g in 1:L) {
              summary.matrix[g, 5] <- mean(Xas[member == 
                                                 g, i])
              summary.matrix[g, 6] <- mean(Xas[member == 
                                                 g, j])
              summary.matrix[g, 7] <- sd(Xas[member == 
                                               g, i])
              summary.matrix[g, 8] <- sd(Xas[member == 
                                               g, j])
            }
            resultsa[[k]][[2]] <- summary.matrix
            ac <- violation.matrix[1:(rvm - 1), 1] <- L - 
              1
            E.1 <- summary.matrix[, 5]
            E.2 <- summary.matrix[, 6]
            d <- E.2 - E.1
            d[d <= minvi] <- 0
            vi <- length(d[d > minvi/2])
            sum.vi <- sum(d)
            t.statistic <- rep(0, L)
            t.pvalue <- rep(1, L)
            if (any(d > 0)) {
              for (gg in 1:L) {
                if (d[gg] > 0) {
                  tt <- t.test(Xas[member == gg, i], Xas[member == 
                                                           gg, j], alternative = "less", paired = TRUE)
                  t.statistic[gg] <- abs(tt$statistic)
                  t.pvalue[gg] <- tt$p.value
                }
              }
            }
            violation.matrix[1, 2:8] <- c(vi, vi/ac, max(d), 
                                          sum.vi, sum.vi/ac, max(t.statistic), sum(t.pvalue < 
                                                                                     alpha))
            vi.matrixA[i, j] <- vi.matrixA[i, j] + sum(t.pvalue < 
                                                         alpha)
            resultsa[[k]][[3]] <- violation.matrix
          }
        }
        vi.matrixA <- vi.matrixA + t(vi.matrixA)
        VIA <- matrix(apply(sign(vi.matrixA), 1, sum))
        dimnames(VIA) <- list(I.labels, paste("step", 1:ncol(VIA)))
      }
        
        if(clusters) {
          resultsC <- list()
          vi.matrixC <- matrix(0,nC,nC)
          dimnames(vi.matrixC) <- list(C.labels,C.labels)
          g <- 0
          gg <- 0
          h <- 0
          cc <- 1
          cd <- 2
          k <- 0
          cluster.mean <- colMeans(Xaggclus)
          
          for (cc in 1:(nC-1)){
            R <- as.matrix(Xaggclus) %*% (matrix(1, nC, nC) - diag(nC)) - Xaggclus[, cc] 
            R[, cc] <- 0
            for (cd in (cc+1):nC){
              k <- k + 1
              rvm <- 2
              violation.matrix <- matrix(0, nrow = 1, ncol = 8)
              dimnames(violation.matrix) <- list(c(t(paste("E(X",cc,")  E(X",cd,")", sep = ""))), 
                                                 c("#ac", "#vi", "#vi/#ac", "maxvi", "sum", "sum/#ac",
                                                   "tmax", "#tsig"))
              resultsC[[k]] <- list()
              resultsC[[k]][[1]] <- list()
              resultsC[[k]][[1]][1] <- C.labels[cc]
              resultsC[[k]][[1]][2] <- C.labels[cd]
              
              sorted.R <- sort(R[,cd])
              group <- max(which(sorted.R==sorted.R[minsize]))
              repeat{
                if(N - max(group) < minsize)break
                group <- c(group,max(which(sorted.R==sorted.R[minsize+max(group)])))
              }
              group <- group[-length(group)]  
              summary.matrix <- matrix(nrow = length(group) + 1, ncol = 8)
              dimnames(summary.matrix)[[2]] <- c("Group", "Lo","Hi", "N", paste("E(X", cc, ")", sep = ""), paste("E(X",cd, ")", sep = ""),paste("SD(X", cc, ")", sep = ""), paste("SD(X",cd, ")", sep = ""))
              summary.matrix[, 1] <- 1:nrow(summary.matrix)
              summary.matrix[, 4] <- c(group, N) - c(0, group)
              group <- c(sorted.R[group],max(sorted.R))
              L <- length(group)
              summary.matrix[, 3] <- group
              summary.matrix[, 2] <- c(min(sorted.R), group[-L] + 1)
              member <- apply(1 - outer(R[,cd], group, "<="),1,sum) + 1
              
              for (g in 1:L){
                summary.matrix[g, 5] <- mean(Xaggclus[member == g, cc])
                summary.matrix[g, 6] <- mean(Xaggclus[member == g, cd])
                summary.matrix[g, 7] <- sd(Xaggclus[member == g, cc])
                summary.matrix[g, 8] <- sd(Xaggclus[member == g, cd])
              }# end g-loop
              resultsC[[k]][[2]] <- summary.matrix
              
              ac <- violation.matrix[1:(rvm - 1), 1] <- L - 1
              E.1 <- summary.matrix[, 5]
              E.2 <- summary.matrix[, 6]
              d <- E.2 - E.1
              d[d <= minvi] <- 0
              vi <- length(d[d > minvi/2])
              sum.vi <- sum(d)
              t.statistic <- rep(0, L)
              t.pvalue <- rep(1, L)
              if (any(d > 0)) {
                for (gg in 1:L) {
                  if (d[gg] > 0) {
                    tt <- t.test(Xaggclus[member==gg,cc],Xaggclus[member==gg,cd],alternative="less", paired = TRUE)
                    t.statistic[gg] <- abs(tt$statistic)
                    t.pvalue[gg] <- tt$p.value
                  }#endif
                }#end gg-loop
              }#endif
              violation.matrix[1, 2:8] <- c(vi,vi/ac, max(d), sum.vi, sum.vi/ac, max(t.statistic),sum(t.pvalue < alpha))
              vi.matrixC[cc,cd] <- vi.matrixC[cc,cd] + sum(t.pvalue < alpha)
              resultsC[[k]][[3]] <- violation.matrix
            }# end d-loop
          }# end c-loop
        }# end if(clusters)
    } # end if(method="MIIO") 
    
    vi.matrix <- vi.matrix + t(vi.matrix)
    vi.matrix.a <- vi.matrix
    VI <- matrix(apply(sign(vi.matrix), 1, sum))
    items.removed <- NULL
    if (item.selection) {
      repeat {
        nvi <- apply(sign(vi.matrix.a), 2, sum)
        if (sum(nvi) == 0) 
          break
        maxvi <- which(nvi == max(nvi))
        if (length(maxvi) > 1) {
          H <- rep(0, length(maxvi))
          for (i in 1:length(maxvi)) H[i] <- coefHTiny(X[, 
                                                         -c(items.removed, maxvi[i])])$H + rnorm(1, 
                                                                                                 0, 1e-05)
          maxvi <- maxvi[which(H == max(H))[1]]
        } # end if
        vi.matrix.a[maxvi, ] <- vi.matrix.a[, maxvi] <- 0
        items.removed <- c(items.removed, maxvi)
        VI. <- matrix(apply(sign(vi.matrix.a), 1, sum))
        VI.[items.removed, ] <- NA
        VI <- cbind(VI, VI.)
      }
    } # end if(item.selection)
    dimnames(VI) <- list(I.labels, paste("step", 1:ncol(VI)))
    if (verbose) 
      print(VI)
    if (length(items.removed) > 0) {
      X <- X[, -items.removed]
      if(clusters) {
        item.cluster.order <- item.cluster.order[-items.removed]
      } 
    } 
    
    if(clusters) {
      vi.matrixC <- vi.matrixC + t(vi.matrixC)
      vi.matrixC.a <- vi.matrixC
      VIC <- matrix(apply(sign(vi.matrixC),1,sum)) 
      dimnames(VIC) <- list(C.labels, paste("step", 1:ncol(VIC)))
      VIclus <- sapply(1:nC, function(x) abs(unique(item.cluster.order)[x] - unique(item.cluster.order)[which(vi.matrixC[x, ] > 0)]))
      VIbetw <- sapply(1:J, function(x) abs(item.cluster.order[x] - item.cluster.order[which(vi.matrix[x, ] > 0)])) # mss hier vi.matrixC of VIC vergelijken met item.cluster.order, van unique > 1, 
      # dan is het between schending.
      if(is.list(VIclus)) {
        MaxDcWICOvi <- matrix(unlist(lapply(VIclus, function(x) ifelse(length(x) > 0, max(x), 0))))
      } else {
        MaxDcWICOvi <- matrix(apply(VIclus, 2, function(x) ifelse(length(x) > 0, max(x), 0)))
      }
      
      
      
      VIC <- cbind(VIC, 
                   MaxDcWICOvi,
                   #tapply(matrix(unlist(lapply(VIbetw, function(x) length(unique(x))))), item.cluster.order.o, sum), # nr of clusters with which the items violate strong ICO, maar zit wel overlap in, bv als item i1 en i2 allebei schenden met item in cluster j, dan is dit aantal 2. Maar als item i1 overlapt met j1 en j2 dan is het aantal 1. Dat is eigenlijk raar! Dus, informatiever om sicovicluster-level: te tellen er any item in dit cluster schendt met any item in ander cluster, en dat sommeren over clusters. EN sicovi-itemlevel: totaal aantal sico schendingen over items heen. Voor eerste is maximale aantal C-1 en voor tweede het maximale aantal J-Jc. Dat hieronder dus gedaan.
                   tapply(matrix(unlist(lapply(VIbetw, function(x) length((x))))), item.cluster.order.o, sum), # nu totaal opgesoms dus echt aantal item-paar schendingen van strong ico.
                   tapply(matrix(unlist(lapply(VIbetw, function(x) length((x[x > 0]))))), item.cluster.order.o, sum), # nu totaal opgesoms dus echt aantal item-paar schendingen van strong ico.
                   unlist(lapply(tapply(lapply(VIbetw, function(x) (unique(x))), item.cluster.order.o, function(y) unlist(y)), function(z) length(unique(z)))),
                   tapply(unlist(lapply(VIbetw, function(x) ifelse(length(x) > 0, max(x), 0))), item.cluster.order.o, max)) # maximale afstand van schending strong ico.
      dimnames(VIC) <- list(C.labels, c("WICOvi", 
                                        "MaxDcWICOvi", 
                                        #"SICOvi",
                                        "IIOtotalvi",
                                        "SICOviItems",
                                        "SICOviClus",
                                        "MaxDcSICOvi"))
      
      # Total number of clusters with which this item violates
      # Maximum cluster distance:
      VIbetw <- cbind("NrBetwClusVio" = matrix(unlist(lapply(VIbetw, function(x) sum((x>0))))), 
                      "NrDiffClusVio" = matrix(unlist(lapply(VIbetw, function(x) length(unique(x[x > 0]))))), 
                      "MaxClusVio" = matrix(unlist(lapply(VIbetw, function(x) ifelse(length(x) > 0, max(x), 0)))))
      VI <- cbind(VI, VIbetw)  
      #dimnames(VIbetw) <- list(I.labels, c("NrBetwClusVio", "NrDiffClusVio", "MaxClusVio"))
      colnames(VI) <- c("NrIIOvio", "NrBetwClusVio", "NrDiffClusVio", "MaxClusVio")
    }
    
    # Computation of HT
    HT <- coefHT(X)
    iio.list <- list(results = results, violations = VI, items.removed = items.removed, 
                     Hi = Hi, HT = HT, method = method, item.mean = item.mean, 
                     m = ncat - 1)
    class(iio.list) <- "iio.class"
    
    if(clusters) {
      Xclust <- X[apply(X, 1, sd) > 0, ] 
      HTB <- coefHTB(Xclust, item.cluster.order)
      
      HT <- HTB
      
      iio.listC <- list(results = resultsC, violations = VIC, items.removed = items.removed, Hi = HiC, HT = HT, method = method, item.mean = cluster.mean, m = ncat-1)
      class(iio.list) <- class(iio.listC) <- "iio.class"
      iio.list <- list("Items" = iio.list, "Clusters" = iio.listC)
      class(iio.list) <- "iio.class"
    } else {
      if (!is.null(level.two.var)) {
        HTA <- coefHT(Xas)
        iio.lista <- list(results = resultsa, violations = VIA, 
                          items.removed = items.removed, Hi = His[, 2], HT = HTA, 
                          method = method, item.mean = item.mean, m = ncat - 
                            1)
        class(iio.lista) <- "iio.class"
        iio.list <- list(Level.1 = iio.list, Level.2 = iio.lista)
        class(iio.list) <- "iio.class"
      }
    }
    
    return(iio.list)
  }


coefHT <- function(Y, largeN = 1000){
   eq.var <- apply(Y, 1, sd)     # item-score variance per respondent 
   Y <- Y[eq.var > 0, ]           # remove respondents without variance, as they do not contribute to the computation of HT 
   N <- nrow(Y)                   # readjust the number of respondents in the analysis (N) 
   tY <- t(Y)                     # take the transpose of Y: used to compute the numerator of HT
   tYm <- apply(tY, 2, sort)      # take the transpose of Y, with all columns sorted: used to compute the denominator of HT

   G <- ceiling(N / largeN)       # Given "largeN:, G is the number of batches in which HT is computed 
                                  # (G - 1) batches of largeN (1000) respondents, and 
                                  # 1 batch of "N - (G-1) * largeN" respondents; i.e., the remainder 
   if (largeN %% N == 0) largeN <- largeN - 1
                                  # A correction for the case that N == largeN (which results in an error "noted by Stefanie Wind, sep 2023" 
   repeat{                        # A correction for the case that "the remainder" consists of one respondent only.  
     if (N %% largeN != 1) break
     largeN <- largeN -1
   }

   sum.S <- 0
   sum.Smax <- 0

   for (i in 1 : G) for (j in 1 : G){
     ni <- ifelse(i < G, largeN, N %% largeN)
     nj <- ifelse(j < G, largeN, N %% largeN)
     S <- var(tY[, (i-1) * largeN + 1 : ni], tY[, (j-1) * largeN + 1 : nj])
     Smax <- var(tYm[, (i-1) * largeN + 1 : ni], tYm[, (j-1) * largeN + 1 : nj])
     if (i == j) diag(S) <- diag(Smax) <- 0
     sum.S <- sum.S + sum(S)
     sum.Smax <- sum.Smax + sum(Smax)
   }
   return(sum.S / sum.Smax)
}

"coefHTB" <- function(X, level.two.var = NULL) {
  # This function can result in biased estimates for unequal cluster sizes (MLcoefH(.) can be used but is very slow.)
  eq.var <- apply(X, 1, sd)
  X <- X[eq.var > 0, ]
  X <- t(X)
  J <- ncol(X)
  C <- length(unique(level.two.var))
  S <- Smax <- 0
  X <- data.frame(X)
  
  r <- table(level.two.var)[1]
  
  rowsi <- rep(1:r, r)[-c(1:r * r - (r-1):0)]
  rowsi <- rep(r * (0:(C - 1)), each = length(rowsi)) + rowsi
  rowsj <- rep(1:r, each = r)[-c(1:r * r - (r-1):0)]
  rowsj <- rep(r * (0:(C - 1)), each = length(rowsj)) + rowsj
  
  Xjall <- X[rowsj, ]
  Xjsort <- apply(Xjall, 2, sort)
  for(i in 1:(J - 1)) {
    Xi <- X[rowsi, i]
    for(j in (i + 1):J) {
      Xj <- Xjall[, j]
      S <- S + var(Xi, Xj)
      Smax <- Smax + var(Xjsort[, c(i, j)])[1,2]
    }
  }
  
  HB <- S/Smax 
  
  S <- var(X)
  Smax <- var(apply(X, 2, sort))
  diag(S) <- 0
  diag(Smax) <- 0
  HW <- sum(S)/sum(Smax)
  
  return(cbind("HT" = HW, "HTB" = HB, "HTB/HT" = HB/HW))
}

# DECREPIT New function above for package version 3.1.1

#coefHT <- function(Y, largeN = 1000){ 
#   eq.var <- apply(Y, 1, sd)
#   Y <- Y[eq.var > 0, ]
#   N <- nrow(Y)
#   tY <- t(Y)
#   tYm <- apply(tY, 2, sort)
#
#   G <- ceiling(N / largeN)
#   sum.S <- 0
#   sum.Smax <- 0
#   repeat{
#     if (N %% largeN != 1) break
#     largeN <- largeN -1
#   }
#   for (i in 1 : G) for (j in 1 : G){
#     ni <- ifelse(i < G, largeN, N %% largeN)
#     nj <- ifelse(j < G, largeN, N %% largeN)
#     S <- var(tY[, (i-1) * largeN + 1 : ni], tY[, (j-1) * largeN + 1 : nj])
#     Smax <- var(tYm[, (i-1) * largeN + 1 : ni], tYm[, (j-1) * largeN + 1 : nj])
#     if (i == j) diag(S) <- diag(Smax) <- 0
#     sum.S <- sum.S + sum(S)
#     sum.Smax <- sum.Smax + sum(Smax)
#   }
#   return(sum.S / sum.Smax)
#}


# DECREPIT New function above 
#coefHT <- function(Y, largeN = 1000){
#    eij <- function(x) {
#        tab <- tabulate(x[ ,1] + 2 * x[ ,2] + 1L, 4)
#        e <- min(tab[3], tab[2])
#        e0 <- ((tab[1] + e) * (e + tab[4])) / nrow(x)
#        return(c(e, e0))
#    }
#
#   eq.var <- apply(Y, 1, sd)
#   Y <- Y[eq.var > 0, ]
#   N <- nrow(Y)
#   tY <- t(Y)
#   if (N < largeN){ 
#      S    <- var(tY)
#      Smax <- var(apply(tY, 2, sort))
#      diag(S) <- diag(Smax) <- 0
#      HT <- sum(S)/sum(Smax)
#   } else {  
#      e1 <- 0L
#      e0 <- 0L
#      for (i in 1 : (ncol(tY) - 1)){
#         for (j in (i + 1) : ncol(tY)){
#            e <- eij(tY[, c(i,j)])
#            e1 <- e1 + e[1]
#            e0 <- e0 + e[2]
#         }
#         cat("\r", "Large sample size, computation may take several minutes. Progress: ", round((i * 100) / (ncol(tY) - 1), 1), "%      ", sep = "") 
#         flush.console()  # Make visible on RGui
#      }
#      cat("\r", "Large sample size, computation may take several minutes. Progress: DONE    ", sep = "") 
#      cat("\n")
#      HT <- 1 - e1/e0
#   }
#   return(HT)
#}
