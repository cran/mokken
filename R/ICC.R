######## Letty Koopman 
######## University of Amsterdam
# ICC() version 16-03-2020
#
# Computation of the intraclass correlation per item and for the total scale
# includes F-test on the null hypothesis that the total scale ICC is 0.

"ICC" <- function(X) { 
  # Computes the intraclass correlation per item and for the total scale of a questionnaire
  #
  # Args:
  #   X: Data matrix with one subject column (the first) and at least one column with item or test scores. 
  # 
  # Returns: 
  #   ICC per item and for the total scale.
  #   F-test for the total scale ICC
  #
  # Error handling:
  #
  if (!is.matrix(X) && !is.data.frame(X)) stop("Data are not matrix or data.frame")
  if (ncol(X) < 2) stop("Data must consist of one groupingscolumn and at least one column with item or test scores")
  if (any(is.na(X))) stop("Missing values are not allowed")
  X <- X[order(X[, 1]), ] 
  Rs <- as.numeric(table(X[, 1]))
  X[, 1] <- rep(1:length(Rs), Rs)
  X <- as.matrix(X)
  if (mode(X)!="numeric") stop("Data must be numeric")
  if(is.null(colnames(X))) colnames(X) <- c("Subs", paste("Item", 1:(ncol(X) - 1)))
  if(any(table(colnames(X)) > 1)){
    warning("At least two items have the same name. The duplicates have received an index.")
    colnames(X) <- make.unique(colnames(X))
  }
  if (any(apply(X, 2, var) < 1e-40)) 
    stop("One or more variables have zero variance.")
  # Ensure each subject has > 1 rater
  if(any(Rs == 1)){ 
    warning('For at least one group there is only 1 respondent. The ICC is computed without this (these) groups(s).') 
    X <- X[!(X[, 1] %in% which(Rs == 1)), ]
    Rs <- as.numeric(table(X[, 1]))
    X[, 1] <- rep(1:length(Rs), Rs)
  }
  subjects <- X[, 1] 
  S <- length(unique(subjects))
  R <- nrow(X)
  hmean <- S / sum( 1 / table(subjects) ) # harmonic mean
  
  iscores <- X[, -1] # item scores
  J <- ncol(X) - 1
  nwiggle <- mean(Rs) - var(Rs) / (S * mean(Rs))
  
  if(J > 1) { 
    if(is.null(colnames(iscores))) colnames(iscores) <- paste("Item", 1:(ncol(X) - 1))
    # Item variances
    ismean <- isvar <- matrix(NA, nrow = S, ncol = J)
    for(i in 1:J) {
      ismean[, i] <- tapply(iscores[, i], subjects, mean)
      isvar[, i] <- tapply(iscores[, i], subjects, var)
    }
    
    iSW2 <- apply(isvar, 2, function(i) 1 / (R - S) * sum((Rs - 1) * i))
    iSB2 <- sapply(1:J, function(i) 1 / (nwiggle * (S - 1)) * sum(Rs * (ismean[, i] - mean(iscores[, i]))^2))
    iTV <- (R - S) / (R - 1) * iSW2 + nwiggle * (S - 1) / (R - 1) * iSB2 
    itau2 <- iSB2 - iSW2 / nwiggle
    isigma2 <- iSW2
    iICC <- data.frame("Item" = colnames(iscores), "ICC" = round(itau2 / (itau2 + isigma2), 3))
    
   
    # Total scores:
    tscores <- rowMeans(iscores)
  } else { 
    iICC <- "Only one column with scores available"
    tscores <- X[, 2]
  }
  # Total scale variances
  smean <- tapply(tscores, subjects, mean)
  svar <- tapply(tscores, subjects, var)
  SW2 <- 1 / (R - S) * sum((Rs - 1) * svar) # within-subject variance
  SB2 <- 1 / (nwiggle * (S - 1)) * sum(Rs * (smean - mean(tscores))^2) # between-object variance 
  TV <- (R - S) / (R - 1) * SW2 + nwiggle * (S - 1) / (R - 1) * SB2 
  tau2 <- SB2 - SW2 / nwiggle
  sigma2 <- SW2
  ICCt <- tau2 / (tau2 + sigma2) # Intraclass correlation
  # Results
  Ftest <- nwiggle * SB2 / SW2
  pval <- pf(Ftest, S - 1, R - S, lower.tail = FALSE)
  ICCF <- data.frame(ICC = round(ICCt, 3), "F" = round(Ftest, 3), "df1" = S - 1, "df2" = R - S, "p.value" = pval)
  rownames(ICCF) <- ""
  results <- list("itemICC" = iICC, "scaleICC" = ICCF)
  return(results)
}

