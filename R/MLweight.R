## MLweight function 18-08-2017. Letty Koopman. 
# Last adjusted 16-04-2018 improve orderings

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
  if(is.null(minx)) {
				minx <- min(X[, -1])
				warning(paste('minx has not been specified. Minimum X value has been determined from the items as value', minx, sep = " "))
			   }
  if(is.null(maxx)) {
				maxx <- max(X[, -1])
				warning(paste('maxx has not been specified. Maximum X value has been determined from the items as value', maxx, sep = " "))
			   }
  X[, -1] <- X[, -1] - minx
  maxx <- maxx - minx
  minx <- 0
  s <- X[, 1] # Subject column
  f1 <- factor(X[, 2], levels = minx:maxx) # First item 
  f2 <- factor(X[, 3], levels = minx:maxx) # Second item
  
  # Compute relative frequencies
  Rel1 <- colSums(table(s, f1) / rowSums(table(s, f1))) 
  names(Rel1) <- paste0(1, minx:maxx) # for item i
  Rel2 <- colSums(table(s, f2) / rowSums(table(s, f2)))
  names(Rel2) <- paste0(2, minx:maxx) # for item j
  
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
          if(any(h >= 1 & h <= maxx) & any(h >= maxx + 1 & h <= maxx * 2)){ # zitten de itemstappen binnen hetzelfde item?
            select[j] <- (all(h[which(h >= 1 & h <= maxx)] == sort(h[which(h >= 1 & h <= maxx)])) & all(h[which(h >= maxx + 1 & h <= maxx * 2)] == sort(h[which(h >= maxx + 1 & h <= maxx * 2)]))) * j
          } else {
            i1 <- ifelse(any(h >= 1 & h <= maxx), all(h[which(h >= 1 & h <= maxx)] == sort(h[which(h >= 1 & h <= maxx)])), 0)
            i2 <- ifelse(any(h >= maxx + 1 & h <= maxx * 2), 
                         all(h[which(h >= maxx + 1 & h <= maxx * 2)] == sort(h[which(h >= maxx + 1 & h <= maxx * 2)])), 0)
            select[j] <- (i1 + i2) * j
          }
        }
        o[[i]] <- matrix(apply(o[[i]][select, ], 1, paste0, collapse = "."))#, ncol = ncol(o[[i]]))
      }
    }
    
    out <- matrix(apply(expand.grid(o), 1, paste0, collapse = "."))
    out <- matrix(unlist(strsplit(out, "[.]")), nrow = nrow(out), byrow = T)
    
    w <- NULL
    Z <- matrix(rep(matrix(all.patterns(2, maxx + 1), nrow = 1), maxx), nrow = maxx, byrow = TRUE)
    Z <- matrix(ifelse(Z < row(Z), 0, 1), ncol = (maxx) * 2, byrow = TRUE)
    for(i in 1:nrow(out)){
      ords <- as.numeric(out[i, ])#as.numeric(unlist(strsplit(out[i, ],NULL)))
      # Compute Z matrix for each item-response pattern
      Z1 <- Z[, (ords)]
      # Compute weights
      w <- rbind(w, apply(Z1, 1, function(x){sum(x * cumsum(abs(x - 1)))}))#
    }
    wr <- matrix(colMeans(w), nrow = 1)
  } else {
    #out <- paste(names(y), collapse = "")
    
    ords <- as.numeric(names(y))#as.numeric(unlist(strsplit(out,NULL)))
    # Compute Z matrix for each item-response pattern
    Z <- matrix(rep(matrix(all.patterns(2, maxx + 1), nrow = 1), maxx), nrow = maxx, byrow = TRUE)
    Z <- matrix(ifelse(Z < row(Z), 0, 1), ncol = (maxx) * 2, byrow = TRUE)
    Z <- Z[, (ords)]
    
    # Compute weights
    wr <- apply(Z, 1, function(x){sum(x * cumsum(abs(x - 1)))})
    
    
  }
  return(wr)
}

