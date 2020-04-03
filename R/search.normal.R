# Aangepast op 9 maart 2020 door Letty Koopman
# - type.se argument: "Z" = Z-score zoals origineel geprogrammeerd, "d1" = delta method one-level, of "d2" delta method two-level
# - test.Hi argument: If FALSE: Test if Hi is significantly larger than zero. If TRUE: Test if Hi is significantly larger than lowerbound c 

# - coefZ kan nu Zij, Zi, Z-score berekenen met null hypothese Hij, Hi, H = lowerbound, met originele berekening (Molenaar & Sijtsma 2000), 
#         delta methode met eenlevel assumptie (Kuijpers et al., 2013) of delta methode met tweelevel assumptie (Koopman et al., 2019)



# Aangepast op 21 oktober 2015, 29 juni 2016
"search.normal" <-
function(X, lowerbound, alpha, StartSet, verbose, type.se, test.Hi, level.two.var){
  #type.se should be "Z" or "delta"

   # Internal functions
   
   any.neg <- function(x){if(any(x < 0)) TRUE else FALSE}

   adjusted.alpha <- function(alpha, K) alpha/(K[1]*(K[1]-1)*.5 + sum(K[-1]))

   fitstring <- function(string.arg,length.arg) substr(paste(string.arg,"                        "),1,length.arg)
   
   newH <- function(j,in.this.set, x, lb, Z.c, type.se, test.Hi, level.two.var){ # can I use the same arguments as in function?

     newX <- cbind(x[, in.this.set == 1], x[, j])
     H.list <- coefH(newX, FALSE)
     if (H.list$Hi[length(H.list$Hi)] < lb) return(-98) # less than lower bound
     Zi <- coefZ(newX, lowerbound = test.Hi * lb, type.se = type.se, level.two.var = level.two.var)$Zi # Test if Hi is significantly larger than zero (test.Hi == FALSE) or lowerbound c (test.Hi == TRUE)
     if (Zi[length(Zi)] < Z.c) return(-97)                      # not significant
     return(H.list$H)
   }
   
   
   if(!is.null(level.two.var)) {
     if (nrow(as.matrix(level.two.var)) != nrow(X)) {
       level.two.var <- NULL
       warning("level.two.var not the same length/nrow as X: level.two.var is ignored.")
     } else if(type.se == "Z") {
       level.two.var <- NULL
       warning("level.two.var is ignored for type.se = 'Z'.")
     } else if(any(is.na(level.two.var))) {
       level.two.var <- NULL
       warning("level.two.var contains missing value(s): level.two.var is ignored.")
     } else {
       X <- X[order(level.two.var), ]
       level.two.var <- sort(level.two.var)
       Rs <- as.numeric(table(level.two.var))
       level.two.var <- rep(1:length(Rs), Rs)
       # Ensure each subject has > 1 rater
       if(any(Rs == 1)){ 
         warning('For at least one group there is only 1 respondent. The aisp is performed without this (these) group(s).') 
         cases <- !(level.two.var %in% which(Rs == 1))
         X <- X[cases, ]
         level.two.var <- level.two.var[cases]
       }
     }
   }
   # initial calculations
   item.label <- dimnames(X)[[2]]
   N <- nrow(X)
   S <- var(X)

   if(any(is.na(StartSet)) | any(is.nan(StartSet)) | (class(StartSet) == "logical" & any(StartSet == TRUE))){
      warning("Start set of items is not properly defined. User-defined StartSet ignored.")
      StartSet <- FALSE
   }   
   if(!(class(StartSet) == "logical" | class(StartSet) == "numeric"| class(StartSet) == "integer")){
      warning("class(StartSet) should be logical or numeric. User-defined StartSet ignored.")
      StartSet <- FALSE
   }   
   StartSet <- unique(StartSet)
   if (class(StartSet) != "logical" & all(StartSet %in% 1 : ncol(X)) == FALSE){
      warning("Start set of items is not properly defined. User-defined StartSet ignored")
      StartSet <- FALSE
   }   
   startset.provided <- FALSE
   if (class(StartSet) != "logical"){
     startset.provided <- TRUE
     cat("Items ", StartSet, " in start set", fill = TRUE) 
   }  
   if(any(is.na(diag(S/S)))) stop("At least one item has no variance")

   Smax <- var(apply(X,2,sort))
   Hij <- S/Smax
   Zij <- coefZ(X, lowerbound = 0, type.se = type.se, level.two.var = level.two.var)$Zij # Test if Hij is significantly larger than zero

   J <- nrow(Hij)
   output <- NULL

   for (lb in lowerbound){
      result <- rep(-99, J);
      j <- 0
      InSet <- rep(0, J)
     
      scale <- 0

      # Start Scaling
      repeat{
        scale <- scale + 1
        step <- 1
        K <- rep(0, J)
      
        if(verbose){ 
          cat("", fill = TRUE)
          cat("SCALE", scale, fill = TRUE)
          cat("", fill = TRUE)
        }  
   
        # Are there at least two items left?
        if(length(InSet[InSet == 0]) < 2){
          if(verbose) cat("Less than two items left. PROCEDURE STOPS", fill = TRUE)
          break
        }
   
        # Compute the critical value for Zij 
        
        K[step] <- length(InSet[InSet == 0])
        Z.c <- abs(qnorm(adjusted.alpha(alpha,K)))
   
        # Determine which items can be selected into the same scale
   
        Hselect <- Hij 
        Hselect[abs(Zij) < Z.c] <- -99
        Hselect[InSet > 0 & InSet < scale,] <- -99
        Hselect[,InSet > 0 & InSet < scale] <- -99
        Hselect[col(Hselect) >= row(Hselect)] <- -99
        eps <- row(Hselect) * 1e-10
        Hselect[Hselect != 99] <- Hselect[Hselect != 99] - eps[Hselect != 99]
        
        # Check if there are any feasible values left
        
        if(max(round(Hselect)) == -99){
          if(verbose) cat("Scale ", scale," could not be formed because all Hj < ",lb," or because no Hij significantly greater than zero. PROCEDURE STOPS", fill = TRUE)
          break
        }
   
        if (class(StartSet) == "logical" | scale > 1){
          first.item <- row(Hselect)[Hselect==max(Hselect)]
          second.item <- col(Hselect)[Hselect==max(Hselect)]
          StartSet <- c(first.item, second.item)
        }  
        if (length(StartSet) == 1 & scale == 1){ 
          first.item <- StartSet 
          max.tmp <- max(c(Hselect[first.item,],Hselect[,first.item]))
          second.item <- which(abs(Hij[first.item,] - max.tmp) < 1e-6)
          StartSet <- c(first.item, second.item)
        }
   
        StartSet <- sort(StartSet)
        # Check whether items meet Mokken's criteria
   
        # (1) All Hij significantly greater than 0 (Check is necessary only if startset has been provided)
        if (startset.provided & scale == 1){
            StartHij <- Hselect[StartSet, StartSet]
            if(any(StartHij[row(StartHij) > col(StartHij)] < 0)) warning("Items in start set do not form a Mokken scale: Some Hij are not significantly greater than zero")
        }
   
        # (2) All Hi greater than lower bound (Check is always necessary)
        if(test.Hi) {
          StartHi <- coefZ(X[, StartSet], lowerbound = lb, type.se = type.se, level.two.var = level.two.var)$Zi
          checkHi <- min(abs(StartHi)) < Z.c
        } else {
          StartHi <- coefH(X[, StartSet], FALSE)[[2]]
          checkHi <- min(StartHi) < lb
        } 
        if(checkHi){
          if(startset.provided & scale == 1) warning("Items in start set do not form a Mokken scale: Some Hj < ",lb) 
          if(!startset.provided | scale > 1) {  
            if(verbose) cat("Scale ", scale," could not be formed due to H < ",lb,". PROCEDURE STOPS", fill = TRUE) 
            break
          }
        }
   
        # Add the StartSet to the scale 
        if(verbose) for (i in 1 : length(StartSet)){
          the.item <- ifelse(StartSet[i] > 9, as.character(StartSet[i]), paste(" ", as.character(StartSet[i]), sep = ""))
          cat("Item ", the.item, ": ", fitstring(item.label[StartSet[i]],20), " Scale", scale," H = ",round(coefH(X[, StartSet], FALSE)[[3]], 2), ifelse(scale == 1 & startset.provided, "StartSet", ""), fill = TRUE)
        } 
        InSet[StartSet] <- scale
   
        # Adding new items
        repeat{
          step <- step + 1
   
          # exclude items from previous scales
          in.this.set <- InSet
          in.this.set <- ifelse(InSet == scale, 1, 0)
          in.this.set <- ifelse(InSet <  scale & InSet > 0, -1, in.this.set)
   
          # exclude items having a negative covariance with the already selected items
          neg1 <- apply(Hij[in.this.set == 1, ], 2, any.neg)
          neg2 <- apply(Hij[, in.this.set == 1], 1, any.neg)
          in.this.set[(neg1 | neg2) & in.this.set == 0] <- -1
   
          # Are there items left after the exclusion?
          available.items <- which(in.this.set == 0)
          if(length(available.items) == 0){
            if(verbose) cat("Scale ", scale," is completed. No items left with Hij => 0",fill = TRUE)
            break
          }
   
         # Compute H and Hi of potentially new items
         result[in.this.set != 0] <- -99  # items already selected in other scales
         K[step] <- length(available.items)
         Z.c <- abs(qnorm(adjusted.alpha(alpha, K)))
         for (j in available.items) result[j] <- newH(j,in.this.set, X, lb, Z.c, type.se, test.Hi, level.two.var)
   
   
         # Is maximum value Hi greater than c?
         if(max(result) < lb){
           if(verbose) cat("Scale ", scale," is completed. No items left such that Hi > ",lb,".", fill = TRUE)
           break
         }
   
         # Add the newly selected item to the scale
         new.item <- row(as.matrix(result))[result == max(result)]
         InSet[new.item] <- scale
         if(verbose){
           the.item <- ifelse(new.item > 9, as.character(new.item), paste(" ", as.character(new.item), sep = ""))
           cat("Item ", the.item, ": ", fitstring(item.label[new.item],20), " Scale", scale, " H = ", round(max(result), 2), fill = TRUE)
         } 
      
      
   
       }
     # BEGIN SEARCH EXTENDED
     #if (search.extended == TRUE) {
        # UNDER CONSTRUCTION
        # in.this.extended.set = in.this.set
        # in.this.extended.set[in.this.extended.set == -1] = 0
        #  repeat{
        #    step <- step + 1
        # 
        #    # exclude items from previous scales
        #  in.this.set <- in.this.extended ????
        #  in.this.set <- ifelse(InSet == scale, 1,0) ????
        #  in.this.set <- ifelse(InSet <  scale & InSet > 0,-1,in.this.set) ????
        # 
        #  # exclude items having a negative covariance with the already selected items
        #  neg1 <- apply(Hij[in.this.set==1,],2,any.neg)
        #  neg2 <- apply(Hij[,in.this.set==1],1,any.neg)
        #  in.this.set[neg1|neg2 & in.this.set==0] <- -1
        # 
        #  # Are there items left after the exclusion?
        #  available.items <- which(in.this.set==0)
        #  if(length(available.items)==0){
        #    if(verbose) cat("Scale ", scale," is completed. No items left with Hij => 0", fill = TRUE)
        #    break
        #  }
        #
        #  # Compute H and Hi of potentially new items
        #  result[in.this.set!=0] <- -99  # items already selected in other scales
        #  K[step] <- length(available.items)
        #  Z.c <- abs(qnorm(adjusted.alpha(alpha,K)))
        #  for (j in available.items) result[j] <- newH(j,in.this.set, X, lb, Z.c)
        #
        #
        #  # Is maximum value Hi greater than c?
        #  if(max(result) < lb){
        #    if(verbose) cat("Scale ", scale," is completed. No items left such that Hi > ",lb,".", fill = TRUE)
        #    break
        #  }
        #
        #  # Add the newly selected item to the scale
        #  new.item <- row(as.matrix(result))[result==max(result)]
        #  InSet[new.item] <- scale
        #  if(verbose) cat("Item: ",fitstring(item.label[new.item],20)," Scale", scale," H = ",round(max(result),2), fill = TRUE)
        #}
        # 
     #}  
        
     # EINDE SEARCH EXTENDED
     # start with next scale
        
     }  
     output <- cbind(output, as.matrix(InSet))
   }
   dimnames(output) <- list(item.label, lowerbound)
   return(output)
}
   
   

"search.normal.old" <-
  function(X, lowerbound, alpha, StartSet, verbose){
    
    # Internal functions
    
    any.neg <- function(x){if(any(x < 0)) TRUE else FALSE}
    
    adjusted.alpha <- function(alpha, K) alpha/(K[1]*(K[1]-1)*.5 + sum(K[-1]))
    
    fitstring <- function(string.arg,length.arg) substr(paste(string.arg,"                        "),1,length.arg)
    
    newH <- function(j,in.this.set, x, lb, Z.c){
      
      newX <- cbind(x[, in.this.set == 1], x[, j])
      H.list <- coefH(newX, FALSE)
      if (H.list$Hi[length(H.list$Hi)] < lb) return(-98) # less than lower bound
      Zi <- coefZ(newX)$Zi
      if (Zi[length(Zi)] < Z.c) return(-97)                      # not significant
      return(H.list$H)
    }
    
    # initial calculations
    
    item.label <- dimnames(X)[[2]]
    N <- nrow(X)
    S <- var(X)
    
    if(any(is.na(StartSet)) | any(is.nan(StartSet)) | (class(StartSet) == "logical" & any(StartSet == TRUE))){
      warning("Start set of items is not properly defined. User-defined StartSet ignored.")
      StartSet <- FALSE
    }   
    if(!(class(StartSet) == "logical" | class(StartSet) == "numeric"| class(StartSet) == "integer")){
      warning("class(StartSet) should be logical or numeric. User-defined StartSet ignored.")
      StartSet <- FALSE
    }   
    StartSet <- unique(StartSet)
    if (class(StartSet) != "logical" & all(StartSet %in% 1 : ncol(X)) == FALSE){
      warning("Start set of items is not properly defined. User-defined StartSet ignored")
      StartSet <- FALSE
    }   
    startset.provided <- FALSE
    if (class(StartSet) != "logical"){
      startset.provided <- TRUE
      cat("Items ", StartSet, " in start set", fill = TRUE) 
    }  
    if(any(is.na(diag(S/S)))) stop("At least one item has no variance")
    
    Smax <- var(apply(X,2,sort))
    Hij <- S/Smax
    Sij <- outer(apply(X,2,var),apply(X,2,var),"*")
    Zij <- (S * sqrt(N-1))/sqrt(Sij)
    
    J <- nrow(Hij)
    output <- NULL
    
    for (lb in lowerbound){
      result <- rep(-99, J);
      j <- 0
      InSet <- rep(0, J)
      
      scale <- 0
      
      # Start Scaling
      repeat{
        scale <- scale + 1
        step <- 1
        K <- rep(0, J)
        
        if(verbose){ 
          cat("", fill = TRUE)
          cat("SCALE", scale, fill = TRUE)
          cat("", fill = TRUE)
        }  
        
        # Are there at least two items left?
        if(length(InSet[InSet == 0]) < 2){
          if(verbose) cat("Less than two items left. PROCEDURE STOPS", fill = TRUE)
          break
        }
        
        # Compute the critical value for Zij 
        
        K[step] <- length(InSet[InSet == 0])
        Z.c <- abs(qnorm(adjusted.alpha(alpha,K)))
        
        # Determine which items can be selected into the same scale
        
        Hselect <- Hij 
        Hselect[abs(Zij) < Z.c] <- -99
        Hselect[InSet > 0 & InSet < scale,] <- -99
        Hselect[,InSet > 0 & InSet < scale] <- -99
        Hselect[col(Hselect) >= row(Hselect)] <- -99
        eps <- row(Hselect) * 1e-10
        Hselect[Hselect != 99] <- Hselect[Hselect != 99] - eps[Hselect != 99]
        
        # Check if there are any feasible values left
        
        if(max(round(Hselect)) == -99){
          if(verbose) cat("Scale ", scale," could not be formed because all Hj < ",lb," or because no Hij significantly greater than zero. PROCEDURE STOPS", fill = TRUE)
          break
        }
        
        if (class(StartSet) == "logical" | scale > 1){
          first.item <- row(Hselect)[Hselect==max(Hselect)]
          second.item <- col(Hselect)[Hselect==max(Hselect)]
          StartSet <- c(first.item, second.item)
        }  
        if (length(StartSet) == 1 & scale == 1){ 
          first.item <- StartSet 
          max.tmp <- max(c(Hselect[first.item,],Hselect[,first.item]))
          second.item <- which(abs(Hij[first.item,] - max.tmp) < 1e-6)
          StartSet <- c(first.item, second.item)
        }
        
        StartSet <- sort(StartSet)
        # Check whether items meet Mokken's criteria
        
        # (1) All Hij significantly greater than 0 (Check is necessary only if startset has been provided)
        if (startset.provided & scale == 1){
          StartHij <- Hselect[StartSet, StartSet]
          if(any(StartHij[row(StartHij) > col(StartHij)] < 0)) warning("Items in start set do not form a Mokken scale: Some Hij are not significantly greater than zero")
        }
        
        # (2) All Hi greater than lower bound (Check is always necessary)
        StartHi  <- coefH(X[, StartSet], FALSE)[[2]]
        if(min(StartHi) < lb){
          if(startset.provided & scale == 1) warning("Items in start set do not form a Mokken scale: Some Hj < ",lb) 
          if(!startset.provided | scale > 1) {  
            if(verbose) cat("Scale ", scale," could not be formed due to H < ",lb,". PROCEDURE STOPS", fill = TRUE) 
            break
          }
        }
        
        # Add the StartSet to the scale 
        if(verbose) for (i in 1 : length(StartSet)){
          the.item <- ifelse(StartSet[i] > 9, as.character(StartSet[i]), paste(" ", as.character(StartSet[i]), sep = ""))
          cat("Item ", the.item, ": ", fitstring(item.label[StartSet[i]],20), " Scale", scale," H = ",round(coefH(X[, StartSet], FALSE)[[3]], 2), ifelse(scale == 1 & startset.provided, "StartSet", ""), fill = TRUE)
        }  
        InSet[StartSet] <- scale
        
        # Adding new items
        repeat{
          step <- step + 1
          
          # exclude items from previous scales
          in.this.set <- InSet
          in.this.set <- ifelse(InSet == scale, 1, 0)
          in.this.set <- ifelse(InSet <  scale & InSet > 0, -1, in.this.set)
          
          # exclude items having a negative covariance with the already selected items
          neg1 <- apply(Hij[in.this.set == 1, ], 2, any.neg)
          neg2 <- apply(Hij[, in.this.set == 1], 1, any.neg)
          in.this.set[(neg1 | neg2) & in.this.set == 0] <- -1
          
          # Are there items left after the exclusion?
          available.items <- which(in.this.set == 0)
          if(length(available.items) == 0){
            if(verbose) cat("Scale ", scale," is completed. No items left with Hij => 0",fill = TRUE)
            break
          }
          
          # Compute H and Hi of potentially new items
          result[in.this.set != 0] <- -99  # items already selected in other scales
          K[step] <- length(available.items)
          Z.c <- abs(qnorm(adjusted.alpha(alpha, K)))
          for (j in available.items) result[j] <- newH(j,in.this.set, X, lb, Z.c)
          
          
          # Is maximum value Hi greater than c?
          if(max(result) < lb){
            if(verbose) cat("Scale ", scale," is completed. No items left such that Hi > ",lb,".", fill = TRUE)
            break
          }
          
          # Add the newly selected item to the scale
          new.item <- row(as.matrix(result))[result == max(result)]
          InSet[new.item] <- scale
          if(verbose){
            the.item <- ifelse(new.item > 9, as.character(new.item), paste(" ", as.character(new.item), sep = ""))
            cat("Item ", the.item, ": ", fitstring(item.label[new.item],20), " Scale", scale, " H = ", round(max(result), 2), fill = TRUE)
          } 
          
          
          
        }
        # BEGIN SEARCH EXTENDED
        #if (search.extended == TRUE) {
        # UNDER CONSTRUCTION
        # in.this.extended.set = in.this.set
        # in.this.extended.set[in.this.extended.set == -1] = 0
        #  repeat{
        #    step <- step + 1
        # 
        #    # exclude items from previous scales
        #  in.this.set <- in.this.extended ????
        #  in.this.set <- ifelse(InSet == scale, 1,0) ????
        #  in.this.set <- ifelse(InSet <  scale & InSet > 0,-1,in.this.set) ????
        # 
        #  # exclude items having a negative covariance with the already selected items
        #  neg1 <- apply(Hij[in.this.set==1,],2,any.neg)
        #  neg2 <- apply(Hij[,in.this.set==1],1,any.neg)
        #  in.this.set[neg1|neg2 & in.this.set==0] <- -1
        # 
        #  # Are there items left after the exclusion?
        #  available.items <- which(in.this.set==0)
        #  if(length(available.items)==0){
        #    if(verbose) cat("Scale ", scale," is completed. No items left with Hij => 0", fill = TRUE)
        #    break
        #  }
        #
        #  # Compute H and Hi of potentially new items
        #  result[in.this.set!=0] <- -99  # items already selected in other scales
        #  K[step] <- length(available.items)
        #  Z.c <- abs(qnorm(adjusted.alpha(alpha,K)))
        #  for (j in available.items) result[j] <- newH(j,in.this.set, X, lb, Z.c)
        #
        #
        #  # Is maximum value Hi greater than c?
        #  if(max(result) < lb){
        #    if(verbose) cat("Scale ", scale," is completed. No items left such that Hi > ",lb,".", fill = TRUE)
        #    break
        #  }
        #
        #  # Add the newly selected item to the scale
        #  new.item <- row(as.matrix(result))[result==max(result)]
        #  InSet[new.item] <- scale
        #  if(verbose) cat("Item: ",fitstring(item.label[new.item],20)," Scale", scale," H = ",round(max(result),2), fill = TRUE)
        #}
        # 
        #}  
        
        # EINDE SEARCH EXTENDED
        # start with next scale
        
      }  
      output <- cbind(output, as.matrix(InSet))
    }
    dimnames(output) <- list(item.label, lowerbound)
    return(output)
  }

