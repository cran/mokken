"aisp" <- function(X, 
                   lowerbound =.3, 
                   search = "normal", 
                   alpha = .05, 
                   StartSet = FALSE,
                   popsize = 20, 
                   maxgens = default.maxgens, 
                   pxover = 0.5, 
                   pmutation = 0.1,
                   verbose = FALSE,
                   type.z = "Z", 
                   test.Hi = FALSE, 
                   level.two.var = NULL)
# default.maxgens <- 10^(log2(ncol(X)/5)) * 1000
# lowerbound =.3; search = "normal"; alpha = .05; StartSet = FALSE; popsize = 20; maxgens = default.maxgens; pxover = 0.5; pmutation = 0.1; verbose = FALSE; type.z = "Z"; test.Hi = FALSE; level.two.var = NULL
{
   X <- check.data(X)
   output <- NULL
   clabels <- as.character(lowerbound)
   rlabels <- dimnames(X)[[2]]
   params <- c(alpha, pxover, pmutation)
   cparams <- c("alpha", "pxover", "pmutation")
   for (i in 1:3){
      if(!is.numeric(params[i])|is.na(params[i])) stop(cparams[i], " is not numeric")
      if (params[i] < 0) {warning(paste("Negative ",cparams[i],". ",cparams[i]," is set to 0")); assign(params[i],0)}
      if (params[i] > 1) {warning(paste(cparams[i]," greater than 1. ",cparams[i]," is set to 1")); assign(params[i],1)}
   }
   for (i in 1:length(lowerbound)) if(!is.numeric(lowerbound[i])|is.na(lowerbound[i])) stop(" lowerbound contains non-numeric values")
   tmp <- coefHTiny(X)$Hij; diag(tmp) <- 0; c.max <- max(tmp); rm(tmp)
   if (any(lowerbound > c.max)) {
      warning(paste("Some lower bounds are greater than max(Hij) rendering all items unscalable. Lower bounds greater than max(Hij) are removed"))
      lowerbound <- lowerbound[lowerbound <= c.max]
      if (length(lowerbound) == 0) stop("no lowerbound provided")
   } 
   default.maxgens <- 10^(log2(ncol(X)/5)) * 1000
   if(is.numeric(popsize)&!is.na(popsize)) popsize <- as.integer(popsize) else stop("popsize is not numeric")
   if(popsize < 1) stop("popsize is nonpositive")
   if(is.numeric(maxgens)&!is.na(maxgens)) maxgens <- as.integer(maxgens) else stop("maxgens is not numeric")
   if(maxgens < 1) stop("maxgens is nonpositive")
   if(search == "ga"){ 
      output <- NULL
      for (lb in lowerbound){ 
         ga <- search.ga(X, popsize, maxgens, alpha, lb, pxover, pmutation)
         output <- cbind(output, ga)
      }
      dimnames(output) <- list(dimnames(X)[[2]], lowerbound)
   } else if(search == "extended") {
      output <- search.extended(verbose) 
   } else {
      if(test.Hi & type.z == "Z") {
         type.z <- "WB"
         warning("type.z has been changed to 'WB' to enable testing Hi > c.")
      }
      if(!is.null(level.two.var)) {
         if (nrow(as.matrix(level.two.var)) != nrow(X)) {
            level.two.var <- NULL
            warning("level.two.var not the same length/nrow as X: level.two.var is ignored.")
         } 
         if(any(is.na(level.two.var))) {
            level.two.var <- NULL
            warning("level.two.var contains missing value(s): level.two.var is ignored.")
         } 
         if(type.z == "Z") {
            type.z <- "WB"
            warning("type.z has been changed to 'WB' to enable testing in multilevel data.")
         }
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
      put <- NULL
      for (lb in lowerbound){ 
         no <- search.normal(X, lb, alpha, StartSet, verbose, type.z, test.Hi, level.two.var)
         output <- cbind(output, no)
      }  
   }
   return(output)   
}



#"aisp.old" <- function(X, 
#                       lowerbound =.3, 
#                       search = "normal", 
#                       alpha = .05, 
#                       StartSet = FALSE,
#                       popsize = 20, 
#                       maxgens = default.maxgens, 
#                       pxover = 0.5, 
#                       pmutation = 0.1,
#                       verbose = FALSE)
#{
#   X <- check.data(X)
#   output <- NULL
#   clabels <- as.character(lowerbound)
#   rlabels <- dimnames(X)[[2]]
#   params <- c(alpha, pxover, pmutation)
#   cparams <- c("alpha", "pxover", "pmutation")
#   for (i in 1:3){
#      if(!is.numeric(params[i])|is.na(params[i])) stop(cparams[i], " is not numeric")
#      if (params[i] < 0) {warning(paste("Negative ",cparams[i],". ",cparams[i]," is set to 0")); assign(params[i],0)}
#      if (params[i] > 1) {warning(paste(cparams[i]," greater than 1. ",cparams[i]," is set to 1")); assign(params[i],1)}
#   }
#   for (i in 1:length(lowerbound)){
#      if(!is.numeric(lowerbound[i])|is.na(lowerbound[i])) stop(" lowerbound contains non-numeric values")
#      if (lowerbound[i] > 1) {warning(paste("Lower bound greater than 1. Lower bound is set to 1")); lowerbound[i] <- 1}
#   }
#   default.maxgens <- 10^(log2(ncol(X)/5)) * 1000
#   if(is.numeric(popsize)&!is.na(popsize)) popsize <- as.integer(popsize) else stop("popsize is not numeric")
#   if(popsize < 1) stop("popsize is nonpositive")
#   if(is.numeric(maxgens)&!is.na(maxgens)) maxgens <- as.integer(maxgens) else stop("maxgens is not numeric")
#   if(maxgens < 1) stop("maxgens is nonpositive")
#   if(search == "ga"){ a
#      output <- NULL
#      for (lb in lowerbound){ 
#         ga <- search.ga(X, popsize, maxgens, alpha, lb, pxover, pmutation)
#         output <- cbind(output, ga)
#      }
#      dimnames(output) <- list(dimnames(X)[[2]], lowerbound)
#   } else if(search == "extended") {
#      output <- search.extended(verbose) 
#   } else {
#      output <- NULL
#      for (lb in lowerbound){ 
#         no <- search.normal(X, lb, alpha, StartSet, verbose)
#         output <- cbind(output, no)
#      }  
#   }  
#   return(output)   
#}
