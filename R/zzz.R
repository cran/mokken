## zzz.R (2008-02-08)

##   Library Loading

## This file is part of the R-package `mokken'.

.First.lib <- function(lib, pkg) {
##    require(nlme, quietly = TRUE)
      library.dynam("mokken", pkg, lib)
}
