\name{mokken-package}
\alias{mokken-package}
\alias{mokken}
\docType{package}

\title{
Mokken Scale Analysis in R
}
\description{
Mokken scale analysis is a scaling procedure for both dichotomous and polytomous items.
It consists of an item selection algorithm to partition a set of items into Mokken scales
 and several methods to check the assumptions of two nonparametric item response theory
models: the monotone homogeneity model and the double monotonicity model.
}
\details{
\tabular{ll}{
Package: \tab mokken\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2009-02-19\cr
License: \tab GPL Version 2 or later\cr
}
The package contains principal functions for Mokken scale analysis.
Version 0 was introduced in Van der Ark (2007). 
Version 1 includes estimation of reliability statistics.
In subversions small bugs were repaired. 
Thanks are due 
 to Daniel van der Palm and J. Hendrik Straat for contributing R code; 
 to Patrick Mair, Rudy Ligtvoet, and J. Hendrik Straat for testing the software; 
 to Michael Dewey, Michael Kubovy, Jue Huang, and Na Yang for reporting bugs;
 to Robert J. Mokken for lending his last name.

}
\author{
L. Andries van der Ark
Maintainer: L. Andries van der Ark <a.vdark@uvt.nl>. 
}

\references{
  Mokken, R. J. (1971) 
  \emph{A Theory and Procedure of Scale Analysis}. 
  Berlin, Germany: De Gruyter.
  
  Molenaar, I.W. and Sijtsma, K. (2000) 
  \emph{User's Manual MSP5 for Windows} [Software manual].
  Groningen, The Netherlands: IEC ProGAMMA.
  
  Sijtsma, K, and Molenaar, I. W. (2002)
  \emph{Introduction to nonparametric item response theory}.
  Thousand Oaks, CA: Sage.
  
  Van der Ark, L. A. (2007). 
  Mokken scale analysis in \pkg{R}. 
  \emph{Journal of Statistical Software}.
  \url{http://www.jstatsoft.org}
}

\keyword{package}

\examples{

# Personality test
data(acl)               

# Select the items of the scale Communality
Communality <- acl[,1:10]

# Compute scalability coefficients
coefH(Communality)

# Investigate the assumption of monotonicity
monotonicity.list <- check.monotonicity(Communality)
summary(monotonicity.list)
plot(monotonicity.list)

# Investigate the assumption of non-intersecting ISRFs using method restscore
restscore.list <- check.restscore(Communality)
summary(restscore.list)
plot(restscore.list)

# Investigate the assumption of non-intersecting ISRFs using method pmatrix
pmatrix.list <- check.pmatrix(Communality)
summary(pmatrix.list)
plot(pmatrix.list)

# Investigate the assumption of IIO using method restscore
check.iio(Communality)

# Compute the reliability of the scale
check.reliability(Communality)

# Partition the the scale into mokken scales
aisp(Communality)

}
