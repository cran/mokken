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
Version: \tab 1.0\cr
Date: \tab 2006-09-29\cr
License: \tab GPL Version 2 or later\cr
}
The package contains principal functions for Mokken scale analysis.
The package is introduced in Van der Ark (2007).
}
\author{
L. Andries van der Ark
Maintainer: L. Andries van der Ark <a.vdark@uvt.nl>
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
\seealso{
}
\examples{

# Personality test
data(acl)               

# Select the items of the scale Communality
acl.com <- acl[,1:10]

# Compute scalability coefficients
coefH(acl.com)

# Investigate the assumption of monotonicity
monotonicity.com <- check.monotonicity(acl.com)
summary(monotonicity.com)
plot(monotonicity.com)

# Investigate the assumption of non-intersecting using method restscore
restscore.com <- check.restscore(acl.com)
summary(restscore.com)
plot(restscore.com)

# Investigate the assumption of non-intersecting using method pmatrix
pmatrix.com <- check.pmatrix(acl.com)
summary(pmatrix.com)
plot(pmatrix.com)

# Partition the entire data set into mokken scales
search.normal(acl.com)
}
