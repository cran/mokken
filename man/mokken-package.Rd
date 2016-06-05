\name{mokken-package}
\alias{mokken-package}
\alias{mokken}
\docType{package}

\title{
Mokken Scale Analysis in R
}
\description{
Mokken scale analysis (Mokken, 1971; Sijtsma and Molenaar, 2002) is a scaling procedure 
for both dichotomous and polytomous items.
It consists of an item selection algorithm to partition a set of items into Mokken scales
 and several methods to check the assumptions of two nonparametric item response theory
models: the monotone homogeneity model and the double monotonicity model.
The output of this R-package resembles the output of the stand-alone program MSP (Molenaar and Sijtsma, 2000).
}
\details{
\tabular{ll}{
Package:    \tab mokken      \cr
Type:       \tab Package     \cr
Version:    \tab 2.8.2       \cr
Date:       \tab 2016-05-21  \cr
License:    \tab GPL Version 2 or later\cr
}
The package contains principal functions for Mokken scale analysis.

The package contains the following data sets
\tabular{ll}{
\code{\link{acl}}                \tab  Scores on a personality checklist. \cr
\code{\link{cavalini}}           \tab  Scores on an inventory on industrial malodor \cr
\code{\link{transreas}}          \tab  Scores on a transitive reasoning test \cr
\code{\link{DS14B}}              \tab  Scores on a Type D test (bootstrap sample) \cr
}

A guide for Mokken scale analysis in R for people who do not know R (Van der Ark, 2010)  is available as a vignette from \url{https://sites.google.com/a/tilburguniversity.edu/avdrark/mokken}.

Thanks are due 
 to Geert H. van Kollenburg, Letty Koopman, Renske E. Kuijpers, Rudy Ligtvoet, Hannah E. M. Oosterhuis, J. Hendrik Straat, and Daniel W. van der Palm for contributing R code; 
 to Geert H. van Kollenburg, Letty Koopman, Rudy Ligtvoet, Patrick Mair, and J. Hendrik Straat for testing the software; 
 to Wijbrandt van Schuur for comments on the vignette;
 to Stephen Cubbellotti, Michael Dewey, Jasmin Durstin, Wilco H. M. Emons, Jue Huang, Michael Kubovy, Ivo Molenaar, Jonathan Rose, Tobias Schlaffer, Klaas Sijtsma, Iris Smits, Jia Jia Syu, Roger Watson, and Na Yang for reporting comments or bugs;
 to Diederick Stoel (ProfitWise) for financial support, 
 to Harrie C. M. Vorst and Pierre Cavalini for permission to use their data;
 to Robert J. Mokken for lending his last name.

Version 0 was introduced in Van der Ark (2007). It included the functions
\tabular{ll}{
\code{\link{coefH}}              \tab Scalability coefficients \cr
\code{\link{coefZ}}              \tab Test statistics for scalability coefficients \cr
\code{\link{check.monotonicity}} \tab Investigate monotonicity assumptions \cr
\code{\link{check.restscore}}    \tab Investigate nonintersection assumption using Method Restscore \cr
\code{\link{check.pmatrix}}      \tab Investigate nonintersection assumption using Method Pmatrix \cr
\code{\link{search.normal}}      \tab Mokken's automated item selection algorithm \cr
}

The following major modifications have been made.
\tabular{ll}{
\code{\link{aisp}}               \tab  More general automated item selection algorithm. \cr
                                 \tab  Function \code{search} has become obsolete (Version 2.0)\cr
\code{\link{check.reliability}}  \tab  Compute reliability coefficients (Version 2.0)\cr
\code{\link{check.iio}}          \tab  Investigate invariant item orderings (Version 2.4) \cr
\code{\link{coefH}}              \tab  Standard errors for scalability coefficients included (Version 2.6)\cr
}

All updates until version 2.7 are described in Van der Ark (2012). The following modifications have been made in Version 2.7 in comparison to previous versions.
\tabular{ll}{
\code{\link{check.errors}}       \tab  Inclusion new function to compute weighted Guttman errors for each person. \cr
\code{\link{check.iio}}          \tab  \code{plot} has been added. \cr
\code{\link{check.monotonicity}} \tab  Computation of number of active pairs for dichotomous items has been corrected. \cr
\code{\link{check.pmatrix}}      \tab  Summary of the results has been corrected. \cr
\code{\link{check.restscore}}    \tab  Code pertaining to IIO has been deleted. The procedure is now equivalent to MSP. \cr
\code{\link{coefH}}              \tab  Option included to compare scalability coefficients across groups\cr
}

The following modifications have been made in Version 2.7.1 in comparison to previous versions.
\tabular{ll}{
\code{\link{mokken}}             \tab  Some legal issues \cr
}

The following modifications have been made in Version 2.7.2 in comparison to previous versions.
\tabular{ll}{
\code{\link{check.iio}}          \tab  Violations of IIO for dichotomous items are now tested using a z-test rather than a t-test. \cr
}

The following modifications have been made in Version 2.7.3 in comparison to previous versions.
\tabular{ll}{
\code{\link{plot.iio.class}}     \tab  Confidence envelopes around estimated response functions \cr
\code{\link{plot.monotonicity.class}}     \tab  Confidence envelopes around estimated response functions \cr
\code{\link{plot.restscore.class}}     \tab  Confidence envelopes around estimated response functions \cr
}

The following modifications have been made in Version 2.8.1 in comparison to previous versions.
\tabular{ll}{
\code{\link{aisp}}     \tab  Startsets have been added \cr
}

The following modifications have been made in Version 2.8.2 in comparison to previous versions.
\tabular{ll}{
\code{\link{recode}}          \tab  New \cr
\code{\link{check.ca}}        \tab  New \cr
\code{\link{check.norms}}     \tab  New \cr
\code{\link{check.errors}}    \tab  Outlier score O+ has been included \cr
}

}
\author{
L. Andries van der Ark
Maintainer: L. Andries van der Ark <L.A.vanderArk@uva.nl>. 
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
  \emph{Journal of Statistical Software, 20}(11), 1-19.
  \url{http://www.jstatsoft.org/v20/i11}

  Van der Ark, L. A. (2010). 
  Getting started with Mokken scale analysis in \pkg{R}. 
  Unpublished manuscript.
  \url{https://sites.google.com/a/tilburguniversity.edu/avdrark/mokken}


  Van der Ark, L. A. (2012). 
  New developments in Mokken scale analysis in \pkg{R}. 
  \emph{Journal of Statistical Software, 48}(5), 1-27.
  \url{http://www.jstatsoft.org/v48/i5}
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

# Investigate the assumption of IIO using method MIIO
iio.list <- check.iio(Communality)
summary(iio.list)
plot(iio.list)

# Compute the reliability of the scale
check.reliability(Communality)

# Partition the the scale into mokken scales
aisp(Communality)

}
