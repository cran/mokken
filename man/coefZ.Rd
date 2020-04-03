\name{coefZ}
\alias{coefZ}

\title{
Computation of Z-Values 
}

\description{
Computes Zij-values of item pairs, Zi-values of items, and Z-value of the entire scale,
which are used to test whether Hij, Hi, and H, respectively, are significantly 
greater than some given lowerbound using the original method Z
(Molenaar and Sijtsma, 2000, pp. 59-62; Sijtsma and Molenaar, p. 40; Van der Ark, 2007; 2010)
or the delta method (Kuijpers, Van der Ark, & Croon, 2013; Koopman, Zijlstra, & Van der Ark, 2020).
The delta method can also handle nested data.

Used in the function \code{aisp}
}
\usage{
coefZ(X, lowerbound = 0, type.se = "delta", level.two.var = NULL)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items. 
          Missing values are not allowed}
  \item{lowerbound}{Value of the null hypothesis to which the scalability are compared to compute the Z-score (see details),
          0 <= \code{lowerbound} < 1. The default is 0.}
  \item{type.se}{Indicates which type of standard errors is used to compute the Z-score:
          "delta": uses standard errors approximated by the delta method (Kuijpers, Van der Ark, Kroon, 2013; Koopman, Zijlstra, Van der Ark, 2020);
          "Z": uses original Z-test (Mokken, 1971; Molenaar and Sijtsma, 2000; Sijtsma and Molenaar, 2002). 
          The default is "delta".}
  \item{level.two.var}{vector of length \code{nrow(X)} or matrix with number of rows equal to \code{nrow(X)} 
          that indicates the level two variable for nested data (Koopman et al., 2020).}
}

\value{
 \item{Zij}{matrix containing the Z-values of the item-pairs}
 \item{Zi}{vector containing Z-values of the items}
 \item{Z}{Z-value of the entire scale}
}

\details{
The Z-score for the estimated item-pair coefficient \eqn{Hij} with standard error \eqn{SE(Hij)} is computed as \deqn{Z = (Hij - lowerbound) / SE(Hij)}
Unlike \code{coefH}, standard errors are not provided.
}

\references{
  Koopman, L. Zijlstra, B. J. H, & Van der Ark, L. A. (2020). 
  \emph{A two-step procedure for scaling multilevel data using Mokken's scalability coefficients}. 
  Manuscript submitted for publication.

  Kuijpers, R. E., Van der Ark, L. A., and Croon, M. A. (2013).
  Standard errors and confidence intervals for scalability coefficients in Mokken scale analysis using marginal models.
  \emph{Sociological Methodology, 43,} 42-69.
 
  Molenaar, I.W. and Sijtsma, K. (2000) 
  \emph{User's Manual MSP5 for Windows} [Software manual].
  Groningen, The Netherlands: IEC ProGAMMA.
  
  Sijtsma, K, and Molenaar, I. W. (2002)
  \emph{Introduction to nonparametric item response theory}.
  Thousand Oaks, CA: Sage.
  
  Van der Ark, L. A. (2007). 
  Mokken scale analysis in \pkg{R}. 
  \emph{Journal of Statistical Software}.
  \url{http://www.jstatsoft.org/v20/i11}

  Van der Ark, L. A. (2010). 
  Getting started with Mokken scale analysis in \pkg{R}. 
  Unpublished manuscript.
  \url{https://sites.google.com/a/tilburguniversity.edu/avdrark/mokken}
 }

\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl}
 L. Koopman 
}


\seealso{
\code{\link{coefH}}, \code{\link{aisp}}
}

\examples{
data(acl)
Communality <- acl[,1:10]

# Compute the Z-score of each coefficient
coefH(Communality)
coefZ(Communality)

# Using lowerbound .3
coefZ(Communality, lowerbound = .3)

# Z-scores for nested data
data(autonomySupport)
scores <- autonomySupport[, -1]
classes <- autonomySupport[, 1]
coefH(scores, level.two.var = classes)
coefZ(scores, level.two.var = classes)


}
\concept{scalability coefficient}
