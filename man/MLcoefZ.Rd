\name{MLcoefZ}
\alias{MLcoefZ}

\title{
Computation of Z-Values for two-level scalability coefficients
}

\description{
Computes Zij-values of item pairs, Zi-values of items, and Z-value of the entire scale,
which are used to test whether Hij, Hi, and H, respectively (within- and between-rater versions), 
are significantly greater a specified lowerbound using the delta method (Koopman, Zijlstra, & Van der Ark, 2020a).
The test uses range-preserving asymptotic theory (Koopman, Zijlstra, & Van der Ark, 2020b).
}
\usage{
MLcoefZ(X, lowerbound = 0)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X) - 1} items. 
	  The first column of \code{X} is assumed to be a subject column, see \code{?MLcoefH()} for details.
          Missing values are not allowed}
  \item{lowerbound}{Value of the null hypothesis to which the scalability are compared to compute the Z-score (see details),
          0 <= \code{lowerbound} < 1. The default is 0.}
}

\value{
 \item{Zij}{matrix containing the Z-values of the item-pairs}
 \item{Zi}{vector containing Z-values of the items}
 \item{Z}{Z-value of the entire scale}
}

\details{
The range-preserving Z-score for the estimated item-pair coefficient \eqn{Hij} with standard error \eqn{SE(Hij)} is computed as 
\deqn{Z = -(log(1 - Hij) - log(1 - lowerbound)) / (SE(Hij) / (1 - Hij)) (Koopman, Zijlstra, Van der Ark, 2020b)}
Unlike \code{MLcoefH}, standard errors are not provided.
}

\references{
  Koopman, L. Zijlstra, B. J. H, & Van der Ark, L. A. (2020a). 
  \emph{A two-step procedure for scaling multilevel data using Mokken's scalability coefficients}. 
  Manuscript submitted for publication.
  
  Koopman, L. Zijlstra, B. J. H, & Van der Ark, L. A. (2020b). 
  \emph{Range-preserving confidence intervals for scalability coefficients in Mokken scale analysis}. 
  Manuscript submitted for publication.
 }

\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl}
 L. Koopman 
}


\seealso{
\code{\link{MLcoefH}}
}

\examples{
data(SWMD)

# Compute the Z-score using lowerbound 0
MLcoefZ(SWMD)

# Using lowerbound .1
MLcoefZ(SWMD, lowerbound = .1)
}
\concept{scalability coefficient}
