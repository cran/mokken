\name{MLcoefZ}
\alias{MLcoefZ}

\title{
Computation of Z-Values for two-level scalability coefficients
}

\description{
Computes Zij-values of item pairs, Zi-values of items, and Z-value of the entire scale,
which are used to test whether Hij, Hi, and H, respectively (within- and between-rater versions), 
are significantly greater a specified lowerbound using the delta method (Koopman et al., in press a).
The test uses either Wald-based (\code{WB}) or range-preserving (\code{RP}) asymptotic theory 
(Koopman et al., in press b).
}

\usage{
MLcoefZ(X, lowerbound = 0, type.z = "WB")
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X) - 1} items. 
	  The first column of \code{X} is assumed to be a subject column, see \code{?MLcoefH()} for details.
          Missing values are not allowed}
  \item{lowerbound}{Value of the null hypothesis to which the scalability are compared to compute the z-score (see details),
          0 <= \code{lowerbound} < 1. The default is 0.}
  \item{type.z}{Indicates which type of z-score is computed:
          "WB": Wald-based z-score based on standard errors as approximated by the delta method 
                (Kuijpers et al., 2013; Koopman et al., in press a);
          "RP": Range-preserving z-score, also based on the delta method (Koopman et al., in press b). 
          The default is "WB".}
}

\value{
 \item{Zij}{matrix containing the Z-values of the item-pairs}
 \item{Zi}{vector containing Z-values of the items}
 \item{Z}{Z-value of the entire scale}
}

\details{
For the estimated item-pair coefficient \eqn{Hij} with standard error \eqn{SE(Hij)}, the Z-score  is computed as \deqn{Zij = (Hij - lowerbound) / SE(Hij)}  
if \code{type.z = "WB"}, and the Z-score  is computed as \deqn{Zij = -(log(1 - Hij) - log(1 - lowerbound)) / (SE(Hij) / (1 - Hij))} if \code{type.z = "RP"}  
(Koopman et al., in press b). For the estimate item-scalability coefficients \eqn{Hi} and total-scalbility coefficients \eqn{H} a similar procedure 
is used. Standard errors of the Z-scores are not provided.
}

\references{
   Koopman, L., Zijlstra, B. J. H., & Van der Ark, L. A. (in press a).
   A two-step, test-guided Mokken scale analysis for nonclustered and clustered data. 
   \emph{Quality of Life Research}. (advanced online publication)
   \doi{10.1007/s11136-021-02840-2}

   Koopman, L., Zijlstra, B. J. H., & Van der Ark, L. A. (in press b).
   Range-preserving confidence intervals and significance tests for scalability coefficients in Mokken scale analysis.
   In M. Wiberg, D. Molenaar, J. Gonzalez, & Kim, J.-S. (Eds.),
   \emph{Quantitative Psychology; The 1st Online Meeting of the Psychometric Society, 2020}.
   Springer.
   \doi{10.1007/978-3-030-74772-5_16}
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
