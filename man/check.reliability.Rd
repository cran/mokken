\name{check.reliability}
\alias{check.reliability}

\title{
Compaction of reliability statistics
}

\description{
Returns a list of reliability statistics: Molenaar Sijtsma (MS, 1984, 1988) statistic (a.k.a rho), 
    Cronbach's (1951) alpha, Guttman's (1945) lambda 2, and the latent class reliability coefficient
    (LCRC; Van der Ark, Van der Palm, and Sijtsma, 2011).
}

\usage{
check.reliability(X, MS = TRUE, alpha = TRUE, lambda.2 = TRUE, LCRC = FALSE, nclass = nclass.default)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items. 
          Missing values are not allowed}
  \item{MS}{Boolean. If \code{TRUE}, The MS statistic is computed.}        
  \item{alpha}{Boolean. If \code{TRUE}, Cronbach's alpha is computed.}        
  \item{lambda.2}{Boolean. If \code{TRUE}, Guttman's Lambda 2 is computed.}        
  \item{LCRC}{Boolean. If \code{TRUE}, the LCRC is computed.} 
  \item{nclass}{Integer. Number of latent classes for the computation of LCRC. 
        By default: half the number of items rounded to the nearest integer.}
}

\details{
The computation of LCRC depends on the package \code{poLCA}, which in its turn depends on
 the packages \code{MASS} and \code{scatterplot3d}. 
Computation of the LCRC may be time consuming if the number of latent classes is large.
The optimal number of latent classes should be determined prior to the computation of the LCRC,
 using software for latent class analysis (e.g., the R-package \code{poLCA}).
}

\value{
  \item{MS}{Molenaar Sijtsma statistic (a.k.a. rho).}
  \item{alpha}{Cronbach's alpha}
  \item{lambda.2}{Guttman's Lambda 2}
  \item{LCRC}{LCRC}
}


\references{
  Cronbach, L. (1951).
  Coefficient alpha and the internal structure of tests.
  \emph{Psychometrika, 16}, 297-334.

  Guttman, L. (1945).
  A basis for analyzing test-retest reliability.
  \emph{Psychometrika, 10},255-282.

  Molenaar, I. W. and K. Sijtsma (1984).
  Internal consistency and reliability in Mokken's nonparametric item response model.
  \emph{Tijdschrift voor onderwijsresearch, 9}, 257--268.

  Molenaar, I. W. and K. Sijtsma (1988).
  Mokken's approach to reliability estimation extended to multicategory items.
  \emph{Kwantitatieve methoden, 9(28)}, 115-126.
  
  Sijtsma, K. and I. W. Molenaar (1987).
  Reliability of test scores in nonparametric item response theory.
  \emph{Psychometrika, 52},79-97.
  
  Van der Ark, L. A. (2007). 
  Mokken scale analysis in \pkg{R}. 
  \emph{Journal of Statistical Software}.
  \url{http://www.jstatsoft.org}

  Van der Ark, L. A. (2010). 
  Computation of the Molenaar Sijtsma statistic. 
  In A. Fink, B. Lausen, W. Seidel, and A. Ultsch (Eds.),
  \emph{Advances in data analysis, data handling and business intelligence} (pp.\ 775-784).
  Berlin: Springer.

  Van der Ark, L. A., Van der Palm, D. W., and Sijtsma (2011).
  A latent class approach to estimating test-score reliability.
  \emph{Applied Psychological Measurement}.
  
  
}

\author{
 L. A. van der Ark \email{a.vdark@uvt.nl}
}

\seealso{
\code{\link{coefH}}, \code{\link{check.iio}}, \code{\link{check.monotonicity}}, \code{\link{check.pmatrix}}, \code{\link{check.restscore}} 
}

\examples{
data(acl)
Communality <- acl[,1:10]
check.reliability(Communality, LCRC = TRUE)
}
\keyword{reliability}
