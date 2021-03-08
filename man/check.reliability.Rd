\name{check.reliability}
\alias{check.reliability}

\title{
Computation of reliability statistics
}

\description{
Returns a list of reliability statistics: Molenaar Sijtsma (MS, 1984, 1988) statistic (a.k.a rho; also see, Sijtsma & Molenaar, 1987; Van der Ark, 2010), 
    Cronbach's (1951) alpha, Guttman's (1945) lambda 2, and the latent class reliability coefficient
    (LCRC; Van der Ark, Van der Palm, & Sijtsma, 2011).
}

\usage{
check.reliability(X, MS = TRUE, alpha = TRUE, lambda.2 = TRUE, 
LCRC = FALSE, nclass = nclass.default, irc = FALSE)
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
        By default: the number of items minus 1.}
  \item{irc}{Boolean.If \code{TRUE}, the item-rest correlation (a.k.a. corrected item-total correlation) is computed.}      
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
  \doi{10.1007/BF02310555}
  
  Guttman, L. (1945).
  A basis for analyzing test-retest reliability.
  \emph{Psychometrika, 10},255-282.
 \doi{10.1007/BF02288892}

  Molenaar, I. W., &  Sijtsma, K. (1984).
  Internal consistency and reliability in Mokken's nonparametric item response model.
  \emph{Tijdschrift voor onderwijsresearch, 9}, 257--268.
  Retrieved from \url{https://pure.uvt.nl/ws/portalfiles/portal/1030704/INTERNAL.PDF}

  Molenaar, I. W., & Sijtsma, K. (1988).
  Mokken's approach to reliability estimation extended to multicategory items.
  \emph{Kwantitatieve methoden, 9(28)}, 115-126.
  Retrieved from \url{https://pure.uvt.nl/ws/portalfiles/portal/1030575/MOKKEN__.PDF}
  
  Sijtsma, K., & Molenaar, I. W. (1987).
  Reliability of test scores in nonparametric item response theory.
  \emph{Psychometrika, 52},79-97.
  \doi{10.1007/BF02293957}
  
  Van der Ark, L. A. (2007). 
  Mokken scale analysis in \pkg{R}. 
  \emph{Journal of Statistical Software}.
  \url{https://www.jstatsoft.org/article/view/v020i11}

  Van der Ark, L. A. (2010). 
  Computation of the Molenaar Sijtsma statistic. 
  In A. Fink, B. Lausen, W. Seidel, & A. Ultsch (Eds.),
  \emph{Advances in data analysis, data handling and business intelligence} (pp. 775-784).
  Springer.
  \doi{10.1007/978-3-642-01044-6_7}

  Van der Ark, L. A., Van der Palm, D. W., & Sijtsma, K. (2011).
  A latent class approach to estimating test-score reliability.
  \emph{Applied Psychological Measurement, 35}, 380-392.
  \doi{10.1177/0146621610392911}
  
  Van der Palm, D. W., Van der Ark, L. A. & Sijtsma, K. (2014).
  A flexible latent class approach to estimating test-score reliability.
  \emph{Journal of Educational Measurement, 51}, 339-357.
  \doi{10.1111/jedm.12053}
}

\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl}
}

\seealso{
\code{\link{check.errors}}, 
\code{\link{check.iio}}, 
\code{\link{check.monotonicity}}, 
\code{\link{check.pmatrix}} 
\code{\link{check.restscore}}, 
\code{\link{coefH}}
}

\examples{
data(acl)
Communality <- acl[,1:10]
check.reliability(Communality, LCRC = TRUE)
}
\concept{reliability}
