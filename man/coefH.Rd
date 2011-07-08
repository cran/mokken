\name{coefH}
\alias{coefH}

\title{
Scalability coefficents H
}

\description{
Computes item-pair scalability coefficents Hij, item scalability coefficents Hi, 
 and scale scalability coefficent H (Loevinger, 1948; Mokken, 1971, pp. 148-153; 
 Molenaar and Sijtsma, 2000, pp. 11-13; Sijtsma and Molenaar, chap. 4; Van der Ark, 2007; 2010),
 plus their standard errors (Kuijpers, Van der Ark, and Croon, 2011; also see Van der Ark, Croon, and Sijtsma, 2008)
}

\usage{
coefH(X, se = TRUE, nice.output = TRUE)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items. 
          Missing values are not allowed}
  \item{se}{Logical: If \code{TRUE}, the standard errors of the scalability coefficients are given}
  \item{nice.output}{Logical: If \code{TRUE}, scalability coefficients and standard errors are combined in an a single object of class \code{noquote}}
} 

\details{
 Does not work if any of the item variances equals zero. 
 Such items should not be used in a test and removed from the data frame.
 
 If \code{nice.output = TRUE} and \code{se = TRUE}, the result is a list of 3 objects of class \code{noquote};
 if \code{nice.output = FALSE} and \code{se = TRUE}, the result is a list of 6 matrices (3 for the scalability coefficients and 3 for the standard errors)+
 if \code{se = FALSE}, the result is a list of 3 matrices (for the scalability coefficients.
 
 Computation of standard errors can be slow for a combination of a large sample size and a large number of items.
}

\value{
 \item{Hij}{scalability coefficients of the item pairs (possibly with standard errors; see details)}
 \item{Hi}{vector containing scalability coefficients of the items (possibly with standard errors; see details)}
 \item{H}{scalability coefficient of the entire scale (possibly with standard error; see details)}
 \item{se.Hij}{standard errors of the scalability coefficients of the item pairs (see details)}
 \item{se.Hi}{standard errors of the scalability coefficients of the items (see details)}
 \item{se.H}{standard error of the scalability coefficient of the entire scale (see details)}
}

\references{

  Kuijpers, R. E., Van der Ark, L. A., and Croon, M. A. (2011). 
  Standard errors and confidence intervals for scalability coefficients in Mokken scale analysis using marginal models.
  Manuscript in preparation.
  
  Loevinger, J. (1948). The technique of homogeneous tests compared with some aspects of `scale analysis' and factor analysis.
  \emph{Psychological Bulletin, 45}, 507-530.
   
  Mokken, R. J. (1971) 
  \emph{A Theory and Procedure of Scale Analysis}. 
  Berlin, Germany: De Gruyter.
  
  Molenaar, I.W. and Sijtsma, K. (2000) 
  \emph{User's Manual MSP5 for Windows} [Software manual].
  Groningen, The Netherlands: IEC ProGAMMA.
  
  Sijtsma, K., and Molenaar, I. W. (2002)
  \emph{Introduction to nonparametric item response theory}.
  Thousand Oaks, CA: Sage.

  Van der Ark, L. A. (2007). 
  Mokken scale analysis in \pkg{R}. 
  \emph{Journal of Statistical Software, 20} (11), 1-19.
  \url{http://www.jstatsoft.org/v20/i11}

  Van der Ark, L. A. (2010). 
  Getting started with Mokken scale analysis in \pkg{R}. 
  Unpublished manuscript.
  \url{http://cran.r-project.org/web/packages/mokken/vignettes/mokken.pdf}

  Van der Ark, L. A., Croon, M. A., and Sijtsma (2008). 
  Mokken scale analysis for dichotomous items using marginal models.
  \emph{Psychometrika, 73}, 183-208.

 }

\author{
 L. A. van der Ark \email{a.vdark@uvt.nl}
}

\seealso{
\code{\link{coefZ}}, \code{\link{search.normal}}
}

\examples{
data(acl)
Communality <- acl[,1:10]
coefH(Communality)

coefH(Communality, se=FALSE)

}
\keyword{scalability coefficients}
