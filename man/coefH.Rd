\name{coefH}
\alias{coefH}

\title{
Scalability coefficents H
}

\description{
Computes item-pair scalability coefficents Hij, item scalability coefficents Hi, 
 and scale scalability coefficent H
}

\usage{
coefH(X)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items. 
          Missing values are not allowed}
}

\details{
 Does not work if any of the item variances equals zero. 
 Such items should not be used in a test and removed from the data frame.
}

\value{
 \item{Hij}{matrix containing scalability coefficients of the item pairs}
 \item{Hi}{vector containing scalability coefficients of the items}
 \item{H}{scalability coefficient of the entire scale}
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

\author{
 L. A. van der Ark \email{a.vdark@uvt.nl}
}

\note{
}

\seealso{
\code{\link{coefZ}}, \code{\link{search.normal}}
}

\examples{
data(acl)
Communality <- acl[,1:10]
coefH(Communality)
}
\keyword{univar}
