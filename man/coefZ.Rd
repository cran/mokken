\name{coefZ}
\alias{coefZ}

\title{
Computation of Z-Values 
}

\description{
Computes Zij-values of item pairs, Zi-values of items, and Z-value of the entire scale,
which are used to test whether Hij, Hi, and H, respectively, are significantly 
greater than zero.

Used in the function \code{search.normal}
}
\usage{
coefZ(X)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items. 
          Missing values are not allowed}
}

\value{
 \item{Zij}{matrix containing the z-values of the item-pairs}
 \item{Zi}{vector containing z-values of the items}
 \item{Z}{real z-value of the entire scale}
}


\references{
  Mokken, R. J. (1971) 
  \emph{A Theory and Procedure of Scale Analysis}. 
  Berlin, Germany: De Gruyter.
  
  Molenaar, I.W. and Sijtsma, K. (2000) 
  \emph{User's Manual MSP5 for Windows} [Software manual].
  Groningen, The Netherlands: IEC ProGAMMA.
  
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
\code{\link{coefH}}, \code{\link{search.normal}}
}

\examples{
data(acl)
acl.com <- acl[,1:10]
coefH(acl.com)
coefZ(acl.com)
}
\keyword{univar}

