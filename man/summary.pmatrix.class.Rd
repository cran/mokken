\name{summary.pmatrix.class}
\alias{summary.pmatrix.class}

\title{
Summarize pmatrix.class objects
}
\description{
S3 Method for \code{summary} of objects of class pmatrix.class.
Summarize checks of nonintersection using method pmatrix.
}
\usage{
\method{summary}{pmatrix.class}(object, ...)
}

\arguments{
  \item{object}{list produced by \code{\link{check.pmatrix}}}
  \item{...}{Optional parameters will be ignored}  
}

\details{
}

\value{
  \item{ppp.summary.matrix}{Matrix with \code{ncol(X)} rows and 7 columns, 
        showing for each item a summary of the violations in the P(++) matrix.}
  \item{pmm.summary.matrix}{ppp.summary.matrix}{Matrix with \code{ncol(X)} rows and 7 columns, 
        showing for each item a summary of the violations in the P(--) matrix.}
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
\code{\link{check.pmatrix}}, \code{\link{plot.pmatrix.class}}
}

\examples{
data(acl)
Communality <- acl[,1:10]
pmatrix.list <- check.pmatrix(Communality)
plot(pmatrix.list)
summary(pmatrix.list)
}

\keyword{print}
