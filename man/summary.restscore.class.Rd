\name{summary.restscore.class}
\alias{summary.restscore.class}

\title{
Summarize restscore.class objects
}
\description{
S3 Method for \code{summary} of objects of class restscore.class.
Summarize checks of nonintersection using method restscore.
}
\usage{
\method{summary}{restscore.class}(object, ...)
}

\arguments{
  \item{object}{list produced by \code{\link{check.restscore}}}
   \item{...}{Optional parameters will be ignored}  

}
\details{
}

\value{
 Matrix with \code{ncol(X)} rows and 9 columns, 
 showing for each item a summary of the violations of nonintersection using method restscore.
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
\code{\link{check.restscore}}, \code{\link{plot.restscore.class}}
}

\examples{
data(acl)
restscore.list <- check.restscore(acl[,1:10])
plot(restscore.list)
summary(restscore.list)
}

\keyword{print}


