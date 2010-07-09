\name{summary.monotonicity.class}
\alias{summary.monotonicity.class}

\title{
Summarize monotonicity.class objects
}

\description{
S3 Method for \code{summary} of objects of class monotonicity.class.
Summarizes checks of monotonicity
}
\usage{
\method{summary}{monotonicity.class}(object, ...)
}

\arguments{
  \item{object}{list produced by \code{\link{check.monotonicity}}}
  \item{...}{Optional parameters will be ignored}  
}


\value{
 Matrix with \code{ncol(X)} rows and 9 columns, 
 showing for each item a summary of the violations of monotonicity 
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

\seealso{
\code{\link{check.monotonicity}}, \code{\link{plot.monotonicity.class}}
}

\examples{
data(acl)
Communality <- acl[,1:10]
monotonicity.list <- check.monotonicity(Communality)
plot(monotonicity.list)
summary(monotonicity.list)
}

\keyword{print}
