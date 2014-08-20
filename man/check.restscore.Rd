\name{check.restscore}
\alias{check.restscore}

\title{
Check of Nonintersection Using Method Restscore
}

\description{
Returns a list (of class \code{restscore.class}) with results from the investigation of nonintersection using method restscore (Mokken, 1971; Molenaar and Sijtsma, 2000; Sijtsma and Molenaar, 2002). 
}
\usage{
check.restscore(X, minvi = 0.03, minsize = default.minsize)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items.  
          Missing values are not allowed}
  \item{minvi}{minimum size of a violation that is reported}
  \item{minsize}{minimum size of a rest score group.  By default
   \code{minsize} = \eqn{N/10} if \eqn{N \ge 500};
   \code{minsize} = \eqn{N/5 if} \eqn{250 \le N < 500}; and
   \code{minsize} = max\eqn{(N/3,50)} if \eqn{N < 250}}
}

\details{
The output is often numerous.  
Procedure may be slow for large data sets. 
Functions \code{plot} and  \code{summary} can be used to summarize the output. 
See Van der Ark (2007) for an example. 
}

\value{
  \item{results}{A list with as many components as there are item pairs.  
  Each component itself is also
  a list containing the results of the check of nonintersection using method restscore.}
  \item{I.labels}{The item labels}
  \item{Hi}{The item scalability coefficients}
  \item{m}{The number of answer categories.}
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
  \url{http://www.jstatsoft.org/v20/i11}
}
\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl}
}

\seealso{
\code{\link{check.errors}}, 
\code{\link{check.iio}}, 
\code{\link{check.monotonicity}}, 
\code{\link{check.pmatrix}}, 
\code{\link{check.reliability}} 
\code{\link{coefH}}, 
\code{\link{plot.restscore.class}}, 
\code{\link{summary.restscore.class}}
}

\examples{
data(acl)
Communality <- acl[,1:10]
restscore.list <- check.restscore(Communality)
plot(restscore.list)
summary(restscore.list)
}
\keyword{restscore}
