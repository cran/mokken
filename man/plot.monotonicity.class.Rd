\name{plot.monotonicity.class}
\alias{plot.monotonicity.class}

\title{
Plot monotonicity.class objects
}

\description{
 S3 Method to \code{plot} objects of class monotonicity.class.
 Graphic display of the checks of monotonicity. 
 One graph for each item plotting the estimated item step response functions.
}
\usage{
\method{plot}{monotonicity.class}(x, items = all.items, ask = TRUE, ...)
}

\arguments{
  \item{x}{Object of class monotonicity.class produced by \code{\link{check.monotonicity}}.}
  \item{items}{vector containing the numbers of the items for which the results are depicted graphically.
                 Default the results for all items are depicted.}
  \item{ask}{Boolean. If \code{TRUE} (default), then \code{par("ask"=TRUE)}; i.e., a hard return between subsequent plots is required.
  If \code{FALSE}, then \code{par("ask"=FALSE)}.}
 \item{...}{Optional graphical parameters will be ignored}  
}
\details{
The default graphic device in R may only display the last graph.
}

\value{
 Returns a graph.
}

\references{
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
\code{\link{check.monotonicity}}, \code{\link{summary.monotonicity.class}}
}

\examples{
data(acl)
Communality <- acl[,1:10]
monotonicity.list <- check.monotonicity(Communality)
plot(monotonicity.list)
summary(monotonicity.list)
}
\keyword{plot}
\keyword{monotonicity}
