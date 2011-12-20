\name{plot.restscore.class}
\alias{plot.restscore.class}

\title{
Plot restscore.class objects
}

\description{
 S3 Method to \code{plot} objects of class restscore.class.
 Graphic display of the checks of restscore. 
 One graph for each item pair plotting the estimated item step response functions.
 Intersection of the lines indicates vaiolation of nonintersection,
}
\usage{
\method{plot}{restscore.class}(x, item.pairs = all.pairs, ask = TRUE, ...)
}

\arguments{
  \item{x}{Object of class restscore.class produced by \code{\link{check.restscore}}.}
  \item{item.pairs}{vector containing the numbers of the item pairs for which the results are depicted graphically.
   For example, \code{item.pairs = 1} prints the results for items 1 and 2,
   \code{item.pairs = 2} prints the results for items 1 and 3,
   \code{item.pairs = J} prints the results for items 1 and $J$, and
   \code{item.pairs = J+1} prints the results for items 2 and 3.
   Default the results for all item pairs are depicted.}
  \item{ask}{Boolean. If \code{TRUE} (default), then \code{par("ask"=TRUE)}; i.e., a hard return between subsequent plots is required.
  If \code{FALSE}, then \code{par("ask"=FALSE)}.}
 \item{...}{Optional graphical parameters will be ignored}  
}

\details{
The default graphic device in R may only display the last graph.
The default number of plots can increase rapidly for large numbers of items.
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
  \url{http://www.jstatsoft.org/v20/i11}
}

\author{
 L. A. van der Ark \email{a.vdark@uvt.nl}
}


\seealso{
\code{\link{check.restscore}}, \code{\link{summary.restscore.class}}
}

\examples{
data(acl)
Communality <- acl[,1:10]
restscore.list <- check.restscore(Communality)
plot(restscore.list)
summary(restscore.list)
}
\keyword{plot}
\keyword{restscore}
