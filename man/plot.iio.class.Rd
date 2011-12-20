\name{plot.iio.class}
\alias{plot.iio.class}

\title{
Plot iio.class objects
}

\description{
 S3 Method to \code{plot} objects of class iio.class.
 Graphic display of the checks of iio. 
 One graph for each item plotting the estimated item step response functions.
}
\usage{
\method{plot}{iio.class}(x, item.pairs = all.pairs, ask = TRUE, ...)
}

\arguments{
  \item{x}{Object of class iio.class produced by \code{\link{check.iio}}.}
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
The plot function corresponds to method MIIO; each graph plots the estimated item response functions (item rest-score functions) 
for two items. The default graphic device in R may only display the last graph.
}

\value{
 Returns a graph.
}


\references{
  Ligtvoet, R., L. A. van der Ark, J. M. te Marvelde, and K. Sijtsma (2010).
  Investigating an invariant item ordering for polytomously scored items. 
  \emph{Educational and Psychological Measurement, 70}, 578-595.

  Ligtvoet, R., L. A. van der Ark, W. P. Bergsma, and K. Sijtsma (2011).
  Polytomous latent scales for the investigation of the ordering of items.
  \emph{Psychometrika, 76}, 200-216.

  Sijtsma, K., R. R. Meijer, and L. A. van der Ark (2011).
  Mokken scale analysis as time goes by: An update for scaling practitioners.
  \emph{Personality and Individual Differences, 50}, 31-37.

  Van der Ark, L. A. (2007). 
  Mokken scale analysis in \pkg{R}. 
  \emph{Journal of Statistical Software, 20} (11), 1-19.
  \url{http://www.jstatsoft.org/v20/i11}
}


\author{
 L. A. van der Ark \email{a.vdark@uvt.nl}
}


\seealso{
\code{\link{check.iio}}, \code{\link{summary.iio.class}}
}

\examples{
data(acl)
Communality <- acl[,1:10]
iio.list <- check.iio(Communality)
summary(iio.list)
plot(iio.list)
}

\keyword{plot}
\keyword{iio}
