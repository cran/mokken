\name{plot.iio.class}
\alias{plot.iio.class}

\title{
Plot iio.class objects
}

\description{
 S3 Method to \code{plot} objects of class iio.class.
 Graphic display of the checks of iio. 
 One graph for each item plotting the estimated item response functions.
}
\usage{
\method{plot}{iio.class}(x, item.pairs = all.pairs, plot.ci = TRUE, color.ci = c("orange","yellow"), alpha.ci = .05, ask = TRUE, ...)
}

\arguments{
  \item{x}{Object of class iio.class produced by \code{\link{check.iio}}.}
  \item{item.pairs}{vector containing the numbers of the item pairs for which the results are depicted graphically.
   For example, \code{item.pairs = 1} prints the results for items 1 and 2,
   \code{item.pairs = 2} prints the results for items 1 and 3,
   \code{item.pairs = J} prints the results for items 1 and $J$, and
   \code{item.pairs = J+1} prints the results for items 2 and 3.
   Default the results for all item pairs are depicted.}
   \item{plot.ci}{Boolean. If \code{TRUE} (default), then confidence envelops are plotted around ISRFs.}              
  \item{color.ci}{Color of the plotted confidence envelops. Defaults are orange for the first item and yellow for the second item.}
  \item{alpha.ci}{Type of plotted (1 - alpha) confidence intervals. By default 95-percent confidence intervals are depicted}
  \item{ask}{Boolean. If \code{TRUE} (default), then \code{par("ask"=TRUE)}; i.e., a hard return between subsequent plots is required.
  If \code{FALSE}, then \code{par("ask"=FALSE)}.}
 \item{...}{Optional graphical parameters will be ignored}  
}
\details{
The plot function corresponds to method MIIO; each graph plots the estimated item response functions (item rest-score functions) for two items. 
For details of the method, see Ligtvoet et al. (2010, 2011); Sijtsma et al. (2012).
For details of the confidence envelopes, see Van der Ark (2012b).
For the implementation in R, see Van der Ark (2012a).
For \code{ask==FALSE}, the default graphic device in R may only display the last graph.
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

  Van der Ark, L. A. (2012a). 
  New developments in Mokken scale analysis in \pkg{R}. 
  \emph{Journal of Statistical Software, 48} (5), 1-27.
  \url{http://www.jstatsoft.org/v20/i11}
  
  Van der Ark, L. A. (2012b). 
  Visualization of nonparametric IRT models. 
  Paper presented at the Annual Meeting of the Psychometric Society, Lincoln NE, July, 2012.

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
