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
\method{plot}{iio.class}(x, item.pairs = all.pairs, ci = TRUE, alpha = .05, 
       color = c("black", "blue"), transparancy = 20, ask = TRUE, ...)
}

\arguments{
  \item{x}{Object of class iio.class produced by \code{\link{check.iio}}.}
  \item{item.pairs}{vector containing the numbers of the item pairs for which the results are depicted graphically.
   For example, \code{item.pairs = 1} prints the results for items 1 and 2,
   \code{item.pairs = 2} prints the results for items 1 and 3,
   \code{item.pairs = J} prints the results for items 1 and $J$, and
   \code{item.pairs = J+1} prints the results for items 2 and 3.
   Default the results for all item pairs are depicted.}
   \item{ci}{Boolean. If \code{TRUE} (default), then confidence envelops are plotted around IRFs.}              
  \item{alpha}{Type of plotted (1 - alpha) confidence intervals. By default 95-percent confidence intervals are depicted}
  \item{color}{Color of the plotted curves and confidence envelops. Defaults are black for the first item and blue for the second item.}
  \item{transparancy}{Transparancy of the confidence intervals. Higher values result in more opaque colors for the confidence intervals.} 
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
  Koopman, L., Zijlstra, B. J. H., & Van der Ark, L. A. (2023a). 
  Assumptions and Properties of Two-Level Nonparametric Item Response 
  Theory Models. Manuscript submitted for publication.
  
  Koopman, L., Zijlstra, B. J. H., & Van der Ark, L. A. (2023b). 
  Evaluating Model Fit in Two-Level Mokken Scale Analysis. 
  Manuscript submitted for publication.

  Ligtvoet, R., L. A. van der Ark, J. M. te Marvelde, & K. Sijtsma (2010).
  Investigating an invariant item ordering for polytomously scored items. 
  \emph{Educational and Psychological Measurement, 70}, 578-595.
  \doi{10.1177/0013164409355697}

  Ligtvoet, R., L. A. van der Ark, W. P. Bergsma, & K. Sijtsma (2011).
  Polytomous latent scales for the investigation of the ordering of items.
  \emph{Psychometrika, 76}, 200-216.
  \doi{10.1007/s11336-010-9199-8}

  Sijtsma, K., R. R. Meijer, & Van der Ark, L. A. (2011).
  Mokken scale analysis as time goes by: An update for scaling practitioners.
  \emph{Personality and Individual Differences, 50}, 31-37.
  \doi{10.1016/j.paid.2010.08.016} 

  Van der Ark, L. A. (2012). 
  New developjements in Mokken scale analysis in \pkg{R}. 
  \emph{Journal of Statistical Software, 48} (5), 1-27.
  \doi{10.18637/jss.v048.i05}

  
  
  Van der Ark, L. A. (2014).
  Visualizing uncertainty of estimated response functions in nonparametric item response theory.
  In R. E. Millsap, L. A. van der Ark, D. Bolt, & C. M. Woods (Eds.),
  \emph{New developments in quantitative psychology} (pp. 59-68). New York: Springer.
  \doi{10.1007/978-1-4614-9348-8_5}
}


\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl}
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

# Compute two-level fit statistics (Koopman et al., 2023a, 2023b)
data("autonomySupport")
dat <- autonomySupport[, -1]
groups <- autonomySupport[, 1]
autonomyMIIO <- check.iio(dat, item.selection = FALSE, level.two.var = groups)
summary(autonomyMIIO)
plot(autonomyMIIO)


}

\concept{iio}
\concept{miio}
\concept{invariant item ordering}
\concept{manifest invariant item ordering}


