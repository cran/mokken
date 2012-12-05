\name{plot.restscore.class}
\alias{plot.restscore.class}

\title{
Plot restscore.class objects
}

\description{
 S3 Method to \code{plot} objects of class restscore.class.
 Graphic display of the checks of restscore. 
 One graph for each item pair plotting the estimated item step response functions (ISRFs); confidence envelopes are optional.
 Intersections of the lines indicate violations of nonintersection,
}
\usage{
\method{plot}{restscore.class}(x, item.pairs = all.pairs, plot.ci = TRUE, color.ci = c("orange","yellow"), alpha.ci = .05, ask = TRUE, ...)
}

\arguments{
  \item{x}{Object of class restscore.class produced by \code{\link{check.restscore}}.}
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
For details of the method, see Molenaar and Sijtsma (2000) and Sijtsma and Molenaar (2002). 
For details of the confidence envelopes, see Van der Ark (2012)
For the implementation in R, see Van der Ark (2007).
For \code{ask==FALSE}, the default graphic device in R may only display the last graph.
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

  Van der Ark, L. A. (2012). 
  Visualization of nonparametric IRT models. 
  Paper presented at the Annual Meeting of the Psychometric Society, Lincoln NE, July, 2012.
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
