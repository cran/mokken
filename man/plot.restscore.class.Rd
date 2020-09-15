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
\method{plot}{restscore.class}(x, item.pairs = all.pairs, ci = TRUE, alpha = .05, 
              color = c("black", "blue"), transparancy = 20, ask = TRUE, ...)
}

\arguments{
  \item{x}{Object of class restscore.class produced by \code{\link{check.restscore}}.}
  \item{item.pairs}{vector containing the numbers of the item pairs for which the results are depicted graphically.
   For example, \code{item.pairs = 1} prints the results for items 1 and 2,
   \code{item.pairs = 2} prints the results for items 1 and 3,
   \code{item.pairs = J} prints the results for items 1 and $J$, and
   \code{item.pairs = J+1} prints the results for items 2 and 3.
   Default the results for all item pairs are depicted.}
   \item{ci}{Boolean. If \code{TRUE} (default), then confidence envelops are plotted around ISRFs.}              
  \item{alpha}{Type of plotted (1 - alpha) confidence intervals. By default 95-percent confidence intervals are depicted}
  \item{color}{Color of the plotted lines and confidence envelops. Defaults are black for the first item and blue for the second item.}
  \item{transparancy}{Transparancy of the confidence intervals. Higher values result in more opaque colors for the confidence intervals.} 
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
   Molenaar, I.W., & Sijtsma, K. (2000) 
  \emph{User's Manual MSP5 for Windows} [Software manual].
  IEC ProGAMMA.
  
  Sijtsma, K., & Molenaar, I. W. (2002)
  \emph{Introduction to nonparametric item response theory}.
  Sage.
  
  Van der Ark, L. A. (2007). 
  Mokken scale analysis in \pkg{R}. 
  \emph{Journal of Statistical Software}.
  \url{https://www.jstatsoft.org/article/view/v020i11}

  Van der Ark, L. A. (2014).
  Visualizing uncertainty of estimated response functions in nonparametric item response theory.
  In R. E. Millsap, L. A. van der Ark, D. Bolt, & C. M. Woods (Eds.),
  \emph{New developments in quantitative psychology} (pp. 59-68). New York: Springer.
  \url{https://doi.org/10.1007/978-1-4614-9348-8_5}


}

\author{
 L. A. van der Ark \email{L.A.vanderArk@uva.nl}
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
\concept{restscore}
