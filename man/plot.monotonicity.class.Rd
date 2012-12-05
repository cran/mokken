\name{plot.monotonicity.class}
\alias{plot.monotonicity.class}

\title{
Plot monotonicity.class objects
}

\description{
 S3 Method to \code{plot} objects of class monotonicity.class.
 Graphic display of the checks of monotonicity. 
 One graph for each item plotting the estimated item step response functions and/or item response function, plus confidence envelopes (Van der Ark, 2012).
}
\usage{
\method{plot}{monotonicity.class}(x, items = all.items, curves = "both", plot.ci = TRUE, color.ci = "orange", alpha.ci = .05, ask = TRUE, ...)
}

\arguments{
  \item{x}{Object of class monotonicity.class produced by \code{\link{check.monotonicity}}.}
  \item{items}{vector containing the numbers of the items for which the results are depicted graphically.
                 Default the results for all items are depicted.}
  \item{curves}{\code{"ISRF"}: Item step response functions (ISRFs) are depicted; \code{"IRF"}: item response functions (IRFs) are depicted;
                \code{"both"}: Both ISRFs and IRFs are depicted. Default is \code{"both"}}               
  \item{plot.ci}{Boolean. If \code{TRUE} (default), then confidence envelops are plotted around IRFs and ISRFs.}              
  \item{color.ci}{Color of the plotted confidence envelops. Default is orange.}
  \item{alpha.ci}{Type of plotted (1 - alpha) confidence intervals. By default 95-percent confidence intervals are depicted}
  \item{ask}{Boolean. If \code{TRUE} (default), then \code{par("ask"=TRUE)}; i.e., a hard return between subsequent plots is required.
  If \code{FALSE}, then \code{par("ask"=FALSE)}.}
 \item{...}{Optional graphical parameters will be ignored}  
}
\details{
For details of the method, see Molenaar and Sijtsma (2000) and Sijtsma and Molenaar (2002). 
For details of the confidence envelopes, see Van der Ark (2012)
For the implementation in R, see Van der Ark (2007).
For \code{curves=="both"}, both plots are plotted simultaneously using \code{layout(matrix(c(1,2)1,2))}.
For \code{ask=="FALSE"}, the default graphic device in R may only display the last graph.
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
