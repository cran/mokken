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
\method{plot}{monotonicity.class}(x, items = all.items, curves = "both", ci = TRUE, 
              alpha = .05, color = "black", transparancy = 20, ask = TRUE, ...)
}

\arguments{
  \item{x}{Object of class monotonicity.class produced by \code{\link{check.monotonicity}}.}
  \item{items}{vector containing the numbers of the items for which the results are depicted graphically.
                 Default the results for all items are depicted.}
  \item{curves}{\code{"ISRF"}: Item step response functions (ISRFs) are depicted; \code{"IRF"}: item response functions (IRFs) are depicted;
                \code{"both"}: Both ISRFs and IRFs are depicted. Default is \code{"both"}}               
  \item{ci}{Boolean. If \code{TRUE} (default), then confidence envelops are plotted around IRFs and ISRFs.}              
  \item{alpha}{Type of plotted (1 - alpha) confidence intervals. By default 95-percent confidence intervals are depicted}
  \item{color}{Color of the plotted curves and confidence envelops. Default is black.}
  \item{transparancy}{Transparancy of the confidence intervals. Higher values result in more opaque colors for the confidence intervals.} 
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
  Koopman, L., Zijlstra, B. J. H., & Van der Ark, L. A. (2023a). 
  Assumptions and Properties of Two-Level Nonparametric Item Response 
  Theory Models. Manuscript submitted for publication.
  
  Koopman, L., Zijlstra, B. J. H., & Van der Ark, L. A. (2023b). 
  Evaluating Model Fit in Two-Level Mokken Scale Analysis. 
  Manuscript submitted for publication.
  
  Molenaar, I.W., & Sijtsma, K. (2000) 
  \emph{User's Manual MSP5 for Windows} [Software manual].
  IEC ProGAMMA.
  
  Sijtsma, K., & Molenaar, I. W. (2002)
  \emph{Introduction to nonparametric item response theory}.
  Sage.
  
  Van der Ark, L. A. (2007). 
  Mokken scale analysis in \pkg{R}. 
  \emph{Journal of Statistical Software}.
  \doi{10.18637/jss.v020.i11}
  

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
\code{\link{check.monotonicity}}, \code{\link{summary.monotonicity.class}}
}

\examples{
data(acl)
Communality <- acl[,1:10]
monotonicity.list <- check.monotonicity(Communality)
plot(monotonicity.list)
summary(monotonicity.list)

# Compute two-level fit statistics (Koopman et al., 2023a, 2023b)
data("autonomySupport")
dat <- autonomySupport[, -1]
groups <- autonomySupport[, 1]
autonomyMM <- check.monotonicity(dat, level.two.var = groups)
summary(autonomyMM)
plot(autonomyMM)
}
\concept{manifest monotonicity}
\concept{monotonicity}
