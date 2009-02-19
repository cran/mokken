\name{plot.pmatrix.class}
\alias{plot.pmatrix.class}

\title{
Plot pmatrix.class objects
}

\description{
 S3 Method to \code{plot} objects of class pmatrix.class.
 Graphic display of the checks of pmatrix. 
 One graph for each item plotting the rows of the P(++) matrix and rows of the P(--) matrix.
 If nonintersection holds the lines in the plots of the P(++) matrix must be nondecreasing
 and the lines in the plots of the P(--) matrix nust be nonincreasing.
}
\usage{
\method{plot}{pmatrix.class}(x, items = all, pmatrix="both", ...)
}

\arguments{
  \item{x}{Object of class pmatrix.class produced by \code{\link{check.pmatrix}}.}
  \item{items}{vector containing the numbers of the item pairs for which the results are depicted graphically.
   Default the results for all items are depicted.}
  \item{pmatrix}{Valid options are \code{"ppp"}, \code{"pmm"}, and \code{"both"}.
   If \code{pmatrix="ppp"}, then the P(++) matrix is plotted,
   if \code{pmatrix="pmm"}, then the P(--) matrix is plotted,
   if \code{pmatrix="both"}, then both the P(++) matrix and P(--) matrix are plotted.}
 \item{...}{Optional graphical parameters will be ignored}  
}

\details{
The default graphic device in R may only display the last graph.

In the plot of the P(++) matrix and the P(--) matrix, the x-axis 
 contains the \eqn{k =(J-1)m} item steps not pertaining to item j in order of popularity (ascending).
Let \eqn{Y_g=1} indicate that the g-th item step has been taken: i.e. \eqn{Xi \geq  y}
and let \eqn{Y_g=0} indicate that the g-th item step has not been taken: i.e. \eqn{Xi <  y}.
The \eqn{m} lines in the plot of the P(++) matrix connect \eqn{P(X_j \geq x, Y_1=1), \dots, P(X_j \geq x, Y_k=1)), x = 1, \dots m}.
The \eqn{m} lines in the plot of the P(--) matrix connect \eqn{P(X_j < x, Y_1=0), \dots, P(X_j < x, Y_k=0)), x = 1, \dots, m}.

If the number of item steps on the x-axis is greater than 10, then the labels are not displayed,
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
\code{\link{check.pmatrix}}, \code{\link{summary.pmatrix.class}}
}

\examples{
data(acl)
Communality <- acl[,1:10]
pmatrix.list <- check.pmatrix(Communality)
plot(pmatrix.list)
summary(pmatrix.list)
}
\keyword{aplot}
