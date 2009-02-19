\name{check.iio}
\alias{check.iio}

\title{
Check of Invariant Item Ordering Using Method Restscore
}

\description{
Returns a list with results from the investigation of invariant item ordering using method restscore. 
For a complete description of this procedure see Ligtvoet et al. (2008).
}
\usage{
check.iio(X, minvi = m * 0.03, minsize = default.minsize, alpha = .05)
}

\arguments{
  \item{X}{matrix or data frame of numeric data 
          containing the responses of \code{nrow(X)} respondents to \code{ncol(X)} items.  
          Missing values are not allowed}
  \item{minvi}{minimum size of a violation that is reported. By default minvi equals .03 times the number of item step response functions (m)}
  \item{minsize}{minimum size of a rest score group.  By default
   \code{minsize} = \eqn{N/10} if \eqn{N \geq 500};
   \code{minsize} = \eqn{N/5 if} \eqn{250 \leq N < 500}; and
   \code{minsize} = \eqn{\max(N/3,50)} if \eqn{N < 250}}
  \item{alpha}{Nominal Type I error for t test. Default alpha - .05} 
}

\value{
  \item{VI}{A matrix: The first column gives the number of violations of IIO according the method restscore. 
  If the number of violations is nonzero, then the item with the largest number of violations is removed.
  If two or more items have the maximum number of violations, then from those items the item producing the lowest value of Loevinger's H is removed.
  The second column shows the number of violations with one item removed, the third column shows the number of violations with two items removed, etc.}
  \item{items.removed}{List of the items removed in chronological order}
  \item{HT}{Coefficient HT for the remaining items}
}

\references{
  Ligtvoet, R., van der Ark, L. A. and Sijtsma (2008) 
  \emph{Investigating an invariant item ordering for polytomously scored items}. 
  Manuscript submitted for publication.
  
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
\code{\link{coefH}}, \code{\link{check.monotonicity}}, \code{\link{check.pmatrix}}, 
\code{\link{check.reliability}}, \code{\link{check.restscore}}
}

\examples{
data(acl)            

Order <- acl[,41:50]               
check.iio(Order)                     

Autonomy <- acl[,91:100]
check.iio(Autonomy)
}
\keyword{univar}
